"""Fetch or package exact raw-source snapshots without changing their provenance.

Supports the retained ChEMBL ``{sources: [...]}`` and RCSB/RuSH list manifests.
Mutable URLs are not a substitute for the recorded bytes: a download with a
different SHA-256 fails rather than silently updating the scientific dataset.
"""

import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys
from urllib.parse import urlparse
from urllib.request import Request, urlopen
import zipfile

from benchmarks.research.protocol import digest, write_json


def _records(data):
    rows = data if isinstance(data, list) else data.get("sources") if isinstance(data, dict) else None
    if not isinstance(rows, list) or not rows:
        raise ValueError("Source manifests require a nonempty list or sources list")
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("Source manifest entries must be objects")
        for key in ("path", "url", "sha256"):
            if not isinstance(row.get(key), str) or not row[key].strip():
                raise ValueError(f"Source manifest entries require {key}")
        if len(row["sha256"]) != 64 or any(c not in "0123456789abcdef" for c in row["sha256"]):
            raise ValueError("Source SHA-256 must be 64 lowercase hexadecimal characters")
        if urlparse(row["url"]).scheme not in {"http", "https"}:
            raise ValueError("Source URLs must use HTTP or HTTPS")
        if "bytes" in row and (
            isinstance(row["bytes"], bool) or not isinstance(row["bytes"], int) or row["bytes"] < 0
        ):
            raise ValueError("Source byte counts must be nonnegative integers")
    return rows


def _local_source(manifest, row, source_root):
    raw = Path(row["path"])
    candidates = [raw] if raw.is_absolute() else [source_root / raw, manifest.parent / raw]
    candidates.append(manifest.parent / raw.name)
    for candidate in candidates:
        if candidate.is_file():
            if digest(candidate) != row["sha256"]:
                raise ValueError(f"Retained source hash mismatch: {candidate}")
            if "bytes" in row and candidate.stat().st_size != row["bytes"]:
                raise ValueError(f"Retained source byte count mismatch: {candidate}")
            return candidate.resolve()
    raise ValueError(f"Retained source absent: {row['path']}; use fetch for the exact recorded URL")


def _download(row, destination, timeout, max_bytes):
    temporary = destination.with_suffix(destination.suffix + ".part")
    size, hasher = 0, hashlib.sha256()
    try:
        request = Request(row["url"], headers={"User-Agent": "ChemBounce-research-snapshot/1"})
        with urlopen(request, timeout=timeout) as response, temporary.open("wb") as stream:
            if urlparse(response.geturl()).scheme not in {"http", "https"}:
                raise ValueError("Source redirect left HTTP/HTTPS")
            while chunk := response.read(1024 * 1024):
                size += len(chunk)
                if size > max_bytes:
                    raise ValueError(f"Source download exceeds max_bytes={max_bytes}")
                hasher.update(chunk)
                stream.write(chunk)
        if hasher.hexdigest() != row["sha256"]:
            raise ValueError(
                f"Downloaded source hash mismatch: {row['url']}. The URL may have changed; "
                "use the retained snapshot instead of treating new bytes as the original dataset."
            )
        if "bytes" in row and size != row["bytes"]:
            raise ValueError("Downloaded source byte count differs from the manifest")
        temporary.replace(destination)
        return size
    finally:
        temporary.unlink(missing_ok=True)


def _safe_artifact(root, relative):
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError("Snapshot artifacts must have relative paths")
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Snapshot artifact escapes its snapshot directory")
    return path


def verify_snapshot(directory):
    root = Path(directory).resolve()
    snapshot = json.loads((root / "snapshot.json").read_text())
    if snapshot.get("schema_version") != 1 or snapshot.get("status") != "complete":
        raise ValueError("Source snapshot is incomplete or unsupported")
    records = snapshot.get("artifacts")
    if not isinstance(records, list) or not records:
        raise ValueError("Source snapshot has no artifacts")
    for row in records:
        path = _safe_artifact(root, row["artifact"])
        if not path.is_file() or digest(path) != row["sha256"]:
            raise ValueError(f"Snapshot artifact hash mismatch or missing file: {row['artifact']}")
        if path.stat().st_size != row["bytes"]:
            raise ValueError(f"Snapshot byte count mismatch: {row['artifact']}")
    return {
        "status": "verified",
        "artifact_count": len(records),
        "total_bytes": sum(row["bytes"] for row in records),
    }


def archive_snapshot(directory, filename):
    root, destination = Path(directory).resolve(), Path(filename).resolve()
    verify_snapshot(root)
    snapshot = json.loads((root / "snapshot.json").read_text())
    # Export only verified artifacts and their manifest, never incidental files
    # subsequently placed in the snapshot directory.
    included = sorted({"snapshot.json", *(row["artifact"] for row in snapshot["artifacts"])})
    if destination.is_relative_to(root):
        raise ValueError("Archive output must be outside the snapshot directory")
    if destination.exists():
        raise ValueError(f"Archive output already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for relative in included:
                path = _safe_artifact(root, relative)
                info = zipfile.ZipInfo(relative, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                archive.writestr(info, path.read_bytes())
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return {"path": str(destination), "sha256": digest(destination), "bytes": destination.stat().st_size}


def acquire_snapshot(
    manifests, output, *, mode="package", source_root=None, timeout=30, max_bytes=1_073_741_824
):
    if mode not in {"package", "fetch"}:
        raise ValueError("Acquisition mode must be package or fetch")
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not 0 < timeout <= 60:
        raise ValueError("Download timeout must be in (0, 60] seconds")
    if isinstance(max_bytes, bool) or not isinstance(max_bytes, int) or max_bytes < 1:
        raise ValueError("max_bytes must be a positive integer per source")
    if not isinstance(manifests, (list, tuple)) or not manifests:
        raise ValueError("At least one source manifest is required")
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Snapshot output directory must be empty")
    source_root = Path(source_root or Path.cwd()).resolve()
    prepared = []
    for index, filename in enumerate(manifests):
        manifest = Path(filename).resolve()
        data = json.loads(manifest.read_text())
        rows = _records(data)
        group = f"{index:02d}-{manifest.parent.name}-{manifest.stem}"
        sources = (
            [_local_source(manifest, row, source_root) for row in rows]
            if mode == "package"
            else [None] * len(rows)
        )
        prepared.append((manifest, data, rows, group, sources))
    output.mkdir(parents=True, exist_ok=True)
    snapshot = {
        "schema_version": 1,
        "status": "incomplete",
        "mode": mode,
        "artifacts": [],
        "source_manifests": [],
        "provenance_policy": "Original manifests and all original metadata/URLs/hashes retained; only portable path fields rewritten",
        "license_policy": "Existing license metadata/files retained verbatim; no license or redistribution permission inferred",
    }

    def record_file(path, role, **metadata):
        relative = str(path.relative_to(output))
        if any(row["artifact"] == relative for row in snapshot["artifacts"]):
            return
        snapshot["artifacts"].append(
            {
                "artifact": relative,
                "role": role,
                "sha256": digest(path),
                "bytes": path.stat().st_size,
                **metadata,
            }
        )

    write_json(output / "snapshot.json", snapshot)
    try:
        for manifest, original, rows, group, sources in prepared:
            directory = output / group
            directory.mkdir()
            original_copy = directory / "original-manifest.json"
            shutil.copyfile(manifest, original_copy)
            record_file(original_copy, "original_manifest", original_path=str(manifest))
            rewritten = copy.deepcopy(original)
            rewritten_rows = _records(rewritten)
            for row, target_row, source in zip(rows, rewritten_rows, sources):
                filename = row["sha256"][:16] + "-" + Path(row["path"]).name
                destination = directory / filename
                if not destination.exists():
                    if mode == "package":
                        shutil.copyfile(source, destination)
                    else:
                        _download(row, destination, timeout, max_bytes)
                if digest(destination) != row["sha256"]:
                    raise ValueError(f"Acquired source changed during acquisition: {destination}")
                target_row["path"] = filename
                role = (
                    "license"
                    if any(term in filename.lower() for term in ("license", "copying"))
                    else "raw_source"
                )
                record_file(destination, role, original_record=copy.deepcopy(row), acquisition=mode)
                write_json(output / "snapshot.json", snapshot)
            portable = directory / "manifest.json"
            write_json(portable, rewritten)
            record_file(portable, "portable_manifest")
            snapshot["source_manifests"].append(
                {
                    "original": str(original_copy.relative_to(output)),
                    "portable": str(portable.relative_to(output)),
                }
            )
        snapshot["status"] = "complete"
        write_json(output / "snapshot.json", snapshot)
        snapshot["verification"] = verify_snapshot(output)
        write_json(output / "snapshot.json", snapshot)
    except Exception as exc:
        snapshot["failure"] = str(exc)
        write_json(output / "snapshot.json", snapshot)
        raise
    return snapshot


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for mode in ("package", "fetch"):
        command = subparsers.add_parser(mode)
        command.add_argument("--manifest", action="append", required=True)
        command.add_argument("--output", required=True)
        command.add_argument("--source-root")
        command.add_argument("--timeout", type=float, default=30)
        command.add_argument("--max-bytes", type=int, default=1_073_741_824)
        command.add_argument("--archive")
    verify = subparsers.add_parser("verify")
    verify.add_argument("snapshot")
    args = parser.parse_args(argv)
    try:
        if args.command == "verify":
            result = verify_snapshot(args.snapshot)
        else:
            snapshot = acquire_snapshot(
                args.manifest,
                args.output,
                mode=args.command,
                source_root=args.source_root,
                timeout=args.timeout,
                max_bytes=args.max_bytes,
            )
            result = {"snapshot": str(Path(args.output).resolve()), **snapshot["verification"]}
            if args.archive:
                result["archive"] = archive_snapshot(args.output, args.archive)
        print(json.dumps(result, sort_keys=True))
        return 0
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
