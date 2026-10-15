#!/usr/bin/env python3
"""Post hoc membership audit of a partial upstream-labelled SMILES sample.

This scratch helper never changes benchmark panels or model/campaign files.
Membership is not proof of exact-checkpoint training; nonmembership is not a
clean holdout. Private row snapshots allow replay after live runs.json changes.
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
import tarfile


def sha(data):
    return hashlib.sha256(data).hexdigest()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def extract(args):
    out = Path(args.output)
    if out.exists() or out.with_suffix(out.suffix + ".provenance.json").exists():
        raise ValueError("Refusing to overwrite extracted sample or provenance")
    archive = Path(args.archive)
    actual = digest(archive)
    if actual != args.archive_sha256:
        raise ValueError("Archive SHA256 mismatch")
    print("Archive SHA256 verified; reading only requested regular member", flush=True)
    data = None
    with tarfile.open(archive, "r|gz") as stream:
        for member in stream:
            if member.name != args.member:
                continue
            if not member.isfile() or member.size > 10_000_000:
                raise ValueError("Requested sample is not a small regular file")
            data = stream.extractfile(member).read()
            if len(data) != member.size:
                raise ValueError("Truncated archive member")
            break
    if data is None:
        raise ValueError("Requested archive member absent")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("xb") as stream:
        stream.write(data)
    provenance = {
        "archive": str(archive.resolve()), "archive_sha256": actual,
        "archive_bytes": archive.stat().st_size, "member": args.member,
        "sample_bytes": len(data), "sample_sha256": sha(data),
        "source_url": args.source_url, "extracted_utc": now(),
        "operation": "regular-member byte read only; no pickle deserialization",
        "scope": "upstream train-labelled evaluation sample; exact-checkpoint membership unverified",
    }
    save(out.with_suffix(out.suffix + ".provenance.json"), provenance)
    print(json.dumps(provenance, indent=2))


def project_sources(spec, base):
    """Read each JSON once, retaining its exact bytes' digest and chemical rows."""
    snapshots, rows, cache = [], [], {}

    def read(relative):
        path = (base / relative).resolve()
        if path not in cache:
            payload = path.read_bytes()
            cache[path] = (json.loads(payload), {
                "path": str(path), "sha256": sha(payload), "bytes": len(payload),
                "read_utc": now(),
            })
        return cache[path]

    def add(source, role, item, smiles, **extra):
        rows.append({"source": source, "role": role, "smiles": smiles, **item, **extra})

    for source in spec:
        name, kind = source["id"], source["kind"]
        data, provenance = read(source["path"])
        meta = {**provenance, "id": name, "kind": kind}
        before = len(rows)
        if kind == "panel":
            lookup = {r["id"]: r for r in data["compounds"]}
            for row in data["compounds"]:
                add(name, "panel_" + row["split"] + "_" + row["label"],
                    {"record_id": row["id"], "target_id": row["target_id"]}, row["smiles"])
            for query in data["queries"]:
                row = lookup[query["compound_id"]]
                add(name, "query", {"query_id": query["id"], "record_id": row["id"],
                    "target_id": row["target_id"]}, row["smiles"])
        elif kind == "controls":
            panel, panel_provenance = read(source["panel"])
            lookup = {r["id"]: r for r in panel["compounds"]}
            queries = {q["id"]: lookup[q["compound_id"]] for q in panel["queries"]}
            meta["linked_panel"] = panel_provenance
            for i, pair in enumerate(data):
                context = {"row_index": i, "query_id": pair["query_id"],
                           "target_id": pair["target_id"]}
                control, query = lookup[pair["compound_id"]], queries[pair["query_id"]]
                if control["target_id"] != pair["target_id"] or query["target_id"] != pair["target_id"]:
                    raise ValueError("Control target lookup mismatch")
                add(name, "control_" + pair["label"], context, control["smiles"],
                    record_id=control["id"])
                add(name, "query", context, query["smiles"], record_id=query["id"])
        elif kind == "matched_pairs":
            for pair in data:
                context = {key: pair[key] for key in
                           ("pair_id", "query_id", "context_id", "assay_id", "target_id")}
                add(name, "control_" + pair["label"], context, pair["control_smiles"],
                    record_id=pair["control_record_id"])
                add(name, "query", context, pair["query_smiles"], record_id=pair["query_record_id"])
        elif kind == "matched_panel":
            lookup = {r["id"]: r for r in data["records"]}
            for query in data["queries"]:
                qrow = lookup[query["record_id"]]
                context = {key: query[key] for key in ("context_id", "assay_id", "target_id")}
                add(name, "query", context, qrow["smiles"], query_id=query["id"], record_id=qrow["id"])
                for rid in query["control_record_ids"]:
                    control = lookup[rid]
                    if control["context_id"] != query["context_id"]:
                        raise ValueError("Matched context mismatch")
                    add(name, "control_" + control["label"], context, control["smiles"],
                        query_id=query["id"], record_id=rid)
        elif kind == "runs":
            meta["run_count"] = len(data)
            meta["run_status_counts"] = dict(Counter(r["status"] for r in data))
            for run in data:
                context = {key: run[key] for key in ("method_id", "query_id", "target_id", "seed", "status")}
                for rank, candidate in enumerate(run.get("response", {}).get("candidates", []), 1):
                    add(name, "returned_candidate", context, candidate.get("smiles"), rank=rank)
                for rank, candidate in enumerate(run.get("evaluation", {}).get("candidates", []), 1):
                    add(name, "evaluated_candidate", context, candidate.get("smiles"), rank=rank)
        elif kind == "generated_scores":
            for row in data:
                context = {k: row.get(k) for k in ("method_id", "query_id", "target_id", "generation_seed", "rank")}
                add(name, "rescored_generated_candidate", context, row.get("smiles"))
        else:
            raise ValueError(f"Unsupported source kind: {kind}")
        meta["projected_rows"] = len(rows) - before
        snapshots.append(meta)
    return {"created_utc": now(), "sources": snapshots, "rows": rows}


def make_normalizer():
    from chembounce_synthesis.common import canonical
    from chembounce_synthesis.identity import chemical_identity

    @lru_cache(maxsize=50000)
    def normalize(smiles):
        out = {"exact": None, "family": None, "exact_error": None, "family_error": None}
        try:
            out["exact"] = canonical(smiles)
        except ValueError as exc:
            out["exact_error"] = str(exc)
            out["family_error"] = "exact_validation_failed"
            return out
        try:
            out["family"] = chemical_identity(smiles).family_smiles
        except ValueError as exc:
            out["family_error"] = str(exc)
        return out
    return normalize


def compare(sample_lines, rows, normalize):
    sample_index = {"exact": defaultdict(list), "family": defaultdict(list)}
    sample_errors = []
    nonempty = []
    for line_index, line in enumerate(sample_lines, 1):
        if not line.strip():
            continue
        # Official .smi member contains one SMILES per line; never guess columns.
        smiles = line.strip()
        nonempty.append(smiles)
        identity = normalize(smiles)
        for level in sample_index:
            if identity[level] is not None:
                sample_index[level][identity[level]].append(line_index)
        if identity["exact_error"] or identity["family_error"]:
            sample_errors.append({"line": line_index, "input_sha256": sha(smiles.encode()), **identity})
    groups = defaultdict(list)
    matches, errors = [], []
    for i, row in enumerate(rows):
        smiles = row.get("smiles")
        if smiles is not None and not isinstance(smiles, str):
            raise ValueError("SMILES values must be strings or null")
        identity = normalize(smiles)
        enriched = {**row, **identity}
        groups[(row["source"], row["role"])].append(enriched)
        groups[("ALL", row["role"])].append(enriched)
        groups[("ALL", "ALL")].append(enriched)
        hit = {level: sample_index[level].get(identity[level], []) for level in sample_index}
        if any(hit.values()):
            matches.append({"row_index": i, **enriched, "sample_line_indices": hit})
        if identity["exact_error"] or identity["family_error"]:
            errors.append({"row_index": i, **enriched})
        if (i + 1) % 2000 == 0:
            print(f"Audited {i + 1}/{len(rows)} projected rows", file=sys.stderr, flush=True)
    summaries = []
    for (source, role), items in sorted(groups.items()):
        summary = {"source": source, "role": role, "rows": len(items),
                   "raw_distinct": len({item["smiles"] for item in items})}
        for level in sample_index:
            valid = [item[level] for item in items if item[level] is not None]
            distinct = set(valid)
            summary[level] = {
                "assessed_rows": len(valid), "unassessed_rows": len(items) - len(valid),
                "distinct": len(distinct), "duplicate_rows": len(valid) - len(distinct),
                "matching_rows": sum(value in sample_index[level] for value in valid),
                "matching_distinct": len(distinct.intersection(sample_index[level])),
            }
        summaries.append(summary)
    sample_summary = {"lines": len(sample_lines), "nonempty_rows": len(nonempty),
                      "empty_lines": len(sample_lines) - len(nonempty), "raw_distinct": len(set(nonempty))}
    for level, index in sample_index.items():
        count = sum(map(len, index.values()))
        sample_summary[level] = {"assessed_rows": count, "unassessed_rows": len(nonempty) - count,
                                 "distinct": len(index), "duplicate_rows": count - len(index)}
    return {"sample": sample_summary, "groups": summaries,
            "matches": matches, "sample_normalization_failures": sample_errors,
            "panel_normalization_failures": errors}


def audit(args):
    out = Path(args.output)
    if out.exists():
        raise ValueError("Refusing to overwrite existing audit directory")
    if bool(args.spec) == bool(args.snapshot):
        raise ValueError("Supply exactly one of --spec or --snapshot")
    if args.snapshot:
        snapshot_bytes = Path(args.snapshot).read_bytes()
        snapshot = json.loads(snapshot_bytes)
    else:
        spec_path = Path(args.spec).resolve()
        spec_bytes = spec_path.read_bytes()
        snapshot = project_sources(json.loads(spec_bytes), spec_path.parent)
        snapshot["spec_sha256"] = sha(spec_bytes)
        snapshot_bytes = (json.dumps(snapshot, indent=2, sort_keys=True) + "\n").encode()
    sample_path = Path(args.sample)
    sample_bytes = sample_path.read_bytes()
    if sha(sample_bytes) != args.sample_sha256:
        raise ValueError("Sample SHA256 mismatch")
    result = compare(sample_bytes.decode().splitlines(), snapshot["rows"], make_normalizer())
    from rdkit import rdBase
    from chembounce_synthesis import identity, common
    result.update({
        "schema": "partial_train_labelled_sample_overlap_v1", "completed_utc": now(),
        "interpretation": [
            "Post hoc audit against an upstream train-labelled evaluation sample only.",
            "Exact SynLlama-1B-2M-91rxns checkpoint training membership remains unverified.",
            "Nonmembership in this small sample does not establish a clean holdout.",
            "No molecular similarity or activity scores influenced source or membership selection.",
            "Groups count appearances; ALL/ALL distinct identities remove repeated roles/targets/seeds/panels.",
            "Live campaign inputs are bounded by each retained source read_utc and SHA256; later runs are not included.",
        ],
        "sample_sha256": sha(sample_bytes), "sample_bytes": len(sample_bytes),
        "snapshot_sha256": sha(snapshot_bytes), "sources": snapshot["sources"],
        "identity_policy": identity.IDENTITY_POLICY, "rdkit_version": rdBase.rdkitVersion,
        "code_sha256": {str(Path(p).resolve()): digest(p) for p in (__file__, identity.__file__, common.__file__)},
        "license_scope": "Private review artifact; no source molecule redistribution permission inferred.",
    })
    out.mkdir(parents=True)
    (out / "input-rows.private.json").write_bytes(snapshot_bytes)
    save(out / "overlap.json", result)
    public = {k: v for k, v in result.items() if k not in
              ("matches", "sample_normalization_failures", "panel_normalization_failures")}
    public["matched_appearance_count"] = len(result["matches"])
    public["sample_normalization_failure_count"] = len(result["sample_normalization_failures"])
    public["panel_normalization_failure_count"] = len(result["panel_normalization_failures"])
    save(out / "overlap-aggregate.json", public)
    print(json.dumps({"output": str(out), "sample": result["sample"],
                      "all": next(x for x in result["groups"] if x["source"] == x["role"] == "ALL")}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    extraction = sub.add_parser("extract")
    for flag in ("archive", "archive-sha256", "member", "source-url", "output"):
        extraction.add_argument("--" + flag, required=True)
    extraction.set_defaults(function=extract)
    auditing = sub.add_parser("audit")
    for flag in ("sample", "sample-sha256", "output"):
        auditing.add_argument("--" + flag, required=True)
    auditing.add_argument("--spec")
    auditing.add_argument("--snapshot")
    auditing.set_defaults(function=audit)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
