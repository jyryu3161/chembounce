"""Bounded, offline, lineage-only migration of one pinned reconstruction index.

The coordinator imports no chemistry. Its children have separate address-space
caps; the validator independently accounts for every inert source occurrence.
An external watchdog is still required around the entire invocation, including
Python startup and receipt writes. A receipt alone never authorizes adoption.
"""

import argparse
import csv
import ctypes
import hashlib
import importlib.metadata
import json
import os
import pickle
import pickletools
from pathlib import Path
import resource
import signal
import stat
import subprocess
import sys
import time

STARTED = time.monotonic()
SOURCE_FILES = ("__init__.py", "assets.py", "catalog.py", "common.py", "reactions.py", "catalog_migration.py")
CEILINGS = {
    "total_seconds": 1800,
    "work_seconds": 1740,
    "startup_seconds": 120,
    "converter_seconds": 600,
    "validator_seconds": 900,
    "bookkeeping_seconds": 120,
    "term_seconds": 15,
    "kill_seconds": 15,
    "receipt_seconds": 30,
    "address_space_bytes": 16 * 1024**3,
}
ENVIRONMENT = {
    "PATH": "/usr/bin:/bin",
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "CUDA_VISIBLE_DEVICES": "",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "HF_HUB_OFFLINE": "1",
    "TRANSFORMERS_OFFLINE": "1",
    "PYTHONDONTWRITEBYTECODE": "1",
    "PYTHONOPTIMIZE": "0",
    "PYTHONNOUSERSITE": "1",
    "PYTHONHASHSEED": "0",
}
ID_POLICY = "Strip input ID when nonblank, otherwise structure-sha256 (not a verified supplier SKU)"
AVAILABILITY_POLICY = "Snapshot membership only; live stock is not checked"


def _digest(path):
    before = Path(path).stat()
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    after = Path(path).stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise ValueError(f"File changed while hashing: {path}")
    return h.hexdigest()


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _decode(data):
    return json.loads(data, object_pairs_hook=_pairs)


def _json(value):
    return json.dumps(value, sort_keys=True, allow_nan=False, separators=(",", ":")).encode()


def _regular(path):
    if not stat.S_ISREG(Path(path).lstat().st_mode):
        raise ValueError(f"Expected nonsymlink regular file: {path}")


def _save(path, value):
    """Exclusive durable records; never replace a previous attempt or receipt."""
    data = _json(value) + b"\n"
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, "wb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    fd = os.open(Path(path).parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _load_plan(path, expected):
    _regular(path)
    with Path(path).open("rb") as stream:
        data = stream.read(1024 * 1024 + 1)
    if len(data) > 1024 * 1024 or hashlib.sha256(data).hexdigest() != expected:
        raise ValueError("Plan size/checksum mismatch")
    plan = _decode(data)
    if type(plan.get("schema_version")) is not int or plan["schema_version"] != 1:
        raise ValueError("Unsupported migration plan")
    limits = plan["limits"]
    if set(limits) != set(CEILINGS) or any(
        type(limits[k]) is not int or not 0 < limits[k] <= cap for k, cap in CEILINGS.items()
    ):
        raise ValueError("Invalid or excessive migration limits")
    if (
        limits["work_seconds"] + sum(limits[k] for k in ("term_seconds", "kill_seconds", "receipt_seconds"))
        > limits["total_seconds"]
    ):
        raise ValueError("Work and cleanup must fit the shared total deadline")
    if set(plan["source"]["sha256"]) != set(SOURCE_FILES):
        raise ValueError("Source pins must cover exactly the migration dependencies")
    for section, key in (("input", "path"), ("source", "root"), ("parent", "path")):
        p = Path(plan[section][key])
        if not p.is_absolute() or p != p.resolve():
            raise ValueError(f"Require a single resolved {section} path")
    output = Path(plan["output_dir"])
    if (
        not output.is_absolute()
        or output.parent != Path(plan["parent"]["path"])
        or output.name in {"", ".", ".."}
    ):
        raise ValueError("Output must be a direct child of the pinned parent")
    if Path(plan["input"]["path"]).suffix not in {".pkl", ".pickle"}:
        raise ValueError("Migration requires one reconstruction pickle")
    for key in ("source", "snapshot", "smiles_column", "id_column"):
        if not isinstance(plan["labels"][key], str) or not plan["labels"][key].strip():
            raise ValueError("Missing exact conversion labels")
    if type(plan["input"]["bytes"]) is not int or plan["input"]["bytes"] <= 0:
        raise ValueError("Invalid input size")
    if set(plan["labels"]) != {"source", "snapshot", "smiles_column", "id_column"}:
        raise ValueError("Unexpected conversion labels")
    if type(plan["baseline"]["unique_records"]) is not int or plan["baseline"]["unique_records"] <= 0:
        raise ValueError("Invalid baseline unique count")
    watchdog = Path(plan["watchdog"]["path"])
    if not watchdog.is_absolute() or watchdog != watchdog.resolve():
        raise ValueError("Watchdog must have a resolved absolute path")
    _regular(watchdog)
    if _digest(watchdog) != plan["watchdog"]["sha256"]:
        raise ValueError("Watchdog bytes changed")
    return plan


def _directory(plan, owner=None):
    parent = Path(plan["parent"]["path"])
    if parent != parent.resolve() or not stat.S_ISDIR(parent.lstat().st_mode):
        raise ValueError("Output parent is not a resolved nonsymlink directory")
    info = parent.stat()
    if [info.st_dev, info.st_ino] != plan["parent"]["device_inode"]:
        raise ValueError("Output parent binding changed")
    output = Path(plan["output_dir"])
    if owner is not None:
        info = output.lstat()
        if not stat.S_ISDIR(info.st_mode) or [info.st_dev, info.st_ino] != owner["directory_device_inode"]:
            raise ValueError("Owned output directory changed")
    return output


def _reserve(plan):
    output = _directory(plan)
    fd = os.open(output.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        parent = os.fstat(fd)
        if [parent.st_dev, parent.st_ino] != plan["parent"]["device_inode"]:
            raise ValueError("Output parent changed before reservation")
        os.mkdir(output.name, mode=0o700, dir_fd=fd)
        child = os.open(output.name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd)
        try:
            info = os.fstat(child)
        finally:
            os.close(child)
    finally:
        os.close(fd)
    binding = [info.st_dev, info.st_ino]
    _directory(plan, {"directory_device_inode": binding})
    return output, binding


def _bindings(plan):
    """All expensive binding checks run inside a bounded child phase."""
    if any(os.environ.get(key) != value for key, value in ENVIRONMENT.items()):
        raise ValueError("Child offline/thread environment differs")
    if len(os.sched_getaffinity(0)) != 1:
        raise ValueError("Child must be bound to one allowed CPU")
    root = Path(plan["source"]["root"])
    if Path(__file__).resolve().parent != root / "chembounce_synthesis":
        raise ValueError("Loaded migration source is outside the pinned checkout")
    for name, expected in plan["source"]["sha256"].items():
        path = root / "chembounce_synthesis" / name
        _regular(path)
        if _digest(path) != expected:
            raise ValueError(f"Source checksum mismatch: {name}")
    for arg, expected in (("HEAD", plan["source"]["commit"]), ("HEAD^{tree}", plan["source"]["tree"])):
        actual = subprocess.check_output(["git", "rev-parse", arg], cwd=root, text=True, timeout=5).strip()
        if actual != expected:
            raise ValueError("Source revision/tree changed")
    if subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True, timeout=5).strip():
        raise ValueError("Source checkout is not clean")
    runtime = plan["runtime"]
    if (
        sys.executable != runtime["executable"]
        or str(Path(sys.executable).resolve()) != runtime["resolved_executable"]
    ):
        raise ValueError("Interpreter path changed")
    if (
        _digest(Path(sys.executable).resolve()) != runtime["executable_sha256"]
        or sys.version != runtime["python_version"]
    ):
        raise ValueError("Interpreter version/bytes changed")
    if {k: importlib.metadata.version(k) for k in ("rdkit", "numpy")} != runtime["packages"]:
        raise ValueError("Normalization dependency versions changed")
    metadata = {
        name: hashlib.sha256(importlib.metadata.distribution(name).read_text("METADATA").encode()).hexdigest()
        for name in ("rdkit", "numpy")
    }
    if metadata != runtime["package_metadata_sha256"]:
        raise ValueError("Normalization distribution metadata changed")
    from rdkit import rdBase

    if rdBase.rdkitVersion != runtime["rdkit_runtime_version"]:
        raise ValueError("RDKit runtime changed")
    path = Path(plan["input"]["path"])
    _regular(path)
    if path.stat().st_size != plan["input"]["bytes"] or _digest(path) != plan["input"]["sha256"]:
        raise ValueError("Input bytes changed")
    return {
        "source_sha256": dict(plan["source"]["sha256"]),
        "input_sha256": plan["input"]["sha256"],
        "runtime": runtime,
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "environment": {key: os.environ[key] for key in ENVIRONMENT},
    }


def _outputs(output):
    result = {}
    for name in ("catalog.tsv", "catalog.tsv.manifest.json"):
        p = output / name
        _regular(p)
        result[name] = _digest(p)
    return result


class _Record:
    def __new__(cls, *args, **kwargs):
        return object.__new__(cls)

    def __init__(self, *args, **kwargs):
        pass

    def __setstate__(self, state):
        self.state = state


class _SourceUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        return type("InertRecord", (_Record,), {"module": module, "name": name})

    def persistent_load(self, pid):
        raise ValueError("Persistent source references are unsupported")


def source_records(path):
    """Independent inert schema extraction, exactly one payload and no tail."""
    with Path(path).open("rb") as stream:
        for op, _, _ in pickletools.genops(stream):
            if op.name in {"EXT1", "EXT2", "EXT4", "PERSID", "BINPERSID"}:
                raise ValueError(f"Unsupported source opcode: {op.name}")
        if stream.read(1):
            raise ValueError("Trailing source bytes after the first pickle payload")
        stream.seek(0)
        data = _SourceUnpickler(stream).load()
        if stream.read(1):
            raise ValueError("Unreturned source tail")
    if not isinstance(data, _Record) or (data.module, data.name) not in {
        ("synllama.chem.fpindex", "FingerprintIndex"),
        ("synllama.chem.matrix", "ReactantReactionMatrix"),
    }:
        raise ValueError("Unsupported source index schema")
    state = getattr(data, "state", None)
    if not isinstance(state, dict):
        raise ValueError("Unsupported source state")
    molecules = state.get("_molecules", state.get("_reactants"))
    if not isinstance(molecules, (tuple, list)):
        raise ValueError("Source molecule sequence is missing")
    for item in molecules:
        if (
            not isinstance(item, _Record)
            or (item.module, item.name) != ("synllama.chem.mol", "Molecule")
            or not isinstance(getattr(item, "state", None), str)
        ):
            raise ValueError("Unsupported source molecule state; input is incomplete")
        yield item.state


def validate_output(plan, output):
    """Full raw-to-row oracle; no converter reader, search or fingerprints."""
    from .common import canonical
    from .catalog import Catalog

    before = _outputs(output)
    manifest_data = (output / "catalog.tsv.manifest.json").read_bytes()
    if hashlib.sha256(manifest_data).hexdigest() != before["catalog.tsv.manifest.json"]:
        raise ValueError("Manifest changed before parsing")
    manifest = _decode(manifest_data)
    rows = {}
    with (output / "catalog.tsv").open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["id", "smiles", "source", "snapshot", "availability"]:
            raise ValueError("Unexpected TSV columns")
        for row in reader:
            if None in row or any(v is None for v in row.values()) or row["id"] in rows:
                raise ValueError("Malformed or duplicate TSV row")
            rows[row["id"]] = row
    labels = plan["labels"]
    origins, rejections, typed_rejections = {}, [], []
    total = accepted = 0
    path = plan["input"]["path"]
    for coordinate, raw in enumerate(source_records(path)):
        total += 1
        try:
            smiles = canonical(raw)
        except ValueError as exc:
            rejections.append({"input_index": 0, "file": path, "record": coordinate, "reason": str(exc)})
            typed_rejections.append(
                {"record_index": coordinate, "error_type": type(exc).__name__, "reason": str(exc)}
            )
            continue
        ident = "structure-sha256:" + hashlib.sha256(smiles.encode()).hexdigest()
        expected = {
            "id": ident,
            "smiles": smiles,
            "source": labels["source"],
            "snapshot": labels["snapshot"],
            "availability": "catalogued",
        }
        if rows.get(ident) != expected:
            raise ValueError(f"Raw source coordinate {coordinate} does not match its exact TSV ID/structure")
        origins.setdefault(ident, []).append(
            {"input_index": 0, "record_index": coordinate, "identifier_origin": "generated_structure_sha256"}
        )
        accepted += 1
    expected_manifest = {
        "schema_version": 2,
        "inputs": [{"path": path, "sha256": plan["input"]["sha256"]}],
        "output_sha256": before["catalog.tsv"],
        "records": len(rows),
        "rejected_records": rejections,
        **labels,
        "record_index_policy": "zero_based_logical_records",
        "record_origins": origins,
        "id_policy": ID_POLICY,
        "availability_policy": AVAILABILITY_POLICY,
    }
    # JSON equality preserves scalar types, unlike Python's True == 1.
    if _json(manifest) != _json(expected_manifest) or set(rows) != set(origins):
        raise ValueError("Manifest is not the complete exact source-coordinate mapping")
    if accepted + len(rejections) != total:
        raise ValueError("Source occurrence partition is incomplete")
    baseline = plan["baseline"]
    if (
        before["catalog.tsv"] != baseline["output_sha256"]
        or len(rows) != baseline["unique_records"]
        or _json([{k: r[k] for k in ("record", "reason")} for r in rejections])
        != _json(baseline["rejected_records"])
    ):
        raise ValueError("Migration differs from the frozen legacy baseline; preserve for review")
    catalog = Catalog(output / "catalog.tsv")
    if set(catalog.blocks) != set(rows):
        raise ValueError("Loaded Catalog ID coverage differs")
    for ident, row in rows.items():
        expected = {
            **row,
            "catalog_sha256": before["catalog.tsv"],
            "identifier_provenance": {
                "status": "recorded",
                "manifest_sha256": before["catalog.tsv.manifest.json"],
                "source_records": [
                    {**origin, "input_path": path, "input_sha256": plan["input"]["sha256"]}
                    for origin in origins[ident]
                ],
            },
            "supplier_sku_verified": False,
            "current_stock_verified": False,
        }
        if _json(catalog.evidence(ident)) != _json(expected):
            raise ValueError("Catalog evidence differs from the complete raw-source mapping")
    if _outputs(output) != before or _digest(path) != plan["input"]["sha256"]:
        raise ValueError("Input or output changed during validation")
    return {
        "complete_source": True,
        "N": total,
        "A": accepted,
        "R": len(rejections),
        "U": len(rows),
        "duplicate_accepted_occurrences": accepted - len(rows),
        "rejections": typed_rejections,
        "output_sha256": before,
        "supplier_sku_verified": False,
        "current_stock_verified": False,
    }


def worker(plan, plan_sha, phase):
    output = Path(plan["output_dir"])
    owner = _decode((output / "owner.json").read_bytes())
    _directory(plan, owner)
    before = _bindings(plan)
    if phase == "startup":
        result = {"bindings": before}
    elif phase == "converter":
        from .assets import convert_catalog

        for name in ("catalog.tsv", "catalog.tsv.manifest.json", "catalog.tsv.manifest.json.tmp"):
            if os.path.lexists(output / name):
                raise ValueError("Conversion output name is already present")
        convert_catalog([plan["input"]["path"]], output / "catalog.tsv", **plan["labels"])
        result = {"output_sha256": _outputs(output)}
    elif phase == "validator":
        converted = _decode((output / "converter.json").read_bytes())
        if _outputs(output) != converted["output_sha256"]:
            raise ValueError("Output changed between converter and validator")
        result = validate_output(plan, output)
    elif phase == "bookkeeping":
        validated = _decode((output / "validator.json").read_bytes())
        if _outputs(output) != validated["output_sha256"]:
            raise ValueError("Output changed after validation")
        result = {"output_sha256": validated["output_sha256"], "bindings": before}
    else:
        raise ValueError("Unsupported migration phase")
    if _bindings(plan) != before:
        raise ValueError("Bindings changed during child phase")
    _directory(plan, owner)
    usage = resource.getrusage(resource.RUSAGE_SELF)
    _save(
        output / f"{phase}.json",
        {
            "phase": phase,
            "status": "complete",
            "plan_sha256": plan_sha,
            "usage": {
                "user_cpu_seconds": usage.ru_utime,
                "system_cpu_seconds": usage.ru_stime,
                "max_rss_kib_linux": usage.ru_maxrss,
            },
            **result,
        },
    )


def _get_subreaper():
    libc = ctypes.CDLL(None, use_errno=True)
    value = ctypes.c_int()
    if libc.prctl(37, ctypes.byref(value), 0, 0, 0) != 0:  # Linux PR_GET_CHILD_SUBREAPER
        raise OSError(ctypes.get_errno(), "Cannot read migration subreaper state")
    return value.value


def _subreaper(value=1):
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(36, value, 0, 0, 0) != 0:  # Linux PR_SET_CHILD_SUBREAPER
        raise OSError(ctypes.get_errno(), "Cannot set migration subreaper state")


def _group_exists(group):
    try:
        os.killpg(group, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def _reap(proc):
    proc.poll()
    if proc.returncode is None:
        return
    while True:
        try:
            pid, _ = os.waitpid(-proc.pid, os.WNOHANG)
        except ChildProcessError:
            break
        if not pid:
            break


def cleanup(proc, limits, total_deadline):
    stages, errors = [], []
    for signum, seconds in (
        (signal.SIGTERM, limits["term_seconds"]),
        (signal.SIGKILL, limits["kill_seconds"]),
    ):
        _reap(proc)
        if not _group_exists(proc.pid):
            break
        try:
            os.killpg(proc.pid, signum)
        except ProcessLookupError:
            pass
        except OSError as exc:
            errors.append(
                {"signal": signal.Signals(signum).name, "error_type": type(exc).__name__, "reason": str(exc)}
            )
        end = min(time.monotonic() + seconds, total_deadline - limits["receipt_seconds"])
        while time.monotonic() < end:
            _reap(proc)
            if not _group_exists(proc.pid):
                break
            time.sleep(min(0.05, max(0, end - time.monotonic())))
        stages.append(signal.Signals(signum).name)
    _reap(proc)
    return {
        "signals": stages,
        "group_gone": not _group_exists(proc.pid),
        "direct_child_reaped": proc.returncode is not None,
        "returncode": proc.returncode,
        "errors": errors,
        "scope": "Owned process group only; an unkillable group is unsafe, never success",
    }


def _child_limits(plan, cpu):
    limit = plan["limits"]["address_space_bytes"]
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    os.sched_setaffinity(0, {cpu})


def run(plan_path, plan_sha, *, started=None):
    start = STARTED if started is None else started
    plan = _load_plan(plan_path, plan_sha)
    limits = plan["limits"]
    total_deadline, work_deadline = start + limits["total_seconds"], start + limits["work_seconds"]
    # Exclusive mkdir anchored to the pinned parent inode, including dangling links.
    output, directory_binding = _reserve(plan)
    owner = {
        "schema_version": 1,
        "plan_path": str(plan_path),
        "plan_sha256": plan_sha,
        "directory_device_inode": directory_binding,
        "pid": os.getpid(),
        "started_monotonic": start,
        "work_deadline_monotonic": work_deadline,
        "total_deadline_monotonic": total_deadline,
        "limits": limits,
    }
    phases, active, result = [], None, "failed"
    interrupted = []
    previous = {}
    previous_affinity = None
    previous_subreaper = None
    try:
        for signum in (signal.SIGINT, signal.SIGTERM):
            previous[signum] = signal.signal(signum, lambda sig, frame: interrupted.append(sig))
        _save(output / "owner.json", owner)
        previous_subreaper = _get_subreaper()
        _subreaper()
        previous_affinity = os.sched_getaffinity(0)
        cpu = min(previous_affinity)
        os.sched_setaffinity(0, {cpu})
        owner["cpu"] = cpu
        env = dict(ENVIRONMENT, PYTHONPATH=plan["source"]["root"])
        for phase in ("startup", "converter", "validator", "bookkeeping"):
            _directory(plan, owner)
            now = time.monotonic()
            if interrupted or now >= work_deadline:
                result = "interrupted" if interrupted else "timeout"
                break
            deadline = min(work_deadline, now + limits[f"{phase}_seconds"])
            if phase == "startup":
                deadline = min(deadline, start + limits["startup_seconds"])
            log = output / f"{phase}.log"
            with log.open("xb") as stream:
                active = subprocess.Popen(
                    [
                        plan["watchdog"]["path"],
                        "--signal=TERM",
                        f"--kill-after={limits['kill_seconds']}s",
                        f"{max(0.001, deadline - time.monotonic())}s",
                        plan["runtime"]["executable"],
                        "-m",
                        "chembounce_synthesis.catalog_migration",
                        "_worker",
                        "--plan",
                        str(plan_path),
                        "--plan-sha256",
                        plan_sha,
                        "--phase",
                        phase,
                    ],
                    cwd=plan["source"]["root"],
                    env=env,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                    preexec_fn=lambda: _child_limits(plan, cpu),
                )
                while active.poll() is None and not interrupted and time.monotonic() < deadline:
                    time.sleep(min(0.05, max(0, deadline - time.monotonic())))
            expired = time.monotonic() >= deadline
            clean = cleanup(active, limits, total_deadline)
            entry = {
                "phase": phase,
                "returncode": active.returncode,
                "deadline_reached": expired,
                "cleanup": clean,
                "elapsed_seconds": time.monotonic() - now,
            }
            phases.append(entry)
            active = None
            if not clean["group_gone"] or not clean["direct_child_reaped"]:
                result = "unsafe_cleanup"
                break
            if interrupted or expired or entry["returncode"] != 0:
                result = "interrupted" if interrupted else "timeout" if expired else "failed"
                break
            evidence = _decode((output / f"{phase}.json").read_bytes())
            if evidence.get("status") != "complete" or evidence.get("plan_sha256") != plan_sha:
                raise ValueError("Missing or mismatched completed phase receipt")
        else:
            result = "validated"
        error = None
    except BaseException as exc:
        error = {"error_type": type(exc).__name__, "reason": str(exc)}
    finally:
        if active is not None:
            clean = cleanup(active, limits, total_deadline)
            phases.append({"phase": "exception_cleanup", "cleanup": clean})
            if not clean["group_gone"] or not clean["direct_child_reaped"]:
                result = "unsafe_cleanup"
        subreaper = {"before": previous_subreaper, "after": None, "restoration": "not_acquired"}
        try:
            if previous_subreaper is not None:
                if result == "unsafe_cleanup":
                    subreaper["restoration"] = "retained_for_unsafe_cleanup"
                else:
                    _subreaper(previous_subreaper)
                    subreaper["restoration"] = "restored_previous"
                subreaper["after"] = _get_subreaper()
                if result != "unsafe_cleanup" and subreaper["after"] != previous_subreaper:
                    raise RuntimeError("Migration subreaper restoration did not retain previous state")
        except (OSError, RuntimeError) as exc:
            subreaper["restoration"] = "failed"
            subreaper["error"] = {"error_type": type(exc).__name__, "reason": str(exc)}
            if result == "validated":
                result = "failed"
        if time.monotonic() >= work_deadline and result == "validated":
            result = "timeout"
        partial = {}
        for name in ("catalog.tsv", "catalog.tsv.manifest.json", "catalog.tsv.manifest.json.tmp"):
            p = output / name
            if os.path.lexists(p):
                s = p.lstat()
                partial[name] = {"bytes": s.st_size, "regular": stat.S_ISREG(s.st_mode)}
        receipt = {
            **owner,
            "status": result,
            "error": error,
            "phases": phases,
            "partial_outputs": partial,
            "elapsed_seconds": time.monotonic() - start,
            "elapsed_seconds_scope": "Before terminal serialization/fsync; the external watchdog covers the complete invocation",
            "received_signals": interrupted,
            "subreaper": subreaper,
            "complete_lifecycle_requires": "pre-write terminal snapshot plus coordinator exit 0 after writing before the original work and total deadlines, plus external watchdog exit 0 within total ceiling; independent result review and adoption remain separate",
            "cost_scope": "Per-child self usage when available in phase receipts; missing counters are unavailable, not zero",
            "final_binding_status": "verified"
            if result == "validated"
            else "unavailable_or_failed; see completed phase receipts",
        }
        if time.monotonic() >= total_deadline:
            receipt["status"] = result = "timeout"
        try:
            _directory(plan, owner)
            _save(output / "terminal.json", receipt)
        finally:
            if previous_affinity is not None:
                os.sched_setaffinity(0, previous_affinity)
            for signum, handler in previous.items():
                signal.signal(signum, handler)
    # The exclusive receipt is a pre-write snapshot. Do not rewrite it to claim
    # knowledge of its own future fsync; the exit decision covers that work.
    finished = time.monotonic()
    return (
        0
        if result == "validated"
        and not interrupted
        and finished < work_deadline
        and finished < total_deadline
        else 3
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "_worker"))
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--phase", choices=("startup", "converter", "validator", "bookkeeping"))
    args = parser.parse_args(argv)
    try:
        if args.action == "run":
            return run(args.plan, args.plan_sha256)
        plan = _load_plan(args.plan, args.plan_sha256)
        try:
            worker(plan, args.plan_sha256, args.phase)
        except BaseException as exc:
            _save(
                Path(plan["output_dir"]) / f"{args.phase}-failure.json",
                {
                    "status": "incomplete",
                    "phase": args.phase,
                    "error_type": type(exc).__name__,
                    "reason": str(exc),
                    "complete_source": False,
                    "plan_sha256": args.plan_sha256,
                },
            )
            raise
        return 0
    except BaseException as exc:
        print(
            json.dumps({"status": "failed", "error_type": type(exc).__name__, "reason": str(exc)}),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
