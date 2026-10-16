"""Tiny inert inputs only: source accounting, failure preservation and ownership."""

import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import pickle
import signal
import subprocess
import sys
import time
import types

import pytest
from rdkit import rdBase

from chembounce_synthesis import catalog_migration as migration
from chembounce_synthesis.assets import convert_catalog
from chembounce_synthesis.common import canonical


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_pickle(path, monkeypatch, strings=("CCO", "OCC", "CCN", "CC.O", "*")):
    names = ("synllama", "synllama.chem", "synllama.chem.mol", "synllama.chem.fpindex")
    with monkeypatch.context() as patch:
        for name in names:
            patch.setitem(sys.modules, name, types.ModuleType(name))
        mol = type("Molecule", (), {"__module__": names[2], "__getstate__": lambda self: self.smiles})
        index = type("FingerprintIndex", (), {"__module__": names[3]})
        sys.modules[names[2]].Molecule, sys.modules[names[3]].FingerprintIndex = mol, index
        data = index()
        data._molecules = []
        for smiles in strings:
            item = mol()
            item.smiles = smiles
            data._molecules.append(item)
        path.write_bytes(pickle.dumps(data))
    return path


@pytest.fixture(scope="module")
def source_repo(tmp_path_factory):
    import shutil

    root = tmp_path_factory.mktemp("migration-source")
    package = root / "chembounce_synthesis"
    package.mkdir()
    original = Path(migration.__file__).parent
    for name in migration.SOURCE_FILES:
        shutil.copyfile(original / name, package / name)
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    subprocess.run(["git", "add", "."], cwd=root, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.name=Fixture",
            "-c",
            "user.email=fixture@example.invalid",
            "commit",
            "-qm",
            "fixture sources",
        ],
        cwd=root,
        check=True,
    )
    return root


@pytest.fixture
def plan(tmp_path, source_repo, monkeypatch):
    path = save_pickle(tmp_path / "index.pkl", monkeypatch)
    labels = {
        "source": "synthetic test index; no supplier",
        "snapshot": "fixture-v1",
        "smiles_column": "smiles",
        "id_column": "id",
    }
    # Independent CSV baseline with known expected two unique structures.
    rows = []
    for smiles in ("CCO", "CCN"):
        ident = "structure-sha256:" + hashlib.sha256(smiles.encode()).hexdigest()
        rows.append(
            {
                "id": ident,
                "smiles": smiles,
                "source": labels["source"],
                "snapshot": labels["snapshot"],
                "availability": "catalogued",
            }
        )
    baseline = tmp_path / "expected.tsv"
    with baseline.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["id", "smiles", "source", "snapshot", "availability"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["id"]))
    rejected = []
    for coordinate, raw in ((3, "CC.O"), (4, "*")):
        with pytest.raises(ValueError) as exc:
            canonical(raw)
        rejected.append({"record": coordinate, "reason": str(exc.value)})
    git = lambda arg: subprocess.check_output(["git", "rev-parse", arg], cwd=source_repo, text=True).strip()
    info = tmp_path.stat()
    limits = dict(migration.CEILINGS)
    limits.update(
        total_seconds=45,
        work_seconds=40,
        startup_seconds=10,
        converter_seconds=10,
        validator_seconds=10,
        bookkeeping_seconds=10,
        term_seconds=1,
        kill_seconds=1,
        receipt_seconds=3,
    )
    executable = Path(sys.executable)
    return {
        "schema_version": 1,
        "source": {
            "root": str(source_repo),
            "commit": git("HEAD"),
            "tree": git("HEAD^{tree}"),
            "sha256": {
                name: digest(source_repo / "chembounce_synthesis" / name) for name in migration.SOURCE_FILES
            },
        },
        "runtime": {
            "executable": str(executable),
            "resolved_executable": str(executable.resolve()),
            "executable_sha256": digest(executable.resolve()),
            "python_version": sys.version,
            "packages": {k: importlib.metadata.version(k) for k in ("rdkit", "numpy")},
            "package_metadata_sha256": {
                k: hashlib.sha256(
                    importlib.metadata.distribution(k).read_text("METADATA").encode()
                ).hexdigest()
                for k in ("rdkit", "numpy")
            },
            "rdkit_runtime_version": rdBase.rdkitVersion,
        },
        "watchdog": {"path": "/usr/bin/timeout", "sha256": digest("/usr/bin/timeout")},
        "input": {"path": str(path), "sha256": digest(path), "bytes": path.stat().st_size},
        "parent": {"path": str(tmp_path), "device_inode": [info.st_dev, info.st_ino]},
        "output_dir": str(tmp_path / "attempt"),
        "labels": labels,
        "baseline": {"output_sha256": digest(baseline), "unique_records": 2, "rejected_records": rejected},
        "limits": limits,
    }


def write_plan(plan):
    path = Path(plan["parent"]["path"]) / "plan.json"
    path.write_text(json.dumps(plan))
    return path, digest(path)


def converted(plan):
    output = Path(plan["output_dir"])
    output.mkdir()
    convert_catalog([plan["input"]["path"]], output / "catalog.tsv", **plan["labels"])
    return output


def test_complete_source_partition_and_no_fingerprints(plan, monkeypatch):
    from chembounce_synthesis.catalog import _FingerprintMapping

    output = converted(plan)
    monkeypatch.setattr(
        _FingerprintMapping, "__getitem__", lambda *args: pytest.fail("fingerprint requested")
    )
    result = migration.validate_output(plan, output)
    assert {k: result[k] for k in ("N", "A", "R", "U")} == {"N": 5, "A": 3, "R": 2, "U": 2}
    assert result["duplicate_accepted_occurrences"] == 1
    assert result["complete_source"] is True
    assert not result["supplier_sku_verified"] and not result["current_stock_verified"]
    assert [r["error_type"] for r in result["rejections"]] == ["ValueError", "ValueError"]


@pytest.mark.parametrize(
    "defect",
    [
        "missing",
        "extra",
        "swapped",
        "repeated",
        "accepted_rejected",
        "bool_index",
        "bool_count",
        "supplied",
        "wrong_input",
        "wrong_source",
        "wrong_snapshot",
        "wrong_column",
        "wrong_rejection",
        "extra_field",
    ],
)
def test_all_lineage_fields_checked_against_full_source(plan, defect):
    output = converted(plan)
    path = output / "catalog.tsv.manifest.json"
    manifest = json.loads(path.read_text())
    ids = sorted(manifest["record_origins"])
    duplicate = next(key for key, values in manifest["record_origins"].items() if len(values) == 2)
    origin = manifest["record_origins"][duplicate][0]
    if defect == "missing":
        manifest["record_origins"][duplicate].pop()
    elif defect == "extra":
        manifest["record_origins"][duplicate].append({**origin, "record_index": 5})
    elif defect == "swapped":
        manifest["record_origins"][ids[0]], manifest["record_origins"][ids[1]] = (
            manifest["record_origins"][ids[1]],
            manifest["record_origins"][ids[0]],
        )
    elif defect == "repeated":
        manifest["record_origins"][duplicate].append(dict(origin))
    elif defect == "accepted_rejected":
        origin["record_index"] = 3
    elif defect == "bool_index":
        origin["input_index"] = False
    elif defect == "bool_count":
        manifest["schema_version"] = True
    elif defect == "supplied":
        origin["identifier_origin"] = "input_provided"
    elif defect == "wrong_input":
        manifest["inputs"][0]["path"] += ".other"
    elif defect in {"wrong_source", "wrong_snapshot"}:
        manifest[defect[6:]] = "other"
    elif defect == "wrong_column":
        manifest["id_column"] = "supplier"
    elif defect == "wrong_rejection":
        manifest["rejected_records"][0]["reason"] = "invented"
    else:
        manifest["unreviewed"] = True
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="complete exact"):
        migration.validate_output(plan, output)


@pytest.mark.parametrize("defect", ["structure", "id", "availability", "source", "extra_row"])
def test_wrong_raw_joins_fail_even_with_rebound_digests(plan, defect):
    output = converted(plan)
    path = output / "catalog.tsv"
    with path.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if defect == "extra_row":
        rows.append({**rows[0], "id": "extra"})
    else:
        key = "smiles" if defect == "structure" else defect
        rows[0][key] = "CO" if defect == "structure" else "wrong"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    mp = output / "catalog.tsv.manifest.json"
    manifest = json.loads(mp.read_text())
    manifest["output_sha256"] = digest(path)
    mp.write_text(json.dumps(manifest))
    plan["baseline"]["output_sha256"] = digest(path)
    with pytest.raises(ValueError, match="source coordinate|complete exact"):
        migration.validate_output(plan, output)


@pytest.mark.parametrize("tail", [b"extra", pickle.dumps(None), b"\x00"])
def test_trailing_source_never_silently_accepted(plan, tail):
    output = converted(plan)
    path = Path(plan["input"]["path"])
    path.write_bytes(path.read_bytes() + tail)
    with pytest.raises(ValueError, match="Trailing"):
        migration.validate_output(plan, output)


@pytest.mark.parametrize(
    "data", [b"\x80\x02\x82\x01.", pickle.dumps({"_molecules": []}), b"\x80\x04", b"Pbad\n."]
)
def test_unsupported_or_truncated_source_is_not_a_rejection(tmp_path, data):
    path = tmp_path / "bad.pkl"
    path.write_bytes(data)
    with pytest.raises((ValueError, EOFError)):
        list(migration.source_records(path))


def test_inert_validator_cannot_execute_reduce(tmp_path):
    marker = tmp_path / "executed"

    class Dangerous:
        def __reduce__(self):
            return os.system, (f"touch {marker}",)

    path = tmp_path / "bad.pkl"
    path.write_bytes(pickle.dumps(Dangerous()))
    with pytest.raises(ValueError, match="schema"):
        list(migration.source_records(path))
    assert not marker.exists()


def test_complete_isolated_coordinator(plan):
    path, pinned = write_plan(plan)
    affinity = os.sched_getaffinity(0)
    assert migration.run(path, pinned, started=time.monotonic()) == 0
    assert os.sched_getaffinity(0) == affinity
    output = Path(plan["output_dir"])
    receipt = json.loads((output / "terminal.json").read_text())
    assert receipt["status"] == "validated"
    assert [p["phase"] for p in receipt["phases"]] == ["startup", "converter", "validator", "bookkeeping"]
    assert all(p["cleanup"]["group_gone"] and p["cleanup"]["direct_child_reaped"] for p in receipt["phases"])
    assert json.loads((output / "validator.json").read_text())["N"] == 5
    assert json.loads((output / "startup.json").read_text())["bindings"]["cpu_affinity"] == [receipt["cpu"]]
    assert (output / "owner.json").stat().st_mode & 0o077 == 0


@pytest.mark.parametrize("defect", ["tail", "source_hash", "runtime", "baseline", "memory"])
def test_owned_failure_preserves_terminal_and_partial_outputs(plan, defect):
    if defect == "tail":
        path = Path(plan["input"]["path"])
        path.write_bytes(path.read_bytes() + b"tail")
        plan["input"].update(sha256=digest(path), bytes=path.stat().st_size)
    elif defect == "source_hash":
        plan["source"]["sha256"]["common.py"] = "0" * 64
    elif defect == "runtime":
        plan["runtime"]["packages"]["rdkit"] = "incorrect"
    elif defect == "baseline":
        plan["baseline"]["output_sha256"] = "0" * 64
    else:
        plan["limits"]["address_space_bytes"] = 1024 * 1024
    path, pinned = write_plan(plan)
    assert migration.run(path, pinned, started=time.monotonic()) == 3
    output = Path(plan["output_dir"])
    receipt = json.loads((output / "terminal.json").read_text())
    assert receipt["status"] == "failed"
    assert all(p["cleanup"]["group_gone"] for p in receipt["phases"])
    if defect in {"tail", "baseline"}:
        assert set(receipt["partial_outputs"]) == {"catalog.tsv", "catalog.tsv.manifest.json"}
        assert json.loads((output / "validator-failure.json").read_text())["complete_source"] is False
    if defect == "memory":
        assert receipt["phases"][0]["returncode"] != 0
        assert not (output / "startup.json").exists()  # Unavailable counters are not fabricated zeros.


@pytest.mark.parametrize("kind", ["directory", "file", "symlink", "dangling"])
def test_existing_output_never_overwritten(plan, kind):
    output = Path(plan["output_dir"])
    if kind == "directory":
        output.mkdir()
    elif kind == "file":
        output.write_text("keep")
    else:
        output.symlink_to(output.parent if kind == "symlink" else output.parent / "missing")
    path, pinned = write_plan(plan)
    with pytest.raises(FileExistsError):
        migration.run(path, pinned, started=time.monotonic())
    assert not (output / "owner.json").exists() if kind == "directory" else os.path.lexists(output)


def test_changed_plan_parent_and_excessive_limits_refused(plan):
    path, pinned = write_plan(plan)
    path.write_text(path.read_text() + " ")
    with pytest.raises(ValueError, match="checksum"):
        migration.run(path, pinned, started=time.monotonic())
    plan["parent"]["device_inode"][1] += 1
    path, pinned = write_plan(plan)
    with pytest.raises(ValueError, match="parent binding"):
        migration.run(path, pinned, started=time.monotonic())
    plan["limits"]["total_seconds"] = 1801
    path, pinned = write_plan(plan)
    with pytest.raises(ValueError, match="excessive"):
        migration.run(path, pinned, started=time.monotonic())


def test_expired_shared_work_deadline_does_not_start_child(plan, monkeypatch):
    path, pinned = write_plan(plan)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **k: pytest.fail("child started after deadline"))
    assert migration.run(path, pinned, started=time.monotonic() - 40) == 3
    receipt = json.loads((Path(plan["output_dir"]) / "terminal.json").read_text())
    assert receipt["status"] == "timeout" and receipt["phases"] == []


def test_owned_term_kill_and_reap_are_bounded():
    migration._subreaper()
    proc = subprocess.Popen(
        [
            sys.executable,
            "-c",
            "import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); print('ready',flush=True); time.sleep(30)",
        ],
        stdout=subprocess.PIPE,
        start_new_session=True,
    )
    assert proc.stdout.readline() == b"ready\n"
    limits = dict(migration.CEILINGS, term_seconds=1, kill_seconds=1, receipt_seconds=1)
    result = migration.cleanup(proc, limits, time.monotonic() + 4)
    assert result["signals"] == ["SIGTERM", "SIGKILL"]
    assert result["group_gone"] and result["direct_child_reaped"]
    assert result["returncode"] == -signal.SIGKILL


def test_unkillable_group_reported_unsafe_not_zero(monkeypatch):
    class Process:
        pid, returncode = 987654321, None

        def poll(self):
            return None

    monkeypatch.setattr(migration, "_group_exists", lambda pid: True)
    monkeypatch.setattr(os, "killpg", lambda *args: None)
    result = migration.cleanup(Process(), dict(migration.CEILINGS), time.monotonic())
    assert not result["group_gone"] and not result["direct_child_reaped"]
    assert result["returncode"] is None


def test_mid_sequence_schema_failure_is_not_counted_as_ordinary_rejection(tmp_path, monkeypatch):
    path = save_pickle(tmp_path / "mid.pkl", monkeypatch, ("CCO", 123, "CCN"))
    stream = migration.source_records(path)
    assert next(stream) == "CCO"
    with pytest.raises(ValueError, match="incomplete"):
        next(stream)


@pytest.mark.parametrize("defect", ["nonzero", "timeout", "interrupt"])
def test_partial_tsv_survives_owned_child_failures(plan, monkeypatch, defect):
    plan["limits"]["converter_seconds"] = 1 if defect == "timeout" else 10
    path, pinned = write_plan(plan)
    original = subprocess.Popen
    output = Path(plan["output_dir"])

    def launch(argv, **kwargs):
        if argv[-1] == "converter":
            code = f"from pathlib import Path; import os,time; Path({str(output / 'catalog.tsv')!r}).write_text('partial TSV'); "
            code += "os._exit(7)" if defect == "nonzero" else "time.sleep(60)"
            process = original([sys.executable, "-c", code], **kwargs)
            if defect == "interrupt":
                # The coordinator's handler is already installed; no subprocess or
                # unrelated process group receives this fixture-only signal.
                limit = time.monotonic() + 3
                while not (output / "catalog.tsv").exists() and time.monotonic() < limit:
                    time.sleep(0.01)
                os.kill(os.getpid(), signal.SIGTERM)
            return process
        return original(argv, **kwargs)

    monkeypatch.setattr(subprocess, "Popen", launch)
    assert migration.run(path, pinned, started=time.monotonic()) == 3
    receipt = json.loads((output / "terminal.json").read_text())
    assert (
        receipt["status"] == {"nonzero": "failed", "timeout": "timeout", "interrupt": "interrupted"}[defect]
    )
    assert set(receipt["partial_outputs"]) == {"catalog.tsv"}
    assert (output / "catalog.tsv").read_text() == "partial TSV"
    assert all(p["cleanup"]["group_gone"] for p in receipt["phases"])
    assert not (output / "validator.json").exists()


@pytest.mark.parametrize("defect", ["input_drift", "output_drift"])
def test_binding_changes_between_children_never_validate(plan, monkeypatch, defect):
    path, pinned = write_plan(plan)
    original = subprocess.Popen

    def launch(argv, **kwargs):
        if defect == "input_drift" and argv[-1] == "validator":
            target = Path(plan["input"]["path"])
            target.write_bytes(target.read_bytes() + b"drift")
        if defect == "output_drift" and argv[-1] == "bookkeeping":
            target = Path(plan["output_dir"]) / "catalog.tsv.manifest.json"
            target.write_text(target.read_text() + " ")
        return original(argv, **kwargs)

    monkeypatch.setattr(subprocess, "Popen", launch)
    assert migration.run(path, pinned, started=time.monotonic()) == 3
    output = Path(plan["output_dir"])
    receipt = json.loads((output / "terminal.json").read_text())
    assert receipt["status"] == "failed" and receipt["partial_outputs"]


def test_orphaned_owned_child_is_reaped_without_touching_unrelated_process():
    migration._subreaper()
    unrelated = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(30)"], start_new_session=True
    )
    proc = subprocess.Popen(
        [
            sys.executable,
            "-c",
            "import os,signal,time; child=os.fork(); signal.signal(signal.SIGTERM,signal.SIG_IGN); print('ready',flush=True) if child==0 else None; time.sleep(30) if child==0 else os._exit(0)",
        ],
        stdout=subprocess.PIPE,
        start_new_session=True,
    )
    try:
        assert proc.stdout.readline() == b"ready\n"
        proc.wait(timeout=3)
        result = migration.cleanup(
            proc,
            dict(migration.CEILINGS, term_seconds=1, kill_seconds=1, receipt_seconds=1),
            time.monotonic() + 4,
        )
        assert result["group_gone"] and result["direct_child_reaped"]
        assert unrelated.poll() is None
    finally:
        unrelated.terminate()
        unrelated.wait(timeout=3)


def test_unsafe_cleanup_cannot_return_validated(plan, monkeypatch):
    path, pinned = write_plan(plan)
    original = migration.cleanup

    def observed_unsafe(*args):
        result = original(*args)
        return {**result, "group_gone": False}

    monkeypatch.setattr(migration, "cleanup", observed_unsafe)
    assert migration.run(path, pinned, started=time.monotonic()) == 3
    receipt = json.loads((Path(plan["output_dir"]) / "terminal.json").read_text())
    assert receipt["status"] == "unsafe_cleanup"


def test_reviewed_external_watchdog_literal_completes_tiny_fixture(plan):
    path, pinned = write_plan(plan)
    completed = subprocess.run(
        [
            "/usr/bin/timeout",
            "--signal=TERM",
            "--kill-after=5s",
            "40s",
            sys.executable,
            "-m",
            "chembounce_synthesis.catalog_migration",
            "run",
            "--plan",
            str(path),
            "--plan-sha256",
            pinned,
        ],
        cwd=plan["source"]["root"],
        env=dict(migration.ENVIRONMENT, PYTHONPATH=plan["source"]["root"]),
        capture_output=True,
        text=True,
        timeout=50,
    )
    assert completed.returncode == 0, completed.stderr
    output = Path(plan["output_dir"])
    assert json.loads((output / "terminal.json").read_text())["status"] == "validated"


def test_reportable_memory_error_has_explicit_incomplete_receipt(plan, monkeypatch, capsys):
    path, pinned = write_plan(plan)
    output = Path(plan["output_dir"])
    output.mkdir()

    def exhausted(*args):
        raise MemoryError("controlled fixture allocation failure")

    monkeypatch.setattr(migration, "worker", exhausted)
    assert (
        migration.main(["_worker", "--plan", str(path), "--plan-sha256", pinned, "--phase", "converter"]) == 2
    )
    failure = json.loads((output / "converter-failure.json").read_text())
    assert failure["error_type"] == "MemoryError" and failure["status"] == "incomplete"
    assert failure["complete_source"] is False
    assert "MemoryError" in capsys.readouterr().err
