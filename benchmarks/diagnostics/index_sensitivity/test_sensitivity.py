"""Synthetic mechanism checks, not scaffold-hopping performance evidence."""

import copy
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
spec = importlib.util.spec_from_file_location("index_sensitivity", HERE / "sensitivity.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return path


def stopped(pid):
    import psutil

    try:
        return psutil.Process(pid).status() == psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess:
        return True


@pytest.fixture
def index_report(tmp_path):
    index = tmp_path / "index"
    index.mkdir()
    for name in ("fingerprints.npy", "smiles.npy"):
        (index / name).write_bytes(b"synthetic")
    metadata = {
        "schema_version": 2,
        "source_sha256": "s" * 64,
        "canonical_sorted_unique": True,
        "rdkit_version": "2025.09.6",
        "rows": 4,
        "bitorder": "big",
        "fp_size": 2048,
        "include_chirality": False,
        "radius": 2,
        "canonicalization_policy": "rdkit-roundtrip-fixed-point-v1",
        "canonicalization_max_roundtrips": 8,
    }
    write(index / "metadata.json", metadata)
    pin = m.directory_pin(index)
    report = {
        "complete": True,
        "input_asset_sha256": pin["files_sha256"],
        "fingerprint_parameters": {
            k: metadata[k] for k in ("bitorder", "fp_size", "include_chirality", "radius")
        },
        "rows": 4,
        "fingerprints_checked": 4,
        "runtime_canonical_unique_strings": 4,
        "raw_strings_strictly_sorted_unique": True,
        "invalid_rows": [],
        "fingerprint_mismatches": [],
        "canonical_string_changed_rows": [],
        "canonical_duplicate_rows": [],
        "database_rdkit_version": "2025.09.6",
        "runtime_rdkit_version": "2025.09.6",
    }
    return pin, metadata, {"rdkit_runtime": "2025.09.6"}, report


def test_independent_report_is_required_and_binds_index(index_report):
    pin, metadata, rt, report = index_report
    assert m.validate_index_report(pin, metadata["source_sha256"], rt, report) == metadata


@pytest.mark.parametrize(
    "field,value",
    [
        ("complete", False),
        ("fingerprints_checked", 3),
        ("rows", 3),
        ("runtime_canonical_unique_strings", 3),
        ("raw_strings_strictly_sorted_unique", False),
        ("invalid_rows", [0]),
        ("fingerprint_mismatches", [0]),
        ("canonical_string_changed_rows", [1, 2, 3, 4, 5]),
        ("canonical_duplicate_rows", [0]),
        ("input_asset_sha256", {}),
        ("runtime_rdkit_version", "2026.09.1"),
    ],
)
def test_index_audit_failures_are_not_approval(index_report, field, value):
    pin, metadata, rt, report = index_report
    report[field] = value
    report["passed_fingerprint_compatibility"] = (
        True  # This alone proves neither fixed-point stability nor completeness.
    )
    with pytest.raises(ValueError):
        m.validate_index_report(pin, metadata["source_sha256"], rt, report)


def test_index_metadata_wrong_source_and_policy_rejected(index_report):
    pin, metadata, rt, report = index_report
    with pytest.raises(ValueError, match="source differs"):
        m.validate_index_report(pin, "wrong", rt, report)
    metadata["canonicalization_policy"] = "one-write"
    write(Path(pin["path"]) / "metadata.json", metadata)
    with pytest.raises(ValueError, match="fixed-point"):
        m.validate_index_report(pin, metadata["source_sha256"], rt, report)


def test_canonicalize_deduplicate_before_cap_and_preserve_order():
    assert m.canonical_selection(["OCC", "CCO", "CCN", "CCC"]) == {
        "canonical_deduplicated": ["CCO", "CCN", "CCC"],
        "downstream_selected": ["CCO", "CCN"],
    }
    assert m.canonical_selection(["CCN", "OCC", "CCC"])["downstream_selected"] == ["CCN", "CCO"]
    assert m.canonical_selection(["[CH3:2][OH:7]", "CO", "CC"]) == {
        "canonical_deduplicated": ["CO", "CC"],
        "downstream_selected": ["CO", "CC"],
    }


def table(path, smiles, newline="\n"):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator=newline)
        writer.writerow(m.TABLE_COLUMNS)
        for s in smiles:
            writer.writerow(
                [
                    "1",
                    "C1CCCCC1",
                    "C1CCCC1",
                    s,
                    s,
                    "0.5000000000000001",
                    "",
                    "0",
                    "0.5",
                    "2",
                    "1",
                    "100",
                    "0",
                    "1",
                ]
            )


def view(path):
    return {
        "table": m.table_view(path),
        "intermediates_sha256": {"fragment_info.tsv": "x"},
        "summary_without_database_path": {"status": "success"},
    }


def test_raw_bytes_full_table_and_selected_order_are_separate(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    table(a, ["CCO", "CCN", "CCC"], "\n")
    table(b, ["CCO", "CCN", "CCC"], "\r\n")
    check = m.compare_views(view(a), view(b))
    assert not check["raw_table_bytes_equal"] and check["full_table_cells_equal"]
    assert check["downstream_first_two_equal"]
    table(b, ["CCO", "CCN", "CCCC"])
    check = m.compare_views(view(a), view(b))
    assert not check["full_table_cells_equal"] and check["downstream_first_two_equal"]
    table(b, ["CCN", "CCO", "CCC"])
    assert not m.compare_views(view(a), view(b))["downstream_first_two_equal"]


def test_invalid_downstream_selection_is_unassessed_not_equal(tmp_path):
    a = tmp_path / "a"
    table(a, ["CC.CC"])
    result = view(a)
    assert result["table"]["downstream_selection_status"] == "error"
    assert m.compare_views(result, result)["downstream_first_two_equal"] is None


def test_four_generator_results_and_bad_headers_rejected(tmp_path):
    path = tmp_path / "a"
    table(path, ["C", "CC", "CCC", "CCCC"])
    with pytest.raises(ValueError, match="three-output"):
        m.table_view(path)
    path.write_text("Final structure\nCC\n")
    with pytest.raises(ValueError, match="schema"):
        m.table_view(path)


def test_failure_and_unattempted_denominators_and_replay_mismatch(tmp_path):
    path = tmp_path / "a"
    table(path, ["CCO", "CCN"])
    v = view(path)
    jobs = [{"id": str(i), "options": {"search_policy": "diverse"}, "seed": 42} for i in range(3)]
    data = {
        "jobs": jobs,
        "references": [
            {
                "job_id": "0",
                "status": "complete_generator_reference",
                "original_campaign_status": "timeout",
                "view": v,
            }
        ],
    }
    reports = [{"job_id": "0", "arm": a, "status": "success", "view": v} for a in ("old", "new")]
    reports += [{"job_id": "1", "arm": "old", "status": "timeout"}]
    result = m.summarize(data, reports)
    assert result["planned_arm_calls"] == 6 and result["unattempted_arm_calls"] == 3
    assert result["condition_statuses"] == {"assessed": 1, "unassessed": 2}
    assert result["conditions"][0]["old_replay_status"] == "matched"
    assert data["references"][0]["original_campaign_status"] == "timeout"
    changed = copy.deepcopy(v)
    changed["table"]["raw_designs"].reverse()
    changed["table"]["downstream_selected"].reverse()
    reports[0]["view"] = changed
    assert (
        m.summarize(data, reports)["conditions"][0]["asset_only_attribution_status"]
        == "replay_mismatch_blocks_attribution"
    )
    with pytest.raises(ValueError, match="Duplicate"):
        m.summarize(data, reports + [reports[0]])


def test_snapshot_is_byte_preserved_and_never_overwrites(tmp_path):
    source, dest = tmp_path / "source", tmp_path / "dest"
    source.write_bytes(b"immutable\r\n")
    pin = m.stable_copy(source, dest)
    assert dest.read_bytes() == source.read_bytes() and m.sha(dest) == pin["sha256"]
    source.write_text("later source changed")
    assert dest.read_bytes() == b"immutable\r\n"
    with pytest.raises(FileExistsError):
        m.stable_copy(source, dest)
    with pytest.raises(ValueError, match="modify"):
        m.protect_output(tmp_path / "campaign" / "new", [tmp_path / "campaign"])


def test_complete_generator_reference_does_not_promote_pipeline_failure(tmp_path, monkeypatch):
    run = tmp_path / "archive" / "run"
    design = run / "designs"
    design.mkdir(parents=True)
    write(design / "run_manifest.json", {"fixture": True})
    write(design / "run_summary.json", {"status": "success"})
    table(design / "overall_result.txt", ["CCO"])
    request = {
        "query": {"id": "q", "smiles": "CC"},
        "seed": 42,
        "kind": "pipeline",
        "method_id": "pipeline",
        "settings": {},
        "budget": {"max_candidates": 3, "wall_seconds": 360},
        "assets": {"fingerprint_db": "/fixture/index", "scaffold_db": "/fixture/scaffolds"},
    }
    write(run / "request.json", request)
    write(run / "evaluation.json", {"status": "timeout"})
    member = {
        "run_directory": str(run),
        "query": request["query"],
        "method_id": "pipeline",
        "settings": {},
        "old_index": {"path": "/fixture/index"},
        "scaffold": {"path": "/fixture/scaffolds"},
    }
    job = {"id": "job", "query_smiles": "CC", "seed": 42, "members": [member]}
    monkeypatch.setattr(m, "inspect_generator", lambda *a: view(design / "overall_result.txt"))
    output = tmp_path / "snapshot"
    refs = m.snapshot_references([job], output, {}, {})
    assert refs[0]["status"] == "complete_generator_reference"
    assert refs[0]["original_campaign_status"] == "timeout"
    frozen = output / refs[0]["snapshot"] / "designs/overall_result.txt"
    assert frozen.read_bytes() == (design / "overall_result.txt").read_bytes()
    (design / "overall_result.txt").write_text("later original mutation")
    assert m.sha(frozen) == refs[0]["files"]["designs/overall_result.txt"]["sha256"]
    (design / "overall_result.txt").unlink()
    assert m.snapshot_references([job], tmp_path / "later", {}, {})[0]["status"] == "incomplete"


@pytest.mark.parametrize("unreaped", [False, True])
def test_bounded_parallel_orchestration_retains_all_252_calls(tmp_path, monkeypatch, unreaped):
    import threading

    path = tmp_path / "table"
    table(path, ["CCO", "CCN"])
    v = view(path)
    jobs = [
        {
            "id": f"{i:064x}",
            "query_smiles": "CC",
            "seed": 42,
            "options": {"search_policy": "diverse"},
            "members": [],
        }
        for i in range(126)
    ]
    plan_path = write(tmp_path / "plan" / "plan.json", {"synthetic": True})
    data = {
        "jobs": jobs,
        "references": [],
        "source_root": str(ROOT),
        "source_sha256": {},
        "runtime": m.runtime(),
        "environment": m.environment_record(m.execution_environment(ROOT)),
        "python_executable": sys.executable,
        "indices": {a: {"path": str(tmp_path / a)} for a in ("old", "new")},
        "scaffold": {},
        "inputs": [],
        "execution": {
            "parallelism": 1 if unreaped else 2,
            "wall_seconds": 360,
            "max_memory_mb": 24000,
            "memory_sample_seconds": 0.05,
        },
    }
    monkeypatch.setattr(m, "verify_plan", lambda *a: data)
    monkeypatch.setattr(m, "inspect_generator", lambda *a: v)
    lock = threading.Lock()
    observed = {"active": 0, "peak": 0, "calls": 0}

    def process(command, cwd, env, directory, **kwargs):
        with lock:
            observed["active"] += 1
            observed["peak"] = max(observed["peak"], observed["active"])
            observed["calls"] += 1
        try:
            time.sleep(0.005)
            if unreaped:
                return {
                    "status": "cleanup_unreaped",
                    "cleanup": {"reaped": False},
                    "resources": {"wall_seconds": 361.5},
                }
            ident = int(m.read_json(directory / "request.json")["job"]["id"], 16)
            if ident == 9:
                return {"status": "timeout", "resources": {"wall_seconds": 360}}
            if ident == 10:
                return {"status": "completed", "resources": {"wall_seconds": 1}}
            if ident == 11:
                (directory / "worker-result.json").write_text("invalid JSON")
            else:
                write(
                    directory / "worker-result.json",
                    {"status": "success", "view": v, "generator_cost": {"fixture": True}},
                )
            return {"status": "completed", "resources": {"wall_seconds": 1}}
        finally:
            with lock:
                observed["active"] -= 1

    monkeypatch.setattr(m, "run_process", process)
    args = type(
        "Args", (), {"output": str(tmp_path / "run"), "plan": str(plan_path), "plan_sha256": m.sha(plan_path)}
    )()
    if unreaped:
        with pytest.raises(RuntimeError, match="could not be reaped"):
            m.run(args)
        assert observed["calls"] == 1
        summary = m.read_json(Path(args.output) / "summary.json")
        assert summary["unattempted_arm_calls"] == 251
        assert summary["arm_statuses"] == {"cleanup_unreaped": 1}
        assert m.read_json(Path(args.output) / "progress.json")["state"] == "interrupted_requires_review"
        return
    result = m.run(args)
    assert result["recorded_arm_calls"] == 252
    assert observed == {"active": 0, "peak": 2, "calls": 252}
    reports = m.read_json(Path(args.output) / "runs.json")
    assert len({(r["job_id"], r["arm"]) for r in reports}) == 252
    summary = m.read_json(Path(args.output) / "summary.json")
    assert summary["planned_conditions"] == 126 and summary["unattempted_arm_calls"] == 0
    assert summary["arm_statuses"] == {
        "success": 246,
        "timeout": 2,
        "missing_worker_result": 2,
        "infrastructure_error": 2,
    }
    assert summary["condition_statuses"] == {"assessed": 123, "unassessed": 3}
    assert m.read_json(Path(args.output) / "progress.json")["state"] == "complete"
    with pytest.raises(ValueError, match="already exists"):
        m.run(args)


def test_process_error_exit_kept(tmp_path):
    result = m.run_process(
        [sys.executable, "-c", "raise SystemExit(7)"], tmp_path, os.environ.copy(), tmp_path, 2, 256
    )
    assert result["status"] == "worker_error" and result["returncode"] == 7
    assert result["resources"]["wall_seconds"] > 0


def test_bounded_reap_reports_uninterruptible_worker(monkeypatch):
    waits, signals = [], []

    class Process:
        pid = 123456789

        def wait(self, timeout):
            waits.append(timeout)
            raise subprocess.TimeoutExpired("uninterruptible-fixture", timeout)

    monkeypatch.setattr(m.os, "killpg", lambda pid, signum: signals.append((pid, signum)))
    result = m.stop_group(Process())
    assert waits == [0.5, 1.0] and len(signals) == 2
    assert result["reaped"] is False and result["status"] == "unreaped_requires_parent_review"


def test_atomic_snapshot_survives_interrupted_write(tmp_path, monkeypatch):
    path = write(tmp_path / "record.json", {"previous": True})

    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt("fixture interruption")

    monkeypatch.setattr(m.json, "dump", interrupted)
    with pytest.raises(KeyboardInterrupt):
        m.atomic_json(path, {"new": True})
    assert m.read_json(path) == {"previous": True}
    assert not path.with_name("record.json.tmp").exists()


def test_real_controller_sigterm_cancels_workers_and_preserves_denominators(tmp_path):
    """Exercise actual OS signals, controller threads and separate worker sessions."""
    import signal

    controller = tmp_path / "controller.py"
    controller.write_text(
        "import importlib.util,json,sys,os,time\nfrom pathlib import Path\nfrom types import SimpleNamespace\n"
        f"spec=importlib.util.spec_from_file_location('helper',{str(HERE / 'sensitivity.py')!r})\n"
        "m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)\n"
        f"base=Path({str(tmp_path)!r});root={str(ROOT)!r}\n"
        "plan=base/'plan'/'plan.json';plan.parent.mkdir();plan.write_text('{}')\n"
        "jobs=[{'id':f'{i:064x}','query_smiles':'CC','seed':42,'options':{'search_policy':'diverse'},'members':[]} for i in range(126)]\n"
        "data={'jobs':jobs,'references':[],'source_root':root,'source_sha256':{},'runtime':{},'environment':{},'python_executable':sys.executable,'indices':{a:{'path':str(base/a)} for a in ('old','new')},'scaffold':{},'inputs':[],'execution':{'parallelism':2,'wall_seconds':360,'max_memory_mb':256,'memory_sample_seconds':0.01}}\n"
        "m.verify_plan=lambda *a:data\nreal_process=m.run_process\n"
        "def process(command,cwd,env,directory,**kwargs):\n"
        " code='import os,signal,time;signal.signal(signal.SIGTERM,signal.SIG_IGN);open('+repr(str(directory/'child.pid'))+\",'w').write(str(os.getpid()));time.sleep(60)\"\n"
        " return real_process([sys.executable,'-c',code],cwd,env,directory,**kwargs)\n"
        "m.run_process=process\n"
        "m.run(SimpleNamespace(output=str(base/'run'),plan=str(plan),plan_sha256=m.sha(plan)))\n"
    )
    with (tmp_path / "controller.log").open("wb") as log:
        proc = subprocess.Popen(
            [sys.executable, str(controller)],
            cwd=ROOT,
            env=m.execution_environment(ROOT),
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            deadline = time.monotonic() + 10
            while time.monotonic() < deadline:
                pidfiles = list((tmp_path / "run").glob("*/*/child.pid"))
                if len(pidfiles) == 2 and all(p.read_text() for p in pidfiles):
                    break
                if proc.poll() is not None:
                    pytest.fail((tmp_path / "controller.log").read_text())
                time.sleep(0.02)
            else:
                pytest.fail("Fixture controller did not start two workers")
            pids = [int(p.read_text()) for p in pidfiles]
            os.kill(proc.pid, signal.SIGTERM)
            assert proc.wait(timeout=5) != 0
            progress = m.read_json(tmp_path / "run/progress.json")
            assert progress["state"] == "interrupted_requires_review" and progress["completed_arm_calls"] == 2
            summary = m.read_json(tmp_path / "run/summary.json")
            assert summary["planned_arm_calls"] == 252 and summary["unattempted_arm_calls"] == 250
            assert summary["arm_statuses"] == {"infrastructure_interruption": 2}
            assert all(stopped(pid) for pid in pids)
        finally:
            if proc.poll() is None:
                m.stop_group(proc)
            for pidfile in (tmp_path / "run").glob("*/*/child.pid"):
                if pidfile.read_text():
                    pid = int(pidfile.read_text())
                    try:
                        os.kill(pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass


def test_timeout_kills_sigterm_ignoring_process_group(tmp_path):
    pidfile = tmp_path / "child.pid"
    child_code = (
        "import os,signal,time;signal.signal(signal.SIGTERM,signal.SIG_IGN);open("
        + repr(str(pidfile))
        + ",'w').write(str(os.getpid()));time.sleep(60)"
    )
    parent_code = (
        "import subprocess,sys,time,signal;signal.signal(signal.SIGTERM,signal.SIG_IGN);subprocess.Popen([sys.executable,'-c',"
        + repr(child_code)
        + "]);time.sleep(60)"
    )
    result = m.run_process(
        [sys.executable, "-c", parent_code], tmp_path, os.environ.copy(), tmp_path, 0.4, 256, 0.01
    )
    assert result["status"] == "timeout"
    assert result["resources"]["wall_seconds"] < 2
    assert pidfile.exists()
    pid = int(pidfile.read_text())
    for _ in range(100):
        if stopped(pid):
            break
        time.sleep(0.01)
    else:
        pytest.fail("Grandchild survived process-group timeout")


def test_memory_limit_recorded(tmp_path):
    result = m.run_process(
        [sys.executable, "-c", "import time;x=bytearray(32*1024**2);time.sleep(3)"],
        tmp_path,
        os.environ.copy(),
        tmp_path,
        2,
        8,
        0.01,
    )
    assert result["status"] == "memory_limit"


@pytest.fixture
def synthetic_campaigns(tmp_path):
    assets = tmp_path / "assets"
    assets.mkdir()
    scaffold = assets / "scaffolds.smi"
    scaffold.write_text("C1CCCCC1\n")
    index = assets / "index"
    index.mkdir()
    for name in ("metadata.json", "fingerprints.npy", "smiles.npy"):
        (index / name).write_text("synthetic fixture")
    config = write(assets / "config.json", {"max_designs": 2})
    asset_specs = {
        "scaffold_db": {"artifact": str(scaffold), "sha256": m.sha(scaffold)},
        "fingerprint_db": {"directory": str(index), "files_sha256": m.directory_pin(index)["files_sha256"]},
        "synthesis_config": {"artifact": str(config), "sha256": m.sha(config)},
    }
    campaigns = []
    for split, n, offset in (("main", 12, 0), ("time", 12, 12), ("target", 6, 24)):
        dataset = write(
            tmp_path / f"{split}-panel.json",
            {
                "queries": [{"id": f"q{offset + i}", "compound_id": f"c{offset + i}"} for i in range(n)],
                "compounds": [
                    {
                        "id": f"c{offset + i}",
                        "smiles": "C" * (offset + i + 1),
                        "target_id": "TARGET",
                        "label": "private",
                    }
                    for i in range(n)
                ],
            },
        )
        for gpu in (0, 1):
            if split == "main":
                methods = (
                    ["pipeline", "pipeline_no_3d", "pipeline_no_refinement"]
                    if gpu == 0
                    else [
                        "synllama_only",
                        "pipeline_legacy_search",
                        "pipeline_encounter_allocation",
                        "pipeline_rule_only",
                    ]
                )
            else:
                methods = ["pipeline"] if gpu == 0 else ["synllama_only"]
            definitions = []
            for name in methods:
                definitions.append(
                    {
                        "id": name,
                        "kind": "synllama_only" if name == "synllama_only" else "pipeline",
                        "assets": asset_specs,
                        "settings": {
                            "generator_options": {
                                **m.OPTIONS,
                                "search_policy": "legacy" if "legacy" in name else "diverse",
                            }
                        },
                    }
                )
            campaign = (
                {"thread_environment": m.THREADS, "parent_protocol_sha256": "1" * 64}
                if split == "main"
                else {"primary_protocol_sha256": "1" * 64}
            )
            protocol = write(
                tmp_path / f"{split}-gpu{gpu}.json",
                {
                    "seeds": [42, 43, 44],
                    "budget": {"wall_seconds": 360, "max_memory_mb": 24000, "max_candidates": 3},
                    "campaign": campaign,
                    "methods": definitions,
                },
            )
            campaigns.append([str(dataset), str(protocol), str(tmp_path / f"archive-{split}-{gpu}")])
    return campaigns


def test_full_grid_dedup_and_private_label_exclusion(synthetic_campaigns):
    jobs, inputs = m.build_grid(synthetic_campaigns)
    assert len(jobs) == 126 and sum(len(j["members"]) for j in jobs) == 270
    assert len(inputs) == 6
    assert sum(i["declared_thread_environment"] is None for i in inputs) == 4
    assert "private" not in json.dumps(jobs)
    assert {j["seed"] for j in jobs} == {42, 43, 44}
    assert sum(j["options"]["search_policy"] == "legacy" for j in jobs) == 36


@pytest.mark.parametrize(
    "change", ["seed", "options", "cap", "missing_campaign", "asset", "parent_link", "threads"]
)
def test_frozen_grid_deviation_rejected(synthetic_campaigns, change):
    campaigns = copy.deepcopy(synthetic_campaigns)
    if change == "missing_campaign":
        campaigns.pop(0)
    else:
        path = Path(campaigns[0][1])
        protocol = m.read_json(path)
        if change == "seed":
            protocol["seeds"] = [42]
        elif change == "options":
            protocol["methods"][0]["settings"]["generator_options"]["calculate_shape"] = True
        elif change == "cap":
            protocol["budget"]["max_candidates"] = 2
        elif change == "asset":
            protocol["methods"][0]["assets"]["scaffold_db"]["sha256"] = "wrong"
        elif change == "parent_link":
            for entry in campaigns[2:]:
                secondary = m.read_json(entry[1])
                secondary["campaign"]["primary_protocol_sha256"] = "2" * 64
                write(Path(entry[1]), secondary)
        elif change == "threads":
            protocol["campaign"]["thread_environment"]["OMP_NUM_THREADS"] = "2"
        write(path, protocol)
    with pytest.raises(ValueError):
        m.build_grid(campaigns)


def test_source_mutation_and_unlisted_source_fail_closed(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    module = root / "chembounce.py"
    module.write_text("# frozen\n")
    expected = {"chembounce.py": m.sha(module)}
    m.source_check(root, expected)
    (root / "unexpected.py").write_text("# extra")
    with pytest.raises(ValueError, match="inventory changed"):
        m.source_check(root, expected)
    (root / "unexpected.py").unlink()
    module.write_text("# changed")
    with pytest.raises(ValueError, match="source mismatch"):
        m.source_check(root, expected)


def test_plan_existing_output_rejected_before_any_sources(tmp_path):
    class Args:
        output = str(tmp_path)

    with pytest.raises(ValueError, match="already exists"):
        m.plan(Args())
    with pytest.raises(ValueError, match="already exists"):
        m.run(Args())


def test_real_frozen_generator_small_fixture_smoke(tmp_path):
    """One tiny synthetic query/index; no public campaign chemistry or model calls."""
    frozen_setting = os.environ.get("CHEMBOUNCE_FROZEN_SOURCE_ROOT")
    if not frozen_setting:
        pytest.skip(
            "Set CHEMBOUNCE_FROZEN_SOURCE_ROOT to the separately prepared original frozen "
            "checkout with its frozen-campaign-code-provenance.json; current develop is not a substitute"
        )
    frozen_root = Path(frozen_setting).resolve()
    frozen_manifest = (
        frozen_root / "benchmark-results/research-validation/frozen-campaign-code-provenance.json"
    )
    if not frozen_manifest.is_file():
        pytest.fail(f"Frozen smoke fixture preparation missing: {frozen_manifest}")
    import numpy as np
    from rdkit import Chem, rdBase
    from rdkit.Chem import rdFingerprintGenerator

    frozen = m.read_json(frozen_manifest)
    assert frozen["implementation_commit"] == "d972124d587e647ef48455cf1a038e6be96e5b04"
    m.source_check(frozen_root, frozen["code_sha256"])
    index = tmp_path / "index"
    index.mkdir()
    scaffold = tmp_path / "scaffolds.smi"
    strings = sorted(["C1CCCCC1", "C1CCCC1", "c1ccccc1", "c1ccncc1"])
    scaffold.write_text("\n".join(strings) + "\n")
    fp = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    np.save(index / "smiles.npy", np.array(strings))
    np.save(
        index / "fingerprints.npy",
        np.array([np.packbits(fp.GetFingerprintAsNumPy(Chem.MolFromSmiles(s))) for s in strings]),
    )
    write(
        index / "metadata.json",
        {
            "schema_version": 2,
            "packed": True,
            "bitorder": "big",
            "fp_size": 2048,
            "radius": 2,
            "include_chirality": False,
            "canonical_sorted_unique": True,
            "rdkit_version": rdBase.rdkitVersion,
            "source_sha256": m.sha(scaffold),
            "rows": 4,
        },
    )
    output = tmp_path / "smoke"
    output.mkdir()
    job = {
        "id": "synthetic-only",
        "query_smiles": "CCOc1ccccc1",
        "seed": 42,
        "options": {**m.OPTIONS, "search_policy": "diverse"},
        "max_designs": 2,
        "scaffold_sha256": m.sha(scaffold),
    }
    env = m.execution_environment(frozen_root)
    request = {
        "output": str(output),
        "source_root": str(frozen_root),
        "source_sha256": frozen["code_sha256"],
        "runtime": m.runtime(),
        "environment": m.environment_record(env),
        "job": job,
        "scaffold": m.file_pin(scaffold),
        "index": m.directory_pin(index),
    }
    path = write(output / "request.json", request)
    command = [sys.executable, str(HERE / "sensitivity.py"), "_worker", "--request", str(path)]
    result = m.run_process(command, frozen_root, env, output, 30, 2048)
    record = m.read_json(output / "worker-result.json")
    assert result["status"] == "completed", (record, (output / "worker.log").read_text())
    assert record["status"] == "success" and record["no_prohibited_modules_imported"]
    assert len(record["view"]["table"]["raw_designs"]) <= 3
    assert len(record["view"]["table"]["downstream_selected"]) <= 2
    request["source_sha256"]["chembounce.py"] = "wrong"
    other = tmp_path / "bad"
    other.mkdir()
    request["output"] = str(other)
    path = write(other / "request.json", request)
    bad = m.run_process(
        [sys.executable, str(HERE / "sensitivity.py"), "_worker", "--request", str(path)],
        frozen_root,
        env,
        other,
        10,
        2048,
    )
    assert bad["status"] == "worker_error"
    assert "Frozen source mismatch" in m.read_json(other / "worker-result.json")["error"]
