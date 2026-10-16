"""Synthetic campaigns only: no panel, chemical runtime, model, GPU or live runner."""

import copy
import fcntl
import multiprocessing
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest

from benchmarks.diagnostics import supervise_comparator_followup as supervisor


FAKE_RUNNER = """
import argparse,hashlib,json,os,pathlib,subprocess,sys,time
p=argparse.ArgumentParser()
for k in ('dataset','protocol','output'): p.add_argument('--'+k,required=True)
a=p.parse_args(); root=pathlib.Path(a.output); root.mkdir(parents=True,exist_ok=False)
load=lambda p:json.loads(pathlib.Path(p).read_text())
def save(p,v):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v))
sha=lambda p:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
protocol=load(a.protocol); panel=load(a.dataset);mode=protocol['fixture_mode']
if mode=='sleep':
 child=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)'])
 save(root/'ready.json',{'pid':os.getpid(),'grandchild':child.pid,'pgid':os.getpgrp()})
 time.sleep(30)
rows=[]
for index,q in enumerate(panel):
 for seed in (42,43,44):
  status=('timeout' if mode=='failed' else 'budget_exceeded') if mode in ('failed','budget') and not rows else 'no_candidates'
  directory=root/protocol['method_id']/f'query-{index:04d}-seed-{seed}'
  response={'candidates':[{'smiles':'rejected-fixture'}]} if status=='budget_exceeded' else {'candidates':[]} if status=='no_candidates' else {}
  if response:save(directory/'response.json',response)
  row={'method_id':protocol['method_id'],'query_id':q['id'],'target_id':q['target_id'],
       'seed':seed,'status':status,'response':response,'evaluation':{'candidates':[]},
       'response_sha256':sha(directory/'response.json') if response else None}
  save(directory/'evaluation.json',row);rows.append(row)
if mode=='incomplete': rows=rows[:-1]
save(root/'runs.json',rows)
save(root/'manifest.json',{'protocol_sha256':sha(a.protocol),
 'dataset':{'dataset_sha256':sha(a.dataset)},'code_sha256':protocol['source_sha256']})
save(root/'campaign.json',{'status':'completed','termination':None,'attempted_runs':len(rows)})
print('fixture completed',len(rows),flush=True)
sys.exit(7 if mode=='unexpected' else 1 if mode in ('failed','budget') else 0)
"""


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(supervisor.encode(value))


@pytest.fixture
def study(tmp_path, monkeypatch):
    prep = supervisor.preparation
    source = tmp_path / "source"
    module = source / "benchmarks/research/runner.py"
    module.parent.mkdir(parents=True)
    module.write_text(FAKE_RUNNER)
    (source / "benchmarks/__init__.py").write_text("")
    (source / "benchmarks/research/__init__.py").write_text("")
    source_info = {"root": str(source), "commit": "synthetic", "files_sha256": {"fixture": "pinned"}}
    runtime_info = {"rdkit_policy": "synthetic", "metadata_probe_crem_root": None}
    monkeypatch.setattr(prep, "source_snapshot", lambda *a: source_info)
    monkeypatch.setattr(prep, "runtime_snapshot", lambda *a: runtime_info)
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    commands, grid, protocols, inputs = [], [], {}, {str(module): prep.digest(module)}
    for split in prep.SPLITS:
        queries = [
            {"id": f"{split}-q{i}", "target_id": f"t{i // 2}", "smiles": "fixture"}
            for i in range(prep.COUNTS[split])
        ]
        dataset = tmp_path / f"{split}-panel.json"
        save(dataset, queries)
        inputs[str(dataset)] = prep.digest(dataset)
        save(prepared / f"public-queries/{split}.json", queries)
        for method in prep.METHODS:
            relative = f"protocols/{split}-{method}.json"
            save(
                prepared / relative,
                {
                    "method_id": method,
                    "fixture_mode": "success",
                    "source_sha256": source_info["files_sha256"],
                },
            )
            protocols[f"{split}/{method}"] = relative
            output = tmp_path / "prepared-runs" / split / method
            commands.append(
                {
                    "split": split,
                    "method_id": method,
                    "expected_attempts": len(queries) * 3,
                    "execution_order": len(commands),
                    "no_retry_or_resume": True,
                    "cwd": str(source),
                    "argv": [
                        sys.executable,
                        "-m",
                        "benchmarks.research.runner",
                        "--dataset",
                        str(dataset),
                        "--protocol",
                        str(prepared / relative),
                        "--output",
                        str(output),
                    ],
                    "environment_overrides": {
                        **prep.THREAD_ENV,
                        "PYTHONPATH": str(source),
                        "CUDA_VISIBLE_DEVICES": "" if method == prep.METHODS[1] else supervisor.GPU_UUID,
                    },
                }
            )
            grid.extend(
                {
                    "split": split,
                    "method_id": method,
                    "query_id": q["id"],
                    "target_id": q["target_id"],
                    "query_index": i,
                    "seed": seed,
                    "output_dir": str(output / method / f"query-{i:04d}-seed-{seed}"),
                }
                for i, q in enumerate(queries)
                for seed in prep.SEEDS
            )
    plan = {
        "expected_attempts": 270,
        "commands": commands,
        "protocols": protocols,
        "source": source_info,
        "runtimes": {sys.executable: runtime_info},
        "input_files_sha256": inputs,
        "input_directory_membership": {},
        "input_python_membership": {},
        "upstream_commits": {},
    }
    save(prepared / "execution-plan.json", plan)
    save(prepared / "grid.json", grid)
    pin(prepared)
    return prepared, tmp_path / "runtime"


def pin(prepared):
    save(
        prepared / "manifest.json",
        {
            "files_sha256": {
                str(p.relative_to(prepared)): supervisor.digest(p)
                for p in prepared.rglob("*")
                if p.is_file() and p.name != "manifest.json"
            }
        },
    )
    return supervisor.digest(prepared / "manifest.json")


def set_mode(prepared, mode):
    command = supervisor.read_json(prepared / "execution-plan.json")["commands"][0]
    path = Path(command["argv"][6])
    value = supervisor.read_json(path)
    value["fixture_mode"] = mode
    save(path, value)
    return pin(prepared)


@pytest.mark.parametrize("mode,failed_status", [("failed", "timeout"), ("budget", "budget_exceeded")])
def test_exact_nine_sequential_campaigns_and_retained_failure(study, mode, failed_status):
    prepared, runtime = study
    sha = set_mode(prepared, mode)
    assert supervisor.launch(prepared / "execution-plan.json", sha, runtime) == 0
    progress = supervisor.read_json(runtime / "progress.json")
    assert (progress["status"], progress["completed_campaigns"], progress["completed_attempts"]) == (
        "completed",
        9,
        270,
    )
    campaigns = progress["campaigns"]
    assert campaigns[0]["exit_code"] == 1 and campaigns[0]["status_counts"][failed_status] == 1
    if mode == "budget":
        rows = supervisor.read_json(Path(campaigns[0]["command"]["argv"][-1]) / "runs.json")
        assert rows[0]["response"]["candidates"] and rows[0]["evaluation"]["candidates"] == []
    assert all(c["child"]["pid"] == c["child"]["process_group_id"] for c in campaigns)
    assert all(a["finished_utc"] <= b["started_utc"] for a, b in zip(campaigns, campaigns[1:]))
    plan = supervisor.read_json(prepared / "execution-plan.json")
    assert [c["command"] for c in campaigns] == plan["commands"]
    for c in campaigns:
        assert c["child"]["create_time"] > 0 and c["child"]["start_ticks"] > 0
        for name, expected in c["logs_sha256"].items():
            assert supervisor.digest(runtime / name) == expected
    with pytest.raises(FileExistsError):
        supervisor.launch(prepared / "execution-plan.json", sha, runtime)


@pytest.mark.parametrize("mode", ["incomplete", "unexpected"])
def test_incomplete_or_unexpected_exit_stops_without_retry(study, mode):
    prepared, runtime = study
    sha = set_mode(prepared, mode)
    assert supervisor.launch(prepared / "execution-plan.json", sha, runtime) == 1
    progress = supervisor.read_json(runtime / "progress.json")
    assert progress["status"] == "failed" and progress["completed_attempts"] == 0
    assert len(progress["campaigns"]) == 1
    assert not (prepared.with_name("prepared-runs") / "scaffold" / supervisor.preparation.METHODS[1]).exists()


@pytest.mark.parametrize("change", ["grid", "gpu", "order"])
def test_grid_and_exact_command_refusal(study, change):
    prepared, runtime = study
    plan = supervisor.read_json(prepared / "execution-plan.json")
    if change == "grid":
        grid = supervisor.read_json(prepared / "grid.json")
        grid[-1] = copy.deepcopy(grid[0])
        save(prepared / "grid.json", grid)
    elif change == "gpu":
        plan["commands"][0]["environment_overrides"]["CUDA_VISIBLE_DEVICES"] = "2"
    else:
        plan["commands"][0], plan["commands"][1] = plan["commands"][1], plan["commands"][0]
    save(prepared / "execution-plan.json", plan)
    sha = pin(prepared)
    with pytest.raises(ValueError):
        supervisor.launch(prepared / "execution-plan.json", sha, runtime)
    assert not runtime.exists()


def test_external_manifest_pin_and_runtime_separation(study):
    prepared, runtime = study
    sha = pin(prepared)
    with pytest.raises(ValueError):
        supervisor.launch(prepared / "execution-plan.json", "0" * 64, runtime)
    for invalid in [prepared / "runtime", prepared.with_name("prepared-runs") / "runtime"]:
        with pytest.raises(ValueError):
            supervisor.launch(prepared / "execution-plan.json", sha, invalid)
    runtime.symlink_to(runtime.with_name("absent"))
    with pytest.raises(FileExistsError):
        supervisor.launch(prepared / "execution-plan.json", sha, runtime)


def test_shared_readonly_manifest_lock_blocks_other_runtime(study):
    prepared, runtime = study
    sha = pin(prepared)
    with (prepared / "manifest.json").open("rb") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(BlockingIOError):
            supervisor.launch(prepared / "execution-plan.json", sha, runtime)
    assert not runtime.exists() and supervisor.digest(prepared / "manifest.json") == sha


def test_detached_launch_passes_lock_and_new_session(study, monkeypatch):
    prepared, runtime = study
    sha = pin(prepared)
    calls = []

    def popen(argv, **kwargs):
        calls.append((argv, kwargs))
        return type("Child", (), {"pid": os.getpid()})()

    monkeypatch.setattr(supervisor.subprocess, "Popen", popen)
    assert supervisor.launch(prepared / "execution-plan.json", sha, runtime, detached=True) == 0
    argv, kwargs = calls[0]
    assert kwargs["start_new_session"] is True and len(kwargs["pass_fds"]) == 1
    assert argv[-2:] == ["--lock-fd", str(kwargs["pass_fds"][0])]
    assert supervisor.read_json(runtime / "detached.json")["argv"] == argv
    assert not (runtime / "started.json").exists()


def test_signal_only_verified_group_and_interruption_is_not_completion(study):
    prepared, runtime = study
    sha = set_mode(prepared, "sleep")
    innocent = subprocess.Popen([sys.executable, "-c", "import time;time.sleep(30)"], start_new_session=True)
    process = multiprocessing.get_context("fork").Process(
        target=supervisor.launch, args=(prepared / "execution-plan.json", sha, runtime)
    )
    process.start()
    try:
        ready = (
            prepared.with_name("prepared-runs")
            / "scaffold"
            / supervisor.preparation.METHODS[0]
            / "ready.json"
        )
        deadline = time.monotonic() + 10
        while not ready.exists() and time.monotonic() < deadline:
            time.sleep(0.02)
        assert ready.exists()
        child = supervisor.read_json(ready)
        os.kill(process.pid, signal.SIGTERM)
        process.join(timeout=10)
        assert not process.is_alive()
        state = supervisor.read_json(runtime / "progress.json")
        assert state["status"] == "interrupted" and state["completed_attempts"] == 0
        assert state["signal"] == signal.SIGTERM and len(state["campaigns"]) == 1
        assert state["campaigns"][0]["child"]["pid"] == child["pid"]
        assert innocent.poll() is None
    finally:
        if process.is_alive():
            process.terminate()
            process.join(5)
        innocent.terminate()
        innocent.wait(timeout=5)


def test_refuse_signal_if_birth_identity_changed(monkeypatch):
    child = type("Child", (), {"pid": 123, "poll": lambda self: None})()
    monkeypatch.setattr(
        supervisor, "birth", lambda pid: {"pid": pid, "create_time": 2, "process_group_id": pid}
    )
    sent = []
    monkeypatch.setattr(supervisor.os, "killpg", lambda *args: sent.append(args))
    with pytest.raises(RuntimeError):
        supervisor.signal_owned(
            child, {"pid": 123, "create_time": 1, "process_group_id": 123}, signal.SIGTERM
        )
    assert sent == []


def test_real_detached_worker_survives_launcher_and_uses_real_verifier(study, monkeypatch):
    """Only synthetic runner code executes; preparation validation is not mocked."""
    prepared, runtime = study
    plan = supervisor.read_json(prepared / "execution-plan.json")
    source = Path(plan["source"]["root"])
    # The reviewed source validator requires exactly 41 committed Python files.
    for index in range(39):
        (source / f"fixture_{index:02d}.py").write_text("# synthetic source pin\n")
    subprocess.run(["git", "init", "-q", str(source)], check=True)
    subprocess.run(["git", "-C", str(source), "add", "."], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(source),
            "-c",
            "user.name=Fixture",
            "-c",
            "user.email=fixture@example.invalid",
            "commit",
            "-qm",
            "synthetic fixture",
        ],
        check=True,
    )
    commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    monkeypatch.undo()
    prep = supervisor.preparation
    pins = prep.Pins()
    plan["source"] = prep.source_snapshot(source, commit, pins)
    plan["runtimes"] = {sys.executable: prep.runtime_snapshot(sys.executable, pins)}
    plan["input_files_sha256"].update(pins.files)
    for command in plan["commands"]:
        path = Path(command["argv"][6])
        protocol = supervisor.read_json(path)
        protocol["source_sha256"] = plan["source"]["files_sha256"]
        save(path, protocol)
    save(prepared / "execution-plan.json", plan)
    sha = pin(prepared)
    assert supervisor.launch(prepared / "execution-plan.json", sha, runtime, detached=True) == 0
    child = supervisor.read_json(runtime / "detached.json")["child"]
    try:
        deadline = time.monotonic() + 40
        state = {}
        while time.monotonic() < deadline:
            if (runtime / "progress.json").exists():
                state = supervisor.read_json(runtime / "progress.json")
                if state["status"] in {"completed", "failed", "interrupted"}:
                    break
            time.sleep(0.05)
        assert state.get("status") == "completed", state
        assert state["completed_attempts"] == 270
        assert state["supervisor"] == child
        assert child["pid"] == child["process_group_id"]
        os.waitpid(child["pid"], 0)
    finally:
        try:
            if supervisor.birth(child["pid"]) == child:
                os.kill(child["pid"], signal.SIGTERM)
                os.waitpid(child["pid"], 0)
        except (ProcessLookupError, ChildProcessError, supervisor.psutil.NoSuchProcess):
            pass
