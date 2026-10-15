"""Real process groups and signals exercise cancellation without GPU workloads."""

from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import psutil
import pytest

from benchmarks.research import runner
from benchmarks.research.protocol import digest, write_json


ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(os.name != "posix", reason="Runner uses POSIX process groups")


@pytest.fixture
def campaign(tmp_path):
    source = tmp_path / "source.txt"
    source.write_text("Synthetic signal test; no experimental claim.\n")
    dataset = {
        "schema_version": 1,
        "dataset_id": "signal-fixture",
        "scope": "synthetic_fixture",
        "sources": [
            {
                "id": "fixture",
                "url": "https://example.invalid/fixture",
                "artifact": str(source),
                "sha256": digest(source),
            }
        ],
        "compounds": [
            {
                "id": ident,
                "smiles": smiles,
                "target_id": "T",
                "label": label,
                "split": split,
                "source_id": "fixture",
                "assay_id": "fixture",
                "label_basis": "synthetic fixture",
            }
            for ident, smiles, label, split in [
                ("q", "CC(=O)Nc1ccccc1", "active", "train"),
                ("a", "CC(=O)Nc1ccncc1", "active", "test"),
                ("n", "CC(=O)Nc1ccsc1", "inactive", "test"),
            ]
        ],
        "queries": [{"id": "query", "compound_id": "q"}],
        "split": {"kind": "scaffold"},
    }
    dataset_path = tmp_path / "panel.json"
    write_json(dataset_path, dataset)
    worker = tmp_path / "worker.py"
    worker.write_text(
        """import json, os, signal, subprocess, sys, time
from pathlib import Path
request = json.loads(Path(sys.argv[1]).read_text())
base = Path(__file__).parent
if request['settings'].get('ignore_term'):
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
ready = base / 'grandchild-ready'
code = "import signal,sys,time; from pathlib import Path; signal.signal(signal.SIGTERM, signal.SIG_IGN); Path(sys.argv[1]).touch(); time.sleep(60)"
grandchild = subprocess.Popen([sys.executable, '-c', code, str(ready)])
while not ready.exists():
    time.sleep(.01)
# This is valid JSON containing a plausible candidate, but the worker is still
# running. It must never be accepted after interruption.
Path(sys.argv[2]).write_text(json.dumps({'candidates': [{'smiles': 'CC(=O)Nc1ccncc1'}], 'evaluations': 1, 'evaluation_unit': 'fixture'}))
(base / ('started-' + str(request['seed']))).write_text(json.dumps({'child': os.getpid(), 'grandchild': grandchild.pid}))
time.sleep(60)
"""
    )
    protocol = {
        "schema_version": 1,
        "seeds": [42, 43],
        "budget": {"wall_seconds": 45, "max_candidates": 5, "max_evaluations": 5},
        "methods": [
            {
                "id": "dummy",
                "kind": "external",
                "settings": {},
                "command": [sys.executable, str(worker), "{request}", "{response}"],
            }
        ],
    }
    protocol_path = tmp_path / "protocol.json"
    write_json(protocol_path, protocol)
    return dataset_path, protocol_path, protocol


def wait_for(predicate, timeout=15):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.02)
    pytest.fail("Timed out waiting for subprocess condition")


def terminal(pid):
    try:
        return psutil.Process(pid).status() == psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess:
        return True


def stop_own_group(process):
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait(timeout=5)


@pytest.mark.parametrize(
    "signum,library,ignore_term",
    [(signal.SIGTERM, False, False), (signal.SIGINT, False, False), (signal.SIGTERM, True, True)],
)
def test_signal_cleans_grandchildren_checkpoints_and_stops_campaign(
    campaign, tmp_path, signum, library, ignore_term
):
    dataset_path, protocol_path, protocol = campaign
    protocol["methods"][0]["settings"]["ignore_term"] = ignore_term
    write_json(protocol_path, protocol)
    output = tmp_path / "runs"
    arguments = ["--dataset", str(dataset_path), "--protocol", str(protocol_path), "--output", str(output)]
    if library:
        # Catching the library's cancellation exception must see the original
        # handlers restored, even after a second signal during forced cleanup.
        program = tmp_path / "library.py"
        program.write_text(
            """import json, signal, sys
from pathlib import Path
from benchmarks.research.runner import BenchmarkInterrupted, run_benchmark
calls = []
def original(signum, frame):
    calls.append(signum)
for signum in (signal.SIGTERM, signal.SIGINT):
    signal.signal(signum, original)
try:
    run_benchmark(sys.argv[1], sys.argv[2], sys.argv[3])
except BenchmarkInterrupted as exc:
    record = {'restored': all(signal.getsignal(s) is original for s in (signal.SIGTERM, signal.SIGINT)), 'original_calls': calls, 'attempts': len(exc.runs)}
    Path(sys.argv[4]).write_text(json.dumps(record))
    raise SystemExit(exc.exit_code)
raise AssertionError('Cancellation was not raised')
"""
        )
        command = [
            sys.executable,
            str(program),
            str(dataset_path),
            str(protocol_path),
            str(output),
            str(tmp_path / "restored.json"),
        ]
    else:
        command = [sys.executable, "-m", "benchmarks.research.runner", *arguments]
    environment = os.environ | {"PYTHONPATH": str(ROOT)}
    unrelated = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(60)"], start_new_session=True
    )
    pids = None
    with (tmp_path / "runner.log").open("w") as log:
        process = subprocess.Popen(
            command, cwd=ROOT, env=environment, stdout=log, stderr=log, start_new_session=True
        )
        try:
            ready = tmp_path / "started-42"
            wait_for(lambda: ready.exists() or process.poll() is not None)
            assert process.poll() is None, (tmp_path / "runner.log").read_text()
            pids = json.loads(ready.read_text())
            # Ensure the active process tree has been sampled, not just launched.
            time.sleep(0.15)
            start = time.monotonic()
            os.kill(process.pid, signum)
            if ignore_term:
                time.sleep(0.15)
                os.kill(process.pid, signal.SIGINT)
            assert process.wait(timeout=8) == 128 + signum
            assert time.monotonic() - start < 5
            wait_for(lambda: all(terminal(pid) for pid in pids.values()), timeout=3)
            assert unrelated.poll() is None
            assert not (tmp_path / "started-43").exists()
            assert not (output / "dummy/query-0000-seed-43").exists()
            runs = json.loads((output / "runs.json").read_text())
            assert len(runs) == 1
            report = runs[0]
            assert report["status"] == "interrupted"
            assert report["termination"]["signal_number"] == signum
            assert report["termination"]["received_phase"] == "worker_monitor"
            assert report["response"] == {}
            assert report["evaluation"]["candidates"] == []
            assert report["evaluation"]["metrics"]["candidate_count"] == 0
            assert report["resources"]["monitor_samples"] >= 1
            assert report["resources"]["sampled_peak_rss_mb"] > 0
            assert report["resources"]["sampled_cpu_seconds_lower_bound"] >= 0
            assert report["cleanup"]["sigkill_sent"]
            assert report["cleanup"]["leader_reaped"]
            attempt = output / "dummy/query-0000-seed-42"
            assert report["response_sha256"] == digest(attempt / "response.json")
            assert json.loads((attempt / "response.json").read_text())["candidates"]
            assert json.loads((attempt / "evaluation.json").read_text()) == report
            summary = json.loads((output / "summary.json").read_text())
            assert summary["dummy"]["status_counts"] == {"interrupted": 1}
            state = json.loads((output / "campaign.json").read_text())
            assert state["status"] == "interrupted" and state["attempted_runs"] == 1
            assert state["termination"]["exit_code"] == 128 + signum
            if library:
                assert json.loads((tmp_path / "restored.json").read_text()) == {
                    "restored": True,
                    "original_calls": [],
                    "attempts": 1,
                }
        finally:
            stop_own_group(process)
            if pids:
                try:
                    os.killpg(pids["child"], signal.SIGKILL)
                except ProcessLookupError:
                    pass
            stop_own_group(unrelated)


@pytest.mark.parametrize("stage", ["load_dataset", "load_protocol", "evidence_verifier", "public_request"])
def test_signal_during_setup_never_starts_attempt(campaign, tmp_path, monkeypatch, stage):
    dataset, protocol, _ = campaign
    original = getattr(runner, stage)

    def interrupt_after_setup(*args, **kwargs):
        result = original(*args, **kwargs)
        os.kill(os.getpid(), signal.SIGTERM)
        return result

    monkeypatch.setattr(runner, stage, interrupt_after_setup)
    output = tmp_path / "runs"
    handlers = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}
    with pytest.raises(runner.BenchmarkInterrupted) as caught:
        runner.run_benchmark(dataset, protocol, output)
    assert caught.value.runs == []
    assert not (output / "dummy").exists()
    assert json.loads((output / "runs.json").read_text()) == []
    state = json.loads((output / "campaign.json").read_text())
    assert state["status"] == "interrupted" and state["attempted_runs"] == 0
    assert all(signal.getsignal(s) == handler for s, handler in handlers.items())


@pytest.mark.parametrize("stage", ["evaluation", "checkpoint"])
def test_signal_after_completed_worker_preserves_outcome_and_stops_next(
    campaign, tmp_path, monkeypatch, stage
):
    dataset, protocol_path, protocol = campaign
    code = "import json,sys; json.dump({'candidates':[{'smiles':'CC(=O)Nc1ccncc1'}], 'evaluations':1, 'evaluation_unit':'fixture'},open(sys.argv[2],'w'))"
    protocol["methods"][0]["command"] = [sys.executable, "-c", code, "{request}", "{response}"]
    write_json(protocol_path, protocol)
    original = runner.evaluate_candidates if stage == "evaluation" else runner._write_checkpoint
    signalled = False

    def interrupt_after_completion(*args, **kwargs):
        nonlocal signalled
        result = original(*args, **kwargs)
        if not signalled and (stage == "evaluation" or Path(args[0]).name == "evaluation.json"):
            signalled = True
            os.kill(os.getpid(), signal.SIGTERM)
        return result

    monkeypatch.setattr(
        runner,
        "evaluate_candidates" if stage == "evaluation" else "_write_checkpoint",
        interrupt_after_completion,
    )
    output = tmp_path / "runs"
    with pytest.raises(runner.BenchmarkInterrupted) as caught:
        runner.run_benchmark(dataset, protocol_path, output)
    assert len(caught.value.runs) == 1
    report = json.loads((output / "runs.json").read_text())[0]
    assert report["status"] == "success" and report["evaluation"]["candidates"]
    assert "termination" not in report
    state = json.loads((output / "campaign.json").read_text())
    assert state["status"] == "interrupted"
    assert state["termination"]["received_phase"] == (
        "evaluation" if stage == "evaluation" else "attempt_checkpoint"
    )
    assert not (output / "dummy/query-0000-seed-43").exists()


@pytest.mark.parametrize("api", ["bounded_method", "run_benchmark"])
def test_library_requires_main_thread_before_work(tmp_path, api):
    output = tmp_path / "unused"
    arguments = ({}, {}, output) if api == "bounded_method" else ("missing", "missing", output)
    with ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(getattr(runner, api), *arguments)
        with pytest.raises(RuntimeError, match="main thread"):
            future.result()
    assert not output.exists()


def test_unrelated_monitor_exception_propagates_and_restores_handlers(tmp_path, monkeypatch):
    handlers = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}
    method = {
        "kind": "external",
        "command": [sys.executable, "-c", "import time; time.sleep(60)", "{request}", "{response}"],
    }
    request = {"budget": {"wall_seconds": 30, "max_candidates": 1, "max_evaluations": 1}}
    pids = []

    def broken_monitor(pid):
        pids.append(pid)
        raise RuntimeError("unrelated monitor failure")

    with monkeypatch.context() as patch:
        patch.setattr(runner.psutil, "Process", broken_monitor)
        with pytest.raises(RuntimeError, match="unrelated monitor failure"):
            runner.bounded_method(method, request, tmp_path / "broken")
    assert len(pids) == 1 and terminal(pids[0])
    assert all(signal.getsignal(s) == handler for s, handler in handlers.items())


@pytest.mark.parametrize("stage", ["monitor", "response_validation", "cleanup"])
def test_standalone_method_signal_scope_and_completion_boundary(tmp_path, monkeypatch, stage):
    code = (
        "import time; time.sleep(60)"
        if stage == "monitor"
        else "import json,sys; json.dump({'candidates':[], 'evaluations':0, 'evaluation_unit':'fixture'},open(sys.argv[2],'w'))"
    )
    method = {"kind": "external", "command": [sys.executable, "-c", code, "{request}", "{response}"]}
    request = {"budget": {"wall_seconds": 30, "max_candidates": 1, "max_evaluations": 1}}
    owner, name = {
        "monitor": (runner.psutil.Process, "children"),
        "response_validation": (runner.json, "loads"),
        "cleanup": (runner, "_stop_group"),
    }[stage]
    original = getattr(owner, name)
    handlers = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}
    signalled = False

    def interrupt(*args, **kwargs):
        nonlocal signalled
        result = original(*args, **kwargs)
        if not signalled:
            signalled = True
            os.kill(os.getpid(), signal.SIGTERM)
        return result

    with monkeypatch.context() as patch:
        patch.setattr(owner, name, interrupt)
        report = runner.bounded_method(method, request, tmp_path / "method")
    assert signalled
    assert report["termination"]["signal_number"] == signal.SIGTERM
    assert report["cleanup"]["leader_reaped"]
    assert all(signal.getsignal(s) == handler for s, handler in handlers.items())
    if stage == "cleanup":
        assert report["status"] == "no_candidates"
        assert report["response"]["evaluation_unit"] == "fixture"
    else:
        assert report["status"] == "interrupted" and report["response"] == {}
    if stage == "monitor":
        assert report["resources"]["monitor_samples"] >= 1


def test_success_restores_custom_handlers(campaign, tmp_path):
    dataset, protocol_path, protocol = campaign
    code = "import json,sys; json.dump({'candidates':[], 'evaluations':0, 'evaluation_unit':'fixture'},open(sys.argv[2],'w'))"
    protocol["methods"][0]["command"] = [sys.executable, "-c", code, "{request}", "{response}"]
    write_json(protocol_path, protocol)
    previous = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}

    def original(signum, frame):
        raise AssertionError("No signal expected")

    try:
        for signum in previous:
            signal.signal(signum, original)
        reports = runner.run_benchmark(dataset, protocol_path, tmp_path / "runs")
        assert len(reports) == 2 and all(r["status"] == "no_candidates" for r in reports)
        assert all(signal.getsignal(s) is original for s in previous)
        state = json.loads((tmp_path / "runs/campaign.json").read_text())
        assert state == {"schema_version": 1, "status": "completed", "attempted_runs": 2, "termination": None}
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)
