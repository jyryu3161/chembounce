"""Synthetic owner authority/timing checks; no real preparation/scientific work."""

import copy
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

from benchmarks.research import polar_capacity_preparation_owner as o

p, d, s = o.p, o.d, o.s
_fixture_spec = importlib.util.spec_from_file_location(
    "preparation_fixtures", Path(__file__).with_name("test_polar_capacity_prepare.py")
)
h = importlib.util.module_from_spec(_fixture_spec)
_fixture_spec.loader.exec_module(h)


def fixture(tmp_path, monkeypatch):
    """Tiny real git/hash/archive authority; constants/module paths are synthetic.

    No load_plan/validate_plan/tracked_inventory/archive validator is replaced.
    """
    spec, args = h.authority(tmp_path, monkeypatch)
    checker = Path(spec["checker"]["path"])
    script = checker / "benchmarks/research/polar_capacity_preparation_owner.py"
    script.write_text("# controlled owner source identity, never executed\n")
    subprocess.run(["git", "-C", str(checker), "add", "."], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(checker),
            "-c",
            "user.name=Synthetic",
            "-c",
            "user.email=test@example.invalid",
            "commit",
            "-qm",
            "owner",
        ],
        check=True,
    )
    spec["checker"]["commit"] = subprocess.check_output(
        ["git", "-C", str(checker), "rev-parse", "HEAD"], text=True
    ).strip()
    spec["checker"]["tree"] = subprocess.check_output(
        ["git", "-C", str(checker), "rev-parse", "HEAD^{tree}"], text=True
    ).strip()
    spec["checker"]["files"][str(script)] = spec["pins"][str(script)] = d.digest(script)
    monkeypatch.setattr(o, "__file__", str(script))
    monkeypatch.setattr(s, "__file__", str(checker / "benchmarks/research/polar_capacity_supervisor.py"))
    monkeypatch.setattr(d, "runtime", lambda: spec["contexts"]["checker"])
    monkeypatch.chdir(checker)
    for key, value in d.environment(checker).items():
        monkeypatch.setenv(key, value)
    auth = d.read(args.authorization)
    auth["external_owner"]["source"] = {"path": str(script), "sha256": d.digest(script)}
    del auth["literal_review"]
    del auth["spec_sha256"]
    payload = o.review_payload(spec, auth, args.spec, args.authorization)
    review_path = Path(args.authorization).with_name("literal-review.json")
    s.write(
        review_path, dict(schema=o.LITERAL_SCHEMA, status=o.REVIEW_STATUS, payload_sha256=o.binding(payload))
    )
    review = dict(path=str(review_path), sha256=d.digest(review_path))
    spec["reviews"].append(review)
    spec["pins"][str(review_path)] = review["sha256"]
    s.write(args.spec, spec)
    auth.update(spec_sha256=d.digest(args.spec), literal_review=review)
    s.write(args.authorization, auth)
    args.spec_sha256, args.authorization_sha256 = d.digest(args.spec), d.digest(args.authorization)
    args.anchor = 10.0
    return spec, auth, args, payload


def test_acyclic_real_authority_roundtrip(tmp_path, monkeypatch):
    spec, auth, args, payload = fixture(tmp_path, monkeypatch)
    assert o.authority_payload(spec, auth, args.spec, args.authorization) == payload
    assert o.authorize(args) == (spec, o.binding(payload))
    assert "spec_sha256" not in payload["authorization"]
    assert "literal_review" not in payload["authorization"]
    assert auth["literal_review"]["path"] not in payload["spec"]["pins"]
    assert not Path(spec["output"]).exists() and not Path(spec["run_output"]).exists()


@pytest.mark.parametrize(
    "mutation",
    [
        "budget",
        "output",
        "clearance",
        "source",
        "review_order",
        "invocation",
        "nan",
        "collision",
        "duplicate",
        "owner_collision",
    ],
)
def test_projection_changes_are_detected(tmp_path, monkeypatch, mutation):
    spec, auth, args, payload = fixture(tmp_path, monkeypatch)
    if mutation == "budget":
        auth["external_owner"]["watchdog_seconds"] += 1
    elif mutation == "output":
        spec["run_output"] += "-changed"
    elif mutation == "clearance":
        auth["clearance"]["sha256"] = "0" * 64
    elif mutation == "source":
        spec["checker"]["commit"] = "changed"
    elif mutation == "review_order":
        spec["reviews"].insert(0, copy.deepcopy(spec["reviews"][0]))
    elif mutation == "invocation":
        args.spec += "-changed"
    elif mutation == "nan":
        spec["x"] = float("nan")
    elif mutation == "duplicate":
        spec["reviews"].append(copy.deepcopy(auth["literal_review"]))
    else:
        auth["literal_review"] = (
            auth["source_review"] if mutation == "collision" else auth["external_owner"]["review"]
        )
    try:
        changed = o.authority_payload(spec, auth, args.spec, args.authorization)
    except ValueError:
        return
    assert o.binding(changed) != o.binding(payload)


@pytest.mark.parametrize(
    "mutation",
    ["source", "env", "review_hash", "review_payload", "dirty", "module", "output", "owner_contract"],
)
def test_authority_rejects_before_work(tmp_path, monkeypatch, mutation):
    spec, auth, args, _ = fixture(tmp_path, monkeypatch)
    if mutation == "source":
        auth["external_owner"]["source"]["sha256"] = "0" * 64
    elif mutation == "env":
        monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    elif mutation == "review_hash":
        Path(auth["literal_review"]["path"]).write_text("{}")
    elif mutation == "review_payload":
        spec["workloads"]["command_markers"].append("different")
    elif mutation == "dirty":
        (Path(spec["checker"]["path"]) / "untracked.py").write_text("# drift")
    elif mutation == "module":
        monkeypatch.setattr(s, "__file__", str(tmp_path / "foreign.py"))
    elif mutation == "output":
        args.output += "-different"
    else:
        auth["external_owner"]["scope"] = "children-only"
    s.write(args.spec, spec)
    auth["spec_sha256"] = args.spec_sha256 = d.digest(args.spec)
    s.write(args.authorization, auth)
    args.authorization_sha256 = d.digest(args.authorization)
    with pytest.raises((ValueError, d.IntegrityError)):
        o.authorize(args)
    assert not Path(spec["output"]).exists()


def complete_receipt(terminal, anchor=10.0):
    return dict(
        status="completed",
        returncode=0,
        cleanup=dict(complete=True, remaining=[]),
        wall_seconds=1.0,
        allocation_seconds=290.0,
        started_monotonic=anchor + 1,
        terminal_sha256=d.digest(terminal),
        identity=dict(pid=123, start_ticks=1),
        resources=dict(samples=1, gaps=[], sampled_tree_peak_rss_bytes=1024),
    )


def saved_prepared(spec, args, *, actual_prepare=False):
    base = Path(spec["output"])
    coord = base / "coordinator"
    handoff = coord / "handoff"
    handoff.mkdir(parents=True)
    queries = p.queries_from_design(spec)
    plan = p.assemble_plan(spec, queries, h.synthetic_discoveries(spec, queries))
    plan_path = coord / "candidate-plan.json"
    s.write(plan_path, plan)
    existing = handoff / "existing-prepare"
    if actual_prepare:
        assert (
            d.main(
                [
                    "prepare",
                    "--plan",
                    str(plan_path),
                    "--plan-sha256",
                    d.digest(plan_path),
                    "--output",
                    str(existing),
                ]
            )
            == 0
        )
    else:
        existing.mkdir()
        s.write(
            existing / "terminal.json",
            dict(
                status="prepared_not_authorized",
                plan_sha256=d.digest(plan_path),
                process_lifetime_ru_maxrss_kib=1,
            ),
        )
        s.write(
            existing / "ledger.json",
            dict(
                plan_sha256=d.digest(plan_path),
                rows=[
                    dict(index=i, status="not_started", request=d.request_for(plan, d.digest(plan_path), i))
                    for i in range(4)
                ],
            ),
        )
    s.write(
        handoff / "terminal.json", dict(status="prepared_not_authorized", plan_sha256=d.digest(plan_path))
    )
    s.write(handoff / "receipt.json", complete_receipt(handoff / "terminal.json", args.anchor))
    terminal = dict(
        status="candidate_requires_independent_review",
        candidate_plan=str(plan_path),
        candidate_sha256=d.digest(plan_path),
        study_launch_authorization_emitted=False,
        prepared_ledger=str(existing / "ledger.json"),
        pins_before_after_matched=True,
        deadline=args.anchor + 300,
    )
    s.write(coord / "terminal.json", terminal)
    receipt = complete_receipt(coord / "terminal.json", args.anchor)
    s.write(coord / "receipt.json", receipt)
    s.write(
        base / "preparation-receipt.json",
        dict(
            status="candidate_requires_independent_review",
            started_monotonic=args.anchor,
            spec_sha256=args.spec_sha256,
            success_allocation_seconds=300,
            watchdog_seconds=330,
            cleanup_seconds=30,
            study_launch_authorization_emitted=False,
            child=receipt,
            candidate=terminal,
            parent_elapsed_before_final_write=2.0,
        ),
    )
    return coord, handoff, existing


def test_actual_existing_prepare_and_deep_saved_artifact_mutations(tmp_path, monkeypatch):
    if os.environ.get("OWNER_COMPOSITION_CHILD") != "1":
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                f"{__file__}::test_actual_existing_prepare_and_deep_saved_artifact_mutations",
                "-q",
                "-p",
                "no:cacheprovider",
            ],
            env={**os.environ, "OWNER_COMPOSITION_CHILD": "1"},
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return
    spec, _, args, _ = fixture(tmp_path, monkeypatch)
    coord, handoff, existing = saved_prepared(spec, args, actual_prepare=True)
    evidence = o.validate_prepared(args, spec, args.anchor)
    assert evidence["planned_slots"] == 4 and len(evidence["files_sha256"]) == 8
    originals = {Path(path): Path(path).read_bytes() for path in evidence["files_sha256"]}
    for mutation in (
        "ledger_status",
        "ledger_request",
        "ledger_order",
        "plan_floor",
        "prep_anchor",
        "prep_late",
        "coordinator_path",
        "cleanup",
        "receipt_late",
        "handoff_terminal",
    ):
        for path, data in originals.items():
            path.write_bytes(data)
        if mutation.startswith("ledger"):
            path = existing / "ledger.json"
            value = d.read(path)
            if mutation == "ledger_status":
                value["rows"][0]["status"] = "completed"
            elif mutation == "ledger_request":
                value["rows"][0]["request"]["arm"] = "capacity"
            else:
                value["rows"].reverse()
            s.write(path, value)
        elif mutation == "plan_floor":
            path = coord / "candidate-plan.json"
            value = d.read(path)
            value["queries"]["attempt-0000"]["config"]["minimum_scores"]["hit_shape"] = 0.69
            s.write(path, value)
            value = d.read(coord / "terminal.json")
            value["candidate_sha256"] = d.digest(path)
            s.write(coord / "terminal.json", value)
        elif mutation.startswith("prep"):
            path = Path(spec["output"]) / "preparation-receipt.json"
            value = d.read(path)
            value[
                "started_monotonic" if mutation == "prep_anchor" else "parent_elapsed_before_final_write"
            ] = 999
            s.write(path, value)
        elif mutation == "coordinator_path":
            value = d.read(coord / "terminal.json")
            value["prepared_ledger"] += ".different"
            s.write(coord / "terminal.json", value)
        elif mutation == "handoff_terminal":
            s.write(handoff / "terminal.json", dict(status="prepared_not_authorized", plan_sha256="0" * 64))
            value = d.read(handoff / "receipt.json")
            value["terminal_sha256"] = d.digest(handoff / "terminal.json")
            s.write(handoff / "receipt.json", value)
        else:
            value = d.read(coord / "receipt.json")
            if mutation == "cleanup":
                value["cleanup"]["remaining"] = [{"pid": 999, "start_ticks": 1}]
            else:
                value["wall_seconds"] = 301
            s.write(coord / "receipt.json", value)
        # Rebind dependent receipts: rejection must include semantic checks,
        # not merely stale copies or a stale exported digest.
        value = d.read(coord / "receipt.json")
        value["terminal_sha256"] = d.digest(coord / "terminal.json")
        s.write(coord / "receipt.json", value)
        prep = d.read(Path(spec["output"]) / "preparation-receipt.json")
        prep.update(child=value, candidate=d.read(coord / "terminal.json"))
        s.write(Path(spec["output"]) / "preparation-receipt.json", prep)
        with pytest.raises((ValueError, d.IntegrityError)):
            o.validate_prepared(args, spec, args.anchor)
    assert not Path(spec["run_output"]).exists()


def adapter_fixture(tmp_path, monkeypatch):
    spec, _, args, payload = fixture(tmp_path, monkeypatch)
    owner, directory = o.paths(spec)
    directory.mkdir(parents=True)
    parent = dict(pid=100, start_ticks=5, ppid=1, state="S")
    child = dict(pid=101, start_ticks=6, ppid=100, state="R")
    monkeypatch.setattr(o.os, "getppid", lambda: 100)
    monkeypatch.setattr(o.os, "getpid", lambda: 101)
    monkeypatch.setattr(s, "identity", lambda pid: parent if pid == 100 else child)
    monkeypatch.setattr(s, "workload_check", lambda x: [])
    monkeypatch.setattr(o.time, "monotonic", lambda: 20)
    argv = o.adapter_argv(args, args.anchor)
    monkeypatch.setattr(sys, "argv", argv[1:])
    monkeypatch.setattr(
        o,
        "process_context",
        lambda pid: {
            "argv": o.owner_argv(args) if pid == 100 else o.adapter_argv(args, args.anchor),
            "cwd": spec["checker"]["path"],
        },
    )
    s.write(
        owner / "owner-receipt.json",
        dict(
            schema=o.OWNER_SCHEMA,
            status="running",
            spec_sha256=args.spec_sha256,
            authorization_sha256=args.authorization_sha256,
            payload_sha256=o.binding(payload),
            started_monotonic=args.anchor,
            adapter_argv=argv,
            owner_argv=o.owner_argv(args),
            environment=o.environment_record(spec["checker"]["path"]),
            environment_sha256=o.binding(dict(os.environ)),
            identity=parent,
        ),
    )
    s.write(
        directory / "receipt.json",
        dict(
            status="running",
            supervisor_identity=parent,
            identity=child,
            started_monotonic=args.anchor,
            allocation_seconds=300,
            watchdog_seconds=330,
            kill_grace_seconds=30,
            argv=argv,
            cwd=spec["checker"]["path"],
        ),
    )
    return spec, args, payload, owner, directory


@pytest.mark.parametrize(
    "mutation",
    [None, "parent", "self", "argv", "anchor", "env", "cwd", "reuse", "overlap", "deadline", "direct"],
)
def test_live_adapter_authority(tmp_path, monkeypatch, mutation):
    spec, args, payload, owner, directory = adapter_fixture(tmp_path, monkeypatch)
    if mutation in {"parent", "self", "argv", "anchor", "cwd"}:
        path = directory / "receipt.json"
        value = d.read(path)
        if mutation in {"parent", "self"}:
            value["supervisor_identity" if mutation == "parent" else "identity"]["start_ticks"] += 1
        elif mutation == "argv":
            value["argv"][-1] = "11.0"
        elif mutation == "anchor":
            value["started_monotonic"] += 1
        else:
            value["cwd"] = str(tmp_path)
        s.write(path, value)
    elif mutation == "env":
        monkeypatch.setenv("UNREVIEWED_EXTRA_ENV", "different")
    elif mutation == "reuse":
        Path(spec["output"]).mkdir()
    elif mutation == "overlap":
        monkeypatch.setattr(s, "workload_check", lambda x: [{"pid": 777}])
    elif mutation == "deadline":
        monkeypatch.setattr(o.time, "monotonic", lambda: 310)
    elif mutation == "direct":
        (owner / "owner-receipt.json").unlink()
    if mutation:
        with pytest.raises((ValueError, FileExistsError, FileNotFoundError)):
            o.adapter_authority(args, spec, o.binding(payload))
    else:
        assert o.adapter_authority(args, spec, o.binding(payload)) == 10.0


def test_clock_includes_supervisor_setup_delay():
    now = [299.0]
    clock = o.OwnerClock(0.0, lambda: now[0])
    now[0] = 301.0  # mkdir/pidfd setup happened before supervisor t0
    assert clock() == 0.0
    assert clock() == 301.0
    now[0] = 330.0
    assert clock() == 330.0
    now[0] = 360.0
    assert clock() == 360.0


@pytest.mark.parametrize("elapsed,write_delay,expected", [(299, 0, 0), (300, 0, 0), (301, 0, 1), (299, 2, 1)])
def test_owner_final_persistence_in_one_allocation(tmp_path, monkeypatch, elapsed, write_delay, expected):
    spec, args, _, owner, directory = adapter_fixture(tmp_path, monkeypatch)
    # Dispatch timing seam only; the genuine artifact validator has separate
    # isolated composition/deep mutation coverage above.
    for path in (directory / "receipt.json", owner / "owner-receipt.json"):
        path.unlink()
    directory.rmdir()
    owner.rmdir()
    monkeypatch.setattr(o, "authorize", lambda a: (spec, "synthetic-payload"))
    monkeypatch.setattr(d, "cpu_guard", lambda: None)
    monkeypatch.setattr(o, "active_identity", lambda pid: dict(pid=100, start_ticks=5))
    monkeypatch.setattr(
        o, "process_context", lambda pid: {"argv": o.owner_argv(args), "cwd": spec["checker"]["path"]}
    )
    now = [0.0]
    monkeypatch.setattr(o.time, "monotonic", lambda: now[0])
    real_write = s.write
    evidence = {"controlled": True}
    monkeypatch.setattr(o, "validate_prepared", lambda *a, **k: evidence)

    def writer(path, value):
        real_write(path, value)
        if Path(path).name == "owner-receipt.json" and value["status"] != "running":
            now[0] += write_delay

    def supervise(argv, **kw):
        assert (kw["allocation"], kw["watchdog"], kw["grace"]) == (300, 330, 30)
        now[0] = 5.0
        assert kw["clock"]() == 0.0 and kw["clock"]() == 5.0
        kw["check_workload"](spec["workloads"])
        kw["output"].mkdir()
        terminal = kw["output"] / "terminal.json"
        real_write(terminal, dict(status="prepared_candidate_validated", evidence=evidence))
        now[0] = elapsed
        return dict(
            status="completed",
            returncode=0,
            cleanup=dict(complete=True, remaining=[]),
            wall_seconds=elapsed,
            allocation_seconds=300,
            started_monotonic=0.0,
            identity=dict(pid=123, start_ticks=1),
            resources=dict(samples=1, gaps=[], sampled_tree_peak_rss_bytes=1024),
            terminal_sha256=d.digest(terminal),
        )

    monkeypatch.setattr(s, "write", writer)
    monkeypatch.setattr(s, "supervise", supervise)
    assert o._run(args, 0.0, {"expired": False}) == expected
    assert s.write is writer
    assert d.read(owner / "owner-receipt.json")["status"] == (
        "prepared_not_authorized" if expected == 0 else "allocation_exceeded"
    )


def test_actual_supervisor_caught_alarm_cannot_begin_final_writer(tmp_path, monkeypatch):
    original = s.write
    calls = []

    def spy(path, value):
        calls.append(value["status"])
        return original(path, value)

    monkeypatch.setattr(s, "write", spy)
    anchor = time.monotonic() - 359.85

    def delayed_check(blocked):
        time.sleep(0.5)  # actual SIGALRM before any child launch
        return []

    with p.administrative_alarm(anchor + 360) as expiry:
        with pytest.raises(o.CeilingExpired):
            with o.bounded_persistence(anchor, expiry):
                s.supervise(
                    [sys.executable, "-c", "raise AssertionError('must not launch')"],
                    cwd=tmp_path,
                    env=os.environ.copy(),
                    output=tmp_path / "supervisor",
                    allocation=300,
                    watchdog=330,
                    grace=30,
                    clock=o.OwnerClock(anchor),
                    check_workload=delayed_check,
                )
        assert expiry["expired"] is True
    assert s.write is spy and calls == ["running"]
    assert d.read(tmp_path / "supervisor/receipt.json")["status"] == "running"


@pytest.mark.parametrize("kind", ["expiry", "write_error"])
def test_persistence_guard_restores_original(tmp_path, kind):
    original = s.write
    expiry = {"expired": kind == "expiry"}
    with pytest.raises((o.CeilingExpired, FileNotFoundError)):
        with o.bounded_persistence(time.monotonic(), expiry):
            s.write(tmp_path / "absent/receipt.json", {})
    assert s.write is original


@pytest.mark.parametrize("elapsed", [300, 301, 361])
def test_exhausted_owner_does_not_read_authority(monkeypatch, elapsed):
    monkeypatch.setattr(o.time, "monotonic", lambda: elapsed)
    monkeypatch.setattr(o, "authorize", lambda args: pytest.fail("post-expiry authority read"))
    original = s.write
    assert o._run(None, 0.0, {"expired": False}) == 1
    assert s.write is original


def test_actual_supervisor_setup_delay_cannot_launch_after_allocation(tmp_path, monkeypatch):
    spec, _, args, _ = fixture(tmp_path, monkeypatch)
    now = [0.0]
    monkeypatch.setattr(o.time, "monotonic", lambda: now[0])
    monkeypatch.setattr(o, "authorize", lambda args: (spec, "controlled"))
    monkeypatch.setattr(d, "cpu_guard", lambda: None)
    monkeypatch.setattr(
        o, "process_context", lambda pid: {"argv": o.owner_argv(args), "cwd": spec["checker"]["path"]}
    )
    monkeypatch.setattr(s, "workload_check", lambda blocked: [])
    original = s.pidfd_open

    def delayed_setup(pid):
        result = original(pid)
        now[0] = 301.0
        return result

    monkeypatch.setattr(s, "pidfd_open", delayed_setup)
    monkeypatch.setattr(s.subprocess, "Popen", lambda *a, **k: pytest.fail("late child launch"))
    assert o._run(args, 0.0, {"expired": False}) == 1
    owner, directory = o.paths(spec)
    receipt = d.read(directory / "receipt.json")
    assert receipt["started_monotonic"] == 0 and receipt["wall_seconds"] == 301
    assert receipt["identity"] is None and receipt["cleanup"]["complete"] is True
    assert d.read(owner / "owner-receipt.json")["status"] == "allocation_exceeded"


def test_actual_supervisor_inherits_absolute_watchdog_and_cleanup(tmp_path):
    ready = tmp_path / "ready"
    script = "import signal,time,pathlib,sys; signal.signal(signal.SIGTERM,signal.SIG_IGN); pathlib.Path(sys.argv[1]).write_text('ready'); time.sleep(10)"
    ticks = iter([329.9, 330.0, 359.9, 359.9, 360.0, 360.0])  # loop reads before AND after /proc sampling
    initial = True

    def actual_ticks():
        nonlocal initial
        if initial:
            initial = False
            deadline = time.monotonic() + 3
            while not ready.exists() and time.monotonic() < deadline:
                time.sleep(0.01)
            assert ready.exists(), "controlled child did not initialize"
        time.sleep(0.03)  # let real TERM/KILL/reap operations progress
        return next(ticks, 360.0)

    receipt = s.supervise(
        [sys.executable, "-c", script, str(ready)],
        cwd=tmp_path,
        env=os.environ.copy(),
        output=tmp_path / "phase",
        allocation=300,
        watchdog=330,
        grace=30,
        clock=o.OwnerClock(0.0, actual_ticks),
        check_workload=lambda blocked: [],
    )
    assert receipt["status"] == "watchdog_timeout"
    assert receipt["started_monotonic"] == 0.0
    assert receipt["cleanup"]["term_monotonic"] == 330.0
    assert receipt["cleanup"]["kill_monotonic"] == 359.9
    assert receipt["wall_seconds"] == 360.0
    assert receipt["cleanup"]["complete"] is True and receipt["cleanup"]["remaining"] == []


def test_live_supervised_adapter_with_genuine_synthetic_prepare(tmp_path, monkeypatch):
    """Actual owner/adapter processes and receipts; synthetic source/runtime authority.

    The fixed adapter backend is replaced only in this temporary script, to call
    genuine existing d.main prepare on tiny pinned inputs instead of importing
    scientific libraries/discovering production runtimes. Live ownership,
    commands, environment, clocks, supervisor and final validators are actual.
    """
    real_owner = Path(o.__file__).resolve()
    real_root = real_owner.parents[2]
    test_file = Path(__file__).resolve()
    spec, auth, args, _ = fixture(tmp_path, monkeypatch)
    script = Path(o.__file__)
    script.write_text(f"""import sys, importlib.util
from pathlib import Path
sys.path.insert(0, {str(real_root)!r})
from benchmarks.research import polar_capacity_preparation_owner as o
p,d,s=o.p,o.d,o.s
spec=d.read(sys.argv[sys.argv.index('--spec')+1])
root=Path(spec['checker']['path'])
o.__file__=__file__
p.__file__=str(root/'benchmarks/research/polar_capacity_prepare.py')
d.__file__=str(root/'benchmarks/research/polar_capacity_development.py')
s.__file__=str(root/'benchmarks/research/polar_capacity_supervisor.py')
d.ARMS={{role: source['commit'] for role,source in spec['sources'].items()}}
d.DESIGN_SHA=spec['design']['sha256'];d.GATE_SHA=spec['gate']['sha256']
d.required_runtime_files=lambda:{{spec['contexts']['checker']['executable'],spec['generation_source']['path']}}
d.runtime=lambda:spec['contexts']['checker']
if sys.argv[1]=='adapter':
    fixture_spec=importlib.util.spec_from_file_location('owner_composition_fixtures',{str(test_file)!r})
    h=importlib.util.module_from_spec(fixture_spec);fixture_spec.loader.exec_module(h)
    def synthetic_prepare(args,anchor,expiry):
        h.saved_prepared(spec,args,actual_prepare=True)
        return 0
    p._run=synthetic_prepare
raise SystemExit(o.main())
""")
    checker = Path(spec["checker"]["path"])
    subprocess.run(["git", "-C", str(checker), "add", "."], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(checker),
            "-c",
            "user.name=Synthetic",
            "-c",
            "user.email=test@example.invalid",
            "commit",
            "-qm",
            "controlled live adapter",
        ],
        check=True,
    )
    spec["checker"]["commit"] = subprocess.check_output(
        ["git", "-C", str(checker), "rev-parse", "HEAD"], text=True
    ).strip()
    spec["checker"]["tree"] = subprocess.check_output(
        ["git", "-C", str(checker), "rev-parse", "HEAD^{tree}"], text=True
    ).strip()
    spec["checker"]["files"][str(script)] = spec["pins"][str(script)] = d.digest(script)
    auth["external_owner"]["source"]["sha256"] = d.digest(script)
    review = auth.pop("literal_review")
    auth.pop("spec_sha256")
    spec["reviews"].remove(review)
    del spec["pins"][review["path"]]
    payload = o.review_payload(spec, auth, args.spec, args.authorization)
    s.write(
        review["path"],
        dict(schema=o.LITERAL_SCHEMA, status=o.REVIEW_STATUS, payload_sha256=o.binding(payload)),
    )
    review["sha256"] = d.digest(review["path"])
    spec["reviews"].append(review)
    spec["pins"][review["path"]] = review["sha256"]
    s.write(args.spec, spec)
    auth.update(literal_review=review, spec_sha256=d.digest(args.spec))
    s.write(args.authorization, auth)
    args.spec_sha256, args.authorization_sha256 = d.digest(args.spec), d.digest(args.authorization)
    result = subprocess.run(
        o.owner_argv(args),
        cwd=checker,
        env=d.environment(checker),
        capture_output=True,
        text=True,
        timeout=30,
    )
    owner, directory = o.paths(spec)
    receipt = d.read(owner / "owner-receipt.json")
    assert result.returncode == 0, (
        result.stdout,
        result.stderr,
        receipt,
        (directory / "stderr.log").read_text(),
    )
    assert receipt["status"] == "prepared_not_authorized"
    assert receipt["child"]["cleanup"]["complete"] is True
    assert receipt["evidence"]["planned_slots"] == 4
    assert not Path(spec["run_output"]).exists()


@pytest.mark.parametrize(
    "stage", ["before", "after", "expiry_before", "expiry_after", "supervisor_swallowed"]
)
def test_owner_interruptions_never_certify_terminal_success(tmp_path, monkeypatch, stage):
    spec, _, args, _ = fixture(tmp_path, monkeypatch)
    owner, directory = o.paths(spec)
    monkeypatch.setattr(o, "authorize", lambda args: (spec, "controlled"))
    monkeypatch.setattr(d, "cpu_guard", lambda: None)
    monkeypatch.setattr(
        o, "process_context", lambda pid: {"argv": o.owner_argv(args), "cwd": spec["checker"]["path"]}
    )
    monkeypatch.setattr(s, "workload_check", lambda blocked: [])
    monkeypatch.setattr(o.time, "monotonic", lambda: 1.0)
    evidence = {"controlled": True}
    monkeypatch.setattr(o, "validate_prepared", lambda *a, **kw: evidence)
    expiry = {"expired": False}
    real_write = s.write
    attempts = []

    def write(path, value):
        assert not expiry["expired"], "new write attempted after sticky expiry"
        attempts.append(value["status"])
        if value["status"] == "prepared_not_authorized":
            if stage.endswith("after"):
                real_write(path, value)
            if stage.startswith("expiry"):
                expiry["expired"] = True
                raise TimeoutError("controlled expiry during persistence")
            raise KeyboardInterrupt("controlled persistence interruption")
        real_write(path, value)

    def supervise(*a, **kw):
        directory.mkdir()
        terminal = directory / "terminal.json"
        real_write(terminal, dict(status="prepared_candidate_validated", evidence=evidence))
        if stage == "supervisor_swallowed":
            expiry["expired"] = True
        return dict(
            status="completed",
            returncode=0,
            cleanup=dict(complete=True, remaining=[]),
            wall_seconds=1.0,
            started_monotonic=0.0,
            allocation_seconds=300,
            identity=dict(pid=123, start_ticks=1),
            resources=dict(samples=1, gaps=[], sampled_tree_peak_rss_bytes=1024),
            terminal_sha256=d.digest(terminal),
        )

    monkeypatch.setattr(s, "write", write)
    monkeypatch.setattr(s, "supervise", supervise)
    if stage in {"before", "after"}:
        with pytest.raises(KeyboardInterrupt):
            o._run(args, 0.0, expiry)
        assert attempts == ["running", "prepared_not_authorized", "persistence_interrupted"]
        assert d.read(owner / "owner-receipt.json")["status"] == "persistence_interrupted"
    else:
        assert o._run(args, 0.0, expiry) == 1
        assert attempts == (
            ["running"] if stage == "supervisor_swallowed" else ["running", "prepared_not_authorized"]
        )
        # Already-persisted bytes are not retroactively cancelled. External exit
        # remains nonzero and is required even if a file says prepared.
        assert d.read(owner / "owner-receipt.json")["status"] == (
            "prepared_not_authorized" if stage == "expiry_after" else "running"
        )
    assert s.write is write


@pytest.mark.parametrize("mutation", ["unsafe", "gap", "samples", "peak", "nonfinite", "boolean_return"])
def test_completed_label_does_not_waive_cleanup_or_resources(tmp_path, mutation):
    terminal = tmp_path / "terminal.json"
    s.write(terminal, {})
    receipt = complete_receipt(terminal)
    assert o.complete(receipt, 290.0)
    if mutation == "unsafe":
        receipt["cleanup"]["complete"] = False
    elif mutation == "gap":
        receipt["resources"]["gaps"] = [{"error": "AccessDenied"}]
    elif mutation == "samples":
        receipt["resources"]["samples"] = 0
    elif mutation == "peak":
        receipt["resources"]["sampled_tree_peak_rss_bytes"] = None
    elif mutation == "nonfinite":
        receipt["wall_seconds"] = float("nan")
    else:
        receipt["returncode"] = False
    assert not o.complete(receipt, 290.0)
