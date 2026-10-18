"""Inert clock/accounting controls; no scientific inputs, models or assets."""

import json

import pytest

from chembounce_synthesis import phase_telemetry as pt


@pytest.fixture
def clock(monkeypatch):
    current = [100, 10]
    monkeypatch.setattr(pt.time, "monotonic_ns", lambda: current[0])
    monkeypatch.setattr(pt.time, "process_time_ns", lambda: current[1])
    return current


def test_fixed_arm_registry_and_no_external_cpu_claim():
    assert len(pt.BUCKETS) == 39 < 64
    assert [len(pt.PHASE_BUCKETS[p]) for p in ("native", "supply", "common")] == [12, 9, 18]
    snapshots = [pt.PhaseTelemetry(p).snapshot() for p in pt.PHASE_BUCKETS]
    assert sum(len(s["buckets"]) for s in snapshots) == 39
    for snapshot in snapshots:
        assert snapshot["identity"]["cpu_scope"].startswith("current process")
        assert "descendant" in snapshot["identity"]["cpu_scope"]
        assert snapshot["identity"]["process_start_ticks"] > 0
        assert snapshot["identity"]["boot_id"]
        assert snapshot["status"] == "closed"
        assert all(r["status"] == "not_reached" for r in snapshot["buckets"].values())
        assert all(r["inclusive_wall_ns"] is None for r in snapshot["buckets"].values())
        json.dumps(snapshot, allow_nan=False)


def test_nested_and_repeated_active_intervals_do_not_use_envelope(clock):
    timing = pt.PhaseTelemetry("native")
    with timing.span("native.worker"):
        clock[:] = [110, 12]
        with timing.span("native.api"):
            clock[:] = [120, 15]
            with timing.span("native.setup"):
                clock[:] = [130, 18]
            clock[:] = [180, 20]  # large disjoint idle gap, not setup duration
            with timing.span("native.setup"):
                clock[:] = [200, 24]
            clock[:] = [210, 25]
        clock[:] = [220, 26]
    result = timing.snapshot()
    worker, api, setup = [result["buckets"][f"native.{n}"] for n in ("worker", "api", "setup")]
    assert setup["inclusive_wall_ns"] == setup["exclusive_wall_ns"] == 30
    assert setup["last_end_wall_ns"] - setup["first_start_wall_ns"] == 80
    assert setup["inclusive_self_cpu_ns"] == setup["exclusive_self_cpu_ns"] == 7
    assert setup["calls"] == setup["completed"] == 2
    assert api["inclusive_wall_ns"] == 100 and api["exclusive_wall_ns"] == 70
    assert api["inclusive_self_cpu_ns"] == 13 and api["exclusive_self_cpu_ns"] == 6
    assert worker["inclusive_wall_ns"] == 120 and worker["exclusive_wall_ns"] == 20
    assert sum(r["exclusive_self_cpu_ns"] for r in (worker, api, setup)) == 16
    assert sum(r["exclusive_wall_ns"] for r in (worker, api, setup)) == 120
    assert setup["parents"] == ["native.api"]
    # Snapshots are copies, not mutable authority over counters.
    setup["parents"].append("forged")
    assert timing.snapshot()["buckets"]["native.setup"]["parents"] == ["native.api"]


@pytest.mark.parametrize(
    "exception", [RuntimeError("failure"), TimeoutError("deadline"), KeyboardInterrupt()]
)
def test_exception_closes_interrupted_scopes_and_propagates_exact_object(clock, exception):
    timing = pt.PhaseTelemetry("common")
    with pytest.raises(type(exception)) as caught:
        with timing.span("common.workflow"):
            with timing.span("common.catalog"):
                clock[:] = [110, 15]
                raise exception
    assert caught.value is exception
    result = timing.snapshot()
    assert result["status"] == "closed" and result["interrupted"]
    for name in ("common.workflow", "common.catalog"):
        assert result["buckets"][name]["status"] == "interrupted"
        assert result["buckets"][name]["interrupted"] == 1
    assert result["buckets"]["common.rules"]["status"] == "not_reached"


def test_cooperative_deadline_uses_existing_observation_not_new_clock(clock):
    timing = pt.PhaseTelemetry("common")
    with timing.span("common.route_search") as span:
        span.interrupt("workflow_deadline_observed")
        clock[:] = [111, 17]
    result = timing.snapshot()["buckets"]["common.route_search"]
    assert result["status"] == "interrupted"
    assert result["reason"] == "workflow_deadline_observed"
    assert result["inclusive_wall_ns"] == 11


def test_skipped_unavailable_not_reached_and_open_are_distinct(clock):
    timing = pt.PhaseTelemetry("native")
    timing.skip("native.reference_preflight")
    timing.unavailable("native.output", "terminated_export")
    with timing.span("native.worker"):
        result = timing.snapshot()
        assert result["status"] == "open" and result["active_depth"] == 1
        assert result["buckets"]["native.worker"]["status"] == "running"
    result = timing.snapshot()
    assert result["buckets"]["native.reference_preflight"]["status"] == "skipped_by_policy"
    assert result["buckets"]["native.output"]["status"] == "unavailable"
    assert result["buckets"]["native.generation"]["status"] == "not_reached"


def test_none_and_snapshot_read_no_clocks(monkeypatch):
    timing = pt.PhaseTelemetry("supply")

    def forbidden():
        pytest.fail("unexpected clock read")

    monkeypatch.setattr(pt.time, "monotonic_ns", forbidden)
    monkeypatch.setattr(pt.time, "process_time_ns", forbidden)
    with pt.phase_span(None, "native.generation"):
        pass
    pt.validate_telemetry(None, "native")
    assert timing.snapshot()["status"] == "closed"


@pytest.mark.parametrize("wrong", ["unknown", "common.catalog"])
def test_unknown_and_other_phase_buckets_fail_closed(wrong):
    timing = pt.PhaseTelemetry("native")
    with pytest.raises(ValueError):
        with timing.span(wrong):
            pass
    assert timing.snapshot()["status"] == "invalid"


def test_bounded_stack_recursion_and_order(clock):
    timing = pt.PhaseTelemetry("native")
    frames = [timing.span(f"native.{n}") for n in ("worker", "api", "setup", "generation", "allocation")]
    for frame in frames[:4]:
        frame.__enter__()
    with pytest.raises(ValueError, match="depth"):
        frames[4].__enter__()
    assert len(timing._stack) == 4
    timing = pt.PhaseTelemetry("native")
    first, second = timing.span("native.worker"), timing.span("native.api")
    first.__enter__()
    second.__enter__()
    with pytest.raises(ValueError, match="stack order"):
        first.__exit__(None, None, None)
    assert timing.snapshot()["status"] == "invalid"
    timing = pt.PhaseTelemetry("native")
    timing.span("native.worker").__enter__()
    with pytest.raises(ValueError, match="recursively"):
        timing.span("native.worker").__enter__()


@pytest.mark.parametrize("kind", ["wall", "cpu", "clock_exception", "noninteger"])
def test_clock_failures_never_clamp_or_publish_valid_durations(clock, monkeypatch, kind):
    timing = pt.PhaseTelemetry("native")
    span = timing.span("native.setup")
    span.__enter__()
    if kind == "wall":
        clock[0] -= 1
    elif kind == "cpu":
        clock[1] -= 1
    elif kind == "noninteger":
        clock[0] = True
    else:

        def broken():
            raise OSError("clock unavailable")

        monkeypatch.setattr(pt.time, "monotonic_ns", broken)
    with pytest.raises((ValueError, OSError)):
        span.__exit__(None, None, None)
    assert timing.snapshot()["status"] == "invalid"


@pytest.mark.parametrize("identity", ["pid", "thread"])
def test_process_thread_reuse_refused(monkeypatch, identity):
    timing = pt.PhaseTelemetry("native")
    if identity == "pid":
        monkeypatch.setattr(pt.os, "getpid", lambda: timing._pid + 1)
    else:
        monkeypatch.setattr(pt.threading, "get_ident", lambda: timing._thread + 1)
    with pytest.raises(ValueError, match="boundaries"):
        timing.span("native.worker").__enter__()


def test_no_custom_callback_context_or_cross_phase_accepted():
    for value in (object(), pt.PhaseTelemetry("supply")):
        with pytest.raises(ValueError, match="exact PhaseTelemetry"):
            pt.validate_telemetry(value, "native")


def test_marking_or_reusing_completed_spans_is_invalid(clock):
    timing = pt.PhaseTelemetry("native")
    span = timing.span("native.setup")
    with span:
        pass
    with pytest.raises(ValueError, match="active span"):
        span.interrupt("late")
    timing = pt.PhaseTelemetry("native")
    timing.skip("native.reference_preflight")
    with pytest.raises(ValueError, match="skipped"):
        timing.span("native.reference_preflight").__enter__()
