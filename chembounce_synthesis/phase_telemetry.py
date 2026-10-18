"""Fixed, opt-in phase counters; external ownership remains the budget authority.

Durations are sums of active intervals, never first/last timestamp envelopes.
Only the current process CPU is sampled. No files are written by this module.
"""

from contextlib import nullcontext
import os
from pathlib import Path
import threading
import time

VERSION = "bounded-phase-telemetry-v1"
MAX_DEPTH = 4
PHASE_BUCKETS = {
    "native": (
        "worker",
        "imports",
        "validation",
        "api",
        "setup",
        "generation",
        "reference_preflight",
        "allocation",
        "output",
        "trace_binding",
        "postflight",
        "publication",
    ),
    "supply": (
        "worker",
        "imports",
        "validation",
        "reference_admission",
        "provider_prepare",
        "provider_sampling",
        "syntax",
        "postflight",
        "publication",
    ),
    "common": (
        "worker",
        "imports",
        "validation",
        "workflow",
        "reference_setup",
        "reference_admission",
        "catalog",
        "rules",
        "evidence_prepare",
        "route_search",
        "forward_evidence",
        "screening",
        "refinement",
        "selection",
        "workflow_finalization",
        "collector",
        "postflight",
        "publication",
    ),
}
BUCKETS = tuple(f"{phase}.{name}" for phase, names in PHASE_BUCKETS.items() for name in names)
assert len(BUCKETS) == len(set(BUCKETS)) < 64


class PhaseTelemetry:
    """One process/phase per instance; at most 39 fixed buckets across an arm.

    Use a single instance in each of the native, supply and common workers.
    A closed snapshot describes coverage, not successful scientific completion.
    Missing snapshots after hard termination are unavailable, not zero-duration.
    """

    def __init__(self, phase):
        if phase not in PHASE_BUCKETS:
            raise ValueError("Unknown telemetry phase")
        self.phase = phase
        self._pid = os.getpid()
        self._thread = threading.get_ident()
        stat = Path("/proc/self/stat").read_text()
        self.identity = {
            "pid": self._pid,
            "process_start_ticks": int(stat[stat.rfind(")") + 2 :].split()[19]),
            "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
            "thread_ident": self._thread,
            "wall_clock": "monotonic_ns",
            "cpu_clock": "process_time_ns",
            "cpu_scope": "current process, all threads; excludes descendant processes",
        }
        self._stack = []
        self._last = None
        self._invalid = False
        self._records = {
            f"{phase}.{name}": {
                "status": "not_reached",
                "calls": 0,
                "completed": 0,
                "interrupted": 0,
                "inclusive_wall_ns": 0,
                "exclusive_wall_ns": 0,
                "inclusive_self_cpu_ns": 0,
                "exclusive_self_cpu_ns": 0,
                "first_start_wall_ns": None,
                "last_end_wall_ns": None,
                "parents": set(),
                "reason": None,
            }
            for name in PHASE_BUCKETS[phase]
        }

    def _fail(self, message):
        self._invalid = True
        raise ValueError(message)

    def _owner(self):
        if os.getpid() != self._pid or threading.get_ident() != self._thread:
            self._fail("Telemetry cannot cross process or thread boundaries")
        if self._invalid:
            raise ValueError("Telemetry context is invalid")

    def _sample(self):
        self._owner()
        try:
            sample = (time.monotonic_ns(), time.process_time_ns())
        except BaseException:
            self._invalid = True
            raise
        if any(type(value) is not int or value < 0 for value in sample):
            self._fail("Telemetry clocks must return nonnegative integer nanoseconds")
        if self._last is not None and any(a < b for a, b in zip(sample, self._last)):
            self._fail("Telemetry clock regressed")
        self._last = sample
        return sample

    def _record(self, name):
        self._owner()
        if name not in self._records:
            self._fail("Unknown or cross-phase telemetry bucket")
        return self._records[name]

    def span(self, name):
        return _Span(self, name)

    def _begin(self, name):
        record = self._record(name)
        if len(self._stack) >= MAX_DEPTH or any(frame["name"] == name for frame in self._stack):
            self._fail("Telemetry depth exceeded or bucket recursively entered")
        if record["status"] in {"skipped_by_policy", "unavailable"}:
            self._fail("Cannot run a skipped or unavailable telemetry bucket")
        wall, cpu = self._sample()
        record["calls"] += 1
        record["status"] = "running"
        if record["first_start_wall_ns"] is None:
            record["first_start_wall_ns"] = wall
        record["parents"].add(self._stack[-1]["name"] if self._stack else None)
        frame = {"name": name, "wall": wall, "cpu": cpu, "child_wall": 0, "child_cpu": 0}
        self._stack.append(frame)
        return frame

    def _end(self, frame, interrupted):
        self._owner()
        if not self._stack or self._stack[-1] is not frame:
            self._fail("Telemetry spans must close in stack order")
        wall, cpu = self._sample()
        elapsed_wall, elapsed_cpu = wall - frame["wall"], cpu - frame["cpu"]
        exclusive_wall = elapsed_wall - frame["child_wall"]
        exclusive_cpu = elapsed_cpu - frame["child_cpu"]
        if exclusive_wall < 0 or exclusive_cpu < 0:
            self._fail("Telemetry child durations exceed parent duration")
        self._stack.pop()
        if self._stack:
            self._stack[-1]["child_wall"] += elapsed_wall
            self._stack[-1]["child_cpu"] += elapsed_cpu
        record = self._records[frame["name"]]
        record["inclusive_wall_ns"] += elapsed_wall
        record["exclusive_wall_ns"] += exclusive_wall
        record["inclusive_self_cpu_ns"] += elapsed_cpu
        record["exclusive_self_cpu_ns"] += exclusive_cpu
        record["last_end_wall_ns"] = wall
        interrupted = interrupted or bool(frame.get("interruption_reason"))
        if frame.get("interruption_reason"):
            record["reason"] = frame["interruption_reason"]
        record["interrupted" if interrupted else "completed"] += 1
        record["status"] = "interrupted" if record["interrupted"] else "complete"

    def _mark(self, name, status, reason):
        record = self._record(name)
        if record["status"] != "not_reached":
            self._fail("Only an unreached telemetry bucket may be explicitly marked")
        if not isinstance(reason, str) or not reason or len(reason) > 160:
            self._fail("Telemetry reason must be a nonempty string of at most 160 characters")
        record.update(status=status, reason=reason)

    def skip(self, name, reason="skipped_by_policy"):
        self._mark(name, "skipped_by_policy", reason)

    def unavailable(self, name, reason):
        self._mark(name, "unavailable", reason)

    def snapshot(self):
        """Copy fixed counters without sampling clocks or serializing/writing files.

        Export/final snapshot serialization is charged by the external owner and
        is not retrospectively included in this snapshot. Active scopes have
        unfinished durations; consumers must not treat this as terminal evidence.
        """
        if os.getpid() != self._pid or threading.get_ident() != self._thread:
            self._fail("Telemetry snapshot belongs to another process or thread")
        records = {}
        for name, record in self._records.items():
            item = dict(record)
            item["parents"] = sorted(record["parents"], key=lambda value: value or "")
            if not record["calls"]:
                for field in (
                    "inclusive_wall_ns",
                    "exclusive_wall_ns",
                    "inclusive_self_cpu_ns",
                    "exclusive_self_cpu_ns",
                ):
                    item[field] = None
            records[name] = item
        return {
            "version": VERSION,
            "phase": self.phase,
            "identity": dict(self.identity),
            "status": "invalid" if self._invalid else "open" if self._stack else "closed",
            "interrupted": any(record["interrupted"] for record in self._records.values()),
            "active_depth": len(self._stack),
            "registry_size_per_arm": len(BUCKETS),
            "max_depth": MAX_DEPTH,
            "duration_semantics": "Sum of active intervals; timestamp envelopes are not durations",
            "budget_authority": "External owned-process terminal measurement; never sum nested inclusive buckets",
            "export_coverage": "Snapshot/export serialization and work outside spans are external-accounted only",
            "buckets": records,
        }


class _Span:
    def __init__(self, telemetry, name):
        self.telemetry, self.name, self.frame = telemetry, name, None

    def __enter__(self):
        if self.frame is not None:
            self.telemetry._fail("Telemetry span objects cannot be reused")
        self.frame = self.telemetry._begin(self.name)
        return self

    def interrupt(self, reason):
        self.telemetry._owner()
        if self.frame is None or not any(frame is self.frame for frame in self.telemetry._stack):
            self.telemetry._fail("Only an active span may be marked interrupted")
        if not isinstance(reason, str) or not reason or len(reason) > 160:
            self.telemetry._fail("Interruption reason must be a nonempty string of at most 160 characters")
        self.frame["interruption_reason"] = reason

    def __exit__(self, exc_type, exc, traceback):
        self.telemetry._end(self.frame, exc_type is not None)
        return False


def validate_telemetry(telemetry, phase):
    if telemetry is not None:
        if type(telemetry) is not PhaseTelemetry or telemetry.phase != phase:
            raise ValueError("Telemetry must be an exact PhaseTelemetry for this process phase")
        telemetry._owner()


def phase_span(telemetry, name):
    """Disabled calls allocate no records and read no clocks."""
    return nullcontext() if telemetry is None else telemetry.span(name)
