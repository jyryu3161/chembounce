"""Optional rejected-pose evidence; never a qualification or rescue path."""

import copy
from dataclasses import asdict, is_dataclass
import hashlib
import json
from pathlib import Path
import sys
import time

VERSION = "rejected-pose-diagnostic-v1"
ARTIFACT = "rejected_pose_diagnostics.json"
LIMITS = {
    "records": 8,
    "run_bytes": 2 * 1024 * 1024,
    "record_bytes": 256 * 1024,
    "atoms": 512,
    "capture_seconds": 1.0,
}
HEADER_RESERVE = 65536


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(encoded(value)).hexdigest()


def source_provenance(check=lambda: None):
    """Small source binding once per collector; no package/asset discovery."""
    root = Path(__file__).parent
    hashes = {}
    for name in ("rejected_pose_diagnostics.py", "selective_scoring.py", "similarity_constraints.py"):
        check()
        hashes[name] = hashlib.sha256((root / name).read_bytes()).hexdigest()
    return {
        "source_sha256": hashes,
        "python": sys.version,
        "rdkit": getattr(sys.modules.get("rdkit"), "__version__", None),
        "scope": "Loaded runtime labels and current source files, not independent runtime authentication.",
    }


class RejectedPoseDiagnostics:
    """One per workflow; private serialized records cannot be mutated by readers.

    All optional bookkeeping, failed capture, serialization and publication work
    is charged to one cooperative allowance. A native or filesystem operation
    already entered can overrun it; there is no hard interruption guarantee.
    """

    def __init__(
        self, *, provenance=None, clock=time.monotonic, started_monotonic=None, deadline_monotonic=None
    ):
        self.clock = clock
        started = clock() if started_monotonic is None else started_monotonic
        self._records = []
        self._bytes = self._calls = self._considered = self._completed = 0
        self._skips = {}
        self._not_applicable = 0
        self._capture_seconds = self._publication_seconds = 0.0
        self._published = False
        self._provenance_error = None
        self._deadline = deadline_monotonic

        def check_setup():
            if deadline_monotonic is not None and clock() >= deadline_monotonic:
                raise TimeoutError("workflow_deadline_exhausted")
            if clock() - started >= LIMITS["capture_seconds"]:
                raise TimeoutError("diagnostic_time_exhausted")

        try:
            check_setup()
            self._provenance = encoded(
                {**source_provenance(check_setup), "workflow_metadata": provenance or {}}
            )
            if len(self._provenance) > HEADER_RESERVE // 2:
                raise ValueError("provenance_limit")
        except Exception as exc:
            self._provenance = b"{}"
            self._provenance_error = "provenance_unavailable:" + type(exc).__name__
        self._elapsed = max(0.0, clock() - started)
        self._setup_seconds = self._elapsed

    def _reason(self, deadline, *, in_flight=0.0):
        deadlines = [value for value in (deadline, self._deadline) if value is not None]
        deadline = min(deadlines) if deadlines else None
        if self._provenance_error:
            return self._provenance_error
        if deadline is not None and self.clock() >= deadline:
            return "workflow_deadline_exhausted"
        if self._elapsed + in_flight >= LIMITS["capture_seconds"]:
            return "diagnostic_time_exhausted"
        return None

    def _skip(self, reason):
        self._skips[reason] = self._skips.get(reason, 0) + 1
        return {"version": VERSION, "status": "unavailable", "reason": reason}

    def charge_overhead(self, elapsed):
        """Charge caller-side opt-in dispatch/import work, without scientific work."""
        self._elapsed += max(0.0, elapsed)
        self._capture_seconds += max(0.0, elapsed)

    def start_call(self, identities, settings, policy, deadline, phase):
        started = self.clock()
        self._calls += 1
        call = _Call(
            self,
            identities,
            asdict(settings) if is_dataclass(settings) else settings,
            policy,
            deadline,
            phase,
            self._calls,
        )
        elapsed = max(0.0, self.clock() - started)
        self._elapsed += elapsed
        self._capture_seconds += elapsed
        return call

    def statistics(self):
        return {
            "version": VERSION,
            "limits": dict(LIMITS),
            "calls": self._calls,
            "call_scope": "Calls reaching valid selective preflight; dispatch overhead is included even if preflight fails.",
            "completed_calls": self._completed,
            "not_applicable_calls": self._not_applicable,
            "rejected_pairs_considered": self._considered,
            "captures": len(self._records),
            "record_bytes": self._bytes,
            "setup_seconds": self._setup_seconds,
            "capture_seconds": self._capture_seconds,
            "publication_seconds": self._publication_seconds,
            "total_observed_seconds": self._elapsed,
            "time_overrun": self._elapsed > LIMITS["capture_seconds"],
            "skips": dict(self._skips),
            "truncated": bool(self._skips),
        }

    def snapshot(self):
        """Defensive copy; complete records occur only in this shared artifact."""
        return {
            "schema_version": 1,
            **self.statistics(),
            "scope": "Diagnostic-only rejected poses, not qualification or exhaustive search. Global pharmacophore assignments were not retained.",
            "publication_completion_within_budget": None,
            "publication_timing_scope": "Final serialization/write completion is not certified by this artifact; consult caller diagnostics.",
            "provenance": json.loads(self._provenance),
            "records": [json.loads(record) for record in self._records],
        }


class _Call:
    def __init__(self, owner, identities, settings, policy, deadline, phase, call_id):
        self.owner = owner
        self.identities, self.settings, self.policy = copy.deepcopy((identities, settings, policy))
        self.deadline, self.phase, self.call_id = deadline, phase, call_id
        self.best, self.unavailable, self.finished = None, None, False

    def consider(
        self,
        pair_index,
        pair,
        aligned,
        witness,
        hit_info,
        candidate_energy,
        reference_energy,
        hit_tanimoto=None,
    ):
        """Retain an existing uniquely owned overlay, with strict first-on-tie."""
        started = self.owner.clock()
        try:
            self.owner._considered += 1
            reason = self.owner._reason(self.deadline)
            if reason:
                # Losing later coverage cannot produce a claimed global best.
                self.unavailable = reason
                self.best = None
                return
            if self.unavailable:
                return
            if len(self.owner._records) >= LIMITS["records"]:
                self.unavailable = "record_count_exhausted"
                return
            if self.best is None or pair["score"] > self.best[1]["score"]:
                if (
                    aligned.GetNumAtoms() > LIMITS["atoms"]
                    or len(witness["reference_pose"]["coordinates"]) > LIMITS["atoms"]
                ):
                    self.unavailable = "atom_limit"
                    self.best = None
                    return
                # Overlay owns its Mol; freeze the small existing audit structures.
                saved_pair, saved_witness, saved_info = copy.deepcopy((pair, witness, hit_info))
                self.best = (
                    pair_index,
                    saved_pair,
                    aligned,
                    saved_witness,
                    saved_info,
                    candidate_energy,
                    reference_energy,
                    hit_tanimoto,
                )
        except Exception as exc:
            self.unavailable = "capture_error:" + type(exc).__name__
            self.best = None
        finally:
            elapsed = max(0.0, self.owner.clock() - started)
            self.owner._elapsed += elapsed
            self.owner._capture_seconds += elapsed

    def finish(self, reason):
        """Serialize only after the scorer has decided final no-joint abstention."""
        started = self.owner.clock()
        try:
            if self.finished:
                return self.owner._skip("call_already_finished")
            self.finished = True
            self.owner._completed += 1
            if reason != "no_jointly_feasible_hit_pose":
                self.owner._not_applicable += 1
                return {"version": VERSION, "status": "not_applicable", "reason": reason}
            unavailable = self.unavailable or self.owner._reason(self.deadline)
            if unavailable:
                return self.owner._skip(unavailable)
            if self.best is None:
                return self.owner._skip("no_completed_rejected_pair")
            index, pair, aligned, witness, hit_info, ce, re, hit_tanimoto = self.best
            if (
                aligned.GetNumAtoms() > LIMITS["atoms"]
                or len(witness["reference_pose"]["coordinates"]) > LIMITS["atoms"]
            ):
                return self.owner._skip("atom_limit")
            from .similarity_constraints import pose_binding

            candidate = pose_binding(aligned, pair["candidate_conformer_id"])
            if candidate["sha256"] != witness["candidate_pose_sha256"]:
                raise ValueError("Already-computed candidate pose binding changed")
            reference = witness["reference_pose"]
            if reference["sha256"] != digest({k: v for k, v in reference.items() if k != "sha256"}):
                raise ValueError("Already-computed reference binding changed")
            if reference["conformer_id"] != pair["reference_conformer_id"]:
                raise ValueError("Reference conformer changed")
            if witness["policy_sha256"] != digest(self.policy):
                raise ValueError("Already-computed constraint policy changed")
            if (
                candidate["smiles"] != self.identities["product"]
                or reference["smiles"] != self.identities["hit"]
            ):
                raise ValueError("Diagnostic molecule identity changed")
            if pair["status"] != "rejected" or not pair["reasons"]:
                raise ValueError("Diagnostic requires a completed rejected pair")
            record = {
                "version": VERSION,
                "status": "diagnostic_only_rejected",
                "call_id": self.call_id,
                "selection": "highest_existing_completed_rejected_scalar_first_on_tie",
                "phase": self.phase,
                "pair_index": index,
                "identities": self.identities,
                "settings": self.settings,
                "policy": self.policy,
                "policy_sha256": digest(self.policy),
                "pair": pair,
                "hit_metrics": hit_info,
                "hit_tanimoto": hit_tanimoto,
                "candidate_relative_energy": ce,
                "reference_relative_energy": re,
                "candidate_pose": candidate,
                "witness": witness,
                "global_pharmacophore_assignment": {
                    "status": "not_retained",
                    "scope": "Existing global scalar/counts only; joint assignments cover polar recall and landmarks only, and may be infeasible/unassessed.",
                },
            }
            raw = encoded(record)
            if len(raw) > LIMITS["record_bytes"]:
                return self.owner._skip("record_bytes_exhausted")
            if self.owner._bytes + len(raw) > LIMITS["run_bytes"] - HEADER_RESERVE:
                return self.owner._skip("run_bytes_exhausted")
            unavailable = self.owner._reason(self.deadline, in_flight=max(0.0, self.owner.clock() - started))
            if unavailable:
                return self.owner._skip(unavailable)
            self.owner._records.append(raw)
            self.owner._bytes += len(raw)
            return {
                "version": VERSION,
                "status": "captured_not_published",
                "availability": "Requires successful shared artifact publication and matching record digest; unresolved otherwise.",
                "artifact": ARTIFACT,
                "record_index": len(self.owner._records) - 1,
                "record_sha256": hashlib.sha256(raw).hexdigest(),
                "call_id": self.call_id,
            }
        except Exception as exc:
            return self.owner._skip("serialization_error:" + type(exc).__name__)
        finally:
            self.best = None
            elapsed = max(0.0, self.owner.clock() - started)
            self.owner._elapsed += elapsed
            self.owner._capture_seconds += elapsed


def _write_bytes(path, data, check):
    """One optional atomic publication; never start an operation after known expiry."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    check()
    try:
        with temporary.open("xb") as stream:
            check()
            stream.write(data)
        check()
        temporary.replace(path)
    finally:
        # Best-effort temporary cleanup is not evidence publication.
        try:
            temporary.unlink(missing_ok=True)
        except OSError:
            pass


def publish(collector, out, *, deadline_monotonic=None, started_monotonic=None):
    """Optional atomic output. Ordinary errors preserve the scientific result.

    Interruption/termination propagates. Reference availability is conditional on
    this result, never predeclared in the scorer audit. The entire emitted byte
    string includes envelope/provenance and is bounded. Final persistence may
    overrun despite cooperative prechecks; no hard wall or durable timing claim.
    """
    started = collector.clock() if started_monotonic is None else started_monotonic
    result = {"version": VERSION, "status": "unavailable"}

    def check():
        reason = collector._reason(deadline_monotonic, in_flight=max(0.0, collector.clock() - started))
        if reason:
            raise TimeoutError(reason)

    try:
        if collector._published:
            raise ValueError("publication_already_attempted")
        collector._published = True
        check()
        data = encoded(collector.snapshot()) + b"\n"
        if len(data) > LIMITS["run_bytes"]:
            raise ValueError("diagnostic_artifact_byte_limit")
        artifact_sha256 = hashlib.sha256(data).hexdigest()
        check()
        _write_bytes(Path(out) / ARTIFACT, data, check)
        result.update(
            status="written",
            artifact=ARTIFACT,
            sha256=artifact_sha256,
            bytes=len(data),
            publication_completion_within_budget=None,
        )
    except Exception as exc:
        result["reason"] = type(exc).__name__ + ": " + str(exc)
    finally:
        elapsed = max(0.0, collector.clock() - started)
        collector._elapsed += elapsed
        collector._publication_seconds += elapsed
    result["publication_seconds"] = elapsed
    result["total_observed_seconds"] = collector._elapsed
    result["time_overrun"] = collector._elapsed > LIMITS["capture_seconds"]
    deadlines = [value for value in (deadline_monotonic, collector._deadline) if value is not None]
    result["deadline_overrun"] = bool(deadlines) and collector.clock() > min(deadlines)
    return result


def complete_summary(collector, out, summary, publication, *, workflow_started, deadline_monotonic=None):
    """One bounded optional update; the initial core summary already exists.

    Its saved elapsed endpoint precedes this final write. Post-return timing is
    returned for caller logging, never certified by a self-referential receipt.
    """
    started = collector.clock()
    result = {"status": "unavailable", "completion_within_budget": None}

    def check():
        reason = collector._reason(deadline_monotonic, in_flight=max(0.0, collector.clock() - started))
        if reason:
            raise TimeoutError(reason)

    try:
        if getattr(collector, "_summary_attempted", False):
            raise ValueError("summary_update_already_attempted")
        collector._summary_attempted = True
        check()
        updated = copy.deepcopy(summary)
        updated["rejected_pose_diagnostics"] = collector.statistics()
        updated["rejected_pose_diagnostics"]["export"] = publication
        updated["rejected_pose_diagnostics"]["summary_update_completion_within_budget"] = None
        updated["elapsed_seconds"] = max(0.0, collector.clock() - workflow_started)
        updated["elapsed_seconds_scope"] = (
            "Through optional sidecar publication; final summary persistence completion is unverified."
        )
        data = encoded(updated) + b"\n"
        check()
        _write_bytes(Path(out) / "summary.json", data, check)
        result["status"] = "written"
    except Exception as exc:
        result["reason"] = type(exc).__name__ + ": " + str(exc)
    finally:
        elapsed = max(0.0, collector.clock() - started)
        collector._elapsed += elapsed
        collector._publication_seconds += elapsed
    result["summary_update_seconds"] = elapsed
    result["total_observed_seconds"] = collector._elapsed
    result["time_overrun"] = collector._elapsed > LIMITS["capture_seconds"]
    deadlines = [value for value in (deadline_monotonic, collector._deadline) if value is not None]
    result["deadline_overrun"] = bool(deadlines) and collector.clock() > min(deadlines)
    return result
