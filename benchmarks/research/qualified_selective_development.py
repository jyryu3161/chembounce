"""Fixed 197-pair development evaluation; no fitting, retries or model calls.

``spec`` prints the literal protocol/declaration requirements. ``prepare`` verifies
an externally hash-frozen declaration and writes all 591 slots before scoring.
``run`` requires the externally frozen plan digest and a new output directory.
The reviewed selective scorer must be present; no legacy-scoring fallback exists.
"""

import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import asdict
from fractions import Fraction
import importlib.abc
import json
import os
from pathlib import Path
import resource
import signal
import sys
import time

from benchmarks.research import qualified_calibration as base

SEEDS = (42, 43, 44)
FLOORS = {"hit_tanimoto": 0.35, "hit_shape": 0.70, "hit_pharmacophore": 0.50}
COUNTS = {"pairs": 197, "active": 111, "inactive": 86, "targets": 6, "queries": 12, "trials": 591}
INPUTS = {
    "panel": "464e71b8bc54ad9002a8ff02e08b249272ae43ce17a4aee3d53c394642439d33",
    "pairs": "bb96b91767da475f313029e2e6a6432e774b1c3f664d7638e7e0df5c988d110f",
    "old_scores": "4b685fda706b9f0ef337f8f2b149a41ea3a89747fac1d6fe8246018b5368fb2b",
    "old_protocol": "6e28e5668d3cd98ebd0caddceec25769407ccca28f53502f6b719fac1318b5a5",
    "identity_manifest": "c6c28f2f00ec0fe4c306ac21179560784bf21ddf89c14131f4a589c95fb7caf9",
    "implementation_decision": "a8a2fb000ca754a84ed1ed509654bd1aad9520811c040ce2e51fa2316b52a987",
    "assignment_decision": "fc0edef0ecc58ae44f5db75f023a88151ffbfc44cddfd05ce4bd11404c7e9a2c",
    "candidate_plan": "71c3bd991fc5713e98c2a6cbf142f59650d4af8945e5990c9336af10d10c73f0",
}
SPEC = {
    "schema": "qualified-selective-development-v2",
    "scope": "previously_inspected_fixed_development_pairs_only",
    "settings": base.SETTINGS,
    "seeds": list(SEEDS),
    "minimum_scores": FLOORS,
    "minimum_polar_reference_recall": 0.75,
    "landmarks": [],
    "core_smiles": None,
    "pose_selection_policy": "joint-feasible-hit-pose-v2",
    "max_hit_pairs_per_trial": 64,
    "max_features": 256,
    "max_landmarks": 32,
    "counts": COUNTS,
    "wall_seconds": 14400,
    "deadline_scope": "Cooperative between bounded operations; individual RDKit C calls are not interruptible. External process wall/memory containment must be frozen before execution.",
    "workers": 1,
    "fit": False,
    "retries": 0,
    "reference_policy": "query_then_seed_then_pair; one prewarm per block; unchanged bounded128 ensemble LRU cleared before/after run; immutable reference digest checked",
    "embedding_bounds": {
        "calls": 627,
        "requested_candidate_conformers": 4728,
        "requested_reference_conformers": 288,
    },
    "old_pose_scope": "saved aligned candidate; deterministic regenerated same-seed reference, matched ID/energy/diagnostics; no original saved reference coordinates or realignment",
    "metrics": "all original pairs including failed/unstarted; equal queries within target then equal targets; old-joint comparison coverage separate",
}
EXTRA_SOURCE = (
    "benchmarks/research/qualified_selective_development.py",
    "benchmarks/research/qualified_calibration.py",
    "benchmarks/research/assay_matched.py",
    "benchmarks/research/assay_metadata.py",
    "benchmarks/research/import_chembl.py",
    "benchmarks/research/protocol.py",
)
FORBIDDEN_IMPORTS = {"torch", "transformers", "safetensors", "tensorflow", "jax", "huggingface_hub"}


def read(path):
    return json.loads(
        Path(path).read_text(), parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value))
    )


def sources():
    root = Path(__file__).resolve().parents[2]
    files = set(EXTRA_SOURCE) | {
        str(p.relative_to(root)) for p in (root / "chembounce_synthesis").glob("*.py")
    }
    files.update(
        {"chembounce_synthesis/selective_scoring.py", "chembounce_synthesis/similarity_constraints.py"}
    )
    return {name: base.file_hash(root / name) for name in sorted(files)}


def verify_declaration(declaration):
    if (
        declaration.get("schema") != "qualified-selective-development-declaration-v2"
        or declaration.get("spec") != SPEC
    ):
        raise ValueError("Declaration changes the fixed engineering protocol")
    if set(declaration.get("inputs", {})) != set(INPUTS):
        raise ValueError("Missing or extra fixed input bindings")
    for key, expected in INPUTS.items():
        item = declaration["inputs"][key]
        if (
            set(item) != {"path", "sha256"}
            or item["sha256"] != expected
            or base.file_hash(item["path"]) != expected
        ):
            raise ValueError("Fixed input changed: " + key)
    if declaration.get("source_sha256") != sources() or declaration.get("runtime") != base.runtime():
        raise ValueError("Source or scientific runtime differs from declaration")
    runtime = declaration["runtime"]
    old_protocol = read(declaration["inputs"]["old_protocol"]["path"])
    for name in ("conformers.py", "common.py", "identity.py"):
        key = "chembounce_synthesis/" + name
        if old_protocol["source_sha256"][key] != declaration["source_sha256"][key]:
            raise ValueError("Old reference numerical kernel differs: " + name)
    if (
        old_protocol["runtime"] != runtime
        or old_protocol["calibration_spec"]["conformers"] != base.SETTINGS
        or old_protocol["calibration_spec"]["seeds"] != list(SEEDS)
    ):
        raise ValueError("Old reference runtime/settings provenance differs")
    if runtime["python"] != "3.10.12" or runtime["packages"] != {
        "rdkit": "2025.9.6",
        "numpy": "1.26.4",
        "scipy": "1.13.1",
    }:
        raise ValueError("Reviewed scientific runtime required")
    if os.environ.get("PYTHONDONTWRITEBYTECODE") != "1":
        raise ValueError("Bytecode writes must be disabled")
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "" or any(
        os.environ.get(k) != "1"
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")
    ):
        raise ValueError("CPU-only, one-thread environment required")


def ordered_slots(pairs):
    if len({p["pair_id"] for p in pairs}) != len(pairs):
        raise ValueError("Duplicate pair IDs")
    groups = defaultdict(list)
    for pair in pairs:
        groups[base.query_key(pair)].append(pair)
    slots = []
    for query in sorted(groups):
        if len({r["query_smiles"] for r in groups[query]}) != 1:
            raise ValueError("Query group has conflicting exact states")
        for seed in SEEDS:
            for pair in sorted(groups[query], key=lambda row: row["pair_id"]):
                slots.append(
                    {"index": len(slots), "pair_id": pair["pair_id"], "seed": seed, "query": list(query)}
                )
    return slots


def prepare(declaration_path, declaration_sha256, output):
    if base.file_hash(declaration_path) != declaration_sha256:
        raise ValueError("Declaration SHA mismatch")
    declaration = read(declaration_path)
    verify_declaration(declaration)
    paths = {key: value["path"] for key, value in declaration["inputs"].items()}
    pairs = read(paths["pairs"])
    if base.counts(pairs) != COUNTS:
        raise ValueError("Fixed 197-pair denominator changed")
    base.verify_inventory(pairs, read(paths["identity_manifest"]))
    # The panel/pairs bytes are fixed independently; verify their projection too.
    from benchmarks.research.assay_matched import control_pairs

    projected = {p["pair_id"]: p for p in control_pairs(read(paths["panel"]))}
    if set(projected) != {p["pair_id"] for p in pairs} or any(
        any(p.get(k) != v for k, v in projected[p["pair_id"]].items()) for p in pairs
    ):
        raise ValueError("Fixed panel and pair identities disagree")
    old = read(paths["old_scores"])
    if len(old) != len(pairs) or len({r["pair_id"] for r in old}) != len(pairs):
        raise ValueError("Old comparison denominator changed")
    by_id = {r["pair_id"]: r for r in old}
    for pair in pairs:
        row = by_id[pair["pair_id"]]
        if any(row.get(k) != v for k, v in pair.items()) or [t["seed"] for t in row["trials"]] != list(SEEDS):
            raise ValueError("Old chosen-pose identity/seed binding mismatch")
    plan = {
        "schema": "qualified-selective-development-plan-v2",
        "declaration": declaration,
        "declaration_sha256": declaration_sha256,
        "declaration_path": str(Path(declaration_path).resolve()),
        "declaration_canonical_sha256": base.json_hash(declaration),
        "pairs": pairs,
        "slots": ordered_slots(pairs),
    }
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    base.save(output / "plan.json", plan)
    base.save(
        output / "preparation.json",
        {
            "status": "plan_frozen_no_scoring",
            "plan_sha256": base.file_hash(output / "plan.json"),
            "counts": COUNTS,
        },
    )
    verify_declaration(declaration)
    return plan


def old_global(row):
    return all(
        t["status"] == "success"
        and all(
            v >= FLOORS[k]
            for k, v in {
                "hit_tanimoto": row["hit_tanimoto"],
                "hit_shape": t["info"]["hit"]["shape"],
                "hit_pharmacophore": t["info"]["hit"]["pharmacophore"],
            }.items()
        )
        for t in row["trials"]
    )


def acceptance_rates(pairs, predictions, label, topology=False):
    by_query = defaultdict(list)
    for pair in pairs:
        if pair["label"] == label and (not topology or pair["topology_change"]):
            by_query[base.query_key(pair)].append(pair)
    queries, targets = [], defaultdict(list)
    for key, rows in sorted(by_query.items()):
        count = sum(predictions.get(r["pair_id"], False) for r in rows)
        value = Fraction(count, len(rows))
        queries.append(
            {
                "query": list(key),
                "accepted": count,
                "denominator": len(rows),
                "rate": float(value),
                "rational": str(value),
            }
        )
        targets[key[0]].append(value)
    values = {key: sum(v, Fraction()) / len(v) for key, v in sorted(targets.items())}
    macro = sum(values.values(), Fraction()) / len(values) if values else None
    return {
        "queries": queries,
        "targets": {key: {"rate": float(v), "rational": str(v)} for key, v in values.items()},
        "macro": float(macro) if macro is not None else None,
        "macro_rational": str(macro) if macro is not None else None,
        "denominator": sum(len(v) for v in by_query.values()),
        "supported_targets": len(values),
    }


def summarize(pairs, slots, results, old_rows):
    if set(results) != set(range(len(slots))) or any(
        row["slot"] != slots[index] for index, row in results.items()
    ):
        raise ValueError("Missing, extra or rebound planned result slots")
    grouped = defaultdict(list)
    for row in results.values():
        grouped[row["slot"]["pair_id"]].append(row)
    new, old_joint, assessed = {}, {}, {}
    for pair in pairs:
        ident = pair["pair_id"]
        rows = grouped[ident]
        if sorted(r["slot"]["seed"] for r in rows) != list(SEEDS):
            raise ValueError("Missing or duplicate required seed")
        new[ident] = all(r["status"] == "success" for r in rows)
        assessed[ident] = all(r.get("old_comparison", {}).get("status") == "assessed" for r in rows)
        old_joint[ident] = assessed[ident] and all(r["old_comparison"]["passed"] for r in rows)
    original = {r["pair_id"]: old_global(r) for r in old_rows}
    rates = {}
    for name, predictions in (
        ("new_joint", new),
        ("old_global_floors_only", original),
        ("old_chosen_pose_joint", old_joint),
    ):
        rates[name] = {
            "active_acceptance": acceptance_rates(pairs, predictions, "active"),
            "inactive_acceptance": acceptance_rates(pairs, predictions, "inactive"),
            "topology_active_acceptance": acceptance_rates(pairs, predictions, "active", True),
            "accepted_pairs": sum(predictions.values()),
        }
    return {
        "rates": rates,
        "trial_status_counts": dict(Counter(r["status"] for r in results.values())),
        "old_joint_assessed_pairs": sum(assessed.values()),
        "old_joint_unassessed_pairs": sum(not v for v in assessed.values()),
        "paired_changes_on_assessed_only": dict(
            Counter(
                ("new_pass" if new[k] else "new_fail") + "/" + ("old_pass" if old_joint[k] else "old_fail")
                for k in new
                if assessed[k]
            )
        ),
        "scope": "All original pair denominators; failed/unstarted zero. Unassessed old-joint comparisons are operational zeros, not measured geometric failures. Assay inactives are not structural-negative truth. No fitting, activity preservation or synthesis claim.",
    }


@contextmanager
def cpu_only():
    """Deny model imports and network/process operations; fail on swallowed attempts."""
    active = [True]
    blocked = []
    if any(name.split(".")[0] in FORBIDDEN_IMPORTS for name in sys.modules):
        raise ValueError("Model runtime already imported")

    class Finder(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            if fullname.split(".")[0] in FORBIDDEN_IMPORTS:
                blocked.append("import:" + fullname)
                raise RuntimeError("Model import forbidden")

    def audit(event, arguments):
        if active[0] and event in {
            "socket.connect",
            "socket.bind",
            "socket.getaddrinfo",
            "socket.gethostbyname",
            "socket.gethostbyaddr",
            "socket.sendto",
            "socket.sendmsg",
            "subprocess.Popen",
            "os.system",
            "os.fork",
            "os.posix_spawn",
        }:
            blocked.append(event)
            raise RuntimeError("Network/process operation forbidden")

    finder = Finder()
    sys.meta_path.insert(0, finder)
    sys.addaudithook(audit)
    try:
        yield blocked
        if blocked:
            raise RuntimeError("Forbidden operation attempted: " + repr(blocked))
    finally:
        active[0] = False
        sys.meta_path.remove(finder)


def _safe_audit(value):
    try:
        base.json_bytes(value)
        return value
    except (ValueError, TypeError):
        return {"malformed_audit": True}


def _execute(plan, old_rows, output, backend, verify, wall_seconds, clock=time.monotonic):
    """Internal injected backend permits synthetic lifecycle tests, never CLI scoring fallback."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    (output / "trials").mkdir()
    (output / "references").mkdir()
    base.save(output / "plan.json", plan)
    pairs = {p["pair_id"]: p for p in plan["pairs"]}
    old = {(r["pair_id"], t["seed"]): t for r in old_rows for t in r["trials"]}
    old_pairs = {r["pair_id"]: r for r in old_rows}
    slots = plan["slots"]
    results = {
        s["index"]: {"slot": s, "status": "not_started", "reason": "planned_not_started"} for s in slots
    }
    started, cpu = clock(), time.process_time()
    deadline = started + wall_seconds
    interrupted = [None]
    handlers = {}
    run_status = "running"
    run_error = None
    reference = None
    group = None
    reference_error = None
    reference_artifact = None

    def checkpoint():
        ledger = [{k: v for k, v in row.items() if k not in {"info", "poses"}} for row in results.values()]
        base.save(
            output / "ledger.json",
            {
                "status": run_status,
                "planned_pairs": len(pairs),
                "planned_trials": len(slots),
                "trials": ledger,
            },
        )
        summary = summarize(plan["pairs"], slots, results, old_rows)
        summary.update(
            run_status=run_status,
            valid_for_interpretation=run_status in {"complete", "deadline", "interrupted"},
            resources={
                "wall_seconds": clock() - started,
                "process_cpu_seconds": time.process_time() - cpu,
                "process_lifetime_peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
                "gpu_seconds": None,
            },
            backend_costs=backend.costs(),
            interruption=interrupted[0],
            run_error=run_error,
        )
        base.save(output / "summary.json", summary)

    def receive(signum, frame):
        interrupted[0] = signum

    # Full denominator is durable before backend.start can embed anything.
    checkpoint()
    try:
        for signum in (signal.SIGINT, signal.SIGTERM):
            handlers[signum] = signal.signal(signum, receive)
        verify()
        backend.start()
        for slot in slots:
            if interrupted[0] is not None or clock() >= deadline:
                run_status = "interrupted" if interrupted[0] is not None else "deadline"
                break
            pair = pairs[slot["pair_id"]]
            current_group = (tuple(slot["query"]), slot["seed"])
            if current_group != group:
                if reference is not None:
                    backend.finish_reference(reference)
                reference, reference_error, reference_artifact, group = None, None, None, current_group
                try:
                    reference, record = backend.reference(pair, slot["seed"], deadline)
                    reference_path = output / "references" / f"{slot['index']:04d}.json"
                    base.save(reference_path, record)
                    reference_artifact = {
                        "path": str(reference_path),
                        "sha256": base.file_hash(reference_path),
                    }
                except Exception as exc:
                    reference_error = {"type": type(exc).__name__, "reason": str(exc)}
            row = {"slot": slot, "status": "running"}
            before, before_cpu = clock(), time.process_time()
            try:
                if interrupted[0] is not None or clock() >= deadline:
                    row.update(status="not_started", reason="deadline_or_signal_before_trial")
                elif reference_error is not None:
                    row.update(status="failed", reason="query_reference_failed", error=reference_error)
                else:
                    row.update(
                        backend.trial(
                            pair,
                            slot["seed"],
                            old[(pair["pair_id"], slot["seed"])],
                            old_pairs[pair["pair_id"]],
                            reference,
                            deadline,
                        )
                    )
                    if row["slot"] != slot or row["status"] not in {"success", "abstained", "failed"}:
                        raise ValueError("Invalid backend trial status")
                    base.json_bytes(row)
            except Exception as exc:
                row = {
                    "slot": slot,
                    "status": "failed",
                    "reason": str(exc),
                    "error_type": type(exc).__name__,
                    "audit": _safe_audit(getattr(exc, "audit", None)),
                }
            row["reference_artifact"] = reference_artifact
            row["costs"] = {
                "wall_seconds": clock() - before,
                "process_cpu_seconds": time.process_time() - before_cpu,
            }
            if (interrupted[0] is not None or clock() >= deadline) and row["status"] == "success":
                row.update(status="failed", reason="completed_after_deadline_or_signal")
            trial_path = output / "trials" / f"{slot['index']:04d}.json"
            base.save(trial_path, row)
            row["artifact"] = {"path": str(trial_path), "sha256": base.file_hash(trial_path)}
            results[slot["index"]] = row
            checkpoint()
        else:
            run_status = (
                "interrupted"
                if interrupted[0] is not None
                else ("deadline" if clock() >= deadline else "complete")
            )
        if reference is not None:
            backend.finish_reference(reference)
        verify()
    except BaseException as exc:
        run_status = "failed"
        run_error = {"type": type(exc).__name__, "reason": str(exc)}
        raise
    finally:
        try:
            backend.close()
            verify()
        except BaseException as exc:
            run_status = "integrity_or_cleanup_failed"
            run_error = {"type": type(exc).__name__, "reason": str(exc)}
            raise
        finally:
            checkpoint()
            for signum, previous in handlers.items():
                signal.signal(signum, previous)
    return results


def reference_record(value, settings, hit):
    return {
        "hit": hit,
        "settings": asdict(settings),
        "force_field": value.force_field,
        "diagnostics": value.failures,
        "energies": {str(k): v for k, v in value.energies.items()},
        "relative_energies": {str(k): v for k, v in value.relative_energies.items()},
        "poses": {str(cid): base.pose_record(value.mol, cid) for cid in value.energies},
    }


class ScientificBackend:
    """Single-process measured wrappers; no persistent/cross-run result cache."""

    def __init__(self):
        self.embeddings, self.requests = [], []
        self.started = False
        self.role = "unassigned"
        self.final_cache = None

    def start(self):
        from chembounce_synthesis import conformers, selective_scoring
        from chembounce_synthesis.similarity_constraints import BoundSimilarity

        self.kernel, self.scorer, self.policy_class = conformers, selective_scoring, BoundSimilarity
        self.original_embed = conformers.AllChem.EmbedMultipleConfs
        self.original_request = selective_scoring.ensemble
        conformers.ensemble.cache_clear()
        self.started = True
        conformers.AllChem.EmbedMultipleConfs = self._embed
        selective_scoring.ensemble = self._ensemble

    def _embed(self, *args, **kwargs):
        requested = kwargs.get("numConfs")
        role_limit = 288 if self.role == "reference_prewarm" else 4728
        role_requested = sum(e["requested_conformers"] for e in self.embeddings if e["role"] == self.role)
        if role_requested + (requested or 0) > role_limit:
            raise ValueError("Embedding role resource bound violated")
        params = kwargs.get("params")
        if (
            len(self.embeddings) >= SPEC["embedding_bounds"]["calls"]
            or requested != 8
            or params.numThreads != 1
        ):
            raise ValueError("Embedding resource bound violated")
        event = {
            "role": self.role,
            "requested_conformers": requested,
            "returned_conformers": None,
            "status": "started",
        }
        self.embeddings.append(event)
        start, cpu = time.monotonic(), time.process_time()
        try:
            ids = self.original_embed(*args, **kwargs)
            event.update(returned_conformers=len(ids), status="returned")
            return ids
        except BaseException as exc:
            event.update(status="failed", error_type=type(exc).__name__, reason=str(exc))
            raise
        finally:
            event.update(wall_seconds=time.monotonic() - start, process_cpu_seconds=time.process_time() - cpu)

    def _ensemble(self, smiles, settings):
        before = self.kernel.ensemble.cache_info()
        event = {"smiles": smiles, "settings": asdict(settings), "context": self.role, "status": "started"}
        self.requests.append(event)
        try:
            result = self.original_request(smiles, settings)
            event.update(status="success", retained_conformers=len(result.energies))
            return result
        except BaseException as exc:
            event.update(status="failed", error_type=type(exc).__name__, reason=str(exc))
            raise
        finally:
            after = self.kernel.ensemble.cache_info()
            event.update(cache_hits=after.hits - before.hits, cache_misses=after.misses - before.misses)

    def reference(self, pair, seed, deadline):
        if time.monotonic() >= deadline:
            raise TimeoutError("Reference deadline exhausted")
        settings = self.kernel.ConformerSettings(seed=seed, **base.SETTINGS)
        policy = self.policy_class(pair["query_smiles"], FLOORS, polar_recall=0.75)
        self.role = "reference_prewarm"
        ensemble = self._ensemble(pair["query_smiles"], settings)
        record = reference_record(ensemble, settings, pair["query_smiles"])
        record.update(policy=policy.policy, policy_sha256=policy.sha256)
        ref = {
            "ensemble": ensemble,
            "settings": settings,
            "policy": policy,
            "record": record,
            "sha256": base.json_hash(record),
        }
        return ref, record

    def finish_reference(self, ref):
        record = reference_record(ref["ensemble"], ref["settings"], ref["policy"].hit)
        record.update(policy=ref["policy"].policy, policy_sha256=ref["policy"].sha256)
        if base.json_hash(record) != ref["sha256"]:
            raise ValueError("Cached query reference mutated")

    def old_comparison(self, pair, old_trial, old_pair, ref):
        from chembounce_synthesis.qualification import validate_result

        scope = "saved_aligned_candidate_with_regenerated_reference_not_original_saved_reference_coordinates"
        try:
            if old_trial["status"] != "success":
                raise ValueError("Old trial did not succeed")
            info = old_trial["info"]
            poses = [base.restore_pose(old_trial["poses"][key]) for key in ("hit", "design")]
            validate_result(
                (info, *poses),
                pair["control_smiles"],
                ref["settings"],
                hit=pair["query_smiles"],
                design=pair["query_smiles"],
            )
            ensemble = ref["ensemble"]
            rid, cid = info["hit"]["reference_conformer_id"], info["candidate_conformer_id"]
            if (
                rid not in ensemble.relative_energies
                or info["hit"]["reference_relative_energy"] != ensemble.relative_energies[rid]
                or info["ensemble_diagnostics"]["hit"] != ensemble.failures
                or info["ensemble_sizes"]["hit"] != len(ensemble.energies)
            ):
                raise ValueError("Old regenerated reference ID/energy/diagnostics binding unavailable")
            witness = ref["policy"].assess(ensemble.mol, poses[0], rid, cid)
            metrics = {
                "hit_tanimoto": old_pair["hit_tanimoto"],
                "hit_shape": info["hit"]["shape"],
                "hit_pharmacophore": info["hit"]["pharmacophore"],
            }
            return {
                "status": "assessed",
                "passed": witness["passed"] and all(metrics[k] >= v for k, v in FLOORS.items()),
                "witness": witness,
                "saved_metrics": metrics,
                "reference_record_sha256": ref["sha256"],
                "scope": scope,
            }
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            return {"status": "unavailable", "reason": str(exc), "scope": scope}

    def trial(self, pair, seed, old_trial, old_pair, ref, deadline):
        old = self.old_comparison(pair, old_trial, old_pair, ref)
        self.role = "candidate_trial"
        try:
            result = self.scorer.score_candidate(
                pair["control_smiles"],
                pair["query_smiles"],
                pair["query_smiles"],
                ref["settings"],
                similarity=ref["policy"],
                deadline_monotonic=deadline,
            )
            self.scorer.validate_result(
                result,
                pair["control_smiles"],
                ref["settings"],
                hit=pair["query_smiles"],
                design=pair["query_smiles"],
                similarity=ref["policy"],
            )
            info, hitpose, designpose = result
            cid, rid = info["candidate_conformer_id"], info["hit"]["reference_conformer_id"]
            restored = ref["policy"].validate_witness(info["similarity_constraints"]["witness"], hitpose, cid)
            if base.pose_record(restored, rid) != base.pose_record(ref["ensemble"].mol, rid):
                # Molblock headers/atom serialization need not be identical: the
                # actual coordinates and canonical atom order are the witness.
                import numpy as np
                from rdkit import Chem

                if Chem.MolToSmiles(restored, canonical=False) != Chem.MolToSmiles(
                    ref["ensemble"].mol, canonical=False
                ) or not np.array_equal(
                    restored.GetConformer(rid).GetPositions(),
                    ref["ensemble"].mol.GetConformer(rid).GetPositions(),
                ):
                    raise ValueError("Selected reference differs from prewarmed reference")
            return {
                "status": "success",
                "info": info,
                "poses": {"hit": base.pose_record(hitpose, cid), "design": base.pose_record(designpose, cid)},
                "old_comparison": old,
                "reference_record_sha256": ref["sha256"],
            }
        except self.scorer.SimilarityAbstention as exc:
            return {"status": "abstained", "reason": exc.reason, "audit": exc.audit, "old_comparison": old}
        except Exception as exc:
            return {
                "status": "failed",
                "reason": str(exc),
                "error_type": type(exc).__name__,
                "old_comparison": old,
            }

    def costs(self):
        return {
            "actual_embedding_calls": len(self.embeddings),
            "actual_embedding_events": self.embeddings,
            "ensemble_requests": self.requests,
            "cache": self.final_cache,
            "scope": "Actual intercepted EmbedMultipleConfs calls; cache misses can fail before embedding. Returned conformers precede optimization/filtering. No assumed unique-reference count.",
        }

    def close(self):
        if self.started:
            self.kernel.AllChem.EmbedMultipleConfs = self.original_embed
            self.scorer.ensemble = self.original_request
            self.final_cache = self.kernel.ensemble.cache_info()._asdict()
            self.kernel.ensemble.cache_clear()
            self.started = False


def run(plan_path, plan_sha256, output):
    if base.file_hash(plan_path) != plan_sha256:
        raise ValueError("Frozen plan digest mismatch")
    plan = read(plan_path)
    declaration = plan["declaration"]
    if (
        plan.get("schema") != "qualified-selective-development-plan-v2"
        or base.json_hash(declaration) != plan["declaration_canonical_sha256"]
        or base.file_hash(plan["declaration_path"]) != plan["declaration_sha256"]
        or read(plan["declaration_path"]) != declaration
        or plan["pairs"] != read(declaration["inputs"]["pairs"]["path"])
        or plan["slots"] != ordered_slots(plan["pairs"])
        or base.counts(plan["pairs"]) != COUNTS
    ):
        raise ValueError("Plan does not match fixed declaration/inventory")

    checks = []

    def verify():
        if blocked:
            raise ValueError("Forbidden operation attempted: " + repr(blocked))
        if base.file_hash(plan_path) != plan_sha256:
            raise ValueError("Plan changed during execution")
        if base.file_hash(plan["declaration_path"]) != plan["declaration_sha256"]:
            raise ValueError("Declaration changed during execution")
        verify_declaration(declaration)
        checks.append({"check": len(checks), "status": "all_source_runtime_data_and_plan_pins_match"})
        if Path(output).exists():
            base.save(
                Path(output) / "integrity.json",
                {
                    "checks": checks,
                    "expected_bindings": declaration,
                    "plan_sha256": plan_sha256,
                    "forbidden_attempts": blocked,
                },
            )

    with cpu_only() as blocked:
        verify()
        old = read(declaration["inputs"]["old_scores"]["path"])
        return _execute(plan, old, output, ScientificBackend(), verify, SPEC["wall_seconds"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("spec")
    for command, item in (("prepare", "declaration"), ("run", "plan")):
        part = sub.add_parser(command)
        part.add_argument("--" + item, required=True)
        part.add_argument("--sha256", required=True)
        part.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "spec":
        print(
            json.dumps(
                {
                    "spec": SPEC,
                    "input_sha256": INPUTS,
                    "source_sha256": sources(),
                    "runtime": base.runtime(),
                    "declaration_schema": "qualified-selective-development-declaration-v2",
                },
                indent=2,
            )
        )
    elif args.command == "prepare":
        prepare(args.declaration, args.sha256, args.output)
    else:
        run(args.plan, args.sha256, args.output)


if __name__ == "__main__":
    main()
