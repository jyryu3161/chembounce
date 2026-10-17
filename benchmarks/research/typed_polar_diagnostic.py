"""Experimental fixed-population pose diagnostic; never a recommendation profile.

The baseline calls the production scorer unchanged. Only failed baseline calls
receive bounded typed-polar rigid alternatives, on the same generated ensembles.
No model, catalog, route search, or synthesis prediction is performed here.
"""

import argparse
import copy
from dataclasses import asdict
import hashlib
import itertools
import json
import math
import sys
from pathlib import Path
import time
from types import FunctionType

import numpy as np
import scipy
from rdkit import rdBase

from chembounce_synthesis import conformers, selective_scoring as scoring
from chembounce_synthesis.common import canonical
from chembounce_synthesis.landmark_proposals import LIMITS, ProposalSearch, declaration
from chembounce_synthesis.similarity_constraints import BoundSimilarity, pose_binding

VERSION = "experimental-typed-polar-rigid-diagnostic-v1"
SEEDS = (42, 43, 44)
FLOORS = {"hit_tanimoto": 0.35, "hit_shape": 0.7, "hit_pharmacophore": 0.5}
PRODUCTION_SOURCES = tuple(
    "chembounce_synthesis/" + name + ".py"
    for name in (
        "__init__",
        "common",
        "conformers",
        "similarity_constraints",
        "selective_scoring",
        "landmark_proposals",
        "qualification",
        "evidence",
        "refinement",
        "profiles",
    )
)


class DiagnosticError(Exception):
    """Parity or input mutation is unknown/error, never a rescued candidate."""


class PolarSearch(ProposalSearch):
    """Reuse existing proper transforms, fixed-pose scoring and bounded accounting.

    Only correspondence enumeration differs from the explicit-landmark helper.
    Inherited proposals are stripped of landmark provenance; their full rigid
    evidence is retained in the separate experimental audit instead.
    """

    def __init__(self, polar_recall=0.75, *, clock=time.monotonic):
        super().__init__(declaration(), clock=clock)
        self.polar_recall = polar_recall
        self.policy = {"version": VERSION, "limits": dict(LIMITS), "polar_recall": polar_recall}
        self.audit["policy"] = copy.deepcopy(self.policy)
        self.audit["rigid_pose_evidence"] = []
        self.audit["accepted_alternatives"] = []

    def _maps(self, ref, cand, requests, row):
        polar = [i for i, f in enumerate(ref) if f["family"] in {"Acceptor", "Donor"}]
        row["polar_reference_features"] = len(polar)
        if not polar:
            row["mapping_reason"] = "no_polar_reference"
            return []
        needed = math.ceil(self.polar_recall * len(polar))
        choices = {
            i: sorted(
                (j for j, f in enumerate(cand) if f["family"] == ref[i]["family"]),
                key=lambda j: (sum((a - b) ** 2 for a, b in zip(ref[i]["position"], cand[j]["position"])), j),
            )
            for i in polar
        }
        capacity = sum(
            min(sum(ref[i]["family"] == kind for i in polar), sum(f["family"] == kind for f in cand))
            for kind in ("Acceptor", "Donor")
        )
        if capacity < needed:
            row["mapping_reason"] = "insufficient_typed_polar_capacity"
            return []
        maps = []

        def walk(required, pairs, used):
            if not self.available():
                return
            if len(maps) >= LIMITS["complete_maps_per_pair"]:
                row["mapping_truncated"] = self.audit["truncated"] = True
                return
            if len(pairs) == len(required):
                maps.append(pairs)
                return
            i = required[len(pairs)]
            for j in choices[i]:
                if row["mapping_edge_visits"] >= LIMITS["mapping_edge_visits_per_pair"]:
                    row["mapping_truncated"] = self.audit["truncated"] = True
                    return
                row["mapping_edge_visits"] += 1
                if j not in used:
                    walk(required, pairs + [(i, j)], used | {j})

        # Subsets allow a missing polar type only when the unchanged recall
        # floor permits it. No zero-weight/dummy correspondences propose poses.
        for required in itertools.combinations(polar, needed):
            if not self.available():
                break
            if len(maps) >= LIMITS["complete_maps_per_pair"]:
                row["mapping_truncated"] = self.audit["truncated"] = True
                break
            if row["helper_edge_visits"] >= LIMITS["helper_edge_visits_per_pair"]:
                row["helper_truncated"] = self.audit["truncated"] = True
                break
            row["helper_edge_visits"] += 1
            walk(required, [], set())
        row["complete_maps"] = len(maps)
        return maps

    def _anchors(self, mapping, ref, cand, row):
        yield mapping  # Only typed polar anchors; no implicit required landmarks.

    def propose(self, *args):
        output = super().propose(*args)
        for info, pose, witness in output:
            evidence = info.pop("landmark_proposal")
            self.audit["rigid_pose_evidence"].append(
                {
                    "rigid_transform": evidence,
                    "candidate_pose": pose_binding(pose, evidence["baseline_candidate_pose"]["conformer_id"]),
                }
            )
        return output


def selected_alternative(info, audit):
    """A rescue must be the exact evaluated, accepted rigid alternative."""
    for index, accepted in enumerate(audit["accepted_alternatives"]):
        evidence = audit["rigid_pose_evidence"][accepted["rigid_evidence_index"]]
        if (
            info["selective_pose_evidence"]["hit_candidate"]
            == accepted["candidate_pose"]
            == evidence["candidate_pose"]
            and info["hit"] == accepted["hit_info"]
            and info["similarity_constraints"]["witness"] == accepted["witness"]
            and evidence["rigid_transform"]["reference_pose_sha256"]
            == accepted["witness"]["reference_pose"]["sha256"]
        ):
            return index
    raise DiagnosticError("Selected feasible pose is not an accepted rigid alternative")


def ensemble_state(ensembles):
    return {
        smiles: {
            "energies": dict(value.energies),
            "relative_energies": dict(value.relative_energies),
            "force_field": value.force_field,
            "poses": {str(cid): pose_binding(value.mol, cid)["sha256"] for cid in value.relative_energies},
        }
        for smiles, value in ensembles.items()
    }


def invoke(function, product, hit, design, settings, similarity):
    try:
        result = function(product, hit, design, settings, similarity=similarity)
        scoring.validate_result(result, product, settings, hit=hit, design=design, similarity=similarity)
        return {"status": "feasible", "info": result[0]}, result
    except scoring.SimilarityAbstention as exc:
        return {"status": "abstained", "reason": exc.reason, "audit": exc.audit}, None


def paired_condition(product, hit, design, seed):
    """One production baseline and, only if needed, one experimental paired call."""
    started, cpu = time.monotonic(), time.process_time()
    product, hit, design = map(canonical, (product, hit, design))
    settings = conformers.ConformerSettings(num_conformers=8, seed=seed, scoring_objective="original-hit-v1")
    similarity = BoundSimilarity(hit, FLOORS, polar_recall=0.75, landmarks=[])
    provider = conformers.ensemble
    if scoring.score_candidate.__globals__["ensemble"] is not provider:
        raise DiagnosticError("Production baseline uses a different ensemble provider")
    construction = time.monotonic()
    ensembles = {
        smiles: conformers.ensemble(smiles, settings) for smiles in dict.fromkeys((product, hit, design))
    }
    construction = time.monotonic() - construction
    before = ensemble_state(ensembles)
    baseline_start = time.monotonic()
    baseline, baseline_result = invoke(scoring.score_candidate, product, hit, design, settings, similarity)
    baseline_seconds = time.monotonic() - baseline_start
    cache_check_start = time.monotonic()
    if (
        scoring.score_candidate.__globals__["ensemble"] is not provider
        or any(provider(smiles, settings) is not value for smiles, value in ensembles.items())
        or ensemble_state(ensembles) != before
    ):
        raise DiagnosticError("Production baseline cache object/identity/coordinates changed")
    cache_check_seconds = time.monotonic() - cache_check_start
    row = {
        "version": VERSION,
        "product": product,
        "hit": hit,
        "design": design,
        "seed": seed,
        "settings": asdict(settings),
        "ensemble_bindings": before,
        "baseline": baseline,
        "cost": {
            "ensemble_construction_seconds": construction,
            "production_baseline_seconds": baseline_seconds,
            "baseline_cache_object_check_seconds": cache_check_seconds,
        },
    }
    if baseline_result is not None:
        row.update(status="baseline_feasible_preserved", prototype=copy.deepcopy(baseline), extra_search=None)
    elif baseline.get("reason") != "no_jointly_feasible_hit_pose" or any(
        pair["status"] != "rejected" for pair in baseline.get("audit", {}).get("pairs", [])
    ):
        row.update(status="baseline_unavailable", prototype={"status": "not_attempted"}, extra_search=None)
    else:
        search = PolarSearch()
        expected = {
            (p["candidate_conformer_id"], p["reference_conformer_id"]): p for p in baseline["audit"]["pairs"]
        }
        pair_count = len(ensembles[product].relative_energies) * len(ensembles[hit].relative_energies)
        calls = {
            "hit_overlay_replays": 0,
            "design_overlays": 0,
            "overlay_replay_seconds": 0.0,
            "callback_seconds": 0.0,
        }

        def overlay(reference, candidate, cid, options):
            is_hit_pair = calls["hit_overlay_replays"] < pair_count
            calls["hit_overlay_replays" if is_hit_pair else "design_overlays"] += 1
            start = time.monotonic()
            try:
                info, aligned = conformers._overlay(reference, candidate, cid, options)
            except (ValueError, RuntimeError) as exc:
                if is_hit_pair:
                    raise DiagnosticError("Baseline overlay replay failed") from exc
                raise
            finally:
                calls["overlay_replay_seconds"] += time.monotonic() - start
            if not is_hit_pair:
                return info, aligned
            callback_start = time.monotonic()
            spent_before = search.spent
            chosen = info, aligned
            chosen_witness = None
            accepted = None
            try:
                rid = info["reference_conformer_id"]
                old = expected.get((cid, rid), {})
                if old.get("status") != "rejected" or any(
                    not math.isclose(info[k], old.get(k, -1), abs_tol=1e-12, rel_tol=1e-12)
                    for k in ("shape", "pharmacophore")
                ):
                    raise DiagnosticError("Baseline overlay replay differs from production audit")
                witness = similarity.assess(reference.mol, aligned, rid, cid)
                current_recall, prior_recall = witness["polar"]["recall"], old["polar_recall"]
                if ((current_recall is None) != (prior_recall is None)) or (
                    current_recall is not None
                    and not math.isclose(current_recall, prior_recall, abs_tol=1e-12, rel_tol=1e-12)
                ):
                    raise DiagnosticError("Baseline polar witness replay mismatch")
                # Scalar closeness cannot excuse a changed acceptance decision.
                reasons = list(witness["reasons"])
                for name, field in (("hit_shape", "shape"), ("hit_pharmacophore", "pharmacophore")):
                    if info[field] < FLOORS[name]:
                        reasons.append(name + "_below_minimum")
                if reasons != old["reasons"] or not reasons:
                    raise DiagnosticError("Baseline replay joint/floor decisions changed")
                if search.available():
                    if "polar_reference_recall_below_minimum" in witness["reasons"]:
                        options = search.propose(
                            reference.mol, aligned, rid, cid, info, witness, similarity, settings
                        )
                        if search.audit["pairs"][-1]["status"] in {
                            "fixed_score_parity_unavailable",
                            "unavailable",
                        }:
                            raise DiagnosticError("Fixed-pose proposal scoring/parity unavailable")
                        feasible = [
                            (extra, pose, decision)
                            for extra, pose, decision in options
                            if decision["passed"]
                            and extra["shape"] >= FLOORS["hit_shape"]
                            and extra["pharmacophore"] >= FLOORS["hit_pharmacophore"]
                        ]
                        if feasible:
                            extra, pose, chosen_witness = max(
                                feasible, key=lambda pair: pair[0]["reference_adjusted_joint"]
                            )
                            chosen = extra, pose
                if chosen_witness is not None:
                    binding = pose_binding(chosen[1], cid)
                    evidence_index = next(
                        (
                            i
                            for i, evidence in enumerate(search.audit["rigid_pose_evidence"])
                            if evidence["candidate_pose"] == binding
                            and evidence["rigid_transform"]["reference_pose_sha256"]
                            == chosen_witness["reference_pose"]["sha256"]
                        ),
                        None,
                    )
                    if evidence_index is None:
                        raise DiagnosticError("Accepted alternative lacks rigid proposal evidence")
                    accepted = {
                        "candidate_pose": binding,
                        "hit_info": copy.deepcopy(chosen[0]),
                        "witness": copy.deepcopy(chosen_witness),
                        "rigid_evidence_index": evidence_index,
                    }
            finally:
                elapsed = time.monotonic() - callback_start
                calls["callback_seconds"] += elapsed
                search.charge(max(0.0, elapsed - (search.spent - spent_before)))
            # Charge callback/witness overhead before accepting an alternative.
            if not search.available():
                return info, aligned
            if accepted is not None:
                search.audit["accepted_alternatives"].append(accepted)
            return chosen

        def fixed_ensemble(smiles, requested):
            if requested != settings or smiles not in ensembles:
                raise DiagnosticError("Unexpected ensemble request")
            return ensembles[smiles]

        # Execute exactly production bytecode with a private globals dictionary.
        # No production module attribute/global or cached ensemble is replaced.
        original = scoring.score_candidate
        namespace = dict(original.__globals__, _overlay=overlay, ensemble=fixed_ensemble)
        experimental = FunctionType(original.__code__, namespace, original.__name__, original.__defaults__)
        experimental.__kwdefaults__ = original.__kwdefaults__.copy()
        prototype_start = time.monotonic()
        try:
            prototype, _ = invoke(experimental, product, hit, design, settings, similarity)
            if prototype["status"] == "feasible":
                row["selected_alternative_index"] = selected_alternative(prototype["info"], search.audit)
            row.update(
                status="rescued" if prototype["status"] == "feasible" else "not_rescued", prototype=prototype
            )
            if prototype["status"] != "feasible" and search.audit["truncated"]:
                row["status"] = "unresolved_allowance_exhausted"
        except Exception as exc:
            row.update(
                status="diagnostic_error",
                prototype={"status": "unknown", "error": f"{type(exc).__name__}: {exc}"},
            )
        row["cost"].update(prototype_replay_seconds=time.monotonic() - prototype_start, **calls)
        row["extra_search"] = search.audit
    if ensemble_state(ensembles) != before:
        row.update(
            status="diagnostic_error", prototype={"status": "unknown", "error": "Cached ensemble mutated"}
        )
    row["cost"].update(total_seconds=time.monotonic() - started, cpu_seconds=time.process_time() - cpu)
    return row


def write(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    raw = args.plan.read_bytes()
    if hashlib.sha256(raw).hexdigest() != args.plan_sha256:
        raise ValueError("Plan byte mismatch")
    plan = json.loads(raw)
    if (
        plan["version"] != VERSION
        or plan["seeds"] != list(SEEDS)
        or len(plan["population"]) != 24
        or plan["source_sha256"] != hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        or plan["limits"] != LIMITS
        or plan["minimum_scores"] != FLOORS
        or plan["polar_recall"] != 0.75
        or [(p["attempt"], p["index"]) for p in plan["population"]]
        != [(attempt, index) for attempt in ("attempt-0000", "attempt-0001") for index in range(12)]
    ):
        raise ValueError("Unsupported prospective diagnostic plan")
    actual_runtime = {
        "python": sys.version.split()[0],
        "rdkit": rdBase.rdkitVersion,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    if plan["runtime_versions"] != actual_runtime:
        raise ValueError("Diagnostic runtime version mismatch")
    if (
        plan["python_executable_sha256"]
        != hashlib.sha256(Path(sys.executable).resolve().read_bytes()).hexdigest()
    ):
        raise ValueError("Diagnostic Python executable mismatch")
    root = Path(__file__).resolve().parents[2]
    for name, expected in plan["production_source_sha256"].items():
        if (
            name not in PRODUCTION_SOURCES
            or hashlib.sha256((root / name).read_bytes()).hexdigest() != expected
        ):
            raise ValueError("Diagnostic production helper source mismatch")
    if set(plan["production_source_sha256"]) != set(PRODUCTION_SOURCES):
        raise ValueError("Incomplete diagnostic production helper pins")
    args.output.mkdir(parents=True, exist_ok=False)
    rows = [
        {"attempt": p["attempt"], "population_index": p["index"], "seed": seed, "status": "not_started"}
        for p in plan["population"]
        for seed in SEEDS
    ]
    report = {
        "version": VERSION,
        "plan_sha256": args.plan_sha256,
        "experimental": True,
        "profile_matching": "not_applicable_experimental",
        "source_sha256": plan["source_sha256"],
        "rdkit": rdBase.rdkitVersion,
        "numpy": np.__version__,
        "runtime_versions": actual_runtime,
        "production_source_sha256": plan["production_source_sha256"],
        "python_executable_sha256": plan["python_executable_sha256"],
        "runtime_pin_scope": "Python executable bytes and package versions, not a complete native-library inventory",
        "rows": rows,
        "scope": "Paired pose-only development diagnostic; no profile matching, synthesis or efficacy claim",
    }
    write(args.output / "summary.json", report)
    for index, (candidate, seed) in enumerate(itertools.product(plan["population"], SEEDS)):
        rows[index]["status"] = "running"
        write(args.output / "summary.json", report)
        outer_started, outer_cpu = time.monotonic(), time.process_time()
        try:
            condition = paired_condition(candidate["product"], candidate["hit"], candidate["design"], seed)
        except Exception as exc:
            condition = {"status": "diagnostic_error", "error": f"{type(exc).__name__}: {exc}"}
        outer_cost = {
            "wall_seconds": time.monotonic() - outer_started,
            "cpu_seconds": time.process_time() - outer_cpu,
        }
        condition["outer_cost_including_errors"] = outer_cost
        name = f"condition-{index:03d}.json"
        write(args.output / name, condition)
        rows[index].update(
            status=condition["status"],
            outer_cost=outer_cost,
            evidence=name,
            evidence_sha256=hashlib.sha256((args.output / name).read_bytes()).hexdigest(),
        )
        for arm in ("baseline", "prototype"):
            observed = condition.get(arm, {})
            metrics = observed.get("info", {})
            rows[index][arm] = {"status": observed.get("status", "unknown")}
            if metrics:
                rows[index][arm].update(
                    shape=metrics["hit"]["shape"],
                    pharmacophore=metrics["hit"]["pharmacophore"],
                    polar_recall=metrics["similarity_constraints"]["witness"]["polar"]["recall"],
                )
        measured = [r["outer_cost"] for r in rows if "outer_cost" in r]
        report["cost_coverage"] = {
            "measured_slots": len(measured),
            "planned_slots": len(rows),
            "unmeasured_slots": len(rows) - len(measured),
            "wall_seconds": sum(c["wall_seconds"] for c in measured),
            "cpu_seconds": sum(c["cpu_seconds"] for c in measured),
            "scope": "All returned/error condition calls; excludes JSON persistence. Not deployable per-input cost.",
        }
        write(args.output / "summary.json", report)
    report["paired_products"] = []
    for index, candidate in enumerate(plan["population"]):
        group = rows[3 * index : 3 * index + 3]
        item = {
            "attempt": candidate["attempt"],
            "population_index": candidate["index"],
            "product": candidate["product"],
            "design": candidate["design"],
        }
        for arm in ("baseline", "prototype"):
            statuses = [r[arm]["status"] for r in group]
            item[arm] = {
                "seed_statuses": statuses,
                "all_three_feasible": statuses == ["feasible"] * 3,
                "bounded_calls_returned": sum(x in {"feasible", "abstained"} for x in statuses),
                "assessment_complete": all(x in {"feasible", "abstained"} for x in statuses)
                and (
                    arm == "baseline"
                    or not any(r["status"] == "unresolved_allowance_exhausted" for r in group)
                ),
            }
        report["paired_products"].append(item)
    write(args.output / "summary.json", report)
    return 1 if any(row["status"] == "diagnostic_error" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
