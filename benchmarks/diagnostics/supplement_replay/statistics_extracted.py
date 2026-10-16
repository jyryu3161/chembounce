from collections import Counter, defaultdict

import statistics

INFRA = "infrastructure_interruption"

RETENTION_FIELDS = (
    "query_scaffold_subgraph_retained",
    "query_generic_ring_skeleton_subgraph_retained",
    "core_retaining_scaffold_expansion",
    "query_scaffold_not_contained",
)

BOUNDED = (
    "attempt_with_hop",
    "attempt_with_verified_route",
    "attempt_with_valid_output",
    "attempt_with_common_3d",
    "common_3d_output_slot_coverage",
    "route_supplied_output_slot_coverage",
    *(f"attempt_with_{name}" for name in RETENTION_FIELDS),
    *(f"{name}_assessed_output_slot_coverage" for name in RETENTION_FIELDS),
)

def key(row):
    return row["split"], row["method_id"], row["query_id"], row["seed"]

def occurrence_key(row):
    return (*key(row), row["rank"], row["smiles"])

def mean(values):
    observed = [v for v in values if v is not None]
    return statistics.mean(observed) if observed else None

def attempt_metrics(attempts, occurrences, scores):
    grouped = defaultdict(list)
    for o in occurrences:
        grouped[key(o)].append(o)
    result = []
    for attempt in attempts:
        rows = grouped[key(attempt)]
        successful = [o for o in rows if scores[occurrence_key(o)]["all_3d_seeds_succeeded"]]
        route_supplied = sum(o["route_supplied"] for o in rows)
        metrics = dict(attempt["evaluation"]["metrics"])
        metrics.update(
            attempt_with_hop=int(metrics["hop_families"] > 0),
            attempt_with_valid_output=int(bool(rows)),
            attempt_with_verified_route=int(metrics["verified_rule_routes"] > 0),
            attempt_with_common_3d=int(bool(successful)),
            common_3d_successful_candidates=len(successful),
            common_3d_failed_candidates=len(rows) - len(successful),
            common_3d_output_slot_coverage=len(successful) / 3,
            route_supplied_occurrences=route_supplied,
            no_route_supplied_occurrences=len(rows) - route_supplied,
            route_supplied_output_slot_coverage=route_supplied / 3,
            route_replay_verified_occurrences=sum(
                o["candidate"]["route_evidence"]["level"] == "catalog_and_rule_verified" for o in rows
            ),
            route_rejected_occurrences=sum(
                o["candidate"]["route_evidence"]["level"] == "route_rejected" for o in rows
            ),
            wall_seconds=attempt["resources"].get("wall_seconds"),
            hit_bit_tanimoto_mean_conditional=mean(o["candidate"]["hit_bit_tanimoto"] for o in rows),
            hit_count_tanimoto_mean_conditional=mean(o["candidate"]["hit_count_tanimoto"] for o in rows),
        )
        for name in ("hit_shape", "hit_pharmacophore", "score"):
            metrics[f"common_{name}_mean_success_conditional"] = mean(
                scores[occurrence_key(o)]["scores"][name] for o in successful
            )
            metrics[f"common_{name}_range_mean_success_conditional"] = mean(
                scores[occurrence_key(o)]["scores"][name + "_stability"]["range"] for o in successful
            )
        for name in RETENTION_FIELDS:
            qualifying = [o for o in rows if o["scaffold_retention"][name] is True]
            assessed = sum(o["scaffold_retention"][name] is not None for o in rows)
            metrics.update(
                {
                    name + "_occurrences": len(qualifying),
                    name + "_unique_families": len({o["candidate"]["family"] for o in qualifying}),
                    name + "_assessed_occurrences": assessed,
                    name + "_unassessed_occurrences": len(rows) - assessed,
                    name + "_assessable_fraction_conditional": assessed / len(rows) if rows else None,
                    name + "_assessed_output_slot_coverage": assessed / 3,
                    "attempt_with_" + name: int(bool(qualifying)),
                }
            )
        declarations = [o["candidate"].get("declared_external_synthesis_evidence") or {} for o in rows]
        steps = [d.get("steps", []) for d in declarations]
        metrics["declared_forward_assessed_occurrences"] = sum(
            bool(s)
            and all(t.get("forward_prediction", {}).get("status") in {"matched", "mismatched"} for t in s)
            for s in steps
        )
        metrics["declared_forward_all_steps_matched_occurrences"] = sum(
            bool(s) and all(t.get("forward_prediction", {}).get("status") == "matched" for t in s)
            for s in steps
        )
        result.append(
            {
                **{k: attempt[k] for k in ("split", "method_id", "query_id", "target_id", "seed", "status")},
                "metrics": metrics,
            }
        )
    return result

def hierarchical(values, slots, interval):
    """Missing conditional scores remain missing; targets always receive equal weight."""
    queries = sorted({s["query_id"] for s in slots})
    targets = sorted({s["target_id"] for s in slots})
    by_query = {
        q: mean(values.get((s["query_id"], s["seed"])) for s in slots if s["query_id"] == q) for q in queries
    }
    per_target = {
        t: mean(
            by_query[q] for q in queries if any(s["query_id"] == q and s["target_id"] == t for s in slots)
        )
        for t in targets
    }
    observed = {
        q: sum(values.get((s["query_id"], s["seed"])) is not None for s in slots if s["query_id"] == q)
        for q in queries
    }
    estimate = interval([v for v in per_target.values() if v is not None], seed=20261016, resamples=10000)
    return {
        "per_query": by_query,
        "per_target": per_target,
        "target_macro": estimate,
        "availability": {
            "planned_slots": len(slots),
            "observed_slots": sum(observed.values()),
            "missing_slots": len(slots) - sum(observed.values()),
            "observed_seeds_per_query": observed,
            "planned_queries": len(queries),
            "observed_queries": sum(v is not None for v in by_query.values()),
            "planned_targets": len(targets),
            "observed_targets": sum(v is not None for v in per_target.values()),
        },
    }

def comparisons(rows, grid, interval):
    output = {}
    for split in preparation.SPLITS:
        planned = [g for g in grid if g["split"] == split and g["method_id"] == preparation.METHODS[0]]
        lookup = {
            m: {(r["query_id"], r["seed"]): r for r in rows if r["split"] == split and r["method_id"] == m}
            for m in preparation.METHODS
        }
        required = {(s["query_id"], s["seed"]) for s in planned}
        if any(set(d) != required for d in lookup.values()):
            raise ValueError("Analysis is missing matched method/query/seed slots")
        names = sorted(next(iter(lookup[preparation.METHODS[0]].values()))["metrics"])
        methods, paired = {}, {}
        for method, local in lookup.items():
            if any(set(r["metrics"]) != set(names) for r in local.values()):
                raise ValueError("Metric schema differs between attempts")
            summary = {
                "attempts": len(local),
                "status_counts": dict(Counter(r["status"] for r in local.values())),
                "metrics": {},
                "complete_case_excluding_infrastructure_only": {},
                "infrastructure_counterfactual_bounds": {},
            }
            for metric in names:
                values = {k: r["metrics"][metric] for k, r in local.items()}
                summary["metrics"][metric] = hierarchical(values, planned, interval)
                complete = {k: r["metrics"][metric] for k, r in local.items() if r["status"] != INFRA}
                summary["complete_case_excluding_infrastructure_only"][metric] = hierarchical(
                    complete, planned, interval
                )
                if metric in BOUNDED:
                    low = {k: 0.0 if r["status"] == INFRA else r["metrics"][metric] for k, r in local.items()}
                    high = {
                        k: 1.0 if r["status"] == INFRA else r["metrics"][metric] for k, r in local.items()
                    }
                    summary["infrastructure_counterfactual_bounds"][metric] = [
                        hierarchical(v, planned, interval)["target_macro"]["mean"] for v in (low, high)
                    ]
            methods[method] = summary
        primary = lookup[preparation.METHODS[0]]
        for method in preparation.METHODS[1:]:
            other = lookup[method]
            paired[method] = {
                "metrics": {},
                "complete_case_excluding_either_infrastructure": {},
                "infrastructure_counterfactual_bounds": {},
            }
            for metric in names:
                paired_values, complete, low, high = {}, {}, {}, {}
                for k in required:
                    a, b = primary[k], other[k]
                    av, bv = a["metrics"][metric], b["metrics"][metric]
                    if av is not None and bv is not None:
                        paired_values[k] = av - bv
                        if a["status"] != INFRA and b["status"] != INFRA:
                            complete[k] = av - bv
                    if metric in BOUNDED:
                        alo, ahi = (0.0, 1.0) if a["status"] == INFRA else (av, av)
                        blo, bhi = (0.0, 1.0) if b["status"] == INFRA else (bv, bv)
                        low[k], high[k] = alo - bhi, ahi - blo
                paired[method]["metrics"][metric] = hierarchical(paired_values, planned, interval)
                paired[method]["complete_case_excluding_either_infrastructure"][metric] = hierarchical(
                    complete, planned, interval
                )
                if metric in BOUNDED:
                    paired[method]["infrastructure_counterfactual_bounds"][metric] = [
                        hierarchical(v, planned, interval)["target_macro"]["mean"] for v in (low, high)
                    ]
        output[split] = {
            "methods": methods,
            "paired_pipeline_current_minus_comparator": paired,
            "counterfactual_bounds_are_confidence_intervals": False,
            "complete_case_is_selection_biased": True,
        }
    return output

import math

import numpy as np

def finite_scores(values):
    if not isinstance(values, list) or not values:
        return None
    if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in values):
        return None
    return [float(v) for v in values]

def score_stability(samples):
    values = finite_scores(samples)
    if values is None:
        return None
    return {
        "n": len(values),
        "mean": statistics.mean(values),
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
        "range": max(values) - min(values),
        "sample_sd": statistics.stdev(values) if len(values) > 1 else None,
    }

def cluster_interval(values, seed=20261015, resamples=2000):
    """Resample target means; candidates and repeat seeds are not independent n."""
    if not values:
        return {"mean": None, "ci95": None, "independent_targets": 0}
    array = np.asarray(values, dtype=float)
    result = {"mean": float(array.mean()), "ci95": None, "independent_targets": len(array)}
    if len(array) < 2:
        result["reason"] = "Fewer than two targets; between-target uncertainty is unidentifiable"
        return result
    rng = np.random.default_rng(seed)
    draws = rng.choice(array, size=(resamples, len(array)), replace=True).mean(axis=1)
    result["ci95"] = [float(x) for x in np.quantile(draws, [0.025, 0.975])]
    return result
