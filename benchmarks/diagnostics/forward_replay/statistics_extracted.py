from collections import Counter, defaultdict

import statistics

import numpy as np

SPLITS = ("scaffold", "time", "target")

PRIMARY = "attempt_with_common_forward_matched_route"

INFRA = "infrastructure_interruption"

COUNTS = (
    "candidate_occurrences",
    "no_route_supplied",
    "route_supplied",
    "route_assessed",
    "route_all_steps_evaluated",
    "route_all_steps_matched",
    "route_assessed_incomplete",
    "route_evaluated_mismatched",
    "route_audit_failed",
)

COVERAGE = ("route_supplied", "route_assessed", "route_all_steps_evaluated", "route_all_steps_matched")

BOUNDED = (
    (PRIMARY,)
    + tuple(f"attempt_with_{name}" for name in COVERAGE)
    + tuple(f"{name}_output_slot_fraction" for name in COVERAGE)
)

CONDITIONAL = "matched_fraction_given_evaluable_routes"

METRICS = COUNTS + BOUNDED + (CONDITIONAL,)

TERMINAL = {
    "success",
    "no_candidates",
    "timeout",
    "memory_limit",
    "budget_exceeded",
    "worker_error",
    "missing_response",
    "invalid_response",
    "launch_error",
    "interrupted",
    INFRA,
}

class ValidationError(ValueError):
    def __init__(self, message, details=None):
        self.details = details or {}
        super().__init__(message)

def require(condition, message, details=None):
    if not condition:
        raise ValidationError(message, details)

def run_key(row):
    return row["method_id"], row["query_id"], row["seed"]

def mean(values):
    values = [v for v in values if v is not None]
    return statistics.mean(values) if values else None

def bootstrap(values, config):
    values = np.array([v for v in values if v is not None], dtype=float)
    result = {
        "mean": float(values.mean()) if len(values) else None,
        "target_clusters": len(values),
        "ci95": None,
    }
    if len(values) > 1:
        draws = (
            np.random.default_rng(config["seed"])
            .choice(values, (config["draws"], len(values)), replace=True)
            .mean(axis=1)
        )
        result["ci95"] = np.quantile(draws, config["quantiles"]).tolist()
    return result

def summarize(rows, targets, names, config, *, bounds=False, planned_attempts=None):
    planned_attempts = len(rows) if planned_attempts is None else planned_attempts
    per_query = {
        q: {name: mean(r["metrics"].get(name) for r in rows if r["query_id"] == q) for name in names}
        for q in targets
    }
    per_target = {
        t: {name: mean(per_query[q][name] for q in targets if targets[q] == t) for name in names}
        for t in sorted(set(targets.values()))
    }
    eligibility = {
        name: {
            "attempts": sum(r["metrics"].get(name) is not None for r in rows),
            "queries": sum(q[name] is not None for q in per_query.values()),
            "targets": sum(t[name] is not None for t in per_target.values()),
        }
        for name in names
    }
    for values in eligibility.values():
        values.update(
            excluded_attempts=planned_attempts - values["attempts"],
            excluded_queries=len(targets) - values["queries"],
            excluded_targets=len(set(targets.values())) - values["targets"],
        )
    return {
        "attempts": len(rows),
        "planned_attempts": planned_attempts,
        "excluded_attempts": planned_attempts - len(rows),
        "metric_eligibility": eligibility,
        "per_query_attempts": {q: sum(r["query_id"] == q for r in rows) for q in targets},
        "per_query": per_query,
        "per_target": per_target,
        "target_macro": {
            name: (
                mean(t[name] for t in per_target.values())
                if bounds
                else bootstrap((t[name] for t in per_target.values()), config)
            )
            for name in names
        },
    }

def paired_rows(left, right, *, complete_case=False):
    a = {(r["query_id"], r["seed"]): r for r in left}
    b = {(r["query_id"], r["seed"]): r for r in right}
    require(
        a.keys() == b.keys() and len(a) == len(left) and len(b) == len(right),
        "Paired query/seed grids differ",
    )
    rows = []
    for key, row in a.items():
        other = b[key]
        if complete_case and INFRA in (row["status"], other["status"]):
            continue
        values = {
            name: (
                row["metrics"][name] - other["metrics"][name]
                if row["metrics"][name] is not None and other["metrics"][name] is not None
                else None
            )
            for name in METRICS
        }
        rows.append({**row, "metrics": values})
    return rows

def aggregate_split(runs, occurrences, targets, caps, config):
    """Independent seed/query/target aggregation; callable on synthetic records."""
    by_attempt = defaultdict(list)
    for row in occurrences:
        by_attempt[row["method_id"], row["query_id"], row["generation_seed"]].append(row)
    methods = defaultdict(list)
    occurrence_methods = defaultdict(list)
    for run in runs:
        output = by_attempt[run_key(run)]
        require(run["status"] == "success" or not output, "Failed attempt has forward occurrences")
        c = Counter(r["status"] for r in output)
        evaluated, matched = sum(r["evaluated"] for r in output), sum(r["matched"] for r in output)
        metrics = {
            "candidate_occurrences": len(output),
            "no_route_supplied": c["no_route_supplied"],
            "route_supplied": len(output) - c["no_route_supplied"],
            "route_assessed": c["assessed"],
            "route_all_steps_evaluated": evaluated,
            "route_all_steps_matched": matched,
            "route_assessed_incomplete": c["assessed"] - evaluated,
            "route_evaluated_mismatched": evaluated - matched,
            "route_audit_failed": c["audit_failed"],
        }
        cap = caps[run["method_id"]]
        require(
            0 <= len(output) <= cap and min(metrics.values()) >= 0, "Invalid coverage counts/candidate cap"
        )
        metrics[PRIMARY] = int(matched > 0)
        metrics.update({f"attempt_with_{name}": int(metrics[name] > 0) for name in COVERAGE})
        metrics.update({f"{name}_output_slot_fraction": metrics[name] / cap for name in COVERAGE})
        metrics[CONDITIONAL] = matched / evaluated if evaluated else None
        methods[run["method_id"]].append(
            {
                "query_id": run["query_id"],
                "target_id": run["target_id"],
                "seed": run["seed"],
                "status": run["status"],
                "metrics": metrics,
            }
        )
        occurrence_methods[run["method_id"]].extend(output)
    result, complete, bounds, paired, paired_complete, paired_bounds = {}, {}, {}, {}, {}, {}
    for method, rows in sorted(methods.items()):
        counts = {name: sum(r["metrics"][name] for r in rows) for name in COUNTS}
        output = occurrence_methods[method]
        input_keys = [key for row in output for key in row["input_keys"]]
        result[method] = {
            "status_counts": dict(Counter(r["status"] for r in rows)),
            "coverage_counts": counts,
            "unique_supplied_route_keys": len({r["route_key"] for r in output if r["route_key"]}),
            "route_step_input_occurrences": len(input_keys),
            "unique_reactant_context_input_keys": len(set(input_keys)),
            "pooled_occurrence_matched_fraction_selected_subset": (
                counts["route_all_steps_matched"] / counts["route_all_steps_evaluated"]
                if counts["route_all_steps_evaluated"]
                else None
            ),
            **summarize(rows, targets, METRICS, config),
        }
        observed = [r for r in rows if r["status"] != INFRA]
        complete[method] = {
            "excluded_infrastructure_attempts": len(rows) - len(observed),
            **summarize(observed, targets, METRICS, config, planned_attempts=len(rows)),
        }
        bound_rows = [
            {
                **r,
                "metrics": {
                    f"{name}_{side}": (value if r["status"] != INFRA else missing)
                    for name in BOUNDED
                    for side, missing in [("lower", 0), ("upper", 1)]
                    for value in [r["metrics"][name]]
                },
            }
            for r in rows
        ]
        bound_names = [f"{name}_{side}" for name in BOUNDED for side in ("lower", "upper")]
        bounds[method] = summarize(bound_rows, targets, bound_names, config, bounds=True)
    if "pipeline" in methods:
        for method, rows in sorted(methods.items()):
            if method == "pipeline":
                continue
            primary = methods["pipeline"]
            paired[method] = summarize(paired_rows(primary, rows), targets, METRICS, config)
            paired_complete[method] = summarize(
                paired_rows(primary, rows, complete_case=True),
                targets,
                METRICS,
                config,
                planned_attempts=len(primary),
            )
            lookup = {(r["query_id"], r["seed"]): r for r in rows}
            differences = []
            for row in primary:
                other = lookup[row["query_id"], row["seed"]]
                values = {}
                for name in BOUNDED:
                    low_a, high_a = (0, 1) if row["status"] == INFRA else (row["metrics"][name],) * 2
                    low_b, high_b = (0, 1) if other["status"] == INFRA else (other["metrics"][name],) * 2
                    values[f"{name}_lower"], values[f"{name}_upper"] = low_a - high_b, high_a - low_b
                differences.append({**row, "metrics": values})
            paired_bounds[method] = summarize(differences, targets, bound_names, config, bounds=True)
    input_targets, input_queries = defaultdict(set), defaultdict(set)
    for occurrence in occurrences:
        for key in occurrence["input_keys"]:
            input_targets[key].add(occurrence["target_id"])
            input_queries[key].add(occurrence["query_id"])
    dependence = {
        "unique_reactant_context_input_keys": len(input_targets),
        "input_keys_shared_across_targets": sum(len(values) > 1 for values in input_targets.values()),
        "input_keys_shared_across_queries": sum(len(values) > 1 for values in input_queries.values()),
        "uncertainty_policy": "Occurrences/routes/reactant reuse are not independent samples; target-cluster intervals are descriptive, with cross-target shared-input dependence and few-target limitations.",
    }
    return {
        "methods": result,
        "paired_pipeline_minus_comparator": paired,
        "dependence_diagnostics": dependence,
        "infrastructure_identification_bounds_not_ci": {
            "methods": bounds,
            "paired_pipeline_minus_comparator": paired_bounds,
        },
        "complete_case_secondary_selection_biased": {
            "methods": complete,
            "paired_pipeline_minus_comparator": paired_complete,
        },
    }
