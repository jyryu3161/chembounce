"""Pure analysis of pinned, collector-validated qualified benchmark attempts.

This module does not qualify chemistry or authenticate arbitrary JSON. Its caller
must trust the pinned collector (or supply a trusted collection digest to the CLI).
No generator acceptance flags, missing attempts, or missing costs are promoted.
"""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path
import statistics

import numpy as np


SCHEMA = "qualified-benchmark-analysis-v1"
STATUSES = {"success", "no_candidates", "timeout", "error", "unsupported", "interrupted", "not_started"}
METRICS = ("hit_tanimoto", "hit_shape", "hit_pharmacophore")
COSTS = ("wall_seconds", "cpu_seconds", "gpu_seconds", "peak_rss_mb", "peak_gpu_mb")
ENDPOINTS = ("at_least_one", "unique_families_at_3_divided_by_3", "at_least_three_families")
BOOTSTRAP_SEED = 20261016
BOOTSTRAP_DRAWS = 2000
SLOT_KEYS = ("method_id", "target_id", "query_id", "seed")


def canonical_sha256(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _text(value):
    return isinstance(value, str) and bool(value.strip())


def _sha(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _number(value):
    return type(value) in (int, float) and math.isfinite(value)


def _slot(row):
    return tuple(row[k] for k in SLOT_KEYS)


def _exact(value):
    return None if value is None else {"value": float(value), "fraction": str(value)}


def _mean(values):
    values = list(values)
    return sum(values, Fraction()) / len(values) if values else None


def validate_plan(plan):
    _require(
        isinstance(plan, dict) and plan.get("schema") == "qualified-benchmark-plan-v1", "Invalid plan schema"
    )
    _require(_text(plan.get("run_id")), "Missing run identity")
    for key in ("protocol_sha256", "collector_source_sha256"):
        _require(_sha(plan.get(key)), f"Invalid {key}")
    methods, queries, seeds = (plan.get(k) for k in ("methods", "queries", "seeds"))
    _require(isinstance(methods, list) and methods and all(_text(x) for x in methods), "Invalid methods")
    _require(len(set(methods)) == len(methods), "Duplicate method")
    _require(
        isinstance(seeds, list) and seeds and all(type(s) is int and s >= 0 for s in seeds), "Invalid seeds"
    )
    _require(len(set(seeds)) == len(seeds), "Duplicate seed")
    _require(isinstance(queries, list) and queries, "Missing queries")
    identities = {}
    for query in queries:
        _require(isinstance(query, dict), "Invalid query")
        _require(all(_text(query.get(k)) for k in ("target_id", "query_id")), "Invalid query identity")
        _require(_sha(query.get("query_sha256")), "Invalid query digest")
        key = query["target_id"], query["query_id"]
        _require(key not in identities, "Duplicate target/query")
        identities[key] = query["query_sha256"]
    slots = [
        (method, target, query, seed) for method, (target, query), seed in product(methods, identities, seeds)
    ]
    if "slot_order" in plan:
        order = plan["slot_order"]
        _require(isinstance(order, list) and len(order) == len(slots), "Invalid planned slot order")
        _require(
            all(isinstance(s, dict) and set(s) == set(SLOT_KEYS) for s in order), "Invalid slot identity"
        )
        _require(all(type(s["seed"]) is int for s in order), "Invalid slot seed")
        _require(set(map(_slot, order)) == set(slots), "Slot order is not the complete planned grid")
        slots = list(map(_slot, order))
    return identities, slots


def _validate_qualified(rows):
    _require(
        isinstance(rows, list) and len(rows) <= 3, "Expected at most three collector-selected recommendations"
    )
    smiles = set()
    for rank, row in enumerate(rows, 1):
        _require(
            isinstance(row, dict) and row.get("rank") == rank and type(row.get("rank")) is int, "Invalid rank"
        )
        _require(
            all(_text(row.get(k)) for k in ("smiles", "family", "scaffold")),
            "Missing recommendation identity",
        )
        _require(row["smiles"] not in smiles, "Duplicate selected molecule")
        smiles.add(row["smiles"])
        if "generic_scaffold" in row:
            _require(
                row["generic_scaffold"] is None or _text(row["generic_scaffold"]), "Invalid generic scaffold"
            )
        bits = row.get("bit_morgan_bits")
        _require(isinstance(bits, list) and bits, "Missing fingerprint bits")
        _require(all(type(b) is int and 0 <= b < 2048 for b in bits), "Invalid fingerprint bits")
        _require(bits == sorted(set(bits)), "Fingerprint bits must be sorted and unique")
        for key in ("metrics", "margins"):
            values = row.get(key)
            _require(isinstance(values, dict) and set(values) == set(METRICS), f"Invalid {key}")
            _require(all(_number(v) and 0 <= v <= 1 for v in values.values()), f"Invalid {key} value")
        _require(
            not ({"accepted", "qualified", "eligible"} & row.keys()),
            "Raw acceptance flags are not analysis evidence",
        )


def _validate_costs(costs):
    _require(
        isinstance(costs, dict) and set(costs) == set(COSTS) | {"native_units", "common_units"},
        "Invalid typed costs",
    )
    for key in COSTS:
        value = costs[key]
        _require(value is None or (_number(value) and value >= 0), f"Invalid {key}")
    for key in ("native_units", "common_units"):
        units = costs[key]
        _require(isinstance(units, dict), "Invalid typed work units")
        _require(
            all(_text(k) and type(v) is int and v >= 0 for k, v in units.items()), "Invalid work unit count"
        )


def validated_attempts(plan, collection):
    """Validate identities and retain every planned denominator; no chemistry inference."""
    identities, slots = validate_plan(plan)
    _require(
        isinstance(collection, dict) and collection.get("schema") == "qualified-benchmark-collection-v1",
        "Invalid collector schema",
    )
    _require(collection.get("plan_sha256") == canonical_sha256(plan), "Stale plan binding")
    _require(
        collection.get("collector_source_sha256") == plan["collector_source_sha256"],
        "Collector source mismatch",
    )
    attempts = collection.get("attempts")
    _require(isinstance(attempts, list), "Invalid attempts")
    expected, indexed = set(slots), {}
    for row in attempts:
        _require(isinstance(row, dict) and all(k in row for k in SLOT_KEYS), "Missing attempt identity")
        _require(
            all(_text(row[k]) for k in SLOT_KEYS[:-1]) and type(row["seed"]) is int,
            "Invalid attempt identity",
        )
        key = _slot(row)
        _require(key in expected, "Extra unplanned attempt")
        _require(key not in indexed, "Duplicate planned attempt")
        _require(row.get("query_sha256") == identities[key[1:3]], "Stale query binding")
        _require(row.get("status") in STATUSES, "Unknown attempt status")
        _validate_qualified(row.get("qualified"))
        _require(
            row["status"] == "success" or not row["qualified"],
            "Failed attempt cannot supply qualified outputs",
        )
        _require(
            not ({"accepted", "qualified_count", "success_rate"} & row.keys()), "Untrusted acceptance summary"
        )
        failures = row.get("failures")
        _require(
            isinstance(failures, list)
            and all(
                isinstance(f, dict) and _text(f.get("stage")) and _text(f.get("reason")) for f in failures
            ),
            "Invalid failure records",
        )
        _validate_costs(row.get("costs"))
        indexed[key] = row
    return [
        indexed.get(key)
        or {
            **dict(zip(SLOT_KEYS, key)),
            "query_sha256": identities[key[1:3]],
            "status": "not_started",
            "qualified": [],
            "failures": [{"stage": "collection", "reason": "missing_planned_attempt"}],
            "costs": {**dict.fromkeys(COSTS), "native_units": {}, "common_units": {}},
        }
        for key in slots
    ]


def attempt_values(row):
    unique = {}
    for item in row["qualified"]:
        unique.setdefault(item["family"], item)
    selected = list(unique.values())
    pairs = list(combinations(selected, 2))
    bits = [set(item["bit_morgan_bits"]) for item in selected]
    bit_pairs = list(combinations(bits, 2))
    generic_available = all(item.get("generic_scaffold") for item in selected)
    return {
        "at_least_one": Fraction(bool(selected)),
        "unique_families_at_3_divided_by_3": Fraction(len(selected), 3),
        "at_least_three_families": Fraction(len(selected) >= 3),
        "bit_morgan_diversity": _mean(1 - Fraction(len(a & b), len(a | b)) for a, b in bit_pairs),
        "scaffold_diversity": _mean(Fraction(a["scaffold"] != b["scaffold"]) for a, b in pairs),
        "generic_scaffold_diversity": _mean(
            Fraction(a["generic_scaffold"] != b["generic_scaffold"]) for a, b in pairs
        )
        if generic_available
        else None,
        "unique_families": len(selected),
        "duplicate_family_outputs": len(row["qualified"]) - len(selected),
    }


def _aggregate(rows, values, key):
    grouped = {}
    for row in rows:
        grouped.setdefault(row["target_id"], {}).setdefault(row["query_id"], []).append(
            values[_slot(row)][key]
        )
    targets, queries = {}, []
    for target, query_groups in sorted(grouped.items()):
        query_means = []
        for query, outcomes in sorted(query_groups.items()):
            observed = [v for v in outcomes if v is not None]
            mean = _mean(observed)
            if mean is not None:
                query_means.append(mean)
            queries.append(
                {
                    "target_id": target,
                    "query_id": query,
                    "estimate": _exact(mean),
                    "supported_attempts": len(observed),
                    "planned_attempts": len(outcomes),
                }
            )
        targets[target] = _mean(query_means)
    supported = [value for value in targets.values() if value is not None]
    return {
        "estimate": _exact(_mean(supported)),
        "targets": {k: _exact(v) for k, v in targets.items()},
        "queries": queries,
        "supported_targets": len(supported),
        "planned_targets": len(targets),
    }, targets


def _cost_report(rows):
    report = {}
    for key in COSTS:
        observed = [row["costs"][key] for row in rows if row["costs"][key] is not None]
        report[key] = {
            "observed_attempts": len(observed),
            "planned_attempts": len(rows),
            "aggregation": "maximum" if key.startswith("peak_") else "sum",
            "observed_value": (max(observed) if key.startswith("peak_") else sum(observed))
            if observed
            else None,
            "complete": len(observed) == len(rows),
        }
    for key in ("native_units", "common_units"):
        names = sorted({name for row in rows for name in row["costs"][key]})
        report[key] = {
            name: {
                "observed_total": sum(row["costs"][key].get(name, 0) for row in rows),
                "observed_attempts": sum(name in row["costs"][key] for row in rows),
                "planned_attempts": len(rows),
            }
            for name in names
        }
    return report


def _distribution(values):
    return {
        "count": len(values),
        "minimum": min(values) if values else None,
        "median": statistics.median(values) if values else None,
        "maximum": max(values) if values else None,
    }


def _interval(values):
    ordered = sorted(values)
    return {
        "lower": _exact(ordered[math.ceil(Fraction(1, 40) * len(ordered)) - 1]),
        "upper": _exact(ordered[math.ceil(Fraction(39, 40) * len(ordered)) - 1]),
    }


def analyze(plan, collection):
    """Return deterministic full-grid estimates and paired target-cluster intervals."""
    rows = validated_attempts(plan, collection)
    values = {_slot(row): attempt_values(row) for row in rows}
    targets = sorted({row["target_id"] for row in rows})
    draws = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED)).integers(
        0, len(targets), size=(BOOTSTRAP_DRAWS, len(targets))
    )
    draws_list = draws.tolist()
    report = {
        "schema": SCHEMA,
        "run_id": plan["run_id"],
        "plan_sha256": canonical_sha256(plan),
        "protocol_sha256": plan["protocol_sha256"],
        "collector_source_sha256": plan["collector_source_sha256"],
        "collection_sha256": canonical_sha256(collection),
        "planned_attempts": len(rows),
        "target_count": len(targets),
        "reported_attempts": len(collection["attempts"]),
        "scope": "Collector-validated computational qualification, not activity, synthetic success, or experimental evidence.",
        "inference_scope": "Descriptive fixed-outcome target intervals only; no confirmatory significance, multiplicity adjustment, or post-hoc winner selection. A two-target development pilot is engineering validation, not comparative inference.",
        "secondary_identity": "Coverage at three requires three unique chemical families, not three exact molecular states.",
        "aggregation": "Equal seeds within query, queries within target, then targets; failures remain zero.",
        "conditional_diversity": "First-ranked representative per family; pairwise mean distance; undefined below two families. Conditional aggregation averages only supported seeds, queries and targets, with support counts.",
        "bootstrap": {
            "unit": "target",
            "target_order": targets,
            "generator": "PCG64",
            "numpy_version": np.__version__,
            "seed": BOOTSTRAP_SEED,
            "draws": BOOTSTRAP_DRAWS,
            "draw_sha256": canonical_sha256(draws_list),
            "interval": "inverse empirical CDF at 1/40 and 39/40; no interpolation",
            "fixed_predictions": True,
            "shared_draws": True,
            "status": "available" if len(targets) >= 2 else "unavailable_single_target",
        },
        "methods": {},
        "paired_differences": [],
        "attempts": [],
    }
    target_estimates, samples = {}, {}
    # Shared draws make an exact ordered target vector sufficient within this call.
    # At most three vectors per method; immutable samples never escape the analysis.
    sample_cache = {}
    for method in plan["methods"]:
        method_rows = [row for row in rows if row["method_id"] == method]
        estimates = {}
        for key in (*ENDPOINTS, "bit_morgan_diversity", "scaffold_diversity", "generic_scaffold_diversity"):
            estimates[key], by_target = _aggregate(method_rows, values, key)
            if key in ENDPOINTS:
                target_estimates[method, key] = [by_target[t] for t in targets]
                if len(targets) >= 2:
                    vector = tuple(target_estimates[method, key])
                    if vector not in sample_cache:
                        sample_cache[vector] = tuple(_mean(vector[i] for i in draw) for draw in draws_list)
                    samples[method, key] = sample_cache[vector]
                else:
                    samples[method, key] = ()
                estimates[key]["interval"] = _interval(samples[method, key]) if len(targets) >= 2 else None
        stages = Counter(f["stage"] for row in method_rows for f in row["failures"])
        stage_attempts = Counter(
            stage for row in method_rows for stage in {f["stage"] for f in row["failures"]}
        )
        reasons = Counter((f["stage"], f["reason"]) for row in method_rows for f in row["failures"])
        report["methods"][method] = {
            "planned_attempts": len(method_rows),
            "estimates": estimates,
            "status_counts": {
                status: sum(row["status"] == status for row in method_rows) for status in sorted(STATUSES)
            },
            "failure_stage_records": dict(sorted(stages.items())),
            "failure_stage_attempts": dict(sorted(stage_attempts.items())),
            "failure_reasons": [
                {"stage": stage, "reason": reason, "count": n}
                for (stage, reason), n in sorted(reasons.items())
            ],
            "costs": _cost_report(method_rows),
            "selected_condition": {
                "scope": "Selected qualified outputs only; not the planned-attempt distribution.",
                **{
                    key: {
                        metric: _distribution(
                            [item[key][metric] for row in method_rows for item in row["qualified"]]
                        )
                        for metric in METRICS
                    }
                    for key in ("metrics", "margins")
                },
            },
            "duplicate_family_outputs": sum(
                values[_slot(row)]["duplicate_family_outputs"] for row in method_rows
            ),
        }
    for left, right in combinations(sorted(plan["methods"]), 2):
        differences = {}
        for key in ENDPOINTS:
            differences[key] = {
                "estimate": _exact(
                    _mean(a - b for a, b in zip(target_estimates[left, key], target_estimates[right, key]))
                ),
                "interval": _interval([a - b for a, b in zip(samples[left, key], samples[right, key])])
                if len(targets) >= 2
                else None,
            }
        report["paired_differences"].append(
            {"left": left, "right": right, "direction": "left_minus_right", "estimates": differences}
        )
    for row in rows:
        report["attempts"].append(
            {
                **{key: row[key] for key in SLOT_KEYS},
                "status": row["status"],
                **{
                    key: _exact(value) if isinstance(value, Fraction) or value is None else value
                    for key, value in values[_slot(row)].items()
                },
            }
        )
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--collection", type=Path, required=True)
    parser.add_argument(
        "--collection-sha256", required=True, help="Trusted externally recorded collector output digest"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    raw = args.collection.read_bytes()
    _require(
        _sha(args.collection_sha256) and hashlib.sha256(raw).hexdigest() == args.collection_sha256,
        "Collection byte digest mismatch",
    )
    result = analyze(json.loads(args.plan.read_bytes()), json.loads(raw))
    result["collection_file_sha256"] = args.collection_sha256
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
