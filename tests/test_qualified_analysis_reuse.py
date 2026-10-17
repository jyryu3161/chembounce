"""Cache-specific synthetic metadata tests; existing independent arithmetic oracles remain unchanged."""

import copy
import json
import sys

import pytest

from benchmarks.research import qualified_analysis as analysis


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def fixture(kind):
    groups = [2, 1, 3] if kind == "unequal_groups" else [1, 1, 1]
    if kind == "single_target":
        groups = [2]
    p = {
        "schema": "qualified-benchmark-plan-v1",
        "run_id": "synthetic-bootstrap-" + kind,
        "protocol_sha256": "a" * 64,
        "collector_source_sha256": "b" * 64,
        "methods": ["a", "b", "c"],
        "seeds": [42, 43, 97],
        "queries": [
            {"target_id": f"t{t}", "query_id": f"q{q}", "query_sha256": f"{10 * t + q + 1:064x}"}
            for t, count in enumerate(groups)
            for q in range(count)
        ],
    }
    rows = []
    for mi, method in enumerate(p["methods"]):
        for query in p["queries"]:
            t, q = int(query["target_id"][1:]), int(query["query_id"][1:])
            for si, seed in enumerate(p["seeds"]):
                if kind == "missing" or (kind == "unequal_groups" and (mi + q + si) % 4 == 0):
                    continue
                n = {
                    "zero": 0,
                    "repeated_nonzero": [1, 2, 3][t],
                    "permuted": ([0, 1, 3], [3, 1, 0], [0, 1, 3])[mi][t],
                    "distinct": (mi + t + si) % 4,
                    "unequal_groups": (mi + t + q + si) % 4,
                    "single_target": (mi + q + si) % 4,
                }[kind]
                qualified = [
                    {
                        "rank": i + 1,
                        "smiles": f"SYNTHETIC_NOT_PARSED_{i}",
                        "family": f"family{i}",
                        "scaffold": f"scaffold{i}",
                        "generic_scaffold": f"generic{i}" if si else None,
                        "bit_morgan_bits": [i + 1, i + 2],
                        "metrics": dict.fromkeys(analysis.METRICS, 0.7),
                        "margins": dict.fromkeys(analysis.METRICS, 0.2),
                    }
                    for i in range(n)
                ]
                rows.append(
                    {
                        "method_id": method,
                        **query,
                        "seed": seed,
                        "status": "success" if n else "timeout",
                        "qualified": qualified,
                        "failures": [] if n else [{"stage": "native", "reason": "synthetic_timeout"}],
                        "costs": {
                            **dict.fromkeys(analysis.COSTS),
                            "native_units": {} if si else {"samples": 0},
                            "common_units": {"calls": si},
                        },
                    }
                )
    c = {
        "schema": "qualified-benchmark-collection-v1",
        "plan_sha256": analysis.canonical_sha256(p),
        "collector_source_sha256": p["collector_source_sha256"],
        "attempts": rows,
    }
    return p, c


def test_ordered_nonzero_vectors_keep_paired_uncertainty():
    p, c = fixture("permuted")
    actual = analysis.analyze(p, c)
    estimates = actual["methods"]
    a = estimates["a"]["estimates"]["at_least_one"]
    b = estimates["b"]["estimates"]["at_least_one"]
    assert a["estimate"] == b["estimate"]
    assert [v["fraction"] for v in a["targets"].values()] == ["0", "1", "1"]
    assert [v["fraction"] for v in b["targets"].values()] == ["1", "1", "0"]
    # A cache that sorts the values would conflate these equal-marginal vectors.
    paired = actual["paired_differences"][0]["estimates"]["at_least_one"]
    assert paired["estimate"]["fraction"] == "0"
    assert paired["interval"]["lower"]["fraction"] != "0"
    assert paired["interval"]["upper"]["fraction"] != "0"


@pytest.mark.parametrize(
    "kind,unique_vectors", [("missing", 1), ("zero", 1), ("repeated_nonzero", 3), ("permuted", 6)]
)
def test_bootstrap_mean_calls_reuse_only_identical_ordered_vectors(monkeypatch, kind, unique_vectors):
    p, c = fixture(kind)
    original = analysis._mean
    calls = []

    def count_sample_means(values):
        # Bootstrap samples call _mean from a comprehension; aggregation and
        # paired point estimates call it from ordinary functions.
        caller = sys._getframe(1)
        if caller.f_code.co_name in {"<listcomp>", "<genexpr>"}:
            calls.append(caller.f_code.co_filename)
        return original(values)

    monkeypatch.setattr(analysis, "_mean", count_sample_means)
    analysis.analyze(p, c)
    assert len(calls) == unique_vectors * analysis.BOOTSTRAP_DRAWS
    assert len(calls) < len(p["methods"]) * len(analysis.ENDPOINTS) * analysis.BOOTSTRAP_DRAWS
    # An unchanged next checkpoint still owns a new cache: no cross-call reuse.
    calls.clear()
    analysis.analyze(p, c)
    assert len(calls) == unique_vectors * analysis.BOOTSTRAP_DRAWS


def test_changed_consecutive_calls_and_mutable_report_values_are_isolated():
    p, c = fixture("repeated_nonzero")
    first = analysis.analyze(p, c)
    expected = canonical(first)
    first["methods"]["a"]["estimates"]["at_least_one"]["interval"]["lower"]["fraction"] = "mutated"
    # Change outcomes without changing the plan, then change the target set.
    changed = copy.deepcopy(c)
    changed["attempts"] = []
    empty = analysis.analyze(p, changed)
    assert empty["methods"]["a"]["estimates"]["at_least_one"]["estimate"]["fraction"] == "0"
    assert empty["planned_attempts"] == first["planned_attempts"]
    single = analysis.analyze(*fixture("single_target"))
    assert single["bootstrap"]["status"] == "unavailable_single_target"
    assert single["methods"]["a"]["estimates"]["at_least_one"]["interval"] is None
    assert canonical(analysis.analyze(p, c)) == expected
