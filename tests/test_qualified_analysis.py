"""Synthetic independent arithmetic and schema adversaries; no chemistry/models."""

import copy
from fractions import Fraction
import hashlib
import json

import numpy as np
import pytest

from benchmarks.research import qualified_analysis as analysis


def plan():
    return {
        "schema": "qualified-benchmark-plan-v1",
        "run_id": "synthetic",
        "protocol_sha256": "a" * 64,
        "collector_source_sha256": "b" * 64,
        "methods": ["a", "b"],
        "seeds": [42, 43],
        "queries": [
            {"target_id": t, "query_id": q, "query_sha256": sha * 64}
            for t, q, sha in [("A", "q1", "1"), ("A", "q2", "2"), ("B", "q3", "3")]
        ],
    }


def qualified(family="f1", rank=1, bits=None, scaffold=None):
    return {
        "rank": rank,
        "smiles": f"molecule_{family}_{rank}",
        "family": family,
        "scaffold": scaffold or f"scaffold_{family}",
        "generic_scaffold": f"generic_{family}",
        "bit_morgan_bits": bits or [1, 2],
        "metrics": dict.fromkeys(analysis.METRICS, 0.7),
        "margins": dict.fromkeys(analysis.METRICS, 0.2),
    }


def row(p, method="a", target="A", query="q1", seed=42, outputs=1, status=None):
    return {
        "method_id": method,
        "target_id": target,
        "query_id": query,
        "seed": seed,
        "query_sha256": next(
            q["query_sha256"] for q in p["queries"] if q["target_id"] == target and q["query_id"] == query
        ),
        "status": status or ("success" if outputs else "no_candidates"),
        "qualified": [qualified(f"f{i}", i + 1) for i in range(outputs)],
        "failures": [],
        "costs": {
            "wall_seconds": 10,
            "cpu_seconds": 8,
            "gpu_seconds": None,
            "peak_rss_mb": 100,
            "peak_gpu_mb": None,
            "native_units": {"samples": 2},
            "common_units": {"forward_calls": 3},
        },
    }


def collection(p, rows):
    return {
        "schema": "qualified-benchmark-collection-v1",
        "plan_sha256": analysis.canonical_sha256(p),
        "collector_source_sha256": p["collector_source_sha256"],
        "attempts": rows,
    }


def test_hierarchy_retains_missing_failures_and_secondary_slots():
    p = plan()
    rows = [
        row(p, outputs=3),
        row(p, query="q2", outputs=1),
        row(p, target="B", query="q3", outputs=0, status="timeout"),
    ]
    report = analysis.analyze(p, collection(p, rows))
    a = report["methods"]["a"]
    # A: ((1+0)/2 + (1+0)/2)/2=1/2; B:0; equal targets=1/4.
    assert a["estimates"]["at_least_one"]["estimate"]["fraction"] == "1/4"
    assert a["estimates"]["unique_families_at_3_divided_by_3"]["estimate"]["fraction"] == "1/6"
    assert a["estimates"]["at_least_three_families"]["estimate"]["fraction"] == "1/8"
    assert a["status_counts"]["not_started"] == 3
    assert a["status_counts"]["timeout"] == 1
    assert report["planned_attempts"] == 12
    assert report["reported_attempts"] == 3
    assert report["methods"]["b"]["estimates"]["at_least_one"]["estimate"]["fraction"] == "0"
    assert a["estimates"]["at_least_one"]["targets"] == {
        "A": {"value": 0.5, "fraction": "1/2"},
        "B": {"value": 0.0, "fraction": "0"},
    }


def test_diversity_deduplicates_families_and_reports_missing_support():
    p = plan()
    r = row(p, outputs=3)
    r["qualified"] = [
        qualified("same", 1, [1, 2]),
        qualified("same", 2, [4, 5]),
        qualified("other", 3, [2, 3]),
    ]
    values = analysis.attempt_values(r)
    assert values["unique_families_at_3_divided_by_3"] == Fraction(2, 3)
    assert values["at_least_three_families"] == 0
    assert values["bit_morgan_diversity"] == Fraction(2, 3)
    assert values["scaffold_diversity"] == 1
    assert values["duplicate_family_outputs"] == 1
    r["qualified"][2].pop("generic_scaffold")
    report = analysis.analyze(p, collection(p, [r]))
    estimates = report["methods"]["a"]["estimates"]
    assert estimates["bit_morgan_diversity"]["estimate"]["fraction"] == "2/3"
    assert estimates["bit_morgan_diversity"]["supported_targets"] == 1
    assert estimates["generic_scaffold_diversity"]["estimate"] is None
    assert report["methods"]["b"]["estimates"]["scaffold_diversity"]["estimate"] is None
    assert analysis.attempt_values(row(p))["bit_morgan_diversity"] is None


def test_costs_unknown_not_zero_and_failure_costs_retained():
    p = plan()
    r = row(p, outputs=0, status="error")
    r["failures"] = [{"stage": "route", "reason": "budget"}, {"stage": "route", "reason": "unsupported"}]
    report = analysis.analyze(p, collection(p, [r]))["methods"]["a"]
    costs = report["costs"]
    assert costs["wall_seconds"] == {
        "observed_attempts": 1,
        "planned_attempts": 6,
        "aggregation": "sum",
        "observed_value": 10,
        "complete": False,
    }
    assert costs["gpu_seconds"]["observed_value"] is None
    assert costs["gpu_seconds"]["observed_attempts"] == 0
    assert costs["peak_rss_mb"]["aggregation"] == "maximum"
    assert costs["native_units"]["samples"] == {
        "observed_total": 2,
        "observed_attempts": 1,
        "planned_attempts": 6,
    }
    assert report["failure_stage_records"]["route"] == 2
    assert report["failure_stage_attempts"]["route"] == 1


def test_bootstrap_exact_shared_draw_oracle_and_no_mutation():
    p = plan()
    rows = [row(p, method=m, target="B", query="q3", seed=s) for m in p["methods"] for s in p["seeds"]]
    rows += [row(p, method="a", query=q, seed=s) for q in ("q1", "q2") for s in p["seeds"]]
    c = collection(p, rows)
    before = copy.deepcopy((p, c))
    report = analysis.analyze(p, c)
    assert (p, c) == before
    assert report == analysis.analyze(p, c)
    draws = np.random.Generator(np.random.PCG64(20261016)).integers(0, 2, size=(2000, 2)).tolist()
    assert report["bootstrap"]["draw_sha256"] == analysis.canonical_sha256(draws)
    # a target outcomes [1,1], b [0,1]: paired differences [1,0].
    oracle = sorted(Fraction(sum(i == 0 for i in draw), 2) for draw in draws)
    paired = report["paired_differences"][0]
    assert paired["direction"] == "left_minus_right"
    estimate = paired["estimates"]["at_least_one"]
    assert estimate["estimate"]["fraction"] == "1/2"
    assert estimate["interval"]["lower"]["fraction"] == str(oracle[49])
    assert estimate["interval"]["upper"]["fraction"] == str(oracle[1949])


def test_identical_methods_have_exact_zero_paired_intervals():
    p = plan()
    report = analysis.analyze(p, collection(p, [row(p, method=m) for m in p["methods"]]))
    for value in report["paired_differences"][0]["estimates"].values():
        assert value["estimate"]["fraction"] == "0"
        assert value["interval"]["lower"]["fraction"] == "0"
        assert value["interval"]["upper"]["fraction"] == "0"


def test_single_target_interval_unavailable_and_empty_collection_zero():
    p = plan()
    p["queries"] = p["queries"][:2]
    report = analysis.analyze(p, collection(p, []))
    assert report["bootstrap"]["status"] == "unavailable_single_target"
    assert report["methods"]["a"]["estimates"]["at_least_one"]["interval"] is None
    assert report["methods"]["a"]["estimates"]["at_least_one"]["estimate"]["fraction"] == "0"


@pytest.mark.parametrize(
    "mutation",
    [
        lambda p, c: c.update(schema="raw-generator-output"),
        lambda p, c: c.update(plan_sha256="f" * 64),
        lambda p, c: c.update(collector_source_sha256="f" * 64),
        lambda p, c: c["attempts"].append(copy.deepcopy(c["attempts"][0])),
        lambda p, c: c["attempts"][0].update(method_id="extra"),
        lambda p, c: c["attempts"][0].update(query_sha256="f" * 64),
        lambda p, c: c["attempts"][0].update(status="timeout"),
        lambda p, c: c["attempts"][0].update(status="unrecognized"),
        lambda p, c: c["attempts"][0].update(seed=True),
        lambda p, c: c["attempts"][0].update(accepted=True),
        lambda p, c: c["attempts"][0]["qualified"][0].update(accepted=True),
        lambda p, c: c["attempts"][0]["qualified"][0].update(rank=True),
        lambda p, c: c["attempts"][0]["qualified"][0].update(rank=2),
        lambda p, c: c["attempts"][0]["qualified"][0].update(bit_morgan_bits=[2, 1]),
        lambda p, c: c["attempts"][0]["qualified"][0].update(bit_morgan_bits=[1, 1]),
        lambda p, c: c["attempts"][0]["qualified"][0].update(bit_morgan_bits=[2048]),
        lambda p, c: c["attempts"][0]["qualified"][0].update(bit_morgan_bits=[True]),
        lambda p, c: c["attempts"][0]["qualified"][0]["metrics"].update(hit_shape=float("nan")),
        lambda p, c: c["attempts"][0]["qualified"][0]["margins"].update(hit_shape=-0.01),
        lambda p, c: c["attempts"][0]["costs"].update(wall_seconds=-1),
        lambda p, c: c["attempts"][0]["costs"].update(gpu_seconds=float("inf")),
        lambda p, c: c["attempts"][0]["costs"]["native_units"].update(samples=True),
        lambda p, c: c["attempts"][0].update(failures=[{"stage": "route"}]),
    ],
)
def test_reject_stale_forged_or_invalid_attempt(mutation):
    p = plan()
    c = collection(p, [row(p)])
    mutation(p, c)
    with pytest.raises(ValueError):
        analysis.analyze(p, c)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda p: p.update(methods=["a", "a"]),
        lambda p: p.update(methods=[]),
        lambda p: p.update(seeds=[42, 42]),
        lambda p: p.update(seeds=[True]),
        lambda p: p.update(queries=[]),
        lambda p: p["queries"].append(copy.deepcopy(p["queries"][0])),
        lambda p: p.pop("protocol_sha256"),
        lambda p: p.update(slot_order=[]),
    ],
)
def test_invalid_plan(mutation):
    p = plan()
    mutation(p)
    with pytest.raises(ValueError):
        analysis.analyze(p, collection(p, []))


def test_duplicate_molecule_rejected_even_different_family():
    p = plan()
    r = row(p, outputs=2)
    r["qualified"][1]["smiles"] = r["qualified"][0]["smiles"]
    with pytest.raises(ValueError, match="Duplicate selected"):
        analysis.analyze(p, collection(p, [r]))


def test_slot_order_preserved_and_bound():
    p = plan()
    _, slots = analysis.validate_plan(p)
    p["slot_order"] = [dict(zip(analysis.SLOT_KEYS, s)) for s in reversed(slots)]
    c = collection(p, [])
    report = analysis.analyze(p, c)
    assert tuple(report["attempts"][0][k] for k in analysis.SLOT_KEYS) == slots[-1]
    p["slot_order"][0] = p["slot_order"][1]
    with pytest.raises(ValueError, match="complete planned grid"):
        analysis.analyze(p, collection(p, []))


def test_cli_requires_trusted_digest_and_refuses_overwrite(tmp_path):
    p = plan()
    plan_path, collection_path, output = (
        tmp_path / name for name in ("plan.json", "collection.json", "report.json")
    )
    plan_path.write_text(json.dumps(p))
    collection_path.write_text(json.dumps(collection(p, [])))
    digest = hashlib.sha256(collection_path.read_bytes()).hexdigest()
    args = [
        "--plan",
        str(plan_path),
        "--collection",
        str(collection_path),
        "--output",
        str(output),
        "--collection-sha256",
    ]
    with pytest.raises(ValueError, match="byte digest"):
        analysis.main(args + ["f" * 64])
    assert not output.exists()
    analysis.main(args + [digest])
    assert json.loads(output.read_text())["collection_file_sha256"] == digest
    with pytest.raises(FileExistsError):
        analysis.main(args + [digest])
