"""Synthetic metadata controls; no saved scientific labels or native geometry."""

import hashlib
import json

import pytest

from benchmarks.research import hit_allocation_comparator as c


def population():
    return [
        dict(
            product=f"p{i}",
            hit="hit",
            design="design",
            hit_tanimoto=0.5,
            first_eligible_route_index=i,
            observations=[],
        )
        for i in range(13)
    ]


def saved(tmp_path, candidate, seed, status, name):
    value = dict(
        product=candidate["product"],
        hit=candidate["hit"],
        design=candidate["design"],
        seed=seed,
        status=status,
    )
    if status == "feasible":
        value["info"] = {"controlled_fixture": True}
    else:
        value.update(reason="no_jointly_feasible_hit_pose", audit={"pairs": [{"status": "rejected"}]})
    path = tmp_path / name
    path.write_text(json.dumps(value))
    return dict(
        path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(), seed=seed, origin="coverage"
    )


def test_ranking_uses_only_original_hit_score_then_first_route():
    rows = population()
    original = [r["product"] for r in rows[:12]]
    rows[-1]["hit_tanimoto"] = 0.9
    expected = ["p12", *original[:11]]
    assert c.select(rows, original)[c.POLICY] == expected
    assert c.select(list(reversed(rows)), original)[c.POLICY] == expected
    for row in rows:
        row.update(design="different", observations=[{"all_three_feasible": row["product"] == "p11"}])
    assert c.select(rows, original)[c.POLICY] == expected
    assert c.select(rows, original)["original_saved"] == original


def test_missing_seed_and_duplicate_disagreement_remain_unknown(tmp_path):
    row = population()[0]
    row["observations"] = [saved(tmp_path, row, s, "feasible", f"s{s}.json") for s in (42, 43)]
    result = c.join_outcomes(row)
    assert result["all_three_feasible"] is None and not result["assessment_complete"]
    row["observations"].append(saved(tmp_path, row, 44, "feasible", "s44.json"))
    assert c.join_outcomes(row)["all_three_feasible"] is True
    row["observations"].append(saved(tmp_path, row, 42, "bounded_infeasible", "duplicate.json"))
    assert c.join_outcomes(row)["all_three_feasible"] is None


def test_changed_or_missing_saved_bytes_are_unknown(tmp_path):
    row = population()[0]
    binding = saved(tmp_path, row, 42, "feasible", "one.json")
    path = tmp_path / "one.json"
    path.write_text("{}")
    assert c.read_observation(binding, row)["status"] == "unknown"
    path.unlink()
    assert c.read_observation(binding, row)["status"] == "unknown"


def test_baseline_only_not_experimental_rescue(tmp_path):
    row = population()[0]
    value = dict(
        product=row["product"],
        hit=row["hit"],
        design=row["design"],
        seed=42,
        baseline={
            "status": "abstained",
            "reason": "no_jointly_feasible_hit_pose",
            "audit": {"pairs": [{"status": "rejected"}]},
        },
        prototype={"status": "feasible", "info": {}},
    )
    path = tmp_path / "paired.json"
    path.write_text(json.dumps(value))
    binding = dict(
        path=str(path),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        seed=42,
        origin="paired_baseline",
    )
    assert c.read_observation(binding, row)["status"] == "bounded_infeasible"


@pytest.mark.parametrize("origin", ["coverage", "paired_baseline"])
@pytest.mark.parametrize(
    "audit",
    [
        None,
        [],
        "invalid",
        {"pairs": None},
        {"pairs": {}},
        {"pairs": [None]},
        {"pairs": [{"status": "rejected"}, None]},
    ],
)
def test_malformed_audit_stays_unknown_without_aborting_comparison(tmp_path, monkeypatch, origin, audit):
    rows = population()
    row = rows[0]
    binding = saved(tmp_path, row, 42, "bounded_infeasible", "malformed.json")
    path = tmp_path / "malformed.json"
    condition = json.loads(path.read_text())
    condition["audit"] = audit
    if origin == "paired_baseline":
        observed = {key: condition.pop(key) for key in ("status", "reason", "audit")}
        condition["baseline"] = observed
    path.write_text(json.dumps(condition))
    binding.update(origin=origin, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    row["observations"] = [
        binding,
        *[saved(tmp_path, row, seed, "feasible", f"s{seed}.json") for seed in (43, 44)],
    ]
    monkeypatch.setattr(c, "describe", lambda _: {"error": "Controlled metadata-only test"})

    assert c.read_observation(binding, row)["status"] == "unknown"
    result = c.compare(rows, [r["product"] for r in rows[:12]])
    assert len(result["outcomes"]) == 13
    assert result["outcomes"][row["product"]]["all_three_feasible"] is None
    assert not result["outcomes"][row["product"]]["assessment_complete"]
    for arm in result["arms"].values():
        assert arm["quota"] == 12 and arm["eligible_products"] == 13
        assert row["product"] in arm["unknown_products"]


def test_full_denominator_and_identity_failures_are_reported(monkeypatch):
    rows = population()
    original = [r["product"] for r in rows[:12]]

    def identity(smiles):
        if smiles == "p0":
            raise ValueError("controlled incomplete normalization")
        return dict(
            family_smiles=smiles,
            scaffold_family_smiles="scaffold",
            ring_skeleton_smiles="ring",
            has_scaffold=True,
        )

    monkeypatch.setattr(c, "describe", identity)
    result = c.compare(rows, original)
    assert len(result["outcomes"]) == 13 and not result["full_population_assessment_complete"]
    for arm in result["arms"].values():
        assert arm["eligible_products"] == 13 and arm["quota"] == 12
        assert arm["exact_product_coverage"] == 12 / 13
        assert len(arm["unknown_products"]) == 12 and not arm["assessment_complete"]
        assert arm["selected_diversity"]["identity_unknown"] == 1
        assert arm["all_three_feasible_observed"] == 0
