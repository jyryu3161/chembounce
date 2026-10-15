"""Independent decision cases for uncertainty, adaptive effort and shared budgets."""

from dataclasses import asdict

import pytest

from chembounce_synthesis.conformers import ConformerSettings
from chembounce_synthesis.refinement import RefinementBudget, RefinementSettings, refine_candidate


def scored(settings, score, shape=None):
    return (
        {
            "score": score,
            "hit": {"shape": score if shape is None else shape, "pharmacophore": score},
            "design": {"shape": score, "pharmacophore": score},
            "candidate_relative_energy": 0.0,
            "settings": asdict(settings),
        },
        f"hit-{settings.seed}",
        f"design-{settings.seed}",
    )


def test_all_seeds_threshold_rejects_lucky_median_and_reports_metric_range(monkeypatch):
    def score(smiles, hit, design, settings, pose):
        return scored(settings, {42: 0.2, 43: 0.8, 44: 0.9}[settings.seed])

    monkeypatch.setattr("chembounce_synthesis.refinement.score_candidate", score)
    args = ("CC", "CC", "CC", ConformerSettings())
    result, report = refine_candidate(*args, RefinementSettings(top_n=1), minimum_scores={"hit_shape": 0.7})
    assert result is None
    assert report["status"] == "threshold_rejected"
    assert report["seed_score_median"] == 0.8
    assert report["metric_ranges"]["hit_shape"] == {"min": 0.2, "max": 0.9}
    result, report = refine_candidate(
        *args,
        RefinementSettings(top_n=1, threshold_policy="representative"),
        minimum_scores={"hit_shape": 0.7},
    )
    assert result[0]["score"] == 0.8


def test_adaptive_stage_uses_final_complete_seed_block_and_keeps_actual_pose(monkeypatch):
    calls = []

    def score(smiles, hit, design, settings, pose):
        calls.append((settings.num_conformers, settings.seed))
        values = (
            {42: 0.3, 43: 0.9, 44: 0.8} if settings.num_conformers == 2 else {42: 0.82, 43: 0.80, 44: 0.81}
        )
        return scored(settings, values[settings.seed])

    monkeypatch.setattr("chembounce_synthesis.refinement.score_candidate", score)
    config = RefinementSettings(
        top_n=1,
        num_conformers=2,
        adaptive_max_conformers=8,
        max_total_conformers=100,
        max_score_range=0.05,
        ranking_statistic="minimum",
    )
    result, report = refine_candidate("CC", "CC", "CC", ConformerSettings(), config)
    assert calls == [(n, seed) for n in (2, 4) for seed in (42, 43, 44)]
    assert report["candidate_conformer_requests"] == 18
    assert report["convergence_status"] == "converged"
    assert report["ranking_score"] == 0.8
    assert result[0]["score"] == 0.81
    assert result[1:] == ("hit-44", "design-44")
    assert report["stages"][0]["convergence_status"] == "unstable"
    assert report["num_conformers"] == 4


def test_shared_budget_reserves_other_family_initial_evaluation(monkeypatch):
    calls = []

    def score(smiles, hit, design, settings, pose):
        calls.append((smiles, settings.num_conformers, settings.seed))
        return scored(settings, (0.1 if settings.seed == 42 else 0.9) if smiles == "CC" else 0.85)

    monkeypatch.setattr("chembounce_synthesis.refinement.score_candidate", score)
    config = RefinementSettings(
        top_n=2, num_conformers=2, adaptive_max_conformers=8, max_total_conformers=18, max_score_range=0.05
    )
    budget = RefinementBudget(config.max_total_conformers, reserved_for_others=6)
    first, first_report = refine_candidate("CC", "CCC", "CCCC", ConformerSettings(), config, budget=budget)
    assert first is None and first_report["status"] == "budget_exhausted"
    assert len(calls) == 3  # 12-request escalation would steal the reserved 6.
    budget.reserved_for_others = 0
    second, second_report = refine_candidate("CCC", "CCC", "CCCC", ConformerSettings(), config, budget=budget)
    assert second is not None
    assert second_report["convergence_status"] == "converged"
    assert budget.requested == 12
    assert budget.score_calls == 6


def test_borderline_criterion_triggers_bounded_more_effort_even_with_stable_score(monkeypatch):
    def score(smiles, hit, design, settings, pose):
        return scored(settings, 0.9, shape=0.71 if settings.num_conformers == 2 else 0.8)

    monkeypatch.setattr("chembounce_synthesis.refinement.score_candidate", score)
    config = RefinementSettings(
        num_conformers=2, adaptive_max_conformers=4, max_total_conformers=18, borderline_margin=0.02
    )
    result, report = refine_candidate(
        "CC", "CC", "CC", ConformerSettings(), config, minimum_scores={"hit_shape": 0.7}
    )
    assert result is not None
    assert report["adaptive_stages"] == 2
    assert result[0]["hit"]["shape"] == 0.8


def test_unstable_at_maximum_fails_even_if_high_median(monkeypatch):
    monkeypatch.setattr(
        "chembounce_synthesis.refinement.score_candidate",
        lambda s, h, d, settings, p: scored(settings, 0.2 if settings.seed == 42 else 0.99),
    )
    result, report = refine_candidate(
        "CC", "CC", "CC", ConformerSettings(), RefinementSettings(max_score_range=0.05)
    )
    assert result is None
    assert report["status"] == "unstable"
    assert report["seed_score_median"] == 0.99


@pytest.mark.parametrize(
    "config",
    [
        {"adaptive_max_conformers": 48},
        {"adaptive_max_conformers": 12},
        {"max_score_range": float("nan")},
        {"max_score_range": 0.1, "seeds": [42]},
        {"max_total_conformers": True},
        {"conformer_growth_factor": 1},
        {"borderline_margin": -1},
        {"threshold_policy": "best"},
        {"ranking_statistic": "best"},
    ],
)
def test_uncertainty_settings_fail_early(config):
    with pytest.raises(ValueError):
        RefinementSettings(**config)


def test_metric_instability_cannot_hide_behind_stable_composite_score(monkeypatch):
    def score(smiles, hit, design, settings, pose):
        return scored(settings, 0.8, shape=0.2 if settings.seed == 42 else 0.99)

    monkeypatch.setattr("chembounce_synthesis.refinement.score_candidate", score)
    result, report = refine_candidate(
        "CC", "CC", "CC", ConformerSettings(), RefinementSettings(max_metric_range=0.1)
    )
    assert result is None and report["status"] == "unstable"
    assert report["seed_score_range"] == 0
    assert report["metric_ranges"]["hit_shape"]["max"] - report["metric_ranges"]["hit_shape"]["min"] == 0.79
