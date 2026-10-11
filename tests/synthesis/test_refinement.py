import pytest

from chembounce_synthesis.conformers import ConformerSettings
from chembounce_synthesis.refinement import RefinementSettings, refine_candidate


def test_representative_is_one_median_pose_not_best_or_mixed(monkeypatch):
    def score(smiles, hit, design, settings, pose):
        values = {42: 0.95, 43: 0.6, 44: 0.8}
        return (
            {"score": values[settings.seed], "settings": {"seed": settings.seed}},
            f"hit-{settings.seed}",
            f"design-{settings.seed}",
        )

    monkeypatch.setattr("chembounce_synthesis.refinement.score_candidate", score)
    result, report = refine_candidate("CC", "CC", "CC", ConformerSettings(), RefinementSettings(top_n=1))
    assert result[1:] == ("hit-44", "design-44")
    assert report["pose_seed"] == 44
    assert report["seed_score_median"] == 0.8
    assert report["seed_score_min"] == 0.6
    assert report["seed_score_max"] == 0.95
    assert len(report["trials"]) == 3


def test_failed_seed_cannot_be_hidden_by_successful_seed(monkeypatch):
    def score(smiles, hit, design, settings, pose):
        if settings.seed == 43:
            raise ValueError("No converged conformers")
        return {"score": 0.9}, None, None

    monkeypatch.setattr("chembounce_synthesis.refinement.score_candidate", score)
    result, report = refine_candidate("CC", "CC", "CC", ConformerSettings(), RefinementSettings(top_n=1))
    assert result is None
    assert report["status"] == "failed"
    assert report["seed_successes"] == 2
    assert report["trials"][1]["reason"] == "No converged conformers"


def test_even_seed_count_keeps_lower_middle_actual_pose(monkeypatch):
    def score(smiles, hit, design, settings, pose):
        return {"score": {42: 0.1, 43: 0.2}[settings.seed]}, settings.seed, settings.seed

    monkeypatch.setattr("chembounce_synthesis.refinement.score_candidate", score)
    result, report = refine_candidate(
        "CC", "CC", "CC", ConformerSettings(), RefinementSettings(top_n=1, seeds=[42, 43])
    )
    assert result[0]["score"] == 0.1
    assert report["pose_seed"] == 42
    assert report["seed_score_median"] == pytest.approx(0.15)


@pytest.mark.parametrize(
    "config",
    [
        {"top_n": -1},
        {"top_n": True},
        {"seeds": []},
        {"seeds": [42, 42]},
        {"seeds": [True]},
        {"seeds": [-1]},
        {"num_conformers": 0},
    ],
)
def test_invalid_refinement_settings(config):
    with pytest.raises(ValueError):
        RefinementSettings(**config)
