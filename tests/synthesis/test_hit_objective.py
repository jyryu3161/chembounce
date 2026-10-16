"""Original-hit scoring mathematics and robust ranking; no calibrated cutoffs."""

import copy
import json
import math
from dataclasses import replace

import numpy as np
import pytest
from rdkit import Chem

from chembounce_synthesis import conformers, workflow
from chembounce_synthesis.cli import main
from chembounce_synthesis.common import canonical
from chembounce_synthesis.conformers import ConformerSettings, score_candidate
from chembounce_synthesis.qualification import validate_result


@pytest.mark.parametrize("seed", [42, 43, 44])
def test_real_original_hit_matches_duplicated_legacy_reference_and_ignores_actual_design(seed):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    hit = "CCCOCC(=O)Nc1ccccc1"
    product = "CCCOCC(=O)Nc1ccncc1"
    settings = ConformerSettings(num_conformers=8, seed=seed)
    legacy = score_candidate(product, hit, hit, settings)
    original = replace(settings, scoring_objective="original-hit-v1")
    for design in (product, "CCNCC", hit):
        result = score_candidate(product, hit, design, original)
        info = result[0]
        assert info["candidate_conformer_id"] == legacy[0]["candidate_conformer_id"]
        assert info["hit"] == legacy[0]["hit"]
        assert info["score"] == pytest.approx(legacy[0]["score"], abs=1e-14)
        assert info["design_smiles"] == canonical(design)
        cid = info["candidate_conformer_id"]
        assert np.array_equal(
            result[1].GetConformer(cid).GetPositions(), legacy[1].GetConformer(cid).GetPositions()
        )
        validate_result(result, canonical(product), original, hit=canonical(hit), design=canonical(design))
        # Design is a rigid overlay of the chosen candidate conformation.
        a, b = (mol.GetConformer(cid).GetPositions() for mol in result[1:])
        assert np.allclose(
            np.linalg.norm(a[:, None] - a[None, :], axis=-1),
            np.linalg.norm(b[:, None] - b[None, :], axis=-1),
            atol=1e-4,
        )


@pytest.mark.parametrize("tie", [False, True])
def test_both_energy_penalties_first_tie_and_single_secondary_overlay(monkeypatch, tie):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    candidate = conformers.Ensemble(
        Chem.MolFromSmiles("CCC"), {0: 0.0, 1: 1.0}, {0: 0.0, 1: 1.0}, "fixture", {}
    )
    hit = conformers.Ensemble(Chem.MolFromSmiles("CCO"), {0: 0.0}, {0: 0.0}, "fixture", {})
    design = conformers.Ensemble(Chem.MolFromSmiles("CCN"), {0: 0.0}, {0: 0.0}, "fixture", {})
    monkeypatch.setattr(
        conformers, "ensemble", lambda s, settings: {"CCC": candidate, "CCO": hit, "CCN": design}[s]
    )
    settings = ConformerSettings(scoring_objective="original-hit-v1", strain_penalty=0.2)
    calls = []

    def overlay(reference, probe, cid, requested):
        calls.append((reference is design, cid))
        raw = 1.0 if cid == 0 else 0.8
        reference_energy = 4.0 if cid == 0 else 0.0
        adjusted = raw * math.exp(-requested.strain_penalty * reference_energy)
        if tie:
            adjusted = math.exp(requested.strain_penalty * candidate.relative_energies[cid]) * 0.5
        return {
            "reference_adjusted_joint": adjusted,
            "reference_relative_energy": reference_energy,
        }, Chem.Mol(probe.mol)

    monkeypatch.setattr(conformers, "_overlay", overlay)
    info, _, _ = score_candidate("CCC", "CCO", "CCN", settings)
    chosen = 0 if tie else 1
    assert info["candidate_conformer_id"] == chosen
    assert info["score"] == pytest.approx(0.5 if tie else 0.8 * math.exp(-0.2))
    assert calls == [(False, 0), (False, 1), (True, chosen)]


@pytest.mark.parametrize(
    "field", ["objective", "score", "reference_energy", "adjusted", "candidate_id", "ensemble_size"]
)
def test_strict_result_rejects_objective_energy_or_ensemble_drift(field):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    product, hit = "CC(=O)Nc1ccncc1", "CC(=O)Nc1ccccc1"
    settings = ConformerSettings(num_conformers=2, scoring_objective="original-hit-v1")
    result = score_candidate(product, hit, product, settings)
    info = copy.deepcopy(result[0])
    if field == "objective":
        info["scoring_objective"] = "joint-v1"
    elif field == "score":
        info["score"] /= 2
    elif field == "reference_energy":
        info["hit"]["reference_relative_energy"] = -1
    elif field == "adjusted":
        info["hit"]["reference_adjusted_joint"] /= 2
    elif field == "candidate_id":
        info["candidate_conformer_id"] = True
    else:
        info["ensemble_sizes"]["candidate"] = 0
    with pytest.raises(ValueError):
        validate_result((info, *result[1:]), product, settings, hit=hit, design=product)


def row(name, hit, design, *, steps=1, scalar=0.5, spread=0.0):
    return {
        "smiles": name,
        "design": name,
        "family_smiles": "same",
        "hit_tanimoto": hit,
        "hit_shape": hit,
        "hit_pharmacophore": hit,
        "design_tanimoto": design,
        "design_shape": design,
        "design_pharmacophore": design,
        "relative_energy": 0.0,
        "reaction_steps": steps,
        "ranking_policy": "original-hit-v1",
        "ranking_score": scalar,
        "score": scalar,
        "seed_score_median": scalar,
        "seed_score_range": spread,
    }


def test_hit_dominance_excludes_design_and_route_and_uses_robust_metrics():
    stronger = row("strong", 0.8, 0.1, steps=8, scalar=0.8)
    weaker = row("weak", 0.7, 1.0, steps=1, scalar=0.9)
    ranked = workflow._rank([weaker, stronger])
    assert [r["smiles"] for r in ranked] == ["strong", "weak"]
    assert [r["pareto_front"] for r in ranked] == [0, 1]
    assert workflow._family_groups(ranked)["same"][0]["smiles"] == "strong"
    stronger["ranking_metrics"] = {"hit_shape": 0.6, "hit_pharmacophore": 0.6}
    stronger["hit_tanimoto"] = 0.6
    assert workflow._rank([stronger, weaker])[0]["smiles"] == "weak"


def test_hit_scalar_then_stability_precede_design_ties_and_mixed_policies_fail():
    a = row("a", 0.8, 0.1, scalar=0.9, spread=0.2)
    b = row("b", 0.8, 1.0, scalar=0.8)
    assert workflow._rank([b, a])[0] is a
    a["ranking_score"] = b["ranking_score"]
    a["seed_score_range"], b["seed_score_range"] = 0.0, 0.2
    assert workflow._rank([b, a])[0] is a
    b["ranking_policy"] = "joint-v1"
    with pytest.raises(ValueError, match="mix ranking"):
        workflow._rank([a, b])
    with pytest.raises(ValueError, match="requires original-hit"):
        workflow._ranking_policy({"ranking_policy": "original-hit-v1"})


@pytest.mark.parametrize("failure", [None, "failed_seed", "objective_drift"])
def test_real_cli_fixed_three_seed_protocol_exports_robust_ranking_and_exact_pose(
    tmp_path, assets, proposal, monkeypatch, failure
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    hit, design = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"
    saved = tmp_path / "saved.json"
    saved.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": design, "response": proposal}]})
    )
    config = {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(saved)},
        "recommendation": {"mode": "qualified"},
        "minimum_scores": {"hit_tanimoto": 0.01, "hit_shape": 0.01, "hit_pharmacophore": 0.01},
        "conformers": {"num_conformers": 8, "seed": 42},
        "refinement": {"top_n": 1, "num_conformers": 8, "seeds": [42, 43, 44], "max_total_conformers": 24},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "max_scored_products": 1,
        "max_results": 1,
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    out = tmp_path / "out"
    if failure:
        from chembounce_synthesis import refinement

        real = refinement.score_candidate

        def score(*args, **kwargs):
            if args[3].seed == 43:
                if failure == "failed_seed":
                    raise ValueError("Required test seed failed")
                result = real(*args, **kwargs)
                result[0]["scoring_objective"] = "joint-v1"
                return result
            return real(*args, **kwargs)

        monkeypatch.setattr(refinement, "score_candidate", score)
    assert main(
        [
            "--hit",
            hit,
            "--design",
            design,
            "--config",
            str(path),
            "--output",
            str(out),
            "--scoring-objective",
            "original-hit-v1",
        ]
    ) == (3 if failure else 0)
    read = lambda name: json.loads((out / name).read_text())
    report = read("refinement.json")["candidates"][0]
    if failure:
        assert not read("results.json")["results"]
        assert report["status"] == "failed"
        assert [trial["seed"] for trial in report["trials"]] == [42, 43, 44]
        assert report["trials"][1]["status"] == "failed"
        return
    selected = read("results.json")["results"][0]
    assert selected["ranking_policy"] == selected["scoring_objective"] == "original-hit-v1"
    assert selected["ranking_score"] == report["seed_score_min"]
    assert [t["settings"]["seed"] for t in report["trials"]] == [42, 43, 44]
    assert all(
        t["settings"]["num_conformers"] == 8 and t["scoring_objective"] == "original-hit-v1"
        for t in report["trials"]
    )
    for metric in ("shape", "pharmacophore"):
        assert selected["ranking_metrics"][f"hit_{metric}"] == min(t["hit"][metric] for t in report["trials"])
        representative = next(t for t in report["trials"] if t["seed"] == selected["pose_seed"])
        assert selected[f"hit_{metric}"] == representative["hit"][metric]
    assert selected["convergence_status"] == "not_assessed"
    assert read("manifest.json")["ranking_metrics"] == [
        "hit_tanimoto",
        "minimum_seed_hit_shape",
        "minimum_seed_hit_pharmacophore",
    ]
