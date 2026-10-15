"""Reject fixed-reference 3D failures before search without changing stereo policy."""

import json
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pytest
from rdkit import Chem

from benchmarks.research import adapters
from chembounce_synthesis import conformers, workflow
from chembounce_synthesis.common import canonical

HIT = "CC(=O)Nc1ccccc1"
DESIGN = "CC(=O)Nc1ccncc1"
UNSPECIFIED = "CC(O)c1ccccc1"
CAGE = "NC12CC3CC(CC(C3)C1)C2"


def config_for(tmp_path, assets, proposal):
    saved = tmp_path / "proposals.json"
    saved.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": DESIGN, "response": proposal}]})
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(saved)},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "conformers": {"num_conformers": 2, "max_iterations": 1000},
        "max_scored_products": 1,
    }


def forbid_expensive_work(monkeypatch):
    calls = {}
    for name in (
        "Catalog",
        "Rules",
        "load_carried_routes",
        "_evidence_evaluator",
        "_provider",
        "resolve_proposal",
        "score_candidate",
    ):
        calls[name] = Mock(side_effect=AssertionError(f"Unexpected expensive work: {name}"))
        monkeypatch.setattr(workflow, name, calls[name])
    embed = Mock(side_effect=AssertionError("Preflight must not embed or resolve stereo"))
    monkeypatch.setattr(conformers.AllChem, "EmbedMultipleConfs", embed)
    calls["embedding"] = embed
    return calls


@pytest.mark.parametrize("hit", [UNSPECIFIED, CAGE])
def test_workflow_rejects_reference_stereo_before_assets_and_preserves_error_artifacts(
    tmp_path, assets, proposal, monkeypatch, hit
):
    config = config_for(tmp_path, assets, proposal)
    monkeypatch.setattr(conformers, "rdShapeAlign", object())
    calls = forbid_expensive_work(monkeypatch)
    output = tmp_path / "run"
    with pytest.raises(ValueError, match="requires specified stereochemistry"):
        workflow.run(hit, [DESIGN], config, output)
    assert not any(call.called for call in calls.values())
    summary = json.loads((output / "summary.json").read_text())
    assert summary["status"] == "error" and summary["error_type"] == "ValueError"
    assert "specified stereochemistry" in summary["error"]
    assert summary["scored_products"] == summary["accepted_products"] == 0
    assert json.loads((output / "routes.json").read_text())["routes"] == []
    assert json.loads((output / "scores.json").read_text())["scores"] == []
    assert not (output / "results.tsv").exists()


def test_workflow_rejects_missing_alignment_before_assets(tmp_path, assets, proposal, monkeypatch):
    config = config_for(tmp_path, assets, proposal)
    monkeypatch.setattr(conformers, "rdShapeAlign", None)
    calls = forbid_expensive_work(monkeypatch)
    with pytest.raises(ValueError, match="requires an RDKit build with rdShapeAlign"):
        workflow.run(HIT, [DESIGN], config, tmp_path / "run")
    assert not any(call.called for call in calls.values())


def test_empty_design_run_keeps_existing_no_scoring_behavior(tmp_path, assets, proposal, monkeypatch):
    config = config_for(tmp_path, assets, proposal)
    monkeypatch.setattr(conformers, "rdShapeAlign", None)
    assert workflow.run(UNSPECIFIED, [], config, tmp_path / "empty") == []
    summary = json.loads((tmp_path / "empty/summary.json").read_text())
    assert summary["status"] == "no_candidates" and summary["scored_products"] == 0


class ReachedExpensiveWork(Exception):
    pass


@pytest.mark.parametrize(
    "hit,design,settings,mode",
    [
        (HIT, DESIGN, {}, "3d"),
        ("C[C@H](O)c1ccccc1", DESIGN, {}, "3d"),
        (UNSPECIFIED, DESIGN, {"reject_unspecified_stereo": False}, "3d"),
        (CAGE, DESIGN, {"reject_unspecified_stereo": False}, "3d"),
        (UNSPECIFIED, DESIGN, {}, "2d"),
        (CAGE, DESIGN, {}, "2d"),
    ],
)
def test_workflow_preserves_valid_permissive_and_2d_paths(
    tmp_path, assets, proposal, monkeypatch, hit, design, settings, mode
):
    config = config_for(tmp_path, assets, proposal)
    config["conformers"].update(settings)
    monkeypatch.setattr(conformers, "rdShapeAlign", None if mode == "2d" else object())
    calls = forbid_expensive_work(monkeypatch)
    calls["load_carried_routes"].side_effect = ReachedExpensiveWork
    with pytest.raises(ReachedExpensiveWork):
        workflow.run(hit, [design], config, tmp_path / "run", scoring_mode=mode)
    calls["load_carried_routes"].assert_called_once()
    assert not calls["embedding"].called


def write_pose(path, smiles, *, is_3d=True):
    # Deliberately collinear coordinates: do not assign stereo or embed a cage.
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(is_3d)
    for index in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(index, (float(index), 0.0, 0.0))
    mol.AddConformer(conf)
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    return path


@pytest.mark.parametrize(
    "pose_smiles,hit,allow,is_3d",
    [
        ("CO", "CO", False, True),
        ("CO", "CC", False, True),
        ("CO", "CO", False, False),
        (UNSPECIFIED, UNSPECIFIED, False, True),
        (UNSPECIFIED, UNSPECIFIED, True, True),
        (CAGE, CAGE, False, True),
        (CAGE, CAGE, True, True),
    ],
)
def test_bound_reference_preflight_matches_existing_loader_semantics(
    tmp_path, monkeypatch, pose_smiles, hit, allow, is_3d
):
    path = write_pose(tmp_path / "bound.sdf", pose_smiles, is_3d=is_3d)
    settings = conformers.ConformerSettings(reject_unspecified_stereo=not allow)
    monkeypatch.setattr(conformers, "rdShapeAlign", object())
    monkeypatch.setattr(
        conformers.AllChem,
        "EmbedMultipleConfs",
        Mock(side_effect=AssertionError("Bound reference must never be embedded")),
    )
    try:
        expected = conformers.reference_from_sdf(path, hit, settings)
    except ValueError as exc:
        with pytest.raises(ValueError) as failure:
            conformers.prepare_reference(hit, settings, path)
        assert str(failure.value) == str(exc)
    else:
        actual = conformers.prepare_reference(hit, settings, path)
        assert actual.bound_pose and actual.energies == expected.energies
        assert canonical(Chem.MolToSmiles(actual.mol)) == canonical(hit)
        assert np.array_equal(
            actual.mol.GetConformer().GetPositions(), expected.mol.GetConformer().GetPositions()
        )


def test_mismatched_bound_reference_is_rejected_before_assets(tmp_path, assets, proposal, monkeypatch):
    config = config_for(tmp_path, assets, proposal)
    config["reference_sdf"] = str(write_pose(tmp_path / "wrong.sdf", "CO"))
    monkeypatch.setattr(conformers, "rdShapeAlign", object())
    calls = forbid_expensive_work(monkeypatch)
    with pytest.raises(ValueError, match="exact chemical state and stereochemistry"):
        workflow.run(HIT, [DESIGN], config, tmp_path / "run")
    assert not any(call.called for call in calls.values())


@pytest.mark.parametrize(
    "hit,allow,mode", [(HIT, False, "3d"), (UNSPECIFIED, True, "3d"), (UNSPECIFIED, False, "2d")]
)
def test_real_scoring_still_returns_candidates_on_supported_paths(
    tmp_path, assets, proposal, hit, allow, mode
):
    if mode == "3d":
        pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)
    config["conformers"]["reject_unspecified_stereo"] = not allow
    rows = workflow.run(hit, [DESIGN], config, tmp_path / "run", scoring_mode=mode)
    assert rows and all(row["score"] > 0 for row in rows)
    assert json.loads((tmp_path / "run/summary.json").read_text())["status"] == "success"
    assert (tmp_path / "run/poses_hit.sdf").exists() == (mode == "3d")


def pipeline_request(tmp_path, config, hit=UNSPECIFIED):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    return {
        "query": {"smiles": hit},
        "seed": 42,
        "settings": {"generator_options": {"calculate_shape": False}},
        "budget": {"max_candidates": 1, "max_evaluations": 1},
        "assets": {"catalog": config["catalog"], "rules": config["rules"], "synthesis_config": str(path)},
        "output_dir": str(tmp_path / "pipeline"),
    }


@pytest.mark.parametrize("failure", ["stereo", "alignment", "bound_state"])
def test_pipeline_rejects_invalid_reference_before_generator(
    tmp_path, assets, proposal, monkeypatch, failure
):
    config = config_for(tmp_path, assets, proposal)
    request = pipeline_request(tmp_path, config, HIT if failure != "stereo" else UNSPECIFIED)
    if failure == "bound_state":
        request["query"]["reference_sdf"] = str(write_pose(tmp_path / "wrong.sdf", "CO"))
    monkeypatch.setattr(conformers, "rdShapeAlign", None if failure == "alignment" else object())
    generator = Mock(side_effect=AssertionError("Invalid reference reached upstream generation"))
    monkeypatch.setattr("chembounce.chembounce", generator)
    calls = forbid_expensive_work(monkeypatch)
    expected = {
        "stereo": "specified stereochemistry",
        "alignment": "rdShapeAlign",
        "bound_state": "exact chemical state",
    }
    with pytest.raises(ValueError, match=expected[failure]):
        adapters.pipeline(request)
    assert not generator.called and not any(call.called for call in calls.values())
    assert not Path(request["output_dir"]).exists()


@pytest.mark.parametrize("variant", ["valid", "allow", "override_allow", "2d", "design_only"])
def test_pipeline_preserves_valid_permissive_and_2d_generator_paths(
    tmp_path, assets, proposal, monkeypatch, variant
):
    config = config_for(tmp_path, assets, proposal)
    if variant == "allow":
        config["conformers"]["reject_unspecified_stereo"] = False
    request = pipeline_request(tmp_path, config, HIT if variant == "valid" else UNSPECIFIED)
    if variant == "override_allow":
        request["settings"]["synthesis_overrides"] = {"conformers": {"reject_unspecified_stereo": False}}
    if variant == "2d":
        request["settings"]["disable_3d"] = True
    if variant == "design_only":
        request["settings"]["disable_synthesis"] = True
        request["assets"].pop("synthesis_config")
    monkeypatch.setattr(conformers, "rdShapeAlign", None if variant in {"2d", "design_only"} else object())
    generator = Mock(side_effect=ReachedExpensiveWork)
    monkeypatch.setattr("chembounce.chembounce", generator)
    with pytest.raises(ReachedExpensiveWork):
        adapters.pipeline(request)
    generator.assert_called_once()
