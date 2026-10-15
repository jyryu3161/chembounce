"""Configured bound-pose checks reach screening and every refinement seed."""

import json

import pytest
from rdkit import Chem

pytest.importorskip("rdkit.Chem.rdShapeAlign")

from chembounce_synthesis import refinement, workflow
from chembounce_synthesis.conformers import ConformerSettings, ensemble

HIT, DESIGN = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"


def config_for(tmp_path, assets, proposal, pattern):
    path = tmp_path / "proposals.json"
    path.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": DESIGN, "response": proposal}]})
    )
    source = tmp_path / "synthetic-reference.sdf"
    reference = ensemble(HIT, ConformerSettings(num_conformers=2))
    with Chem.SDWriter(str(source)) as writer:
        writer.write(reference.mol, confId=min(reference.relative_energies))
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "conformers": {"num_conformers": 2},
        "reference_sdf": str(source),
        "pocket": {
            "anchors": [
                {
                    "id": "specified-anchor",
                    "kind": "atom",
                    "reference_atom_index": 0,
                    "candidate_smarts": pattern,
                    "candidate_atom_index": 0,
                    "max_distance": 100,
                    "required": True,
                }
            ]
        },
    }


def test_impossible_required_anchor_rejects_screened_products(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, proposal, "[S]")
    out = tmp_path / "reject"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    records = json.loads((out / "routes.json").read_text())["routes"]
    assert records[0]["selection_status"] == "3d_failed"
    assert "pocket" in records[0]["reason"].lower() or "anchor" in records[0]["reason"].lower()
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["pocket"]["reference_sdf_sha256"]


def test_pocket_constraints_are_applied_to_every_refinement_seed(tmp_path, assets, proposal, monkeypatch):
    config = config_for(tmp_path, assets, proposal, "[C]")
    config["refinement"] = {"top_n": 1, "num_conformers": 2, "seeds": [42, 43, 44]}
    real_score = refinement.score_candidate
    received = []

    def checked_score(*args, **kwargs):
        assert kwargs["pose_constraints"] is not None
        received.append(args[3].seed)
        return real_score(*args, **kwargs)

    monkeypatch.setattr(refinement, "score_candidate", checked_score)
    out = tmp_path / "refined"
    rows = workflow.run(HIT, [DESIGN], config, out)
    assert rows and received == [42, 43, 44]
    report = json.loads((out / "refinement.json").read_text())["candidates"][0]
    assert all(trial["pocket"]["passed"] for trial in report["trials"])


def test_pocket_requires_explicit_reference_before_loading_models(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, proposal, "[C]")
    del config["reference_sdf"]
    with pytest.raises(ValueError, match="reference_sdf"):
        workflow._config(config)
