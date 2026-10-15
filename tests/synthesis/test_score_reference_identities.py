"""Score identities stay attributable when one product has multiple references."""

import copy
import csv
import json

import pytest
from rdkit import Chem

from chembounce_synthesis import workflow
from chembounce_synthesis.common import canonical


HIT_INPUT = "O=C(C)Nc1ccccc1"
DESIGN_INPUTS = ["O=C(C)Nc1ccncc1", "O=C(C)Nc1cccnc1"]


@pytest.mark.parametrize("mode,refine", [("2d", False), ("3d", False), ("3d", True)])
def test_score_exports_identify_exact_references_and_preserve_metrics(
    tmp_path, assets, proposal, monkeypatch, mode, refine
):
    if mode == "3d":
        pytest.importorskip("rdkit.Chem.rdShapeAlign")
    proposals = tmp_path / "proposals.json"
    proposals.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "proposals": [{"design": design, "response": proposal} for design in DESIGN_INPUTS],
            }
        )
    )
    config = {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(proposals)},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "conformers": {"num_conformers": 2, "max_iterations": 1000},
        "refinement": {"top_n": int(refine), "num_conformers": 2, "seeds": [42, 43, 44]},
        "max_scored_products": 2,
    }
    original_metrics = {}
    real_score, real_refine = workflow.score_candidate, workflow.refine_candidate

    def score(product, hit, design, *args, **kwargs):
        result = real_score(product, hit, design, *args, **kwargs)
        original_metrics["screening", product, design] = copy.deepcopy(result[0])
        return result

    def refined(product, hit, design, *args, **kwargs):
        result, report = real_refine(product, hit, design, *args, **kwargs)
        if result is not None:
            original_metrics["refined", product, design] = copy.deepcopy(result[0])
        return result, report

    monkeypatch.setattr(workflow, "score_candidate", score)
    monkeypatch.setattr(workflow, "refine_candidate", refined)
    out = tmp_path / "output"
    results = workflow.run(HIT_INPUT, DESIGN_INPUTS, config, out, scoring_mode=mode)
    manifest = json.loads((out / "manifest.json").read_text())
    scores = json.loads((out / "scores.json").read_text())["scores"]
    routes = {row["route_id"]: row for row in json.loads((out / "routes.json").read_text())["routes"]}
    screening = [row for row in scores if row["scoring_stage"] != "refined"]
    expected_designs = {canonical(design) for design in DESIGN_INPUTS}
    assert len(screening) == 2
    assert {row["product"] for row in screening} == {canonical(DESIGN_INPUTS[0])}
    assert {row["design_smiles"] for row in screening} == expected_designs
    assert manifest["hit"] == canonical(HIT_INPUT) != HIT_INPUT
    assert set(manifest["designs"]) == expected_designs
    assert {row["scoring_stage"] for row in scores} == (
        {"screening", "refined"} if refine else {"2d_screening" if mode == "2d" else "screening"}
    )
    for row in scores:
        assert row["hit_smiles"] == manifest["hit"]
        assert row["design_smiles"] == routes[row["route_id"]]["design"]
        assert row["design_smiles"] in manifest["designs"]
        assert row["product"] == routes[row["route_id"]]["verification"]["product"]
        if mode == "3d":
            original = original_metrics[row["scoring_stage"], row["product"], row["design_smiles"]]
            assert isinstance(row["hit"], dict) and isinstance(row["design"], dict)
            for key, value in original.items():
                assert row[key] == value
        else:
            assert row["design"] == row["design_smiles"]  # Preserve legacy 2D representation.
            assert row["scoring_mode"] == "2d"
    with (out / "results.tsv").open() as stream:
        table = list(csv.DictReader(stream, delimiter="\t"))
    assert len(results) == len(table) == 1  # One chemical family, two design assessments.
    for result, exported in zip(results, table):
        score_row = next(
            row
            for row in scores
            if row["route_id"] == result["route_id"] and row["scoring_stage"] == result["scoring_stage"]
        )
        assert score_row["design_smiles"] == result["design"] == exported["design"]
        assert score_row["product"] == result["smiles"] == exported["smiles"]
        assert score_row["score"] == result["score"] == float(exported["score"])
        if mode == "3d":
            for role in ("hit", "design"):
                for metric in ("shape", "pharmacophore"):
                    assert score_row[role][metric] == result[f"{role}_{metric}"]
                reference = next(iter(Chem.SDMolSupplier(str(out / f"references_{role}.sdf"))))
                assert canonical(Chem.MolToSmiles(reference)) == score_row[f"{role}_smiles"]
