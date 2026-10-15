"""Independent evidence applies to each route before choosing product provenance."""

import json

import pytest

pytest.importorskip("rdkit.Chem.rdShapeAlign")

from chembounce_synthesis import workflow
from chembounce_synthesis.common import canonical
from chembounce_synthesis.evidence import SynthesisEvidenceEvaluator

HIT, DESIGN = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"


def config_for(tmp_path, assets, proposals):
    path = tmp_path / "proposals.json"
    path.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": DESIGN, "response": p} for p in proposals]})
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "conformers": {"num_conformers": 2},
        "synthesis_evidence": {"policy": {"require_forward_agreement": True}},
    }


def test_no_predictor_cannot_pass_requested_independent_agreement(tmp_path, assets, proposal):
    out = tmp_path / "unknown"
    assert workflow.run(HIT, [DESIGN], config_for(tmp_path, assets, [proposal]), out) == []
    records = json.loads((out / "routes.json").read_text())["routes"]
    assert records[0]["selection_status"] == "evidence_rejected"
    assert records[0]["synthesis_evidence"]["steps"][0]["forward_prediction"]["status"] == "unknown"
    summary = json.loads((out / "summary.json").read_text())
    assert summary["scored_products"] == 0
    assert summary["synthesis_evidence_usage"]["predictor_calls"] == 0


def test_different_reactant_route_can_rescue_same_exact_product(tmp_path, assets, proposal, monkeypatch):
    # Two real RDKit replay routes make the identical amide from different acyl
    # donors. The independent test predictor rejects the first donor only.
    catalog, rules = assets
    with catalog.path.open("a") as stream:
        stream.write("acetic-anhydride\tCC(=O)OC(C)=O\ttest fixture\t2026-10-15\n")
    rule_data = json.loads(rules.path.read_text())
    rule_data["rules"].append(
        {
            "id": "anhydride",
            "smarts": "[C:1](=[O:2])[O:3][C:4](=[O:5])[C:6].[N;H1,H2:7]>>[C:1](=[O:2])[N:7]",
            "source": "test alternative acyl-donor rule",
        }
    )
    rules.path.write_text(json.dumps(rule_data))
    alternative = {
        "reactions": [
            {"reaction_template": "anhydride", "reactants": ["CC(=O)OC(C)=O", "Nc1ccncc1"], "product": DESIGN}
        ],
        "building_blocks": ["CC(=O)OC(C)=O", "Nc1ccncc1"],
    }

    class Predictor:
        provenance = {"provider": "controlled_independent_route_predictor", "uses_conditions": False}

        def token_budget(self, reactants, reagents, conditions):
            return {"input_tokens": 8, "max_output_tokens": 24}

        def predict(self, reactants, reagents, conditions):
            product = DESIGN if canonical("CC(=O)OC(C)=O") in reactants else "CCC(=O)Nc1ccncc1"
            return {"products": [product], "complete": [True], "input_tokens": 8, "output_tokens": 12}

    engine = SynthesisEvidenceEvaluator(Predictor(), {"require_forward_agreement": True})
    monkeypatch.setattr(workflow, "_evidence_evaluator", lambda config: engine)
    out = tmp_path / "rescued"
    rows = workflow.run(HIT, [DESIGN], config_for(tmp_path, assets, [proposal, alternative]), out)
    assert len(rows) == 1 and rows[0]["smiles"] == DESIGN
    records = json.loads((out / "routes.json").read_text())["routes"]
    assert records[0]["selection_status"] == "evidence_rejected"
    selected = next(r for r in records if r["route_id"] == rows[0]["route_id"])
    assert selected["route"]["nodes"][-1]["rule_id"] == "anhydride"
    assert selected["synthesis_evidence"]["eligibility"]["passed"]
    assert records[0]["route_id"] not in rows[0]["alternative_route_ids"]
    assert engine.usage["predictor_calls"] == 2


def test_evidence_model_path_resolves_from_config_directory(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, [proposal])
    config["synthesis_evidence"]["predictor"] = {"kind": "reactiont5", "model_path": "models/forward"}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    result = workflow.load_config(path)
    assert result["synthesis_evidence"]["predictor"]["model_path"] == str(tmp_path / "models/forward")


def test_independent_prediction_budget_is_shared_before_alternate_routes(
    tmp_path, assets, proposal, monkeypatch
):
    catalog, rules = assets
    with catalog.path.open("a") as stream:
        stream.write("acetic-anhydride\tCC(=O)OC(C)=O\ttest fixture\t2026-10-15\n")
    rule_data = json.loads(rules.path.read_text())
    rule_data["rules"].append(
        {
            "id": "anhydride",
            "smarts": "[C:1](=[O:2])[O:3][C:4](=[O:5])[C:6].[N;H1,H2:7]>>[C:1](=[O:2])[N:7]",
            "source": "test alternative acyl-donor rule",
        }
    )
    rules.path.write_text(json.dumps(rule_data))
    alternative = {
        "reactions": [
            {"reaction_template": "anhydride", "reactants": ["CC(=O)OC(C)=O", "Nc1ccncc1"], "product": DESIGN}
        ],
        "building_blocks": ["CC(=O)OC(C)=O", "Nc1ccncc1"],
    }
    other = "CCC(=O)Nc1ccncc1"
    second = {
        "reactions": [
            {"reaction_template": "amide", "reactants": ["CCC(=O)O", "Nc1ccncc1"], "product": other}
        ],
        "building_blocks": ["CCC(=O)O", "Nc1ccncc1"],
    }
    config = config_for(tmp_path, assets, [proposal, alternative])
    path = tmp_path / "proposals.json"
    data = json.loads(path.read_text())
    data["proposals"].append({"design": other, "response": second})
    path.write_text(json.dumps(data))

    class Predictor:
        provenance = {"provider": "controlled_budget_predictor", "uses_conditions": False}

        def token_budget(self, reactants, reagents, conditions):
            return {"input_tokens": 8, "max_output_tokens": 24}

        def predict(self, reactants, reagents, conditions):
            product = other if "CCC(=O)O" in reactants else DESIGN
            return {"products": [product], "complete": [True], "input_tokens": 8, "output_tokens": 12}

    engine = SynthesisEvidenceEvaluator(
        Predictor(), {"require_forward_agreement": True}, max_predictor_calls=2
    )
    monkeypatch.setattr(workflow, "_evidence_evaluator", lambda config: engine)
    out = tmp_path / "shared"
    rows = workflow.run(HIT, [DESIGN, other], config, out)
    assert {row["smiles"] for row in rows} == {DESIGN, other}
    records = json.loads((out / "routes.json").read_text())["routes"]
    unfunded = next(record for record in records if record["route"]["nodes"][-1]["rule_id"] == "anhydride")
    prediction = unfunded["synthesis_evidence"]["steps"][0]["forward_prediction"]
    assert prediction["status"] == "unknown"
    assert "budget" in prediction["reason"]
    assert engine.usage["predictor_calls"] == 2
