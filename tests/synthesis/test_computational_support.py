"""Strict route evidence integration with controlled raw model beams, never inference."""

import copy
import csv
import json
import time

import pytest
from rdkit import Chem

from chembounce_synthesis import workflow
from chembounce_synthesis.computational_support import ComputationalRouteSupport
from chembounce_synthesis.evidence import evaluator_from_config, validate_evidence_config
from chembounce_synthesis.forward_support import PROFILE, PinnedReactionT5Forward, strict_smiles

HIT, DESIGN = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"


@pytest.fixture
def concrete_origins():
    from rdkit.Chem import rdChemReactions

    rxn = rdChemReactions.ReactionFromSmarts("[C:1]-[O:2]>>[C:1]")
    product = rxn.RunReactants((Chem.MolFromSmiles("CO"),))[0][0]
    if not product.GetAtomWithIdx(0).HasProp("react_idx"):
        pytest.skip(
            "Installed RDKit lacks concrete reaction source-slot metadata; strict support fails closed"
        )


def config(tmp_path, assets, proposal):
    path = tmp_path / "proposals.json"
    path.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": DESIGN, "response": proposal}]})
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "recommendation": {"mode": "qualified"},
        "minimum_scores": {"hit_tanimoto": 0.1, "hit_shape": 0.01, "hit_pharmacophore": 0.01},
        "conformers": {"num_conformers": 2},
        "synthesis_evidence": {
            "profile": PROFILE,
            "policy": {"require_forward_agreement": True},
            "predictor": {"kind": "reactiont5", "model_path": str(tmp_path / "controlled-model")},
        },
        "search": {"beam_width": 1, "building_block_top_k": 1},
        "max_scored_products": 1,
    }


def test_real_route_audit_cache_and_fresh_final_comparison(assets, route, strict_predictor, concrete_origins):
    evaluator = ComputationalRouteSupport(strict_predictor)
    first = evaluator.evaluate_route(route, *assets)
    assert first["eligibility"]["passed"]
    assert first["route_audit"]["strict_eligible"]
    assert evaluator.validate_route(first, route, *assets)
    assert len(strict_predictor.calls) == 1
    assert evaluator.evaluate_route(route, *assets)["eligibility"]["passed"]
    assert len(strict_predictor.calls) == 1 and evaluator.usage["cache_hits"] == 1
    for mutate in (
        lambda e: e["steps"][0]["forward_prediction"]["raw"]["predictions"][0].update(raw_product="CC"),
        lambda e: e.update(route_binding_sha256="supplied-claim"),
        lambda e: e["route_audit"]["steps"][0]["reactant_slots"].reverse(),
        lambda e: e["steps"][0]["forward_prediction"].update(matched_rank=2),
    ):
        forged = copy.deepcopy(first)
        mutate(forged)
        with pytest.raises(ValueError):
            evaluator.validate_route(forged, route, *assets)


@pytest.mark.parametrize(
    "products,status",
    [
        (["CC", "CO", "CN"], "mismatched"),
        (["CC", [0, 70], "CN"], "partially_evaluated"),
        (["CC", "invalid", "CN"], "partially_evaluated"),
        ([DESIGN + ".", "CC", "CO"], "partially_evaluated"),
        ([DESIGN + ".O", "CC", "CO"], "mismatched"),
        ([DESIGN, [0, 2, 1], [0, 70]], "matched"),
    ],
)
def test_all_beams_retained_no_partial_disagreement_or_component_repair(
    assets, route, strict_predictor, products, status, concrete_origins
):
    strict_predictor.products = products
    evidence = ComputationalRouteSupport(strict_predictor).evaluate_route(route, *assets)
    prediction = evidence["steps"][0]["forward_prediction"]
    assert prediction["status"] == status
    assert len(prediction["predictions"]) == 3
    assert evidence["eligibility"]["passed"] is (status == "matched")
    assert all("original_decoded_text" in p and "token_ids" in p for p in prediction["predictions"])


def test_payload_multiset_cache_ignores_conditions_but_not_molecular_reagents(
    assets, route, strict_predictor
):
    evaluator = ComputationalRouteSupport(strict_predictor)
    context = {"conditions": {"temperature_C": 80}, "source": "not transmitted"}
    first = evaluator._raw_prediction(["CC(=O)O", "Nc1ccncc1"], context, None)
    second = evaluator._raw_prediction(["Nc1ccncc1", "CC(=O)O"], {"conditions": {"time_h": 5}}, None)
    assert first == second and len(strict_predictor.calls) == 1
    evaluator._raw_prediction(["CC(=O)O", "Nc1ccncc1"], {"reagents": ["O"]}, None)
    assert len(strict_predictor.calls) == 2
    assert set(strict_predictor.calls[0]) == {
        "reactants",
        "reagents",
        "prompt",
        "input_token_ids",
        "input_tokens",
    }
    duplicate = strict_predictor.prepare(["CC", "CC"], [])
    assert duplicate["reactants"] == ["CC", "CC"]


@pytest.mark.parametrize(
    "limits", [{"max_predictor_calls": 0}, {"max_input_tokens": 0}, {"max_output_tokens": 767}]
)
def test_global_budgets_fail_closed_before_call(assets, route, strict_predictor, limits):
    evaluator = ComputationalRouteSupport(strict_predictor, **limits)
    evidence = evaluator.evaluate_route(route, *assets)
    assert not evidence["eligibility"]["passed"] and not strict_predictor.calls


def test_deadline_and_no_predictor_do_not_claim_evidence(assets, route, strict_predictor):
    evaluator = ComputationalRouteSupport(strict_predictor)
    evidence = evaluator.evaluate_route(route, *assets, deadline_monotonic=time.monotonic() - 1)
    assert not evidence["eligibility"]["passed"] and not strict_predictor.calls
    assert not ComputationalRouteSupport().evaluate_route(route, *assets)["eligibility"]["passed"]


@pytest.mark.parametrize(
    "patch",
    [
        {"predictor": None},
        {"policy": {"require_forward_agreement": False}},
        {"policy": {"require_forward_agreement": True, "require_observed_success": True}},
        {"product_matching": "component"},
        {"predictor": {"kind": "reactiont5", "model_path": "/missing", "top_k": 4}},
    ],
)
def test_strict_config_cannot_silently_downgrade(patch):
    values = {"profile": PROFILE, "predictor": {"kind": "reactiont5", "model_path": "/missing"}}
    values.update(patch)
    with pytest.raises(ValueError):
        validate_evidence_config(values)


def test_missing_assets_fail_without_loading_any_model(tmp_path):
    with pytest.raises(ValueError, match="directory is missing"):
        PinnedReactionT5Forward(str(tmp_path / "missing"))
    with pytest.raises(ValueError, match="asset mismatch"):
        PinnedReactionT5Forward(str(tmp_path))


@pytest.mark.parametrize(
    "text", ["CC.", ".CC", "C..C", "CC name", "CC\n", "CC |atomProp:0.foo.bar|", "<unk>"]
)
def test_strict_parser_no_name_cxsmiles_or_empty_component_repair(text):
    with pytest.raises(ValueError):
        strict_smiles(text)


def test_actual_workflow_strict_support_exports_and_failure(tmp_path, assets, proposal, strict_predictor):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    configured = config(tmp_path, assets, proposal)
    out = tmp_path / "accepted"
    rows = workflow.run(HIT, [DESIGN], configured, out)
    assert len(rows) == 1 and rows[0]["evidence_profile"] == PROFILE
    assert rows[0]["route_support_binding"]
    with (out / "results.tsv").open() as stream:
        table = list(csv.DictReader(stream, delimiter="\t"))
    assert table[0]["route_support_binding"] == rows[0]["route_support_binding"]
    for filename in ("poses_hit.sdf", "poses_design.sdf"):
        pose = next(iter(Chem.SDMolSupplier(str(out / filename))))
        assert pose.GetProp("evidence_profile") == PROFILE
        assert pose.GetProp("route_support_binding") == rows[0]["route_support_binding"]
    strict_predictor.products = ["CC", "CO", "CN"]
    assert not workflow.run(HIT, [DESIGN], configured, tmp_path / "rejected")
    saved = json.loads((tmp_path / "rejected/routes.json").read_text())
    assert saved


def test_named_profile_requires_local_model_and_blocks_legacy(tmp_path, assets, proposal):
    values = config(tmp_path, assets, proposal)
    values["recommendation_profile"] = "qualified-development-v1"
    values.pop("synthesis_evidence")
    with pytest.raises(ValueError, match="explicit local"):
        workflow._config(values)
    values["synthesis_evidence"] = {"profile": "legacy-conditional-v1"}
    with pytest.raises(ValueError, match="requires computational"):
        workflow._config(values)
    # Historical explicit engineering mode remains available and labelled.
    values.pop("recommendation_profile")
    legacy = evaluator_from_config(values["synthesis_evidence"])
    assert legacy.predictor is None


def test_registered_input_markers_allowed_but_unknown_token_rejected(strict_predictor, monkeypatch):
    import numpy as np

    class Tokenizer:
        eos_token_id = 1
        all_special_ids = [0, 1, 2, 266, 267]
        unknown = False

        def __call__(self, prompt, **kwargs):
            assert prompt == "REACTANT:CCREAGENT:"
            return {"input_ids": np.array([[266, 2 if self.unknown else 40, 267, 1]])}

    strict_predictor.tokenizer = Tokenizer()
    strict_predictor.max_input_tokens = 512
    prepared = PinnedReactionT5Forward.prepare(strict_predictor, ["CC"], [])
    assert prepared["input_token_ids"] == [266, 40, 267, 1]
    strict_predictor.tokenizer.unknown = True
    with pytest.raises(ValueError, match="unknown"):
        PinnedReactionT5Forward.prepare(strict_predictor, ["CC"], [])


def test_failed_calls_charge_reserved_tokens_and_do_not_retry(assets, route, strict_predictor, monkeypatch):
    def fail(prepared):
        raise RuntimeError("controlled generation failure")

    monkeypatch.setattr(strict_predictor, "predict_raw", fail)
    evaluator = ComputationalRouteSupport(strict_predictor)
    first = evaluator.evaluate_route(route, *assets)
    assert not first["eligibility"]["passed"]
    assert evaluator.usage["predictor_calls"] == 1
    assert evaluator.usage["output_tokens"] == 768
    evaluator.evaluate_route(route, *assets)
    assert evaluator.usage["predictor_calls"] == 1
    with pytest.raises(ValueError, match="finite"):
        evaluator.evaluate_route(route, *assets, deadline_monotonic=float("nan"))
    with pytest.raises(ValueError, match="actual route"):
        evaluator.evaluate(first["route_audit"]["verification"])


def test_multistep_requires_each_exact_intermediate_and_shared_budget(
    assets, strict_predictor, concrete_origins
):
    intermediate = "CC(=O)NCC(=O)O"
    product = "CC(=O)NCC(=O)Nc1ccncc1"
    route = {
        "schema_version": 1,
        "root": "r1",
        "nodes": [
            {"id": "a", "building_block_id": "acetic"},
            {"id": "g", "building_block_id": "glycine"},
            {"id": "n", "building_block_id": "pyridyl"},
            {"id": "r0", "rule_id": "amide", "reactants": ["a", "g"], "product": intermediate},
            {"id": "r1", "rule_id": "amide", "reactants": ["r0", "n"], "product": product},
        ],
    }
    strict_predictor.products = lambda prepared: [
        product if intermediate in prepared["reactants"] else intermediate,
        "CC",
        "CO",
    ]
    evaluator = ComputationalRouteSupport(strict_predictor)
    report = evaluator.evaluate_route(route, *assets)
    assert report["eligibility"]["passed"] and evaluator.usage["predictor_calls"] == 2
    assert (
        report["route_audit"]["steps"][1]["reactant_slots"][0]["parent_binding_sha256"]
        == report["route_audit"]["steps"][0]["step_binding_sha256"]
    )
    assert evaluator.validate_route(report, route, *assets)
    limited = ComputationalRouteSupport(strict_predictor, max_predictor_calls=1).evaluate_route(
        route, *assets
    )
    assert not limited["eligibility"]["passed"]
    assert limited["steps"][0]["forward_prediction"]["status"] == "matched"
    assert limited["steps"][1]["forward_prediction"]["status"] == "unknown"


def test_final_boundary_rejects_forged_report_after_scoring(
    tmp_path, assets, proposal, strict_predictor, monkeypatch
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    original = ComputationalRouteSupport.evaluate_route

    def corrupt(self, *args, **kwargs):
        report = original(self, *args, **kwargs)
        assert report["eligibility"]["passed"]
        report["steps"][0]["forward_prediction"]["raw"]["predictions"][0]["raw_product"] = "CC"
        return report

    monkeypatch.setattr(ComputationalRouteSupport, "evaluate_route", corrupt)
    out = tmp_path / "forged"
    assert workflow.run(HIT, [DESIGN], config(tmp_path, assets, proposal), out) == []
    record = json.loads((out / "routes.json").read_text())["routes"][0]
    assert record["selection_status"] == "qualification_rejected"
    assert "fresh route/input/raw-prediction" in " ".join(record["qualification"]["reasons"])


def test_alternate_route_rescues_without_borrowing_and_continuation_rechecks(
    tmp_path, assets, proposal, strict_predictor
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from chembounce_synthesis.continuation import continuation_inputs

    catalog, rules = assets
    with catalog.path.open("a") as stream:
        stream.write("acetic-anhydride\tCC(=O)OC(C)=O\ttest fixture\t2026-10-15\n")
    data = json.loads(rules.path.read_text())
    data["rules"].append(
        {
            "id": "anhydride",
            "smarts": "[C:1](=[O:2])[O:3][C:4](=[O:5])[C:6].[N;H1,H2:7]>>[C:1](=[O:2])[N:7]",
            "source": "controlled alternative donor",
        }
    )
    rules.path.write_text(json.dumps(data))
    configured = config(tmp_path, assets, proposal)
    alternative = {
        "reactions": [
            {"reaction_template": "anhydride", "reactants": ["CC(=O)OC(C)=O", "Nc1ccncc1"], "product": DESIGN}
        ],
        "building_blocks": ["CC(=O)OC(C)=O", "Nc1ccncc1"],
    }
    from pathlib import Path

    path = Path(configured["provider"]["path"])
    data = json.loads(path.read_text())
    data["proposals"].append({"design": DESIGN, "response": alternative})
    path.write_text(json.dumps(data))
    strict_predictor.products = lambda prepared: [
        DESIGN if "CC(=O)OC(C)=O" in prepared["reactants"] else "CC",
        "CO",
        "CN",
    ]
    out = tmp_path / "alternative"
    rows = workflow.run(HIT, [DESIGN], configured, out)
    assert len(rows) == 1
    records = json.loads((out / "routes.json").read_text())["routes"]
    selected = next(r for r in records if r["route_id"] == rows[0]["route_id"])
    assert selected["route"]["nodes"][-1]["rule_id"] == "anhydride"
    rejected = [r for r in records if not r["synthesis_evidence"]["eligibility"]["passed"]]
    assert rejected and all(r["route_id"] not in rows[0]["alternative_route_ids"] for r in rejected)
    assert len(strict_predictor.calls) == 2
    hit, designs, carried = continuation_inputs(out)
    strict_predictor.products = ["CC", "CO", "CN"]
    assert workflow.run(hit, designs, carried, tmp_path / "resumed") == []
    assert len(strict_predictor.calls) == 4


def test_incomplete_pad_tokens_consume_full_output_budget(assets, route, strict_predictor):
    strict_predictor.products = [[0] + [0] * 256 for _ in range(3)]
    evaluator = ComputationalRouteSupport(strict_predictor, max_output_tokens=768)
    report = evaluator.evaluate_route(route, *assets)
    prediction = report["steps"][0]["forward_prediction"]
    assert prediction["status"] == "partially_evaluated"
    assert not report["eligibility"]["passed"]
    assert prediction["raw"]["output_tokens"] == evaluator.usage["output_tokens"] == 768
    assert all(not p["complete"] and len(p["content_token_ids"]) == 256 for p in prediction["predictions"])
    second = evaluator._raw_prediction(["CCC(=O)O", "Nc1ccncc1"], {}, None)
    assert second["status"] == "unknown" and "budget" in second["reason"]
    assert evaluator.usage["predictor_calls"] == len(strict_predictor.calls) == 1
    assert evaluator.usage["output_tokens"] == 768
