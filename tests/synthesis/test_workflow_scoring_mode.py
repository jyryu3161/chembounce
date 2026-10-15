"""Paired scoring ablation: actual routes, evidence queues and allocation agree."""

import copy
import json

import pytest

from chembounce_synthesis import workflow
from chembounce_synthesis.evidence import SynthesisEvidenceEvaluator
from chembounce_synthesis.identity import scaffold_family_key


HIT = "CC(=O)Nc1ccccc1"
DESIGNS = ["CC(=O)Nc1ccncc1", "CC(=O)NC1CCCCC1"]


def inputs(tmp_path, assets):
    catalog, rules = assets
    with catalog.path.open("a") as stream:
        stream.write("cyclohexyl\tNC1CCCCC1\tunit test fixture\t2026-10-16\n")
    proposals = []
    for design, amine in zip(DESIGNS, ("Nc1ccncc1", "NC1CCCCC1")):
        proposal = {
            "design": design,
            "response": {
                "reactions": [
                    {"reaction_template": "amide", "reactants": ["CC(=O)O", amine], "product": design}
                ],
                "building_blocks": ["CC(=O)O", amine],
            },
        }
        # Repeated proposals exercise alternate-route evidence-cache handling.
        proposals.extend([proposal, copy.deepcopy(proposal)])
    path = tmp_path / "proposals.json"
    path.write_text(json.dumps({"schema_version": 1, "proposals": proposals}))
    return {
        "schema_version": 1,
        "catalog": str(catalog.path),
        "rules": str(rules.path),
        "provider": {"kind": "saved", "path": str(path)},
        "search": {"building_block_top_k": 2, "beam_width": 4, "max_combinations": 16},
        "expansion": {"rounds": 1, "max_routes": 6},
        "allocation_policy": "balanced",
        "conformers": {"num_conformers": 2, "max_iterations": 1000},
        "refinement": {"top_n": 0},
        "max_scored_products": 2,
    }


def test_same_saved_proposals_share_route_evidence_and_balanced_allocation(tmp_path, assets, monkeypatch):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = inputs(tmp_path, assets)
    traces, current = {}, None

    class FixturePredictor:
        # This test oracle replays a rule from reactants only. It is deliberately
        # not independent scientific evidence and never leaves this unit test.
        provenance = {"provider": "unit-test reaction oracle; not scientific evidence"}

        def token_budget(self, reactants, reagents, conditions):
            return {"input_tokens": 5, "max_output_tokens": 20}

        def predict(self, reactants, reagents, conditions):
            traces[current]["predictor_inputs"].append([list(reactants), list(reagents), dict(conditions)])
            products, _ = assets[1].apply("amide", tuple(reactants))
            return {"products": list(products), "complete": [True] * len(products), "output_tokens": 8}

    class TracedEvaluator(SynthesisEvidenceEvaluator):
        def evaluate(self, verification, **kwargs):
            traces[current]["evidence_order"].append(copy.deepcopy(verification))
            return super().evaluate(verification, **kwargs)

    monkeypatch.setattr(
        workflow,
        "_evidence_evaluator",
        lambda config: TracedEvaluator(
            FixturePredictor(),
            {"require_forward_agreement": True},
            max_predictor_calls=4,
        ),
    )
    results = {}
    for mode in ("3d", "2d"):
        current = mode
        traces[mode] = {"evidence_order": [], "predictor_inputs": []}
        out = tmp_path / mode
        rows = workflow.run(HIT, DESIGNS, config, out, scoring_mode=mode)
        assert rows
        results[mode] = {"rows": rows}
        for filename in ("summary", "routes", "proposals", "manifest"):
            results[mode][filename] = json.loads((out / f"{filename}.json").read_text())

    left, right = results["3d"], results["2d"]
    assert left["proposals"] == right["proposals"]
    assert traces["3d"] == traces["2d"]
    assert len(traces["3d"]["predictor_inputs"]) == 4
    assert len(traces["3d"]["evidence_order"]) > 4
    for key in (
        "search_reports",
        "verified_routes",
        "scored_products",
        "screening_allocation",
        "synthesis_evidence_usage",
    ):
        assert left["summary"][key] == right["summary"][key]
    assert left["summary"]["synthesis_evidence_usage"]["cache_hits"] > 0
    assert left["summary"]["scored_products"] == 2
    allocation = left["summary"]["screening_allocation"]
    assert len({scaffold_family_key(row["product"]) for row in allocation}) == 2
    compare_keys = (
        "route_id",
        "design",
        "route",
        "eligible_for_evidence",
        "eligible_for_scoring",
        "synthesis_evidence",
    )
    assert [{k: row.get(k) for k in compare_keys} for row in left["routes"]["routes"]] == [
        {k: row.get(k) for k in compare_keys} for row in right["routes"]["routes"]
    ]
    assert (tmp_path / "3d/poses_hit.sdf").is_file()
    assert not list((tmp_path / "2d").glob("*.sdf"))
    assert all(row["hit_shape"] is not None for row in left["rows"])
    assert all(row["hit_shape"] is None and row["seed_count"] == 0 for row in right["rows"])


def test_2d_has_no_alignment_dependency_or_conformer_work(tmp_path, assets, monkeypatch):
    config = inputs(tmp_path, assets)
    config["refinement"] = {"top_n": 2}

    def forbidden(*args, **kwargs):
        raise AssertionError("2D ablation invoked a geometric operation")

    for name in (
        "prepare_reference",
        "score_candidate",
        "ensemble",
        "reference_from_sdf",
        "refine_candidate",
    ):
        monkeypatch.setattr(workflow, name, forbidden)
    rows = workflow.run(HIT, DESIGNS, config, tmp_path / "2d", scoring_mode="2d")
    assert rows and all(row["refinement_conformers"] == 0 for row in rows)
    manifest = json.loads((tmp_path / "2d/manifest.json").read_text())
    assert manifest["scoring_backend"]["requested_refinement"]["top_n"] == 2
    assert manifest["scoring_backend"]["effective_refinement"]["top_n"] == 0


def test_2d_rejects_geometric_thresholds_and_pocket_constraints(tmp_path, assets):
    config = inputs(tmp_path, assets)
    config["minimum_scores"] = {"hit_shape": 0.5}
    with pytest.raises(ValueError, match="2D scoring cannot enforce"):
        workflow.run(HIT, DESIGNS, config, tmp_path / "threshold", scoring_mode="2d")
    config.pop("minimum_scores")
    # A valid pocket config is checked by the general configuration validator
    # before the scoring-mode compatibility check; use a minimal valid form.
    config["reference_sdf"] = str(tmp_path / "reference.sdf")
    config["pocket"] = {
        "receptor_pdb": str(tmp_path / "receptor.pdb"),
        "reference_ligand": {"chain": "A", "residue_number": 1, "residue_name": "LIG"},
    }
    with pytest.raises(ValueError, match="2D scoring cannot enforce"):
        workflow.run(HIT, DESIGNS, config, tmp_path / "pocket", scoring_mode="2d")
