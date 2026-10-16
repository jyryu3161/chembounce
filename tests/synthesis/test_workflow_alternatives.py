"""Route/design fallback and fair allocation regressions."""

import json

import pytest

pytest.importorskip("rdkit.Chem.rdShapeAlign")

from chembounce_synthesis import workflow
from chembounce_synthesis.workflow import _balanced_order, _candidate_order, _pareto, run

HIT = "CC(=O)Nc1ccccc1"
DESIGN = "CC(=O)Nc1ccncc1"
OTHER = "CCC(=O)Nc1ccncc1"


def fallback_config(tmp_path, assets, proposal):
    path = tmp_path / "proposals.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "proposals": [{"design": design, "response": proposal} for design in (DESIGN, OTHER)],
            }
        )
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "conformers": {"num_conformers": 2},
        "refinement": {"top_n": 1, "num_conformers": 2, "seeds": [42, 43, 44]},
        "max_scored_products": 4,
    }


def test_failed_primary_design_cannot_erase_refinable_alternative(tmp_path, assets, proposal, monkeypatch):
    original_screen, original_refine = workflow.score_candidate, workflow.refine_candidate
    seen = []

    def screen(product, hit, design, *args, **kwargs):
        result = original_screen(product, hit, design, *args, **kwargs)
        result[0]["score"] = 0.95 if design == DESIGN else 0.8
        return result

    def refine(product, hit, design, *args, **kwargs):
        seen.append(design)
        if design == DESIGN:
            return None, {"status": "failed", "trials": [], "reason": "controlled design-specific failure"}
        return original_refine(product, hit, design, *args, **kwargs)

    monkeypatch.setattr(workflow, "score_candidate", screen)
    monkeypatch.setattr(workflow, "refine_candidate", refine)
    out = tmp_path / "out"
    rows = run(HIT, [DESIGN, OTHER], fallback_config(tmp_path, assets, proposal), out)
    assert seen == [DESIGN, OTHER]
    assert len(rows) == 1 and rows[0]["design"] == OTHER
    records = json.loads((out / "routes.json").read_text())["routes"]
    assert [r["selection_status"] for r in records] == ["refinement_failed", "accepted"]
    assert rows[0]["alternative_route_ids"] == [records[0]["route_id"]]
    comparison = rows[0]["hit_identity_comparison"]
    assert comparison["status"] == "recorded"
    assert comparison["original_hit_smiles"] == HIT
    assert comparison["final_product_smiles"] == rows[0]["smiles"] == DESIGN
    assert comparison["final_product"]["exact_smiles"] != rows[0]["design"]
    assert comparison["change_class"] == "scaffold"
    report = json.loads((out / "summary.json").read_text())
    assert report["screened_products"] == 1
    assert report["screened_product_design_pairs"] == 2
    assert report["refinement_attempted"] == 2
    assert len((out / "screening_alternatives.tsv").read_text().splitlines()) == 3


def test_balanced_allocation_visits_scaffolds_and_designs_without_sorting_away_priority():
    items = [
        {"product": "c1ccncc1", "design": "A", "id": "A1"},
        {"product": "Cc1ccncc1", "design": "A", "id": "A2"},
        {"product": "CCc1ccncc1", "design": "B", "id": "B1"},
        {"product": "c1ccoc1", "design": "C", "id": "C1"},
        {"product": "Cc1ccoc1", "design": "C", "id": "C2"},
    ]
    result = _balanced_order(items, product=lambda r: r["product"], design=lambda r: r["design"])
    assert [r["id"] for r in result] == ["A1", "C1", "B1", "C2", "A2"]
    for item in items:
        item["records"] = [{"carried_from": "old"}] if item["id"] == "A2" else [{}]
    assert _candidate_order(items, "balanced")[0]["id"] == "A2"


def test_robust_pareto_cannot_reward_lucky_representative_metrics():
    common = {"hit_tanimoto": 0.5, "design_tanimoto": 0.5, "relative_energy": 0, "reaction_steps": 1}
    unlucky = {
        **common,
        **{key: 0.99 for key in ("hit_shape", "hit_pharmacophore", "design_shape", "design_pharmacophore")},
        "ranking_metrics": {
            key: 0.2 for key in ("hit_shape", "hit_pharmacophore", "design_shape", "design_pharmacophore")
        },
    }
    stable = {
        **common,
        **{key: 0.8 for key in ("hit_shape", "hit_pharmacophore", "design_shape", "design_pharmacophore")},
    }
    _pareto([unlucky, stable])
    assert stable["pareto_front"] == 0 and unlucky["pareto_front"] == 1


def test_global_refinement_budget_protects_later_distinct_scaffold(tmp_path, assets, proposal, monkeypatch):
    import chembounce_synthesis.refinement as refinement_module

    other = "CCNC(C)=O"
    config = fallback_config(tmp_path, assets, proposal)
    config["refinement"] = {
        "top_n": 2,
        "num_conformers": 2,
        "seeds": [42, 43, 44],
        "adaptive_max_conformers": 8,
        "max_total_conformers": 18,
        "max_score_range": 0.05,
    }
    path = tmp_path / "proposals.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "proposals": [
                    {"design": DESIGN, "response": proposal},
                    {
                        "design": other,
                        "response": {
                            "reactions": [
                                {
                                    "reaction_template": "amide",
                                    "reactants": ["CC(=O)O", "CCN"],
                                    "product": other,
                                }
                            ],
                            "building_blocks": ["CC(=O)O", "CCN"],
                        },
                    },
                ],
            }
        )
    )
    original = workflow.score_candidate
    calls = []

    def score(product, hit, design, settings, pose):
        result = original(product, hit, design, settings, pose)
        value = (0.2 if settings.seed == 42 else 0.99) if product == DESIGN else 0.8
        result[0]["score"] = value
        for role in ("hit", "design"):
            for metric in ("shape", "pharmacophore"):
                result[0][role][metric] = value
        calls.append((product, settings.num_conformers, settings.seed))
        return result

    def screen(product, hit, design, *args, **kwargs):
        result = original(product, hit, design, *args, **kwargs)
        value = 0.99 if product == DESIGN else 0.8
        result[0]["score"] = value
        result[0]["candidate_relative_energy"] = 0
        for role in ("hit", "design"):
            for metric in ("shape", "pharmacophore"):
                result[0][role][metric] = value
        return result

    monkeypatch.setattr(workflow, "score_candidate", screen)
    monkeypatch.setattr(refinement_module, "score_candidate", score)
    out = tmp_path / "budget"
    rows = run(HIT, [DESIGN, other], config, out)
    reports = json.loads((out / "refinement.json").read_text())["candidates"]
    assert reports[0]["product"] == DESIGN
    assert reports[0]["status"] == "budget_exhausted"
    assert reports[0]["global_budget"]["reserved_for_others"] == 6
    assert reports[1]["status"] == "success"
    assert len(rows) == 1 and rows[0]["smiles"] == other
    assert [call[1] for call in calls] == [2] * 6
    summary = json.loads((out / "summary.json").read_text())
    assert summary["refinement_budget"]["candidate_conformer_requests"] == 12
    assert summary["screened_scaffold_families"] == 2


def test_family_dedup_never_labels_other_stereoisomer_as_same_product_route(tmp_path, assets, proposal):
    first, second = "CC(=O)N[C@@H](C)C(=O)O", "CC(=O)N[C@H](C)C(=O)O"
    config = fallback_config(tmp_path, assets, proposal)
    config["refinement"] = {"top_n": 0}
    path = tmp_path / "proposals.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "proposals": [
                    {
                        "design": product,
                        "response": {
                            "reactions": [
                                {
                                    "reaction_template": "amide",
                                    "reactants": ["CC(=O)O", amine],
                                    "product": product,
                                }
                            ],
                            "building_blocks": ["CC(=O)O", amine],
                        },
                    }
                    for product, amine in ((first, "N[C@@H](C)C(=O)O"), (second, "N[C@H](C)C(=O)O"))
                ],
            }
        )
    )
    out = tmp_path / "stereo"
    rows = run(HIT, [first, second], config, out)
    assert len(rows) == 1
    row = rows[0]
    routes = {r["route_id"]: r for r in json.loads((out / "routes.json").read_text())["routes"]}
    assert all(
        routes[rid]["verification"]["product"] == row["smiles"]
        for rid in [row["route_id"], *row["alternative_route_ids"]]
    )
    assert len(row["alternative_family_candidates"]) == 1
    other = row["alternative_family_candidates"][0]
    assert other["smiles"] != row["smiles"]
    assert routes[other["route_id"]]["verification"]["product"] == other["smiles"]
