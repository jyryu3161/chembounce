"""Root-only route guidance under unchanged configured search budgets.

The 0.5 Tanimoto boundary is an engineering fixture, not a calibrated cutoff.
"""

import copy
import json

import pytest

from chembounce_synthesis import workflow
from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.cli import main
from chembounce_synthesis.common import canonical
from chembounce_synthesis.identity import classify_change
from chembounce_synthesis.route_objective import RootProductObjective, search_policy
from chembounce_synthesis.routes import SearchLimits, expand_route, resolve_proposal, verify_route

HIT = "CC(=O)Nc1ccccc1"
GOOD = "CC(=O)Nc1ccncc1"
DESIGN = "CCC(=O)Nc1ccncc1"


def narrow_catalog(assets):
    path = assets[0].path
    lines = path.read_text().splitlines()
    keep = {"acetic", "propionic", "pyridyl", "phenyl"}
    path.write_text(
        "\n".join([lines[0], *(line for line in lines[1:] if line.split("\t")[0] in keep)]) + "\n"
    )
    return Catalog(path), assets[1]


def proposal_for(product=DESIGN, acid="CCC(=O)O", amine="Nc1ccncc1"):
    return {
        "reactions": [{"reaction_template": "amide", "reactants": [acid, amine], "product": product}],
        "building_blocks": [acid, amine],
    }


def objective(hit=HIT, design=DESIGN, minimum=0.5, *, rounds=0, cache_size=4096, core=None):
    config = {
        "recommendation": {"mode": "qualified"},
        "expansion": {"rounds": rounds},
        "route_search_cache_size": cache_size,
    }
    return RootProductObjective(
        hit,
        design,
        {"hit_tanimoto": minimum},
        lambda product: workflow.check_constraints(product, hit, design, core_smiles=core),
        search_policy(config),
    )


def config_for(tmp_path, assets, proposal, **extra):
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
        "minimum_scores": {"hit_tanimoto": 0.5, "hit_shape": 0.01, "hit_pharmacophore": 0.01},
        "conformers": {"num_conformers": 2},
        "search": {"building_block_top_k": 2, "beam_width": 1, "max_combinations": 4},
        "max_scored_products": 1,
        "max_results": 1,
        **extra,
    }


def read(out, name):
    return json.loads((out / name).read_text())


def test_beam_one_retains_hit_similar_hop_under_equal_actual_combination_count(assets):
    assets = narrow_catalog(assets)
    limits = SearchLimits(building_block_top_k=2, beam_width=1, max_combinations=4)
    legacy, before = resolve_proposal(proposal_for(), *assets, limits)
    guided, after = resolve_proposal(proposal_for(), *assets, limits, root_objective=objective())
    assert legacy[0]["verification"]["product"] == DESIGN
    assert guided[0]["verification"]["product"] == GOOD
    assert before.combinations == after.combinations == 4
    assert len(legacy) == len(guided) == 1
    assert guided[0]["search_eligibility"]["hit_tanimoto"] == pytest.approx(0.6)
    assert guided[0]["search_eligibility"]["design_tanimoto"] == pytest.approx(0.6071428571428571)
    assert verify_route(guided[0]["route"], *assets)["product"] == GOOD
    assert after.root_ineligible_state_visits and after.root_objective_cache_hits


def test_identical_hit_cannot_crowd_out_a_real_hop(assets):
    assets = narrow_catalog(assets)
    limits = SearchLimits(building_block_top_k=2, beam_width=1, max_combinations=4)
    roots, _ = resolve_proposal(
        proposal_for(HIT, "CC(=O)O", "Nc1ccccc1"), *assets, limits, root_objective=objective(design=HIT)
    )
    assert len(roots) == 1 and roots[0]["verification"]["product"] == GOOD
    assert classify_change(HIT, roots[0]["verification"]["product"]) == "scaffold"


def test_model_claim_is_not_substituted_for_actual_design_reference(assets):
    assets = narrow_catalog(assets)
    obj = objective(design=GOOD)
    roots, _ = resolve_proposal(
        proposal_for(), *assets, SearchLimits(building_block_top_k=2, beam_width=1), root_objective=obj
    )
    assert roots[0]["model_claimed_product"] == DESIGN
    assert roots[0]["search_eligibility"]["design_tanimoto"] == 1.0


def test_low_hit_similarity_intermediate_and_incomplete_core_are_not_globally_filtered(assets):
    hit = "CC(=O)NCC(=O)Nc1ccccc1"
    final = "CC(=O)NCC(=O)Nc1ccncc1"
    intermediate = "CC(=O)NCC(=O)O"
    proposal = {
        "reactions": [
            {"reaction_template": "amide", "reactants": ["CC(=O)O", "NCC(=O)O"], "product": intermediate},
            {"reaction_template": "amide", "reactants": [intermediate, "Nc1ccncc1"], "product": final},
        ],
        "building_blocks": ["CC(=O)O", "NCC(=O)O", "Nc1ccncc1"],
    }
    obj = objective(hit, final, core="NCC(=O)N")
    roots, report = resolve_proposal(
        proposal,
        *assets,
        SearchLimits(building_block_top_k=1, beam_width=1, max_combinations=16),
        root_objective=obj,
    )
    assert roots and roots[0]["verification"]["product"] == final
    assert all(item["product"] != intermediate for item in obj.snapshot()["decision_sample"])
    assert len(verify_route(roots[0]["route"], *assets)["steps"]) == 2
    assert report.combinations <= 16


def test_expansion_rebuild_and_frontier_prefer_original_hit(assets, route):
    assets = narrow_catalog(assets)
    route["nodes"][0]["building_block_id"] = "propionic"
    route["nodes"][-1]["product"] = DESIGN
    obj = objective(rounds=1)
    items, report = expand_route(
        route,
        *assets,
        SearchLimits(beam_width=1, max_combinations=2),
        rounds=1,
        max_routes=1,
        building_block_top_k=1,
        target_smiles=DESIGN,
        root_objective=obj,
    )
    assert len(items) == 1 and items[0]["verification"]["product"] == GOOD
    assert report.combinations <= 2 and report.root_prune_calls >= 2
    assert items[0]["search_eligibility"]["eligible"]


def test_multiround_seed_only_recovery_does_not_spend_eligible_output_slot(assets, route):
    assets = narrow_catalog(assets)
    initial = "CCC(=O)Nc1ccccc1"
    route["nodes"][0]["building_block_id"] = "propionic"
    route["nodes"][1]["building_block_id"] = "phenyl"
    route["nodes"][-1]["product"] = initial
    obj = objective(rounds=2)
    items, report = expand_route(
        route,
        *assets,
        SearchLimits(beam_width=1, max_combinations=4),
        rounds=2,
        max_routes=1,
        building_block_top_k=1,
        target_smiles=DESIGN,
        root_objective=obj,
    )
    assert len(items) == 1 and items[0]["verification"]["product"] == GOOD
    assert items[0]["expansion_round"] == 2
    assert report.round_reports[0]["eligible_routes"] == 0
    assert report.round_reports[0]["seed_only_frontier_routes"] == 1
    assert report.round_reports[1]["eligible_routes"] == 1
    assert report.combinations == 4 and report.rounds_completed == 2
    assert all(item["search_eligibility"]["eligible"] for item in items)


def test_same_hit_seed_recovers_in_workflow_without_evidence_or_3d_for_seed(tmp_path, assets):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    assets = narrow_catalog(assets)
    config = config_for(
        tmp_path,
        assets,
        proposal_for(HIT, "CC(=O)O", "Nc1ccccc1"),
        search={"building_block_top_k": 1, "beam_width": 1, "max_combinations": 2},
        expansion={"rounds": 1, "max_routes": 1, "building_block_top_k": 1},
    )
    out = tmp_path / "seed"
    rows = workflow.run(HIT, [DESIGN], config, out)
    assert rows and rows[0]["smiles"] == GOOD
    records = read(out, "routes.json")["routes"]
    seed = next(r for r in records if r["verification"]["product"] == HIT)
    assert seed["selection_status"] == "search_seed_only" and not seed["eligible_for_scoring"]
    assert "synthesis_evidence" not in seed
    assert read(out, "summary.json")["scored_products"] == 1


def test_cache_is_bounded_scoped_and_explicitly_sampled(assets):
    assets = narrow_catalog(assets)
    obj = objective(cache_size=1)
    limits = SearchLimits(building_block_top_k=2, beam_width=1, max_combinations=4)
    first, _ = resolve_proposal(proposal_for(), *assets, limits, root_objective=obj)
    second, _ = resolve_proposal(proposal_for(), *assets, limits, root_objective=obj)
    assert first == second
    audit = obj.snapshot()
    assert audit["current_cached_unique_products"] == len(audit["decision_sample"]) == 1
    assert audit["counts"]["evictions"] > 0
    assert audit["counts"]["evaluations"] + audit["counts"]["cache_hits"] == audit["counts"]["requests"]
    assert audit["costs"]["constraint_checks"] == audit["counts"]["evaluations"]
    assert objective(design=GOOD).assess(GOOD)["design_tanimoto"] == 1.0
    assert obj.assess(GOOD)["design_tanimoto"] < 1.0


def test_explicit_limits_and_deterministic_replay_are_unchanged(assets):
    assets = narrow_catalog(assets)
    for budget in (1, 2, 4):
        limits = SearchLimits(building_block_top_k=2, beam_width=1, max_combinations=budget)
        left, a = resolve_proposal(proposal_for(), *assets, limits, root_objective=objective())
        right, b = resolve_proposal(
            copy.deepcopy(proposal_for()), *assets, limits, root_objective=objective()
        )
        assert left == right and a == b and a.combinations <= budget and len(left) <= 1


def test_streaming_root_prune_and_exact_stereo_replay(assets):
    hit = "N[C@@H](C)C(=O)Nc1ccccc1"
    design = "N[C@@H](C)C(=O)Nc1ccncc1"
    roots, report = resolve_proposal(
        proposal_for(design, "N[C@@H](C)C(=O)O"),
        *assets,
        SearchLimits(building_block_top_k=4, beam_width=1, max_combinations=32),
        root_objective=objective(hit, design, minimum=0.01),
    )
    assert report.root_prune_calls > 1 and report.combinations <= 32
    assert roots and roots[0]["verification"]["product"] == canonical(design)
    assert verify_route(roots[0]["route"], *assets)["product"] == canonical(design)
    wrong = copy.deepcopy(roots[0]["route"])
    next(node for node in wrong["nodes"] if node["id"] == wrong["root"])["product"] = (
        "N[C@H](C)C(=O)Nc1ccncc1"
    )
    with pytest.raises(ValueError, match="stereochemistry"):
        verify_route(wrong, *assets)


def test_qualified_cli_defaults_to_hit_search_and_design_ablation_keeps_final_gate(tmp_path, assets):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    assets = narrow_catalog(assets)
    config = config_for(tmp_path, assets, proposal_for())
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    out = tmp_path / "qualified"
    args = ["--hit", HIT, "--design", DESIGN, "--config", str(path), "--output", str(out)]
    assert main(args) == 0
    assert read(out, "results.json")["results"][0]["smiles"] == GOOD
    assert read(out, "manifest.json")["route_search_policy"]["version"] == "original-hit-v1"
    audit = read(out, "search_objectives.json")["design_objectives"][0]
    assert audit["counts"]["rejected_evaluations"] > 0
    out = tmp_path / "ablation"
    args[-1] = str(out)
    assert main(args + ["--route-search-policy", "design-v1"]) == 3
    assert read(out, "manifest.json")["route_search_policy"]["version"] == "design-v1"
    assert read(out, "results.json")["policy"]["minimum_scores"]["hit_tanimoto"] == 0.5
    assert read(out, "summary.json")["scored_products"] == 0
    policy = read(out, "search_objectives.json")
    assert not policy["design_objectives"]
    assert policy["policy"]["cache_size_per_design"] == 0
    assert "Legacy" in policy["policy"]["ranking"]


@pytest.mark.parametrize("requirement", ["similarity", "anchor", "property"])
def test_default_workflow_records_early_root_failures_without_scoring(tmp_path, assets, requirement):
    assets = narrow_catalog(assets)
    config = config_for(tmp_path, assets, proposal_for())
    if requirement == "similarity":
        config["minimum_scores"]["hit_tanimoto"] = 0.99
    elif requirement == "anchor":
        config["core_smiles"] = "Nc1ccccc1"
    else:
        config["property_thresholds"] = {"mw": [1, 2]}
    out = tmp_path / requirement
    assert workflow.run(HIT, [DESIGN], config, out) == []
    summary = read(out, "summary.json")
    assert summary["scored_products"] == 0
    assert not read(out, "routes.json")["routes"]
    audit = read(out, "search_objectives.json")["design_objectives"][0]
    assert audit["counts"]["rejected_evaluations"] > 0
    assert all(not item["eligible"] and item["reasons"] for item in audit["decision_sample"])
    assert summary["search_reports"][0]["root_filtered_state_visits"] > 0


def test_qualified_root_rejects_unspecified_product_before_scoring(tmp_path, assets):
    with assets[0].path.open("a") as stream:
        stream.write("lactic-unknown\tCC(O)C(=O)O\ttest fixture\tfixture-date\n")
    config = config_for(
        tmp_path,
        assets,
        proposal_for("CC(O)C(=O)Nc1ccncc1", "CC(O)C(=O)O"),
        search={"building_block_top_k": 1, "beam_width": 1, "max_combinations": 2},
    )
    config["minimum_scores"]["hit_tanimoto"] = 0.01
    out = tmp_path / "stereo"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    audit = read(out, "search_objectives.json")["design_objectives"][0]
    assert any(
        "specified stereochemistry" in reason
        for item in audit["decision_sample"]
        for reason in item["reasons"]
    )
    assert read(out, "summary.json")["scored_products"] == 0
