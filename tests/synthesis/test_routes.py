import copy
import json

import pytest

from chembounce_synthesis.common import canonical
from chembounce_synthesis.routes import (
    SearchLimits,
    expand_route,
    parse_response,
    resolve_proposal,
    verify_route,
)


def test_replay_requires_actual_blocks_and_exact_product(assets, route):
    catalog, rules = assets
    evidence = verify_route(route, catalog, rules)
    assert evidence["product"] == canonical("CC(=O)Nc1ccncc1")
    assert evidence["experimental_validation"] is False
    assert evidence["building_blocks"][0]["catalog_sha256"] == catalog.sha256
    forged = copy.deepcopy(route)
    forged["nodes"][-1]["product"] = "CCC(=O)Nc1ccncc1"
    with pytest.raises(ValueError, match="regenerate"):
        verify_route(forged, catalog, rules)
    forged["nodes"][0]["building_block_id"] = "LLM invented SKU"
    with pytest.raises(ValueError, match="absent"):
        verify_route(forged, catalog, rules)


def test_unknown_rules_are_not_fuzzily_repaired(assets, proposal):
    proposal["reactions"][0]["reaction_template"] += "X"
    with pytest.raises(ValueError, match="allowlist"):
        resolve_proposal(proposal, *assets)


def test_resolver_reorders_slots_and_reconstructs_inexact_claim(assets, proposal):
    limits = SearchLimits(building_block_top_k=2, beam_width=8)
    routes, report = resolve_proposal(proposal, *assets, limits)
    assert routes and report.combinations > 0
    assert any(r["exact_model_product"] for r in routes)
    proposal["reactions"][0]["product"] = "CC(=O)NCCc1ccncc1"
    routes, _ = resolve_proposal(proposal, *assets, limits)
    assert routes and all(not r["exact_model_product"] for r in routes)
    for item in routes:
        assert item["verification"]["product"] == verify_route(item["route"], *assets)["product"]


def test_multistep_dag_preserves_actual_intermediates(assets):
    proposal = {
        "reactions": [
            {
                "reaction_template": "amide",
                "reactants": ["Nc1ccncc1", "CC(=O)NCC(=O)O"],
                "product": "CC(=O)NCC(=O)Nc1ccncc1",
            },
            {"reaction_template": "amide", "reactants": ["NCC(=O)O", "CC(=O)O"], "product": "CC(=O)NCC(=O)O"},
        ],
        "building_blocks": ["CC(=O)O", "NCC(=O)O", "Nc1ccncc1"],
    }
    routes, _ = resolve_proposal(proposal, *assets, SearchLimits(building_block_top_k=1, beam_width=8))
    assert routes
    assert any(item["exact_model_product"] for item in routes)
    for item in routes:
        assert len(verify_route(item["route"], *assets)["steps"]) == 2


@pytest.mark.parametrize(
    "change,match",
    [
        (lambda r: r["nodes"].append({"id": "unused", "building_block_id": "ethylamine"}), "unused"),
        (lambda r: r["nodes"][-1].update(reactants=["r0", "b1"]), "Cyclic"),
        (lambda r: r["nodes"][-1].update(reactants=["missing", "b1"]), "Missing"),
        (lambda r: r["nodes"].append(r["nodes"][0].copy()), "Duplicate"),
        (lambda r: r["nodes"][0].update(product="C"), "Unexpected"),
    ],
)
def test_invalid_topology_is_rejected(assets, route, change, match):
    change(route)
    with pytest.raises(ValueError, match=match):
        verify_route(route, *assets)


def test_stereoisomer_mismatch_not_hidden_by_fingerprint(assets, route):
    route["nodes"][1]["building_block_id"] = "alanine_s"
    products, _ = assets[1].apply("amide", ("CC(=O)O", assets[0].get("alanine_s").smiles))
    route["nodes"][-1]["product"] = products[0]
    verify_route(route, *assets)
    route["nodes"][-1]["product"] = (
        products[0].replace("@@", "@") if "@@" in products[0] else products[0].replace("@", "@@")
    )
    with pytest.raises(ValueError, match="stereochemistry"):
        verify_route(route, *assets)


def test_expansion_replays_and_respects_locked_leaves(assets, route):
    items, report = expand_route(
        route, *assets, SearchLimits(building_block_top_k=8), rounds=2, max_routes=20, locked_nodes=["b1"]
    )
    assert items and report.combinations
    assert any(i["verification"]["product"] == canonical("CCC(=O)Nc1ccncc1") for i in items)
    for item in items:
        assert next(n for n in item["route"]["nodes"] if n["id"] == "b1")["building_block_id"] == "pyridyl"
        verify_route(item["route"], *assets)
    assert expand_route(route, *assets, locked_nodes=["b0", "b1"])[0] == []


def test_caps_are_reported(assets, proposal):
    _, report = resolve_proposal(proposal, *assets, SearchLimits(max_combinations=1, building_block_top_k=8))
    assert report.combination_truncations > 0


def test_model_schema_never_executes_code(proposal):
    assert parse_response("```json\n" + json.dumps(proposal) + "\n```") == proposal
    with pytest.raises(ValueError):
        parse_response('__import__("os").system("false")')
    with pytest.raises(ValueError):
        parse_response({"reactions": "not a list", "building_blocks": []})


def test_catalog_locks_do_not_depend_on_route_node_names(assets, route):
    items, _ = expand_route(
        route, *assets, SearchLimits(building_block_top_k=8), locked_building_block_ids=["pyridyl"]
    )
    assert items
    for item in items:
        assert next(n for n in item["route"]["nodes"] if n["id"] == "b1")["building_block_id"] == "pyridyl"


def test_expansion_reserves_later_rounds_and_records_substitution_similarity(assets, route):
    from rdkit import DataStructs
    from chembounce_synthesis.common import FP, molecule

    items, report = expand_route(
        route,
        *assets,
        SearchLimits(building_block_top_k=1, max_combinations=100),
        building_block_top_k=4,
        rounds=2,
        max_routes=12,
    )
    assert len(items) <= 12
    assert report.rounds_completed == 2
    assert all(r["retained_routes"] for r in report.round_reports)
    assert any({s["node_id"] for s in item["substitutions"]} == {"b0", "b1"} for item in items), (
        "A two-round expansion must be able to replace both reactants"
    )
    for item in items:
        verify_route(item["route"], *assets)
        for substitution in item["substitutions"]:
            expected = DataStructs.TanimotoSimilarity(
                FP.GetFingerprint(molecule(substitution["from_smiles"])),
                FP.GetFingerprint(molecule(substitution["to_smiles"])),
            )
            assert substitution["tanimoto"] == pytest.approx(expected)


def test_expansion_neighbour_count_and_optional_tanimoto_cutoff(assets, route):
    small, _ = expand_route(route, *assets, building_block_top_k=1, locked_nodes=["b1"])
    large, _ = expand_route(route, *assets, building_block_top_k=4, locked_nodes=["b1"])
    assert len(large) > len(small)
    assert expand_route(route, *assets, min_tanimoto=1.0)[0] == []


def test_catalog_all_reaction_roles_are_applied_before_shortlisting(assets):
    blocks = assets[0].nearest("CCN", "[N]", k=1, required_roles=("C(=O)O",))
    assert blocks
    assert blocks[0].id in {"glycine", "alanine_s", "alanine_r"}


def test_implicit_reagent_atoms_are_exposed_in_evidence(assets, tmp_path):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "implicit.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "implicit-methylation",
                        "smarts": "[N;H2:1]>>[N:1]C",
                        "source": "fixture with an implicit carbon source",
                    }
                ],
            }
        )
    )
    rules = Rules(path)
    catalog = assets[0]
    products, _ = rules.apply("implicit-methylation", (catalog.get("ethylamine").smiles,))
    route = {
        "schema_version": 1,
        "root": "r",
        "nodes": [
            {"id": "b", "building_block_id": "ethylamine"},
            {"id": "r", "rule_id": "implicit-methylation", "reactants": ["b"], "product": products[0]},
        ],
    }
    evidence = verify_route(route, catalog, rules)
    assert evidence["has_implicit_atom_sources"] is True
    assert evidence["reagents_and_conditions_validated"] is False
    assert evidence["steps"][0]["unmapped_product_template_atoms"][0]["element"] == "C"


def test_alkene_stereochemistry_survives_reaction_replay(assets):
    from chembounce_synthesis.catalog import Catalog

    path = assets[0].path
    with path.open("a") as stream:
        stream.write("crotonic-E\tC/C=C/C(=O)O\ttest fixture\tfixture-date\n")
    catalog = Catalog(path)
    rules = assets[1]
    reactants = (catalog.get("crotonic-E").smiles, catalog.get("ethylamine").smiles)
    products, _ = rules.apply("amide", reactants)
    assert products
    expected = canonical("C/C=C/C(=O)NCC")
    assert expected in products
    route = {
        "schema_version": 1,
        "root": "r",
        "nodes": [
            {"id": "a", "building_block_id": "crotonic-E"},
            {"id": "b", "building_block_id": "ethylamine"},
            {"id": "r", "rule_id": "amide", "reactants": ["a", "b"], "product": "C/C=C\\C(=O)NCC"},
        ],
    }
    with pytest.raises(ValueError, match="stereochemistry"):
        verify_route(route, catalog, rules)
