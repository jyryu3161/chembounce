"""Concrete deterministic fixtures, not empirical reaction success records."""

import copy
import json

import pytest
from rdkit import Chem

from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.common import canonical
from chembounce_synthesis.reactions import Rules
from chembounce_synthesis.route_audit import atom_origins, audit_route


def has_concrete_origin_metadata():
    """Probe mapped/off-template atoms and two distinct slots, not a version string."""
    from rdkit.Chem import rdChemReactions

    reaction = rdChemReactions.ReactionFromSmarts("[C:1].[C:2]>>[C:1][C:2]")
    product = reaction.RunReactants((Chem.MolFromSmiles("CC"), Chem.MolFromSmiles("C")))[0][0]
    atoms = list(product.GetAtoms())
    return (
        all(a.HasProp("react_idx") and a.HasProp("react_atom_idx") for a in atoms)
        and {a.GetIntProp("react_idx") for a in atoms} == {0, 1}
        and any(not a.HasProp("old_mapno") for a in atoms)
    )


@pytest.fixture(scope="session")
def concrete_origin_metadata():
    if not has_concrete_origin_metadata():
        pytest.skip(
            "RDKit lacks concrete slot/atom metadata; positive witness checks require that capability"
        )


def fixture(tmp_path, smiles, smarts, *, product=None, repeated=False):
    catalog_path = tmp_path / "audit-catalog.tsv"
    catalog_path.write_text(
        "id\tsmiles\tsource\tsnapshot\n"
        + "".join(f"b{i}\t{s}\tsynthetic test\t2026-10-16\n" for i, s in enumerate(smiles))
    )
    rules_path = tmp_path / "audit-rules.json"
    rules_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {"id": "r", "smarts": smarts, "source": "synthetic fixture; not validated chemistry"}
                ],
            }
        )
    )
    catalog, rules = Catalog(catalog_path), Rules(rules_path)
    reactants = [f"b{i}" for i in range(len(smiles))]
    if repeated:
        reactants = ["b0", "b0"]
    products, _ = rules.apply("r", tuple(catalog.get(i).smiles for i in reactants))
    assert products, "Fixture needs a valid principal product"
    route = {
        "schema_version": 1,
        "root": "r",
        "nodes": [
            *[{"id": f"b{i}", "building_block_id": f"b{i}"} for i in range(len(smiles))],
            {"id": "r", "rule_id": "r", "reactants": reactants, "product": product or products[0]},
        ],
    }
    return route, catalog, rules


def audit_fixture(tmp_path, smiles, smarts, **kwargs):
    route, catalog, rules = fixture(tmp_path, smiles, smarts, **kwargs)
    return audit_route(route, catalog, rules)


def changes(result, kind="tetrahedral"):
    return [
        c
        for s in result["steps"]
        for w in s["witnesses"]
        for c in w["stereo"]["changes"]
        if c["type"] == kind
    ]


@pytest.mark.usefixtures("concrete_origin_metadata")
def test_amide_audits_off_template_atoms_and_lost_oxygen(assets, route):
    result = audit_route(route, *assets)
    assert result["strict_eligible"] is True
    assert result["verification"]["experimental_validation"] is False
    step = result["steps"][0]
    atoms = step["witnesses"][0]["atom_origins"]
    assert atoms["status"] == "complete"
    assert any(a["old_map_number"] is None for a in atoms["atoms"])
    assert len(atoms["atoms"]) == Chem.MolFromSmiles(step["product"]).GetNumHeavyAtoms()
    assert {a["element"] for a in atoms["lost_input_atoms"]} == {"O"}
    assert sorted(a["canonical_product_atom_index"] for a in atoms["atoms"]) == list(
        range(len(atoms["atoms"]))
    )
    assert result == audit_route(copy.deepcopy(route), *assets)


def test_missing_template_donor_is_unsupported_even_with_same_element_in_input(tmp_path):
    result = audit_fixture(tmp_path, ["CCO"], "[C:1]>>[C:1]O")
    assert result["strict_eligible"] is False
    assert result["atom_origin_status"] == "unsupported"
    assert result["steps"][0]["template_audit"]["unbound_product_atoms"][0]["element"] == "O"
    assert any(
        p["reason"] == "missing_concrete_reactant_origin"
        for w in result["steps"][0]["witnesses"]
        for p in w["atom_origins"]["problems"]
    )


@pytest.mark.usefixtures("concrete_origin_metadata")
def test_identical_reactant_instances_are_distinct_equivalents(tmp_path):
    result = audit_fixture(tmp_path, ["C"], "[C:1].[C:2]>>[C:1][C:2]", repeated=True)
    assert result["strict_eligible"]
    step = result["steps"][0]
    assert [s["node_id"] for s in step["reactant_slots"]] == ["b0", "b0"]
    assert [
        (a["reactant_slot"], a["reactant_atom_index"]) for a in step["witnesses"][0]["atom_origins"]["atoms"]
    ] == [(0, 0), (1, 0)]


@pytest.mark.usefixtures("concrete_origin_metadata")
def test_shared_dag_intermediate_reused_in_two_slots_keeps_parent_binding(tmp_path):
    route, catalog, rules = fixture(tmp_path, ["C"], "[C:1].[C:2]>>[C:1][C:2]", repeated=True)
    route["nodes"][-1]["id"] = "middle"
    products, _ = rules.apply("r", ("CC", "CC"))
    route["nodes"].append(
        {"id": "r", "rule_id": "r", "reactants": ["middle", "middle"], "product": products[0]}
    )
    result = audit_route(route, catalog, rules)
    assert result["strict_eligible"]
    first, second = result["steps"]
    assert [s["parent_binding_sha256"] for s in second["reactant_slots"]] == [
        first["step_binding_sha256"]
    ] * 2
    assert {a["reactant_slot"] for a in second["witnesses"][0]["atom_origins"]["atoms"]} == {0, 1}


@pytest.mark.usefixtures("concrete_origin_metadata")
def test_symmetric_origin_ambiguity_and_bounded_truncation_retained(tmp_path):
    route, catalog, rules = fixture(tmp_path, ["CCC"], "[C:1]>>[C:1]")
    result = audit_route(route, catalog, rules)
    step = result["steps"][0]
    assert len(step["witnesses"]) == 3
    assert step["origin_mapping_ambiguous"]
    assert step["origin_mapping_distinct_witnesses"] > 1
    capped = audit_route(route, catalog, rules, product_limit=1)
    assert capped["enumeration_truncated"]
    assert capped["steps"][0]["enumeration"]["examined_entries"] == 1
    assert capped["route_binding_sha256"] != result["route_binding_sha256"]


@pytest.mark.parametrize("smiles", ["[13CH3]C", "[2H]C([2H])C", "[3H]CC"])
@pytest.mark.usefixtures("concrete_origin_metadata")
def test_isotopes_and_explicit_isotope_hydrogens_have_concrete_origins(tmp_path, smiles):
    result = audit_fixture(tmp_path, [smiles], "[C:1]>>[C:1]")
    assert result["strict_eligible"]
    atoms = result["steps"][0]["witnesses"][0]["atom_origins"]["atoms"]
    assert any(a["isotope"] for a in atoms)
    assert all(a["isotope"] == a["source_isotope"] for a in atoms)


@pytest.mark.usefixtures("concrete_origin_metadata")
def test_isotope_transmutation_is_not_a_valid_donor(tmp_path):
    result = audit_fixture(tmp_path, ["C"], "[C:1]>>[13C:1]")
    assert not result["strict_eligible"]
    assert (
        result["steps"][0]["witnesses"][0]["atom_origins"]["problems"][0]["reason"]
        == "element_or_isotope_transmutation"
    )


def test_element_transmutation_is_not_a_valid_donor(tmp_path):
    result = audit_fixture(tmp_path, ["C"], "[C:1]>>[N:1]")
    assert not result["strict_eligible"]
    assert result["atom_origin_status"] == "unsupported"


def test_duplicate_donor_and_lost_metadata_fail_closed():
    product, source = Chem.MolFromSmiles("CC"), Chem.MolFromSmiles("CC")
    for atom in product.GetAtoms():
        atom.SetIntProp("react_idx", 0)
        atom.SetIntProp("react_atom_idx", 0)
    _, result = atom_origins(product, [source], {})
    assert result["status"] == "unsupported"
    assert any(p["reason"] == "duplicate_donor_atom_within_slot" for p in result["problems"])
    product.GetAtomWithIdx(0).ClearProp("react_idx")
    _, result = atom_origins(product, [source], {})
    assert any(p["reason"] == "missing_concrete_reactant_origin" for p in result["problems"])


@pytest.mark.parametrize(
    "smiles,smarts,kind",
    [
        ("F[C@H](Cl)Br", "[C:1]>>[C:1]", "tetrahedral"),
        ("F/C=C/Cl", "[C:1]=[C:2]>>[C:1]=[C:2]", "double_bond"),
        ("OC[C@H](CN)F", "[CH2:1][OH:2]>>[CH3:1]", "tetrahedral"),
        ("OC/C(CN)=C(/Cl)Br", "[CH2:1][OH:2]>>[CH3:1]", "double_bond"),
    ],
)
@pytest.mark.usefixtures("concrete_origin_metadata")
def test_retained_geometry_uses_mapped_neighbors_not_priority_labels(tmp_path, smiles, smarts, kind):
    result = audit_fixture(tmp_path, [smiles], smarts)
    assert result["strict_eligible"]
    assert changes(result, kind)
    assert {c["classification"] for c in changes(result, kind)} == {"retained"}


@pytest.mark.usefixtures("concrete_origin_metadata")
def test_retained_tetra_geometry_can_change_cip_label(tmp_path):
    result = audit_fixture(tmp_path, ["OC[C@H](CN)F"], "[CH2:1][OH:2]>>[CH3:1]")
    before = Chem.FindMolChiralCenters(Chem.MolFromSmiles(result["building_blocks"][0]["smiles"]))
    after = Chem.FindMolChiralCenters(Chem.MolFromSmiles(result["product"]))
    assert before[0][1] != after[0][1]
    assert changes(result)[0]["classification"] == "retained"


@pytest.mark.usefixtures("concrete_origin_metadata")
def test_template_specified_inversion_requires_external_stereo_support(tmp_path):
    result = audit_fixture(
        tmp_path, ["F[C@H](Cl)Br"], "[C@:1]([F:2])([Cl:3])[Br:4]>>[C@@:1]([F:2])([Cl:3])[Br:4]"
    )
    assert result["strict_eligible"]
    assert result["requires_stereo_prediction"]
    assert {c["classification"] for c in changes(result)} == {"inverted"}


@pytest.mark.parametrize(
    "smiles,smarts",
    [
        ("FC(Cl)Br", "[C:1]>>[C:1]"),
        ("FC=CCl", "[C:1]=[C:2]>>[C:1]=[C:2]"),
    ],
)
def test_unspecified_route_state_stereo_is_not_rescued_by_exact_replay(tmp_path, smiles, smarts):
    result = audit_fixture(tmp_path, [smiles], smarts)
    assert not result["strict_eligible"]
    assert result["stereo_status"] == "unsupported"


@pytest.mark.usefixtures("concrete_origin_metadata")
def test_new_specified_center_remains_template_asserted(tmp_path):
    # Alkene reactant is not stereogenic (two identical F); product tetrahedral
    # isotope substitution has explicit stereo but cannot pass atom-origin audit.
    result = audit_fixture(
        tmp_path, ["FC(F)=CCl"], "[C:1]([F:2])([F:3])=[C:4][Cl:5]>>[C@:1]([F:2])([18F:3])[C:4][Cl:5]"
    )
    assert changes(result)[0]["classification"] == "new_or_template_asserted"
    assert result["requires_stereo_prediction"]
    assert not result["strict_eligible"]  # isotope creation is independently unsupported


def test_bindings_include_node_ids_order_and_asset_hashes(assets, route):
    original = audit_route(route, *assets)
    changed = copy.deepcopy(route)
    changed["nodes"][0]["id"] = "renamed"
    changed["nodes"][-1]["reactants"][0] = "renamed"
    assert audit_route(changed, *assets)["route_binding_sha256"] != original["route_binding_sha256"]
    changed = copy.deepcopy(route)
    changed["nodes"][-1]["product"] = canonical("CCC(=O)Nc1ccncc1")
    with pytest.raises(ValueError, match="regenerate"):
        audit_route(changed, *assets)
    assets[1].path.write_text(assets[1].path.read_text() + " ")
    with pytest.raises(ValueError, match="changed"):
        audit_route(route, *assets)


@pytest.mark.parametrize("limit", [0, -1, True, 65])
def test_product_limit_is_strictly_bounded(assets, route, limit):
    with pytest.raises(ValueError):
        audit_route(route, *assets, product_limit=limit)


@pytest.mark.parametrize(
    "smiles,smarts,kind",
    [
        ("O=C(F)Cl", "[C:1](=[O:2])([F:3])[Cl:4]>>[C@H:1]([OH:2])([F:3])[Cl:4]", "tetrahedral"),
        ("FC#CCl", "[F:3][C:1]#[C:2][Cl:4]>>[F:3]/[C:1]=[C:2]/[Cl:4]", "double_bond"),
    ],
)
@pytest.mark.usefixtures("concrete_origin_metadata")
def test_new_stereo_with_complete_atom_origins_is_conditional_computational_support(
    tmp_path, smiles, smarts, kind
):
    result = audit_fixture(tmp_path, [smiles], smarts)
    assert result["strict_eligible"]
    assert result["requires_stereo_prediction"]
    assert {c["classification"] for c in changes(result, kind)} == {"new_or_template_asserted"}


@pytest.mark.usefixtures("concrete_origin_metadata")
def test_double_bond_inversion_is_mapped_geometry_not_label_assertion(tmp_path):
    result = audit_fixture(tmp_path, ["F/C=C/Cl"], r"[F:3]/[C:1]=[C:2]/[Cl:4]>>[F:3]/[C:1]=[C:2]\[Cl:4]")
    assert result["strict_eligible"]
    assert result["requires_stereo_prediction"]
    assert {c["classification"] for c in changes(result, "double_bond")} == {"inverted"}


def test_map_collision_diagnostics_are_not_silently_collapsed():
    from rdkit.Chem import rdChemReactions
    from chembounce_synthesis.route_audit import template_maps

    _, result = template_maps(rdChemReactions.ReactionFromSmarts("[C:1].[C:1]>>[C:1][C:1]"))
    assert result["reactant_map_collisions"] == [1]
    assert result["product_map_collisions"] == [1]


def test_verification_status_dictionary_is_not_an_authoritative_route(assets, route):
    from chembounce_synthesis.routes import verify_route

    with pytest.raises(ValueError, match="schema_version"):
        audit_route(verify_route(route, *assets), *assets)


def test_ordered_slots_with_identical_structures_have_distinct_bound_parent_ids(tmp_path):
    route, catalog, rules = fixture(tmp_path, ["C", "C"], "[C:1].[C:2]>>[C:1][C:2]")
    first = audit_route(route, catalog, rules)
    route["nodes"][-1]["reactants"].reverse()
    second = audit_route(route, catalog, rules)
    assert first["product"] == second["product"]
    assert first["steps"][0]["step_binding_sha256"] != second["steps"][0]["step_binding_sha256"]
    assert first["route_binding_sha256"] != second["route_binding_sha256"]


def test_runtime_without_concrete_origin_metadata_fails_closed(tmp_path):
    result = audit_fixture(tmp_path, ["CC"], "[C:1]>>[C:1]")
    if has_concrete_origin_metadata():
        assert result["strict_eligible"]
        assert result["atom_origin_status"] == "complete"
    else:
        assert result["strict_eligible"] is False
        assert result["atom_origin_status"] == "unsupported"
        assert all(
            p["reason"] == "missing_concrete_reactant_origin"
            for step in result["steps"]
            for witness in step["witnesses"]
            for p in witness["atom_origins"]["problems"]
        )
        assert all(
            witness["atom_origins"]["problems"] for step in result["steps"] for witness in step["witnesses"]
        )
