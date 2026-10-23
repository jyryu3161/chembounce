"""Synthetic catalog/rule controls for explicit route-first generation."""

import copy
import json
import time
from types import SimpleNamespace

import pytest

from chembounce_synthesis.route_audit import audit_route
from chembounce_synthesis.common import canonical, molecule
from chembounce_synthesis.route_first import (
    HitStereoAnchors,
    ReverseHintError,
    RouteFirstLimits,
    RouteFirstSearch,
    check_hit_stereo_anchors,
    generate,
    normalize_reverse_hint,
    reverse_partition_issue,
    reverse_role_issue,
    typed_reverse_hint,
)
from chembounce_synthesis.routes import verify_route
from chembounce_synthesis.workflow import check_constraints, run

HIT = "CC(=O)Nc1ccccc1"
TWO_STEP_HIT = "CC(=O)NCC(=O)Nc1ccccc1"


def limits(**overrides):
    return {
        "max_rule_inspections": 8,
        "max_rule_scans_per_target": 8,
        "max_rules": 1,
        "fallback_rules": 1,
        "max_catalog_rows_per_query": 20,
        "max_catalog_row_visits": 200,
        "max_fingerprint_calls": 200,
        "max_combinations": 100,
        "max_products": 200,
        "max_routes": 100,
        "max_steps": 2,
        **overrides,
    }


def test_actual_one_and_two_step_routes_have_complete_forward_and_atom_stereo_replay(assets):
    catalog, rules = assets
    routes, report = generate(TWO_STEP_HIT, catalog, rules, limits())
    assert report["depth_counts"]["1"] > 0
    assert report["depth_counts"]["2"] > 0
    assert report["catalog_row_visits"] <= 200
    assert report["fingerprint_calls"] <= 200
    assert report["total_fingerprint_calls"] <= 200
    assert report["combinations"] <= 100
    assert (
        0
        < report["inner_role_compatible"]
        <= report["inner_role_checks"]
        <= report["products_generated"]
        <= 200
    )
    assert report["route_audit_attempts"] <= report["limits"]["max_routes"]
    assert all(item["generation"]["reference"] == "self_product_design" for item in routes)
    for item in routes:
        verified = verify_route(item["route"], catalog, rules, 16)
        audit = audit_route(item["route"], catalog, rules, 16)
        assert verified["product"] == audit["product"] == item["verification"]["product"]
        assert audit["strict_eligible"]
        assert item["generation"]["route_audit_binding_sha256"] == audit["route_binding_sha256"]
    chiral = next(item for item in routes if "@" in item["verification"]["product"])
    chiral_audit = audit_route(chiral["route"], catalog, rules, 16)
    assert any(
        change["classification"] == "retained"
        for step in chiral_audit["steps"]
        for witness in step["witnesses"]
        for change in witness["stereo"]["changes"]
    )
    two = next(item for item in routes if len(item["verification"]["steps"]) == 2)
    forged = copy.deepcopy(two["route"])
    inner = next(node for node in forged["nodes"] if node["id"] == "ir")
    inner["product"] = "C"
    with pytest.raises(ValueError, match="regenerate"):
        verify_route(forged, catalog, rules, 16)


def test_product_anchor_and_stereo_are_not_invented_and_search_is_deterministic(assets):
    catalog, rules = assets

    def accept(product):
        check_constraints(product, HIT, product, "hop", core_smiles="Nc1ccccn1")

    first, report = generate(HIT, catalog, rules, limits(max_steps=1), accept=accept)
    second, again = generate(HIT, catalog, rules, limits(max_steps=1), accept=accept)
    assert first == second
    assert report == again
    assert first and all("n" in item["verification"]["product"] for item in first)
    assert any("Reconstructed" in reason for reason in report["product_rejections"])
    with pytest.raises(ValueError, match="Unknown route_first option"):
        generate(HIT, catalog, rules, {"unbounded": True})


def test_global_retrieval_combination_and_deadline_caps(assets):
    catalog, rules = assets
    _, report = generate(
        TWO_STEP_HIT,
        catalog,
        rules,
        limits(max_catalog_rows_per_query=1, max_catalog_row_visits=2, max_fingerprint_calls=1),
    )
    assert report["catalog_row_visits"] <= 2
    assert report["fingerprint_calls"] <= 1
    _, expired = generate(HIT, catalog, rules, limits(), deadline_monotonic=time.monotonic() - 1)
    assert expired["routes_verified"] == 0
    assert "deadline_monotonic" in expired["limits_reached"]
    _, combinations = generate(TWO_STEP_HIT, catalog, rules, limits(max_combinations=2))
    assert combinations["combinations"] <= 2
    _, capped = generate(TWO_STEP_HIT, catalog, rules, limits(max_routes=2))
    assert capped["routes_verified"] <= 2
    assert capped["route_audit_attempts"] <= 2
    assert capped["depth_counts"]["2"] >= 1


def test_reverse_guides_require_disjoint_complete_parent_atom_partition(tmp_path):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "ether",
                        "smarts": "[#6;$([#6]~[#6]);!$([#6]=O):2][#8;H1:3]."
                        "[Cl,Br,I][#6;H2;$([#6]~[#6]):4]>>[CH2:4][O:3][#6:2]",
                        "source": "test",
                    },
                    {"id": "split", "smarts": "[C:1].[C:2]>>[C:1][C:2]", "source": "test"},
                    {"id": "unary", "smarts": "[C:1][O:2]>>[C:1][O:2]", "source": "test"},
                ],
            }
        )
    )
    rules = Rules(path)
    catalog = SimpleNamespace(sha256="test")
    settings = RouteFirstLimits(**limits(max_steps=1, max_reverse_outcomes=1))
    overlap = RouteFirstSearch("C1CCOCC1", catalog, rules, settings)
    assert overlap._reverse(rules.get("ether"), "C1CCOCC1", False) == []
    assert overlap.report["reverse_partition_checks"] == 1
    assert overlap.report["reverse_partition_rejections"] == {"reused_parent_atom": 1}
    assert overlap.report["reverse_fallbacks"] == 0
    assert overlap._reverse(rules.get("ether"), "C1CCOCC1", True)
    assert overlap.report["reverse_partition_checks"] == 1  # explicit fallback skips reverse

    linear = RouteFirstSearch("CC", catalog, rules, settings)
    guides = linear._reverse(rules.get("split"), "CC", False)
    assert guides and len(guides[0]) == 2
    assert linear.report["reverse_partition_checks"] == 1
    assert linear.report["reverse_partition_rejections"] == {}
    assert linear.report["reverse_fallbacks"] == 0

    unary = RouteFirstSearch("CCO", catalog, rules, settings)
    assert unary._reverse(rules.get("unary"), "CCO", False)
    assert unary.report["reverse_partition_checks"] == 1
    assert unary.report["reverse_partition_rejections"] == {}


def test_reverse_partition_metadata_rejects_missing_reused_and_invalid_origins():
    from rdkit import Chem

    def outcome(records):
        mol = Chem.MolFromSmiles("C" * len(records))
        for atom, record in zip(mol.GetAtoms(), records):
            for key, value in record.items():
                atom.SetIntProp(key, value)
        return (mol,)

    def origin(index):
        return {"react_idx": 0, "react_atom_idx": index}

    assert reverse_partition_issue(outcome([origin(0), {}]), 1) is None
    assert reverse_partition_issue(outcome([origin(0), origin(0)]), 2) == "reused_parent_atom"
    assert reverse_partition_issue(outcome([origin(0), {}]), 2) == "missing_parent_atoms"
    assert reverse_partition_issue(outcome([{"react_idx": 0}]), 1) == "incomplete_origin_metadata"
    assert reverse_partition_issue(outcome([{"react_atom_idx": 0}]), 1) == "incomplete_origin_metadata"
    assert reverse_partition_issue(outcome([{"react_idx": 1, "react_atom_idx": 0}]), 1) == "wrong_parent_slot"
    assert reverse_partition_issue(outcome([origin(1)]), 1) == "origin_out_of_range"


def test_reverse_hint_normalization_is_copy_only_and_bounded(tmp_path):
    from rdkit import Chem
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "join", "smarts": "[C:1].[N:2]>>[C:1][N:2]", "source": "test"}],
            }
        )
    )
    rule = Rules(path).get("join")
    parent = molecule("CN1CCCCC1")
    raw = molecule("N1CCCCC1")
    nitrogen = raw.GetAtomWithIdx(0)
    nitrogen.SetIntProp("react_idx", 0)
    nitrogen.SetIntProp("react_atom_idx", 1)
    nitrogen.SetIntProp("old_mapno", 2)
    nitrogen.SetNumExplicitHs(2)
    nitrogen.SetNoImplicit(True)
    hint, edits = normalize_reverse_hint(raw, parent, rule)
    assert hint == canonical("N1CCCCC1")
    assert edits["reset_template_n_hydrogens"] == 1
    assert raw.GetAtomWithIdx(0).GetNumExplicitHs() == 2  # raw reverse witness untouched

    raw_aromatic = Chem.MolFromSmiles("*c1ccc([C@H](F)Cl)cc1")
    normalized, edits = normalize_reverse_hint(raw_aromatic, parent, rule)
    assert "@" in normalized and "*" not in normalized
    assert edits["removed_originless_terminal_dummy"] == 1
    assert raw_aromatic.GetAtomWithIdx(0).GetAtomicNum() == 0

    bridge = Chem.MolFromSmiles("C*C")
    with pytest.raises(ReverseHintError, match="unsafe_dummy"):
        normalize_reverse_hint(bridge, parent, rule)
    origin_dummy = Chem.MolFromSmiles("*C")
    origin_dummy.GetAtomWithIdx(0).SetIntProp("react_atom_idx", 0)
    with pytest.raises(ReverseHintError, match="unsafe_dummy"):
        normalize_reverse_hint(origin_dummy, parent, rule)
    stereo_reference_dummy = Chem.MolFromSmiles("*/C=C/F")
    with pytest.raises(ReverseHintError, match="unsafe_dummy_stereo"):
        normalize_reverse_hint(stereo_reference_dummy, parent, rule)
    stereo_reference_dummy.GetAtomWithIdx(0).SetIntProp("old_mapno", 3)
    with pytest.raises(ReverseHintError, match="unsupported_placeholder_stereo"):
        typed_reverse_hint(stereo_reference_dummy, parent, rule, {(0, 0): 3}, 0, (9, 0, 0, "F"))


def test_reverse_aryl_amine_hints_reach_retrieval_with_chiral_anchor(tmp_path):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {"id": "aryl-h2", "smarts": "[Cl,Br,I][c:1].[N;H2:2]>>[c:1][N:2]", "source": "test"},
                    {"id": "aryl-n", "smarts": "[Cl,F][c:1].[N:2]>>[c:1][N:2]", "source": "test"},
                ],
            }
        )
    )
    rules = Rules(path)
    hit = "c1ccc(N2CCC[C@H](C)C2)cc1"
    search = RouteFirstSearch(
        hit,
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits(max_steps=1)),
    )
    for ident in ("aryl-h2", "aryl-n"):
        guides = search._reverse(rules.get(ident), hit, False)
        assert len(guides) == 1
        assert len(guides[0]) == 2
        assert all("*" not in part for part in guides[0])
        assert "@" in guides[0][1]
        assert "NH2" not in guides[0][1]
    assert search.report["reverse_partition_rejections"] == {}
    assert search.report["reverse_hint_rejections"] == {}
    assert search.report["reverse_hint_normalizations"]["removed_originless_terminal_dummy"] == 1
    assert search.report["reverse_hint_normalizations"]["reset_template_n_hydrogens"] == 1
    assert search.report["reverse_hint_normalizations"]["typed_originless_placeholder_instantiated"] == 1
    assert search.report["reverse_hint_quality_counts"]["typed"] == 1
    assert search.report["reverse_hint_quality_counts"]["approximate"] == 1
    assert (
        search._guide_quality[("aryl-h2", hit, search._reverse(rules.get("aryl-h2"), hit, False)[0])]
        == "approximate"
    )
    assert (
        search._guide_quality[("aryl-n", hit, search._reverse(rules.get("aryl-n"), hit, False)[0])] == "typed"
    )
    assert search.report["reverse_role_trials"] <= 16 * 2
    assert search.report["reverse_role_work"] <= 1024


def test_typed_reverse_role_correspondence_rejects_unrelated_site_and_truncation():
    from rdkit import Chem

    role = Chem.MolFromSmarts("[N:1]")
    unrelated = Chem.MolFromSmiles("NCCN")
    unrelated.GetAtomWithIdx(1).SetIntProp("old_mapno", 1)
    assert reverse_role_issue(unrelated, role, {}, 0) == "role_mismatch_at_mapped_positions"
    assert reverse_role_issue(Chem.MolFromSmiles("NCCN"), role, {}, 0) == "missing_or_ambiguous_role_origin"

    many = Chem.MolFromSmiles("C" * 70)
    many.GetAtomWithIdx(0).SetIntProp("old_mapno", 1)
    assert reverse_role_issue(many, Chem.MolFromSmarts("[C:1]"), {}, 0) == "role_match_cap"


def test_typed_reverse_attempt_caps_share_all_precursor_parts_and_preserve_approximate_hint(
    tmp_path, monkeypatch
):
    from rdkit import Chem

    import chembounce_synthesis.route_first as route_first
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "aryl", "smarts": "[Cl,F][c:1].[N:2]>>[c:1][N:2]", "source": "test"}],
            }
        )
    )
    rules = Rules(path)
    hit = "c1ccc(N2CCCCC2)cc1"
    search = RouteFirstSearch(hit, SimpleNamespace(sha256="test"), rules, RouteFirstLimits(**limits()))

    def late_choice(_raw, _parent, _rule, _tags, slot, representative=None, count_match=None):
        if slot < 2 and representative[3] == "O-":
            return "IC", 0
        raise ReverseHintError("role_mismatch_at_mapped_positions")

    monkeypatch.setattr(route_first, "typed_reverse_hint", late_choice)
    raw = Chem.MolFromSmiles("*C")
    assert search._typed_reverse_outcome((raw, raw, raw), molecule(hit), rules.get("aryl"), {}) is None
    assert search.report["reverse_role_trials"] == search.report["reverse_role_work"] == 16
    assert search.report["reverse_role_rejections"]["outcome_trial_budget"] == 1

    # Saturation rejects typed construction, then retains a partition-valid
    # approximation from the real reverse result without adding outcomes.
    fresh = RouteFirstSearch(hit, SimpleNamespace(sha256="test"), rules, RouteFirstLimits(**limits()))
    fresh.report["reverse_role_work"] = route_first.MAX_REVERSE_ROLE_WORK_PER_SEARCH
    guides = fresh._reverse(rules.get("aryl"), hit, False)
    assert len(guides) == 1
    assert fresh.report["reverse_role_trials"] == 0
    assert fresh.report["reverse_hint_quality_counts"]["approximate"] == 1
    assert fresh.report["reverse_hint_quality_counts"]["implicit_fallback"] == 0
    assert fresh.report["reverse_role_work"] == route_first.MAX_REVERSE_ROLE_WORK_PER_SEARCH

    expired = RouteFirstSearch(
        hit,
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        deadline_monotonic=time.monotonic() - 1,
    )
    assert expired._reverse(rules.get("aryl"), hit, False) == []
    assert expired.report["reverse_role_work"] == 0


def test_inner_beam_uses_only_intermediates_matching_outer_reactant_role(tmp_path):
    from chembounce_synthesis.catalog import Catalog
    from chembounce_synthesis.reactions import Rules

    catalog_path = tmp_path / "catalog.tsv"
    catalog_path.write_text(
        "id\tsmiles\tsource\tsnapshot\n"
        "amide\tCC(=O)N\ttest\ttest\n"
        "acid\tCC(=O)O\ttest\ttest\n"
        "aniline\tNc1ccccc1\ttest\ttest\n"
    )
    rules_path = tmp_path / "rules.json"
    rules_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {"id": "inner", "smarts": "[C:1]>>[C:1]", "source": "test"},
                    {
                        "id": "outer",
                        "smarts": "[C:1](=[O:2])[O;H1:3].[N;H1,H2:4]>>[C:1](=[O:2])[N:4]",
                        "source": "test",
                    },
                ],
            }
        )
    )
    catalog, rules = Catalog(catalog_path), Rules(rules_path)
    search = RouteFirstSearch(
        HIT,
        catalog,
        rules,
        RouteFirstLimits(**limits(beam_width=1, building_block_top_k=2)),
    )
    search._phase = 2
    target, outer_guide = "CC(=O)N", ("CC(=O)N", "Nc1ccccc1")
    inner = search._stage(rules.get("inner"), (target,), "i")
    ranked = sorted(inner, key=lambda state: (search._rank(state["product"], target), state["product"]))
    assert len(ranked) == 2
    assert ranked[0]["product"] == canonical(target)
    assert search._stage(rules.get("outer"), outer_guide, "o", replace=(0, ranked[:1])) == []

    compatible = search._compatible_inner(inner, rules.get("outer"), 0)
    selected = sorted(
        compatible, key=lambda state: (search._rank(state["product"], target), state["product"])
    )[: search.limits.beam_width]
    assert [state["product"] for state in selected] == [canonical("CC(=O)O")]
    products = search._stage(rules.get("outer"), outer_guide, "o", replace=(0, selected))
    assert canonical(HIT) in {state["product"] for state in products}
    assert search.report["inner_role_checks"] == 2
    assert search.report["inner_role_compatible"] == 1
    assert search.report["inner_role_rejected"] == 1


def test_feature_index_targets_guide_size_without_unbounded_per_query_scan(assets):
    catalog, rules = assets
    search = RouteFirstSearch(
        HIT, catalog, rules, RouteFirstLimits(**limits(max_steps=1, max_catalog_rows_per_query=2))
    )
    small = search._pool("CC(=O)O")
    large = search._pool("O=C(O)c1ccccc1")
    assert len(small) == len(large) == 2
    assert small != large
    assert search.report["catalog_index_rows_scanned"] == len(catalog.mols)
    assert search.report["catalog_bucket_rankings"] == 2 * search.report["catalog_index_buckets"]
    assert search.report["catalog_row_visits"] <= search.limits.max_catalog_row_visits


def test_template_filter_and_inner_early_stop_preserve_full_scan_priority(tmp_path):
    from chembounce_synthesis.reactions import Rules

    specifications = [
        ("a-small", "[C:1]>>[C:1]"),
        ("b-medium", "[C:1][C:2]>>[C:1][C:2]"),
        ("c-invalid", "[C:1][C:2][C:3]>>[C:1][C:2][C:3]N"),
        ("d-large", "[C:1][C:2][C:3]>>[C:1][C:2][C:3]"),
        ("e-tie", "[C:1][C:2][C:3]>>[C:1]([C:2])[C:3]"),
        ("f-fallback", "[C:1].[O:2]>>[C:1]=[O:2]"),
        ("z-outside", "[C:1][C:2][C:3][C:4]>>[C:1][C:2][C:3][C:4]"),
    ]
    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {"id": ident, "smarts": smarts, "source": "test"} for ident, smarts in specifications
                ],
            }
        )
    )
    rules = Rules(path)
    catalog = SimpleNamespace(sha256="test")
    search = RouteFirstSearch(
        "CCCC",
        catalog,
        rules,
        RouteFirstLimits(
            **limits(
                max_steps=1,
                max_rule_inspections=12,
                max_rule_scans_per_target=6,
                max_inner_rules=2,
                max_rules=3,
            )
        ),
    )
    inner = search._select_rules("CCCC", 2)
    assert [rule.id for rule, _ in inner] == ["d-large", "e-tie"]
    assert search.report["rule_template_checks"] == 6
    assert [row["rule_id"] for row in search.report["rule_template_exclusions"]] == ["c-invalid"]
    assert search.report["rule_inspections"] == 8  # 6 static + 2 priority matches
    assert [rule.id for rule, _ in search._select_rules("CCCC", 2)] == ["d-large", "e-tie"]
    assert search.report["rule_inspections"] == 8  # cached selection
    assert "z-outside" not in [rule.id for rule, _ in inner]

    # The outer selection still scans the bounded ID prefix, then applies its
    # original product-first and reactant-role fallback policy.
    outer_search = RouteFirstSearch(
        "CCCC",
        catalog,
        rules,
        RouteFirstLimits(
            **limits(
                max_steps=1,
                max_rule_inspections=12,
                max_rule_scans_per_target=6,
                max_rules=5,
                fallback_rules=1,
            )
        ),
    )
    outer = outer_search._select_rules("CCCC", 5, fallback=1)
    assert [rule.id for rule, fallback in outer if not fallback] == [
        "d-large",
        "e-tie",
        "b-medium",
        "a-small",
    ]
    assert [(rule.id, fallback) for rule, fallback in outer if fallback] == [("f-fallback", True)]
    assert outer_search.report["rule_inspections"] == 11
    assert outer_search.report["rule_inspections"] <= outer_search.limits.max_rule_inspections


def test_template_checks_share_global_rule_budget(tmp_path):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {"id": f"r{i}", "smarts": f"[C:{i + 1}]>>[C:{i + 1}]", "source": "test"} for i in range(5)
                ],
            }
        )
    )
    rules = Rules(path)
    search = RouteFirstSearch(
        "C",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits(max_steps=1, max_rule_inspections=5, max_rule_scans_per_target=5)),
    )
    assert search._select_rules("C", 1) == []
    assert search.report["rule_template_checks"] == search.report["rule_inspections"] == 5
    assert "max_rule_inspections" in search.report["limits_reached"]


def test_exact_catalog_target_replaces_pool_row_within_query_cap(assets):
    catalog, rules = assets
    search = RouteFirstSearch(
        HIT, catalog, rules, RouteFirstLimits(**limits(max_steps=1, max_catalog_rows_per_query=1))
    )
    search._pool = lambda _target: ("CC(=O)O",)
    blocks = search._nearest("Nc1ccccc1", "[N;H1,H2]")
    assert [block.smiles for block in blocks] == ["Nc1ccccc1"]
    assert search.report["catalog_row_visits"] == 1


def test_stereo_anchor_posting_is_used_only_for_matching_guide_and_role():
    hit = canonical("C[C@H]1CCCNC1")
    inverse = canonical("C[C@@H]1CCCNC1")
    ordinary = canonical("N1CCCCC1")
    values = (hit, inverse, ordinary)
    catalog = SimpleNamespace(
        sha256="synthetic",
        mols={smiles: molecule(smiles) for smiles in values},
        by_smiles={smiles: [smiles] for smiles in values},
        blocks={smiles: SimpleNamespace(id=smiles, smiles=smiles) for smiles in values},
    )
    rules = SimpleNamespace(sha256="synthetic")
    search = RouteFirstSearch(
        hit,
        catalog,
        rules,
        RouteFirstLimits(**limits(max_steps=1, max_catalog_rows_per_query=1)),
        hit_anchors=HitStereoAnchors(hit),
    )
    search._index()
    search._pool = lambda _target: (ordinary,)
    guide = canonical("C[C@H]1CCCN(C)C1")
    assert [block.smiles for block in search._nearest(guide, "[N;H1,H2]")] == [hit]
    assert search.report["anchor_index_matches"] == [1]
    assert search.report["anchor_role_compatible"] == 1
    assert search.report["guide_anchor_checks"] == 1
    assert [block.smiles for block in search._nearest(ordinary, "[N;H1,H2]")] == [ordinary]
    assert search.report["anchor_role_compatible"] == 1
    assert search.report["catalog_row_visits"] <= search.limits.max_catalog_row_visits


def test_original_hit_stereo_anchor_preserves_local_geometry_and_rejects_inversion():
    hit = "N[C@@H](C)C(=O)Nc1ccccc1"
    retained = "C[C@H](N)C(=O)Nc1ccncc1"
    inverted = "C[C@@H](N)C(=O)Nc1ccncc1"
    check_hit_stereo_anchors(hit, retained)
    with pytest.raises(ValueError, match="inverted"):
        check_hit_stereo_anchors(hit, inverted)
    with pytest.raises(ValueError, match="missing or has ambiguous"):
        check_hit_stereo_anchors(hit, "CC(=O)Nc1ccncc1")
    aromatic = "N[C@H](C)c1ccccc1"
    check_hit_stereo_anchors(aromatic, aromatic)
    alkene = "C/C=C/c1ccccc1"
    check_hit_stereo_anchors(alkene, alkene)
    with pytest.raises(ValueError, match="inverted"):
        check_hit_stereo_anchors(alkene, "C/C=C\\c1ccccc1")
    with pytest.raises(ValueError, match="inverted"):
        check_hit_stereo_anchors("N[C@H](CCCCl)CCCBr", "N[C@@H](CCCCl)CCCBr")
    check_hit_stereo_anchors("N[C@H](CCCCl)CCCBr", "N[C@H](CCCCl)CCCBr")
    long_alkene = "CCC/C(CCCO)=C(/F)Cl"
    check_hit_stereo_anchors(long_alkene, long_alkene)
    with pytest.raises(ValueError, match="inverted"):
        check_hit_stereo_anchors(long_alkene, "CCC/C(CCCO)=C(\\F)Cl")
    complex_hit = "O=C(O)COC[C@H]1CC[C@@H](Cn2nc(-c3ccccc3)c(-c3ccccc3)c2SCCO)CC1"
    check_hit_stereo_anchors(complex_hit, complex_hit)


def test_unmapped_product_atom_and_unspecified_stereo_fail_strict_route_audit(tmp_path):
    from chembounce_synthesis.catalog import Catalog
    from chembounce_synthesis.reactions import Rules

    catalog_path = tmp_path / "catalog.tsv"
    catalog_path.write_text("id\tsmiles\tsource\tsnapshot\na\tC\ttest\t2026\nb\tCC\ttest\t2026\n")
    rules_path = tmp_path / "rules.json"
    rules_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "bad", "smarts": "[C:1].[C:2]>>[C:1][C:2]C", "source": "synthetic"}],
            }
        )
    )
    routes, report = generate("CCC", Catalog(catalog_path), Rules(rules_path), limits(max_steps=1))
    assert routes == []
    assert report["rule_template_exclusions"][0]["rule_id"] == "bad"
    assert report["rule_template_exclusions"][0]["template_audit"]["unbound_product_atoms"]
    assert report["route_audit_attempts"] == 0

    catalog_path.write_text(
        "id\tsmiles\tsource\tsnapshot\n"
        "acid\tCC(=O)O\ttest\t2026\n"
        "alanine\tNC(C)C(=O)O\ttest\t2026\n"
        "phenyl\tNc1ccccc1\ttest\t2026\n"
        "pyridyl\tNc1ccccn1\ttest\t2026\n"
    )
    rules_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "amide",
                        "smarts": "[C:1](=[O:2])[O;H1:3].[N;H1,H2:4]>>[C:1](=[O:2])[N:4]",
                        "source": "synthetic",
                    }
                ],
            }
        )
    )
    routes, report = generate(
        "CC(=O)N[C@H](C)C(=O)Nc1ccccc1",
        Catalog(catalog_path),
        Rules(rules_path),
        limits(),
    )
    assert routes and all("@" not in item["verification"]["product"] for item in routes)
    assert any("atom/stereo audit" in reason for reason in report["product_rejections"])


def test_real_workflow_uses_generated_route_without_proposal_provider_or_native_hopping(
    assets, monkeypatch, tmp_path
):
    from chembounce_synthesis import workflow
    import chembounce

    def forbidden(*args, **kwargs):
        pytest.fail("SynLlama/native hopping must not run in route-first mode")

    monkeypatch.setattr(workflow, "SynLlama", forbidden)
    monkeypatch.setattr(chembounce, "chembounce", forbidden)
    catalog, rules = assets
    config = {
        "schema_version": 1,
        "catalog": str(catalog.path),
        "rules": str(rules.path),
        "route_first": limits(max_steps=1),
        "max_designs": 3,
        "max_scored_products": 3,
        "conformers": {"num_conformers": 2},
    }
    output = tmp_path / "route-first"
    rows = run(HIT, (), config, output)
    assert rows
    assert all(row["design"] == row["smiles"] for row in rows)
    assert json.loads((output / "proposals.json").read_text())["proposals"] == []
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["provider"] is None
    assert manifest["route_first"]["source_sha256"]["route_first.py"]
    generated = json.loads((output / "route_first_candidates.json").read_text())
    assert generated["candidates"] and len(generated["selected_designs"]) <= 3
    summary = json.loads((output / "summary.json").read_text())
    assert summary["route_first_downstream"]["workflow_reverified_route_records"] > 0
    assert summary["route_first_downstream"]["qualified_products"] == 0
    routes = json.loads((output / "routes.json").read_text())["routes"]
    assert all(record["generation"]["policy"] == "catalog-rule-route-first-v2" for record in routes)


def test_preintake_exact_hit_similarity_gate_does_not_spend_design_or_route_slots(assets, tmp_path):
    catalog, rules = assets
    output = tmp_path / "high-threshold"
    rows = run(
        HIT,
        (),
        {
            "schema_version": 1,
            "catalog": str(catalog.path),
            "rules": str(rules.path),
            "route_first": limits(max_steps=1),
            "max_designs": 3,
            "minimum_scores": {"hit_tanimoto": 0.99},
        },
        output,
    )
    summary = json.loads((output / "summary.json").read_text())
    generated = json.loads((output / "route_first_candidates.json").read_text())
    assert rows == []
    assert generated["selected_designs"] == []
    assert summary["route_first"]["candidate_designs"] == 0
    assert summary["route_first"]["route_audit_attempts"] == 0
    assert summary["route_first"]["admission_fingerprint_calls"] > 0
    assert (
        summary["route_first"]["total_fingerprint_calls"]
        <= summary["route_first"]["limits"]["max_fingerprint_calls"]
    )
    assert any(
        "2D original-hit/design similarity" in reason
        for reason in summary["route_first"]["product_rejections"]
    )


def test_workflow_passes_compiled_chiral_hit_anchors_to_generator(assets, monkeypatch, tmp_path):
    from chembounce_synthesis import workflow

    catalog, rules = assets
    original = workflow.generate_route_first
    seen = []

    def wrapped(*args, **kwargs):
        anchors = kwargs["hit_anchors"]
        assert isinstance(anchors, HitStereoAnchors)
        assert len(anchors.queries) == 1
        seen.append(anchors)
        return original(*args, **kwargs)

    monkeypatch.setattr(workflow, "generate_route_first", wrapped)
    run(
        "N[C@@H](C)C(=O)Nc1ccccc1",
        (),
        {
            "schema_version": 1,
            "catalog": str(catalog.path),
            "rules": str(rules.path),
            "route_first": limits(max_steps=1),
        },
        tmp_path / "chiral-seam",
        scoring_mode="2d",
    )
    assert len(seen) == 1


def test_two_step_product_reaches_bounded_design_intake(assets, tmp_path):
    catalog, rules = assets
    config = {
        "schema_version": 1,
        "catalog": str(catalog.path),
        "rules": str(rules.path),
        "route_first": limits(max_steps=2),
        "max_designs": 3,
        "max_scored_products": 3,
        "conformers": {"num_conformers": 2},
    }
    output = tmp_path / "two-step"
    run(TWO_STEP_HIT, (), config, output)
    generated = json.loads((output / "route_first_candidates.json").read_text())
    selected = set(generated["selected_designs"])
    assert len(selected) <= 3
    routes = json.loads((output / "routes.json").read_text())["routes"]
    assert any(
        record["design"] in selected and len(record["verification"]["steps"]) == 2 for record in routes
    )


def test_public_route_first_api_and_main_cli_skip_hopping_but_retain_strict_property_intersection(
    assets, monkeypatch, tmp_path
):
    import chembounce
    from chembounce_synthesis import recommend_routes, workflow

    catalog, rules = assets
    model = tmp_path / "model"
    model.mkdir()
    config = {
        "schema_version": 1,
        "catalog": str(catalog.path),
        "rules": str(rules.path),
        "route_first": limits(max_steps=1),
        "synthesis_evidence": {"predictor": {"kind": "reactiont5", "model_path": str(model)}},
    }
    calls = []

    def fake_run(hit, designs, effective, output_dir, **kwargs):
        assert designs == ()
        assert "provider" not in effective
        calls.append((effective, kwargs))
        return []

    def forbidden(*args, **kwargs):
        pytest.fail("Native hopping/SynLlama must not run")

    monkeypatch.setattr(workflow, "run", fake_run)
    monkeypatch.setattr(chembounce, "chembounce", forbidden)
    assert recommend_routes(HIT, config, tmp_path / "api") == []
    assert calls[-1][0]["recommendation_profile"] == "qualified-selective-engineering-v7"
    assert calls[-1][1]["scoring_mode"] == "3d"
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    assert (
        chembounce.main(
            [
                "-i",
                HIT,
                "-o",
                str(tmp_path / "main"),
                "--synthesis-config",
                str(path),
                "--route-first",
                "--recommend",
                "--mw_max",
                "400",
            ]
        )
        == 3
    )
    assert calls[-1][0]["property_thresholds"]["mw"][1] == 400
    manifest = json.loads((tmp_path / "main" / "run_manifest.json").read_text())
    assert manifest["native_hopping"] == "not_run"
    assert manifest["design_supply"] == "route_first_catalog_rule_v2"
    assert (
        chembounce.main(
            [
                "-i",
                HIT,
                "-o",
                str(tmp_path / "forbidden"),
                "--synthesis-config",
                str(path),
                "--route-first",
                "--recommend",
                "--no-synllama",
            ]
        )
        == 2
    )


def test_strict_route_first_evidence_requires_each_two_step_forward_match(assets, strict_predictor, tmp_path):
    from chembounce_synthesis import recommend_routes

    catalog, rules = assets
    model = tmp_path / "model"
    model.mkdir()
    config = {
        "schema_version": 1,
        "catalog": str(catalog.path),
        "rules": str(rules.path),
        "route_first": limits(max_steps=2),
        "synthesis_evidence": {"predictor": {"kind": "reactiont5", "model_path": str(model)}},
        "max_designs": 3,
        "max_scored_products": 2,
    }

    def forward(prepared):
        products, _ = rules.apply("amide", tuple(prepared["reactants"]), 16)
        return [products[0], "CC", "CO"] if products else ["CC", "CO", "CCC"]

    strict_predictor.products = forward
    good = tmp_path / "good"
    recommend_routes(TWO_STEP_HIT, config, good)
    good_routes = json.loads((good / "routes.json").read_text())["routes"]
    assert any(
        len(row["verification"]["steps"]) == 2 and row["synthesis_evidence"]["eligibility"]["passed"]
        for row in good_routes
    )
    assert json.loads((good / "manifest.json").read_text())["provider"] is None
    assert (
        json.loads((good / "summary.json").read_text())["route_first_downstream"][
            "producer_supported_route_records"
        ]
        > 0
    )

    strict_predictor.products = lambda prepared: (
        ["CC", "CO", "CCC"] if "NCC(=O)O" in prepared["reactants"] else forward(prepared)
    )
    bad = tmp_path / "bad"
    recommend_routes(TWO_STEP_HIT, config, bad)
    bad_routes = json.loads((bad / "routes.json").read_text())["routes"]
    assert all(
        not row["synthesis_evidence"]["eligibility"]["passed"]
        for row in bad_routes
        if len(row["verification"]["steps"]) == 2
    )
