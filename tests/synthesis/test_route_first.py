"""Synthetic catalog/rule controls for explicit route-first generation."""

import copy
import json
import time
from types import SimpleNamespace

import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from chembounce_synthesis import route_first as route_first_module
from chembounce_synthesis.route_audit import audit_route
from chembounce_synthesis.common import canonical, molecule
from chembounce_synthesis.route_first import (
    HitStereoAnchors,
    ReverseHintError,
    RouteFirstLimits,
    RouteFirstSearch,
    check_hit_stereo_anchors,
    continuation_role_origin,
    generate,
    normalize_reverse_hint,
    reverse_partition_issue,
    reverse_role_issue,
    topology_rule_specification,
    typed_reverse_hint,
)
from chembounce_synthesis.route_first_nh import candidate_role_n, prove_parent_n, prove_raw_n
from chembounce_synthesis.routes import verify_route
from chembounce_synthesis.workflow import check_constraints, run

HIT = "CC(=O)Nc1ccccc1"
TWO_STEP_HIT = "CC(=O)NCC(=O)Nc1ccccc1"


def test_mapped_nh_proof_and_second_source_use_existing_ranked_pool(tmp_path, monkeypatch):
    from chembounce_synthesis.catalog import BuildingBlock
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "couple",
                        "smarts": "[Cl:3][c:1].[N:2]>>[c:1][N:2]",
                        "source": "test",
                    }
                ],
            }
        )
    )
    rules = Rules(path)
    rule = rules.get("couple")
    names = ("CNC", "CN", "N", "Clc1ccccc1", "unknown")
    blocks = {name: BuildingBlock(name, name, "test", "2026") for name in names}
    catalog = SimpleNamespace(
        sha256="test",
        mols={name: molecule("[2H]NC" if name == "unknown" else name) for name in names},
        blocks=blocks,
        by_smiles={name: {name} for name in names},
    )
    settings = RouteFirstLimits(**limits())
    off = RouteFirstSearch("CNc1ccccc1", catalog, rules, settings)
    on = RouteFirstSearch("CNc1ccccc1", catalog, rules, settings, prefer_mapped_nh_second=True)
    guide = off._reverse(rule, off.hit, False, record_centers=True)[0]
    assert on._reverse(rule, on.hit, False, record_centers=True) == [guide]
    assert off.report["reverse_attempts"] == on.report["reverse_attempts"] == 1
    assert off.report["fingerprint_calls"] == on.report["fingerprint_calls"] == 1
    assert off.report["mapped_nh_preference"]["proof_contexts"] == 0
    proof = on._mapped_nh_proof(rule, guide, False)
    assert proof is not None and proof[2] == 1
    assert on._mapped_nh_proof(rule, guide, False) == proof
    assert on.report["mapped_nh_preference"]["proof_contexts"] == 1

    key = (guide[1], rule.roles[1], False)
    legacy = (blocks["CNC"], blocks["N"])
    on._nearest_cache[key] = legacy
    on._mapped_nh_ranked[key] = (
        ((-1.0, 0, "CNC"), (-0.9, 0, "CN"), (-0.8, 0, "N")),
        True,
    )
    selected = on._nearest(guide[1], rule.roles[1], mapped_nh_context=(rule, proof))
    assert selected == (blocks["CNC"], blocks["CN"])
    assert on._nearest_cache[key] == legacy
    assert on.report["catalog_row_visits"] == 0
    assert on.report["fingerprint_calls"] == 1
    assert on.report["mapped_nh_preference"]["replacements"] == 1
    assert on._nearest(guide[1], rule.roles[1]) == legacy

    on.deadline = time.monotonic() - 1
    assert on._nearest(guide[1], rule.roles[1], mapped_nh_context=(rule, proof)) == legacy
    on.deadline = None
    on._mapped_nh_selected.clear()
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_CHECKS_PER_QUERY", 0)
    assert on._nearest(guide[1], rule.roles[1], mapped_nh_context=(rule, proof)) == legacy
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_CHECKS_PER_QUERY", 8192)
    on._mapped_nh_ranked[key] = (((-1.0, 0, "CNC"), (-0.9, 0, "unknown"), (-0.8, 0, "CN")), True)
    assert on._nearest(guide[1], rule.roles[1], mapped_nh_context=(rule, proof)) == legacy
    on._mapped_nh_ranked[key] = (((-1.0, 0, "CNC"), (-0.8, 0, "CN")), False)
    assert on._nearest(guide[1], rule.roles[1], mapped_nh_context=(rule, proof)) == legacy
    on._mapped_nh_ranked.clear()
    assert on._nearest(guide[1], rule.roles[1], mapped_nh_context=(rule, proof)) == legacy


def test_mapped_nh_unknown_proof_and_strict_role_checks_fail_closed(tmp_path, monkeypatch):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "couple", "smarts": "[Cl:3][c:1].[N:2]>>[c:1][N:2]", "source": "test"}],
            }
        )
    )
    rules = Rules(path)
    rule = rules.get("couple")
    search = RouteFirstSearch(
        "CNc1ccccc1",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        prefer_mapped_nh_second=True,
    )
    guide = search._reverse(rule, search.hit, False)[0]
    registry = search._mapped_nh_raw[(rule.id, search.hit, False)]
    registry["raw_cap_saturated"] = True
    assert search._mapped_nh_proof(rule, guide, False) is None
    registry["raw_cap_saturated"] = False
    registry["incomplete"] = "snapshot_budget"
    assert search._mapped_nh_proof(rule, guide, False) is None
    registry["incomplete"] = None
    monkeypatch.setattr(route_first_module, "MAX_MAPPED_NH_PROOF_UNITS", 0)
    assert search._mapped_nh_proof(rule, guide, False) is None
    monkeypatch.setattr(route_first_module, "MAX_MAPPED_NH_PROOF_UNITS", 32768)
    monkeypatch.setattr(route_first_module, "MAX_MAPPED_NH_PROOF_MATCH_REQUESTS", 0)
    assert search._mapped_nh_proof(rule, guide, False) is None
    monkeypatch.setattr(route_first_module, "MAX_MAPPED_NH_PROOF_MATCH_REQUESTS", 512)
    assert search._mapped_nh_proof(rule, guide, False) is not None

    conflicting = RouteFirstSearch(
        "CNc1ccccc1",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        prefer_mapped_nh_second=True,
    )
    assert conflicting._reverse(rule, conflicting.hit, False) == [guide]
    alternative = conflicting._mapped_nh_raw[(rule.id, conflicting.hit, False)]["by_guide"][guide][0]
    raw_parts, clean_parts = alternative
    modified = tuple(Chem.Mol(part) for part in clean_parts)
    site_atom = next(atom for atom in modified[1].GetAtoms() if atom.GetAtomicNum() == 7)
    site_atom.SetIntProp("react_atom_idx", 999)
    conflicting._mapped_nh_raw[(rule.id, conflicting.hit, False)]["by_guide"][guide].append(
        (raw_parts, modified)
    )
    assert conflicting._mapped_nh_proof(rule, guide, False) is None

    parent = molecule("CNc1ccccc1")
    site, reason = prove_parent_n(rule.reaction, parent)
    assert reason is None
    raw = registry["by_guide"][guide][0][0]
    assert prove_raw_n((raw,), parent, site) == (True, None)
    assert candidate_role_n(molecule("CN"), rule.reaction.GetReactants()[1], 2) == (True, None)
    assert candidate_role_n(molecule("CNC"), rule.reaction.GetReactants()[1], 2) == (False, None)
    assert candidate_role_n(molecule("[2H]NC"), rule.reaction.GetReactants()[1], 2)[0] is None
    assert prove_parent_n(rule.reaction, molecule("CN(C)c1ccccc1"))[0] is None
    assert prove_parent_n(rule.reaction, molecule("CNc1ccc(NC)cc1"))[0] is None
    assert RouteFirstLimits(**limits(prefer_mapped_nh_second=True)).prefer_mapped_nh_second
    assert RouteFirstLimits(**limits(prefer_continuation_second=True)).prefer_continuation_second
    with pytest.raises(ValueError, match="must be a boolean"):
        RouteFirstLimits(**limits(prefer_mapped_nh_second=1))
    with pytest.raises(ValueError, match="must be a boolean"):
        RouteFirstLimits(**limits(prefer_continuation_second="yes"))


def test_mapped_nh_approximate_guide_can_prove_raw_n_without_role_promotion(tmp_path):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "couple",
                        "smarts": "[Cl:3][c:1].[N;!$(NS(=O)=O):2]>>[c:1][N:2]",
                        "source": "test",
                    }
                ],
            }
        )
    )
    rules = Rules(path)
    rule = rules.get("couple")
    search = RouteFirstSearch(
        "CS(=O)(=O)Nc1ccccc1",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        prefer_mapped_nh_second=True,
        prefer_continuation_second=True,
    )
    guide = search._reverse(rule, search.hit, False, record_centers=True)[0]
    assert search._guide_quality[(rule.id, search.hit, guide)] == "approximate"
    assert not molecule(guide[1]).HasSubstructMatch(rule.reaction.GetReactants()[1])
    proof = search._mapped_nh_proof(rule, guide, False)
    assert proof is not None and proof[2] == 1
    assert search.report["reverse_attempts"] == 1
    assert search.report["mapped_nh_preference"]["raw_snapshot_attempts"] == 1
    assert search._continuation_snapshots == 0
    search._phase = 2
    assert search._reverse(rule, "CNc1ccccc1", False)
    assert search._continuation_snapshots == 1
    assert search._mapped_nh_snapshots == 1


def test_public_route_first_mapped_nh_opt_in_changes_only_second_source(tmp_path):
    from chembounce_synthesis.catalog import Catalog
    from chembounce_synthesis.reactions import Rules

    catalog_path = tmp_path / "catalog.tsv"
    names = (
        "Clc1ccccc1",
        "CNC1CCCCC1",
        "CCNC1CCCCC1",
        "NC1CCCCC1",
        "CN",
        "CCN",
        "CNC",
        "CC(C)N",
    )
    catalog_path.write_text(
        "id\tsmiles\tsource\tsnapshot\n"
        + "".join(f"{index}\t{name}\ttest\t2026\n" for index, name in enumerate(names))
    )
    rules_path = tmp_path / "rules.json"
    rules_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "couple",
                        "smarts": "[Cl:3][c:1].[N;!$(NS(=O)=O):2]>>[c:1][N:2]",
                        "source": "test",
                    }
                ],
            }
        )
    )
    catalog, rules = Catalog(catalog_path), Rules(rules_path)
    hit = "O=S(=O)(C1CCCCC1)Nc1ccccc1"
    settings = limits(max_steps=1, max_outer_guides=1)
    default, default_report = generate(hit, catalog, rules, settings)
    explicit_off, off_report = generate(hit, catalog, rules, {**settings, "prefer_mapped_nh_second": False})
    on, on_report = generate(hit, catalog, rules, {**settings, "prefer_mapped_nh_second": True})

    def signature(rows):
        return [
            (
                item["verification"]["product"],
                tuple(node.get("building_block_id") for node in item["route"]["nodes"]),
            )
            for item in rows
        ]

    assert signature(default) == signature(explicit_off)
    assert default_report["mapped_nh_preference"]["proof_contexts"] == 0
    assert off_report["mapped_nh_preference"]["proof_contexts"] == 0
    assert on_report["mapped_nh_preference"]["replacements"] == 1
    assert on_report["reverse_attempts"] == default_report["reverse_attempts"]
    assert on_report["fingerprint_calls"] == default_report["fingerprint_calls"]
    assert signature(default)[0] == signature(on)[0]
    assert signature(on)[1][0] == "CNc1ccccc1"
    assert signature(on)[1][1][:2] == ("0", "4")
    assert on_report["routes_verified"] == 2

    both = RouteFirstSearch(
        hit,
        catalog,
        rules,
        RouteFirstLimits(**settings),
        prefer_continuation_second=True,
        prefer_mapped_nh_second=True,
    )
    both.run()
    assert both.report["mapped_nh_preference"]["ranked_rows_retained"] > 0
    assert both._continuation_ranked == {}
    assert both._continuation_selected == {}

    reversed_rule = rdChemReactions.ReactionFromSmarts("[N:2].[Cl:3][c:1]>>[c:1][N:2]")
    reversed_rule.Initialize()
    reversed_site, reason = prove_parent_n(reversed_rule, molecule("CNc1ccccc1"))
    assert reason is None and reversed_site["source_slot"] == 0


def test_continuation_origin_distinguishes_intact_split_and_edited_raw_sites():
    parent = molecule("CCO")
    reverse = rdChemReactions.ReactionFromSmarts("[C:1][C:2][O:3]>>[C:1][C:2].[O:3]")
    reverse.Initialize()
    outcome = reverse.RunReactants((parent,), maxProducts=4)[0]
    assert reverse_partition_issue(outcome, parent.GetNumAtoms()) is None

    intact = continuation_role_origin(parent, outcome, Chem.MolFromSmarts("[CH3]"))
    assert intact["classification"] == "single_source_local_intact"
    assert intact["source_slot"] == 0
    assert intact["matched_parent_atoms"] == [0]
    assert not intact["full_smarts_preservation_proven"]
    assert not intact["concrete_leaf_requirement_proven"]

    edited = continuation_role_origin(parent, outcome, Chem.MolFromSmarts("[CH2]"))
    assert edited["classification"] == "single_source_edited"
    assert edited["source_slot"] == 0
    assert edited["bond_changes"] == [1]

    split = continuation_role_origin(parent, outcome, Chem.MolFromSmarts("[CH2][OH]"))
    assert split["classification"] == "split_across_sources"
    assert split["owner_slots"] == [0, 1]
    assert split["bond_changes"] == [1, 2]


def test_continuation_origin_ambiguity_truncation_and_invalid_origin_stay_unknown():
    parent = molecule("CCO")
    reverse = rdChemReactions.ReactionFromSmarts("[C:1][C:2][O:3]>>[C:1][C:2].[O:3]")
    reverse.Initialize()
    outcome = reverse.RunReactants((parent,), maxProducts=4)[0]
    assert continuation_role_origin(parent, outcome, Chem.MolFromSmarts("[C]"))["reason"] == (
        "ambiguous_role_embeddings"
    )
    assert continuation_role_origin(parent, outcome, Chem.MolFromSmarts("[F]"))["reason"] == (
        "role_not_in_parent"
    )
    assert (
        continuation_role_origin(
            parent, outcome, Chem.MolFromSmarts("[CH3]"), deadline_monotonic=time.monotonic() - 1
        )["reason"]
        == "deadline"
    )

    damaged = tuple(Chem.Mol(part) for part in outcome)
    damaged[0].GetAtomWithIdx(0).ClearProp("react_atom_idx")
    assert (
        continuation_role_origin(parent, damaged, Chem.MolFromSmarts("[CH3]"))["reason"]
        == "partition_incomplete_origin_metadata"
    )

    chain = molecule("C" * 65)
    identity = rdChemReactions.ReactionFromSmarts("[C:1]>>[C:1]")
    identity.Initialize()
    # The bounded query sentinel is independent of reverse outcome length.
    copied = Chem.Mol(chain)
    for atom in copied.GetAtoms():
        atom.SetIntProp("react_idx", 0)
        atom.SetIntProp("react_atom_idx", atom.GetIdx())
    report = continuation_role_origin(chain, (copied,), Chem.MolFromSmarts("[C]"))
    assert report["reason"] == "role_match_truncated"
    assert report["match_count"] == 65


def test_continuation_origin_changes_recursive_stereo_and_bounds_are_conservative(monkeypatch):
    parent = molecule("CO")

    def sourced_copy(mol):
        raw = Chem.Mol(mol)
        for atom in raw.GetAtoms():
            atom.SetIntProp("react_idx", 0)
            atom.SetIntProp("react_atom_idx", atom.GetIdx())
        return raw

    role = Chem.MolFromSmarts("[C;$([C]-[O])]")
    raw = sourced_copy(parent)
    raw.GetAtomWithIdx(0).SetFormalCharge(1)
    raw.GetAtomWithIdx(0).SetNumExplicitHs(1)
    raw.GetAtomWithIdx(0).SetIsAromatic(True)
    changed = continuation_role_origin(parent, (raw,), role)
    assert changed["classification"] == "single_source_edited"
    assert changed["atom_changes"] == [0]
    assert not changed["full_smarts_preservation_proven"]

    split = (Chem.MolFromSmiles("C"), Chem.MolFromSmiles("O"))
    split[0].GetAtomWithIdx(0).SetIntProp("react_idx", 0)
    split[0].GetAtomWithIdx(0).SetIntProp("react_atom_idx", 0)
    split[1].GetAtomWithIdx(0).SetIntProp("react_idx", 0)
    split[1].GetAtomWithIdx(0).SetIntProp("react_atom_idx", 1)
    recursive = continuation_role_origin(parent, split, role)
    assert recursive["classification"] == "single_source_edited"
    assert recursive["bond_changes"] == [0]
    assert not recursive["full_smarts_preservation_proven"]

    symmetric = molecule("CC")
    ambiguous = continuation_role_origin(symmetric, (sourced_copy(symmetric),), Chem.MolFromSmarts("[C][C]"))
    assert ambiguous["reason"] == "ambiguous_role_embeddings"
    assert ambiguous["match_count"] == 2

    chiral = molecule("C[C@H](F)Cl")
    chiral_role = Chem.MolFromSmarts("C[C@H](F)Cl")
    stereo = continuation_role_origin(chiral, (sourced_copy(chiral),), chiral_role)
    assert stereo["reason"] == "specified_stereo_mapping_unproven"

    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_ATOMS", 1)
    assert continuation_role_origin(parent, (sourced_copy(parent),), role)["reason"] == ("parent_size_limit")
    assert (
        continuation_role_origin(molecule("C"), (sourced_copy(molecule("C")),), Chem.MolFromSmarts("[C][O]"))[
            "reason"
        ]
        == "parent_size_limit"
    )
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_ATOMS", 256)
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_PRECURSORS", 1)
    assert continuation_role_origin(parent, split, role)["reason"] == "precursor_size_limit"
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_PRECURSORS", 4)
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_ATOMS", 1)
    assert (
        continuation_role_origin(molecule("C"), (sourced_copy(parent),), Chem.MolFromSmarts("[C]"))["reason"]
        == "precursor_size_limit"
    )
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_ATOMS", 256)
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_BONDS", 0)
    assert continuation_role_origin(parent, (sourced_copy(parent),), role)["reason"] == ("parent_size_limit")
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_BONDS", 512)
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ORIGIN_VISITS", 1)
    assert continuation_role_origin(parent, (sourced_copy(parent),), role)["reason"] == ("origin_visit_limit")


def test_opt_in_continuation_bridge_preserves_reverse_guides_and_cache(tmp_path, monkeypatch):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "continuation_rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "inner",
                        "smarts": "[N:1][C:2][C:3].[O:4]>>[N:1][C:2][C:3][O:4]",
                        "source": "test",
                    }
                ],
            }
        )
    )
    rules = Rules(path)
    catalog = SimpleNamespace(sha256="test")
    settings = RouteFirstLimits(**limits())
    legacy = RouteFirstSearch("NCCO", catalog, rules, settings)
    observed = RouteFirstSearch("NCCO", catalog, rules, settings, diagnose_continuation=True)
    inner = rules.get("inner")
    old_guides = legacy._reverse(inner, "NCCO", False)
    new_guides = observed._reverse(inner, "NCCO", False)
    assert old_guides == new_guides == [("CCN", "O")]
    assert legacy.report == observed.report
    assert legacy._reverse_cache == observed._reverse_cache
    assert legacy._guide_quality == observed._guide_quality
    assert legacy._guide_centers == observed._guide_centers
    assert observed._continuation_snapshots == 1
    assert legacy._reverse(inner, "NCCO", False, record_centers=True) == observed._reverse(
        inner, "NCCO", False, record_centers=True
    )
    assert legacy._guide_centers == observed._guide_centers
    assert legacy.report == observed.report

    outer = SimpleNamespace(
        id="outer", reaction=rdChemReactions.ReactionFromSmarts("[N:1].[C:2]>>[N:1][C:2]")
    )
    guide = new_guides[0]
    positive = observed.continuation_guide_origin(inner, "NCCO", False, guide, outer, 0)
    assert positive["classification"] == "single_source_local_intact_normalized_role"
    assert positive["source_slot"] == 0
    assert positive["full_role_on_normalized_guide_same_site"]
    assert not positive["concrete_leaf_requirement_proven"]
    count = observed._continuation_role_contexts
    assert observed.continuation_guide_origin(inner, "NCCO", False, guide, outer, 0) == positive
    assert observed._continuation_role_contexts == count

    registry = observed._continuation_raw[(inner.id, "NCCO", False)]
    registry["raw_cap_saturated"] = True
    observed._continuation_context_cache.clear()
    assert observed.continuation_guide_origin(inner, "NCCO", False, guide, outer, 0)["reason"] == (
        "raw_outcome_cap_saturated"
    )
    attempts = observed.report["reverse_attempts"]
    observed._continuation_raw.clear()
    assert observed._reverse(inner, "NCCO", False) == new_guides
    assert observed.report["reverse_attempts"] == attempts
    assert observed.continuation_guide_origin(inner, "NCCO", False, guide, outer, 0)["reason"] == (
        "raw_provenance_missing_or_cached_before_capture"
    )

    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_CONTEXTS", 0)
    fresh = RouteFirstSearch("NCCO", catalog, rules, settings, diagnose_continuation=True)
    fresh._reverse(inner, "NCCO", False)
    assert fresh.continuation_guide_origin(inner, "NCCO", False, guide, outer, 0)["reason"] == (
        "role_context_budget"
    )


def test_continuation_bridge_ambiguous_normalized_site_and_snapshot_budget(tmp_path, monkeypatch):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "inner", "smarts": "[C:1][C:2].[O:3]>>[C:1][C:2][O:3]", "source": "test"}],
            }
        )
    )
    rules = Rules(path)
    inner = rules.get("inner")
    outer = SimpleNamespace(
        id="outer", reaction=rdChemReactions.ReactionFromSmarts("[CH3:1].[N:2]>>[C:1][N:2]")
    )

    def search():
        return RouteFirstSearch(
            "CCO",
            SimpleNamespace(sha256="test"),
            rules,
            RouteFirstLimits(**limits()),
            diagnose_continuation=True,
        )

    observed = search()
    guide = observed._reverse(inner, "CCO", False)[0]
    assert observed.continuation_guide_origin(inner, "CCO", False, guide, outer, 0)["reason"] == (
        "normalized_role_embedding_ambiguous_or_absent"
    )
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_RAW_SNAPSHOTS", 0)
    saturated = search()
    assert saturated._reverse(inner, "CCO", False) == [guide]
    assert saturated.continuation_guide_origin(inner, "CCO", False, guide, outer, 0)["reason"] == (
        "snapshot_budget"
    )

    raw = Chem.MolFromSmiles("CO")
    clean = Chem.Mol(raw)
    for mol in (raw, clean):
        for atom in mol.GetAtoms():
            atom.SetIntProp("react_idx", 0)
            atom.SetIntProp("react_atom_idx", atom.GetIdx())
    clean.GetAtomWithIdx(0).SetNumRadicalElectrons(1)
    assert not route_first_module._normalized_origin_site_unchanged(raw, clean, [0])


def test_continuation_failed_snapshot_copy_is_charged_before_allocation(tmp_path, monkeypatch):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "r", "smarts": "[C:1].[O:2]>>[C:1][O:2]", "source": "test"}],
            }
        )
    )
    search = RouteFirstSearch(
        "CO",
        SimpleNamespace(sha256="test"),
        Rules(path),
        RouteFirstLimits(**limits()),
        diagnose_continuation=True,
    )
    registry = {"by_guide": {}, "incomplete": None, "raw_cap_saturated": False}
    parts = (molecule("C"), molecule("O"))
    with monkeypatch.context() as patched:
        patched.setattr(route_first_module.Chem, "Mol", lambda _: (_ for _ in ()).throw(RuntimeError("copy")))
        search._capture_continuation_raw(registry, ("C", "O"), parts, parts)
    assert registry["incomplete"] == "snapshot_copy_failed"
    assert search._continuation_snapshots == 1
    assert search._continuation_snapshot_units > 0
    assert not registry["by_guide"]


def test_continuation_second_uses_same_ranked_pool_and_can_enable_outer_role(tmp_path, monkeypatch):
    from chembounce_synthesis.catalog import BuildingBlock
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "inner", "smarts": "[C:1][C:2].[O:3]>>[C:1][C:2][O:3]", "source": "test"}],
            }
        )
    )
    rules = Rules(path)
    inner = rules.get("inner")
    outer = SimpleNamespace(
        id="outer", reaction=rdChemReactions.ReactionFromSmarts("[Cl:1][C:2].[N:3]>>[C:2][N:3]")
    )
    choices = ("CCC", "CCCC", "ClCCC", "O")
    blocks = {item: BuildingBlock(item, item, "fixture", "2026-10-24") for item in choices}
    catalog = SimpleNamespace(
        sha256="test",
        mols={item: molecule(item) for item in choices},
        blocks=blocks,
        by_smiles={item: {item} for item in choices},
    )
    search = RouteFirstSearch(
        "ClCCO", catalog, rules, RouteFirstLimits(**limits()), prefer_continuation_second=True
    )
    search._phase = 2
    guide = search._reverse(inner, "ClCCO", False)[0]
    origin = search.continuation_guide_origin(inner, "ClCCO", False, guide, outer, 0)
    assert origin["classification"] == "single_source_local_intact_normalized_role"
    assert origin["source_slot"] == 0

    phase_one = RouteFirstSearch(
        "ClCCO", catalog, rules, RouteFirstLimits(**limits()), prefer_continuation_second=True
    )
    assert phase_one._reverse(inner, "ClCCO", False) == [guide]
    assert phase_one._continuation_snapshots == 0
    phase_one._phase = 2
    assert phase_one._reverse(inner, "ClCCO", False) == [guide]
    assert phase_one.continuation_guide_origin(inner, "ClCCO", False, guide, outer, 0)["reason"] == (
        "raw_provenance_missing_or_cached_before_capture"
    )

    base_key = (guide[0], inner.roles[0], False)
    legacy = (blocks["CCC"], blocks["CCCC"])
    search._nearest_cache[base_key] = legacy
    search._continuation_ranked[base_key] = (
        ((-1.0, 0, "CCC"), (-0.8, 0, "CCCC"), (-0.7, 0, "ClCCC")),
        True,
    )
    context = (outer, 0, 0, inner.id, "ClCCO", guide)
    chosen = search._nearest(guide[0], inner.roles[0], continuation_context=context)
    assert chosen[0] == legacy[0]
    assert chosen[1] == blocks["ClCCC"]
    assert search._nearest_cache[base_key] == legacy
    assert search.report["continuation_preference"]["role_checks"] == 2
    assert search.report["catalog_row_visits"] == 0

    assert search.report["fingerprint_calls"] == 1

    future_role = outer.reaction.GetReactants()[0]
    old_products, _ = rules.apply(inner.id, (legacy[1].smiles, "O"), 16)
    new_products, _ = rules.apply(inner.id, (chosen[1].smiles, "O"), 16)
    assert old_products and not any(molecule(x).HasSubstructMatch(future_role) for x in old_products)
    assert any(molecule(x).HasSubstructMatch(future_role) for x in new_products)

    other = SimpleNamespace(
        id="other", reaction=rdChemReactions.ReactionFromSmarts("[C:1].[N:2]>>[C:1][N:2]")
    )
    other_context = (other, 0, 0, inner.id, "ClCCO", guide)
    assert search._nearest(guide[0], inner.roles[0], continuation_context=other_context) == legacy
    assert search._nearest(guide[0], inner.roles[0], continuation_context=context) == chosen
    assert search.report["catalog_row_visits"] == 0
    search.deadline = time.monotonic() - 1
    assert search._nearest(guide[0], inner.roles[0], continuation_context=context) == legacy
    search.deadline = None

    search._continuation_selected.clear()
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_SELECTED_CONTEXTS", 0)
    assert search._nearest(guide[0], inner.roles[0], continuation_context=context) == legacy
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_SELECTED_CONTEXTS", 128)
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_CHECKS_PER_QUERY", 1)
    assert search._nearest(guide[0], inner.roles[0], continuation_context=context) == legacy
    assert not search._continuation_selected
    monkeypatch.setattr(route_first_module, "MAX_CONTINUATION_ROLE_CHECKS_PER_QUERY", 8192)
    monkeypatch.setattr(
        search,
        "_catalog_block_stereo_eligible",
        lambda smiles: None if smiles == "ClCCC" else True,
    )
    assert search._nearest(guide[0], inner.roles[0], continuation_context=context) == legacy
    assert not search._continuation_selected

    search.deadline = time.monotonic() - 1
    assert search._nearest(guide[0], inner.roles[0], continuation_context=context) == legacy
    search.deadline = None
    search._continuation_ranked.clear()
    assert search._nearest(guide[0], inner.roles[0], continuation_context=context) == legacy
    assert search.report["catalog_row_visits"] == 0


def test_optional_continuation_preference_preserves_strict_one_and_two_step_routes(assets):
    catalog, rules = assets
    search = RouteFirstSearch(
        TWO_STEP_HIT,
        catalog,
        rules,
        RouteFirstLimits(**limits()),
        prefer_continuation_second=True,
    )
    routes, report = search.run()
    assert routes
    assert report["depth_counts"]["1"] > 0
    assert report["depth_counts"]["2"] > 0
    assert report["catalog_row_visits"] <= search.limits.max_catalog_row_visits
    assert report["continuation_preference"]["ranked_rows_retained"] <= 65536
    assert report["continuation_preference"]["role_checks"] <= 65536
    for item in routes:
        audit = audit_route(item["route"], catalog, rules, search.limits.reaction_product_limit)
        assert audit["strict_eligible"]


def test_continuation_rank_retention_is_phase_two_only_and_never_repools_for_context(tmp_path, monkeypatch):
    from chembounce_synthesis.catalog import BuildingBlock
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "r", "smarts": "[C:1][C:2].[O:3]>>[C:1][C:2][O:3]", "source": "test"}],
            }
        )
    )
    rules = Rules(path)
    rows = ("CCC", "CCCC", "ClCCC")
    catalog = SimpleNamespace(
        sha256="test",
        mols={s: molecule(s) for s in rows},
        blocks={s: BuildingBlock(s, s, "fixture", "2026-10-24") for s in rows},
        by_smiles={s: {s} for s in rows},
    )
    outer = SimpleNamespace(
        id="outer", reaction=rdChemReactions.ReactionFromSmarts("[Cl:1][C:2].[N:3]>>[C:2][N:3]")
    )
    context = (outer, 0, 0, "r", "ClCCO", ("CCCl", "O"))
    role = rules.get("r").roles[0]

    def make_search(phase):
        search = RouteFirstSearch(
            "ClCCO", catalog, rules, RouteFirstLimits(**limits()), prefer_continuation_second=True
        )
        search._phase = phase
        calls = []

        def pool(*_args, **_kwargs):
            calls.append(1)
            return rows

        monkeypatch.setattr(search, "_pool", pool)
        return search, calls

    first, calls = make_search(1)
    old = first._nearest("CCCl", role, continuation_context=context)
    assert calls == [1]
    assert first.report["continuation_preference"]["ranked_rows_retained"] == 0
    first._phase = 2
    assert first._nearest("CCCl", role, continuation_context=context) == old
    assert calls == [1]

    with monkeypatch.context() as patched:
        patched.setattr(route_first_module, "MAX_CONTINUATION_RANKED_ROWS", 1)
        small, calls = make_search(2)
        small._nearest("CCCl", role, continuation_context=context)
        assert calls == [1]
        assert small.report["continuation_preference"]["ranked_rows_retained"] == 0
        assert small.report["continuation_preference"]["role_checks"] == 0

    second, calls = make_search(2)
    second._nearest("CCCl", role, continuation_context=context)
    assert calls == [1]
    assert 0 < second.report["continuation_preference"]["ranked_rows_retained"] <= len(rows)
    before = second.report["catalog_row_visits"]
    second._nearest("CCCl", role, continuation_context=(outer, 1, 0, "r", "ClCCO", ("CCCl", "O")))
    assert calls == [1]
    assert second.report["catalog_row_visits"] == before


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


def test_empty_stage_retrieval_opt_in_preserves_legacy_calls_and_skips_unneeded_roles():
    rule = SimpleNamespace(id="join", roles=("role0", "role1"))
    applied = []
    rules = SimpleNamespace(
        sha256="test",
        apply=lambda _rule, reactants, _limit: (applied.append(reactants) or ("CN",), False),
    )
    catalog = SimpleNamespace(sha256="test")

    def search(enabled):
        return RouteFirstSearch(
            "CN",
            catalog,
            rules,
            RouteFirstLimits(**limits(skip_empty_stage_retrieval=enabled)),
        )

    off, on = search(False), search(True)
    off_calls, on_calls = [], []
    off._nearest = lambda target, _role, **_kw: (
        off_calls.append(target) or (() if target == "empty" else (SimpleNamespace(id="b", smiles="N"),))
    )
    on._nearest = lambda target, _role, **_kw: (
        on_calls.append(target) or (() if target == "empty" else (SimpleNamespace(id="b", smiles="N"),))
    )
    assert off._stage(rule, ("empty", "later"), "x") == []
    assert on._stage(rule, ("empty", "later"), "x") == []
    assert off_calls == ["empty", "later"] and on_calls == ["empty"]
    assert off.report["empty_stage_short_circuits"] == 0
    assert on.report["empty_stage_short_circuits"] == 1
    assert on.report["empty_stage_unvisited_role_queries"] == 1
    assert off.report["combinations"] == on.report["combinations"] == 0
    assert off.report["products_generated"] == on.report["products_generated"] == 0
    assert not applied

    # A known-empty replacement at slot 1 prevents even the earlier slot-0 lookup.
    off_calls.clear()
    on_calls.clear()
    assert off._stage(rule, ("earlier", "replaced"), "x", replace=(1, [])) == []
    assert on._stage(rule, ("earlier", "replaced"), "x", replace=(1, [])) == []
    assert off_calls == ["earlier"] and on_calls == []
    assert on.report["empty_replacement_stage_short_circuits"] == 1
    assert on.report["empty_stage_unvisited_role_queries"] == 2
    assert on.report["combinations"] == 0

    # The same known-empty replacement at slot 0 cannot trigger either role lookup.
    on_calls.clear()
    assert on._stage(rule, ("replaced", "later"), "x", replace=(0, ())) == []
    assert on_calls == []
    assert on.report["empty_replacement_stage_short_circuits"] == 2
    assert on.report["empty_stage_unvisited_role_queries"] == 3

    # A slot that legacy zip would never replace must retain the legacy route.
    on_calls.clear()
    assert len(on._stage(rule, ("first", "second"), "x", replace=(0.5, []))) == 1
    assert on_calls == ["first", "second"]
    on_calls.clear()
    assert on._stage(rule, ("earlier", "replaced"), "x", replace=(1, iter(()))) == []
    assert on_calls == ["earlier"]
    assert on.report["empty_replacement_stage_short_circuits"] == 2

    on_calls.clear()
    assert on._stage(rule, ("earlier", "empty"), "x") == []
    assert on_calls == ["earlier", "empty"]
    assert on.report["empty_stage_unvisited_role_queries"] == 3
    assert len(off._stage(rule, (x for x in ("first", "second")), "x")) == 1


def test_nonempty_stage_order_and_config_bool_are_unchanged():
    rule = SimpleNamespace(id="join", roles=("role0", "role1"))
    applied = []
    rules = SimpleNamespace(
        sha256="test",
        apply=lambda _rule, reactants, _limit: (applied.append(reactants) or ("CN",), False),
    )
    catalog = SimpleNamespace(sha256="test")

    def run(enabled):
        search = RouteFirstSearch(
            "CN",
            catalog,
            rules,
            RouteFirstLimits(**limits(skip_empty_stage_retrieval=enabled)),
        )
        calls = []
        search._nearest = lambda target, _role, **_kw: (
            calls.append(target)
            or (
                SimpleNamespace(id=f"{target}-a", smiles="C"),
                SimpleNamespace(id=f"{target}-b", smiles="N"),
            )
        )
        states = search._stage(rule, ("first", "second"), "x")
        return calls, states, search.report

    off_calls, off_states, off_report = run(False)
    on_calls, on_states, on_report = run(True)
    assert off_calls == on_calls == ["first", "second"]
    assert off_states == on_states
    assert off_report["combinations"] == on_report["combinations"] == 4
    assert off_report["products_generated"] == on_report["products_generated"] == 4
    assert on_report["empty_stage_short_circuits"] == 0
    assert applied == [("C", "C"), ("C", "N"), ("N", "C"), ("N", "N")] * 2
    assert RouteFirstLimits(**limits()).skip_empty_stage_retrieval is False
    with pytest.raises(ValueError, match="skip_empty_stage_retrieval must be a boolean"):
        RouteFirstLimits(**limits(skip_empty_stage_retrieval=1))


def test_actual_one_and_two_step_routes_have_complete_forward_and_atom_stereo_replay(assets):
    catalog, rules = assets
    routes, report = generate(TWO_STEP_HIT, catalog, rules, limits())
    assert report["depth_counts"]["1"] > 0
    assert report["depth_counts"]["2"] > 0
    assert report["catalog_row_visits"] <= 200
    assert report["fingerprint_calls"] <= 200
    assert report["total_fingerprint_calls"] <= 200
    assert report["total_fingerprint_calls"] == (
        report["hit_fingerprint_calls"]
        + report["generation_fingerprint_calls"]
        + report["admission_rank_fingerprint_calls"]
        + report["admission_fingerprint_calls"]
    )
    assert report["admission_attempted_products"] <= report["admission_ranked_products"]
    assert (
        report["admission_ranked_products"] + report["admission_fp_unassessed_products"]
        == report["admission_pending_unique_products"]
    )
    assert report["combinations"] <= 100
    assert (
        0
        < report["inner_role_compatible"]
        <= report["inner_role_checks"]
        <= report["products_generated"]
        <= 200
    )
    assert report["route_audit_attempts"] <= report["limits"]["max_routes"]
    assert len(report["phase2_plan_trace"]) <= 8
    assert any(slot["inner_stage_calls"] for plan in report["phase2_plan_trace"] for slot in plan["slots"])
    assert all(
        slot["inner_role_compatible_states"] <= slot["inner_stage_states"]
        and slot["pool_selected_blocks"] <= slot["pool_role_compatible"] + slot["pool_cache_hits"] * 2
        for plan in report["phase2_plan_trace"]
        for slot in plan["slots"]
    )
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


def test_phase_two_executes_typed_guides_one_per_outer_plan_before_later_guides(tmp_path, monkeypatch):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "round_robin_rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {"id": "outer-a", "smarts": "[C:1].[O:2]>>[C:1][O:2]", "source": "test"},
                    {"id": "outer-b", "smarts": "[C:1].[N:2]>>[C:1][N:2]", "source": "test"},
                    {"id": "inner-a", "smarts": "[C:1].[F:2]>>[C:1][F:2]", "source": "test"},
                    {"id": "inner-b", "smarts": "[C:1].[Cl:2]>>[C:1][Cl:2]", "source": "test"},
                ],
            }
        )
    )
    rules = Rules(path)
    settings = RouteFirstLimits(
        **limits(max_catalog_row_visits=10, max_outer_guides=2, max_rules=2, fallback_rules=1)
    )

    def make_search():
        search = RouteFirstSearch("CCO", SimpleNamespace(sha256="test"), rules, settings)
        attempts = []

        def select(target, maximum, fallback=0):
            search.report["rule_inspections"] += 1
            if search._phase == 1:
                return [(rules.get("outer-a"), False), (rules.get("outer-b"), False)]
            return [(rules.get("inner-a"), False), (rules.get("inner-b"), False)]

        def reverse(rule, target, fallback, *, record_centers=False):
            if search._phase == 1:
                guide = ("CC", "O") if rule.id == "outer-a" else ("CCC", "N")
                search._guide_quality[(rule.id, target, guide)] = "typed"
                search._guide_centers[(rule.id, target, guide)] = (None,)
                return [guide]
            guides = [("C", "O"), ("C", "N")] if rule.id == "inner-a" else [("CC", "N")]
            for guide in guides:
                if guide == ("C", "N"):
                    quality = "implicit_fallback"
                elif rule.id == "inner-a":
                    quality = "approximate"
                else:
                    quality = "typed"
                search._guide_quality[(rule.id, target, guide)] = quality
            return guides

        def stage(rule, guide, prefix, **kwargs):
            if prefix == "i" and search.report["catalog_row_visits"] < settings.max_catalog_row_visits:
                search.report["catalog_row_visits"] += 5
                attempts.append((rule.id, tuple(guide)))
            return []

        monkeypatch.setattr(search, "_select_rules", select)
        monkeypatch.setattr(search, "_reverse", reverse)
        monkeypatch.setattr(search, "_stage", stage)
        monkeypatch.setattr(search, "_admit", lambda: None)
        return search, attempts

    search, attempts = make_search()
    _, report = search.run()
    assert attempts == [("inner-b", ("CC", "N")), ("inner-b", ("CC", "N"))]
    assert report["catalog_row_visits"] == settings.max_catalog_row_visits
    plans = report["phase2_plan_trace"]
    assert [row["catalog_rows_execution"] for row in plans] == [5, 5]
    assert [row["rule_inspections_preparation"] for row in plans] == [2, 2]
    assert [row["typed_prepared"] for row in plans] == [2, 2]
    assert [row["approximate_prepared"] for row in plans] == [2, 2]
    assert [row["fallback_prepared"] for row in plans] == [2, 2]
    assert all(row["typed_attempted"] >= 1 for row in plans)
    assert all(row["fallback_attempted"] >= 1 for row in plans)
    assert all(row["catalog_exhausted_attempts"] > 0 for row in plans)

    expired, attempted = make_search()
    original_stage = expired._stage

    def expire_after_first(rule, guide, prefix, **kwargs):
        result = original_stage(rule, guide, prefix, **kwargs)
        if prefix == "i" and attempted:
            expired.deadline = time.monotonic() - 1
        return result

    monkeypatch.setattr(expired, "_stage", expire_after_first)
    _, expired_report = expired.run()
    assert len(attempted) == 1
    assert expired_report["catalog_row_visits"] <= settings.max_catalog_row_visits
    assert any(
        row["work_items_deferred"] and row["defer_reason"] == "deadline_monotonic"
        for row in expired_report["phase2_plan_trace"]
    )

    preparation_expired, preparation_attempts = make_search()
    original_reverse = preparation_expired._reverse

    def expire_during_preparation(rule, target, fallback, *, record_centers=False):
        guides = original_reverse(rule, target, fallback, record_centers=record_centers)
        if preparation_expired._phase == 2:
            preparation_expired.deadline = time.monotonic() - 1
        return guides

    monkeypatch.setattr(preparation_expired, "_reverse", expire_during_preparation)
    _, preparation_report = preparation_expired.run()
    first_plan = preparation_report["phase2_plan_trace"][0]
    assert first_plan["prepared_work_items"] > 0
    assert first_plan["preparation_limit_observed"] == "deadline_monotonic"
    assert first_plan["work_items_deferred"] == first_plan["prepared_work_items"]
    assert not preparation_attempts


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


def test_generation_fingerprint_ceiling_preserves_pending_route_admission(assets, route):
    catalog, rules = assets

    def accept(_product):
        accept.fingerprint_calls += 1

    accept.fingerprint_calls = 0
    search = RouteFirstSearch(
        HIT,
        catalog,
        rules,
        RouteFirstLimits(**limits(max_steps=1, max_fingerprint_calls=8, max_products=1)),
        accept=accept,
    )
    search._pending = [
        {"root": route["root"], "product": route["nodes"][-1]["product"], "nodes": route["nodes"]}
    ]
    assert search.report["admission_fingerprint_reserve"] == 2
    assert search.report["generation_fingerprint_ceiling"] == 6
    assert search.report["depth_one_fingerprint_ceiling"] == 4
    search.report["fingerprint_calls"] = 4
    assert search._fingerprint("CCCC", purpose="generation") is not None
    assert search._fingerprint("CCCCC", purpose="generation") is not None
    assert search._fingerprint("CCCCCC", purpose="generation") is None
    assert "generation_fingerprint_budget" in search.report["limits_reached"]
    assert "depth_one_fingerprint_budget" not in search.report["limits_reached"]
    search._admit()
    assert search.report["admission_pending_unique_products"] == 1
    assert search.report["admission_ranked_products"] == search.report["admission_attempted_products"] == 1
    assert search.report["admission_fp_unassessed_products"] == 0
    assert search.report["admission_callback_unattempted_products"] == 0
    assert search.report["routes_verified"] == 1
    assert search.routes[0]["verification"]["product"] == canonical("CC(=O)Nc1ccncc1")
    assert search.report["fingerprint_calls"] + search.report["admission_fingerprint_calls"] <= 8

    accept.fingerprint_calls = 0
    two_step = RouteFirstSearch(
        HIT,
        catalog,
        rules,
        RouteFirstLimits(**limits(max_steps=2, max_fingerprint_calls=8, max_products=1)),
        accept=accept,
    )
    two_step.report["fingerprint_calls"] = 4
    assert two_step._fingerprint("CCCC", purpose="generation") is None
    assert "depth_one_fingerprint_budget" in two_step.report["limits_reached"]
    two_step._phase = 2
    assert two_step._fingerprint("CCCC", purpose="generation") is not None


def test_tiny_fingerprint_reservation_handles_cached_and_uncached_products(assets):
    catalog, rules = assets
    hit_route = {
        "root": "r0",
        "product": canonical(HIT),
        "nodes": [
            {"id": "b0", "building_block_id": "acetic"},
            {"id": "b1", "building_block_id": "phenyl"},
            {"id": "r0", "rule_id": "amide", "reactants": ["b0", "b1"], "product": canonical(HIT)},
        ],
    }
    uncached = {
        "root": "r0",
        "product": canonical("CC(=O)Nc1ccncc1"),
        "nodes": [
            {"id": "b0", "building_block_id": "acetic"},
            {"id": "b1", "building_block_id": "pyridyl"},
            {
                "id": "r0",
                "rule_id": "amide",
                "reactants": ["b0", "b1"],
                "product": canonical("CC(=O)Nc1ccncc1"),
            },
        ],
    }

    def accept(_product):
        accept.fingerprint_calls += 1

    accept.fingerprint_calls = 0
    search = RouteFirstSearch(
        HIT,
        catalog,
        rules,
        RouteFirstLimits(**limits(max_steps=1, max_fingerprint_calls=2, max_products=2)),
        accept=accept,
    )
    search._pending = [hit_route, uncached]
    search._admit()
    assert search.report["admission_fingerprint_reserve"] == 1
    assert search.report["admission_pending_unique_products"] == 2
    assert search.report["admission_ranked_products"] == search.report["admission_attempted_products"] == 1
    assert search.report["admission_fp_unassessed_products"] == 1
    assert search.report["routes_verified"] == 1
    assert search.routes[0]["verification"]["product"] == canonical(HIT)
    assert search.report["fingerprint_calls"] + search.report["admission_fingerprint_calls"] == 2

    accept.fingerprint_calls = 0
    miss = RouteFirstSearch(
        HIT,
        catalog,
        rules,
        RouteFirstLimits(**limits(max_steps=1, max_fingerprint_calls=2, max_products=2)),
        accept=accept,
    )
    miss._pending = [uncached]
    miss._admit()
    assert miss.report["admission_ranked_products"] == miss.report["admission_attempted_products"] == 0
    assert miss.report["admission_fp_unassessed_products"] == 1
    assert miss.report["routes_verified"] == 0
    assert not miss.report["product_rejections"]  # Unassessed is not a chemical failure.

    accept.fingerprint_calls = 0
    spare_generation = RouteFirstSearch(
        HIT,
        catalog,
        rules,
        RouteFirstLimits(**limits(max_steps=1, max_fingerprint_calls=3, max_products=2)),
        accept=accept,
    )
    spare_generation._pending = [uncached]
    spare_generation._admit()
    assert spare_generation.report["admission_fingerprint_reserve"] == 1
    assert spare_generation.report["admission_reserved_cost_units"] == 2
    assert spare_generation.report["admission_ranked_products"] == 1
    assert spare_generation.report["routes_verified"] == 1
    assert (
        spare_generation.report["fingerprint_calls"] + spare_generation.report["admission_fingerprint_calls"]
        == 3
    )


def test_admission_reservation_deadline_and_callback_counter_contract(assets):
    catalog, rules = assets
    state = {"root": "r0", "product": canonical(HIT), "nodes": []}
    expired = RouteFirstSearch(
        HIT,
        catalog,
        rules,
        RouteFirstLimits(**limits(max_steps=1, max_fingerprint_calls=2)),
        deadline_monotonic=time.monotonic() - 1,
    )
    expired._pending = [state]
    expired._admit()
    assert expired.report["admission_pending_unique_products"] == 1
    assert expired.report["admission_fp_unassessed_products"] == 1
    assert expired.report["admission_attempted_products"] == 0
    assert expired.report["admission_callback_unattempted_products"] == 0
    assert "deadline_monotonic" in expired.report["limits_reached"]

    def deadline_after_first(_product):
        deadline_after_first.fingerprint_calls += 1
        partially_ranked.deadline = time.monotonic() - 1

    deadline_after_first.fingerprint_calls = 0
    partially_ranked = RouteFirstSearch(
        HIT,
        catalog,
        rules,
        RouteFirstLimits(**limits(max_steps=1, max_fingerprint_calls=4)),
        accept=deadline_after_first,
    )
    partially_ranked._pending = [state, {**state, "product": canonical("CC(=O)Nc1ccncc1")}]
    partially_ranked._admit()
    assert partially_ranked.report["admission_ranked_products"] == 2
    assert partially_ranked.report["admission_attempted_products"] == 1
    assert partially_ranked.report["admission_callback_unattempted_products"] == 1
    assert partially_ranked.report["admission_fp_unassessed_products"] == 0

    def broken_accept(_product):
        broken_accept.fingerprint_calls += 2

    broken_accept.fingerprint_calls = 0
    broken = RouteFirstSearch(
        HIT,
        catalog,
        rules,
        RouteFirstLimits(**limits(max_steps=1, max_fingerprint_calls=4)),
        accept=broken_accept,
    )
    broken._pending = [state]
    with pytest.raises(RuntimeError, match="one-fingerprint contract"):
        broken._admit()
    _, capped = generate(TWO_STEP_HIT, catalog, rules, limits(max_routes=2))
    assert capped["routes_verified"] <= 2
    assert capped["route_audit_attempts"] <= 2
    assert capped["depth_counts"]["2"] >= 1


def test_catalog_source_stereo_screen_fills_top_k_and_unknown_keeps_legacy_rank(tmp_path, monkeypatch):
    from chembounce_synthesis import route_first
    from chembounce_synthesis.catalog import Catalog

    path = tmp_path / "stereo-catalog.tsv"
    path.write_text(
        "id\tsmiles\tsource\tsnapshot\n"
        "bad\tCC(F)Br\ttest\t2026\n"
        "specified\tC[C@H](F)Br\ttest\t2026\n"
        "achiral\tCCBr\ttest\t2026\n"
    )
    catalog = Catalog(path)
    settings = RouteFirstLimits(
        **limits(
            max_steps=1,
            max_catalog_index_rows=10,
            max_catalog_rows_per_query=20,
            max_catalog_row_visits=200,
        )
    )

    def search():
        return RouteFirstSearch("CC(F)Br", catalog, SimpleNamespace(sha256="test"), settings)

    screened = search()
    assert [block.id for block in screened._nearest("CC(F)Br", "[C]", rule_id="r", role_slot=0)] == [
        "specified",
        "achiral",
    ]
    assert screened.report["catalog_stereo_exclusions"] == 1
    assert screened.report["catalog_stereo_exclusion_samples"][0]["building_block_id"] == "bad"
    assert screened.report["catalog_stereo_perceptions"] == 3
    assert screened.report["catalog_stereo_atom_bond_charges"] <= (
        route_first.MAX_CATALOG_STEREO_ATOM_BOND_CHARGES
    )
    assert len(screened._catalog_stereo_cache) <= route_first.MAX_CATALOG_STEREO_CACHE_ENTRIES
    before = screened.report["catalog_stereo_perceptions"]
    assert screened._nearest("CC(F)Br", "[C]", rule_id="different", role_slot=1)
    assert screened.report["catalog_stereo_perceptions"] == before  # unchanged nearest key/cache

    monkeypatch.setattr(route_first, "MAX_CATALOG_STEREO_PERCEPTIONS", 0)
    capped = search()
    assert [block.id for block in capped._nearest("CC(F)Br", "[C]")] == ["bad", "specified"]
    assert capped.report["catalog_stereo_unknowns"]["perception_budget"] >= 1
    assert capped._nearest_cache == {}  # no incomplete selection cached
    monkeypatch.setattr(route_first, "MAX_CATALOG_STEREO_PERCEPTIONS", 2048)

    monkeypatch.setattr(route_first, "MAX_CATALOG_STEREO_ATOM_BOND_CHARGES", 0)
    size_capped = search()
    assert [block.id for block in size_capped._nearest("CC(F)Br", "[C]")] == [
        "bad",
        "specified",
    ]
    assert size_capped.report["catalog_stereo_unknowns"]["atom_bond_budget"] >= 1
    assert size_capped.report["catalog_stereo_perceptions"] == 0
    assert size_capped.report["catalog_stereo_atom_bond_charges"] == 0
    assert size_capped._nearest_cache == {}
    monkeypatch.setattr(route_first, "MAX_CATALOG_STEREO_ATOM_BOND_CHARGES", 131072)

    monkeypatch.setattr(route_first, "MAX_CATALOG_STEREO_CACHE_ENTRIES", 0)
    cache_capped = search()
    assert [block.id for block in cache_capped._nearest("CC(F)Br", "[C]")] == [
        "bad",
        "specified",
    ]
    assert cache_capped.report["catalog_stereo_unknowns"]["cache_budget"] >= 1
    assert cache_capped.report["catalog_stereo_perceptions"] == 0
    assert cache_capped._nearest_cache == {}
    monkeypatch.setattr(route_first, "MAX_CATALOG_STEREO_CACHE_ENTRIES", 2048)

    def failed_perception(_mol):
        raise ValueError("synthetic perception failure")

    monkeypatch.setattr(route_first, "stereo_state", failed_perception)
    failed = search()
    assert [block.id for block in failed._nearest("CC(F)Br", "[C]")] == ["bad", "specified"]
    assert failed.report["catalog_stereo_unknowns"]["perception_failed"] >= 1
    assert failed._nearest_cache == {}
    assert failed.report["catalog_stereo_perceptions"] <= route_first.MAX_CATALOG_STEREO_PERCEPTIONS

    expired = search()
    expired.deadline = time.monotonic() - 1
    assert expired._catalog_block_stereo_eligible("CC(F)Br") is None
    assert expired.report["catalog_stereo_perceptions"] == 0

    source_sulfur = "O=S(=O)(Cl)N=S(=O)(c1ccc(F)cc1)C(F)(F)F"
    sulfur_catalog = SimpleNamespace(sha256="test", mols={source_sulfur: molecule(source_sulfur)})
    sulfur = RouteFirstSearch(source_sulfur, sulfur_catalog, SimpleNamespace(sha256="test"), settings)
    monkeypatch.undo()
    assert sulfur._catalog_block_stereo_eligible(source_sulfur) is False


def test_stereo_screen_preserves_one_anchor_one_ordinary_reservation(tmp_path):
    from rdkit import Chem

    from chembounce_synthesis.catalog import Catalog

    path = tmp_path / "anchor-ordinary-stereo.tsv"
    path.write_text(
        "id\tsmiles\tsource\tsnapshot\n"
        "a0\tCCBr\ttest\t2026\n"
        "a1\tCCCBr\ttest\t2026\n"
        "o0\tCC(F)Br\ttest\t2026\n"
        "o1\tC[C@H](F)Br\ttest\t2026\n"
    )
    catalog = Catalog(path)
    settings = RouteFirstLimits(
        **limits(
            max_steps=1,
            max_catalog_index_rows=10,
            max_catalog_rows_per_query=20,
            max_catalog_row_visits=200,
        )
    )
    search = RouteFirstSearch("CC(F)Br", catalog, SimpleNamespace(sha256="test"), settings)
    search._index()
    anchors = (canonical("CCBr"), canonical("CCCBr"))
    search._anchor_queries = (Chem.MolFromSmarts("[C]"),)
    search._anchor_postings = (anchors,)
    search._anchor_scores = {smiles: index for index, smiles in enumerate(anchors)}
    selected = search._nearest("CC(F)Br", "[C]")
    assert [block.id for block in selected] == ["a0", "o1"]
    assert search.report["catalog_stereo_exclusions"] == 1


def _fallback_second_search(tmp_path, *, hit="OC1CCCCC1"):
    from chembounce_synthesis.catalog import Catalog

    path = tmp_path / "fallback-second-catalog.tsv"
    path.write_text(
        "id\tsmiles\tsource\tsnapshot\n"
        "alkyne\tC#CCO\ttest\ttest\n"
        "branched\tCC(C)O\ttest\ttest\n"
        "linear\tCCCO\ttest\ttest\n"
        "ring\tOC1CCCCC1\ttest\ttest\n"
        "unspecified\tCC(F)O\ttest\ttest\n"
    )
    catalog = Catalog(path)
    return RouteFirstSearch(
        hit,
        catalog,
        SimpleNamespace(sha256="test"),
        RouteFirstLimits(**limits(max_steps=1, max_catalog_index_rows=10)),
    )


def test_explicit_missing_role_second_uses_hit_context_without_changing_first_or_fp_work(tmp_path):
    search = _fallback_second_search(tmp_path)
    old = search._nearest("C", "[C][O]")
    fp_before = search.report["fingerprint_calls"]
    assert [block.id for block in old] == ["alkyne", "branched"]
    changed = search._nearest("C", "[C][O]", fallback_missing_role=True)
    assert [block.id for block in changed] == ["alkyne", "ring"]
    assert search.report["fingerprint_calls"] == fp_before
    assert search.report["fallback_hit_second_replacements"] == 1
    assert search.report["fallback_hit_second_tanimotos"] <= 5
    assert (
        search.report["fallback_hit_second_tanimotos"] <= route_first_module.MAX_FALLBACK_HIT_SECOND_TANIMOTOS
    )
    sample = search.report["fallback_hit_second_samples"][0]
    assert (sample["first"], sample["old_second"], sample["new_second"]) == (
        canonical("C#CCO"),
        canonical("CC(C)O"),
        canonical("OC1CCCCC1"),
    )
    assert sample["comparison_count"] == search.report["fallback_hit_second_tanimotos"]
    assert sample["new_second_hit_tanimoto"] > sample["old_second_hit_tanimoto"]
    assert search._nearest("C", "[C][O]") == old
    assert search._nearest("C", "[C][O]", fallback_missing_role=True) == changed
    assert len(search._nearest_cache) == 2  # active and ordinary contexts cannot alias

    anchored = _fallback_second_search(tmp_path)
    anchored._index()
    anchored._anchor_queries = (Chem.MolFromSmarts("[C]"),)
    anchored._anchor_postings = ((canonical("C#CCO"),),)
    anchored._anchor_scores = {canonical("C#CCO"): 0}
    assert [block.id for block in anchored._nearest("C", "[C][O]", fallback_missing_role=True)] == [
        "alkyne",
        "ring",
    ]
    assert anchored.report["fallback_hit_second_replacements"] == 1

    unspecified = _fallback_second_search(tmp_path, hit="CC(F)O")
    assert [block.id for block in unspecified._nearest("C", "[C][O]", fallback_missing_role=True)] == [
        "alkyne",
        "branched",
    ]
    assert unspecified.report["catalog_stereo_exclusions"] >= 1


def test_explicit_fallback_provenance_marks_only_missing_role_and_not_implicit(tmp_path):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "fallback-rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "join", "smarts": "[C:1].[N:2]>>[C:1][N:2]", "source": "test"}],
            }
        )
    )
    rules = Rules(path)
    hit = "CCO"
    explicit = RouteFirstSearch(
        hit, SimpleNamespace(sha256="test"), rules, RouteFirstLimits(**limits(max_steps=1))
    )
    guide = explicit._reverse(rules.get("join"), hit, True)[0]
    assert guide == (hit, "C")
    assert explicit._guide_missing_role_slots[("join", hit, guide, True)] == (1,)
    calls = []
    explicit._nearest = lambda target, role, **kwargs: calls.append((target, kwargs)) or ()
    explicit._stage(rules.get("join"), guide, "o", fallback_missing_slots=(1,))
    assert [call[1]["fallback_missing_role"] for call in calls] == [False, True]
    calls.clear()
    explicit._stage(rules.get("join"), guide, "i")
    assert [call[1]["fallback_missing_role"] for call in calls] == [False, False]

    implicit = RouteFirstSearch(
        hit, SimpleNamespace(sha256="test"), rules, RouteFirstLimits(**limits(max_steps=1))
    )
    implicit._topology_spec = lambda _rule: (None, "unsupported_for_test")
    assert implicit._reverse(rules.get("join"), hit, False)[0] == guide
    assert implicit._guide_missing_role_slots == {}
    assert implicit._guide_quality[("join", hit, guide)] == "implicit_fallback"

    plan_only = RouteFirstSearch(
        hit, SimpleNamespace(sha256="test"), rules, RouteFirstLimits(**limits(max_steps=1))
    )
    propagated = []
    plan_only._stage = lambda *_args, **kwargs: propagated.append(kwargs["fallback_missing_slots"]) or []
    _routes, report = plan_only.run()
    assert report["outer_plans_selected"][0]["hint_quality"] == "explicit_fallback"
    assert propagated == [(1,)]


def test_hit_context_second_budget_missing_cache_unknown_and_deadline_restore_old_second(
    tmp_path, monkeypatch
):
    from chembounce_synthesis import route_first

    baseline = _fallback_second_search(tmp_path)
    old = baseline._nearest("C", "[C][O]")

    monkeypatch.setattr(route_first, "MAX_FALLBACK_HIT_SECOND_TANIMOTOS", 0)
    bounded = _fallback_second_search(tmp_path)
    assert bounded._nearest("C", "[C][O]", fallback_missing_role=True) == old
    assert bounded.report["fallback_hit_second_tanimotos"] == 0
    assert bounded.report["fallback_hit_second_legacy_reasons"]["tanimoto_budget"] == 1
    assert ("C", "[C][O]", True) not in bounded._nearest_cache
    monkeypatch.setattr(route_first, "MAX_FALLBACK_HIT_SECOND_TANIMOTOS", 65536)

    partial = _fallback_second_search(tmp_path)
    original_pool = partial._pool

    def incomplete_pool(*args, **kwargs):
        pool = original_pool(*args, **kwargs)
        partial._pool_last_complete = False
        return pool

    partial._pool = incomplete_pool
    assert partial._nearest("C", "[C][O]", fallback_missing_role=True) == old
    assert partial.report["fallback_hit_second_legacy_reasons"]["incomplete_original_selection"] == 1
    assert ("C", "[C][O]", True) not in partial._nearest_cache

    monkeypatch.setattr(route_first, "MAX_CATALOG_STEREO_PERCEPTIONS", 2)
    stereo_capped = _fallback_second_search(tmp_path)
    assert stereo_capped._nearest("C", "[C][O]", fallback_missing_role=True) == old
    assert stereo_capped.report["fallback_hit_second_legacy_reasons"]["stereo_budget"] == 1
    assert ("C", "[C][O]", True) not in stereo_capped._nearest_cache
    monkeypatch.setattr(route_first, "MAX_CATALOG_STEREO_PERCEPTIONS", 2048)

    missing = _fallback_second_search(tmp_path)
    source_screen = missing._catalog_block_stereo_eligible

    def evict_after_second(smiles):
        result = source_screen(smiles)
        if smiles == canonical("CC(C)O"):
            missing._fp_cache.pop(canonical("OC1CCCCC1"), None)
        return result

    missing._catalog_block_stereo_eligible = evict_after_second
    assert missing._nearest("C", "[C][O]", fallback_missing_role=True) == old
    assert missing.report["fallback_hit_second_legacy_reasons"]["missing_cached_fingerprint"] == 1
    assert ("C", "[C][O]", True) not in missing._nearest_cache

    unknown = _fallback_second_search(tmp_path)
    original = unknown._catalog_block_stereo_eligible

    def uncertain_top(smiles):
        if smiles == canonical("OC1CCCCC1"):
            return None
        return original(smiles)

    unknown._catalog_block_stereo_eligible = uncertain_top
    assert unknown._nearest("C", "[C][O]", fallback_missing_role=True) == old
    assert unknown.report["fallback_hit_second_legacy_reasons"]["unknown_source_stereo"] == 1
    assert ("C", "[C][O]", True) not in unknown._nearest_cache

    deadline = _fallback_second_search(tmp_path)
    similarity = route_first.DataStructs.TanimotoSimilarity

    def expire_during_new_scores(reference, candidate):
        score = similarity(reference, candidate)
        if reference is deadline.hit_fp:
            deadline.deadline = time.monotonic() - 1
        return score

    monkeypatch.setattr(route_first.DataStructs, "TanimotoSimilarity", expire_during_new_scores)
    assert deadline._nearest("C", "[C][O]", fallback_missing_role=True) == old
    assert deadline.report["fallback_hit_second_legacy_reasons"]["deadline_during_scoring"] == 1
    assert ("C", "[C][O]", True) not in deadline._nearest_cache

    monkeypatch.setattr(route_first, "MAX_FALLBACK_HIT_SECOND_SAMPLES", 0)
    no_samples = _fallback_second_search(tmp_path)
    assert [block.id for block in no_samples._nearest("C", "[C][O]", fallback_missing_role=True)] == [
        "alkyne",
        "ring",
    ]
    assert no_samples.report["fallback_hit_second_samples"] == []
    assert no_samples.report["fallback_hit_second_samples_omitted"] == 1


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
    # Exercise the legacy raw partition guard when the optional local query
    # construction is unavailable; the default nonring query filters this cut.
    overlap._nonring_unavailable_rules.add("ether")
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


def test_topology_support_is_conservative_about_lost_atoms_and_multi_bond_cuts():
    from rdkit.Chem import rdChemReactions

    cases = (
        ("[C:1][O:2].[N:3]>>[C:1][N:3]", 1, None),
        ("[C:1][O:2].[N:3]>>[C:1][N:3][O:2]", None, "retained_reactant_edge_deleted"),
        ("[C:1][C:2].[N:3]>>[C:1]1[C:2][N:3]1", 2, None),
        ("[C:1]>>[C:1]", None, "not_two_reactants_one_product"),
        ("[C:1].[N:1]>>[C:1]", None, "ambiguous_or_unbound_maps"),
    )
    for smarts, expected_cross, expected_reason in cases:
        specification, reason = topology_rule_specification(rdChemReactions.ReactionFromSmarts(smarts))
        assert reason == expected_reason
        assert (len(specification[2]) if specification else None) == expected_cross


def test_inner_topology_prunes_only_complete_impossibility_and_retains_unknowns(tmp_path, monkeypatch):
    from chembounce_synthesis import route_first
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "topology_rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {"id": "split", "smarts": "[C:1].[C:2]>>[C:1][C:2]", "source": "test"},
                    {"id": "unary", "smarts": "[C:1]>>[C:1]", "source": "test"},
                ],
            }
        )
    )
    rules = Rules(path)
    catalog = SimpleNamespace(sha256="test")
    settings = RouteFirstLimits(**limits(max_steps=1))
    ring = RouteFirstSearch("C1CC1", catalog, rules, settings)
    ring._phase = 2
    assert [rule.id for rule, _ in ring._select_rules("C1CC1", 2)] == ["unary"]
    assert ring.report["topology_pruned_rules"] == 1
    assert ring.report["topology_match_requests"] == 1
    assert ring.report["rule_inspections"] <= settings.max_rule_inspections
    charged = ring.report["topology_spec_graph_visits"]
    assert charged <= route_first.MAX_TOPOLOGY_SPEC_GRAPH_VISITS
    assert ring._topology_assess("C1CC1", molecule("C1CC1"), rules.get("split")) == "impossible"
    assert ring.report["topology_match_requests"] == 1
    assert ring.report["topology_spec_graph_visits"] == charged
    assert ring.report["topology_assessment_cache_hits"] == 1
    assert len(ring._topology_assessment_cache) <= route_first.MAX_TOPOLOGY_MATCH_REQUESTS

    # The 65th match is a sentinel; a truncated enumeration cannot prove all
    # cuts impossible even when early ring matches are impossible.
    long_chain = "C" * 70
    truncated = RouteFirstSearch(long_chain, catalog, rules, settings)
    assert truncated._topology_assess(long_chain, molecule(long_chain), rules.get("split")) == "unknown"
    assert truncated.report["topology_unknown_reasons"] == {"match_truncation": 1}
    assert truncated.report["topology_match_requests"] == 1

    monkeypatch.setattr(route_first, "MAX_TOPOLOGY_MATCH_REQUESTS", 0)
    capped = RouteFirstSearch("C1CC1", catalog, rules, settings)
    capped._phase = 2
    assert [rule.id for rule, _ in capped._select_rules("C1CC1", 2)] == ["split", "unary"]
    assert capped.report["topology_unknown_reasons"]["match_request_budget"] == 1
    assert capped.report["topology_match_requests"] == 0
    assert len(capped._topology_assessment_cache) == 0
    monkeypatch.setattr(route_first, "MAX_TOPOLOGY_MATCH_REQUESTS", 128)
    monkeypatch.setattr(route_first, "MAX_TOPOLOGY_SPEC_GRAPH_VISITS", 0)
    graph_capped = RouteFirstSearch("C1CC1", catalog, rules, settings)
    graph_capped._phase = 2
    assert [rule.id for rule, _ in graph_capped._select_rules("C1CC1", 2)] == ["split", "unary"]
    assert graph_capped.report["topology_unknown_reasons"]["spec_budget_or_deadline"] == 2
    assert graph_capped.report["topology_spec_graph_visits"] == 0
    expired = RouteFirstSearch("C1CC1", catalog, rules, settings, deadline_monotonic=time.monotonic() - 1)
    assert expired._topology_assess("C1CC1", molecule("C1CC1"), rules.get("split")) == "unknown"
    assert expired.report["topology_match_requests"] == 0


def test_single_cross_edge_local_reverse_suppresses_implicit_resurrection(tmp_path):
    from rdkit.Chem import rdChemReactions

    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "split.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "split", "smarts": "[C:1].[C:2]>>[C:1][C:2]", "source": "test"}],
            }
        )
    )
    rules = Rules(path)
    rule = rules.get("split")
    original = rdChemReactions.ReactionToSmarts(rule.reaction)
    settings = RouteFirstLimits(**limits(max_steps=2))
    ring = RouteFirstSearch("C1CC1", SimpleNamespace(sha256="test"), rules, settings)
    ring._phase = 2
    assert ring._reverse(rule, "C1CC1", False) == []
    assert ring.report["reverse_attempts"] == ring.report["reverse_nonring_query_calls"] == 1
    assert ring.report["reverse_nonring_empty"] == 1
    assert ring.report["reverse_fallbacks"] == 0
    assert ring._reverse(rule, "C1CC1", True)
    assert ring.report["reverse_hint_quality_counts"]["explicit_fallback"] == 1
    assert rdChemReactions.ReactionToSmarts(rule.reaction) == original

    multi_path = tmp_path / "multi.json"
    multi_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "multi",
                        "smarts": "[C:1][C:2].[N:3]>>[C:1]1[C:2][N:3]1",
                        "source": "test",
                    }
                ],
            }
        )
    )
    multi_rules = Rules(multi_path)
    multi = RouteFirstSearch("C1CN1", SimpleNamespace(sha256="test"), multi_rules, settings)
    multi._phase = 2
    multi._reverse(multi_rules.get("multi"), "C1CN1", False)
    assert multi.report["reverse_nonring_query_calls"] == 0

    linear = RouteFirstSearch("CC", SimpleNamespace(sha256="test"), rules, settings)
    linear._phase = 2
    guides = linear._reverse(rule, "CC", False)
    assert guides and 1 <= linear.report["reverse_partition_checks"] <= settings.max_reverse_outcomes
    before = linear.report["reverse_attempts"]
    linear._phase = 1
    assert linear._reverse(rule, "CC", False) == guides
    assert linear.report["reverse_attempts"] == before
    assert linear.report["reverse_nonring_query_compiles"] == 1
    assert linear.report["reverse_nonring_query_calls"] == 1
    assert linear.report["reverse_outcomes"] <= settings.max_reverse_outcomes
    assert rdChemReactions.ReactionToSmarts(rule.reaction) == original

    outer_first = RouteFirstSearch("CC", SimpleNamespace(sha256="test"), rules, settings)
    first_guides = outer_first._reverse(rule, "CC", False, record_centers=True)
    assert first_guides == guides
    assert outer_first.report["reverse_nonring_query_calls"] == 1
    outer_first._phase = 2
    assert outer_first._reverse(rule, "CC", False) == first_guides
    assert outer_first.report["reverse_attempts"] == 1
    assert outer_first.report["reverse_nonring_query_calls"] == 1


def test_reverse_deadline_expiry_during_static_check_or_local_compile_skips_chemistry(tmp_path, monkeypatch):
    from chembounce_synthesis import route_first
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "split", "smarts": "[C:1].[C:2]>>[C:1][C:2]", "source": "test"}],
            }
        )
    )
    rules = Rules(path)
    rule = rules.get("split")
    settings = RouteFirstLimits(**limits(max_steps=2))
    static = RouteFirstSearch("CC", SimpleNamespace(sha256="test"), rules, settings, time.monotonic() + 30)

    def expire_static(_rule):
        static.deadline = time.monotonic() - 1
        return None, "spec_budget_or_deadline"

    monkeypatch.setattr(static, "_topology_spec", expire_static)
    assert static._reverse(rule, "CC", False) == []
    assert static.report["reverse_attempts"] == 1
    assert static.report["reverse_partition_checks"] == 0
    assert static.report["reverse_fallbacks"] == 0
    assert static.report["reverse_local_template_compiles"] == 0

    compiled = RouteFirstSearch("CC", SimpleNamespace(sha256="test"), rules, settings, time.monotonic() + 30)
    original_compile = route_first.tagged_reverse_reaction

    def expire_compile(rule):
        result = original_compile(rule)
        compiled.deadline = time.monotonic() - 1
        return result

    monkeypatch.setattr(route_first, "tagged_reverse_reaction", expire_compile)
    assert compiled._reverse(rule, "CC", False) == []
    assert compiled.report["reverse_attempts"] == 1
    assert compiled.report["reverse_nonring_query_calls"] == 0
    assert compiled.report["reverse_partition_checks"] == 0
    assert compiled.report["reverse_fallbacks"] == 0


def test_optional_outer_center_metadata_never_repeats_cached_reverse_work(tmp_path):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [{"id": "split", "smarts": "[C:1].[C:2]>>[C:1][C:2]", "source": "test"}],
            }
        )
    )
    rules = Rules(path)
    rule = rules.get("split")
    search = RouteFirstSearch("CC", SimpleNamespace(sha256="test"), rules, RouteFirstLimits(**limits()))
    guides = search._reverse(rule, search.hit, False, record_centers=True)
    before = (search.report["reverse_attempts"], search.report["reverse_role_trials"])
    assert search._reverse(rule, search.hit, False) == guides
    assert (search.report["reverse_attempts"], search.report["reverse_role_trials"]) == before
    assert search._guide_centers[(rule.id, search.hit, guides[0])] == (((0, 1),),)

    inner_first = RouteFirstSearch("CC", SimpleNamespace(sha256="test"), rules, RouteFirstLimits(**limits()))
    guides = inner_first._reverse(rule, inner_first.hit, False)
    before = (inner_first.report["reverse_attempts"], inner_first.report["reverse_role_trials"])
    assert inner_first._reverse(rule, inner_first.hit, False, record_centers=True) == guides
    assert (inner_first.report["reverse_attempts"], inner_first.report["reverse_role_trials"]) == before
    assert inner_first._guide_centers[(rule.id, inner_first.hit, guides[0])] == (None,)


def test_outer_plan_cut_diversity_and_approximate_fallback_reachability(tmp_path, monkeypatch):
    import chembounce_synthesis.route_first as route_first
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {"id": "split", "smarts": "[C:1].[C:2]>>[C:1][C:2]", "source": "test"},
                    {"id": "unary", "smarts": "[C:1][C:2]>>[C:1][C:2]", "source": "test"},
                    {"id": "fallback", "smarts": "[C:1].[O:2]>>[C:1]=[O:2]", "source": "test"},
                ],
            }
        )
    )
    rules = Rules(path)

    def plan_only():
        search = RouteFirstSearch(
            "CC",
            SimpleNamespace(sha256="test"),
            rules,
            RouteFirstLimits(**limits(max_steps=1, max_rules=3, fallback_rules=1, max_outer_guides=3)),
        )
        search._stage = lambda *_args, **_kwargs: []
        return search.run()[1]

    report = plan_only()
    assert [(row["rule_id"], row["center_kind"]) for row in report["outer_plans_selected"]] == [
        ("split", "cut"),
        ("unary", "no_cut"),
        ("fallback", "explicit_fallback"),
    ]
    assert len({tuple(row["guide"]) for row in report["outer_plans_selected"]}) == 3
    assert report["outer_plans_selected"][0]["cut_bond_centers"] == [[[0, 1]]]

    monkeypatch.setattr(route_first, "MAX_OUTER_CENTER_METADATA_CHECKS", 0)
    unknown = plan_only()
    assert any(row["center_kind"] == "unknown" for row in unknown["outer_plans_selected"])
    assert unknown["outer_center_metadata_incomplete"] >= 1
    assert len(unknown["outer_plans_selected"]) == 3

    approximate_path = tmp_path / "approximate.json"
    approximate_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {"id": "aryl-h2", "smarts": "[Cl,Br,I][c:1].[N;H2:2]>>[c:1][N:2]", "source": "test"}
                ],
            }
        )
    )
    approximate_rules = Rules(approximate_path)
    approximate = RouteFirstSearch(
        "c1ccc(N2CCC[C@H](C)C2)cc1",
        SimpleNamespace(sha256="test"),
        approximate_rules,
        RouteFirstLimits(**limits(max_steps=1, max_outer_guides=1)),
    )
    approximate._stage = lambda *_args, **_kwargs: []
    _, retained = approximate.run()
    assert retained["outer_plans_selected"][0]["hint_quality"] == "approximate"


def test_two_distinct_cut_centers_from_one_rule_keep_two_outer_slots(tmp_path):
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
    search = RouteFirstSearch(
        "Nc1ccc(Nc2ccccc2)cc1",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits(max_steps=1, max_outer_guides=2)),
    )
    search._stage = lambda *_args, **_kwargs: []
    _, report = search.run()
    plans = report["outer_plans_selected"]
    assert len(plans) == 2
    assert all(plan["rule_id"] == "aryl" and plan["hint_quality"] == "typed" for plan in plans)
    assert plans[0]["cut_bond_centers"] != plans[1]["cut_bond_centers"]
    assert len({tuple(plan["guide"]) for plan in plans}) == 2


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


def test_reverse_hint_restores_only_proven_cleaved_ring_query_flags(tmp_path):
    from chembounce_synthesis.reactions import Rules
    from chembounce_synthesis.route_first import tagged_reverse_reaction

    path = tmp_path / "ring-rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "close-ring",
                        "smarts": "[C:1][C:2][C:3][C:4][C:5][C:6]>>[c:1]1[c:2][c:3][c:4][c:5][c:6]1",
                        "source": "test",
                    }
                ],
            }
        )
    )
    rule = Rules(path).get("close-ring")
    parent = molecule("c1ccccc1")
    raw = molecule("CCCCCC")
    for index, atom in enumerate(raw.GetAtoms()):
        atom.SetIsAromatic(True)  # RDKit reverse-query metadata after opening the ring.
        atom.SetIntProp("old_mapno", index + 1)
        atom.SetIntProp("react_idx", 0)
        atom.SetIntProp("react_atom_idx", index)
    original = Chem.Mol(raw).ToBinary(Chem.PropertyPickleOptions.AllProps)
    normalized, edits = normalize_reverse_hint(raw, parent, rule)
    assert normalized == "CCCCCC"
    assert edits["reset_product_aromatic_query_flags"] == 6
    _, tags = tagged_reverse_reaction(rule)
    assert typed_reverse_hint(raw, parent, rule, tags, 0)[0] == normalized
    assert raw.ToBinary(Chem.PropertyPickleOptions.AllProps) == original

    # An intact aromatic ring keeps its original flags and bonds.
    closed = molecule("c1ccccc1")
    closed.GetAtomWithIdx(0).SetIntProp("old_mapno", 1)
    closed.GetAtomWithIdx(0).SetIntProp("react_idx", 0)
    closed.GetAtomWithIdx(0).SetIntProp("react_atom_idx", 0)
    kept, edits = normalize_reverse_hint(closed, parent, rule)
    assert kept == canonical("c1ccccc1")
    assert edits["reset_product_aromatic_query_flags"] == 0

    # Missing origin metadata and a surviving aromatic bond are not repairs.
    ambiguous = Chem.Mol(raw)
    ambiguous.GetAtomWithIdx(0).ClearProp("react_atom_idx")
    with pytest.raises(ReverseHintError, match="hint_sanitization_failed"):
        normalize_reverse_hint(ambiguous, parent, rule)
    aromatic_bond = Chem.Mol(raw)
    bond = aromatic_bond.GetBondBetweenAtoms(0, 1)
    bond.SetBondType(Chem.BondType.AROMATIC)
    bond.SetIsAromatic(False)
    with pytest.raises(ReverseHintError, match="hint_sanitization_failed"):
        normalize_reverse_hint(aromatic_bond, parent, rule)


def test_reverse_hint_query_proof_rejects_generic_and_mixed_atomic_numbers():
    from chembounce_synthesis.route_first import _role_forces_aliphatic_element

    def proves(smarts):
        return _role_forces_aliphatic_element(Chem.MolFromSmarts(smarts).GetAtomWithIdx(0), 6)

    assert proves("[C]")
    assert proves("[C;!N]")
    assert not proves("[#6]")  # Atomic number also admits aromatic carbon.
    assert not proves("[#6,C]")
    assert not proves("[!C]")
    assert not proves("[c]")


def test_reverse_hint_cleaved_pyrazole_has_typed_roles_only_when_both_match(tmp_path):
    from rdkit import rdBase

    from chembounce_synthesis.reactions import Rules
    from chembounce_synthesis.route_first import tagged_reverse_reaction

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "pyrazole",
                        "smarts": (
                            "[NH2,NH3+:3]-[N$([NH](N)([CX4])):2]."
                            "[C$([CH](C)(C)([CX4])),C$([CH2](C)(C)):6]"
                            "(-[C$(C(=O)(C)([CX4])),C$([CH](=O)(C)):5]=[O:9])"
                            "-[C$(C(=O)(C)([CX4])),C$([CH](=O)(C)):7]=[O:10]"
                            ">>[c:7]1[n:3][n:2][c:5][c:6]1"
                        ),
                        "source": "test",
                    }
                ],
            }
        )
    )
    rules = Rules(path)
    rule = rules.get("pyrazole")
    with rdBase.BlockLogs():
        products, truncated = rules.apply("pyrazole", ("CNN", "CC(=O)CC(C)=O"), product_limit=4)
    assert products and not truncated
    parent = molecule(products[0])
    reverse, tags = tagged_reverse_reaction(rule)
    with rdBase.BlockLogs():
        outcomes = reverse.RunReactants((parent,), maxProducts=4)
    assert outcomes and len(outcomes) <= 4
    typed_pairs = []
    for outcome in outcomes:
        if reverse_partition_issue(outcome, parent.GetNumAtoms()) is not None:
            continue
        try:
            typed_pairs.append(
                tuple(
                    typed_reverse_hint(part, parent, rule, tags, slot)[0] for slot, part in enumerate(outcome)
                )
            )
        except ReverseHintError:
            continue
    assert ("CN[NH2+]", "CC(=O)CC(C)=O") in typed_pairs
    search = RouteFirstSearch(
        products[0], SimpleNamespace(sha256="test"), rules, RouteFirstLimits(**limits())
    )
    counted = [search._typed_reverse_outcome(outcome, parent, rule, tags) for outcome in outcomes]
    assert any(
        item is not None and item[1]["typed_reset_product_aromatic_query_flags"] == 4 for item in counted
    )


def test_reverse_hint_actual_role_still_rejects_wrong_mapped_site(tmp_path):
    from chembounce_synthesis.reactions import Rules
    from chembounce_synthesis.route_first import tagged_reverse_reaction

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "mapped-ring",
                        "smarts": "[N:1][C:2][C:3][C:4][C:5][C:6]>>[n:1]1[c:2][c:3][c:4][c:5][c:6]1",
                        "source": "test",
                    }
                ],
            }
        )
    )
    rule = Rules(path).get("mapped-ring")
    parent = molecule("c1ccncc1")
    raw = molecule("CCCCCN")
    for index, atom in enumerate(raw.GetAtoms()):
        atom.SetIsAromatic(True)
        atom.SetIntProp("old_mapno", index + 1)
        atom.SetIntProp("react_idx", 0)
        atom.SetIntProp("react_atom_idx", index)
    # Source/map element mismatch cannot be cleared into a false typed role.
    with pytest.raises(ReverseHintError, match="hint_sanitization_failed"):
        normalize_reverse_hint(raw, parent, rule)
    _, tags = tagged_reverse_reaction(rule)
    with pytest.raises(ReverseHintError, match="hint_sanitization_failed"):
        typed_reverse_hint(raw, parent, rule, tags, 0)


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

    def late_choice(
        _raw, _parent, _rule, _tags, slot, representative=None, count_match=None, aromatic_changes=None
    ):
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


def test_optional_typed_prefix_keeps_only_validated_acid_or_halide_hint(tmp_path):
    from rdkit import Chem

    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "amide",
                        "smarts": (
                            "[Cl,OH,O-:3][C$(C(=O)([CX4,c])),C$([CH](=O)):2]=[O:4]."
                            "[N$([NH2,NH3+1]([CX4,c])),N$([NH]([CX4,c])([CX4,c])):6]"
                            ">>[N+0:6]-[C:2]=[O:4]"
                        ),
                        "source": "test",
                    },
                    {
                        "id": "aryl",
                        "smarts": "[Cl,Br,I][c:1].[N;H2;$(NC):2]>>[c:1][N:2]",
                        "source": "test",
                    },
                ],
            }
        )
    )
    rules = Rules(path)
    for ident, hit, old, mixed in (
        ("amide", "CC(N)=O", ("CC=O", "N"), ("CC(=O)O", "N")),
        ("aryl", "Nc1ccccc1", ("c1ccccc1", "[NH2]"), ("Clc1ccccc1", "[NH2]")),
    ):
        searches = [
            RouteFirstSearch(
                hit,
                SimpleNamespace(sha256="test"),
                rules,
                RouteFirstLimits(**limits()),
                diagnose_continuation=True,
                preserve_typed_prefix=enabled,
            )
            for enabled in (False, True)
        ]
        assert searches[0]._reverse(rules.get(ident), hit, False) == [old]
        assert searches[1]._reverse(rules.get(ident), hit, False) == [mixed]
        for field in (
            "reverse_attempts",
            "reverse_role_trials",
            "reverse_role_work",
            "reverse_role_match_checks",
        ):
            assert searches[0].report[field] == searches[1].report[field]
        assert searches[1].report["mixed_typed_prefix_hints"] == 1
        assert searches[1].report["mixed_typed_prefix_parts"] == 1
        assert searches[1].report["reverse_hint_quality_counts"]["approximate"] == 1
        assert searches[1].report["reverse_hint_quality_counts"]["typed"] == 0
        assert searches[0].report["reverse_hint_normalizations"] == {"removed_originless_terminal_dummy": 1}
        assert searches[1].report["reverse_hint_normalizations"] == {
            "typed_originless_placeholder_instantiated": 1
        }
        registry = searches[1]._continuation_raw[(ident, hit, False)]
        assert tuple(registry["by_guide"]) == (mixed,)
        raw, clean = registry["by_guide"][mixed][0]
        assert len(raw) == len(clean) == 2
        assert tuple(Chem.MolToSmiles(part, isomericSmiles=True) for part in clean) == mixed
        assert searches[1]._reverse(rules.get(ident), hit, False) == [mixed]
        assert searches[1].report["reverse_attempts"] == 1

    fully_typed = RouteFirstSearch(
        "CC(=O)NC",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        preserve_typed_prefix=True,
    )
    assert fully_typed._reverse(rules.get("amide"), "CC(=O)NC", False) == [("CC(=O)O", "CN")]
    assert fully_typed.report["reverse_hint_quality_counts"]["typed"] == 1
    assert fully_typed.report["mixed_typed_prefix_hints"] == 0


def test_optional_typed_prefix_preserves_failure_and_budget_fallbacks(tmp_path, monkeypatch):
    import chembounce_synthesis.route_first as route_first
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "amide",
                        "smarts": (
                            "[Cl,OH,O-:3][C$(C(=O)([CX4,c])),C$([CH](=O)):2]=[O:4]."
                            "[N$([NH2,NH3+1]([CX4,c])),N$([NH]([CX4,c])([CX4,c])):6]"
                            ">>[N+0:6]-[C:2]=[O:4]"
                        ),
                        "source": "test",
                    },
                    {
                        "id": "reversed-amide",
                        "smarts": (
                            "[N$([NH2,NH3+1]([CX4,c])),N$([NH]([CX4,c])([CX4,c])):6]."
                            "[Cl,OH,O-:3][C$(C(=O)([CX4,c])),C$([CH](=O)):2]=[O:4]"
                            ">>[N+0:6]-[C:2]=[O:4]"
                        ),
                        "source": "test",
                    },
                ],
            }
        )
    )
    rules = Rules(path)
    rule = rules.get("amide")
    first_failed = RouteFirstSearch(
        "CC(N)=O",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        preserve_typed_prefix=True,
    )
    assert first_failed._reverse(rules.get("reversed-amide"), "CC(N)=O", False) == [("N", "CC=O")]
    assert first_failed.report["reverse_role_trials"] == 1
    assert first_failed.report["mixed_typed_prefix_hints"] == 0
    capped = RouteFirstSearch(
        "CC(N)=O",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        preserve_typed_prefix=True,
    )
    capped.report["reverse_role_work"] = route_first.MAX_REVERSE_ROLE_WORK_PER_SEARCH
    assert capped._reverse(rule, "CC(N)=O", False) == [("CC=O", "N")]
    assert capped.report["mixed_typed_prefix_hints"] == 0
    assert capped.report["reverse_role_trials"] == 0
    expired = RouteFirstSearch(
        "CC(N)=O",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        deadline_monotonic=time.monotonic() - 1,
        preserve_typed_prefix=True,
    )
    assert expired._reverse(rule, "CC(N)=O", False) == []
    assert expired.report["reverse_attempts"] == 0
    broken = RouteFirstSearch(
        "CC(N)=O",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        preserve_typed_prefix=True,
    )

    def invalid_approx(*args, **kwargs):
        raise route_first.ReverseHintError("invalid_remaining_precursor")

    monkeypatch.setattr(route_first, "normalize_reverse_hint", invalid_approx)
    assert broken._reverse(rule, "CC(N)=O", False) == []
    assert broken.report["mixed_typed_prefix_hints"] == 0
    assert broken.report["reverse_hint_rejections"]["invalid_remaining_precursor"] == 1


def test_typed_representative_size_proof_for_mixed_and_full_guides(tmp_path):
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "amide",
                        "smarts": (
                            "[Cl,OH,O-:3][C$(C(=O)([CX4,c])),C$([CH](=O)):2]=[O:4]."
                            "[N$([NH2,NH3+1]([CX4,c])),N$([NH]([CX4,c])([CX4,c])):6]"
                            ">>[N+0:6]-[C:2]=[O:4]"
                        ),
                        "source": "test",
                    }
                ],
            }
        )
    )
    rules = Rules(path)
    for hit, expected, quality in (
        ("CC(N)=O", ("CC(=O)O", "N"), "approximate"),
        ("CC(=O)NC", ("CC(=O)O", "CN"), "typed"),
    ):
        search = RouteFirstSearch(
            hit,
            SimpleNamespace(sha256="test"),
            rules,
            RouteFirstLimits(**limits()),
            preserve_typed_prefix=True,
            correct_typed_prefix_size=True,
        )
        assert search._reverse(rules.get("amide"), hit, False, record_centers=True) == [expected]
        assert search._guide_quality[("amide", hit, expected)] == quality
        assert search._guide_representative_heavy[("amide", hit, expected)] == 1
        assert search.report["representative_origin_known_guides"] == 1
        assert search.report["representative_origin_structural_units"] > 0
        assert search._reverse(rules.get("amide"), hit, False, record_centers=True) == [expected]
        assert search.report["reverse_attempts"] == 1
        assert search.report["representative_origin_checks"] == 1


def test_representative_size_proof_fails_closed_and_uses_inside_absolute_value(tmp_path):
    from rdkit import Chem

    import chembounce_synthesis.route_first as route_first
    from chembounce_synthesis.reactions import Rules

    path = tmp_path / "rules.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "amide",
                        "smarts": (
                            "[Cl,OH,O-:3][C$(C(=O)([CX4,c])),C$([CH](=O)):2]=[O:4]."
                            "[N$([NH2,NH3+1]([CX4,c])),N$([NH]([CX4,c])([CX4,c])):6]"
                            ">>[N+0:6]-[C:2]=[O:4]"
                        ),
                        "source": "test",
                    }
                ],
            }
        )
    )
    rules = Rules(path)
    rule = rules.get("amide")
    search = RouteFirstSearch(
        "CC(N)=O",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        preserve_typed_prefix=True,
        correct_typed_prefix_size=True,
        diagnose_continuation=True,
    )
    guide = search._reverse(rule, search.hit, False, record_centers=True)[0]
    raw, clean = search._continuation_raw[(rule.id, search.hit, False)]["by_guide"][guide][0]
    assert search._representative_origin_count(raw, clean, 1, search.hit_mol.GetNumAtoms()) == (
        1,
        None,
    )
    assert search._outer_size_delta(rule, guide, guide, 0) == 0
    search._guide_representative_heavy.clear()
    assert search._outer_size_delta(rule, guide, guide, 0) == 1
    search._guide_representative_heavy[(rule.id, search.hit, guide)] = 1
    assert search._outer_size_delta(rule, guide, guide, 1) == 1
    assert search.report["representative_origin_legacy_reasons"]["retrieval_target_changed"] == 1
    modified = tuple(Chem.Mol(part) for part in raw)
    parent = next(atom for atom in modified[0].GetAtoms() if atom.HasProp("react_atom_idx"))
    parent.SetAtomicNum(7)
    assert (
        search._representative_origin_count(modified, clean, 1, search.hit_mol.GetNumAtoms())[1]
        == "parent_atom_changed"
    )
    marked = tuple(Chem.Mol(part) for part in raw)
    dummy = next(atom for atom in marked[0].GetAtoms() if atom.GetAtomicNum() == 0)
    dummy.SetIntProp("react_atom_idx", 0)
    dummy.SetIntProp("react_idx", 0)
    assert search._representative_origin_count(marked, clean, 1, search.hit_mol.GetNumAtoms())[1] is not None

    counts, unknown = {}, {}
    search._merge_representative_origin(counts, unknown, guide, 1, None)
    search._merge_representative_origin(counts, unknown, guide, 0, None)
    assert guide not in counts
    assert unknown[guide] == "raw_alternatives_disagree"
    search._merge_representative_origin(counts, unknown, guide, 1, None)
    assert guide not in counts
    # Subtract inside abs: an already undersized normalized guide becomes
    # farther from the hit when its proven representative is discounted.
    assert search._outer_size_delta(rule, guide, ("C", "N"), 0) == 3

    capped = RouteFirstSearch(
        "CC(N)=O",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        preserve_typed_prefix=True,
        correct_typed_prefix_size=True,
    )
    capped.report["representative_origin_structural_units"] = route_first.MAX_REPRESENTATIVE_ORIGIN_VISITS
    assert capped._reverse(rule, capped.hit, False, record_centers=True) == [guide]
    assert capped._guide_representative_heavy == {}
    assert capped.report["representative_origin_legacy_reasons"]["structural_unit_limit"] == 1
    assert (
        capped.report["representative_origin_structural_units"]
        == route_first.MAX_REPRESENTATIVE_ORIGIN_VISITS
    )
    entry_capped = RouteFirstSearch(
        "CC(N)=O",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        preserve_typed_prefix=True,
        correct_typed_prefix_size=True,
    )
    entry_capped.report["representative_origin_guide_entries"] = route_first.MAX_REPRESENTATIVE_ORIGIN_GUIDES
    assert entry_capped._reverse(rule, entry_capped.hit, False, record_centers=True) == [guide]
    assert entry_capped._guide_representative_heavy == {}
    assert entry_capped.report["representative_origin_checks"] == 0
    assert entry_capped.report["representative_origin_legacy_reasons"]["guide_entry_limit"] == 1
    legacy = RouteFirstSearch(
        "CC(N)=O",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        preserve_typed_prefix=True,
    )
    assert legacy._reverse(rule, legacy.hit, False, record_centers=True) == [guide]
    assert legacy._outer_size_delta(rule, guide, guide, 0) == 1
    assert legacy.report["representative_origin_checks"] == 0
    expired = RouteFirstSearch(
        "CC(N)=O",
        SimpleNamespace(sha256="test"),
        rules,
        RouteFirstLimits(**limits()),
        deadline_monotonic=time.monotonic() - 1,
        preserve_typed_prefix=True,
        correct_typed_prefix_size=True,
    )
    assert expired._representative_origin_count(raw, clean, 1, search.hit_mol.GetNumAtoms()) == (
        None,
        "deadline",
    )
    assert expired.report["representative_origin_structural_units"] == 0


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


def test_reaction_pool_selected_roles_are_frequency_ranked_and_map_free():
    catalog = SimpleNamespace(sha256="synthetic")
    rules = SimpleNamespace(sha256="synthetic")
    search = RouteFirstSearch("Brc1ccccc1", catalog, rules, RouteFirstLimits(**limits()))
    aryl = SimpleNamespace(roles=("[Cl,Br,I][c:1]", "[N;$(NC):2]"))
    aliphatic = SimpleNamespace(roles=("[C@H:1](F)Br", "[N;$(NC):2]"))
    plans = [(rule, (), False, (), "cut", "typed") for rule in (aryl, aliphatic, aryl)]
    search._select_role_indexes(plans)
    assert [row["occurrences"] for row in search.report["reaction_pool_selected_roles"]] == [3, 2]
    assert search._selected_role_indexes == ("[N&$(NC)]", "[Cl,Br,I]c")
    assert all(":1" not in role and ":2" not in role for role in search._selected_role_indexes)
    chiral = route_first_module._normalized_role_smarts("[C@H:1](F)Br")
    assert "@" in chiral and ":1" not in chiral
    assert "$(NC)" in search._selected_role_indexes[0]


def test_reaction_pool_midpoint_and_same_site_reuse_are_context_bound():
    search = RouteFirstSearch(
        "Brc1ccccc1",
        SimpleNamespace(sha256="synthetic"),
        SimpleNamespace(sha256="synthetic"),
        RouteFirstLimits(**limits()),
    )
    roles = ("[Cl,Br,I]c", "[N;H2]c")
    search._selected_role_indexes = roles
    search._role_index_queries = {role: Chem.MolFromSmarts(role) for role in roles}
    aryl_features, aryl_index, complete = search._pool_metadata("Brc1ccccc1", "Brc1ccccc1", "[#6][Br]")
    assert complete and aryl_index == roles[0]
    assert aryl_features == route_first_module._pool_bins((7, 1, 6, 1))
    aliphatic_features, aliphatic_index, complete = search._pool_metadata("CCBr", "CCBr", "[#6][Br]")
    assert complete and aliphatic_index is None
    assert aliphatic_features == route_first_module._pool_bins((3, 1, 0, 0))
    assert search.report["reaction_pool_index_choices"]["same_site"] == 1
    assert search.report["reaction_pool_index_choices"]["legacy"] == 1
    assert search.report["reaction_pool_caller_only_sites"] >= 1

    # Retain some substituent size without changing the actual guide graph.
    larger = "Clc1ccc(CC(=O)NCC)cc1"
    midpoint, selected, complete = search._pool_metadata(larger, larger, "[Cl]c")
    assert complete and selected == roles[0]
    assert search._features(molecule(larger)) == (3, 1, 1, 1)
    assert midpoint == (2, 1, 1, 1)


def test_reaction_role_index_pool_is_complete_contextual_and_total_capped(monkeypatch):
    values = tuple(map(canonical, ("Brc1ccccc1", "Clc1ccccc1", "Nc1ccccc1", "CC")))
    catalog = SimpleNamespace(
        sha256="synthetic",
        mols={smiles: molecule(smiles) for smiles in values},
    )
    rules = SimpleNamespace(sha256="synthetic")
    settings = RouteFirstLimits(**limits(max_catalog_rows_per_query=2))
    aryl, amine = "[Cl,Br,I]c", "[N;H2]c"

    def make_search():
        search = RouteFirstSearch("Brc1ccccc1", catalog, rules, settings)
        search._selected_role_indexes = (aryl, amine)
        search._role_index_queries = {role: Chem.MolFromSmarts(role) for role in (aryl, amine)}
        return search

    search = make_search()
    aryl_pool = search._pool("Brc1ccccc1", guide="Brc1ccccc1", role=aryl)
    amine_pool = search._pool("Brc1ccccc1", guide="Brc1ccccc1", role=amine)
    assert search.report["reaction_pool_role_index_complete"]
    assert search.report["reaction_pool_role_index_checks"] == 2 * len(values)
    assert len(aryl_pool) == len(amine_pool) == len(set(aryl_pool)) == 2
    assert set(aryl_pool) == {canonical("Brc1ccccc1"), canonical("Clc1ccccc1")}
    assert amine_pool[0] == canonical("Nc1ccccc1")
    assert len(search._catalog_pools) == 2  # Same guide, distinct effective index context.
    assert search.report["reaction_pool_indexed_rows_selected"] <= 4

    monkeypatch.setattr(route_first_module, "MAX_ROLE_INDEX_MATCH_CHECKS", 1)
    capped = make_search()
    capped_pool = capped._pool("Brc1ccccc1", guide="Brc1ccccc1", role=aryl)
    assert capped.report["reaction_pool_role_index_checks"] == 1
    assert not capped.report["reaction_pool_role_index_complete"]
    assert capped.report["reaction_pool_role_index_failure"] == "match_check_cap"
    legacy = RouteFirstSearch("Brc1ccccc1", catalog, rules, settings)
    assert capped_pool == legacy._pool("Brc1ccccc1")
    assert capped.report["reaction_pool_index_disabled_incomplete"] == 1

    monkeypatch.setattr(route_first_module, "MAX_ROLE_INDEX_MATCH_CHECKS", 500000)
    failed = make_search()
    failed._role_index_queries[aryl] = None
    assert failed._pool("Brc1ccccc1", guide="Brc1ccccc1", role=aryl) == legacy._pool("Brc1ccccc1")
    assert not failed.report["reaction_pool_role_index_complete"]
    assert failed.report["reaction_pool_role_index_failure"] == "role_match_error"


def test_reaction_pool_unknown_match_budget_and_expired_deadline_use_legacy(monkeypatch):
    catalog = SimpleNamespace(
        sha256="synthetic",
        mols={"CCBr": molecule("CCBr")},
        by_smiles={"CCBr": ["block"]},
        blocks={"block": SimpleNamespace(id="block", smiles="CCBr")},
    )
    rules = SimpleNamespace(sha256="synthetic")
    search = RouteFirstSearch("CCBr", catalog, rules, RouteFirstLimits(**limits()))
    search._selected_role_indexes = ("[Cl,Br,I]c",)
    search._role_index_queries = {"[Cl,Br,I]c": Chem.MolFromSmarts("[Cl,Br,I]c")}
    monkeypatch.setattr(route_first_module, "MAX_POOL_SITE_MATCH_REQUESTS", 0)
    full, selected, complete = search._pool_metadata("CCBr", "CCBr", "[#6][Br]")
    assert full == search._features(molecule("CCBr")) and selected is None and not complete
    assert search.report["reaction_pool_site_match_requests"] == 0
    assert search.report["reaction_pool_metadata_unknowns"]["site_match_request_cap"] == 1
    assert not search._pool_metadata_cache
    assert [block.smiles for block in search._nearest("CCBr", "[#6][Br]")] == ["CCBr"]
    assert not search._nearest_cache

    larger = "Clc1ccc(CC(=O)NCC)cc1"
    different = RouteFirstSearch(larger, catalog, rules, RouteFirstLimits(**limits()))
    different._selected_role_indexes = ("[Cl,Br,I]c", "[N;H2]c")
    different._role_index_queries = {
        role: Chem.MolFromSmarts(role) for role in different._selected_role_indexes
    }
    features, index, complete = different._pool_metadata(larger, larger, "[Cl]c")
    assert features == (3, 1, 1, 1) and index is None and not complete
    assert not different._pool_metadata_cache

    monkeypatch.setattr(route_first_module, "MAX_POOL_GRAPH_VISITS", 0)
    graph_capped = RouteFirstSearch("CCBr", catalog, rules, RouteFirstLimits(**limits()))
    assert graph_capped._pool("CCBr", guide="CCBr", role="[#6][Br]") == ("CCBr",)
    assert graph_capped.report["reaction_pool_graph_visits"] == 0
    assert not graph_capped._catalog_pools

    expired = RouteFirstSearch(
        "CCBr", catalog, rules, RouteFirstLimits(**limits()), deadline_monotonic=time.monotonic() - 1
    )
    expired._selected_role_indexes = search._selected_role_indexes
    expired._role_index_queries = search._role_index_queries
    assert expired._pool("CCBr", guide="CCBr", role="[#6][Br]") == ()
    assert expired.report["reaction_pool_role_index_checks"] == 0


def test_reaction_pool_truncated_site_match_does_not_cache_partial_metadata():
    search = RouteFirstSearch(
        "CC",
        SimpleNamespace(sha256="synthetic"),
        SimpleNamespace(sha256="synthetic"),
        RouteFirstLimits(**limits()),
    )
    search._selected_role_indexes = ("[C]", "[N]")
    search._role_index_queries = {role: Chem.MolFromSmarts(role) for role in search._selected_role_indexes}
    long_chain = "C" * 70
    features, index, complete = search._pool_metadata(long_chain, long_chain, "[#6]")
    assert features == search._features(molecule(long_chain))
    assert index is None and not complete
    assert search.report["reaction_pool_metadata_unknowns"]["window_match_truncated"] == 1
    assert search.report["reaction_pool_site_match_requests"] == 0
    assert not search._pool_metadata_cache


def test_reaction_pool_oversized_guide_skips_added_raw_count_work(monkeypatch):
    search = RouteFirstSearch(
        "CC",
        SimpleNamespace(sha256="synthetic"),
        SimpleNamespace(sha256="synthetic"),
        RouteFirstLimits(**limits()),
    )

    def unexpected_raw_counts(*_args, **_kwargs):
        raise AssertionError("Oversized guide must use legacy full-guide features")

    monkeypatch.setattr(route_first_module, "_pool_raw_counts", unexpected_raw_counts)
    features, index, complete = search._pool_metadata("CC", "C" * 257, "[C]")
    assert features == search._features(molecule("CC")) and index is None and not complete
    assert search.report["reaction_pool_metadata_unknowns"]["query_size_or_role_parse"] == 1
    assert search.report["reaction_pool_graph_visits"] == 0


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
    search._pool = lambda _target, **_kwargs: ("CC(=O)O",)
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
    search._pool = lambda _target, **_kwargs: (ordinary,)
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
    assert report["catalog_stereo_exclusions"] > 0
    assert any(item["building_block_id"] == "alanine" for item in report["catalog_stereo_exclusion_samples"])


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
    assert manifest["route_first"]["source_sha256"]["route_first_nh.py"]
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
