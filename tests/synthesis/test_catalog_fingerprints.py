"""Lazy lookup must preserve the validated catalog and eager search results."""

import gc
import weakref
from collections.abc import MutableMapping
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
from rdkit import Chem, DataStructs

import chembounce_synthesis.catalog as catalog_module
from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.common import FP, canonical, molecule
from chembounce_synthesis.routes import SearchLimits, expand_route, resolve_proposal, verify_route


def fingerprint_spy(monkeypatch):
    spy = Mock(wraps=FP.GetFingerprint)
    monkeypatch.setattr(catalog_module, "FP", SimpleNamespace(GetFingerprint=spy))
    return spy


def write_catalog(tmp_path, rows):
    path = tmp_path / "catalog.tsv"
    path.write_text(
        "id\tsmiles\tsource\tsnapshot\tavailability\n" + "".join("\t".join(row) + "\n" for row in rows)
    )
    return path


def test_construction_role_matching_and_route_replay_need_no_catalog_fingerprints(assets, route, monkeypatch):
    original, rules = assets
    expected = verify_route(route, original, rules)
    spy = fingerprint_spy(monkeypatch)
    catalog = Catalog(original.path)

    assert catalog.get("acetic") == original.get("acetic")
    assert catalog.evidence("acetic") == original.evidence("acetic")
    assert catalog.compatible(rules.get("amide").roles[0]) == original.compatible(rules.get("amide").roles[0])
    assert verify_route(route, catalog, rules) == expected
    spy.assert_not_called()


def test_keys_membership_and_unknown_keys_do_not_calculate_fingerprints(assets, monkeypatch):
    spy = fingerprint_spy(monkeypatch)
    catalog = Catalog(assets[0].path)
    fps = catalog.fps
    assert isinstance(fps, MutableMapping)
    assert len(fps) == len(catalog.mols)
    assert list(fps) == list(fps.keys()) == list(catalog.mols)
    assert list(reversed(fps)) == list(reversed(catalog.mols))
    assert all(smiles in fps for smiles in catalog.mols)
    assert "missing" not in fps
    assert fps.get("missing") is None
    marker = object()
    assert fps.get("missing", marker) is marker
    assert fps.pop("missing", marker) is marker
    with pytest.raises(KeyError):
        fps["missing"]
    with pytest.raises(KeyError):
        del fps["missing"]
    with pytest.raises(KeyError):
        fps.pop("missing")
    spy.assert_not_called()


def test_lookup_uses_canonical_molecule_and_captured_generator_once(assets, monkeypatch):
    spy = fingerprint_spy(monkeypatch)
    catalog = Catalog(assets[0].path)
    replacement = Mock(side_effect=AssertionError("Generator policy changed after construction"))
    monkeypatch.setattr(catalog_module, "FP", SimpleNamespace(GetFingerprint=replacement))
    key = catalog.get("alanine_s").smiles
    value = catalog.fps[key]
    assert value.ToBitString() == FP.GetFingerprint(molecule(key)).ToBitString()
    assert catalog.fps[key] is value
    assert catalog.fps.get(key) is value
    assert catalog.fps.setdefault(key, None) is value
    assert spy.call_count == 1
    assert spy.call_args.args[0] is catalog.mols[key]
    replacement.assert_not_called()


@pytest.mark.parametrize("view", ["items", "values", "copy"])
def test_value_views_and_copy_return_real_fingerprints(assets, monkeypatch, view):
    spy = fingerprint_spy(monkeypatch)
    catalog = Catalog(assets[0].path)
    fps = catalog.fps
    if view == "copy":
        result = fps.copy()
        assert isinstance(result, dict)
        values = list(result.values())
        assert list(result) == list(fps)
        result.clear()
        assert len(fps) == len(catalog.mols)
    else:
        values_view = getattr(fps, view)()
        spy.assert_not_called()
        iterator = iter(values_view)
        first = next(iterator)
        assert spy.call_count == 1
        result = [first, *iterator]
        values = [value for _, value in result] if view == "items" else result
        if view == "items":
            assert [key for key, _ in result] == list(fps)
    assert [value.ToBitString() for value in values] == [
        FP.GetFingerprint(mol).ToBitString() for mol in catalog.mols.values()
    ]
    assert spy.call_count == len(fps)
    assert all(fps[key] is value for key, value in zip(fps, values))
    assert spy.call_count == len(fps)


def test_mapping_mutations_preserve_order_and_do_not_force_discarded_values(assets, monkeypatch):
    spy = fingerprint_spy(monkeypatch)
    catalog = Catalog(assets[0].path)
    fps = catalog.fps
    keys = list(fps)
    override = FP.GetFingerprint(molecule("C"))
    fps[keys[0]] = override
    assert fps[keys[0]] is override
    del fps[keys[1]]
    assert keys[1] not in fps and keys[1] in catalog.mols
    fps.update({"added": override})
    assert list(fps) == [key for key in keys if key != keys[1]] + ["added"]
    assert fps.popitem() == ("added", override)
    spy.assert_not_called()

    # Existing pending values are returned, not replaced by setdefault.
    assert (
        fps.setdefault(keys[2], override).ToBitString()
        == FP.GetFingerprint(catalog.mols[keys[2]]).ToBitString()
    )
    assert spy.call_count == 1
    last_key, last_value = fps.popitem()
    assert last_key == keys[-1]
    assert last_value.ToBitString() == FP.GetFingerprint(catalog.mols[last_key]).ToBitString()
    assert spy.call_count == 2
    assert fps.setdefault("new", override) is override
    assert fps.pop("new") is override
    fps.clear()
    assert not fps and len(catalog.mols) == len(keys)
    assert spy.call_count == 2
    with pytest.raises(KeyError):
        fps.popitem()


@pytest.mark.parametrize(
    "bad_row,error",
    [
        (("bad", "notSMILES", "test", "date", "catalogued"), "Invalid SMILES"),
        (("bad", "notSMILES", "test", "date", "unavailable"), "Invalid SMILES"),
        (("bad", "CCO.Cl", "test", "date", "unavailable"), "connected molecule"),
        (("bad", "CC*", "test", "date", "catalogued"), "without dummy atoms"),
        (("bad", "CCO", "", "date", "catalogued"), "Empty catalog provenance"),
        (("valid", "CCN", "test", "date", "catalogued"), "Duplicate building-block ID"),
        (("bad", "CCO", "test", "date", "unknown"), "Unknown availability"),
    ],
)
def test_all_rows_still_validate_before_any_fingerprint(tmp_path, monkeypatch, bad_row, error):
    spy = fingerprint_spy(monkeypatch)
    path = write_catalog(tmp_path, [("valid", "CCO", "test", "date", "catalogued"), bad_row])
    with pytest.raises(ValueError, match=error):
        Catalog(path)
    spy.assert_not_called()


def test_unavailable_catalog_is_rejected_without_fingerprints(tmp_path, monkeypatch):
    spy = fingerprint_spy(monkeypatch)
    path = write_catalog(tmp_path, [("absent", "CCO", "test", "date", "unavailable")])
    with pytest.raises(ValueError, match="no eligible building blocks"):
        Catalog(path)
    spy.assert_not_called()


def test_duplicate_skus_availability_and_canonical_atom_order_are_preserved(tmp_path, monkeypatch):
    spy = fingerprint_spy(monkeypatch)
    path = write_catalog(
        tmp_path,
        [
            ("z-sku", "OCC", "test", "date", "catalogued"),
            ("a-sku", "CCO", "test", "date", "in_stock"),
            ("unavailable-sku", "CCO", "test", "date", "unavailable"),
            ("absent-only", "CCCl", "test", "date", "unavailable"),
            ("stereo", "N[C@@H](C)C(=O)O", "test", "date", "catalogued"),
            ("isotope-charge", "[13CH3]C[NH3+]", "test", "date", "catalogued"),
        ],
    )
    catalog = Catalog(path)
    assert len(catalog.blocks) == 6 and len(catalog.fps) == 3
    assert catalog.by_smiles["CCO"] == ["z-sku", "a-sku"]
    assert "CCCl" not in catalog.fps
    with pytest.raises(ValueError, match="unavailable"):
        catalog.get("unavailable-sku")
    spy.assert_not_called()
    for smiles, mol in catalog.mols.items():
        expected_mol = molecule(smiles)
        assert Chem.MolToSmiles(mol, canonical=False) == Chem.MolToSmiles(expected_mol, canonical=False)
        assert catalog.fps[smiles].ToBitString() == FP.GetFingerprint(expected_mol).ToBitString()
    assert spy.call_count == 3
    assert catalog.nearest("OCC", "[O;H1]", k=1)[0].id == "a-sku"


@pytest.mark.parametrize(
    "query,role,k,minimum,required",
    [
        ("CC(=O)O", "[C](=O)[O;H1]", 8, 0.0, ()),
        ("Nc1ccncc1", "[N;H1,H2]", 4, 0.0, ()),
        ("CCN", "[N;H1,H2]", 5, 0.2, ("[C](=O)[O;H1]",)),
        ("N[C@@H](C)C(=O)O", "[N;H1,H2]", 8, 0.0, ()),
        ("CCO", "[Cl]", 3, 0.0, ()),
    ],
)
def test_nearest_matches_eager_fingerprint_dictionary(assets, query, role, k, minimum, required):
    lazy = Catalog(assets[0].path)
    eager = Catalog(assets[0].path)
    eager.fps = {smiles: FP.GetFingerprint(mol) for smiles, mol in eager.mols.items()}
    assert lazy.nearest(query, role, k, minimum, required) == eager.nearest(query, role, k, minimum, required)


def test_exact_match_and_similarity_diversity_ties_keep_deterministic_order(tmp_path, monkeypatch):
    # Force equal fingerprints to exercise both nearest and diversity tie rules.
    def constant_fingerprint(mol):
        result = DataStructs.ExplicitBitVect(2048)
        result.SetBit(1)
        return result

    monkeypatch.setattr(catalog_module, "FP", SimpleNamespace(GetFingerprint=constant_fingerprint))
    rows = [
        (f"sku-{smiles}", smiles, "test", "date", "catalogued")
        for smiles in ("CO", "CCO", "C", "CC", "CCC", "CCCC")
    ]
    path = write_catalog(tmp_path, rows)
    lazy, eager = Catalog(path), Catalog(path)
    eager.fps = {smiles: constant_fingerprint(mol) for smiles, mol in eager.mols.items()}
    expected = ["CCC", "C", "CC", "CO"]
    assert [b.smiles for b in lazy.nearest("CCC", "[#6]", k=4)] == expected
    assert lazy.nearest("CCC", "[#6]", k=4) == eager.nearest("CCC", "[#6]", k=4)
    lazy.nearest.cache_clear()
    assert [b.smiles for b in lazy.nearest("CCC", "[#6]", k=4)] == expected


def test_route_resolution_expansion_and_evidence_match_eager_catalog(assets, route, proposal):
    lazy, rules = assets
    eager = Catalog(lazy.path)
    eager.fps = {smiles: FP.GetFingerprint(mol) for smiles, mol in eager.mols.items()}
    limits = SearchLimits(building_block_top_k=4, beam_width=8)
    assert resolve_proposal(proposal, lazy, rules, limits) == resolve_proposal(proposal, eager, rules, limits)
    assert expand_route(route, lazy, rules, limits, rounds=2, max_routes=10) == expand_route(
        route, eager, rules, limits, rounds=2, max_routes=10
    )
    assert verify_route(route, lazy, rules) == verify_route(route, eager, rules)


def test_catalog_instances_do_not_share_fingerprint_values_and_are_released(assets, monkeypatch):
    spy = fingerprint_spy(monkeypatch)
    first, second = Catalog(assets[0].path), Catalog(assets[0].path)
    key = canonical("CC(=O)O")
    assert first.fps[key] is not second.fps[key]
    assert first.fps[key].ToBitString() == second.fps[key].ToBitString()
    assert spy.call_count == 2
    first.fps[key] = None
    assert second.fps[key] is not None
    # Also exercise the existing per-instance nearest/role cache cycle.
    second.nearest("CC(=O)O", "[C](=O)[O;H1]", k=2)
    refs = [weakref.ref(first), weakref.ref(first.fps), weakref.ref(second), weakref.ref(second.fps)]
    del first, second
    gc.collect()
    assert all(ref() is None for ref in refs)
