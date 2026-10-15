"""Chemical-family equivalence must not mutate the structures being evaluated."""

import pytest
from rdkit import Chem

from chembounce_synthesis import identity
from chembounce_synthesis.common import canonical, molecule
from chembounce_synthesis.identity import (
    chemical_identity,
    classify_change,
    family_key,
    scaffold_family_key,
)


@pytest.mark.parametrize(
    "original,changed",
    [
        ("O=c1cccc[nH]1", "Oc1ccccn1"),
        ("CC(=O)Nc1ccncc1", "CC(=O)Nc1cc[nH+]cc1"),
        ("CC(=O)O", "CC(=O)[O-]"),
        ("NCC(=O)O", "[NH3+]CC(=O)[O-]"),
        ("C[N+](C)(C)CC(=O)O", "C[N+](C)(C)CC(=O)[O-]"),
        ("CC[CH+]CN", "CC[CH+]C[NH3+]"),
        ("CC(=O)Nc1ccccc1", "CC(=O)Nc1cc[13cH]cc1"),
        ("C1CC[C@H]2CCCC[C@@H]2C1", "C1CC[C@@H]2CCCC[C@@H]2C1"),
        ("OC", "[2H]OC"),
    ],
)
def test_state_changes_share_family_but_retain_exact_state(original, changed):
    before, after = chemical_identity(original), chemical_identity(changed)
    assert before.family_smiles == after.family_smiles
    assert before.scaffold_family_smiles == after.scaffold_family_smiles
    assert classify_change(original, changed) == "state_only"
    assert before.exact_smiles == canonical(original)
    assert after.exact_smiles == canonical(changed)
    assert before.exact_smiles != after.exact_smiles


@pytest.mark.parametrize(
    "charged,other",
    [("CC[CH2+]", "CCC"), ("[BH4-]", "B"), ("C[N+](C)(C)C", "CN(C)C")],
)
def test_family_identity_does_not_remove_irreversible_charges(charged, other):
    parent = molecule(family_key(charged))
    assert Chem.GetFormalCharge(parent) == Chem.GetFormalCharge(molecule(charged))
    assert family_key(charged) != family_key(other)


@pytest.mark.parametrize(
    "original,changed,classification",
    [
        ("Oc1ccccc1", "c1ccccc1O", "unchanged"),
        ("Cc1ccccc1", "CCc1ccccc1", "decoration"),
        ("CC", "CCC", "decoration"),
        ("c1ccccc1", "c1ccncc1", "scaffold"),
        ("C1CCCCC1", "C1CC1", "scaffold"),
        ("CC", "C1CC1", "scaffold"),
        ("C1CC1", "CC", "scaffold"),
    ],
)
def test_change_classes_separate_state_decoration_and_scaffold(original, changed, classification):
    assert classify_change(original, changed) == classification


def test_atom_aware_scaffold_differs_from_ring_linker_topology():
    benzene, pyridine = chemical_identity("c1ccccc1"), chemical_identity("c1ccncc1")
    assert benzene.scaffold_family_smiles != pyridine.scaffold_family_smiles
    assert benzene.ring_skeleton_smiles == pyridine.ring_skeleton_smiles
    pyridone = chemical_identity("O=c1cccc[nH]1")
    assert pyridone.ring_skeleton_smiles == pyridine.ring_skeleton_smiles
    assert not chemical_identity("CCC").has_scaffold
    assert scaffold_family_key("CCC") == ""


def test_identity_is_canonical_and_policy_versioned():
    mol = molecule("CC(=O)Nc1cc[nH+]cc1")
    expected = chemical_identity(Chem.MolToSmiles(mol))
    for _ in range(12):
        assert chemical_identity(Chem.MolToSmiles(mol, doRandom=True)) == expected
    assert chemical_identity("[CH3:1][OH:2]") == chemical_identity("CO")
    assert expected.policy_version == identity.IDENTITY_POLICY
    assert expected.rdkit_version


def test_disconnected_salts_are_rejected_instead_of_silently_stripped():
    with pytest.raises(ValueError, match="connected"):
        family_key("C[NH3+].[Cl-]")


def test_incomplete_tautomer_normalization_fails_closed(monkeypatch):
    identity._identity.cache_clear()
    with monkeypatch.context() as patch:
        patch.setattr(identity, "TAUTOMER_LIMIT", 1)
        with pytest.raises(ValueError, match="enumeration incomplete"):
            family_key("CC(=O)CC(=O)C")
    identity._identity.cache_clear()
    assert family_key("CC(=O)CC(=O)C")


@pytest.mark.parametrize(
    "original,changed",
    [("O=c1cccc[nH]1", "Oc1ccccn1"), ("CC(=O)Nc1ccncc1", "CC(=O)Nc1cc[nH+]cc1")],
)
def test_hopping_rejects_state_only_false_positive(original, changed):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from chembounce_synthesis.workflow import check_constraints

    with pytest.raises(ValueError, match="only chemical state"):
        check_constraints(changed, original, changed)


def test_exact_core_and_design_preservation_are_not_standardized():
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from chembounce_synthesis.workflow import check_constraints

    neutral = "CC(=O)Nc1ccncc1"
    protonated = "CC(=O)Nc1cc[nH+]cc1"
    hit = "CC(=O)Nc1ccccc1"
    check_constraints(protonated, hit, protonated)
    with pytest.raises(ValueError, match="preserve the selected design"):
        check_constraints(protonated, hit, neutral, preserve_design_scaffold=True)
    with pytest.raises(ValueError, match="anchor"):
        check_constraints(neutral, hit, neutral, core_smiles="[nH+]1ccccc1")


def test_expansion_does_not_count_sidechain_protonation_as_a_new_analogue():
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from chembounce_synthesis.workflow import check_constraints

    with pytest.raises(ValueError, match="unchanged design chemical family"):
        check_constraints("O=C([O-])c1ccccc1", "c1ccncc1", "O=C(O)c1ccccc1", mode="expand")


def test_failed_normalization_is_cached_but_never_becomes_an_identity(monkeypatch):
    identity.clear_identity_caches()
    original = identity._canonical_tautomer
    calls = []

    def counting(mol):
        calls.append(Chem.MolToSmiles(mol))
        return original(mol)

    monkeypatch.setattr(identity, "_canonical_tautomer", counting)
    with monkeypatch.context() as patch:
        patch.setattr(identity, "TAUTOMER_LIMIT", 1)
        for _ in range(2):
            with pytest.raises(ValueError, match="enumeration incomplete"):
                family_key("CC(=O)CC(=O)C")
        assert len(calls) == 1
        assert identity.identity_cache_info()["identities_including_failures"]["hits"] == 1
    # Changing an enumeration bound changes the cache key: the cached failure
    # cannot contaminate a fully completed normalization under the normal bound.
    assert family_key("CC(=O)CC(=O)C")
    assert len(calls) == 2
    identity.clear_identity_caches()


def test_shared_scaffold_normalization_is_reused_without_changing_families():
    identity.clear_identity_caches()
    first, second = chemical_identity("Cc1ccncc1"), chemical_identity("CCc1ccncc1")
    assert first.family_smiles != second.family_smiles
    assert first.scaffold_family_smiles == second.scaffold_family_smiles
    assert identity.identity_cache_info()["scaffolds"]["hits"] == 1
    assert chemical_identity("Cc1ccncc1") == first
    assert identity.identity_cache_info()["exact_inputs"]["hits"] == 1
    identity.clear_identity_caches()
