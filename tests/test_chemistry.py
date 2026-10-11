import pytest
from rdkit import Chem
import chemistry
import utils
import scaffoldgraph_fragmenter as sg
from conftest import ADP


def mol(s):
    return Chem.MolFromSmiles(s)


def test_adp_all_scaffolds_no_longer_crash():
    target = mol(ADP)
    fragments = sg.get_all_murcko_fragments(target, False, 1)
    assert len(fragments) == 3
    for fragment in fragments:
        results = chemistry.replace(target, fragment, mol("c1ccncc1"), 5, 40)
        assert results
        for score, smiles in results:
            product = mol(smiles)
            assert product is not None
            assert len(Chem.GetMolFrags(product)) == 1
            assert all(a.GetAtomicNum() > 0 for a in product.GetAtoms())
            assert 0 <= score <= 1


def test_retained_adp_ribose_stereo_and_phosphate():
    target = mol(ADP)
    pattern = mol("c1ncc2nc[nH]c2n1")
    retained_core = mol("O=P(O)(O)OP(=O)(O)OC[C@H]1O[C@H](N)[C@H](O)[C@@H]1O")
    # Use a retained ribose/phosphate fragment without defining its new substituent.
    retained_core = mol("O=P(O)(O)OP(=O)(O)OC[C@H]1OC[C@H](O)[C@@H]1O")
    assert target.HasSubstructMatch(retained_core, useChirality=True)
    for _, smiles in chemistry.replace(target, pattern, mol("c1ccncc1"), 8, 40):
        assert mol(smiles).HasSubstructMatch(retained_core, useChirality=True)


@pytest.mark.parametrize(
    "smiles,pattern",
    [
        ("C[C@H](O)c1ccccc1", "c1ccccc1"),
        ("N[C@@H](CO)c1ccccc1", "c1ccccc1"),
        ("COc1cc(Cl)ccc1N", "c1ccccc1"),
        ("Cc1c(C)cccc1", "c1ccccc1"),
        ("C1CC2CCC1C2", "C1CCCCC1"),
    ],
)
def test_identity_replacement_can_reconstruct_input(smiles, pattern):
    target, query = mol(smiles), mol(pattern)
    result = chemistry.replace(target, query, query, 500, 500)
    assert Chem.MolToSmiles(target) in {s for _, s in result}


def test_deterministic_and_bounded():
    target, query, replacement = mol(ADP), mol("C1CCOC1"), mol("C1CCNC1")
    first = chemistry.replace(target, query, replacement, 3, 10)
    assert first == chemistry.replace(target, query, replacement, 3, 10)
    assert 0 < len(first) <= 3


@pytest.mark.parametrize("invalid", ["", "garbage", "*CC", "C1CCC", "[R]"])
def test_invalid_smiles(invalid):
    with pytest.raises(utils.SmilesInputError):
        utils.init_mol(invalid)


def test_acyclic_fragmentation_and_depth():
    assert sg.get_all_murcko_fragments(mol("CCO")) == []
    target = mol("c1ccccc1CC1CCCCC1")
    assert len(sg.get_all_murcko_fragments(target, False, 0)) == 1
    assert len(sg.get_all_murcko_fragments(target, False, 1)) > 1


@pytest.mark.parametrize(
    "smiles,pattern",
    [
        ("C=C1CCCC1", "C1CCCC1"),
        ("C/C=C/c1ccccc1", "c1ccccc1"),
        ("[13CH3][C@H](O)c1ccccc1", "c1ccccc1"),
    ],
)
def test_bond_order_double_bond_stereo_and_isotopes_survive(smiles, pattern):
    target, query = mol(smiles), mol(pattern)
    candidates = chemistry.replace(target, query, query, 100, 100)
    assert Chem.MolToSmiles(target) in {s for _, s in candidates}


def test_partial_fused_aromatic_cut_has_actionable_error():
    with pytest.raises(ValueError, match="entire aromatic ring system"):
        chemistry.replace(mol("Cc1cccc2ccccc12"), mol("c1ccccc1"), mol("c1ccncc1"))


def test_impossible_attachment_capacity_stops_before_enumeration(monkeypatch):
    # Twelve retained methyl groups cannot attach to six aromatic carbon atoms.
    target = mol("CC1(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C1(C)C")

    def unexpected(*args):
        pytest.fail("Impossible replacement entered expensive attachment enumeration")

    monkeypatch.setattr(chemistry, "_attach", unexpected)
    assert chemistry.replace(target, mol("C1CCCCC1"), mol("c1ccccc1")) == []


@pytest.mark.parametrize(
    "smiles,pattern,replacement",
    [
        ("Cc1c(C)c(O)c(N)c(C)c1F", "c1ccccc1", "C1CCNCC1"),
        ("CC1(C)CCCCC1", "C1CCCCC1", "C1CCCC1"),
        ("Cn1cncc1C", "c1c[nH]cn1", "c1cc[nH]c1"),
        ("COc1ccccc1", "c1ccccc1", "[NH2+]1CCCCC1"),
        ("COc1ccccc1", "c1ccccc1", "C1CSCC1"),
    ],
)
def test_capacity_pruning_preserves_feasible_products(monkeypatch, smiles, pattern, replacement):
    args = (mol(smiles), mol(pattern), mol(replacement), 20, 50)
    bounded = chemistry.replace(*args)
    assert bounded
    monkeypatch.setattr(chemistry, "_insufficient_capacity", lambda *args: False)
    assert chemistry.replace(*args) == bounded
