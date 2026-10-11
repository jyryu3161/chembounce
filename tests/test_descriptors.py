import numpy as np
import pytest
from rdkit import Chem
from cats_module import CATS2D
import utils


def test_cats_all_smarts_and_distinct_pairs():
    cats = CATS2D()
    labels = cats.getPcoreGroups(Chem.MolFromSmiles("[NH3+]CCN"))
    assert "P" in labels[0] and "P" in labels[-1]
    assert len(labels[-1]) == len(set(labels[-1]))
    assert len(cats.getCATs2D(Chem.MolFromSmiles("NCCO"))) == 210


def test_cats_occupancy_scales_every_row():
    mol = Chem.MolFromSmiles("NCCO")
    raw = np.array(CATS2D().getCATs2D(mol)).reshape(21, 10)
    cats = CATS2D(scale_type="occ")
    scaled = np.array(cats.getCATs2D(mol)).reshape(21, 10)
    labels = cats.getPcoreGroups(mol)
    for i, pair in enumerate(cats.cats_desc):
        divisor = sum(label.count(pair[0]) + label.count(pair[1]) for label in labels)
        if divisor:
            np.testing.assert_allclose(scaled[i], raw[i] / divisor)


def test_properties_accept_numpy_numbers_and_reject_nonfinite():
    assert utils.determ_props({"mw": np.float64(50)}, {"mw": (0, 100)})[0]
    for value in (np.nan, np.inf, -np.inf, None):
        assert not utils.determ_props({"mw": value}, {"mw": (0, 100)})[0]


def test_shape_is_repeatable_without_cache():
    first = utils.shape_descriptor("CCOc1ccccc1", 42).copy()
    utils.shape_descriptor.cache_clear()
    second = utils.shape_descriptor("CCOc1ccccc1", 42)
    np.testing.assert_allclose(first, second, atol=1e-12)
    assert utils.calc_electron_shape("CCOc1ccccc1", "CCOc1ccccc1") == pytest.approx(1)


def test_electroshape_agrees_with_oddt_reference():
    pytest.importorskip("oddt")
    from rdkit.Chem import AllChem
    from oddt.toolkits.rdk import Molecule
    from oddt.shape import electroshape as reference
    from shape_descriptors import electroshape

    mol = Chem.AddHs(Chem.MolFromSmiles("CCOc1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)
    np.testing.assert_allclose(electroshape(mol), reference(Molecule(mol)), rtol=2e-5, atol=2e-5)


@pytest.mark.parametrize("smiles", ["NCCO", "[NH3+]CC(=O)[O-]", "Clc1ccc(O)c(N)c1", "CC(C)(C)C1CCCCC1"])
def test_vectorized_cats_matches_ordered_pair_definition(smiles):
    mol = Chem.MolFromSmiles(smiles)
    cats = CATS2D()
    labels = cats.getPcoreGroups(mol)
    distances = Chem.GetDistanceMatrix(mol)
    reference = np.zeros((21, 10))
    for row, pair in enumerate(cats.cats_desc):
        a, b = pair
        for i in range(mol.GetNumAtoms()):
            for j in range(mol.GetNumAtoms()):
                distance = int(distances[i, j])
                if distance <= 9 and (
                    (a in labels[i] and b in labels[j]) or (b in labels[i] and a in labels[j])
                ):
                    reference[row, distance] += 1
    np.testing.assert_array_equal(cats.getCATs2D(mol), reference.ravel())
    reverse = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    np.testing.assert_array_equal(cats.getCATs2D(mol), cats.getCATs2D(reverse))


def test_graph_helpers_exclude_origin_and_preserve_single_atom_isotope():
    molecule = Chem.MolFromSmiles("[13CH3]CCO")
    assert utils.get_neighbor_atoms(molecule, 0, 0) == []
    assert utils.get_neighbor_atoms(molecule, 0, 2) == [1, 2]
    assert utils.get_submol_from_atom_list(molecule, [0]).GetAtomWithIdx(0).GetIsotope() == 13
