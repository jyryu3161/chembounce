"""Graph traversal must retain strict validation, atom order and map normalization."""

from types import SimpleNamespace
from unittest.mock import Mock

import pytest
from rdkit import Chem

from chembounce_synthesis import common


@pytest.mark.parametrize(
    "smiles",
    [
        "C",
        "OCC",
        "c1cc[nH]c1",
        "[13CH3][C@H](F)Cl",
        "[CH3:4][C@@H:2]([OH:3])[NH3+:1]",
        "N[C@@H]([2H])C(=O)O",
        "C/C=C/C",
        "C/C=C\\C",
        "C[C@H](F)[C@@H](O)Cl |a:1,3|",
        "[H][C@](F)(Cl)Br |a:1|",
    ],
)
def test_normalization_retains_atom_order_properties_and_stereo(smiles):
    expected = common.parse_absolute_smiles(smiles)
    for atom in expected.GetAtoms():
        atom.SetAtomMapNum(0)
    actual = common.molecule(smiles)
    assert actual.ToBinary(Chem.PropertyPickleOptions.AllProps) == expected.ToBinary(
        Chem.PropertyPickleOptions.AllProps
    )
    assert Chem.MolToSmiles(actual, canonical=False) == Chem.MolToSmiles(expected, canonical=False)
    assert common.FP.GetFingerprint(actual).ToBitString() == common.FP.GetFingerprint(expected).ToBitString()
    assert Chem.MolToCXSmiles(actual) == Chem.MolToCXSmiles(expected)


@pytest.mark.parametrize("smiles", ["[CH3:9]*", "*[CH3:9]", "[CH3:9]Cl.[Na+]", "C.*"])
def test_rejected_graph_is_not_partially_map_normalized(smiles, monkeypatch):
    mol = common.parse_absolute_smiles(smiles)
    before = mol.ToBinary(Chem.PropertyPickleOptions.AllProps)
    monkeypatch.setattr(common, "parse_absolute_smiles", lambda _: mol)
    with pytest.raises(ValueError, match="connected molecule without dummy atoms"):
        common.molecule(smiles)
    assert mol.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


def test_connectivity_failure_precedes_atom_access(monkeypatch):
    mol = SimpleNamespace(GetNumAtoms=Mock(side_effect=AssertionError("atoms accessed")))
    monkeypatch.setattr(common, "parse_absolute_smiles", lambda _: mol)
    monkeypatch.setattr(Chem, "GetMolFrags", lambda _: ((0,), (1,)))
    with pytest.raises(ValueError, match="connected molecule"):
        common.molecule("disconnected")
    mol.GetNumAtoms.assert_not_called()


@pytest.mark.parametrize("dummy_index", [None, 0, 2])
def test_atom_access_once_and_all_validation_before_map_writes(monkeypatch, dummy_index):
    events = []

    def atom(index):
        return SimpleNamespace(
            GetAtomicNum=lambda: events.append(("check", index)) or (0 if index == dummy_index else 6),
            SetAtomMapNum=lambda value: events.append(("map", index, value)),
        )

    atoms = [atom(index) for index in range(3)]
    mol = SimpleNamespace(
        GetNumAtoms=lambda: len(atoms),
        GetAtomWithIdx=lambda index: events.append(("get", index)) or atoms[index],
    )
    monkeypatch.setattr(common, "parse_absolute_smiles", lambda _: mol)
    monkeypatch.setattr(Chem, "GetMolFrags", lambda _: ((0, 1, 2),))
    if dummy_index is None:
        assert common.molecule("fixture") is mol
        count = 3
    else:
        with pytest.raises(ValueError, match="without dummy atoms"):
            common.molecule("fixture")
        count = dummy_index + 1
    expected = [event for index in range(count) for event in (("get", index), ("check", index))]
    if dummy_index is None:
        expected += [("map", index, 0) for index in range(3)]
    assert events == expected


@pytest.mark.parametrize("property_name,value", [("_MolFileBondStereo", 1), ("_MolFileBondCfg", 3)])
def test_unresolved_wedge_properties_are_checked_without_mutation(property_name, value):
    mol = Chem.MolFromSmiles("CC")
    mol.GetBondWithIdx(0).SetIntProp(property_name, value)
    before = mol.ToBinary(Chem.PropertyPickleOptions.AllProps)
    with pytest.raises(ValueError, match="unresolved.*wedging"):
        common.require_absolute_stereo(mol)
    assert mol.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


@pytest.mark.parametrize("direction", [Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH])
def test_wedge_on_specified_tetrahedral_origin_remains_supported(direction):
    mol = Chem.MolFromSmiles("F[C@H](Cl)Br")
    # Bond 1 starts at the specified center; bond 0 starts at fluorine.
    bond = mol.GetBondWithIdx(1)
    bond.SetBondDir(direction)
    before = mol.ToBinary(Chem.PropertyPickleOptions.AllProps)
    common.require_absolute_stereo(mol)
    assert mol.ToBinary(Chem.PropertyPickleOptions.AllProps) == before


@pytest.mark.parametrize("stereo", ["STEREOATROPCW", "STEREOATROPCCW"])
def test_atrop_rejection_precedes_wedge_access_and_later_bonds(stereo):
    # Covers named enum states even on older RDKit builds without atrop enums.
    first = SimpleNamespace(GetStereo=lambda: stereo)
    get_bond = Mock(side_effect=[first, AssertionError("later bond accessed")])
    mol = SimpleNamespace(GetStereoGroups=lambda: (), GetNumBonds=lambda: 2, GetBondWithIdx=get_bond)
    with pytest.raises(ValueError, match="atropisomeric"):
        common.require_absolute_stereo(mol)
    get_bond.assert_called_once_with(0)


def test_nonabsolute_group_rejection_precedes_bond_access():
    group = SimpleNamespace(GetGroupType=lambda: Chem.StereoGroupType.STEREO_OR)
    mol = SimpleNamespace(
        GetStereoGroups=lambda: (group,),
        GetNumBonds=Mock(side_effect=AssertionError("bonds accessed before group check")),
    )
    with pytest.raises(ValueError, match="nonabsolute enhanced stereochemistry"):
        common.require_absolute_stereo(mol)
    mol.GetNumBonds.assert_not_called()


@pytest.mark.parametrize("mol_type", [Chem.Mol, Chem.RWMol])
def test_public_guard_rechecks_current_bond_and_atom_properties(mol_type):
    mol = mol_type(Chem.MolFromSmiles("CC"))
    common.require_absolute_stereo(mol)
    bond = mol.GetBondWithIdx(0)
    bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
    with pytest.raises(ValueError, match="unresolved.*wedging"):
        common.require_absolute_stereo(mol)
    bond.GetBeginAtom().SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    common.require_absolute_stereo(mol)
    bond.GetBeginAtom().SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    with pytest.raises(ValueError, match="unresolved.*wedging"):
        common.require_absolute_stereo(mol)
    bond.SetBondDir(Chem.BondDir.NONE)
    common.require_absolute_stereo(mol)
    if hasattr(Chem.BondStereo, "STEREOATROPCW"):
        bond.SetStereo(Chem.BondStereo.STEREOATROPCW)
        with pytest.raises(ValueError, match="atropisomeric"):
            common.require_absolute_stereo(mol)
        bond.SetStereo(Chem.BondStereo.STEREONONE)
        common.require_absolute_stereo(mol)


def test_public_guard_rechecks_mutable_graph_size_and_stereo_groups():
    mol = Chem.RWMol(Chem.MolFromSmiles("F[C@H](Cl)Br"))
    common.require_absolute_stereo(mol)
    group = Chem.CreateStereoGroup(Chem.StereoGroupType.STEREO_OR, mol, [1])
    mol.SetStereoGroups([group])
    with pytest.raises(ValueError, match="nonabsolute enhanced stereochemistry"):
        common.require_absolute_stereo(mol)
    mol.SetStereoGroups([])
    common.require_absolute_stereo(mol)
    new_index = mol.AddAtom(Chem.Atom(6))
    mol.AddBond(0, new_index, Chem.BondType.SINGLE)
    mol.GetBondBetweenAtoms(0, new_index).SetBondDir(Chem.BondDir.BEGINDASH)
    with pytest.raises(ValueError, match="unresolved.*wedging"):
        common.require_absolute_stereo(mol)
    mol.RemoveAtom(new_index)
    common.require_absolute_stereo(mol)
