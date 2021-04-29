import copy
import oddt
from oddt import shape
import itertools
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import scaffoldgraph as sg
from rdkit import RDLogger
import os
import warnings

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
warnings.filterwarnings('ignore')

def calc_shape(smiles1, smiles2):
    query_mol = oddt.toolkit.readstring('smi', smiles1)
    query_mol2 = oddt.toolkit.readstring('smi', smiles2)
    query_mol.make3D()
    query_mol2.make3D()
    
    query_shape = shape.usr_cat(query_mol)
    query_shape2 = shape.usr_cat(query_mol2)

    similarity = shape.usr_similarity(query_shape, query_shape2)
    return similarity

def calc_electron_shape(smiles1, smiles2):
    query_mol = oddt.toolkit.readstring('smi', smiles1)
    query_mol2 = oddt.toolkit.readstring('smi', smiles2)
    query_mol.make3D()
    query_mol2.make3D()
    query_shape = shape.electroshape(query_mol)
    query_shape2 = shape.electroshape(query_mol2)

    similarity = shape.usr_similarity(query_shape, query_shape2)
    return similarity

def calc_tanimoto_sim(mol1, mol2):
    fp1 = AllChem.GetMorganFingerprint(mol1,2)
    fp2 = AllChem.GetMorganFingerprint(mol2,2)
    sim = DataStructs.TanimotoSimilarity(fp1, fp2)
    return sim

def get_atom_object_from_idx(mol, atom_idx_list):
    atom_obj_list = []
    for each_idx in atom_idx_list:
        atom_obj_list.append(mol.GetAtomWithIdx(each_idx))
    return atom_obj_list

def get_neighbor_atoms(mol, idx, max_len=1):
    neighbors = mol.GetAtomWithIdx(idx).GetNeighbors()
    idx_list = []
    for each_neighbor in neighbors:
        idx_list.append(each_neighbor.GetIdx())
    
    tmp_idx_list = copy.deepcopy(idx_list)
    for i in range(max_len-1):
        for each_idx in tmp_idx_list:
            neighbors = mol.GetAtomWithIdx(each_idx).GetNeighbors()
            for each_neighbor in neighbors:
                idx_list.append(each_neighbor.GetIdx())
                
        idx_list = list(set(idx_list))
        tmp_idx_list = copy.deepcopy(idx_list)
    
    return idx_list

def get_atoms_from_bond(mol, bond_idx):
    begin_atom = None
    end_atom = None
    for bond in mol.GetBonds():
        if bond.GetIdx() == bond_idx:
            begin_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
    return begin_atom, end_atom

def get_bonds_from_atom_list(mol, atom_list):
    bonds = []
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()
        if begin_atom in atom_list and end_atom in atom_list:
            bonds.append(bond.GetIdx())
    return bonds

def get_submol_from_atom_list(mol, atom_list):
    if len(atom_list) == 1:
        atom_idx = atom_list[0]
        symbol = get_symbol_from_atom(mol, atom_idx)
        tmp_submol = Chem.MolFromSmiles(symbol)
    else:
        bonds = get_bonds_from_atom_list(mol, atom_list)
        tmp_submol = Chem.PathToSubmol(mol, bonds)
    return tmp_submol

def get_symbol_from_atom(mol, atom_idx):
    symbol = mol.GetAtomWithIdx(atom_idx).GetSymbol()
    return symbol

def get_substructure_info(mol, pattern):
    substructure_candidate_atom_list = []
    for each_substr in mol.GetSubstructMatches(pattern):
        substructure_candidate_atom_list.append(each_substr)
    return substructure_candidate_atom_list

def init_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    submolblock = Chem.MolToMolBlock(mol)
    mol = Chem.MolFromMolBlock(submolblock)
    return mol, Chem.MolToSmiles(mol)