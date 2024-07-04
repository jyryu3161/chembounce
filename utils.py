#!/usr/bin/env python3
"""
util codes for ChemBounce
"""
import copy
import oddt
from oddt import shape
import itertools
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from rdkit import DataStructs
import scaffoldgraph as sg
from rdkit import RDLogger
import os
import sys
import warnings
import cats_module
from scipy.spatial.distance import euclidean, cosine
import pubchempy as pcp
import pickle
import math
from rdkit.Chem import Descriptors
# SA score calcualtion function - dependency of importing path by the rdkit version
import rdkit.RDPaths as RDPaths
import rdkit.RDConfig as RDConfig
sys.path.append(RDPaths.RDContribDir)
sys.path.append(os.environ['RDBASE'])
try:
    from Contrib.SA_Score.sascorer import calculateScore as calc_SA
except:
    try:
        from SA_Score.sascorer import calculateScore as calc_SA
    except:
        from moses.metrics import SA as calc_SA
import molvs

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
warnings.filterwarnings('ignore')

#### Reference data ####

PLF_LOC=os.path.split(os.path.abspath(__file__))[0]

# Calling reference data
def call_frag_db():
    fragment_file = os.path.join(PLF_LOC,'data','Scaffolds_processed.txt')
    fragment_pkl_file = os.path.join(PLF_LOC,'data','fragment_data.pickle')
    with open(fragment_pkl_file, 'rb') as f:
        fragments_DB = pickle.load(f)
    return fragment_file, fragment_pkl_file, fragments_DB


# Calling reference data
def _call_frag_db_smi_(fragment_file:str=''):
    if not os.path.isfile(fragment_file):
        fragment_file = os.path.join(PLF_LOC,'data','fragment_data.smi.txt')
        if not os.path.isfile(fragment_file):
            fragment_file = os.path.join(PLF_LOC,'data','Scaffolds_processed.txt')
    with open(fragment_file, 'rb') as f:
        fragments_DB = f.read().decode().splitlines()
    return fragment_file, fragments_DB


#### Sub functions ####

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

def calculate_cats_des(mol1, mol2):
    cats = cats_module.CATS2D()
    cats1 = cats.getCATs2D(mol1)
    cats2 = cats.getCATs2D(mol2)
    
    return euclidean(cats1, cats2)

def get_cid_by_structure(smiles):
    try:
        compounds = pcp.get_compounds(smiles, 'smiles')
        if compounds:
            return compounds[0].cid
        else:
            return 'N/A'
    except Exception as e:
        print(f"Error: {e}")
        return 'N/A' 

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

def get_standard_smiles(smiles):
    std_smi = molvs.standardize_smiles(smiles)
    smi_val = molvs.validate_smiles(smiles)
    return std_smi, smi_val


def _get_molecular_prop_(mol):
    sa = calc_SA(mol)
    qed = Chem.QED.default(mol)
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    return sa, qed, mw, logp


def _get_molecular_prop_extd_(mol):
    sa, qed, mw, logp = _get_molecular_prop_(mol)
    props = {'sa':sa,'qed':qed,'mw':mw,'logp':logp}
    props['n_hdonor'] = Lipinski.NumHDonors(mol)
    props['n_hacceptor'] = Lipinski.NumHAcceptors(mol)
    return props


#### Counts and numbers for limit on the main fuction ####

# frag_max_n = overall_max_n/frag_n*3
# scaffold_top_n = int(frag_max_n/cand_max_n__rplc*2)
# (cand_max_n__rplc ~ 10)
# # TODO - scaffold_top_n numbers : with Tanimoto threshold?
def _scaffold_no_reassign_(overall_max_n:int=None, # recommend: ~ 10000
                           frag_max_n:int=None, # recommend: ~ 1000
                           scaffold_top_n:int=None, # recommend: ~ 200
                           cand_max_n__rplc:int=None, # recommend: <= 10
                           _merge_structure_top_n_:int=100,
                           frag_n:int=1,
                          ):
    if not _merge_structure_top_n_:
        _merge_structure_top_n_ = 100
    # define overall_max_n and frag_max_n
    # if frag_max_n is defined
    if not overall_max_n and frag_max_n:
        overall_max_n = int(frag_max_n*frag_n*1.5)
    # if overall_max_n is defined
    elif not frag_max_n and overall_max_n:
        frag_max_n = int(overall_max_n/frag_n*1.5)
    # if both overall_max_n and frag_max_n are not defined
    elif not overall_max_n and not frag_max_n:
        if scaffold_top_n and cand_max_n__rplc:
            frag_max_n = int(scaffold_top_n*cand_max_n__rplc)
            overall_max_n = int(frag_max_n*frag_n)
        elif not scaffold_top_n and cand_max_n__rplc:
            overall_max_n = 10000
            frag_max_n = int(overall_max_n/frag_n*1.5)
        elif scaffold_top_n and not cand_max_n__rplc:
            cand_max_n__rplc = 10
            frag_max_n = int(scaffold_top_n*cand_max_n__rplc*2)
            overall_max_n = int(frag_max_n*frag_n*1.5)
        else: # all has not been defined
            overall_max_n = 10000
            frag_max_n = int(overall_max_n/frag_n*1.5)
    if not cand_max_n__rplc and not scaffold_top_n:
        cand_max_n__rplc = 10
        scaffold_top_n = int(frag_max_n/cand_max_n__rplc*3)
    elif not cand_max_n__rplc and scaffold_top_n:
        cand_max_n__rplc = int(frag_max_n/scaffold_top_n*3)
    elif cand_max_n__rplc and not scaffold_top_n:
        scaffold_top_n = int(frag_max_n/cand_max_n__rplc*3)
    return overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_


#### Determination of properties for candidates ####

# thresholds: {METRIC:(MIN_VAL,MAX_VAL)}
def _default_thrs_(lipinski:bool=True):
    thresholds = dict()
    if lipinski:
        thresholds.update({
            'mw':(0.0,500.0),
            'logp':(float('-inf'),5),
            'n_hdonor':(0,5),
            'n_hacceptor':(0,10),
        })
    return thresholds


# props (dict) : {METRIC:VAL}
# thresholds (dict) : {METRTIC:(min,max)}
def determ_props(props:dict,thresholds:dict,w_all_failed_reason:bool=False,ignore_null_val=False):
    determ_results = []
    for metric, (min_val, max_val) in thresholds.items():
        # Not defined
        if metric not in props:
            determ_results.append(f"{metric}:undefined")
            if not w_all_failed_reason:
                return False, determ_results
            continue
        # property value is invalid: nan, None, or else
        target_val = props[metric]
        if type(target_val) not in [float,int]:
            if ignore_null_val:
                continue
            else:
                determ_results.append(f"{metric}:invalid_value:{target_val}")
                if not w_all_failed_reason:
                    return False, determ_results
                continue
        elif math.isnan(target_val):
            if ignore_null_val:
                continue
            else:
                determ_results.append(f"{metric}.nan_value:{target_val}")
                if not w_all_failed_reason:
                    return False, determ_results
                continue
        if type(min_val) in [int,float]:
            if target_val<min_val:
                determ_results.append(f"{metric}.under_minimum:{target_val}<{min_val}")
                if not w_all_failed_reason:
                    return False, determ_results
        if type(max_val) in [int,float]:
            if target_val>max_val:
                determ_results.append(f"{metric}.upper_maximum:{target_val}>{max_val}")
                if not w_all_failed_reason:
                    return False, determ_results
    if determ_results:
        return False, determ_results
    else:
        return True, determ_results
            
