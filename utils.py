#!/usr/bin/env python3
"""
util codes for ChemBounce
"""
import copy
import oddt
from oddt import shape
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from rdkit import DataStructs
from rdkit import RDLogger
import os
import sys
import warnings
import cats_module
from scipy.spatial.distance import euclidean
import pubchempy as pcp
import pickle
import math
from rdkit.Chem import Descriptors
import numpy as np
import pandas as pd
# SA score calcualtion function - dependency of importing path by the rdkit version
import rdkit.RDPaths as RDPaths
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
from multiprocessing import Pool, cpu_count
from functools import partial

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
            

class SmilesInputError(Exception):
    def __init__(self, value, message='Invalid SMILES input for rdkit'):
        self.value = value
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message}\nReceived SMILES input: {self.value}'


#### Fingerprint-based fast similarity search functions ####

def load_fingerprint_db(fp_file=None):
    """
    Load pre-computed fingerprint database
    
    Args:
        fp_file: Path to fingerprint file. If None, use default location.
        
    Returns:
        dict containing fingerprints, smiles, and metadata
    """
    if fp_file is None:
        fp_file = os.path.join(PLF_LOC, 'data', 'scaffold_fingerprints.npz')
    
    if not os.path.exists(fp_file):
        raise FileNotFoundError(f"Fingerprint file not found: {fp_file}")
    
    print(f"Loading fingerprint database from {fp_file}...")
    data = np.load(fp_file, allow_pickle=True)
    
    # Convert to proper format
    fp_data = {
        'fingerprints': data['fingerprints'],
        'smiles': data['smiles'].tolist() if isinstance(data['smiles'], np.ndarray) else data['smiles'],
        'fp_size': int(data['fp_size']),
        'radius': int(data['radius'])
    }
    
    print(f"Loaded {len(fp_data['fingerprints'])} fingerprints")
    return fp_data


def calculate_bulk_tanimoto(query_fp, db_fps, top_n=None, threshold=0.0):
    """
    Calculate Tanimoto similarity between query and database fingerprints using numpy
    
    Args:
        query_fp: Query fingerprint (numpy array)
        db_fps: Database fingerprints (numpy array)
        top_n: Return only top N results
        threshold: Minimum similarity threshold
        
    Returns:
        indices and similarities of matching compounds
    """
    # Ensure query_fp is 1D
    if query_fp.ndim > 1:
        query_fp = query_fp.flatten()
    
    # Calculate Tanimoto similarity using numpy operations
    # Tanimoto = intersection / union = (A & B) / (A | B)
    intersection = np.sum(db_fps & query_fp, axis=1)
    union = np.sum(db_fps | query_fp, axis=1)
    
    # Avoid division by zero
    similarities = np.zeros(len(db_fps))
    non_zero_union = union > 0
    similarities[non_zero_union] = intersection[non_zero_union] / union[non_zero_union]
    
    # Apply threshold
    valid_mask = similarities >= threshold
    valid_indices = np.where(valid_mask)[0]
    valid_similarities = similarities[valid_mask]
    
    # Sort by similarity (descending)
    sorted_indices = np.argsort(valid_similarities)[::-1]
    
    # Apply top_n limit if specified
    if top_n is not None and len(sorted_indices) > top_n:
        sorted_indices = sorted_indices[:top_n]
    
    # Return indices and similarities
    result_indices = valid_indices[sorted_indices]
    result_similarities = valid_similarities[sorted_indices]
    
    return result_indices, result_similarities


def search_similar_scaffolds_fast(query_mol, fp_db, scaffold_top_n=None, threshold=0.3, use_multiprocessing=False, n_jobs=None):
    """
    Fast similarity search using pre-computed fingerprints
    
    Args:
        query_mol: Query molecule (RDKit mol object)
        fp_db: Fingerprint database (from load_fingerprint_db)
        scaffold_top_n: Maximum number of results
        threshold: Minimum similarity threshold
        use_multiprocessing: Whether to use multiprocessing (default: False for backward compatibility)
        n_jobs: Number of processes (None = use all CPUs, -1 = use all CPUs)
        
    Returns:
        pandas Series with SMILES as index and similarity as values
    """
    # Generate fingerprint for query molecule
    query_fp = AllChem.GetMorganFingerprintAsBitVect(
        query_mol, 
        fp_db['radius'], 
        nBits=fp_db['fp_size']
    )
    
    # Convert to numpy array
    query_fp_arr = np.zeros((fp_db['fp_size'],), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(query_fp, query_fp_arr)
    
    # Get query SMILES for filtering
    query_smiles = Chem.MolToSmiles(query_mol)
    query_len = len(query_smiles)
    
    # Pre-filter by SMILES length (same logic as original)
    db_fps = fp_db['fingerprints']
    db_smiles = fp_db['smiles']
    
    # Create length filter mask
    length_mask = np.ones(len(db_smiles), dtype=bool)
    for i, smiles in enumerate(db_smiles):
        smiles_len = len(smiles)
        if smiles_len > query_len * 2 or query_len > smiles_len * 2:
            length_mask[i] = False
    
    # Apply length filter
    filtered_indices = np.where(length_mask)[0]
    if len(filtered_indices) == 0:
        return pd.Series(dtype=float)
    
    filtered_fps = db_fps[filtered_indices]
    
    # Calculate similarities
    if use_multiprocessing:
        indices, similarities = calculate_bulk_tanimoto_multiprocessing(
            query_fp_arr, 
            filtered_fps, 
            top_n=scaffold_top_n,
            threshold=max(threshold, 0.1),  # At least 0.1 as in original
            n_jobs=n_jobs
        )
    else:
        indices, similarities = calculate_bulk_tanimoto(
            query_fp_arr, 
            filtered_fps, 
            top_n=scaffold_top_n,
            threshold=max(threshold, 0.1)  # At least 0.1 as in original
        )
    
    # Map back to original indices
    original_indices = filtered_indices[indices]
    
    # Create result series
    result_smiles = [db_smiles[i] for i in original_indices]
    result_series = pd.Series(similarities, index=result_smiles)
    
    # Remove exact match (similarity = 1.0)
    result_series = result_series[result_series < 1.0]
    
    # Apply final threshold
    if threshold:
        result_series = result_series[result_series > threshold]
    
    # Apply top_n limit
    if scaffold_top_n and len(result_series) > scaffold_top_n:
        result_series = result_series.iloc[:scaffold_top_n]
    
    result_series.name = 'Tanimoto Similarity'
    return result_series


def _calculate_tanimoto_chunk(args):
    """Helper function for multiprocessing: calculates Tanimoto similarity for a chunk"""
    query_fp, db_fps_chunk, start_idx, threshold = args
    
    # Calculate Tanimoto similarity using numpy operations
    intersection = np.sum(db_fps_chunk & query_fp, axis=1)
    union = np.sum(db_fps_chunk | query_fp, axis=1)
    
    # Avoid division by zero
    similarities = np.zeros(len(db_fps_chunk))
    non_zero_union = union > 0
    similarities[non_zero_union] = intersection[non_zero_union] / union[non_zero_union]
    
    # Apply threshold and return valid results with adjusted indices
    valid_mask = similarities >= threshold
    valid_indices = np.where(valid_mask)[0]
    
    # Adjust indices to original position
    adjusted_indices = valid_indices + start_idx
    valid_similarities = similarities[valid_mask]
    
    return adjusted_indices, valid_similarities


def calculate_bulk_tanimoto_multiprocessing(query_fp, db_fps, top_n=None, threshold=0.0, n_jobs=None):
    """
    Calculate Tanimoto similarity using multiprocessing
    
    Args:
        query_fp: Query fingerprint (numpy array)
        db_fps: Database fingerprints (numpy array)
        top_n: Return only top N results
        threshold: Minimum similarity threshold
        n_jobs: Number of processes (None = use all CPUs)
        
    Returns:
        indices and similarities of matching compounds
    """
    # Ensure query_fp is 1D
    if query_fp.ndim > 1:
        query_fp = query_fp.flatten()
    
    # Determine number of processes
    if n_jobs is None:
        n_jobs = cpu_count()
    elif n_jobs == -1:
        n_jobs = cpu_count()
    else:
        n_jobs = min(n_jobs, cpu_count())
    
    # For small datasets, use single process
    if len(db_fps) < 1000 or n_jobs == 1:
        return calculate_bulk_tanimoto(query_fp, db_fps, top_n, threshold)
    
    # Split data into chunks
    chunk_size = max(1, len(db_fps) // n_jobs)
    chunks = []
    for i in range(0, len(db_fps), chunk_size):
        end_idx = min(i + chunk_size, len(db_fps))
        chunks.append((query_fp, db_fps[i:end_idx], i, threshold))
    
    # Process chunks in parallel
    with Pool(n_jobs) as pool:
        results = pool.map(_calculate_tanimoto_chunk, chunks)
    
    # Combine results
    all_indices = []
    all_similarities = []
    for indices, similarities in results:
        all_indices.extend(indices)
        all_similarities.extend(similarities)
    
    if not all_indices:
        return np.array([]), np.array([])
    
    # Convert to numpy arrays
    all_indices = np.array(all_indices)
    all_similarities = np.array(all_similarities)
    
    # Sort by similarity (descending)
    sorted_idx = np.argsort(all_similarities)[::-1]
    
    # Apply top_n limit if specified
    if top_n is not None and len(sorted_idx) > top_n:
        sorted_idx = sorted_idx[:top_n]
    
    # Return sorted results
    result_indices = all_indices[sorted_idx]
    result_similarities = all_similarities[sorted_idx]
    
    return result_indices, result_similarities


def search_similar_scaffolds_fast_mp(query_mol, fp_db, scaffold_top_n=None, threshold=0.3, use_multiprocessing=True, n_jobs=None):
    """
    Fast similarity search using pre-computed fingerprints with optional multiprocessing
    
    Args:
        query_mol: Query molecule (RDKit mol object)
        fp_db: Fingerprint database (from load_fingerprint_db)
        scaffold_top_n: Maximum number of results
        threshold: Minimum similarity threshold
        use_multiprocessing: Whether to use multiprocessing (default: True)
        n_jobs: Number of processes (None = use all CPUs, -1 = use all CPUs)
        
    Returns:
        pandas Series with SMILES as index and similarity as values
    """
    # Generate fingerprint for query molecule
    query_fp = AllChem.GetMorganFingerprintAsBitVect(
        query_mol, 
        fp_db['radius'], 
        nBits=fp_db['fp_size']
    )
    
    # Convert to numpy array
    query_fp_arr = np.zeros((fp_db['fp_size'],), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(query_fp, query_fp_arr)
    
    # Get query SMILES for filtering
    query_smiles = Chem.MolToSmiles(query_mol)
    query_len = len(query_smiles)
    
    # Pre-filter by SMILES length (same logic as original)
    db_fps = fp_db['fingerprints']
    db_smiles = fp_db['smiles']
    
    # Create length filter mask (can be parallelized for very large datasets)
    length_mask = np.ones(len(db_smiles), dtype=bool)
    for i, smiles in enumerate(db_smiles):
        smiles_len = len(smiles)
        if smiles_len > query_len * 2 or query_len > smiles_len * 2:
            length_mask[i] = False
    
    # Apply length filter
    filtered_indices = np.where(length_mask)[0]
    if len(filtered_indices) == 0:
        return pd.Series(dtype=float)
    
    filtered_fps = db_fps[filtered_indices]
    
    # Calculate similarities with or without multiprocessing
    if use_multiprocessing:
        indices, similarities = calculate_bulk_tanimoto_multiprocessing(
            query_fp_arr, 
            filtered_fps, 
            top_n=scaffold_top_n,
            threshold=max(threshold, 0.1),  # At least 0.1 as in original
            n_jobs=n_jobs
        )
    else:
        indices, similarities = calculate_bulk_tanimoto(
            query_fp_arr, 
            filtered_fps, 
            top_n=scaffold_top_n,
            threshold=max(threshold, 0.1)
        )
    
    # Map back to original indices
    original_indices = filtered_indices[indices]
    
    # Create result series
    result_smiles = [db_smiles[i] for i in original_indices]
    result_series = pd.Series(similarities, index=result_smiles)
    
    # Remove exact match (similarity = 1.0)
    result_series = result_series[result_series < 1.0]
    
    # Apply final threshold
    if threshold:
        result_series = result_series[result_series > threshold]
    
    # Apply top_n limit
    if scaffold_top_n and len(result_series) > scaffold_top_n:
        result_series = result_series.iloc[:scaffold_top_n]
    
    result_series.name = 'Tanimoto Similarity'
    return result_series


def generate_query_fingerprint(mol, radius=2, fp_size=2048):
    """
    Generate fingerprint for a query molecule
    
    Args:
        mol: RDKit mol object
        radius: Morgan fingerprint radius
        fp_size: Size of fingerprint
        
    Returns:
        numpy array of fingerprint
    """
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=fp_size)
    arr = np.zeros((fp_size,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr
