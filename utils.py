#!/usr/bin/env python3
"""
util codes for ChemBounce
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from rdkit import DataStructs
import os
import cats_module
from scipy.spatial.distance import euclidean
import math
from rdkit.Chem import Descriptors
import numpy as np
from rdkit.Contrib.SA_Score.sascorer import calculateScore as calc_SA
from rdkit.Chem.MolStandardize import rdMolStandardize
from functools import lru_cache
from numbers import Real
from chemistry import count_fingerprint

#### Reference data ####

PLF_LOC = os.path.split(os.path.abspath(__file__))[0]


# Calling reference data
def call_frag_db():
    """
    Load fragment database from SMILES (legacy return tuple retained).

    Returns:
        tuple: (fragment_file_path, fragment_pkl_file_path, fragments_DB)
            - fragment_file_path: Path to text file with scaffolds
            - fragment_pkl_file_path: Path to pickle file
            - fragments_DB: List of loaded fragment molecules
    """
    fragment_file = os.path.join(PLF_LOC, "data", "Scaffolds_processed.txt")
    fragment_pkl_file = os.path.join(PLF_LOC, "data", "fragment_data.pickle")
    # Legacy API keeps its tuple shape but loads transparent SMILES, never pickle.
    with open(fragment_file) as stream:
        fragments_DB = [init_mol(line.strip())[0] for line in stream if line.strip()]
    return fragment_file, fragment_pkl_file, fragments_DB


# Calling reference data
def _call_frag_db_smi_(fragment_file: str = ""):
    """
    Load fragment database as SMILES strings for low memory mode.

    Args:
        fragment_file: Path to fragment file (optional, uses default if empty)

    Returns:
        tuple: (fragment_file_path, fragments_DB)
            - fragment_file_path: Path to the fragment file used
            - fragments_DB: List of SMILES strings
    """
    if not os.path.isfile(fragment_file):
        fragment_file = os.path.join(PLF_LOC, "data", "fragment_data.smi.txt")
        if not os.path.isfile(fragment_file):
            fragment_file = os.path.join(PLF_LOC, "data", "Scaffolds_processed.txt")
    with open(fragment_file, "rb") as f:
        fragments_DB = f.read().decode().splitlines()
    return fragment_file, fragments_DB


#### Sub functions ####


@lru_cache(maxsize=256)
def shape_descriptor(smiles, seed=42, method="electroshape"):
    """Seeded ETKDGv3/UFF conformer; heavy-atom shape descriptor.

    A single conformer is an approximate ranking descriptor, not evidence of
    binding affinity. Failures are reported by the caller, never scored as zero.
    """
    import shape_descriptors as shape

    mol, canonical = init_mol(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    if AllChem.EmbedMolecule(mol, params) < 0:
        raise ValueError("3D embedding failed")
    if not AllChem.UFFHasAllMoleculeParams(mol):
        raise ValueError("UFF parameters unavailable")
    status = AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
    if status != 0:
        raise ValueError("UFF optimization did not converge")
    heavy = Chem.RemoveHs(mol)
    descriptor = shape.electroshape(heavy) if method == "electroshape" else shape.usrcat(heavy)
    if not np.isfinite(descriptor).all():
        raise ValueError("Non-finite shape descriptor")
    descriptor.setflags(write=False)
    return descriptor


def calc_shape(smiles1, smiles2, seed=42):
    import shape_descriptors as shape

    return float(
        shape.similarity(shape_descriptor(smiles1, seed, "usrcat"), shape_descriptor(smiles2, seed, "usrcat"))
    )


def calc_electron_shape(smiles1, smiles2, seed=42):
    import shape_descriptors as shape

    return float(shape.similarity(shape_descriptor(smiles1, seed), shape_descriptor(smiles2, seed)))


def calc_tanimoto_sim(mol1, mol2):
    """
    Calculate Tanimoto similarity between two molecules using Morgan fingerprints.

    Args:
        mol1: RDKit molecule object
        mol2: RDKit molecule object

    Returns:
        float: Tanimoto similarity score (0-1)
    """
    fp1 = count_fingerprint(mol1)
    fp2 = count_fingerprint(mol2)
    sim = DataStructs.TanimotoSimilarity(fp1, fp2)
    return sim


def calculate_cats_des(mol1, mol2):
    """
    Calculate CATS (Chemically Advanced Template Search) descriptor distance.

    Args:
        mol1: RDKit molecule object
        mol2: RDKit molecule object

    Returns:
        float: Euclidean distance between CATS descriptors
    """
    cats = cats_module.CATS2D()
    cats1 = cats.getCATs2D(mol1)
    cats2 = cats.getCATs2D(mol2)

    return euclidean(cats1, cats2)


def get_cid_by_structure(smiles):
    """
    Get PubChem CID (Compound ID) from SMILES structure.

    Args:
        smiles: SMILES string of molecule

    Returns:
        str: PubChem CID or 'N/A' if not found
    """
    import pubchempy as pcp

    try:
        compounds = pcp.get_compounds(smiles, "smiles")
        if compounds:
            return compounds[0].cid
        else:
            return "N/A"
    except Exception as e:
        print(f"Error: {e}")
        return "N/A"


def get_atom_object_from_idx(mol, atom_idx_list):
    """
    Get atom objects from a list of atom indices.

    Args:
        mol: RDKit molecule object
        atom_idx_list: List of atom indices

    Returns:
        list: List of RDKit atom objects
    """
    atom_obj_list = []
    for each_idx in atom_idx_list:
        atom_obj_list.append(mol.GetAtomWithIdx(each_idx))
    return atom_obj_list


def get_neighbor_atoms(mol, idx, max_len=1):
    """
    Get neighbor atoms within specified distance from a given atom.

    Args:
        mol: RDKit molecule object
        idx: Atom index to start from
        max_len: Maximum distance in bonds (default: 1)

    Returns:
        list: List of atom indices within max_len bonds
    """
    if not isinstance(max_len, int) or max_len < 0:
        raise ValueError("max_len must be a nonnegative integer")
    mol.GetAtomWithIdx(idx)  # validate the starting index
    visited, frontier = {idx}, {idx}
    for _ in range(max_len):
        frontier = {
            neighbor.GetIdx() for index in frontier for neighbor in mol.GetAtomWithIdx(index).GetNeighbors()
        } - visited
        visited.update(frontier)
    return sorted(visited - {idx})


def get_atoms_from_bond(mol, bond_idx):
    """
    Get the begin and end atom indices of a bond.

    Args:
        mol: RDKit molecule object
        bond_idx: Bond index

    Returns:
        tuple: (begin_atom_idx, end_atom_idx) or (None, None) if not found
    """
    begin_atom = None
    end_atom = None
    for bond in mol.GetBonds():
        if bond.GetIdx() == bond_idx:
            begin_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
    return begin_atom, end_atom


def get_bonds_from_atom_list(mol, atom_list):
    """
    Get all bonds between atoms in the given list.

    Args:
        mol: RDKit molecule object
        atom_list: List of atom indices

    Returns:
        list: List of bond indices
    """
    bonds = []
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()
        if begin_atom in atom_list and end_atom in atom_list:
            bonds.append(bond.GetIdx())
    return bonds


def get_submol_from_atom_list(mol, atom_list):
    """
    Extract a submolecule from a list of atoms.

    Args:
        mol: RDKit molecule object
        atom_list: List of atom indices

    Returns:
        RDKit molecule object: Submolecule containing specified atoms
    """
    keep = set(atom_list)
    if any(i < 0 or i >= mol.GetNumAtoms() for i in keep):
        raise IndexError("Atom index out of range")
    edit = Chem.RWMol(mol)
    for index in reversed(range(mol.GetNumAtoms())):
        if index not in keep:
            edit.RemoveAtom(index)
    return edit.GetMol()


def get_symbol_from_atom(mol, atom_idx):
    """
    Get atomic symbol from atom index.

    Args:
        mol: RDKit molecule object
        atom_idx: Atom index

    Returns:
        str: Atomic symbol (e.g., 'C', 'N', 'O')
    """
    symbol = mol.GetAtomWithIdx(atom_idx).GetSymbol()
    return symbol


def get_substructure_info(mol, pattern):
    """
    Find all substructure matches of pattern in molecule.

    Args:
        mol: RDKit molecule object to search in
        pattern: RDKit molecule object representing the pattern

    Returns:
        list: List of tuples containing atom indices for each match
    """
    substructure_candidate_atom_list = []
    for each_substr in mol.GetSubstructMatches(pattern):
        substructure_candidate_atom_list.append(each_substr)
    return substructure_candidate_atom_list


def init_mol(smiles):
    """
    Initialize molecule from SMILES with standardization.

    Args:
        smiles: SMILES string

    Returns:
        tuple: (RDKit molecule object, canonical SMILES string)
    """
    if not isinstance(smiles, str) or not smiles.strip():
        raise SmilesInputError(smiles, "SMILES must be a nonempty string")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None or not mol.GetNumAtoms():
        raise SmilesInputError(smiles)
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise SmilesInputError(smiles, "Dummy atoms are not supported in input structures")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol, Chem.MolToSmiles(mol, isomericSmiles=True)


def get_standard_smiles(smiles):
    mol, _ = init_mol(smiles)
    standardized = rdMolStandardize.Cleanup(mol)
    return Chem.MolToSmiles(standardized, isomericSmiles=True), []


def _get_molecular_prop_(mol):
    """
    Calculate basic molecular properties.

    Args:
        mol: RDKit molecule object

    Returns:
        tuple: (SA_score, QED, molecular_weight, logP)
    """
    sa = calc_SA(mol)
    qed = Chem.QED.default(mol)
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    return sa, qed, mw, logp


def _get_molecular_prop_extd_(mol):
    """
    Calculate extended molecular properties including H-bond donors/acceptors.

    Args:
        mol: RDKit molecule object

    Returns:
        dict: Dictionary with keys 'sa', 'qed', 'mw', 'logp', 'n_hdonor', 'n_hacceptor'
    """
    sa, qed, mw, logp = _get_molecular_prop_(mol)
    props = {"sa": sa, "qed": qed, "mw": mw, "logp": logp}
    props["n_hdonor"] = Lipinski.NumHDonors(mol)
    props["n_hacceptor"] = Lipinski.NumHAcceptors(mol)
    return props


#### Counts and numbers for limit on the main fuction ####


def _scaffold_no_reassign_(
    overall_max_n: int = None,  # recommend: ~ 10000
    frag_max_n: int = None,  # recommend: ~ 1000
    scaffold_top_n: int = None,  # recommend: ~ 200
    cand_max_n__rplc: int = None,  # recommend: <= 10
    _merge_structure_top_n_: int = 100,
    frag_n: int = 1,
):
    """
    Reassign iteration limits for scaffold search based on provided parameters.

    Args:
        overall_max_n: Maximum total candidates (recommend: ~10000)
        frag_max_n: Maximum candidates per fragment (recommend: ~1000)
        scaffold_top_n: Top N scaffolds to test (recommend: ~200)
        cand_max_n__rplc: Maximum candidates per replacement (recommend: <=10)
        _merge_structure_top_n_: Top N structures for merging (default: 100)
        frag_n: Number of fragments

    Returns:
        tuple: (overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_)
    """
    for name, val in locals().copy().items():
        if val is not None and (
            not isinstance(val, int) or isinstance(val, bool) or val < (0 if name == "frag_n" else 1)
        ):
            raise ValueError(f"{name} must be a positive integer (frag_n may be zero)")
    frag_max_n = frag_max_n or 1000
    overall_max_n = overall_max_n or frag_max_n * max(1, frag_n)
    cand_max_n__rplc = cand_max_n__rplc or 10
    scaffold_top_n = scaffold_top_n or max(1, math.ceil(frag_max_n / cand_max_n__rplc) * 3)
    return overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_ or 100


#### Determination of properties for candidates ####


# thresholds: {METRIC:(MIN_VAL,MAX_VAL)}
def _default_thrs_(lipinski: bool = True):
    """
    Get default molecular property thresholds.

    Args:
        lipinski: If True, apply Lipinski's rule of five thresholds

    Returns:
        dict: Dictionary of thresholds {metric: (min_value, max_value)}
    """
    thresholds = dict()
    if lipinski:
        thresholds.update(
            {
                "mw": (0.0, 500.0),
                "logp": (float("-inf"), 5),
                "n_hdonor": (0, 5),
                "n_hacceptor": (0, 10),
            }
        )
    return thresholds


# props (dict) : {METRIC:VAL}
# thresholds (dict) : {METRTIC:(min,max)}
def determ_props(props: dict, thresholds: dict, w_all_failed_reason: bool = False, ignore_null_val=False):
    """
    Determine if molecular properties pass given thresholds.

    Args:
        props: Dictionary of molecular properties {metric: value}
        thresholds: Dictionary of thresholds {metric: (min, max)}
        w_all_failed_reason: If True, return all failed reasons; if False, stop at first failure
        ignore_null_val: If True, ignore invalid/NaN values

    Returns:
        tuple: (passed: bool, failure_reasons: list)
    """
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
        if not isinstance(target_val, Real):
            if ignore_null_val:
                continue
            else:
                determ_results.append(f"{metric}:invalid_value:{target_val}")
                if not w_all_failed_reason:
                    return False, determ_results
                continue
        elif not math.isfinite(target_val):
            if ignore_null_val:
                continue
            else:
                determ_results.append(f"{metric}.nan_value:{target_val}")
                if not w_all_failed_reason:
                    return False, determ_results
                continue
        if isinstance(min_val, Real):
            if target_val < min_val:
                determ_results.append(f"{metric}.under_minimum:{target_val}<{min_val}")
                if not w_all_failed_reason:
                    return False, determ_results
        if isinstance(max_val, Real):
            if target_val > max_val:
                determ_results.append(f"{metric}.upper_maximum:{target_val}>{max_val}")
                if not w_all_failed_reason:
                    return False, determ_results
    if determ_results:
        return False, determ_results
    else:
        return True, determ_results


class SmilesInputError(Exception):
    def __init__(self, value, message="Invalid SMILES input for rdkit"):
        self.value = value
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}\nReceived SMILES input: {self.value}"


# Public compatibility imports; implementation lives in the dedicated module.
from fingerprint_db import (
    load_fingerprint_db as load_fingerprint_db,
    calculate_bulk_tanimoto as calculate_bulk_tanimoto,
    search_similar_scaffolds_fast as search_similar_scaffolds_fast,
    calculate_bulk_tanimoto_multiprocessing as calculate_bulk_tanimoto_multiprocessing,
    search_similar_scaffolds_fast_mp as search_similar_scaffolds_fast_mp,
    generate_query_fingerprint as generate_query_fingerprint,
)
