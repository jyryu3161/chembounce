import os
import glob
import warnings
import time
import logging
import numpy as np
import math
import time
import argparse
from multiprocessing import Process, Queue
import utils

from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Descriptors
from moses.metrics import SA
from rdkit.Chem.Descriptors import qed
from tqdm import tqdm

import numpy as np
from numpy import dot
from numpy.linalg import norm
from scipy.spatial import distance
from moses.metrics import mol_passes_filters, QED, SA, logP
from molvs.fragment import LargestFragmentChooser
from molvs import standardize_smiles
from molvs import validate_smiles
from rdkit.Chem.Scaffolds import MurckoScaffold

def get_rules(mol):
    choices = ['#6', '#7', '#8', '#9', '#16', '#17', '#35']
    rxn_smart_list = []
    for X in choices:
        for Y in choices:
            if mol.HasSubstructMatch(Chem.MolFromSmarts('[' + X + ']')):
                if mol.HasSubstructMatch(Chem.MolFromSmarts('[' + Y + ']')):
                    if X != Y:
                        tmp_smart = '[X:1]>>[Y:1]'.replace('X', X).replace('Y', Y)
                        rxn_smart_list.append(tmp_smart)
    return rxn_smart_list

def run_change_atoms(query_smiles):
    candidate_smiles_list = []
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol != None:
        rxn_smart_list = get_rules(query_mol)

        for rxn_smart in rxn_smart_list:
            rxn = AllChem.ReactionFromSmarts(rxn_smart)
            new_mol_trial = rxn.RunReactants((query_mol,))
            for results in new_mol_trial:
                mol = results[0]
                if mol != None:
                    tmp_smiles = Chem.MolToSmiles(mol)
                    tmp_mol = Chem.MolFromSmiles(tmp_smiles)
                    if tmp_mol != None:
                        validate_results = validate_smiles(tmp_smiles)
                        if len(validate_results) == 0:
                            candidate_smiles_list.append(tmp_smiles)

        candidate_smiles_list = list(set(candidate_smiles_list))
    return candidate_smiles_list

def run_transformation(query_smiles, cycle_n=2):
    candidate_smiles_list_set1 = [query_smiles]
    for cycle in range(cycle_n):
        for i in range(len(candidate_smiles_list_set1)):
            tmp_query = candidate_smiles_list_set1[i]
            tmp_set = run_change_atoms(tmp_query)
            candidate_smiles_list_set1+=tmp_set
        candidate_smiles_list_set1 = list(set(candidate_smiles_list_set1))
    candidate_smiles_list = list(set(candidate_smiles_list_set1))
    final_candidates = check_isosteres(query_smiles, candidate_smiles_list)
    return final_candidates
    
def check_isosteres(query_smiles, smiles_list):
    final_candidates = []
    
    m = Chem.MolFromSmiles(query_smiles)
    scaffold = MurckoScaffold.GetScaffoldForMol(m)
    scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    scaffold_smiles = Chem.MolToSmiles(scaffold)

    for smiles in smiles_list:
        tmp_mol = Chem.MolFromSmiles(smiles)
        tmp_scaffold = MurckoScaffold.GetScaffoldForMol(tmp_mol)
        tmp_scaffold = MurckoScaffold.MakeScaffoldGeneric(tmp_scaffold)
        tmp_scaffold_smiles = Chem.MolToSmiles(tmp_scaffold)
        if scaffold_smiles == tmp_scaffold_smiles:
            final_candidates.append(smiles)
    return final_candidates

    
    # total_smiles = run_transformation(target_smiles, cycle_num)


