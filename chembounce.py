#!/usr/bin/env python3
"""
ChemBounce main functions and subfuctions
"""
import logging
import utils
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import AllChem
import tqdm
import gc
import numpy as np
import pandas as pd

import scaffoldgraph as sg
import itertools
import copy
import os
import sys
import warnings
import pickle
import copy
import datetime


# Removal of substructures for given molecule
def remove_substructures(mol, pattern_mol):
    substructure_list = utils.get_substructure_info(mol, pattern_mol)
    candidate_structures = []
    
    for each_structure_atom_list in substructure_list:
        temp_mol = copy.deepcopy(mol)
        temp_mol2 = copy.deepcopy(mol)
        
        edit_mol = Chem.EditableMol(temp_mol)
        edit_mol2 = Chem.EditableMol(temp_mol2)
        
        original_atom_info = {}
        rm_atom_idx_list = []
        for each_atom_idx in sorted(each_structure_atom_list, reverse=True):
            neighbor_atoms = utils.get_neighbor_atoms(mol, each_atom_idx)
            if set(neighbor_atoms).issubset(set(each_structure_atom_list)) == False:
                linker_atom = list(set(neighbor_atoms).difference(set(each_structure_atom_list)))[0]
                original_atom_info[linker_atom] = each_atom_idx
            else:
                rm_atom_idx_list.append(each_atom_idx)
        
        for i in each_structure_atom_list:
            for j in each_structure_atom_list:
                if i != j:
                    if temp_mol.GetBondBetweenAtoms(i, j) != None:
                        edit_mol.RemoveBond(i, j)
        
        for each_atom_idx in rm_atom_idx_list:
            edit_mol.RemoveAtom(each_atom_idx)
        
        temp_mol = edit_mol.GetMol()
        Chem.Kekulize(temp_mol)
        smiles = Chem.MolToSmiles(temp_mol, False, kekuleSmiles=True)
        smiles_list = smiles.split('.')
        for each_frag in smiles_list:
            temp_pattern_mol = Chem.MolFromSmiles(each_frag)
            
            temp_substructure_list = utils.get_substructure_info(mol, temp_pattern_mol)
            temp_substructure_list = list(temp_substructure_list[0])
            target_idx = None
            for each_atom in temp_substructure_list:
                if each_atom in original_atom_info:
                    temp_substructure_list.append(original_atom_info[each_atom])
                    target_idx = original_atom_info[each_atom]
            
            edit_mol = Chem.EditableMol(mol)
            edit_mol.ReplaceAtom(target_idx, Chem.Atom(0))
            submol = edit_mol.GetMol()

            submol = utils.get_submol_from_atom_list(submol, temp_substructure_list)  
            frag_smiles = Chem.MolToSmiles(submol)
            candidate_structures.append(frag_smiles)
            
        break
    return candidate_structures


def merge_structures(template_mol, mol, frag_mol, top_n=3):
    template_smiles = Chem.MolToSmiles(template_mol)
    substitution_smiles = '[*]'
    substitution_mol = Chem.MolFromSmiles(substitution_smiles)
    
    combo = Chem.CombineMols(mol, frag_mol)
    
    substructure_list = utils.get_substructure_info(combo, substitution_mol)
    target_pos_atom = substructure_list[0][0] 
    
    neighbor_atoms = utils.get_neighbor_atoms(combo, target_pos_atom, max_len=1)
    neighbor_atom = neighbor_atoms[0]
    
    bonds = utils.get_bonds_from_atom_list(combo, [target_pos_atom, neighbor_atom])
    target_bond = bonds[0]
    target_bond_obj = combo.GetBondWithIdx(target_bond)
    target_bond_type = target_bond_obj.GetBondType()
    
    structure_candidates = []
    substructure_list = utils.get_substructure_info(combo, mol)
#     for each_atom in substructure_list[0]:
    for each_atom in substructure_list[0]:
        symbol = combo.GetAtomWithIdx(each_atom).GetSymbol()
        
        combo = copy.deepcopy(combo)
        edit_combo_mol = Chem.EditableMol(combo)        
        edit_combo_mol.AddBond(neighbor_atom, each_atom, order=target_bond_type)
        edit_combo_mol.RemoveAtom(target_pos_atom)
        submol = edit_combo_mol.GetMol()
        
        submol_smiles = Chem.MolToSmiles(submol)
        tmp_submol = Chem.MolFromSmiles(submol_smiles)

        if tmp_submol == None:
            try:
                submolblock = Chem.MolToMolBlock(submol)
            except:
                continue
            submol = Chem.MolFromMolBlock(submolblock)
        else:
            submol = tmp_submol

        if submol != None:
            submol_smiles = Chem.MolToSmiles(submol)
            
            sim = utils.calc_tanimoto_sim(template_mol, submol)
               # sim = utils.calc_electron_shape(template_smiles, submol_smiles)
            structure_candidates.append([sim, submol_smiles])
    
    structure_candidates.sort(reverse=True)
    return structure_candidates

def replace_molecule(target_mol, pattern_mol, replace_mol, max_n, _merge_structure_top_n_):
    candidate_structures = remove_substructures(target_mol, pattern_mol)
    
    final_candidates = []
    saved_structures = {}
    # Permutation for merge order for fragment remnants
    candidate_sets = itertools.permutations(candidate_structures)

    # code_ rearr. : unlimit ub/ 
    for each_set in candidate_sets: # permut round
        candidate_set = list(each_set)
        start_replace_scaffold_list = [(1.,replace_mol)] # seed mol : reset for each_set
        for each_structure in candidate_set: # frag round
            frag_mol = Chem.MolFromSmiles(each_structure)
            tmp_replace_scaffold_list = []
            for _score_,each_replace_mol in start_replace_scaffold_list[:max_n]:
                structure_candidates = merge_structures(
                    target_mol, each_replace_mol, frag_mol, _merge_structure_top_n_)
                for each_candidate in structure_candidates:
                    each_score = each_candidate[0]
                    each_str = each_candidate[1]
                    if each_str not in saved_structures:
                        saved_structures[each_str] = True
                        each_str_mol, each_str = utils.init_mol(each_str)
                        tmp_replace_scaffold_list.append(tuple([each_score,each_str_mol]))
            tmp_replace_scaffold_list = list(set(tmp_replace_scaffold_list))
            tmp_replace_scaffold_list.sort(reverse=True,key=lambda x: x[0])
            start_replace_scaffold_list = tmp_replace_scaffold_list # update seed mole
        final_candidates.extend(
            [tuple([score,Chem.MolToSmiles(mol)]) for (score, mol) in start_replace_scaffold_list])
    
    final_candidates.sort(reverse=True)
    final_candidates = final_candidates[0:max_n]
    return final_candidates


def read_fragments(target_smiles, fragment_file):
    fragments = []           
    num_lines = sum(1 for line in open(fragment_file,'r'))
    with open(fragment_file,'r') as f:
        for line in tqdm.tqdm(f, total=num_lines):
            line = line.strip()
            if '.' not in line.strip():
                if len(line.strip()) < len(target_smiles):
                    smiles = line.strip()
                    mol = Chem.MolFromSmiles(smiles)
                    fragments.append(mol)
    return fragments

def search_similar_scaffolds(original_scaffold, fragments_DB,
#                              scaffold_top_n:int, threshold:int=None,
                             low_mem:bool=False, tqdm_quiet:bool=False):
    original_scaffold_smiles = Chem.MolToSmiles(original_scaffold)
    
    scaffold_scores = dict()
    for each_frag_candidate in tqdm.tqdm(fragments_DB,desc=original_scaffold_smiles,disable=tqdm_quiet):
        if low_mem: # candidate of SMILES obj.
            candidate_smiles = each_frag_candidate
            candidate_mol = Chem.MolFromSmiles(candidate_smiles)
        else: # Candidate of mol obj.
            candidate_mol = each_frag_candidate
            candidate_smiles = Chem.MolToSmiles(each_frag_candidate)
        if len(candidate_smiles) > len(original_scaffold_smiles)*2:
            continue
        if len(original_scaffold_smiles) > len(candidate_smiles)*2:
            continue
        sim = utils.calc_tanimoto_sim(original_scaffold, candidate_mol)
        if sim > 0.1:
            if sim != 1.0:
                scaffold_scores[candidate_smiles] = sim
    scaffold_scores = pd.Series(scaffold_scores)
    scaffold_scores.sort_values(ascending=False)
    return scaffold_scores

    
# Main function
def chembounce(target_smiles:str,
               fragments_DB:list,
#                core_smiles:str='C', # TODO - TBD
               threshold:float=0.5,
               overall_max_n:int=None,
               frag_max_n:int=None,
               scaffold_top_n:int=None,
               cand_max_n__rplc:int=None,
               _merge_structure_top_n_:int=100,
               output_dir:str='./output',
               low_mem:bool=False,
               tqdm_quiet:bool=False,
               fragments:list=[],
               replace_scaffold_files:list=[],
              ):
    if type(output_dir)==str:
        os.makedirs(output_dir,exist_ok=True)
    
    # TODO - move calling DB to the section of replace scaffold
    if not fragments_DB:
        print("Calling fragment DB..")
        if not low_mem:
            fragments_DB = utils.call_frag_db()[2]
        else:
            fragments_DB = utils._call_frag_db_smi_()[1]
    target_mol, target_smiles = utils.init_mol(target_smiles)
    
    result_features = [
        'Fragment_no','Original scaffold',
        'Replaced scaffold','Final structure','Standardized final structure',
        'Tanimoto Similarity','Electron shape Similarity',
        'CATS2D dist', 'QED', 'SAscore', 'logP']
    fp = open(os.path.join(output_dir,'overall_result.txt'), 'wb')
    fp.write('\t'.join(result_features).encode())
    fp.write('\n'.encode())
    # Fragment info
    if fragments:
        frags = [Chem.MolFromSmiles(i) for i in fragments]
    else:
        frags = sg.get_all_murcko_fragments(target_mol, break_fused_rings=False)
    frags_n = len(frags)
    saved_results = {target_smiles:1}
    
    frag_info = pd.DataFrame([[Chem.MolToSmiles(i) for i in frags]],index=['SMILES']).T
    frag_info.index.name = 'Fragment_no'
    frag_info.to_csv(os.path.join(output_dir,'fragment_info.tsv'),sep='\t')
    print(f"Fragments found\t: {len(frags)}")
    
    # Getting replace_mol_list
    print("Finding scaffolds for replacement...")
#     scaffold_len_n_l = []
    for frag_no, pattern_mol in tqdm.tqdm(enumerate(frags), desc='Fragments'):
        # Pre-defined replace_scaffold_list for given fragment
        _scf_f_n_ = os.path.join(output_dir,f"replace_scaffold_list.fragment_{frag_no:04d}.tsv")
        predefined=False
        # pre-defined scaffold list
        if fragments and replace_scaffold_files:
            curr_frag_f = replace_scaffold_files[frag_no]
            print(f"Finding scaffolds for fragment {Chem.MolToSmiles(pattern_mol)} at {curr_frag_f}")
            if os.path.isfile(curr_frag_f):
                try:
                    # pre-defined one
                    _scf_ser = pd.read_csv(curr_frag_f,sep='\t',index_col=0)
                    # In case that replace_scaffold_file is SMILES-only without data (without header)
                    if len(_scf_ser.columns) == 0:
                        with open(curr_frag_f,'rb') as f:
                            _scf_ser = f.read().decode().splitlines()
                        _scf_ser = pd.Series(list(range(len(_scf_ser),0,-1)),index=_scf_ser)
                        _scf_ser.name='Order'
                    predefined=True
                    _scf_ser.to_csv(_scf_f_n_,sep='\t')
                except:
                    predefined=False
        if not predefined:
            print(f"Finding alternatives for\t\t{Chem.MolToSmiles(pattern_mol)}")
            replace_scaffold_ser = search_similar_scaffolds(
                original_scaffold=pattern_mol,
                fragments_DB=fragments_DB,
                low_mem=low_mem,
                tqdm_quiet=tqdm_quiet,
            )
            replace_scaffold_ser.name='Tanimoto Similarity'
            replace_scaffold_ser.to_csv(_scf_f_n_,sep='\t')
            
    overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_ = utils._scaffold_no_reassign_(
        overall_max_n=overall_max_n,
        frag_max_n=frag_max_n,
        scaffold_top_n=scaffold_top_n,
        cand_max_n__rplc=cand_max_n__rplc,
        _merge_structure_top_n_=_merge_structure_top_n_,
        frag_n=len(frags),
    )
    # TODO remove here
    print(overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_)
    #
    # Finding candidates
    _overall_cand_cnt_ = 0
    
#     # TODO-remove here
#     return frags,overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_
#     #
    
    for frag_no, pattern_mol in enumerate(frags):
        gc.collect()
        original_scaffold = Chem.MolToSmiles(pattern_mol)
        pattern_mol, original_scaffold = utils.init_mol(original_scaffold)
        _frag_start = datetime.datetime.now()
        _frag_cand_cnt_ = 0
        # Fragment IO
        f_frag = open(os.path.join(output_dir,f'fragment_result_{frag_no:04d}.txt'), 'wb')
        f_frag.write('\t'.join(result_features).encode())
        f_frag.write('\n'.encode())
        frag_saved_results = {target_smiles:1}
        
        replace_scaffold_df = pd.read_csv(
            os.path.join(output_dir,f"replace_scaffold_list.fragment_{frag_no:04d}.tsv"),
            sep='\t',index_col=0)
        replace_scaffold_df = replace_scaffold_df.sort_values(
            by=replace_scaffold_df.columns[0],ascending=False)
        
        # TODO - application of scaffold_top_n when replace_scaffold_files is imposed
        if fragments and replace_scaffold_files:
            scaffold_top_n__frg = scaffold_top_n
#             scaffold_top_n__frg = len(replace_scaffold_df.index) # all th SMILES of replace_scaffold_list will be proceeded at reach to overall_max_n
        else:
            scaffold_top_n__frg = scaffold_top_n
        
        # TODO - limit iteration of screening to satisfy overall_max_n
        replc_scf_it_obj = tqdm.tqdm(enumerate(list(replace_scaffold_df.index)[:scaffold_top_n__frg]), desc='Scaffold')
        for replc_scf_cnt, replace_scaffold in replc_scf_it_obj:
            gc.collect()
            replace_mol, replace_scaffold = utils.init_mol(replace_scaffold)
            try:
                final_candidates = replace_molecule(
                    target_mol, pattern_mol, replace_mol,
                    max_n=cand_max_n__rplc,
                    _merge_structure_top_n_=_merge_structure_top_n_)
            except Exception as e:
                continue
            for cand_cnt, [_, each_candidate] in enumerate(final_candidates):
                if each_candidate in frag_saved_results:
                    continue
                frag_saved_results[each_candidate] = 1
                mol = Chem.MolFromSmiles(each_candidate)
                # Tanimoto similarity cutoff for selection
                try: # Possible error in calculation of Tanimoto similarity
                    tanimoto_sim = utils.calc_tanimoto_sim(target_mol, mol)
                except:
                    continue
                # Thresholds
                if tanimoto_sim >= threshold:
                    # calculation of electron similarity and subscores
                    try: # Possible error in calculation of electron shape
                        electron_shape_sim = utils.calc_electron_shape(target_smiles, each_candidate)
                    except:
                        electron_shape_sim = float('nan')
                    sa_score, qed_score, _mw_, logp_score = utils._get_molecular_prop_(mol)
                    
                    cats_des_dist = utils.calculate_cats_des(target_mol, mol)
                    # Standardization of smiles structure
                    try:
                        _fin_structure, valid_info = utils.get_standard_smiles(each_candidate)
                    except Exception as e:
                        print(f'Failed to find the standard SMILES of {each_candidate}',e)
                        _fin_structure = ''
                    f_frag.write('\t'.join([str(i) for i in [
                        frag_no, original_scaffold,
                        replace_scaffold, each_candidate, _fin_structure,
                        tanimoto_sim, electron_shape_sim,
                        cats_des_dist, qed_score, sa_score, logp_score,
                    ]]).encode())
                    f_frag.write('\n'.encode())
                    _frag_cand_cnt_ += 1
                    if each_candidate not in saved_results:
                        fp.write('\t'.join([str(i) for i in [
                            frag_no, original_scaffold,
                            replace_scaffold, each_candidate, _fin_structure,
                            tanimoto_sim, electron_shape_sim,
                            cats_des_dist, qed_score, sa_score, logp_score,
                        ]]).encode())
                        saved_results[each_candidate] = 1
                        fp.write('\n'.encode())
                        _overall_cand_cnt_ += 1
            # Early stopping
            # 2-fold more iterations in case that number of candidates is less than max_n
#             if replc_scf_cnt%10==9:
#                 print(replc_scf_cnt,scaffold_top_n__frg)
            if _frag_cand_cnt_ > frag_max_n or replc_scf_cnt > scaffold_top_n__frg:
                print(_frag_cand_cnt_, "candidates have found")
                replc_scf_it_obj.close()
                print("Time cost for fragment ",frag_no)
                print(datetime.datetime.now() - _frag_start)
                f_frag.close()
                break
        print("Time cost for fragment ",frag_no)
        print(datetime.datetime.now() - _frag_start)
        f_frag.close()
    fp.close()
    result_df = pd.read_csv(os.path.join(output_dir,'overall_result.txt'),sep='\t')
    print(f"Found hopped structures\t: {result_df.shape[0]}")
    return result_df


# Main function
def main():
    start = datetime.datetime.now()
    print("Started at \t",start)
    
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    warnings.filterwarnings('ignore')
    RDLogger.DisableLog('rdApp.*')
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    mylogger = logging.getLogger("log")
    mylogger.setLevel(logging.INFO)

    stream_hander = logging.StreamHandler()
    stream_hander.setFormatter(formatter)
    mylogger.addHandler(stream_hander)
    
    file_handler = logging.FileHandler('log.log')
    file_handler.setFormatter(formatter)
    mylogger.addHandler(file_handler)
    
    parser = utils.argument_parser()
    
    options = parser.parse_args()    
    target_smiles = options.input_smiles
#     core_smiles = options.core_smiles
    output_dir = os.path.abspath(options.output_dir)
    
    os.makedirs(output_dir,exist_ok=True)
    
    if options.low_mem:
        _, fragments_DB = utils._call_frag_db_smi_()
    else:
        _, _, fragments_DB = utils.call_frag_db()

    result_df = chembounce(
        target_smiles=target_smiles,
        fragments_DB=fragments_DB,
#         core_smiles=options.core_smiles, # TODO-further dev
        overall_max_n=options.overall_max_n,
        frag_max_n=options.frag_max_n,
        scaffold_top_n=options.scaffold_top_n,
        cand_max_n__rplc=options.cand_max_n__rplc,
        _merge_structure_top_n_=100,
        threshold=options.tanimoto_threshold,
        output_dir=output_dir,
        low_mem=options.low_mem,
        fragments=options.fragments,
        replace_scaffold_files=options.replace_scaffold_files,
    )
    end = datetime.datetime.now()
    print("Finished at\t",end)
    cost = end-start
    print("Time cost: \t",end-start)

if __name__ == '__main__':
    main()
