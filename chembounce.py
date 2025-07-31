#!/usr/bin/env python3
"""
ChemBounce main functions and subfuctions
"""

import os
import sys
import tqdm
import gc
import numpy as np
import pandas as pd
import itertools
import copy
import json
import warnings
import pickle
import copy
import datetime
import logging
sys.path.append(os.path.split(os.path.abspath(__file__))[0])

from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import AllChem

import scaffoldgraph_fragmenter as sg
import utils
import cli
import cost_estimation


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
                             low_mem:bool=False, tqdm_quiet:bool=False,
                             scaffold_top_n:int=None, threshold:float=0.3, # In case to limit results
                            ):
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
        # Similarity threshold : original_scaffold vs frag_candidate
        sim = utils.calc_tanimoto_sim(original_scaffold, candidate_mol)
        if sim > 0.1:
            if sim != 1.0:
                scaffold_scores[candidate_smiles] = sim
    scaffold_scores = pd.Series(scaffold_scores)
    scaffold_scores = scaffold_scores.sort_values(ascending=False)
    # Limit
    if threshold:
        scaffold_scores = scaffold_scores.loc[scaffold_scores>threshold]
    if scaffold_top_n and type(scaffold_top_n)==int:
        if scaffold_scores.shape[0]>scaffold_top_n:
            scaffold_scores = scaffold_scores.iloc[:scaffold_top_n]
        
    return scaffold_scores

    
def get_frags_cands(target_smiles:str,
                    target_mol,
                    fragments_DB:list=[],
                    overall_max_n:int=None,
                    frag_max_n:int=None,
                    scaffold_top_n:int=None,
                    cand_max_n__rplc:int=None,
                    _merge_structure_top_n_:int=100,
                    murcko_frag_itr_rnd:int=1,
                    _search_scf_thr_:float=0.3,
                    fragments:list=[],
                    replace_scaffold_files:list=[],
                    output_dir:str='./output',
                    low_mem:bool=False,
                    tqdm_quiet:bool=False,
                   ):
    # Fragment info
    if fragments:
        frags = [Chem.MolFromSmiles(i) for i in fragments]
    else:
        # Modified function
        frags = sg.get_all_murcko_fragments(
            target_mol,
            break_fused_rings=False,
            iteration_round=murcko_frag_itr_rnd)
    
    frag_info = pd.DataFrame([[Chem.MolToSmiles(i) for i in frags]],index=['SMILES']).T
    frag_info.index.name = 'Fragment_no'
    frag_info.to_csv(os.path.join(output_dir,'fragment_info.tsv'),sep='\t')
    print(f"Fragments found\t: {len(frags)}")
    
    # Organize how much tests and candidates to find
    overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_ = utils._scaffold_no_reassign_(
        overall_max_n=overall_max_n,
        frag_max_n=frag_max_n,
        scaffold_top_n=scaffold_top_n,
        cand_max_n__rplc=cand_max_n__rplc,
        _merge_structure_top_n_=_merge_structure_top_n_,
        frag_n=len(frags),
    )
    print(
        "<Applied iteration limits>\n\t",
        f"overall_max_n:\t{overall_max_n}\n\t",
        f"frag_max_n:\t{frag_max_n}\n\t",
        f"scaffold_top_n:\t{scaffold_top_n}\n\t",
        f"cand_max_n__rplc:\t{cand_max_n__rplc}\n\t",
        f"_merge_structure_top_n_:\t{_merge_structure_top_n_}\n\t",)
    _search_scf_max_n_ = None
    if low_mem:
        _search_scf_max_n_ = scaffold_top_n
    # Getting replace_mol_list
    print("Finding scaffolds for replacement...")
    for frag_no, pattern_mol in tqdm.tqdm(enumerate(frags), desc='Fragments'):
        # Pre-defined replace_scaffold_list for given fragment
        _scf_f_n_ = os.path.join(output_dir,f"replace_scaffold_list.fragment_{frag_no:04d}.tsv")
        predefined=False
        # pre-defined scaffold list
        if fragments and replace_scaffold_files:
            if type(replace_scaffold_files)==str:
                curr_frag_f = replace_scaffold_files
            elif len(replace_scaffold_files)==1:
                curr_frag_f = replace_scaffold_files[0]
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
            else:
                print(f"File not found: {curr_frag_f}")
        if not predefined:
            print(f"Finding alternatives for\t\t{Chem.MolToSmiles(pattern_mol)}")
            if not fragments_DB:
                print("Calling fragment DB..")
                if not low_mem:
                    fragments_DB = utils.call_frag_db()[2]
                else:
                    fragments_DB = utils._call_frag_db_smi_()[1]
            replace_scaffold_ser = search_similar_scaffolds(
                original_scaffold=pattern_mol,
                fragments_DB=fragments_DB,
                low_mem=low_mem,
                tqdm_quiet=tqdm_quiet,
                scaffold_top_n=_search_scf_max_n_,
                threshold=_search_scf_thr_,
            )
            replace_scaffold_ser.name='Tanimoto Similarity'
            replace_scaffold_ser.to_csv(_scf_f_n_,sep='\t')
    
    return frags, overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_
    
    
def make_scaffold_hopping(target_smiles:str,
                          target_mol,
                          frags:list,
                          core_mol=None,
                          tanimoto_threshold:float=0.5,
                          overall_max_n:int=None,
                          frag_max_n:int=None,
                          scaffold_top_n:int=None,
                          cand_max_n__rplc:int=None,
                          _merge_structure_top_n_:int=100,
                          output_dir:str='./output',
                          candidate_thresholds:dict=dict(),
                         ):
    result_features = [
        'Fragment_no','Original scaffold',
        'Replaced scaffold','Final structure','Standardized final structure',
        'Tanimoto Similarity','Electron shape Similarity',
        'CATS2D dist', 'QED', 'SAscore', 'logP', 'MW', 'H_Donors', 'H_Acceptors']
    fp = open(os.path.join(output_dir,'overall_result.txt'), 'wb')
    fp.write('\t'.join(result_features).encode())
    fp.write('\n'.encode())
    # Fragment info
    saved_results = {target_smiles:1}
    # Finding candidates
    _overall_cand_cnt_ = 0
    for frag_no, pattern_mol in enumerate(frags):
        print(f"Screening started for fragment_{frag_no:04d}: {Chem.MolToSmiles(pattern_mol)}")
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
        
        # TODO - limit iteration of screening to satisfy overall_max_n
        replc_scf_it_obj = tqdm.tqdm(enumerate(list(replace_scaffold_df.index)[:scaffold_top_n]), desc='Scaffold')
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
                cand_mol = Chem.MolFromSmiles(each_candidate)
                # Tanimoto similarity cutoff for selection
                try: # Possible error in calculation of Tanimoto similarity
                    tanimoto_sim = utils.calc_tanimoto_sim(target_mol, cand_mol)
                except:
                    continue
                # Thresholds
                if tanimoto_sim >= tanimoto_threshold:
                    if core_mol:
                        if not target_mol.HasSubstructMatch(core_mol):
                            continue
                    # calculation of electron similarity and subscores
                    try: # Possible error in calculation of electron shape
                        electron_shape_sim = utils.calc_electron_shape(target_smiles, each_candidate)
                    except:
                        electron_shape_sim = float('nan')
                    props = utils._get_molecular_prop_extd_(mol=cand_mol)
                    # Property determination
                    passed, determ_results = utils.determ_props(
                        props=props,thresholds=candidate_thresholds,
                        w_all_failed_reason=False,ignore_null_val=False)
                    if not passed:
                        continue
                    cats_des_dist = utils.calculate_cats_des(target_mol, cand_mol)
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
                        cats_des_dist, props['qed'],props['sa'],props['logp'],
                        props['mw'],props['n_hdonor'],props['n_hacceptor'],
                    ]]).encode())
                    f_frag.write('\n'.encode())
                    _frag_cand_cnt_ += 1
                    if each_candidate not in saved_results:
                        fp.write('\t'.join([str(i) for i in [
                            frag_no, original_scaffold,
                            replace_scaffold, each_candidate, _fin_structure,
                            tanimoto_sim, electron_shape_sim,
                            cats_des_dist, props['qed'],props['sa'],props['logp'],
                            props['mw'],props['n_hdonor'],props['n_hacceptor'],
                        ]]).encode())
                        saved_results[each_candidate] = 1
                        fp.write('\n'.encode())
                        _overall_cand_cnt_ += 1
            # Early stopping
            # 2-fold more iterations in case that number of candidates is less than max_n
            if _frag_cand_cnt_ > frag_max_n or replc_scf_cnt > scaffold_top_n:
                print(_frag_cand_cnt_, "candidates have found")
                replc_scf_it_obj.close()
                break
        print(_frag_cand_cnt_, "candidates have found")
        print("Time cost for fragment ",frag_no)
        print(datetime.datetime.now() - _frag_start)
        f_frag.close()
    fp.close()
    return os.path.join(output_dir,'overall_result.txt')


# Main function
@cost_estimation.EstDecoMem
def chembounce(target_smiles:str,
               fragments_DB:list=[],
               core_smiles:str=None,
               tanimoto_threshold:float=0.5,
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
               candidate_thresholds:dict=dict(),
               murcko_frag_itr_rnd:int=1, # Iteration round limiation
               _search_scf_thr_:float=0.3, # Tanimoto similarity threshold for search_scaffold function
              ):
    # IO
    if type(output_dir)==str:
        os.makedirs(output_dir,exist_ok=True)
    try:
        target_mol, target_smiles = utils.init_mol(target_smiles)
    except Exception as e:
        raise utils.SmilesInputError(value=target_smiles)
    if core_smiles:
        try:
            core_mol, core_smiles = utils.init_mol(core_smiles)
        except:
            core_mol=None
            core_smiles=None
    else:
        core_mol=None
        core_smiles=None
    if core_mol:
        if not(target_mol.HasSubstructMatch(core_mol)):
            raise utils.SmilesInputError(value=core_smiles, message='Core SMILES is not contained in the input SMILES')
    
    frags, overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_ = get_frags_cands(
        target_smiles=target_smiles,
        target_mol=target_mol,
        fragments_DB=fragments_DB,
        overall_max_n=overall_max_n,
        frag_max_n=frag_max_n,
        scaffold_top_n=scaffold_top_n,
        cand_max_n__rplc=cand_max_n__rplc,
        _merge_structure_top_n_=_merge_structure_top_n_,
        output_dir=output_dir,
        low_mem=low_mem,
        tqdm_quiet=tqdm_quiet,
        fragments=fragments,
        replace_scaffold_files=replace_scaffold_files,
        murcko_frag_itr_rnd=murcko_frag_itr_rnd,
        _search_scf_thr_=_search_scf_thr_,
    )
    
    o_f = make_scaffold_hopping(
        target_smiles=target_smiles,
        target_mol=target_mol,
        frags=frags,
        core_mol=core_mol,
        tanimoto_threshold=tanimoto_threshold,
        overall_max_n=overall_max_n,
        frag_max_n=frag_max_n,
        scaffold_top_n=scaffold_top_n,
        cand_max_n__rplc=cand_max_n__rplc,
        _merge_structure_top_n_=_merge_structure_top_n_,
        output_dir=output_dir,
        candidate_thresholds=candidate_thresholds,
    )
    result_df = pd.read_csv(o_f,sep='\t')
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
    
    parser = cli.argument_parser()
    options = parser.parse_args()
    target_smiles = options.input_smiles
    output_dir = os.path.abspath(options.output_dir)
    
    os.makedirs(output_dir,exist_ok=True)
    
    # candidate threshold
    candidate_thresholds = utils._default_thrs_(lipinski=not options.wo_lipinski)
    user_thr = {
        'qed':(options.qed_min,options.qed_max),
        'sa':(options.sa_min,options.sa_max),
        'logp':(options.logp_min,options.logp_max),
        'mw':(options.mw_min,options.mw_max),
        'n_hdonor':(options.h_donor_min,options.h_donor_max),
        'n_hacceptor':(options.h_acceptor_min,options.h_acceptor_max),
    }
    for metric, (min_val,max_val) in user_thr.items():
        if type(min_val) in [float,int] or type(max_val) in [float,int]:
            if metric in candidate_thresholds:
                o_min,o_max = candidate_thresholds[metric]
                if type(min_val) in [float,int]:
                    o_min = min_val
                if type(max_val) in [float,int]:
                    o_max = max_val
                candidate_thresholds[metric] = (o_min,o_max)    
            else:
                candidate_thresholds[metric] = (min_val,max_val)
    print(f"Applied candidate thresholds (Lipinski\'s rule of five : {not options.wo_lipinski}):")
    print('\n'.join([
        f"\t{metric} :\tMin.:{min_val}\tMax.:{max_val}" for metric, (min_val,max_val) in candidate_thresholds.items()]))
    print(f'Tanimoto Similarity :\t{options.tanimoto_threshold}')
    
    result_df, resource_cost = chembounce(
        target_smiles=target_smiles,
#         fragments_DB=fragments_DB, # Requirement has been changed
        overall_max_n=options.overall_max_n,
        core_smiles=options.core_smiles,
        frag_max_n=options.frag_max_n,
        scaffold_top_n=options.scaffold_top_n,
        cand_max_n__rplc=options.cand_max_n__rplc,
        _merge_structure_top_n_=100,
        tanimoto_threshold=options.tanimoto_threshold,
        output_dir=output_dir,
        low_mem=options.low_mem,
        fragments=options.fragments,
        replace_scaffold_files=options.replace_scaffold_files,
        candidate_thresholds=candidate_thresholds,
        murcko_frag_itr_rnd=1,
        _search_scf_thr_=0.3,
#     )
        #murcko_frag_itr_rnd:int=options.murcko_frag_itr_rnd,
        #_search_scf_thr_:float=options.search_scf_thr,
    )
    # if options.estimate_cost:
    with open(os.path.join(output_dir,'resource_cost.json'),'wb') as f:
        f.write(json.dumps(resource_cost).encode())
    print(f"""#### Cost estimatation ####
            Elapsed time:              {resource_cost['elapsed_time']}
            Total CPU time:            {resource_cost['total_cpu_time']}
            CPU usage percent:         {resource_cost['cpu_usage_percent']}
            Maximal memory usage (MB): {resource_cost['max_memory_mb']}
            Average memory usage (MB): {resource_cost['avg_memory_mb']}
        """)
    end = datetime.datetime.now()
    print("Finished at\t",end)
    cost = end-start
    print("Time cost: \t",end-start)

    
if __name__ == '__main__':
    main()
