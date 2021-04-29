import logging
import utils
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import AllChem
import tqdm

import numpy as np

import scaffoldgraph as sg
import itertools
import copy
import os
import warnings

import bioisostere
import feasibility_utils

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
    structure_candidates = structure_candidates[0:3]
    return structure_candidates

def replace_molecule(target_mol, pattern_mol, replace_mol, top_n):
    start_replace_scaffold_list = [replace_mol]
    candidate_structures = remove_substructures(target_mol, pattern_mol)
    
    final_candidates = []
    saved_structures = {}
    candidate_sets = itertools.permutations(candidate_structures)
    
    for each_set in candidate_sets:
        candidate_set = list(each_set)
        for each_structure in candidate_set:
            frag_mol = Chem.MolFromSmiles(each_structure)
            tmp_replace_scaffold_list = []

            for each_replace_mol in start_replace_scaffold_list:
                structure_candidates = merge_structures(target_mol, each_replace_mol, frag_mol, top_n)

                for each_candidate in structure_candidates:
                    each_score = each_candidate[0]
                    each_str = each_candidate[1]
                    if each_str not in saved_structures:
                        saved_structures[each_str] = True
                        final_candidates.append([each_score, each_str])
                        each_str_mol, each_str = utils.init_mol(each_str)
                        tmp_replace_scaffold_list.append(each_str_mol)
            tmp_replace_scaffold_list = list(set(tmp_replace_scaffold_list))
            start_replace_scaffold_list = tmp_replace_scaffold_list
            start_replace_scaffold_list = start_replace_scaffold_list[0:3]
    
    final_candidates.sort(reverse=True)
    final_candidates = final_candidates[0:top_n]
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

def search_similar_scaffolds(original_scaffold, fragments_DB, scaffold_top_n):
    scaffold_scores = []
    threshold = 0.7
    for each_frag_candidate in tqdm.tqdm(fragments_DB):
        sim = utils.calc_tanimoto_sim(original_scaffold, each_frag_candidate)
        if sim > threshold:
            if sim != 1.0:
                scaffold_scores.append([sim, Chem.MolToSmiles(each_frag_candidate)])
    
    replace_scaffold_list = []
    scaffold_scores.sort(reverse=True)
    scaffold_scores = scaffold_scores[0:scaffold_top_n]
    for each_data in scaffold_scores:
        replace_scaffold_list.append(each_data[1])
    return replace_scaffold_list

def main():
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    warnings.filterwarnings('ignore')
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    mylogger = logging.getLogger("log")
    mylogger.setLevel(logging.INFO)

    stream_hander = logging.StreamHandler()
    stream_hander.setFormatter(formatter)
    mylogger.addHandler(stream_hander)
    
    file_handler = logging.FileHandler('log.log')
    file_handler.setFormatter(formatter)
    mylogger.addHandler(file_handler)
    
    
    target_smiles = 'C1=CC(=CC=C1C(=O)OC2=CC3=C(C=C2)C=C(C=C3)C(=N)N)N=C(N)N'
    core_smiles = 'C1=CC2=C(C=C1)C=CC=C2'
    output_dir = './output/'
    
    scaffold_top_n = 1000
    final_top_n = 10
    bioisostere_flag = True
    
    try:
        os.mkdir(output_dir)
    except:
        pass
    
    fragment_file = './data/Scaffolds_processed.txt'
    fragments_DB = read_fragments(target_smiles, fragment_file)
    
    
    target_mol, target_smiles = utils.init_mol(target_smiles)
    
    fp = open(output_dir+'/result.txt', 'w')
    fp.write('%s\t%s\t%s\t%s\t%s\n'%('Original scaffold', 'Replacement scaffold', 'Final product', 'TC Similarity', 'Electron Similarity'))
    frags = sg.get_all_murcko_fragments(target_mol, break_fused_rings=False)
    
    for pattern_mol in tqdm.tqdm(frags):
        original_scaffold = Chem.MolToSmiles(pattern_mol)
        pattern_mol, original_scaffold = utils.init_mol(original_scaffold)
        
        replace_scaffold_list = search_similar_scaffolds(pattern_mol, fragments_DB, scaffold_top_n)
        
        for replace_scaffold in tqdm.tqdm(replace_scaffold_list):
            replace_mol, replace_scaffold = utils.init_mol(replace_scaffold)
            try:
                final_candidates = replace_molecule(target_mol, pattern_mol, replace_mol, final_top_n)
            except:
                continue
                
            for each_score, each_candidate in final_candidates:
                try:
                    sim = utils.calc_electron_shape(target_smiles, each_candidate)
                    if sim >= 0.7:
                        fp.write('%s\t%s\t%s\t%s\t%s\n'%(original_scaffold, replace_scaffold, each_candidate, each_score, sim))
                except:
                    continue
    fp.close()
    
    candidates = []
    candidate_smiles_list = []
    with open(output_dir+'/result.txt', 'r') as fp:
        fp.readline()
        for line in fp:
            sptlist = line.strip().split('\t')
            smiles = sptlist[2].strip()
            candidates.append(smiles)
            candidate_smiles_list.append(smiles)
    candidates = list(set(candidates))
    fp2 = open(output_dir+'/result_bioisostere.txt', 'w')
    fp2.write('%s\t%s\n'%('Original Smiles', 'Bioisostere'))
    for smiles in tqdm.tqdm(candidates):
        try:
            total_smiles_list = bioisostere.run_transformation(smiles, 1)
        except:
            continue
        for each_smiles in total_smiles_list:
            fp2.write('%s\t%s\n'%(smiles, each_smiles))
            candidate_smiles_list.append(each_smiles)
    fp2.close()
    
    candidate_smiles_list = list(set(candidate_smiles_list))
    feasible_smiles_list = feasibility_utils.check_feasibility(candidate_smiles_list)
    fp = open(output_dir+'/result_final_structures.txt', 'w')
    fp.write('%s\n'%('Smiles'))
    for smiles in feasible_smiles_list:
        fp.write('%s\n'%(smiles))
    fp.close()
    return

if __name__ == '__main__':
    main()
    