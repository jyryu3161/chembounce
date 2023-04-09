import os
import glob
import warnings
import time
import logging
import numpy as np
import math
import time
import argparse

from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem.Descriptors import qed
from tqdm import tqdm

import numpy as np
from keras.models import model_from_json 
from numpy import dot
from numpy.linalg import norm
from scipy.spatial import distance
import re
from rdkit.Chem.Descriptors import ExactMolWt
import utils
import scaffoldgraph as sg

def calculate_fingerprints(smiles_list, r=2): 
    ecfp_features = []
    smi_list = []
    for i in range(len(smiles_list)):
        #try:
        smi = smiles_list[i]
        mol = Chem.MolFromSmiles(smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=r, nBits=2048)
        bits = fp.ToBitString()
        merged_feature = []
        feature = []
        for f in bits:
            feature.append(int(f))
            merged_feature.append(int(f))

        smi_list.append(smi)
        ecfp_features.append(merged_feature)
        #except:
        #    pass
    return np.asarray(ecfp_features), smi_list

def check_feasibility(input_smiles_list):
    compound_model = './data/compound_model.json'
    compound_model_weight = './data/compound_model.h5'
    
    json_file = open(compound_model, "r")
    compound_loaded_model_json = json_file.read() 
    json_file.close()
    compound_loaded_encoder_model = model_from_json(compound_loaded_model_json)
    compound_loaded_encoder_model.load_weights(compound_model_weight)

    compound_features, smiles_candidates = calculate_fingerprints(input_smiles_list)
    compound_results = compound_loaded_encoder_model.predict(compound_features)
    compound_results = np.where(compound_results > 0.5, 1, 0)
    
    feasible_compound_smiles = []
    for i in tqdm(range(len(compound_features))):
        original_feat = compound_features[i]
        pred_feat = compound_results[i]
        smiles = smiles_candidates[i]
        dist = np.linalg.norm(original_feat-pred_feat)    
        cos_sim = dot(original_feat, pred_feat)/(norm(original_feat)*norm(pred_feat))
        if cos_sim > 0.9:
            feasible_compound_smiles.append(smiles)
    
    return feasible_compound_smiles
