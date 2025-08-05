#!/usr/bin/env python3
"""
Generate pre-computed fingerprints for scaffold database
This significantly speeds up similarity search
"""

import os
import sys
import numpy as np
import pickle
import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.SimDivFilters import MaxMinPicker
import argparse
import time

def sample_by_diversity(fps, percentages, seed=42):
    """
    Selects diverse subsets of fingerprints using the MaxMinPicker algorithm.
    Ensures that smaller percentage subsets are contained within larger ones.

    Args:
        fps: A list of RDKit fingerprint objects.
        percentages: A list of percentages (e.g., [10, 25, 50]) for the subsets.
        seed: Random seed for reproducibility.

    Returns:
        A dictionary mapping each percentage to a list of selected indices.
    """
    if not fps:
        return {p: [] for p in percentages}

    num_fps = len(fps)
    picker = MaxMinPicker()
    
    def dist_func(i, j, fps=fps):
        return 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])

    results = {}
    # Sort percentages to pick for the largest first, ensuring subsets are nested
    sorted_percentages = sorted(percentages, reverse=True)
    
    # Find the largest percentage that results in at least one pick
    p_max = 0
    num_to_pick_max = 0
    for p in sorted_percentages:
        num_to_pick = int(num_fps * p / 100)
        if num_to_pick > 0:
            p_max = p
            num_to_pick_max = num_to_pick
            break
    
    if num_to_pick_max > 0:
        # Pick indices for the largest valid percentage
        last_picked_indices = list(picker.LazyPick(dist_func, num_fps, num_to_pick_max, seed=seed))
        results[p_max] = last_picked_indices
        
        # Create subsets for smaller percentages from the largest set
        for p in sorted_percentages:
            if p == p_max:
                continue
            num_to_pick = int(num_fps * p / 100)
            results[p] = last_picked_indices[:num_to_pick]

    # Ensure all requested percentages are in the final dict
    for p in percentages:
        if p not in results:
            results[p] = []
            
    return results

def generate_fingerprints(input_file, output_file, fp_size=1024, radius=2, batch_size=10000):
    """
    Generate Morgan fingerprints for all molecules in the input file,
    and create diversity-based subsets.
    
    Args:
        input_file: Path to scaffold SMILES file
        output_file: Path to save fingerprint data. This will be used as a base name.
        fp_size: Size of fingerprint bit vector
        radius: Morgan fingerprint radius
        batch_size: Number of molecules to process at once
    """
    
    print(f"Reading scaffolds from {input_file}")
    
    # Count total lines for progress bar
    with open(input_file, 'r') as f:
        total_lines = sum(1 for _ in f)
    
    print(f"Total scaffolds: {total_lines}")
    
    # Process in batches to manage memory
    all_fps = []
    all_smiles = []
    valid_indices = []
    
    with open(input_file, 'r') as f:
        batch_smiles = []
        batch_fps = []
        
        for i, line in enumerate(tqdm.tqdm(f, total=total_lines, desc="Generating fingerprints")):
            smiles = line.strip()
            if not smiles:
                continue
                
            batch_smiles.append(smiles)
            
            # Generate fingerprint
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=fp_size)
                    batch_fps.append(fp)
                    valid_indices.append(i)
                else:
                    # Skip invalid molecules
                    batch_smiles.pop()
            except:
                # Skip molecules that cause errors
                batch_smiles.pop()
                continue
            
            # Save batch when full
            if len(batch_smiles) >= batch_size:
                all_smiles.extend(batch_smiles)
                all_fps.extend(batch_fps)
                batch_smiles = []
                batch_fps = []
        
        # Save remaining batch
        if batch_smiles:
            all_smiles.extend(batch_smiles)
            all_fps.extend(batch_fps)
    
    print(f"Generated {len(all_fps)} valid fingerprints out of {total_lines} scaffolds")
    
    if not all_fps:
        print("No valid fingerprints were generated. Exiting.")
        return 0

    # Convert fingerprints to numpy array for faster operations
    print("Converting to numpy array...")
    fp_array = np.zeros((len(all_fps), fp_size), dtype=np.uint8)
    
    for i, fp in enumerate(tqdm.tqdm(all_fps, desc="Converting to array")):
        arr = np.zeros((fp_size,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fp_array[i] = arr
    
    # Save original full data (100%)
    print(f"\nSaving full dataset (100%) to {output_file}")
    data = {
        'fingerprints': fp_array,
        'smiles': np.array(all_smiles, dtype=object),
        'valid_indices': np.array(valid_indices, dtype=int),
        'fp_size': fp_size,
        'radius': radius,
        'version': '1.0'
    }
    np.savez_compressed(output_file, **data)
    
    return len(all_fps)

def verify_fingerprints(fp_file, sample_size=100):
    """Verify the generated fingerprint file"""
    print(f"\nVerifying {fp_file}")
    
    data = np.load(fp_file, allow_pickle=True)
    fps = data['fingerprints']
    smiles = data['smiles']
    
    print(f"Total fingerprints: {len(fps)}")
    print(f"Fingerprint shape: {fps.shape}")
    
    # Test a few random similarities
    if len(fps) > sample_size:
        indices = np.random.choice(len(fps), sample_size, replace=False)
        
        print(f"\nTesting {sample_size} random similarity calculations...")
        start_time = time.time()
        
        for i in range(10):
            idx1, idx2 = np.random.choice(indices, 2, replace=False)
            
            # Calculate Tanimoto similarity using numpy operations
            fp1 = fps[idx1]
            fp2 = fps[idx2]
            
            intersection = np.sum(fp1 & fp2)
            union = np.sum(fp1 | fp2)
            similarity = intersection / union if union > 0 else 0.0
            
            # Verify with RDKit
            mol1 = Chem.MolFromSmiles(smiles[idx1])
            mol2 = Chem.MolFromSmiles(smiles[idx2])
            fp1_rdkit = AllChem.GetMorganFingerprintAsBitVect(mol1, data['radius'].item(), nBits=data['fp_size'].item())
            fp2_rdkit = AllChem.GetMorganFingerprintAsBitVect(mol2, data['radius'].item(), nBits=data['fp_size'].item())
            similarity_rdkit = DataStructs.TanimotoSimilarity(fp1_rdkit, fp2_rdkit)
            
            print(f"Pair {i+1}: {similarity:.4f} (numpy) vs {similarity_rdkit:.4f} (RDKit)")
        
        elapsed = time.time() - start_time
        print(f"\nTime for {10} similarity calculations: {elapsed:.4f} seconds")

def main():
    parser = argparse.ArgumentParser(description='Generate fingerprints for scaffold database')
    parser.add_argument('-i', '--input', default='data/Scaffolds_processed.txt',
                        help='Input SMILES file')
    parser.add_argument('-o', '--output', default='data/scaffold_fingerprints.npz',
                        help='Output fingerprint file')
    parser.add_argument('--fp-size', type=int, default=2048,
                        help='Fingerprint size (default: 2048)')
    parser.add_argument('--radius', type=int, default=2,
                        help='Morgan fingerprint radius (default: 2)')
    parser.add_argument('--batch-size', type=int, default=10000,
                        help='Batch size for processing (default: 10000)')
    parser.add_argument('--verify', action='store_true',
                        help='Verify the generated fingerprints')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate fingerprints
    start_time = time.time()
    num_fps = generate_fingerprints(
        args.input,
        args.output,
        fp_size=args.fp_size,
        radius=args.radius,
        batch_size=args.batch_size
    )
    
    elapsed = time.time() - start_time
    print(f"\nTotal time: {elapsed:.2f} seconds")
    print(f"Speed: {num_fps / elapsed:.0f} molecules/second")
    
    # Verify if requested
    if args.verify:
        verify_fingerprints(args.output)

if __name__ == '__main__':
    main()