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
import argparse
import time

def generate_fingerprints(input_file, output_file, fp_size=2048, radius=2, batch_size=10000):
    """
    Generate Morgan fingerprints for all molecules in the input file
    
    Args:
        input_file: Path to scaffold SMILES file
        output_file: Path to save fingerprint data
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
    
    # Convert fingerprints to numpy array for faster operations
    print("Converting to numpy array...")
    fp_array = np.zeros((len(all_fps), fp_size), dtype=np.uint8)
    
    for i, fp in enumerate(tqdm.tqdm(all_fps, desc="Converting to array")):
        arr = np.zeros((fp_size,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fp_array[i] = arr
    
    # Save data
    print(f"Saving to {output_file}")
    data = {
        'fingerprints': fp_array,
        'smiles': all_smiles,
        'valid_indices': valid_indices,
        'fp_size': fp_size,
        'radius': radius,
        'version': '1.0'
    }
    
    # Use numpy's compressed format for efficient storage
    np.savez_compressed(output_file, **data)
    
    # Print statistics
    file_size = os.path.getsize(output_file) / (1024 * 1024)  # MB
    print(f"Fingerprint file size: {file_size:.2f} MB")
    print(f"Compression ratio: {file_size / (len(all_fps) * fp_size / 8 / 1024 / 1024):.2%}")
    
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