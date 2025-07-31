# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## ChemBounce Overview

ChemBounce is a computational chemistry tool for scaffold hopping in drug discovery. It systematically replaces molecular fragments to generate novel chemical structures while maintaining desired properties.

## Development Commands

### Environment Setup
```bash
# Create conda environment from environment.yml
conda env create -f environment.yml

# Activate environment
conda activate chembounce

# Download required reference data (first time only)
bash install.sh
```

### Running ChemBounce
```bash
# Basic run
python chembounce.py -o OUTPUT_DIR -i "INPUT_SMILES" -n 1000 -t 0.5

# Example run (see run.sh)
python chembounce.py \
    -o ./test_output \
    -i "CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)CO)Cl" \
    -n 1000 -t 0.5 \
    --cand_max_n__rplc 10
```

### Running in Google Colab
ChemBounce provides a Jupyter notebook for Google Colab usage: `chembounce_colab.ipynb`

## Architecture Overview

### Core Components

1. **chembounce.py** - Main entry point and orchestration
   - `chembounce()` function decorated with `@cost_estimation.EstDecoMem` for performance tracking
   - Handles fragment decomposition, scaffold replacement, and structure recombination
   - Manages filtering based on molecular property thresholds

2. **scaffoldgraph_fragmenter.py** - Fragment generation
   - Modified from ScaffoldGraph library
   - `get_all_murcko_fragments()` - Generates Murcko fragments with iteration control
   - Supports breaking or preserving fused ring systems

3. **utils.py** - Core molecular operations
   - Molecular similarity calculations (Tanimoto, shape, electroshape, CATS descriptors)
   - Property calculations (QED, SA score, logP, MW, H-donors/acceptors)
   - Fragment database loading from `data/` directory
   - Structure manipulation and standardization
   - **New**: Fast fingerprint-based similarity search functions
     - `load_fingerprint_db()` - Loads pre-computed fingerprints
     - `search_similar_scaffolds_fast()` - Bulk similarity search using numpy
     - `calculate_bulk_tanimoto()` - Optimized Tanimoto calculations

4. **cli.py** - Command line interface
   - Argument parsing with extensive threshold and iteration control options
   - Supports Lipinski's rule of five filtering (can be disabled)

5. **cost_estimation.py** - Performance monitoring
   - Tracks CPU time, memory usage, and elapsed time
   - Decorator-based implementation for clean integration

6. **cats_module.py** - CATS (Chemically Advanced Template Search) descriptors
   - Used for molecular similarity assessment

7. **generate_fingerprints.py** - Pre-compute fingerprints for fast search
   - Generates Morgan fingerprints for all scaffolds
   - Saves as compressed numpy array
   - Enables 100-1000x faster similarity search

### Data Flow

1. Input SMILES → Fragment decomposition (scaffoldgraph_fragmenter)
2. Each fragment → Scaffold database search (utils)
3. Candidate scaffolds → Structure recombination (chembounce)
4. Combined structures → Property filtering (utils)
5. Valid candidates → Similarity scoring and ranking
6. Results → Output files with detailed metrics

### Key Algorithms

- **Fragment Decomposition**: Uses Murcko scaffold decomposition with configurable iteration rounds
- **Scaffold Replacement**: Searches fragment database using similarity metrics
- **Structure Merging**: Combines original molecule with replacement scaffolds at attachment points
- **Filtering Pipeline**: 
  - Molecular property thresholds (MW, logP, QED, SA score, H-donors/acceptors)
  - Tanimoto similarity threshold
  - Optional Lipinski's rule of five

### Output Structure

- `overall_result.txt` - Main results with final structures and standardized SMILES
- `fragment_info.tsv` - Decomposed fragments from input
- `replace_scaffold_list.fragment_XXXX.tsv` - Candidate scaffolds per fragment
- `fragment_result_XXXX.txt` - Detailed results per fragment
- `resource_cost.json` - Performance metrics

## Important Dependencies

- RDKit (2020.09.5) - Core cheminformatics library
- ScaffoldGraph - Fragment decomposition
- ODDT - Shape similarity calculations
- PubChemPy - Chemical data retrieval
- MolVS - Molecule validation and standardization

## Notes on Memory and Performance

### Fingerprint-based Fast Search (New)
- **Major Performance Improvement**: Pre-computed fingerprints enable 100-1000x faster similarity search
- Generate fingerprints once: `bash prepare_fingerprints.sh` (takes several hours for 5M scaffolds)
- Uses optimized numpy operations instead of iterating through 5 million scaffolds one by one
- Enabled by default, disable with `--no-fingerprint-db` if needed

### General Performance Tips
- Use `-l` flag for low memory mode when processing large datasets
- Fragment iteration parameters control search space:
  - `--overall_max_n` - Total candidates across all fragments
  - `--frag_max_n` - Candidates per fragment
  - `--scaffold_top_n` - Scaffolds to test per fragment
  - `--cand_max_n__rplc` - Candidates per scaffold replacement