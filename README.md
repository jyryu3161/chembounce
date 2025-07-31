# ChemBounce

<!--<img src="./assets/graphical_abstract.png" width="450px"></img>-->


## Using in Google Colaboratory

ChemBounce is also available in [Google Colaboratory](https://colab.research.google.com/github/jyryu3161/chembounce/blob/main/chembounce_colab.ipynb)

## Dependencies

Install the conda environment as follows: 

```bash
conda env create -f environment.yml
```

For a simpler installation with more flexible dependencies:
```bash
conda env create -f environment_simple.yml
```

**Difference between environment files:**
- `environment.yml`: Full environment with exact versions (strict, ensures reproducibility)
- `environment_simple.yml`: Core dependencies only with flexible versions (easier installation, fewer conflicts)

To download reference files:
```bash
bash install.sh
```

### Performance Optimization

ChemBounce uses pre-computed fingerprints for **100-1000x faster** similarity search. 

**The default scaffold database already includes pre-computed fingerprints** (`data/scaffold_fingerprints.npz`), so you don't need to generate them.

#### Using Custom Scaffold Databases

If you want to use your own scaffold database, you need to generate fingerprints for it:

```bash
# Generate fingerprint database for custom scaffold database
bash prepare_fingerprints.sh /path/to/custom_scaffolds.txt /path/to/output_fingerprints.npz
```

This creates a fingerprint database that dramatically speeds up scaffold searching. The original implementation searches through 5 million scaffolds one by one, while the fingerprint-based search uses optimized bulk operations.

## How to run

Before running, make sure the conda environment is activated: `conda activate chembounce`

```bash
python chembounce.py -o OUTPUT_DIRECTORY -i INPUT_SMILES -n FRAGMENT_MAX_N -t TANIMOTO_THRESHOLD
```
Example usage:

```
python chembounce.py \
    -o ./output_test \
    -i "CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)CO)Cl" \
    -n 100 -t 0.5 \
    --cand_max_n__rplc 10
```

Example with custom scaffold database:

```
# First, generate fingerprints for your custom scaffold database
bash prepare_fingerprints.sh my_scaffolds.txt my_scaffolds_fingerprints.npz

# Then use it with ChemBounce
python chembounce.py \
    -o ./output_test \
    -i "CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)CO)Cl" \
    -n 100 -t 0.5 \
    --scaffold-db my_scaffolds.txt \
    --fingerprint-db my_scaffolds_fingerprints.npz
```

### Parameters


- `-i INPUT_SMILES` : Input SMILES, the target molecular structure.
- `-o OUTPUT_DIRECTORY` : Output location
- `-n FRAGMENT_MAX_N` : (`int`) Maximal number of scaffold-hopped candidates for a fragment
- `-t TANIMOTO_THRESHOLD` : (`float`; `0` to `1`) Threshold for Tanimoto similarity between `INPUT_SMILES` and generated SMILES. (default is `0.5`)
- `--core_smiles CORE_SMILES` : SMILES of core structure for the `INPUT_SMILES`, which should not be changed. If the `CORE_SMILES` is not in the `INPUT_SMILES`, it would raise an error.


#### Options for iteration numbers

Limit the number of iterations

- `--overall_max_n` : Maximal number of scaffold-hopped candidates for overall fragments
- `--scaffold_top_n` : Number of scaffolds to test for a fragment
- `--cand_max_n__rplc` : Maximal number of candidates for a replaced scaffold


#### Options for Thresholds

Thresholds for candidates

- `-t TANIMOTO_THRESHOLD` : (`float`; `0` to `1`) Threshold for Tanimoto similarity between `INPUT_SMILES` and generated SMILES. (default is `0.5`)
- `--qed_min QED_MIN` : (`float`) Minimal QED score
- `--qed_max QED_MAX` : (`float`) Maximal QED score
- `--sa_min SA_MIN` : (`float`; `0.0` to `10.0`) Minimal SAscore
- `--sa_max SA_MAX` : (`float`) Maximal SAscore
- `--logp_min LOGP_MIN` : (`float`) Minimal logP limit
- `--logp_max LOGP_MAX` : (`float`) Maximal logP limit
- `--mw_min MW_MIN` : (`float`; above `0.0`) Minimal molecular weight
- `--mw_max MW_MAX` : (`float`; above `0.0`) Maximal molecular weight
- `--h_donor_min H_DONOR_MIN` : (above `0`) Minimal H donor limit
- `--h_donor_max H_DONOR_MAX` : (above `0`) Maximal H donor limit
- `--h_acceptor_min H_ACCEPTOR_MIN` : (above `0`) Minimal H acceptor limit
- `--h_acceptor_max H_ACCEPTOR_MAX` : (above `0`) Maximal H acceptor limit
- `--wo_lipinski` : Turn off *Lipinski's rule of five*, which is `logp_max=5`, `mw_max=500`, `h_donor_max=5`, `h_acceptor_max=10`
This option will turn off the limitation of candidates by Lipinski's rule of five.
If one of the threshold categories for this rule is additionally defined, the threshold of Lipinski's rule is ignored and replaced by the user-defined threshold.


#### Optional parameters

- `-l` : (Optional) Low memory mode
- `--use-fingerprint-db` : Use pre-computed fingerprint database for fast similarity search (default: enabled)
- `--no-fingerprint-db` : Disable fingerprint database and use the original slow similarity search
- `--scaffold-db` : Path to custom scaffold database file (SMILES format, one per line)
- `--fingerprint-db` : Path to custom fingerprint database file (NPZ format) corresponding to the custom scaffold database
- `--fragments` : Specific fragment SMILES of the input structure to replace. For multiple fragments, repeatedly use this option: e.g. `--fragments SMILES_A --fragments SMILES_B`
- `--replace_scaffold_files` : User-defined library (if solely imposed). OR Replace scaffold file for the `fragments` option with TSV format, with or without priority score (the higher, the more priority).
  * NOTE: One of the result files, `replace_scaffold_list.fragment_XXXX.tsv`, can be used for `replace_scaffold_files`.
  * NOTE: For multiple `fragments` and their corresponding `replace_scaffold_files`, the order of `fragments` and `replace_scaffold_files` should be matched:
      e.g.: `--fragments SMILES_A --replace_scaffold_files FILE_A --fragments SMILES_B --replace_scaffold_files FILE_B --fragments SMILES_C --replace_scaffold_files FILE_C`


## Results

The result file is saved as `OUTPUT_DIRECTORY/overall_result.txt` in tab-delimited format. The `Final structure` column contains the final result SMILES of scaffold hopping and `Standardized final structure` is the standardized SMILES format for the corresponding result. For sub-data:

- Fragments of the input molecule are saved as `~/fragment_info.tsv`
- Considered replacement scaffold lists for each fragment are saved as `~/replace_scaffold_list.fragment_XXXX.tsv`
- Results for each fragment are saved as `~/fragment_result_XXXX.txt`

## Citation

<!--```bibtex
@article{,
  title    = "",
  author   = "",
  journal  = "",
  month    = "",
  year     =  2024
}
```
-->

## Contact

For more information, check [CSB_Lab](https://www.csb-lab.net/)

