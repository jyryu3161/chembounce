# ChemBounce

<!--<img src="./assets/graphical_abstract.png" width="450px"></img>-->

## Dependencies

Install conda environment as follows: 

```bash
conda env create -n ENVIRONMENT_NAME -f environment.yml
```

For download reference files,
```bash
bash install.sh
```


## How to run

Before running, make sure conda environment is activate `conda activate ENVIRONMENT_NAME`

```bash
python chembounce.py -o OUTPUT_DIRECTORY -i INPUT_SMILES -n FRAGMENT_MAX_N -t TANIMOTO_THRESHOLD
```
See the example of `run.sh`


### Parameters


- `-i INPUT_SMILES` : Input SMILES, the target molecular structure
- `-o OUTPUT_DIRECTORY` : Output location
- `-n FRAGMENT_MAX_N` : (`int`) "Maximal number of scaffold-hopped candidates for a fragment
- `-t TANIMOTO_THRESHOLD` : (`float`; `0` to `1`) Threshold for Tanimoto simiarity between `INPUT_SMILES` and generated SMILES. (default is `0.5`)


#### Options for iteration numbers

Limit numbers of iterations

- `--overall_max_n` : Maximal number of scaffold-hopped candidates for overall fragments
- `--scaffold_top_n` : Number of scaffolds to test for a fragment
- `--cand_max_n__rplc` : Maximal number of candidates for a replaced scaffold


#### Options for Thresholds

Thresholds for candidates

- `-t TANIMOTO_THRESHOLD` : (`float`; `0` to `1`) Threshold for Tanimoto simiarity between `INPUT_SMILES` and generated SMILES. (default is `0.5`)
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
- `--wo_lipinski` : Turn off  *Lipinski's rule of five*, which is `logp_max=5`, `qed_max=500`, `h_donor_max=5`, `h_acceptor_max=10`
This option will turn off the limitation of candidates by lipinski's rule of five
If one of the threhold category for this rule is additionally defined, the threshold of Lipinski's rule is ignored and replaced by the user-defined threshold.


#### Optional parameters

- `-l` : (Optional) Low memory mode
- `--fragments` : Specific fragment SMILES of the input structure to replace. For multiple fragment, repeatedly impose this option: e.g. `--fragments SMILES_A --fragments SMILES_B`
- `--replace_scaffold_files` :  Replace scaffold file, for the `fragments` option with tsv format, with or without priority score (the higher, with more priority).
  * NOTE: One of the result files, `replace_scaffold_list.fragment_XXXX.tsv`, can be used for `replace_scaffold_files`.
  * NOTE: for multiple `fragments` and its corresponding `replace_scaffold_files`, the order of `fragments` and `replace_scaffold_files` should be matched:
      e.g.: `--fragments SMILES_A --replace_scaffold_files FILE_A --fragments SMILES_B --replace_scaffold_files FILE_B --fragments SMILES_C --replace_scaffold_files FILE_C`
  * NOTE: For multiple 


## Results

Result file will is saved as `OUTPUT_DIRECTORY/overall_result.txt` with tap-delimitted format. `Final structure` column is the final result SMILES of scaffold hopping and `Standardized final structure` is standardized SMILES format for the corresponding result. For sub-data,

- Fragments of the input molecule is saved as `~/fragment_info.tsv`
- Considered replace scaffold list for each fragment is saved as `~/replace_scaffold_list.fragment_XXXX.tsv`.
- Results for each fragment is saved `~/fragment_result_XXXX.txt`.

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

For more information : check [CSB_Lab](https://www.csb-lab.net/)

