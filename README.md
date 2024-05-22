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
python chembounce.py -o OUTPUT_DIRECTORY -i INPUT_SMILES -c CORE_SMILES -n NUMBER_OF_FRAGMENTS -t SIMILARITY_THRESHOLD
```
See the example of `run.sh`


### Parameters

- `-o OUTPUT_DIRECTORY` : Output location
- `-i INPUT_SMILES` : Input SMILES, the target molecular structure
- `-c CORE_SMILES` : (Optional) Core SMILES which should not be altered while scaffold hopping
- `-n NUMBER_OF_FRAGMENTS` : (`int`) Number of fragments to test for scaffold hopping
- `-t SIMILARITY_THRESHOLD` : (`float`; `0` to `1`) Threshold for Tanimoto simiarity between `INPUT_SMILES` and generated SMILES. (default is `0.5`)
- `-l` : (Optional) Low memory mode.


## Result

Result file will be generated as follows: `OUTPUT_DIRECTORY/result.txt` with tap-delimitted format. `Final structure` column is the final result SMILES of scaffold hopping and `Standardized final structure` is standardized SMILES format for the corresponding result.


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

