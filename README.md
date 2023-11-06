# scaffoldhopper
scaffold hopping

# Installation

For download reference files, run `bash install.sh`
Install conda environment as follows: `conda create env -f environment.yml`


# How to run

Before running, make sure conda environment is activate `conda activate scaffold`

`python run_scaffold_hopper.py -o OUTPUT_DIRECTORY -i INPUT_SMILES -n NUMBER_OF_STRUCTURES -t SIMILARITY_THRESHOLD`

`OUTPUT_DIRECTORY` Location of output directory
`INPUT_SMILES` SMILES structure for scaffold hopping
`NUMBER_OF_STRUCTURES` (int) Number of structure to generate by scaffold hopping
`SIMILARITY_THRESHOLD` (float; 0 to 1) Threshold for Tanimoto simiarity between `INPUT_SMILES` and generated SMILES. (default is 0.7)


# Result

Result file will be generated as follows: `OUTPUT_DIRECTORY/result.txt` with tap-delimitted format. `Final structure` column is the final result SMILES of scaffold hopping and `Standardized final structure` is standardized SMILES format for the corresponding result.





