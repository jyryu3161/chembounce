#!/bin/bash
# Script to generate fingerprint database for ChemBounce

echo "ChemBounce Fingerprint Database Generator"
echo "========================================"
echo ""
echo "This script will generate a pre-computed fingerprint database"
echo "for fast similarity search in ChemBounce."
echo ""
echo "Usage: bash prepare_fingerprints.sh [input_file] [output_file]"
echo "  input_file: Path to scaffold SMILES file (default: data/Scaffolds_processed.txt)"
echo "  output_file: Path to save fingerprints (default: data/scaffold_fingerprints.npz)"
echo ""

# Parse command line arguments
INPUT_FILE="${1:-data/Scaffolds_processed.txt}"
OUTPUT_FILE="${2:-data/scaffold_fingerprints.npz}"

# If only one argument is provided, assume it's the input file and generate output name
if [ $# -eq 1 ]; then
    # Generate output filename based on input filename
    OUTPUT_DIR=$(dirname "$INPUT_FILE")
    INPUT_BASENAME=$(basename "$INPUT_FILE" | sed 's/\.[^.]*$//')
    OUTPUT_FILE="${OUTPUT_DIR}/${INPUT_BASENAME}_fingerprints.npz"
fi

echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo ""

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Scaffold database not found at $INPUT_FILE!"
    if [ "$INPUT_FILE" = "data/Scaffolds_processed.txt" ]; then
        echo "Please run 'bash install.sh' first to download the default scaffold database."
    else
        echo "Please provide a valid input file path."
    fi
    exit 1
fi

# Check if fingerprint already exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "Fingerprint database already exists at $OUTPUT_FILE"
    read -p "Do you want to regenerate it? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 0
    fi
fi

# Activate conda environment if not already active
if [[ "$CONDA_DEFAULT_ENV" != "chembounce" ]]; then
    echo "Activating chembounce conda environment..."
    conda activate chembounce
fi

# Run fingerprint generation
echo "Starting fingerprint generation..."
echo "This may take several hours. You can monitor progress below:"
echo ""

python generate_fingerprints.py \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_FILE" \
    --fp-size 2048 \
    --radius 2 \
    --batch-size 10000 \
    --verify

echo ""
echo "Fingerprint generation complete!"
echo "Output saved to: $OUTPUT_FILE"
echo ""