#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
input_file="${1:-${script_dir}/data/Scaffolds_processed.txt}"
output_file="${2:-${input_file%.*}_fingerprints.npz}"
exec "${PYTHON:-python}" "${script_dir}/generate_fingerprints.py" \
  --input "$input_file" --output "$output_file" --fp-size 2048 --radius 2 --verify
