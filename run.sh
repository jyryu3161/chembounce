#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "${PYTHON:-python}" "${script_dir}/chembounce.py" \
  -o "${1:-${script_dir}/output_example}" \
  -i 'CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)CO)Cl' \
  -n 100 -t 0.5 --cand_max_n__rplc 10
