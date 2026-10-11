#!/usr/bin/env python3
"""Describe generated candidates; these statistics do not measure biological utility.

Scaffold novelty means a different canonical, stereo-free Murcko scaffold from
that of the input, not novelty against patents, known compounds or the library.
"""

import argparse
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid generated SMILES: {smiles}")
    core = MurckoScaffold.GetScaffoldForMol(mol)
    Chem.RemoveStereochemistry(core)
    return Chem.MolToSmiles(core)


def summarize(root):
    rows = []
    for directory in sorted(Path(root).iterdir()):
        if not directory.is_dir():
            continue
        manifest = directory / "run_manifest.json"
        evaluation = directory / "evaluation_status.json"
        if evaluation.exists():
            status = json.loads(evaluation.read_text())["status"]
            if status not in ("success", "no_candidates"):
                rows.append({"input": directory.name, "status": status, "candidates": None})
                continue
        if not manifest.exists():
            continue
        info = json.loads(manifest.read_text())
        results = pd.read_csv(manifest.parent / "overall_result.txt", sep="\t")
        original = info["parameters"]["target_smiles"]
        original_core = scaffold(original)
        cores = [scaffold(s) for s in results["Final structure"]]
        rows.append(
            {
                "input": manifest.parent.name,
                "status": json.loads((directory / "run_summary.json").read_text())["status"],
                "candidates": len(results),
                "unique_murcko_scaffolds": len(set(cores)),
                "fraction_changed_murcko": sum(c != original_core for c in cores) / len(cores)
                if cores
                else None,
                "median_tanimoto": results["Tanimoto Similarity"].median(),
                "median_qed": results["QED"].median(),
                "median_sa": results["SAscore"].median(),
                "min_mw": results["MW"].min(),
                "max_mw": results["MW"].max(),
            }
        )
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory")
    args = parser.parse_args()
    table = summarize(args.directory)
    table.to_csv(Path(args.directory) / "candidate_statistics.tsv", sep="\t", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
