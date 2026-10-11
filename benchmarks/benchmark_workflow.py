#!/usr/bin/env python3
"""Reproducible multi-input workflow evaluation; records zero-result cases too."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
from rdkit import Chem
from chembounce import chembounce


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fingerprint-db")
    parser.add_argument("--scaffold-db")
    parser.add_argument("--trust-legacy-db", action="store_true")
    parser.add_argument("--output", default="benchmark-results/workflow")
    parser.add_argument("--shape", action="store_true")
    args = parser.parse_args()
    root = Path(args.output)
    root.mkdir(parents=True, exist_ok=True)
    panel = json.loads(Path(__file__).with_name("molecules.json").read_text())
    reports = []
    for record in panel:
        out = root / record["name"]
        result, cost = chembounce(
            record["smiles"],
            output_dir=out,
            fingerprint_db_file=args.fingerprint_db,
            scaffold_db_file=args.scaffold_db,
            trust_legacy_db=args.trust_legacy_db,
            candidate_thresholds={},
            frag_max_n=10,
            overall_max_n=20,
            scaffold_top_n=20,
            calculate_shape=args.shape,
            seed=42,
            _merge_structure_top_n_=50,
        )
        molecules = [Chem.MolFromSmiles(s) for s in result["Final structure"]]
        assert all(m is not None and len(Chem.GetMolFrags(m)) == 1 for m in molecules)
        assert result["Final structure"].is_unique
        assert len(result) <= 20
        assert result["Tanimoto Similarity"].ge(0.5).all()
        (out / "resource_cost.json").write_text(json.dumps(cost, indent=2) + "\n")
        summary = json.loads((out / "run_summary.json").read_text())
        reports.append(
            {
                "name": record["name"],
                "status": summary["status"],
                "candidates": len(result),
                "seconds": cost["elapsed_time"],
                "sampled_peak_rss_mb": cost["max_memory_mb"],
                "shape_failures": summary["counts"].get("shape_failed", 0),
                "manifest": str(out / "run_manifest.json"),
            }
        )
    pd.DataFrame(reports).to_csv(root / "summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(reports).to_string(index=False))


if __name__ == "__main__":
    main()
