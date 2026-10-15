#!/usr/bin/env python3
"""Audit known-pair admission and 2D gating, not blind generation performance.

Input: JSON list with reference_smiles, hop_smiles and target identifiers. Each
known hop is supplied explicitly as a replacement of the full reference. This
isolates candidate exclusion by the two policies; it does not establish that
ChemBounce can generate an independently held-out active structure.
"""

import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from chembounce import _search_diverse_scaffolds, chembounce, search_similar_scaffolds


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)
    pairs = json.loads(args.pairs.read_text())
    library = sorted(
        {
            Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(pair[key])))
            for pair in pairs
            for key in ("reference_smiles", "hop_smiles")
        }
    )
    rows = []
    for pair in pairs:
        folder = args.output / pair["target"]
        folder.mkdir(parents=True, exist_ok=True)
        replacements = folder / "provided-hop.smi"
        replacements.write_text(pair["hop_smiles"] + "\n")
        reference_scaffold = MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(pair["reference_smiles"]))
        hop_scaffold = Chem.MolToSmiles(
            MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(pair["hop_smiles"]))
        )
        legacy_scaffolds = search_similar_scaffolds(
            reference_scaffold, library, scaffold_top_n=3, threshold=0.3
        )
        diverse_scaffolds, selection = _search_diverse_scaffolds(
            reference_scaffold, library, 3, 0.3, args.seed, 0.5
        )
        record = {
            "target": pair["target"],
            "reference_smiles": pair["reference_smiles"],
            "hop_smiles": pair["hop_smiles"],
            "legacy_scaffold_admitted": hop_scaffold in legacy_scaffolds.index,
            "diverse_scaffold_admitted": hop_scaffold in diverse_scaffolds.index,
            "scaffold_pool_size": len(library),
            "diverse_scaffold_selection": selection,
        }
        for policy in ("legacy", "diverse"):
            results, _ = chembounce(
                pair["reference_smiles"],
                output_dir=folder / policy,
                fragments=[pair["reference_smiles"]],
                replace_scaffold_files=[replacements],
                candidate_thresholds={},
                calculate_shape=False,
                overall_max_n=1,
                frag_max_n=1,
                scaffold_top_n=1,
                cand_max_n__rplc=1,
                _merge_structure_top_n_=1,
                tanimoto_threshold=0.5,
                search_policy=policy,
                exploration_fraction=0.5,
                search_candidate_budget=1,
                seed=args.seed,
            )
            expected = Chem.MolToSmiles(Chem.MolFromSmiles(pair["hop_smiles"]))
            record[f"{policy}_provided_candidate_admitted"] = expected in set(results["Final structure"])
        rows.append(record)
    report = {
        "schema_version": 1,
        "scope": "Known-pair admissibility/filter diagnostic with explicitly supplied candidates and a six-scaffold pair-derived library. Not independent generation, activity prediction, generalization, or scaffold recovery performance.",
        "parameters": {
            "candidate_threshold": 0.5,
            "scaffold_threshold": 0.3,
            "scaffold_top_n": 3,
            "exploration_fraction": 0.5,
            "seed": args.seed,
        },
        "pairs": rows,
    }
    (args.output / "admission-audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                "pair_count": len(rows),
                **{
                    f"{policy}_provided_candidate_admitted": sum(
                        row[f"{policy}_provided_candidate_admitted"] for row in rows
                    )
                    for policy in ("legacy", "diverse")
                },
            }
        )
    )


if __name__ == "__main__":
    main()
