#!/usr/bin/env python3
"""Compare small/broader BB expansion using identical saved proposals and 3D settings.

Use a saved provider so differences are attributable to expansion, not resampling.
Counts and scores are diagnostics for the supplied case, not a success-rate study.
"""

import argparse
import copy
import gc
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem
from chembounce_synthesis.common import canonical, digest, write_json
from chembounce_synthesis.workflow import load_config, run


def isotope_parent(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--hit", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-scored", type=int, default=64)
    args = parser.parse_args()
    config = load_config(args.config)
    if config["provider"]["kind"] != "saved":
        parser.error("Use saved proposals to isolate expansion from model sampling")
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        parser.error("Output must be a fresh directory")
    out.mkdir(parents=True, exist_ok=True)
    result = {
        "scope": "Single-case expansion ablation with fixed saved proposals; not a held-out benchmark",
        "input_config_sha256": digest(args.config),
        "hit": canonical(args.hit),
        "design": canonical(args.design),
        "cases": {},
    }
    for name, expansion in (
        ("small", {"building_block_top_k": 3, "rounds": 1, "max_routes": 12}),
        ("broader", {"building_block_top_k": 16, "rounds": 2, "max_routes": 128}),
    ):
        settings = copy.deepcopy(config)
        settings.setdefault("expansion", {}).update(expansion)
        settings["minimum_scores"] = {}
        settings["max_scored_products"] = args.max_scored
        settings["max_results"] = args.max_scored
        write_json(out / f"{name}-config.json", settings)
        rows = run(args.hit, [args.design], settings, out / name)
        summary = json.loads((out / name / "summary.json").read_text())
        manifest = json.loads((out / name / "manifest.json").read_text())
        distinct = [r for r in rows if isotope_parent(r["smiles"]) != isotope_parent(args.design)]
        core = Chem.MolFromSmiles(settings["core_smiles"]) if settings.get("core_smiles") else None
        result["cases"][name] = {
            "expansion": expansion,
            "code_revision": manifest["git_revision"],
            "code_sha256": manifest["code_sha256"],
            "catalog_sha256": manifest["catalog_sha256"],
            "rules_sha256": manifest["rules_sha256"],
            "proposals_sha256": digest(settings["provider"]["path"]),
            "conformers": manifest["effective_conformers"],
            "returned_products": len(rows),
            "nonisotopic_analogues": len(distinct),
            "all_preserve_core": core is None
            or all(Chem.MolFromSmiles(r["smiles"]).HasSubstructMatch(core, useChirality=True) for r in rows),
            "scored_products": summary["scored_products"],
            "limits_reached": summary["limits_reached"],
            "elapsed_seconds": summary["elapsed_seconds"],
            "top_rows": rows[:5],
            "best_design_shape_row": max(distinct, key=lambda r: r["design_shape"], default=None),
            "best_joint_3d_row": max(distinct, key=lambda r: r["score"], default=None),
        }
        write_json(out / "comparison.json", result)
        print(
            json.dumps(
                {
                    "case": name,
                    "returned": len(rows),
                    "nonisotopic": len(distinct),
                    "elapsed_seconds": summary["elapsed_seconds"],
                }
            ),
            flush=True,
        )
        gc.collect()


if __name__ == "__main__":
    main()
