#!/usr/bin/env python3
"""Audit multi-seed ranking and accumulated candidates using saved model responses."""

import argparse
import csv
import gc
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem, DataStructs

from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.common import FP, canonical, digest, molecule, write_json
from chembounce_synthesis.conformers import pharmacophore_overlap
from chembounce_synthesis.continuation import continuation_inputs
from chembounce_synthesis.reactions import Rules
from chembounce_synthesis.routes import verify_route
from chembounce_synthesis.workflow import load_config, run


def read_table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Config with a saved provider and any desired core")
    parser.add_argument("--hit", required=True)
    parser.add_argument("--design", required=True)
    parser.add_argument("--include-run", help="Optional earlier pool to preserve from the first stage")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    config = load_config(args.config)
    if config["provider"].get("kind") != "saved":
        parser.error("Use saved proposals so the comparison does not resample the model")
    out = Path(args.output).resolve()
    if out.exists() and any(out.iterdir()):
        parser.error("Use a fresh output directory")
    out.mkdir(parents=True, exist_ok=True)
    config.setdefault("expansion", {}).update(building_block_top_k=16, rounds=1, max_routes=128)
    config.setdefault("conformers", {}).update(num_conformers=12, seed=42)
    config["refinement"] = {"top_n": 10, "num_conformers": 24, "seeds": [42, 43, 44]}
    config["max_scored_products"] = config["max_results"] = 64
    if args.include_run:
        config.setdefault("carry_forward", []).append(str(Path(args.include_run).resolve()))
    write_json(out / "first-config.json", config)
    first = run(args.hit, [args.design], config, out / "first")
    before = {p.name: digest(p) for p in (out / "first").iterdir()}
    gc.collect()
    hit, designs, continued = continuation_inputs(out / "first")
    continued["expansion"]["rounds"] = 2
    second = run(hit, designs, continued, out / "continued")
    assert before == {p.name: digest(p) for p in (out / "first").iterdir()}
    first_pool = {r["smiles"] for r in read_table(out / "first/screening.tsv")}
    second_pool = {r["smiles"] for r in read_table(out / "continued/screening.tsv")}
    assert first_pool <= second_pool, "A previously screened product was lost"
    old_pool = set()
    if args.include_run:
        old = Path(args.include_run)
        table = old / "screening.tsv" if (old / "screening.tsv").exists() else old / "results.tsv"
        old_pool = {r["smiles"] for r in read_table(table)}
        assert old_pool <= first_pool, "An earlier product was lost during first-stage carryover"

    # Build a small audit catalog by looking up selected IDs in the original TSV.
    # This rechecks membership without fingerprinting the entire catalog again.
    routes, ids = {}, set()
    for name, rows in (("first", first), ("continued", second)):
        records = json.loads((out / name / "routes.json").read_text())["routes"]
        lookup = {r["route_id"]: r for r in records}
        routes[name] = [lookup[r["route_id"]]["route"] for r in rows]
        ids.update(
            n["building_block_id"]
            for route in routes[name]
            for n in route["nodes"]
            if "building_block_id" in n
        )
    with (
        open(config["catalog"], newline="") as stream,
        (out / "audit-catalog.tsv").open("w", newline="") as dest,
    ):
        reader = csv.DictReader(stream, delimiter="\t")
        writer = csv.DictWriter(dest, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(r for r in reader if r["id"] in ids)
    catalog, rules = Catalog(out / "audit-catalog.tsv"), Rules(config["rules"])
    report = {
        "scope": "Single supplied case; saved model outputs; not held-out or experimental validation",
        "hit": canonical(args.hit),
        "design": canonical(args.design),
        "input_config_sha256": digest(args.config),
        "earlier_pool_size": len(old_pool),
        "earlier_pool_retained": old_pool <= first_pool,
        "first_pool_size": len(first_pool),
        "continued_pool_size": len(second_pool),
        "lost_products": sorted(first_pool - second_pool),
        "first_files_unchanged": True,
        "cases": {},
    }
    core = molecule(config["core_smiles"]) if config.get("core_smiles") else None
    for name, rows in (("first", first), ("continued", second)):
        errors = []
        for row, route in zip(rows, routes[name]):
            assert verify_route(route, catalog, rules)["product"] == row["smiles"]
            mol = molecule(row["smiles"])
            assert core is None or mol.HasSubstructMatch(core, useChirality=True)
            assert row["seed_count"] == row["seed_successes"] == 3
            assert row["seed_score_min"] <= row["score"] == row["seed_score_median"] <= row["seed_score_max"]
            for role, reference in (("hit", args.hit), ("design", args.design)):
                value = DataStructs.TanimotoSimilarity(
                    FP.GetFingerprint(mol), FP.GetFingerprint(molecule(reference))
                )
                assert abs(value - row[f"{role}_tanimoto"]) < 1e-12
        for role in ("hit", "design"):
            poses = list(Chem.SDMolSupplier(str(out / name / f"poses_{role}.sdf")))
            refs = list(Chem.SDMolSupplier(str(out / name / f"references_{role}.sdf")))
            assert len(poses) == len(refs) == len(rows)
            for row, pose, reference in zip(rows, poses, refs):
                assert pose.GetIntProp("pose_seed") == reference.GetIntProp("pose_seed") == row["pose_seed"]
                overlap, _ = pharmacophore_overlap(
                    reference, pose, 0, 0, sigma=config.get("conformers", {}).get("pharmacophore_sigma", 1.0)
                )
                errors.append(abs(overlap - row[f"{role}_pharmacophore"]))
        assert max(errors, default=0) < 1e-4
        manifest = json.loads((out / name / "manifest.json").read_text())
        summary = json.loads((out / name / "summary.json").read_text())
        report["cases"][name] = {
            "returned_products": len(rows),
            "best_row": rows[0] if rows else None,
            "all_selected_routes_replayed": True,
            "all_selected_cores_preserved": True,
            "tanimoto_recomputed": True,
            "max_sdf_pharmacophore_error": max(errors, default=0),
            "catalog_sha256": manifest["catalog_sha256"],
            "rules_sha256": manifest["rules_sha256"],
            "code_sha256": manifest["code_sha256"],
            "git_revision": manifest["git_revision"],
            "provider": manifest["provider"]["provider"],
            "summary": {k: v for k, v in summary.items() if k not in ("search_reports", "failures")},
        }
    write_json(out / "audit.json", report)
    print(json.dumps({k: v for k, v in report.items() if k != "cases"}, indent=2))


if __name__ == "__main__":
    main()
