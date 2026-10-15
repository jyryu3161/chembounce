"""Build portable, hashed method configurations from local pinned assets.

Variants are explicit JSON settings overrides, not hidden changes in the runner.
A saved proposal provider is permitted for diagnostic replay but cannot count as
a live-inference wall-time comparison; its provenance remains in run artifacts.
"""

import argparse
import copy
import json
import os
from pathlib import Path

from benchmarks.research.protocol import digest, write_json, synthesis_asset_paths


def artifact(path, base):
    path = Path(path).resolve()
    if path.is_dir():
        return {
            "directory": os.path.relpath(path, base),
            "files_sha256": {str(p.relative_to(path)): digest(p) for p in path.rglob("*") if p.is_file()},
        }
    return {"artifact": os.path.relpath(path, base), "sha256": digest(path)}


def build(args):
    out = args.output.resolve()
    if out.exists():
        raise ValueError("Protocol already exists; preserve previous experiment configuration")
    out.parent.mkdir(parents=True, exist_ok=True)
    shared = {"catalog": artifact(args.catalog, out.parent)}
    if args.rules:
        shared["rules"] = artifact(args.rules, out.parent)
    methods = []
    caps = {
        "nearest_catalog": args.max_evaluations,
        "nearest_3d": 256,
        "reaction_enumeration": 4096,
        "synllama_only": 4096,
        "pipeline": 256,
    }
    if args.method_evaluation_caps:
        caps.update(json.loads(args.method_evaluation_caps.read_text()))
    for kind in args.methods:
        method = {
            "id": kind,
            "kind": kind,
            "assets": copy.deepcopy(shared),
            "settings": {},
            "budget_overrides": {"max_evaluations": min(caps[kind], args.max_evaluations)},
        }
        if kind in {"nearest_catalog", "nearest_3d"}:
            method["assets"].pop("rules", None)
            if args.retrieval_library:
                method["assets"]["catalog"] = artifact(args.retrieval_library, out.parent)
        if kind in {"synllama_only", "pipeline"}:
            method["inventories"] = [
                {
                    "kind": kind,
                    "status": "unknown",
                    "reason": "Upstream complete chemical inventory unavailable; no leakage-free pretraining claim",
                }
                for kind in ("training", "prompt")
            ]
            if args.prompt_inventory:
                method["inventories"][1] = artifact(args.prompt_inventory, out.parent) | {
                    "kind": "prompt",
                    "status": "available",
                    "coverage": "complete",
                }
        if kind == "synllama_only":
            if not args.model_directory or not args.prompt:
                raise ValueError("synllama_only requires model directory and prompt")
            directory = args.model_directory.resolve()
            method["assets"]["model_directory"] = {
                "directory": os.path.relpath(directory, out.parent),
                "files_sha256": {
                    str(p.relative_to(directory)): digest(p) for p in directory.rglob("*") if p.is_file()
                },
            }
            method["assets"]["prompt"] = artifact(args.prompt, out.parent)
            method["settings"] = {"device": args.device, "samples": args.samples}
        if kind == "pipeline":
            if not args.synthesis_config or not args.scaffold_db:
                raise ValueError("pipeline requires synthesis config and scaffold DB")
            method["assets"].update(
                synthesis_config=artifact(args.synthesis_config, out.parent),
                scaffold_db=artifact(args.scaffold_db, out.parent),
            )
            if args.fingerprint_db:
                method["assets"]["fingerprint_db"] = artifact(args.fingerprint_db, out.parent)
            for name, path in synthesis_asset_paths(args.synthesis_config).items():
                method["assets"][name] = artifact(path, out.parent)
            method["settings"] = {
                "generator_options": {
                    "calculate_shape": False,
                    "candidate_thresholds": {},
                    "search_policy": "diverse",
                    "exploration_fraction": 0.25,
                    "frag_max_n": args.max_candidates,
                    "scaffold_top_n": args.max_candidates,
                }
            }
            if args.generator_options:
                method["settings"]["generator_options"].update(json.loads(args.generator_options.read_text()))
        methods.append(method)
    if args.variants:
        variants = json.loads(args.variants.read_text())
        lookup = {m["id"]: m for m in methods}
        for variant in variants:
            variant_kind = variant.get("variant_kind", "single_factor_ablation")
            if variant_kind not in {"single_factor_ablation", "compound_baseline"}:
                raise ValueError("Unknown comparison variant_kind")
            parent = lookup[variant["parent"]]
            method = copy.deepcopy(parent)
            method["id"] = variant["id"]
            path_overrides = variant["settings"].get("synthesis_overrides", {})
            if set(path_overrides) & {"provider", "catalog", "rules", "reference_sdf", "carry_forward"}:
                raise ValueError(
                    "Variants cannot replace pinned external assets; create a separately pinned parent"
                )
            if path_overrides.get("synthesis_evidence", {}).get("predictor"):
                raise ValueError("Variant predictor replacement requires a separately pinned parent")
            for name, value in variant["settings"].items():
                if isinstance(value, dict):
                    method["settings"].setdefault(name, {}).update(value)
                else:
                    method["settings"][name] = value
            method["ablation" if variant_kind == "single_factor_ablation" else "compound_baseline"] = {
                "kind": variant_kind,
                "parent": parent["id"],
                "description": variant["description"],
                "variant_file_sha256": digest(args.variants),
            }
            methods.append(method)
    protocol = {
        "schema_version": 1,
        "seeds": args.seeds,
        "budget": {
            "wall_seconds": args.wall_seconds,
            "max_candidates": args.max_candidates,
            "max_evaluations": args.max_evaluations,
            "max_memory_mb": args.max_memory_mb,
        },
        "search_space_audit_max_records": args.search_space_audit_max_records,
        "methods": methods,
        "planned_unexecuted_methods": [
            {
                "id": "rush",
                "source_url": "https://github.com/molML/RUSH",
                "reason": "Published implementation requires its configured environment and ROCS license; use external adapter only when actually available",
            }
        ],
    }
    write_json(out, protocol)
    return protocol


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--rules", type=Path)
    parser.add_argument(
        "--retrieval-library",
        type=Path,
        help="Train-only end-product library for nearest2D/3D; reaction methods retain building-block catalog",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=("nearest_catalog", "nearest_3d", "reaction_enumeration", "synllama_only", "pipeline"),
        default=["nearest_catalog", "reaction_enumeration"],
    )
    parser.add_argument("--scaffold-db", type=Path)
    parser.add_argument("--fingerprint-db", type=Path)
    parser.add_argument("--synthesis-config", type=Path)
    parser.add_argument("--model-directory", type=Path)
    parser.add_argument("--prompt", type=Path)
    parser.add_argument("--prompt-inventory", type=Path)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--samples", type=int, default=4)
    parser.add_argument("--generator-options", type=Path)
    parser.add_argument(
        "--variants", type=Path, help="JSON list of named parent/settings/description ablations"
    )
    parser.add_argument("--seeds", type=int, nargs="+", default=[42, 43, 44])
    parser.add_argument("--wall-seconds", type=float, default=300)
    parser.add_argument("--max-candidates", type=int, default=20)
    parser.add_argument("--max-evaluations", type=int, default=1_000_000)
    parser.add_argument(
        "--method-evaluation-caps", type=Path, help="JSON kind-to-count caps; wall/output caps remain common"
    )
    parser.add_argument("--max-memory-mb", type=float, default=24000)
    parser.add_argument("--search-space-audit-max-records", type=int, default=5000)
    args = parser.parse_args(argv)
    if any(k not in {"nearest_catalog", "nearest_3d"} for k in args.methods) and not args.rules:
        parser.error("Reaction and synthesis methods require --rules")
    build(args)


if __name__ == "__main__":
    main()
