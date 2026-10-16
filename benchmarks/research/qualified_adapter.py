"""Bounded native design supply and common qualified assessment for fresh runs.

Invoked only by qualified_runner; historical adapters and studies stay unchanged.
"""

import argparse
import copy
import csv
import json
from pathlib import Path
import time

from benchmarks.research.protocol import digest, write_json

ARMS = {"full", "baseline", "crem", "mol2mol", "synllama"}


def read(path):
    def invalid(value):
        raise ValueError(f"Nonfinite JSON: {value}")

    path = Path(path)
    if path.stat().st_size > 100_000_000:
        raise ValueError("JSON artifact exceeds the 100 MB integrity-read limit")
    return json.loads(path.read_text(), parse_constant=invalid)


def native_designs(request):
    from chembounce import chembounce

    options = copy.deepcopy(request["settings"].get("generator_options", {}))
    forbidden = {"output_dir", "seed", "overall_max_n", "scaffold_db_file", "fingerprint_db_file"}
    if forbidden & options.keys():
        raise ValueError("Native generator resources must come from the frozen request")
    options.update(output_dir=Path(request["output_dir"]) / "designs", seed=request["seed"], overall_max_n=12)
    for asset, option in (("scaffold_db", "scaffold_db_file"), ("fingerprint_db", "fingerprint_db_file")):
        if asset in request["assets"]:
            options[option] = request["assets"][asset]
    frame, cost = chembounce(request["query"]["smiles"], **options)
    supply = list(frame["Final structure"])
    if len(supply) > 12:
        raise ValueError("Native generator exceeded the declared pool")
    write_json(Path(request["output_dir"]) / "design-pool.json", {"designs": supply, "generator_cost": cost})
    return {
        "candidates": [{"smiles": s} for s in supply],
        "evaluations": len(supply),
        "evaluation_unit": "returned native designs",
        "metadata": {"generator_cost": cost},
    }


def native_supply(request, output):
    """Read only a completed current-attempt ledger, never a historical top-three list."""
    arm = request["settings"]["arm"]
    if arm == "synllama":
        return [request["query"]["smiles"]], None, {"raw_count": 1, "policy": "direct original hit"}
    native = Path(request["native"]["directory"])
    manifest = request["native"]["files_sha256"]
    if not manifest or any(
        not (native / name).resolve().is_relative_to(native.resolve()) or digest(native / name) != value
        for name, value in manifest.items()
    ):
        raise ValueError("Native artifacts changed after bounded execution")
    prior = read(native / "request.json")
    if prior["query"] != request["query"] or prior["seed"] != request["seed"]:
        raise ValueError("Native supply query/generation-seed mismatch")
    response = read(native / "response.json")
    if arm == "full":
        supply = read(native / "design-pool.json")["designs"]
        if supply != [r["smiles"] for r in response["candidates"]] or len(supply) > 12:
            raise ValueError("Native pool disagrees with response")
        return supply, None, {"raw_count": len(supply), "policy": "bounded native order"}
    if arm == "baseline":
        return baseline_bridge(native, response, request, output)
    if arm == "crem":
        from benchmarks.external.crem_baseline import select_products

        meta = response["metadata"]
        if meta.get("completed_within_budget") is not True:
            raise ValueError("Incomplete CReM ledger cannot supply designs")
        rows = [json.loads(line) for line in (native / "response.raw.jsonl").read_text().splitlines()]
        if (
            len(rows) > 256
            or [r["index"] for r in rows] != list(range(len(rows)))
            or len(rows) != response["evaluations"]
        ):
            raise ValueError("CReM raw ledger is incomplete")
        selected, decisions = select_products(request["query"]["smiles"], [r["raw_smiles"] for r in rows], 12)
    else:
        from benchmarks.external.mol2mol.selection import select

        raw = read(native / "mol2mol-raw.json")
        from chembounce_synthesis.common import canonical

        if (
            raw["input_smiles"] != canonical(request["query"]["smiles"])
            or raw["seed"] != request["seed"]
            or digest(native / "mol2mol-raw.json") != response["metadata"]["raw_sha256"]
        ):
            raise ValueError("Mol2Mol raw ledger identity mismatch")
        rows = raw["rows"]
        if len(rows) > 256 or len(rows) != response["evaluations"]:
            raise ValueError("Mol2Mol raw ledger is incomplete")
        selected, decisions = select(
            request["query"]["smiles"], rows, maximum=12, deadline=request["deadline_monotonic"]
        )
    if len(rows) > 256:
        raise ValueError("Native raw cap exceeded")
    write_json(output / "native-selection.json", {"raw_count": len(rows), "decisions": decisions})
    return (
        [r["smiles"] for r in selected],
        None,
        {
            "raw_count": len(rows),
            "selection_count": len(selected),
            "policy": "complete fresh raw ledger; native ranking to twelve",
        },
    )


def baseline_bridge(native, response, request, output):
    """Only the historical method's selected primary routes may be reassessed."""
    with (native / "synthesis/results.tsv").open() as stream:
        results = list(csv.DictReader(stream, delimiter="\t"))
    records = read(native / "synthesis/routes.json")["routes"]
    if len(results) > 3 or [r["smiles"] for r in results] != [r["smiles"] for r in response["candidates"]]:
        raise ValueError("Baseline selected products disagree with native response")
    routes = []
    for row in results:
        matches = [r for r in records if r["route_id"] == row["route_id"]]
        if len(matches) != 1:
            raise ValueError("Baseline primary route is missing or ambiguous")
        record = matches[0]
        response_row = response["candidates"][len(routes)]
        if response_row.get("design") != row["design"] or response_row.get("route") != record["route"]:
            raise ValueError("Baseline response differs from selected primary graph")
        if record["design"] != row["design"] or record["verification"]["product"] != row["smiles"]:
            raise ValueError("Baseline selected route/product/design mismatch")
        routes.append({"route_id": row["route_id"], "design": row["design"], "route": record["route"]})
    designs = list(dict.fromkeys(row["design"] for row in results))
    bridge = output / "selected-primary-bridge"
    bridge.mkdir()
    write_json(
        bridge / "manifest.json",
        {"schema_version": 1, "hit": request["query"]["smiles"], "designs": designs, "config": {}},
    )
    write_json(bridge / "routes.json", {"schema_version": 1, "routes": routes})
    write_json(bridge / "empty-proposals.json", {"schema_version": 1, "proposals": []})
    return designs, bridge, {"raw_count": len(results), "policy": "selected primary three only; no rescue"}


def common_config(request, bridge=None):
    from chembounce_synthesis.workflow import load_config
    from chembounce_synthesis.profiles import (
        FIXED_CONFORMERS,
        FIXED_REFINEMENT,
        PROFILE,
        load_profile,
        profile_provenance,
    )

    config = load_config(request["assets"]["synthesis_config"])
    if config.get("recommendation_profile") != PROFILE:
        raise ValueError("Common assessment requires the named reviewed development profile")
    if (
        config.get("minimum_scores") != load_profile(PROFILE)["minimum_scores"]
        or config.get("ranking_policy") != "original-hit-v1"
        or config.get("route_search_policy") != "original-hit-v1"
        or profile_provenance(config)["status"] != "development_protocol_matched"
    ):
        raise ValueError("Common numerical kernel/runtime/minima/objective must match frozen calibration")
    if config["conformers"] != FIXED_CONFORMERS:
        raise ValueError("Common screening must retain fixed seed42/eight-conformer protocol")
    for key, value in FIXED_REFINEMENT.items():
        if key not in {"top_n", "max_total_conformers"} and config["refinement"].get(key) != value:
            raise ValueError("Common final scoring protocol changed")
    expected = {"max_design_attempts": 12, "max_designs": 3, "max_results": 3, "max_scored_products": 12}
    if any(config.get(k) != v for k, v in expected.items()):
        raise ValueError("Common design/product budgets differ from frozen design")
    if config["refinement"]["top_n"] != 6 or config["refinement"]["max_total_conformers"] != 144:
        raise ValueError("Common refinement budget differs from frozen design")
    search = {
        "building_block_top_k": 3,
        "beam_width": 6,
        "max_combinations": 256,
        "max_steps": 12,
        "reaction_product_limit": 64,
    }
    expansion = {"rounds": 1, "max_routes": 12, "building_block_top_k": 6, "min_tanimoto": 0.0}
    if config.get("search") != search or config.get("expansion") != expansion:
        raise ValueError("Common route-search resource settings changed")
    evidence = config["synthesis_evidence"]
    if evidence["profile"] != "computational-route-support-v1" or any(
        evidence.get(key) != value
        for key, value in {
            "max_predictor_calls": 24,
            "max_input_tokens": 12288,
            "max_output_tokens": 18432,
        }.items()
    ):
        raise ValueError("Common strict forward budgets changed")
    if config["provider"].get("kind") != "synllama" or config["provider"].get("samples") != 2:
        raise ValueError("Common provider must declare two SynLlama samples per eligible design")
    if config["provider"].get("max_new_tokens") != 1200:
        raise ValueError("Common SynLlama output cap changed")
    config["provider"]["seed"] = request["seed"]  # generation only; never touch scoring seeds
    arm = request["settings"]["arm"]
    if arm == "synllama":
        config["provider"]["samples"] = 6
    if bridge is not None:
        config["provider"] = {"kind": "saved", "path": str(bridge / "empty-proposals.json")}
        config["carry_forward"] = [str(bridge)]
        config["expansion"]["rounds"] = 0
    elif config.get("carry_forward"):
        raise ValueError("No cross-attempt continuation or cache is allowed")
    ablation = request["settings"].get("ablation")
    if ablation == "no_extra_refill_capacity":
        config["max_design_attempts"] = 3
    elif ablation == "design_only_route_search":
        config["route_search_policy"] = "design-v1"
    elif ablation is not None:
        raise ValueError("Unknown benchmark ablation")
    if ablation and arm != "full":
        raise ValueError("Ablations apply only to the full method")
    return config


def common(request):
    from chembounce_synthesis.workflow import run

    output = Path(request["output_dir"])
    designs, bridge, supply = native_supply(request, output)
    config = common_config(request, bridge)
    write_json(output / "design-supply.json", {"designs": designs, **supply})
    rows = run(
        request["query"]["smiles"],
        designs,
        config,
        output / "synthesis",
        deadline_monotonic=request["deadline_monotonic"] - request["finalization_margin_seconds"],
    )
    summary = read(output / "synthesis/summary.json")
    files = {str(p.relative_to(output)): digest(p) for p in (output / "synthesis").iterdir() if p.is_file()}
    return {
        "candidates": [{"smiles": r["smiles"], "route_id": r["route_id"]} for r in rows],
        "evaluations": summary.get("scored_products", 0),
        "evaluation_unit": "common candidate 3D screens",
        "metadata": {"files_sha256": files, "supply": supply, "completed_monotonic": time.monotonic()},
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", required=True)
    parser.add_argument("--response", required=True)
    parser.add_argument("--phase", choices=("native", "common"), required=True)
    args = parser.parse_args(argv)
    request = read(args.request)
    if request["settings"].get("arm") not in ARMS:
        raise ValueError("Unknown qualified benchmark arm")
    result = native_designs(request) if args.phase == "native" else common(request)
    write_json(args.response, result)


if __name__ == "__main__":
    main()
