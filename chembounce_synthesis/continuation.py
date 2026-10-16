"""Carry reaction evidence forward; never trust previous scores or verification."""

import copy
import json
from pathlib import Path

from .common import canonical, digest


def read_manifest(directory):
    path = Path(directory).resolve() / "manifest.json"
    manifest = json.loads(path.read_text())
    if manifest.get("schema_version") != 1 or not isinstance(manifest.get("config"), dict):
        raise ValueError(f"Invalid continuation manifest: {path}")
    return manifest


def continuation_inputs(directory):
    """Resume using saved model responses without loading or calling the model."""
    directory = Path(directory).resolve()
    manifest = read_manifest(directory)
    config = copy.deepcopy(manifest["config"])
    config["provider"] = {"kind": "saved", "path": str(directory / "proposals.json")}
    # The latest routes already contain the accumulated evidence from older runs.
    config["carry_forward"] = [str(directory)]
    preflight = manifest.get("design_preflight", {})
    designs = (
        preflight.get("saved_design_supply", preflight["requested_designs"])
        if preflight.get("policy", {}).get("version") == "bounded-refill-v1"
        else manifest["designs"]
    )
    return manifest["hit"], designs, config


def load_carried_routes(directories, hit, designs, *, include_excluded=False):
    """Read and deduplicate evidence, leaving chemical verification to workflow.run."""
    candidates = {design: [] for design in designs}
    sources, seen = [], set()
    for directory in directories:
        directory = Path(directory).resolve()
        manifest = read_manifest(directory)
        if canonical(manifest["hit"]) != hit:
            raise ValueError(f"Carried run uses a different original hit: {directory}")
        path = directory / "routes.json"
        data = json.loads(path.read_text())
        if data.get("schema_version") != 1 or not isinstance(data.get("routes"), list):
            raise ValueError(f"Invalid carried routes: {path}")
        source = {
            "directory": str(directory),
            "manifest_sha256": digest(directory / "manifest.json"),
            "routes_sha256": digest(path),
            "loaded_routes": 0,
            "duplicate_routes": 0,
            "excluded_design_routes": 0,
        }
        if include_excluded:
            source["prior_design_preflight"] = manifest.get("design_preflight")
            source["unselected_design_routes_retained"] = 0
        for record in data["routes"]:
            design = canonical(record["design"])
            if design not in candidates:
                if not include_excluded:
                    source["excluded_design_routes"] += 1
                    continue
                candidates[design] = []
            key = (design, json.dumps(record["route"], sort_keys=True))
            if key in seen:
                source["duplicate_routes"] += 1
                continue
            seen.add(key)
            if include_excluded and design not in designs:
                source["unselected_design_routes_retained"] += 1
            candidates[design].append(
                {
                    "route": record["route"],
                    "carried_from": {"run": str(directory), "route_id": record["route_id"]},
                }
            )
            source["loaded_routes"] += 1
        sources.append(source)
    return candidates, sources
