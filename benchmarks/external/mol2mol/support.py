"""Small, model-free provenance and protocol helpers."""

import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import sys

UPSTREAM_COMMIT = "eb61443a1a93a8cbd06510548b25cc77cfd9d3a6"
MODEL_SHA256 = "a0f40a2fb4baedaf51ae199bb38a61a8c8d26f03062977eaf06e34079e9f22ca"
MODEL_SHA512 = (
    "6c7a1153aecfeef7fafe0774a66f305f1e1a7bebd4ef8781cd3bd1ef113ff5629"
    "ec502c23384e24b28b4796b4e3449cf90df0135dec18048250cf9c8b1bb7c41"
)
ORIGINAL_FILES = {
    "pyproject.toml": "21e01bf499c0e2923b8ff5c329353e9501ef90425c98405f968d5de04d58b525",
    "requirements-linux-64.lock": "0dae6cef52c3c5018cf1ffa31dbeaf93212b145f3ea912bbd3c90d193028e7a8",
    "priors/sha512.sum": "3132dd89b4a786e6f7a4f4b5e73137ec660c9160e721183d7e39a04739b4689a",
    "LICENSE": "9f0f6abba6fe77a76fdf38bb3fa08904325c0508e888ece3adde95fb64a02bd5",
}
COMMON_RDKIT = "2025.09.6"
SETTINGS = {
    "protocol": "mol2mol-pretrained-v1-proposal",
    "raw_sample_cap": 256,
    "batch_size": 32,
    "temperature": 1.0,
    "strategy": "multinomial",
    "input_randomization": False,
    "device": "cuda:0",
    "filter_margin_seconds": 15.0,
    "batch_reserve_multiplier": 1.25,
}


ADAPTER_FILES = (
    "__init__.py",
    "adapter.py",
    "observation.py",
    "prepare.py",
    "selection.py",
    "support.py",
    "worker.py",
)
COMMON_FILES = (
    "chembounce_synthesis/__init__.py",
    "chembounce_synthesis/common.py",
    "chembounce_synthesis/identity.py",
)


def asset_paths(bundle):
    return {
        "mol2mol_model": bundle["model"],
        **{
            "adapter_source_" + Path(name).stem: str(Path(bundle["adapter_root"]) / name)
            for name in ADAPTER_FILES
        },
        **{
            "common_source_" + Path(name).stem: str(Path(bundle["common_root"]) / name)
            for name in COMMON_FILES
        },
    }


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    def reject(value):
        raise ValueError(f"Nonfinite JSON: {value}")

    return json.loads(Path(path).read_text(), parse_constant=reject)


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def versions():
    return dict(
        sorted(
            (d.metadata["Name"].lower().replace("_", "-"), d.version)
            for d in importlib.metadata.distributions()
        )
    )


def check_hashes(root, entries):
    root = Path(root).resolve()
    for relative, expected in entries.items():
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or digest(path) != expected:
            raise ValueError(f"Pinned file changed: {relative}")


def verify_bundle(bundle, *, sampler=False):
    if bundle["upstream_commit"] != UPSTREAM_COMMIT:
        raise ValueError("Wrong official upstream commit")
    model = Path(bundle["model"])
    if model.stat().st_size != 80198341 or digest(model) != MODEL_SHA256:
        raise ValueError("Official checkpoint size/SHA256 mismatch; refusing pickle load")
    if digest(model, "sha512") != MODEL_SHA512:
        raise ValueError("Official checkpoint SHA512 mismatch; refusing pickle load")
    check_hashes(bundle["upstream"], bundle["source_sha256"])
    check_hashes(bundle["adapter_root"], bundle["adapter_sha256"])
    check_hashes(bundle["common_root"], bundle["common_sha256"])
    if sampler:
        if Path(sys.executable).absolute() != Path(bundle["sampler_python"]).absolute():
            raise ValueError("Wrong sampler interpreter")
        if versions() != bundle["sampler_versions"]:
            raise ValueError("Sampler environment versions changed")
        if digest(sys.executable) != bundle["sampler_python_sha256"]:
            raise ValueError("Sampler Python binary changed")
    else:
        if versions() != bundle["common_versions"]:
            raise ValueError("Common-filter environment versions changed")
        if digest(sys.executable) != bundle["common_python_sha256"]:
            raise ValueError("Common-filter Python binary changed")
        if Path(__file__).resolve().parent != Path(bundle["adapter_root"]).resolve():
            raise ValueError("Adapter imported outside pinned source root")


def upstream_origins(bundle):
    """Every imported upstream Python module must be an unchanged tracked file."""
    root = Path(bundle["upstream"]).resolve()
    result = {}
    for name, module in sorted(sys.modules.items()):
        if name == "reinvent" or name.startswith("reinvent."):
            origin = getattr(module, "__file__", None)
            if not origin:
                # Namespace packages have no source, but all search paths must be pinned-tree paths.
                if not all(Path(p).resolve().is_relative_to(root) for p in module.__path__):
                    raise ValueError(f"Unpinned upstream namespace: {name}")
                continue
            path = Path(origin).resolve()
            if not path.is_relative_to(root):
                raise ValueError(f"Upstream import outside pinned checkout: {name}")
            relative = path.relative_to(root).as_posix()
            actual = digest(path)
            if bundle["source_sha256"].get(relative) != actual:
                raise ValueError(f"Upstream import is not a pinned tracked file: {name}")
            result[name] = {"relative_path": relative, "sha256": actual}
    return result


def validate_request(request):
    allowed = {
        "schema_version",
        "kind",
        "method_id",
        "query",
        "seed",
        "budget",
        "assets",
        "settings",
        "output_dir",
        "deadline_monotonic",
        "finalization_margin_seconds",
    }
    if set(request) - allowed or request.get("kind") != "external" or request.get("schema_version") != 1:
        raise ValueError("Expected whitelisted external public_request schema")
    # target_id is public bookkeeping only. No reference structures or labels are accepted.
    if set(request["query"]) - {"id", "smiles", "target_id"}:
        raise ValueError("Only public query SMILES and identifiers are permitted")
    if request["settings"] != SETTINGS:
        raise ValueError("Settings must exactly match the fixed proposal; no target/check_nll overrides")
    if type(request["seed"]) is not int or request["seed"] not in (42, 43, 44):
        raise ValueError("Expected proposal seed 42, 43 or 44")
    budget = request["budget"]
    if budget != {"wall_seconds": 360, "max_memory_mb": 24000, "max_candidates": 3, "max_evaluations": 256}:
        raise ValueError("Expected fixed 360s/24000MB/3-output/256-raw budget")
    expected_assets = {"mol2mol_bundle", "mol2mol_model"}
    expected_assets.update("adapter_source_" + Path(name).stem for name in ADAPTER_FILES)
    expected_assets.update("common_source_" + Path(name).stem for name in COMMON_FILES)
    if set(request["assets"]) != expected_assets:
        raise ValueError("Only the exact pinned Mol2Mol bundle/model/source assets are accepted")
    if (
        isinstance(request["deadline_monotonic"], bool)
        or not isinstance(request["deadline_monotonic"], (int, float))
        or not math.isfinite(request["deadline_monotonic"])
    ):
        raise ValueError("Finite outer runner deadline required")
    if request["finalization_margin_seconds"] != 2.0:
        raise ValueError("Expected runner finalization margin of 2 seconds")
