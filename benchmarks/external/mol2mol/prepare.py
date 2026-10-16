"""Verify a prepared isolated official environment; never execute model inference.

Download/install commands are documented in README.md. This tool records a new
bundle only after tracked-source, checkpoint, dependency and import checks pass.
"""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from benchmarks.external.mol2mol.support import (
    COMMON_RDKIT,
    COMMON_FILES,
    ADAPTER_FILES,
    asset_paths,
    MODEL_SHA256,
    MODEL_SHA512,
    ORIGINAL_FILES,
    SETTINGS,
    UPSTREAM_COMMIT,
    check_hashes,
    digest,
    read_json,
    versions,
    write_json,
)


def command(argv, *, cwd=None):
    return subprocess.check_output(argv, cwd=cwd, text=True).strip()


def source_inventory(upstream):
    if command(["git", "rev-parse", "HEAD"], cwd=upstream) != UPSTREAM_COMMIT:
        raise ValueError("Wrong official checkout commit")
    check_hashes(upstream, ORIGINAL_FILES)
    inventory = {}
    # Sparse prior exclusion is allowed; every Python source in the Git tree must be present.
    records = command(["git", "ls-tree", "-r", UPSTREAM_COMMIT], cwd=upstream)
    for record in records.splitlines():
        details, relative = record.split("\t", 1)
        _, kind, object_id = details.split()
        if kind != "blob" or not (relative.endswith(".py") or relative in ORIGINAL_FILES):
            continue
        path = upstream / relative
        # Compare exact bytes to their Git blob, not merely a mutable working-tree status.
        actual_object = command(["git", "hash-object", "--no-filters", str(path)], cwd=upstream)
        if actual_object != object_id:
            raise ValueError(f"Official tracked source differs: {relative}")
        inventory[relative] = digest(path)
    return inventory


def prepare(upstream, sampler_python, model, output, installation_report):
    import rdkit
    from rdkit import rdBase
    from chembounce_synthesis import identity

    if rdBase.rdkitVersion != COMMON_RDKIT:
        raise ValueError(f"Run preparation in the common RDKit {COMMON_RDKIT} environment")
    if output.exists():
        raise ValueError("Use a new preparation directory; existing evidence is immutable")
    output.mkdir(parents=True)
    inventory = source_inventory(upstream)
    if (
        model.stat().st_size != 80198341
        or digest(model) != MODEL_SHA256
        or digest(model, "sha512") != MODEL_SHA512
    ):
        raise ValueError("Official model size/checksum mismatch")
    checksum_lines = (upstream / "priors/sha512.sum").read_text().splitlines()
    if not any(line.split() == [MODEL_SHA512, "mol2mol_medium_similarity.prior"] for line in checksum_lines):
        raise ValueError("Model not bound to official sha512.sum entry")
    # Query only environment metadata, never torch/model loading.
    code = "import importlib.metadata as m,json; print(json.dumps(dict(sorted((d.metadata['Name'].lower().replace('_','-'),d.version) for d in m.distributions()))))"
    sampler_versions = json.loads(command([str(sampler_python), "-c", code]))
    expected = {
        "torch": "2.3.1+cu121",
        "torchvision": "0.18.1+cu121",
        "rdkit": "2023.9.5",
        "reinvent": "4.4.22",
    }
    for package, version in expected.items():
        if sampler_versions.get(package) != version:
            raise ValueError(f"Sampler package {package} must be {version}")
    pip_check = command([str(sampler_python), "-m", "pip", "check"])
    (output / "pip-check.txt").write_text(pip_check + "\n")
    freeze = command([str(sampler_python), "-m", "pip", "freeze", "--all"])
    (output / "resolved-environment.txt").write_text(freeze + "\n")
    # Validate all declared pyproject requirements, including exact torch pins.
    check = """import importlib.metadata as m,json,sys,tomli
from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
p=tomli.load(open(sys.argv[1],'rb'))['project']
assert SpecifierSet(p['requires-python']).contains('.'.join(map(str,sys.version_info[:3])))
report=[]
for raw in p['dependencies']:
 r=Requirement(raw)
 if r.marker and not r.marker.evaluate(): continue
 actual=m.version(r.name)
 assert r.specifier.contains(actual), (raw,actual)
 report.append({'requirement':raw,'installed':actual})
print(json.dumps(report))"""
    requirements = json.loads(command([str(sampler_python), "-c", check, str(upstream / "pyproject.toml")]))
    root = Path(__file__).resolve().parents[3]
    adapter = Path(__file__).resolve().parent
    report = read_json(installation_report)
    installed = {
        entry["metadata"]["name"].lower().replace("_", "-"): entry["metadata"]["version"]
        for entry in report["install"]
    }
    if installed.get("reinvent") != "4.4.22":
        raise ValueError("Installation report must include official reinvent4.4.22")
    for name, version in installed.items():
        if sampler_versions.get(name) != version:
            raise ValueError(f"Installation report/version mismatch: {name}")
    bundle = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "upstream_commit": UPSTREAM_COMMIT,
        "upstream": str(upstream),
        "source_sha256": inventory,
        "model": str(model),
        "model_sha256": MODEL_SHA256,
        "model_sha512": MODEL_SHA512,
        "model_bytes": model.stat().st_size,
        "sampler_python": str(sampler_python),
        "sampler_python_sha256": digest(sampler_python),
        "sampler_versions": sampler_versions,
        "common_python": sys.executable,
        "common_python_sha256": digest(sys.executable),
        "common_versions": versions(),
        "common_rdkit": rdBase.rdkitVersion,
        "common_import_origins": {
            "rdkit": str(Path(rdkit.__file__).resolve()),
            "identity": str(Path(identity.__file__).resolve()),
        },
        "adapter_root": str(adapter),
        "adapter_sha256": {name: digest(adapter / name) for name in ADAPTER_FILES},
        "common_root": str(root),
        "common_sha256": {p: digest(root / p) for p in COMMON_FILES},
        "pyproject_requirements_validation": requirements,
        "upstream_original_files_sha256": ORIGINAL_FILES,
        "resolved_environment_sha256": digest(output / "resolved-environment.txt"),
        "pip_check_sha256": digest(output / "pip-check.txt"),
        "installation_report_sha256": digest(installation_report),
        "dependency_scope": "pyproject-compliant new resolution; NOT reproduction of contradictory upstream lock",
        "official_code_license": "Apache-2.0",
        "official_checkpoint_license": "Selected Git asset is distributed under same-commit Apache-2.0 LICENSE; registered Zenodo priors also declare Apache-2.0",
        "zenodo_byte_equivalence_to_git_checkpoint": "not established",
        "license_metadata": {
            "doi": "10.5281/zenodo.15641297",
            "api": "https://api.datacite.org/dois/10.5281/zenodo.15641297",
            "parent_saved_response_sha256": "da379ed3254a24d987feb4e83d9c03b0e9c16e2de062b4b849a1de8d78cd3e9f",
            "concept_doi": "10.5281/zenodo.15641296",
            "rights": "Apache-2.0; Copyright 2025 AstraZeneca",
        },
        "checkpoint_source": f"https://github.com/MolecularAI/REINVENT4/blob/{UPSTREAM_COMMIT}/priors/mol2mol_medium_similarity.prior",
        "checkpoint_release": "https://zenodo.org/records/15641296",
        "inference_executed": False,
        "actual_observer_equivalence": "pending separate authorized preflight",
    }
    bundle_path = output / "bundle.json"
    write_json(bundle_path, bundle)
    # CUDA is hidden for import-only preparation; official trusted pickle is NOT loaded.
    import os

    subprocess.run(
        [
            str(sampler_python),
            str(adapter / "worker.py"),
            "--bundle",
            str(bundle_path),
            "--inspect",
            "--output",
            str(output / "import-inspection.json"),
        ],
        check=True,
        env=dict(os.environ, CUDA_VISIBLE_DEVICES="", PYTHONNOUSERSITE="1", OMP_NUM_THREADS="1"),
    )
    imports = read_json(output / "import-inspection.json")
    method = {
        "id": "mol2mol_pretrained",
        "kind": "external",
        "command": [
            sys.executable,
            str(adapter / "adapter.py"),
            "--request",
            "{request}",
            "--response",
            "{response}",
        ],
        "settings": SETTINGS,
        "budget_overrides": {"max_evaluations": 256},
        "assets": {
            "mol2mol_bundle": {"artifact": str(bundle_path), "sha256": digest(bundle_path)},
            **{
                name: {"artifact": path, "sha256": digest(path)} for name, path in asset_paths(bundle).items()
            },
        },
    }
    write_json(output / "method-fragment.json", method)
    write_json(
        output / "validation.json",
        {
            "status": "prepared_not_approved_for_inference",
            "bundle_sha256": digest(bundle_path),
            "method_fragment_sha256": digest(output / "method-fragment.json"),
            "import_inspection_sha256": digest(output / "import-inspection.json"),
            "upstream_imported_modules": len(imports["upstream_imports"]),
            "model_loaded": False,
            "inference_executed": False,
            "pyproject_requirements_checked": len(requirements),
            "pending": [
                "parent review",
                "official-checkpoint observer ON/OFF equality preflight",
                "fixed GPU resource and process-group preflight",
                "final protocol asset pinning",
            ],
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("upstream", "sampler-python", "model", "output", "installation-report"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    prepare(
        args.upstream.resolve(),
        args.sampler_python.absolute(),
        args.model.resolve(),
        args.output.resolve(),
        args.installation_report.resolve(),
    )


if __name__ == "__main__":
    main()
