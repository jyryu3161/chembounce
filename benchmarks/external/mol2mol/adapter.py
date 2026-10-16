"""External runner adapter: common identity process plus inherited-group sampler."""

import argparse
from collections import Counter
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from benchmarks.external.mol2mol.selection import select
from benchmarks.external.mol2mol.support import (
    COMMON_RDKIT,
    asset_paths,
    digest,
    read_json,
    validate_request,
    verify_bundle,
    write_json,
)


OWNED_FILES = (
    "mol2mol-job.json",
    "mol2mol-raw.json",
    "mol2mol-selection.json",
    "sampler.stdout.log",
    "sampler.stderr.log",
    "response.json",
    "mol2mol-response-provenance.json",
)


def fresh_output(output_dir, response_path):
    output = Path(output_dir).resolve()
    response = Path(response_path)
    if not output.is_dir() or response.parent.resolve() != output or response.name != "response.json":
        raise ValueError("Expected runner-created attempt directory and response.json")
    for name in OWNED_FILES:
        path = output / name
        if path.exists() or path.is_symlink() or any(output.glob(name + ".*.tmp")):
            raise FileExistsError(f"Refusing to overwrite existing attempt evidence: {name}")
    return output


def run_child(command, *, deadline, log_dir, env=None):
    """Inherit the runner's owned process group and deadline, including grandchildren.

    Never use start_new_session here: bounded_method owns group-wide hard cleanup.
    SIGTERM propagates as failure, never as a successful partial result.
    """
    child = None
    previous = {}

    def interrupted(signum, frame):
        raise InterruptedError(f"Sampler adapter received signal {signum}")

    for signum in (signal.SIGTERM, signal.SIGINT):
        previous[signum] = signal.signal(signum, interrupted)
    try:
        with (
            (Path(log_dir) / "sampler.stdout.log").open("w") as stdout,
            (Path(log_dir) / "sampler.stderr.log").open("w") as stderr,
        ):
            child = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=env)
            while child.poll() is None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError("Sampler exceeded shared outer deadline")
                time.sleep(min(0.05, remaining))
            if child.returncode:
                raise RuntimeError(f"Official sampler exited {child.returncode}")
            if time.monotonic() >= deadline:
                raise TimeoutError("Sampler completed after shared deadline")
    finally:
        if child is not None and child.poll() is None:
            child.terminate()
            try:
                child.wait(timeout=0.2)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait(timeout=0.2)
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def run(request, response_path):
    validate_request(request)
    if os.getpgrp() != os.getpid():
        raise RuntimeError("Run through bounded_method: adapter must lead the owned process group")
    output = fresh_output(request["output_dir"], response_path)
    bundle_path = request["assets"]["mol2mol_bundle"]
    bundle = read_json(bundle_path)
    verify_bundle(bundle)
    if any(
        Path(request["assets"][name]).resolve() != Path(path).resolve()
        for name, path in asset_paths(bundle).items()
    ):
        raise ValueError("Protocol asset paths disagree with preparation bundle")
    from rdkit import rdBase
    from chembounce_synthesis.identity import chemical_identity, IDENTITY_POLICY

    if rdBase.rdkitVersion != COMMON_RDKIT:
        raise ValueError(f"Common identity requires RDKit {COMMON_RDKIT}")
    import chembounce_synthesis.identity as identity_module
    import chembounce_synthesis.common as common_module

    for module in (identity_module, common_module):
        if not Path(module.__file__).resolve().is_relative_to(Path(bundle["common_root"]).resolve()):
            raise ValueError("Common chemistry imported outside pinned checkout")
    reference = chemical_identity(request["query"]["smiles"])
    deadline = request["deadline_monotonic"] - request["finalization_margin_seconds"]
    if time.monotonic() >= deadline:
        raise TimeoutError("Deadline before sampling preparation")
    job = {
        "smiles": reference.exact_smiles,
        "seed": request["seed"],
        "settings": request["settings"],
        "deadline_monotonic": deadline,
    }
    job_path, raw_path = output / "mol2mol-job.json", output / "mol2mol-raw.json"
    write_json(job_path, job)
    worker = Path(__file__).with_name("worker.py")
    command = [
        bundle["sampler_python"],
        str(worker),
        "--bundle",
        bundle_path,
        "--job",
        str(job_path),
        "--output",
        str(raw_path),
    ]
    env = dict(
        os.environ,
        PYTHONNOUSERSITE="1",
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        PYTHONHASHSEED=str(request["seed"]),
    )
    # No child-specific CUDA reassignment: activation chooses one visible GPU outside the adapter.
    run_child(command, deadline=deadline, log_dir=output, env=env)
    raw = read_json(raw_path)
    if raw["process_group_id"] != os.getpgrp():
        raise ValueError("Official sampler escaped runner process group")
    if raw["input_smiles"] != reference.exact_smiles or raw["seed"] != request["seed"]:
        raise ValueError("Raw result does not belong to this query/seed")
    if len(raw["rows"]) > 256:
        raise ValueError("Raw cap exceeded")
    candidates, annotated = select(reference.exact_smiles, raw["rows"], deadline=deadline)
    selection_path = output / "mol2mol-selection.json"
    write_json(selection_path, annotated)
    if time.monotonic() >= deadline:
        raise TimeoutError("Deadline after common selection")
    response = {
        "candidates": candidates,
        "evaluations": len(raw["rows"]),
        "evaluation_unit": "official Mol2Mol raw sampled sequences",
        "metadata": {
            "baseline": "plain pretrained Mol2Mol plus declared selection",
            "bundle_sha256": digest(bundle_path),
            "raw_sha256": digest(raw_path),
            "selection_sha256": digest(selection_path),
            "selection_counts": dict(Counter(row["selection_status"] for row in annotated)),
            "identity_policy": IDENTITY_POLICY,
            "common_rdkit": rdBase.rdkitVersion,
            "sampler_versions": bundle["sampler_versions"],
            "cuda": raw["cuda"],
            "stop_reason": raw["stop_reason"],
            "batches": len(raw["batches"]),
            "max_sequence_length": raw["max_sequence_length"],
            "input_token_capacity": raw["input_token_capacity"],
            "training_inventory": "unknown",
            "route_coverage": "no supplied route",
            "adapter_process_group_id": os.getpgrp(),
            "sampler_command": command,
        },
    }
    write_json(response_path, response)
    write_json(
        output / "mol2mol-response-provenance.json",
        {
            "response_sha256": digest(response_path),
            "request_job_sha256": digest(job_path),
            "bundle_sha256": digest(bundle_path),
            "raw_sha256": digest(raw_path),
            "selection_sha256": digest(selection_path),
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", required=True)
    parser.add_argument("--response", required=True)
    args = parser.parse_args()
    run(read_json(args.request), args.response)


if __name__ == "__main__":
    main()
