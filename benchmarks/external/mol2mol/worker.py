"""Dedicated official runtime. --inspect imports only; it never loads a model."""

import argparse
import math
import os
from pathlib import Path
import random
import sys
import time

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from benchmarks.external.mol2mol.observation import bind_outputs, observe_decode, validate_tokens
from benchmarks.external.mol2mol.support import digest, read_json, upstream_origins, verify_bundle, write_json


def official_api(bundle):
    verify_bundle(bundle, sampler=True)
    sys.path.insert(0, bundle["upstream"])
    import torch
    import numpy as np
    from reinvent.models.transformer.mol2mol.mol2mol import Mol2MolModel
    from reinvent.models.transformer.core.dataset.dataset import Dataset

    origins = upstream_origins(bundle)
    return torch, np, Mol2MolModel, Dataset, origins


def sample(job, bundle, destination):
    deadline = job["deadline_monotonic"]
    settings = job["settings"]
    if time.monotonic() >= deadline:
        raise TimeoutError("Deadline before official imports")
    torch, np, model_type, dataset_type, origins = official_api(bundle)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    random.seed(job["seed"])
    np.random.seed(job["seed"])
    if not torch.cuda.is_available():
        raise RuntimeError("Declared CUDA sampler unavailable; no silent CPU fallback")
    if torch.cuda.device_count() != 1:
        raise RuntimeError("Expose exactly one authorized GPU for this sampler")
    device = torch.device(settings["device"])
    torch.cuda.set_device(device)
    torch.manual_seed(job["seed"])
    torch.cuda.manual_seed_all(job["seed"])
    torch.cuda.reset_peak_memory_stats(device)
    if time.monotonic() >= deadline:
        raise TimeoutError("Deadline before verified checkpoint load")
    model = model_type.load_from_file(bundle["model"], mode="inference", device=device)
    model.set_temperature(settings["temperature"])
    input_token_capacity = int(model.network.src_embed[1].pe.shape[1])
    validate_tokens(job["smiles"], model.tokenizer, model.vocabulary, input_token_capacity)
    size = settings["batch_size"]
    dataset = dataset_type([job["smiles"]] * size, model.vocabulary, model.tokenizer)
    if len(dataset) != size:
        raise ValueError("Official dataset silently lost an input")
    src, mask = dataset_type.collate_fn([dataset[i] for i in range(size)])
    source_rows = src.tolist()
    src, mask = src.to(device), mask.to(device)
    rows, batches, longest = [], [], 0.0
    result = {
        "schema_version": 1,
        "rows": rows,
        "batches": batches,
        "pid": os.getpid(),
        "process_group_id": os.getpgrp(),
        "seed": job["seed"],
        "input_smiles": job["smiles"],
        "max_sequence_length": model.max_sequence_length,
        "input_token_capacity": input_token_capacity,
        "termination_semantics": "BOS plus at most max_sequence_length-1 sampled tokens; EOS then zero padding",
        "upstream_imports": origins,
        "stop_reason": "raw_cap",
    }
    for offset in range(0, settings["raw_sample_cap"], size):
        reserve = settings["filter_margin_seconds"] + settings["batch_reserve_multiplier"] * longest
        if time.monotonic() + reserve >= deadline:
            result["stop_reason"] = "early_finish_reserve"
            break
        start = time.monotonic()
        with observe_decode(model.vocabulary) as captured:
            inputs, outputs, nlls = model.sample(src, mask, settings["strategy"])
        if any(text != job["smiles"] for text in inputs):
            raise ValueError("Official input decoding changed canonical public input")
        batch = bind_outputs(
            captured,
            source_rows,
            inputs,
            outputs,
            nlls,
            eos=model.vocabulary["$"],
            bos=model.vocabulary["^"],
            max_length=model.max_sequence_length,
            tokenizer=model.tokenizer,
            vocabulary=model.vocabulary,
        )
        for row in batch:
            if not math.isfinite(row["native_nll"]):
                row["native_nll"] = str(row["native_nll"])
            row["batch_index"] = len(batches)
            rows.append(row)
        duration = time.monotonic() - start
        longest = max(longest, duration)
        batches.append({"offset": offset, "count": size, "wall_seconds": duration})
        # Local diagnostics persist even if a later batch/outer run fails.
        write_json(destination, result)
        if time.monotonic() >= deadline:
            raise TimeoutError("Batch crossed outer deadline; no accepted partial response")
    result["upstream_imports"] = upstream_origins(bundle)
    result["cuda"] = {
        "device": str(device),
        "name": torch.cuda.get_device_name(device),
        "torch_version": torch.__version__,
        "compiled_cuda": torch.version.cuda,
        "visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "peak_allocated_bytes": torch.cuda.max_memory_allocated(device),
        "peak_reserved_bytes": torch.cuda.max_memory_reserved(device),
        "scope": "PyTorch allocator in sampler only; not full device/process memory",
    }
    write_json(destination, result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True)
    parser.add_argument("--inspect", action="store_true")
    parser.add_argument("--job")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    bundle = read_json(args.bundle)
    if args.inspect:
        torch, np, _, _, origins = official_api(bundle)
        write_json(
            args.output,
            {
                "model_loaded": False,
                "inference_executed": False,
                "torch": torch.__version__,
                "numpy": np.__version__,
                "compiled_cuda": torch.version.cuda,
                "upstream_imports": origins,
                "runtime_imports": {
                    name: {
                        "path": str(Path(sys.modules[name].__file__).resolve()),
                        "sha256": digest(sys.modules[name].__file__),
                    }
                    for name in ("torch", "torch._C", "numpy", "rdkit", "rdkit.Chem.rdchem")
                    if name in sys.modules
                },
            },
        )
    else:
        if not args.job:
            parser.error("--job is required for authorized future sampling")
        sample(read_json(args.job), bundle, args.output)


if __name__ == "__main__":
    main()
