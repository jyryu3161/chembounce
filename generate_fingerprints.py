#!/usr/bin/env python3
"""Build pickle-free packed Morgan databases (NPZ or memory-mappable directory)."""

import argparse
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np
from rdkit import Chem, rdBase
from rdkit.SimDivFilters import MaxMinPicker
from rdkit import DataStructs

from fingerprint_db import generate_query_fingerprint, load_fingerprint_db


def sample_by_diversity(fps, percentages, seed=42):
    if any(not 0 <= p <= 100 for p in percentages):
        raise ValueError("Percentages must lie in [0,100]")
    if not fps or not percentages:
        return {p: [] for p in percentages}
    count = int(len(fps) * max(percentages) / 100)
    picks = (
        list(
            MaxMinPicker().LazyPick(
                lambda i, j: 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j]), len(fps), count, seed=seed
            )
        )
        if count
        else []
    )
    return {p: picks[: int(len(fps) * p / 100)] for p in percentages}


def generate_fingerprints(input_file, output_file, fp_size=2048, radius=2, batch_size=8192):
    if fp_size < 1 or radius < 0 or batch_size < 1:
        raise ValueError("Fingerprint size/batch size must be positive; radius must be nonnegative")
    source, destination = Path(input_file), Path(output_file)
    if destination.exists():
        raise FileExistsError(f"Database already exists: {destination}; choose a new output path")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open() as stream:
        n_lines = sum(1 for _ in stream)
    if n_lines == 0:
        raise ValueError("Scaffold input is empty")
    digest = hashlib.sha256()
    with source.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".fingerprints-") as scratch:
        scratch = Path(scratch)
        width = (fp_size + 7) // 8
        raw = np.lib.format.open_memmap(
            scratch / "raw.npy", mode="w+", dtype=np.uint8, shape=(n_lines, width)
        )
        rows = {}
        rejected = 0
        duplicates = 0
        blanks = 0
        rejection_counts = {}
        rejection_examples = []
        with source.open() as stream:
            for line_number, line in enumerate(stream, 1):
                smiles = line.strip()
                if not smiles:
                    blanks += 1
                    continue
                with rdBase.BlockLogs():
                    mol = Chem.MolFromSmiles(smiles)
                reason = None
                if mol is None or not mol.GetNumAtoms():
                    reason = "invalid_smiles"
                elif len(Chem.GetMolFrags(mol)) != 1:
                    reason = "multiple_components"
                elif any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
                    reason = "dummy_atoms"
                if reason:
                    rejected += 1
                    rejection_counts[reason] = rejection_counts.get(reason, 0) + 1
                    if len(rejection_examples) < 10:
                        rejection_examples.append({"line": line_number, "smiles": smiles, "reason": reason})
                    continue
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(0)
                smiles = Chem.MolToSmiles(mol)
                if smiles in rows:
                    duplicates += 1
                    continue
                index = len(rows)
                raw[index] = np.packbits(generate_query_fingerprint(mol, radius, fp_size))
                rows[smiles] = index
                if len(rows) % 100000 == 0:
                    print(
                        f"Processed {line_number}/{n_lines} lines; {len(rows)} unique scaffolds", flush=True
                    )
        if not rows:
            raise ValueError("No valid connected scaffolds in input")
        smiles = np.asarray(sorted(rows), dtype=str)
        order = np.asarray([rows[s] for s in smiles])
        stage = scratch / "database"
        stage.mkdir()
        fps = np.lib.format.open_memmap(
            stage / "fingerprints.npy", mode="w+", dtype=np.uint8, shape=(len(rows), width)
        )
        for start in range(0, len(rows), batch_size):
            fps[start : start + batch_size] = raw[order[start : start + batch_size]]
        fps.flush()
        np.save(stage / "smiles.npy", smiles, allow_pickle=False)
        metadata = {
            "schema_version": 2,
            "packed": True,
            "bitorder": "big",
            "fp_size": fp_size,
            "radius": radius,
            "include_chirality": False,
            "canonical_sorted_unique": True,
            "source_sha256": digest.hexdigest(),
            "rdkit_version": rdBase.rdkitVersion,
            "input_lines": n_lines,
            "rejected_lines": rejected,
            "duplicate_lines": duplicates,
            "blank_lines": blanks,
            "rejection_counts_json": json.dumps(rejection_counts, sort_keys=True),
            "rejection_examples_json": json.dumps(rejection_examples),
            "rows": len(rows),
        }
        (stage / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        if destination.suffix == ".npz":
            stage_npz = scratch / "database.npz"
            np.savez_compressed(stage_npz, fingerprints=fps, smiles=smiles, **metadata)
            stage_npz.replace(destination)
        else:
            stage.replace(destination)
    print(f"Generated {len(rows)} unique scaffolds; rejected {rejected} lines; saved {destination}")
    return len(rows)


def verify_fingerprints(fp_file, sample_size=100, seed=42, trust_legacy=False):
    if sample_size < 1:
        raise ValueError("sample_size must be positive")
    data = load_fingerprint_db(fp_file, trust_legacy)
    count = len(data["smiles"])
    indices = np.random.default_rng(seed).choice(count, min(sample_size, count), replace=False)
    for index in indices:
        mol = Chem.MolFromSmiles(data["smiles"][index])
        if mol is None:
            raise ValueError(f"Invalid database SMILES at row {index}")
        expected = generate_query_fingerprint(mol, data["radius"], data["fp_size"])
        if data["packed"]:
            expected = np.packbits(expected)
        if not np.array_equal(expected, data["fingerprints"][index]):
            raise ValueError(f"Fingerprint does not match SMILES at row {index}")
    return len(indices)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", default="data/Scaffolds_processed.txt")
    parser.add_argument("-o", "--output", default="data/scaffold_fingerprints.npz")
    parser.add_argument("--fp-size", type=int, default=2048)
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--batch-size", type=int, default=8192)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    try:
        generate_fingerprints(args.input, args.output, args.fp_size, args.radius, args.batch_size)
        if args.verify:
            print(f"Verified {verify_fingerprints(args.output)} rows against RDKit")
    except (ValueError, OSError) as exc:
        parser.exit(2, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
