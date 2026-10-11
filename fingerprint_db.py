"""Validated Morgan-bit databases and bounded-memory deterministic search."""

import json
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator

POPCOUNT = np.array([i.bit_count() for i in range(256)], dtype=np.uint8)


@lru_cache(maxsize=32)
def _generator(radius, fp_size):
    return rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size)


def generate_query_fingerprint(mol, radius=2, fp_size=2048):
    return _generator(radius, fp_size).GetFingerprintAsNumPy(mol)


def load_fingerprint_db(fp_file=None, trust_legacy=False):
    if fp_file is None:
        root = Path(__file__).parent / "data"
        preferred = root / "fingerprints_mw250_v2"
        path = preferred if preferred.is_dir() else root / "scaffold_fingerprints_mw250.npz"
    else:
        path = Path(fp_file)
    if path.is_dir():
        metadata = json.loads((path / "metadata.json").read_text())
        data = {
            **metadata,
            "fingerprints": np.load(path / "fingerprints.npy", mmap_mode="r"),
            "smiles": np.load(path / "smiles.npy", mmap_mode="r"),
        }
    else:
        with np.load(path, allow_pickle=trust_legacy) as archive:
            try:
                data = {key: archive[key] for key in archive.files}
            except ValueError as exc:
                if "Object arrays cannot be loaded" not in str(exc):
                    raise
                raise ValueError(
                    "Legacy object-array database requires explicit trust: use --trust-legacy-db only for trusted files, then regenerate in the 2.0 format"
                ) from exc
    if int(data.get("schema_version", 1)) not in (1, 2):
        raise ValueError("Unsupported fingerprint database schema")
    if str(data.get("bitorder", "big")) != "big" or bool(data.get("include_chirality", False)):
        raise ValueError("Unsupported fingerprint bit order or chirality setting")
    for key in ("fingerprints", "smiles", "fp_size", "radius"):
        if key not in data:
            raise ValueError(f"Fingerprint database missing {key}")
    for key in ("fp_size", "radius"):
        value = np.asarray(data[key])
        if value.ndim != 0 or value.dtype.kind not in "iu":
            raise ValueError(f"{key} must be an integer scalar")
    data["fp_size"], data["radius"] = int(data["fp_size"]), int(data["radius"])
    data["packed"] = bool(data.get("packed", False))
    if data["fp_size"] < 1 or data["radius"] < 0:
        raise ValueError("Invalid fingerprint parameters")
    fps, smiles = data["fingerprints"], data["smiles"]
    width = (data["fp_size"] + 7) // 8 if data["packed"] else data["fp_size"]
    if fps.ndim != 2 or fps.shape != (len(smiles), width) or fps.dtype != np.uint8:
        raise ValueError("Fingerprint dimensions or dtype do not match metadata")
    if np.ndim(smiles) != 1 or any(not isinstance(s, str) or not s for s in smiles):
        raise ValueError("Database SMILES must be nonempty strings")
    for start in range(0, len(fps), 8192):
        chunk = fps[start : start + 8192]
        if not data["packed"] and np.any(chunk > 1):
            raise ValueError("Unpacked fingerprints must contain only 0 and 1")
        remainder = data["fp_size"] % 8
        if data["packed"] and remainder and np.any(chunk[:, -1] & ((1 << (8 - remainder)) - 1)):
            raise ValueError("Nonzero padding bits in packed database")
    data["path"] = str(path.resolve())
    return data


def _rank(indices, scores, top_n):
    if top_n is not None and len(indices) > top_n:
        # Retain the lowest row indices when the cutoff has tied scores.
        cutoff = np.partition(scores, len(scores) - top_n)[len(scores) - top_n]
        better = np.flatnonzero(scores > cutoff)
        tied = np.flatnonzero(scores == cutoff)
        tied = tied[np.argsort(indices[tied], kind="stable")][: top_n - len(better)]
        selected = np.concatenate((better, tied))
        indices, scores = indices[selected], scores[selected]
    order = np.lexsort((indices, -scores))
    return indices[order], scores[order]


def calculate_bulk_tanimoto(
    query_fp, db_fps, top_n=None, threshold=0.0, packed=False, chunk_size=8192, exclude=None
):
    if not 0 <= threshold <= 1 or chunk_size < 1 or (top_n is not None and top_n < 1):
        raise ValueError("Invalid similarity threshold or search limit")
    query = np.asarray(query_fp, dtype=np.uint8).reshape(-1)
    if db_fps.ndim != 2 or db_fps.shape[1] != len(query):
        raise ValueError("Query and database fingerprint widths differ")
    indices, scores = np.empty(0, dtype=np.int64), np.empty(0, dtype=float)
    index_chunks, score_chunks = [], []
    query_sum = POPCOUNT[query].sum() if packed else query.sum()
    for start in range(0, len(db_fps), chunk_size):
        chunk = db_fps[start : start + chunk_size]
        if packed:
            intersection = POPCOUNT[np.bitwise_and(chunk, query)].sum(axis=1)
            counts = POPCOUNT[chunk].sum(axis=1)
        else:
            intersection = np.bitwise_and(chunk, query).sum(axis=1)
            counts = chunk.sum(axis=1)
        union = counts + query_sum - intersection
        # Match RDKit binary Tanimoto: two empty bit vectors score zero.
        similarities = np.divide(intersection, union, out=np.zeros(len(chunk)), where=union != 0)
        valid = similarities >= threshold
        if exclude is not None:
            valid &= ~exclude[start : start + len(chunk)]
        local = np.flatnonzero(valid)
        if top_n is None:
            index_chunks.append(local + start)
            score_chunks.append(similarities[local])
        else:
            indices = np.concatenate((indices, local + start))
            scores = np.concatenate((scores, similarities[local]))
            indices, scores = _rank(indices, scores, top_n)
    if top_n is None and index_chunks:
        indices, scores = np.concatenate(index_chunks), np.concatenate(score_chunks)
    return _rank(indices, scores, top_n)


def search_similar_scaffolds_fast(
    query_mol, fp_db, scaffold_top_n=None, threshold=0.3, use_multiprocessing=False, n_jobs=None
):
    query = generate_query_fingerprint(query_mol, fp_db["radius"], fp_db["fp_size"])
    packed = fp_db.get("packed", False)
    if packed:
        query = np.packbits(query)
    smiles = np.asarray(fp_db["smiles"])
    exclude = smiles == Chem.MolToSmiles(query_mol)
    # Duplicate rows must not consume top-k slots. New databases are canonical.
    if "_duplicate_mask" not in fp_db:
        if bool(fp_db.get("canonical_sorted_unique", False)):
            duplicate = np.zeros(len(smiles), dtype=bool)
            if len(smiles) > 1 and np.any(smiles[:-1] >= smiles[1:]):
                raise ValueError("Database claims sorted unique SMILES but rows disagree")
        else:
            _, first = np.unique(smiles, return_index=True)
            duplicate = np.ones(len(smiles), dtype=bool)
            duplicate[first] = False
        fp_db["_duplicate_mask"] = duplicate
    duplicate = fp_db["_duplicate_mask"]
    canonical_db = bool(fp_db.get("canonical_sorted_unique", False))
    fetch_n = scaffold_top_n
    while True:
        indices, scores = calculate_bulk_tanimoto(
            query,
            fp_db["fingerprints"],
            top_n=fetch_n,
            threshold=threshold,
            packed=packed,
            exclude=exclude | duplicate,
        )
        if canonical_db:
            return pd.Series(scores, index=smiles[indices], name="Tanimoto Similarity", dtype=float)
        # Legacy libraries contain pseudoatoms, disconnected entries and equivalent
        # SMILES. Validate selected hits and refill top-k after those exclusions.
        valid = {}
        rejected = fp_db.setdefault("_legacy_rejected_rows", set())
        query_smiles = Chem.MolToSmiles(query_mol)
        for index, score in zip(indices, scores):
            with rdBase.BlockLogs():
                mol = Chem.MolFromSmiles(str(smiles[index]))
            if (
                mol is None
                or not mol.GetNumAtoms()
                or len(Chem.GetMolFrags(mol)) != 1
                or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())
            ):
                rejected.add(int(index))
                continue
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
            canonical = Chem.MolToSmiles(mol)
            if canonical != query_smiles and canonical not in valid:
                valid[canonical] = score
        if (
            scaffold_top_n is None
            or len(valid) >= scaffold_top_n
            or len(indices) < fetch_n
            or fetch_n >= len(smiles)
        ):
            rows = list(valid.items())[:scaffold_top_n]
            return pd.Series(dict(rows), name="Tanimoto Similarity", dtype=float)
        fetch_n = min(len(smiles), max(fetch_n * 2, 1))


def calculate_bulk_tanimoto_multiprocessing(query_fp, db_fps, top_n=None, threshold=0.0, n_jobs=None):
    """Compatibility entry point; chunking avoids copying the DB into workers."""
    if n_jobs is not None and n_jobs != -1 and n_jobs < 1:
        raise ValueError("n_jobs must be positive or -1")
    return calculate_bulk_tanimoto(query_fp, db_fps, top_n, threshold)


def search_similar_scaffolds_fast_mp(*args, **kwargs):
    return search_similar_scaffolds_fast(*args, **kwargs)
