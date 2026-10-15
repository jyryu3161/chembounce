"""Score a fixed same-assay control panel with the existing common scorer.

Selection and labels come exclusively from the validated panel. This evaluator
does not select controls, tune thresholds, or assess generated compounds.
"""

import argparse
import concurrent.futures
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import platform
import sys
import time

import numpy as np
from rdkit import rdBase
from rdkit.Chem import RDConfig
import scipy

from benchmarks.research.assay_matched import (
    aggregate_scores,
    control_pairs,
    public_queries,
    validate_panel,
)
from benchmarks.research.controls import similarity_scores
from benchmarks.research.protocol import digest
from chembounce_synthesis.conformers import ConformerSettings

METRICS = ("bit_tanimoto", "count_tanimoto", "hit_shape", "hit_pharmacophore", "score")
THREE_D_METRICS = METRICS[2:]
CODE_FILES = (
    "benchmarks/research/assay_control_scoring.py",
    "benchmarks/research/assay_matched.py",
    "benchmarks/research/controls.py",
    "benchmarks/research/evaluation.py",
    "benchmarks/research/protocol.py",
    "benchmarks/research/assay_metadata.py",
    "benchmarks/research/import_chembl.py",
    "chembounce_synthesis/conformers.py",
    "chembounce_synthesis/identity.py",
    "chembounce_synthesis/common.py",
    "chembounce_synthesis/exploration.py",
    "chemistry.py",
)


def _save(path, data):
    """Replace each snapshot atomically; a new run owns its output directory."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def _json_hash(data):
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def _unreported_trial(seed, reason):
    # An unexpected exception can occur after some work. Do not claim that a
    # seed failed chemically, or was never attempted, without a returned record.
    return {"seed": seed, "status": "unreported", "reason": reason}


def _error_result(pair, seeds, exc, *, seconds=None, stage="scorer"):
    return {
        **pair,
        "scores": {},
        "trials": [_unreported_trial(seed, f"No seed result returned: {stage} error") for seed in seeds],
        "status": "error",
        "all_3d_seeds_succeeded": False,
        "error": f"{type(exc).__name__}: {exc}",
        "error_stage": stage,
        "seconds": seconds,
    }


def score_pair(payload):
    """Spawn-safe worker; expected seed failures retain the scorer's 2D scores."""
    pair, conformers, seeds = payload
    started = time.monotonic()
    try:
        scores, trials = similarity_scores(pair["query_smiles"], pair["control_smiles"], seeds, conformers)
        scores = dict(scores)
        by_seed = {}
        for trial in trials:
            seed = trial["seed"]
            if seed not in seeds or seed in by_seed:
                raise ValueError("Scorer returned an unknown or duplicate seed")
            by_seed[seed] = trial
        trials = [
            by_seed.get(seed, _unreported_trial(seed, "Scorer returned no record for this requested seed"))
            for seed in seeds
        ]
        complete = all(trial["status"] == "success" for trial in trials)
        if not complete:
            # The common scorer already does this. Guard the export contract as
            # well: an incomplete seed ensemble cannot become a 3D estimate.
            for metric in THREE_D_METRICS:
                scores.pop(metric, None)
                scores.pop(metric + "_stability", None)
        for metric in METRICS:
            value = scores.get(metric)
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
            ):
                raise ValueError(f"Nonfinite or nonnumeric scorer metric: {metric}")
        row = {
            **pair,
            "scores": scores,
            "trials": trials,
            "status": "success",
            "all_3d_seeds_succeeded": complete,
            "seconds": time.monotonic() - started,
        }
        json.dumps(row, allow_nan=False)
        return row
    except Exception as exc:
        return _error_result(pair, seeds, exc, seconds=time.monotonic() - started)


def _pair_results(pairs, conformers, seeds, workers):
    if workers == 1:
        for pair in pairs:
            yield score_pair((pair, conformers, seeds))
        return
    try:
        executor = concurrent.futures.ProcessPoolExecutor(
            max_workers=workers, mp_context=multiprocessing.get_context("spawn")
        )
    except Exception as exc:
        for pair in pairs:
            yield _error_result(pair, seeds, exc, stage="worker_initialization")
        return
    with executor:
        pending = {}
        for pair in pairs:
            try:
                pending[executor.submit(score_pair, (pair, conformers, seeds))] = pair
            except Exception as exc:
                yield _error_result(pair, seeds, exc, stage="worker_submission")
        for future in concurrent.futures.as_completed(pending):
            try:
                yield future.result()
            except Exception as exc:
                yield _error_result(pending[future], seeds, exc, stage="worker_execution")


def run_scoring(
    panel_path,
    output,
    conformers_path,
    seeds=(42, 43, 44),
    workers=3,
    *,
    bootstrap_seed=20261016,
    resamples=10000,
):
    """Validate sources, score every selected pair and save metric-wise summaries."""
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Output already exists; choose a new directory: {output}")
    for name, value in (("workers", workers), ("resamples", resamples)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    seeds = list(seeds)
    if not seeds or any(
        isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed <= 2**31 - 1 for seed in seeds
    ):
        raise ValueError("Scoring seeds must be nonempty integers in [0, 2**31-1]")
    if len(set(seeds)) != len(seeds):
        raise ValueError("Scoring seeds must be unique")
    if isinstance(bootstrap_seed, bool) or not isinstance(bootstrap_seed, int) or bootstrap_seed < 0:
        raise ValueError("bootstrap_seed must be a nonnegative integer")
    conformer_bytes = Path(conformers_path).read_bytes()
    conformers = json.loads(conformer_bytes)
    if not isinstance(conformers, dict):
        raise ValueError("Conformers JSON must be an object of ConformerSettings fields")
    settings = ConformerSettings(**conformers)
    panel_bytes = Path(panel_path).read_bytes()
    panel = json.loads(panel_bytes)
    validate_panel(panel, verify_sources=True)
    pairs = list(control_pairs(panel))
    if not pairs:
        raise ValueError("Selected panel contains no control pairs")
    pair_ids = {pair["pair_id"] for pair in pairs}
    if len(pair_ids) != len(pairs):
        raise ValueError("Selected panel contains duplicate pair IDs")
    root = Path(__file__).resolve().parents[2]
    code_hashes = {name: digest(root / name) for name in CODE_FILES}
    effective_conformers = asdict(settings)
    effective_conformers.pop("seed")  # Every requested trial overrides this field.
    effective_settings = {
        "conformers": effective_conformers,
        "seeds": seeds,
        "workers": workers,
        "metrics": list(METRICS),
        "bootstrap": {"seed": bootstrap_seed, "resamples": resamples},
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "panel.json").write_bytes(panel_bytes)
    (output / "conformers.json").write_bytes(conformer_bytes)
    _save(output / "control_pairs.json", pairs)
    _save(output / "public_queries.json", public_queries(panel))
    _save(
        output / "protocol.json",
        {
            "schema": "assay-matched-scoring-v1",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "panel_sha256": hashlib.sha256(panel_bytes).hexdigest(),
            "module_sha256": code_hashes["benchmarks/research/assay_matched.py"],
            "pairs_sha256": digest(output / "control_pairs.json"),
            "conformers_sha256": hashlib.sha256(conformer_bytes).hexdigest(),
            "settings_sha256": _json_hash(effective_settings),
            "effective_settings": effective_settings,
            "scoring_code_sha256": code_hashes,
            "source_payloads_verified": True,
            "source_sha256": {source["id"]: source["sha256"] for source in panel["sources"]},
            "conformers": conformers,
            "seeds": seeds,
            "workers": workers,
            "runtime": {
                "python": platform.python_version(),
                "rdkit": rdBase.rdkitVersion,
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "system": platform.system(),
                "machine": platform.machine(),
                "rdkit_base_features_sha256": digest(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"),
            },
            "scope": "Fixed same-assay similarity controls; not a generation benchmark or activity probability",
            "status_semantics": (
                "Pair status success means the scorer returned, not that all 3D seeds succeeded. "
                "all_3d_seeds_succeeded=false and missing 3D metrics retain seed failures; "
                "valid 2D scores remain available. Unexpected scorer/worker errors have status error; "
                "unreported seed records mean execution/outcome is unknown, not a measured chemical failure."
            ),
            "summary_semantics": {
                "framework-summary.json": "Each metric uses its available scores; all selected pairs enter bounds",
                "framework-paired-summary.json": (
                    "All metrics use only pairs with every 3D seed successful; "
                    "all selected pairs still enter bounds"
                ),
            },
        },
    )
    started = time.monotonic()
    results = {}

    def snapshot(stage, error=None):
        ordered = [results[pair["pair_id"]] for pair in pairs if pair["pair_id"] in results]
        _save(output / "scores.json", ordered)
        progress = {
            "expected_pairs": len(pairs),
            "scored_pairs": len(ordered),
            "failed_3d_pairs": sum(not row["all_3d_seeds_succeeded"] for row in ordered),
            "unexpected_errors": sum(row["status"] == "error" for row in ordered),
            "complete": stage == "complete",
            "stage": stage,
            "seconds": time.monotonic() - started,
        }
        if error:
            progress["error"] = error
        _save(output / "progress.json", progress)
        print(json.dumps(progress, sort_keys=True), file=sys.stderr, flush=True)
        return ordered

    snapshot("scoring")
    try:
        for row in _pair_results(pairs, conformers, seeds, workers):
            if row["pair_id"] not in pair_ids or row["pair_id"] in results:
                raise ValueError("Worker returned an unknown or duplicate pair ID")
            results[row["pair_id"]] = row
            snapshot("scoring")
        if len(results) != len(pairs):
            raise RuntimeError("Scoring stopped without a result for every selected pair")
        rows = snapshot("summarizing")
        _save(
            output / "framework-summary.json",
            {metric: aggregate_scores(panel, rows, metric, bootstrap_seed, resamples) for metric in METRICS},
        )
        paired = [{**row, "scores": row["scores"] if row["all_3d_seeds_succeeded"] else {}} for row in rows]
        _save(
            output / "framework-paired-summary.json",
            {
                metric: aggregate_scores(panel, paired, metric, bootstrap_seed, resamples)
                for metric in METRICS
            },
        )
    except BaseException as exc:
        snapshot(
            "interrupted" if isinstance(exc, KeyboardInterrupt) else "error", f"{type(exc).__name__}: {exc}"
        )
        raise
    snapshot("complete")
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--conformers", required=True, type=Path, help="JSON file of ConformerSettings fields"
    )
    parser.add_argument("--seeds", nargs="+", type=int, default=[42, 43, 44])
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--bootstrap-seed", type=int, default=20261016)
    parser.add_argument("--resamples", type=int, default=10000)
    args = parser.parse_args(argv)
    try:
        run_scoring(
            args.panel,
            args.output,
            args.conformers,
            args.seeds,
            args.workers,
            bootstrap_seed=args.bootstrap_seed,
            resamples=args.resamples,
        )
    except (OSError, ValueError, TypeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
