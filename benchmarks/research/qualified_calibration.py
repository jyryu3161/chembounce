"""Frozen development calibration: CPU rescoring and a separately invoked fit.

This tool never selects new controls, changes a scoring protocol, or runs a model.
An externally frozen procedure digest is mandatory for either data operation.
Previously inspected controls remain development data, including LOTO results.
"""

import argparse
import concurrent.futures
from dataclasses import asdict
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import importlib.metadata
import itertools
import json
import math
import multiprocessing
import os
from pathlib import Path
import platform
import statistics
import sys


SEEDS = (42, 43, 44)
METRICS = ("hit_tanimoto", "hit_shape", "hit_pharmacophore")
QUERY_FIELDS = ("target_id", "assay_id", "context_id", "query_id")
LEVELS = (Fraction(1, 100), Fraction(1, 20), Fraction(1, 10), Fraction(3, 20), Fraction(1, 5))
SETTINGS = {
    "scoring_objective": "original-hit-v1",
    "num_conformers": 8,
    "energy_window": 5.0,
    "strain_penalty": 0.05,
    "max_iterations": 500,
    "pharmacophore_sigma": 1.0,
    "reject_unspecified_stereo": True,
}
SPEC = {
    "version": "qualified-development-calibration-v1",
    "scope": "previously_inspected_development_only",
    "seeds": list(SEEDS),
    "conformers": SETTINGS,
    "fingerprint": {"kind": "binary_morgan", "radius": 2, "bits": 2048, "chirality": True},
    "features": ["hit_tanimoto", "all_seed_min_hit_shape", "all_seed_min_hit_pharmacophore"],
    "query_group": list(QUERY_FIELDS),
    "primary_rates": "equal_supported_queries_within_target_then_equal_supported_targets",
    "performance_denominator": "all_original_labelled_pairs_including_computational_failures",
    "grid_support": "complete_validated_active_common_support_training_only",
    "grid_weights": "1/(supported_targets*queries_in_target*scaffolds_in_query*rows_in_scaffold)",
    "quantile_levels": [str(q) for q in LEVELS],
    "quantile_rule": "inverse_cdf_observed_values_exact_rational_weights_grouped_ties_no_rounding",
    "grid": "cartesian_product_positive_deduplicated_quantiles_in_feature_order",
    "acceptance": "eligible_and_every_feature_greater_than_or_equal_to_cutoff",
    "constraints": {"active_macro": "4/5", "each_target_active": "3/5", "topology_active_macro": "4/5"},
    "selection": ["maximize_inactive_rejection", "maximize_active_retention", "ascending_cutoff_tuple"],
    "arithmetic": "exact_rational_weights_rates_constraints_and_objective_ties",
    "folds": "leave_one_target_out_sorted_target_ids",
    "infeasible": "abstain_all_pairs_geometric_discrimination_not_assessed_no_goal_relaxation",
    "scaffold_sensitivity": "equal_original_scaffolds_then_queries_then_targets_reporting_only",
    "bootstrap": {
        "generator": "numpy.random.Generator(numpy.random.PCG64(20261016))",
        "draws": 2000,
        "unit": "sorted_target_ids_with_replacement_same_number_as_targets",
        "quantiles": ["1/40", "39/40"],
        "interval": "inverse_empirical_cdf_no_interpolation",
        "predictions": "fixed_final_fit_and_locked_oof_separately_no_refitting",
    },
    "stability": "report_only_no_hard_range_gate_convergence_not_assessed",
}
CODE_FILES = (
    "benchmarks/research/qualified_calibration.py",
    "benchmarks/research/assay_matched.py",
    "benchmarks/research/assay_metadata.py",
    "benchmarks/research/import_chembl.py",
    "benchmarks/research/protocol.py",
    "chembounce_synthesis/conformers.py",
    "chembounce_synthesis/qualification.py",
    "chembounce_synthesis/common.py",
    "chembounce_synthesis/identity.py",
)


def json_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def json_hash(value):
    return hashlib.sha256(json_bytes(value)).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def load_procedure(path, expected_sha256):
    """Digest is supplied by the external stage gate, not inferred from the file."""
    if file_hash(path) != expected_sha256:
        raise ValueError("Frozen procedure SHA-256 mismatch")
    procedure = json.loads(Path(path).read_text())
    if procedure.get("calibration_spec") != SPEC:
        raise ValueError("Frozen calibration_spec does not match this implementation")
    if not isinstance(procedure.get("panel_sha256"), str) or len(procedure["panel_sha256"]) != 64:
        raise ValueError("Frozen procedure requires panel_sha256")
    counts = procedure.get("expected_counts", {})
    if counts != {"pairs": 197, "active": 111, "inactive": 86, "targets": 6, "queries": 12, "trials": 591}:
        raise ValueError("Frozen development panel counts differ from the reviewed 197-pair scope")
    for key in (
        "panel",
        "identity_manifest",
        "historical_scores",
        "historical_protocol",
        "independent_preregistration_review",
        "machine_spec",
        "planning_record",
    ):
        verified_reference(path, procedure, key)
    if procedure["panel"]["sha256"] != procedure["panel_sha256"]:
        raise ValueError("Frozen panel digest fields disagree")
    if json.loads(verified_reference(path, procedure, "machine_spec").read_text()) != SPEC:
        raise ValueError("Machine spec artifact differs from frozen procedure")
    return procedure


def verified_reference(procedure_path, procedure, key):
    record = procedure.get(key, {})
    if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
        raise ValueError(f"Frozen procedure requires {key} path and SHA-256")
    path = (Path(procedure_path).resolve().parent / record["path"]).resolve()
    if file_hash(path) != record["sha256"]:
        raise ValueError(f"Frozen {key} digest mismatch")
    return path


def verify_inventory(pairs, manifest):
    frozen = manifest.get("pairs", [])
    if (
        manifest.get("expected_pairs") != len(pairs)
        or manifest.get("expected_requested_trials") != len(pairs) * len(SEEDS)
        or len(frozen) != len(pairs)
        or len({p["pair_id"] for p in frozen}) != len(pairs)
    ):
        raise ValueError("Frozen identity inventory denominator mismatch")
    current = {p["pair_id"]: p for p in pairs}
    if any(
        p["pair_id"] not in current or any(current[p["pair_id"]].get(k) != v for k, v in p.items())
        for p in frozen
    ):
        raise ValueError("Frozen identity inventory differs from selected panel")


def historical_parity(rows, historical):
    """Diagnostic old/new scalar comparison; never supplies new eligibility."""
    old = {r["pair_id"]: r for r in historical}
    if len(old) != len(historical) or set(old) != {r["pair_id"] for r in rows}:
        raise ValueError("Historical parity pair inventory differs")
    comparisons = []
    for row in rows:
        previous = old[row["pair_id"]]
        if any(previous[k] != row[k] for k in ("query_smiles", "control_smiles", "label", *QUERY_FIELDS)):
            raise ValueError("Historical parity input identity differs")
        old_trials = {t["seed"]: t for t in previous["trials"]}
        if len(old_trials) != len(SEEDS) or set(old_trials) != set(SEEDS):
            raise ValueError("Historical parity seed inventory differs")
        for trial in row["trials"]:
            prior = old_trials[trial["seed"]]
            deltas = None
            if trial["status"] == prior["status"] == "success":
                info = trial["info"]
                values = {
                    "score": info["score"],
                    "hit_shape": info["hit"]["shape"],
                    "hit_pharmacophore": info["hit"]["pharmacophore"],
                }
                deltas = {key: value - prior["metrics"][key] for key, value in values.items()}
            comparisons.append(
                {
                    "pair_id": row["pair_id"],
                    "seed": trial["seed"],
                    "old_status": prior["status"],
                    "new_status": trial["status"],
                    "deltas": deltas,
                }
            )
    return {
        "scope": "Historical scalar parity only; not new qualification evidence or validation",
        "trials": comparisons,
        "compared_successes": sum(c["deltas"] is not None for c in comparisons),
        "status_changes": sum(c["old_status"] != c["new_status"] for c in comparisons),
        "max_absolute_delta": max(
            (abs(v) for c in comparisons if c["deltas"] is not None for v in c["deltas"].values()),
            default=None,
        ),
    }


def runtime():
    from rdkit.Chem import RDConfig

    feature_path = Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"
    return {
        "python": platform.python_version(),
        "executable": sys.executable,
        "packages": {name: importlib.metadata.version(name) for name in ("rdkit", "numpy", "scipy")},
        "pharmacophore_feature_definition": {"path": str(feature_path), "sha256": file_hash(feature_path)},
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        "threads": {
            key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
        },
    }


def verify_runtime(procedure):
    current = runtime()
    expected = procedure["runtime_for_scoring"]
    if (
        current["python"] != expected["python"]
        or current["packages"]["rdkit"] != expected["rdkit"]
        or current["packages"]["numpy"] != expected["numpy"]
        or Path(current["executable"]).absolute() != Path(expected["interpreter"]).absolute()
    ):
        raise ValueError("Runtime differs from frozen scoring protocol")


def source_hashes():
    root = Path(__file__).resolve().parents[2]
    return {name: file_hash(root / name) for name in CODE_FILES}


def panel_pairs(panel):
    from benchmarks.research.assay_matched import control_pairs
    from chembounce_synthesis.identity import chemical_identity

    rows, scaffolds = [], {}
    for pair in control_pairs(panel):
        query, control = (chemical_identity(pair[k]) for k in ("query_smiles", "control_smiles"))
        if (
            query.exact_smiles != pair["query_smiles"]
            or control.exact_smiles != pair["control_smiles"]
            or query.scaffold_family_smiles != pair["query_scaffold_family"]
            or control.scaffold_family_smiles != pair["control_scaffold_family"]
            or query.family_smiles != pair["query_family"]
            or control.family_smiles != pair["control_family"]
            or not query.has_scaffold
            or not control.has_scaffold
            or query.scaffold_family_smiles == control.scaffold_family_smiles
        ):
            raise ValueError("Panel exact-state/family/scaffold identity mismatch or non-hop")
        for identity in (query, control):
            scaffolds.setdefault(identity.scaffold_family_smiles, set()).add(pair["target_id"])
        rows.append(
            {
                **pair,
                "topology_change": query.ring_skeleton_smiles != control.ring_skeleton_smiles,
                "identity_policy": query.policy_version,
                "identity_rdkit_version": query.rdkit_version,
            }
        )
    if any(len(targets) > 1 for targets in scaffolds.values()):
        raise ValueError("Query/control scaffold family crosses target folds; frozen scope no longer holds")
    return rows


def counts(rows):
    return {
        "pairs": len(rows),
        "active": sum(r["label"] == "active" for r in rows),
        "inactive": sum(r["label"] == "inactive" for r in rows),
        "targets": len({r["target_id"] for r in rows}),
        "queries": len({query_key(r) for r in rows}),
        "trials": len(rows) * len(SEEDS),
    }


def query_key(row):
    return tuple(row[key] for key in QUERY_FIELDS)


def pose_record(pose, conformer_id):
    from rdkit import Chem

    conformer = pose.GetConformer(conformer_id)
    return {
        "molblock": Chem.MolToMolBlock(pose, confId=conformer_id, forceV3000=True),
        "conformer_id": conformer_id,
        "is_3d": conformer.Is3D(),
        "coordinates": conformer.GetPositions().tolist(),
    }


def restore_pose(record):
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    pose = Chem.MolFromMolBlock(record["molblock"], removeHs=False)
    coords = record["coordinates"]
    if pose is None or len(coords) != pose.GetNumAtoms() or record["is_3d"] is not True:
        raise ValueError("Malformed saved pose")
    cid = record["conformer_id"]
    if isinstance(cid, bool) or not isinstance(cid, int) or cid < 0:
        raise ValueError("Malformed saved conformer identity")
    conf = pose.GetConformer()
    conf.SetId(cid)
    conf.Set3D(True)
    for i, xyz in enumerate(coords):
        if len(xyz) != 3 or not all(finite(v) for v in xyz):
            raise ValueError("Nonfinite saved pose coordinates")
        conf.SetAtomPosition(i, Point3D(*xyz))
    return pose


def failed_row(pair, reason, stage):
    return {
        **pair,
        "status": "error",
        "all_3d_seeds_succeeded": False,
        "features": None,
        "trials": [{"seed": s, "status": "unreported", "reason": reason, "stage": stage} for s in SEEDS],
    }


def score_pair(pair):
    """One fixed block; seed failures remain explicit and never delete a pair."""
    from rdkit import DataStructs
    from chembounce_synthesis.common import FP, molecule
    from chembounce_synthesis.conformers import ConformerSettings, score_candidate
    from chembounce_synthesis.qualification import validate_result

    row = {**pair, "status": "success", "trials": []}
    row["hit_tanimoto"] = DataStructs.TanimotoSimilarity(
        FP.GetFingerprint(molecule(pair["query_smiles"])),
        FP.GetFingerprint(molecule(pair["control_smiles"])),
    )
    for seed in SEEDS:
        settings = ConformerSettings(**SETTINGS, seed=seed)
        trial = {"seed": seed, "settings": asdict(settings)}
        try:
            result = score_candidate(
                pair["control_smiles"], pair["query_smiles"], pair["query_smiles"], settings
            )
            validate_result(
                result,
                pair["control_smiles"],
                settings,
                hit=pair["query_smiles"],
                design=pair["query_smiles"],
            )
            info, hit_pose, design_pose = result
            trial.update(
                status="success",
                info=info,
                poses={
                    "hit": pose_record(hit_pose, info["candidate_conformer_id"]),
                    "design": pose_record(design_pose, info["candidate_conformer_id"]),
                },
            )
            json_bytes(trial)
        except Exception as exc:
            trial = {**trial, "status": "failed", "error_type": type(exc).__name__, "reason": str(exc)}
            trial.pop("info", None)
            trial.pop("poses", None)
        row["trials"].append(trial)
    row["all_3d_seeds_succeeded"] = all(t["status"] == "success" for t in row["trials"])
    row["status"] = "success" if row["all_3d_seeds_succeeded"] else "failed"
    row["features"] = features(row) if row["all_3d_seeds_succeeded"] else None
    if row["features"] is not None:
        row["dispersion"] = {
            metric: {
                "min": min(values),
                "max": max(values),
                "range": max(values) - min(values),
                "median": statistics.median(values),
                "sample_sd": statistics.stdev(values),
            }
            for metric in ("shape", "pharmacophore")
            for values in [[t["info"]["hit"][metric] for t in row["trials"]]]
        }
        row["convergence_status"] = "not_assessed"
    return row


def features(row):
    trials = row["trials"]
    if [t["seed"] for t in trials] != list(SEEDS) or any(t["status"] != "success" for t in trials):
        return None
    result = [row["hit_tanimoto"]] + [
        min(t["info"]["hit"][m] for t in trials) for m in ("shape", "pharmacophore")
    ]
    if not all(finite(v) and 0 <= v <= 1 for v in result):
        raise ValueError("Nonfinite or out-of-range calibration feature")
    return result


def validate_scored_row(row, pair):
    """Do not trust cached eligibility, summaries, input identities or poses."""
    from rdkit import DataStructs
    from chembounce_synthesis.common import FP, molecule
    from chembounce_synthesis.conformers import ConformerSettings
    from chembounce_synthesis.qualification import validate_result

    if any(row.get(k) != v for k, v in pair.items()):
        raise ValueError("Scored row differs from frozen pair")
    trials = row.get("trials", [])
    if [t.get("seed") for t in trials] != list(SEEDS) or any(type(t.get("seed")) is not int for t in trials):
        raise ValueError("Missing, duplicate or reordered required seeds")
    for trial in trials:
        if trial.get("status") == "success":
            settings = ConformerSettings(**SETTINGS, seed=trial["seed"])
            if trial.get("settings") != asdict(settings):
                raise ValueError("Trial settings differ from frozen protocol")
            info = trial["info"]
            poses = tuple(restore_pose(trial["poses"][role]) for role in ("hit", "design"))
            validate_result(
                (info, *poses),
                pair["control_smiles"],
                settings,
                hit=pair["query_smiles"],
                design=pair["query_smiles"],
            )
        elif trial.get("status") not in {"failed", "unreported"} or not trial.get("reason"):
            raise ValueError("Missing structured seed failure")
        elif trial["status"] == "failed" and trial.get("settings") != {
            **SETTINGS,
            "seed": trial["seed"],
        }:
            raise ValueError("Failed trial settings differ from frozen protocol")
    complete = all(t["status"] == "success" for t in trials)
    if row.get("all_3d_seeds_succeeded") is not complete:
        raise ValueError("Seed-completeness summary mismatch")
    if row.get("status") not in ({"success"} if complete else {"failed", "error"}):
        raise ValueError("Pair status disagrees with seed completion")
    if "hit_tanimoto" in row:
        expected = DataStructs.TanimotoSimilarity(
            FP.GetFingerprint(molecule(pair["query_smiles"])),
            FP.GetFingerprint(molecule(pair["control_smiles"])),
        )
        if row["hit_tanimoto"] != expected:
            raise ValueError("Binary Morgan similarity mismatch")
    if row.get("features") != (features(row) if complete else None):
        raise ValueError("Stored feature summary mismatch")


def rescore(panel_path, procedure_path, procedure_sha256, output, workers=1):
    from benchmarks.research.assay_matched import validate_panel
    from chembounce_synthesis.conformers import ConformerSettings

    if isinstance(workers, bool) or not isinstance(workers, int) or not 1 <= workers <= 4:
        raise ValueError("CPU workers must be in [1,4]")
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Rescoring requires explicit CUDA_VISIBLE_DEVICES='' before invocation")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        if os.environ.get(name) != "1":
            raise ValueError(f"Rescoring requires {name}=1 before invocation")
    procedure = load_procedure(procedure_path, procedure_sha256)
    verify_runtime(procedure)
    if asdict(ConformerSettings(**SETTINGS, seed=SEEDS[0])) != {**SETTINGS, "seed": SEEDS[0]}:
        raise ValueError("Scorer settings class differs from frozen protocol")
    if file_hash(panel_path) != procedure["panel_sha256"]:
        raise ValueError("Panel digest differs from frozen procedure")
    panel = json.loads(Path(panel_path).read_text())
    validate_panel(panel, verify_sources=True)
    pairs = panel_pairs(panel)
    if counts(pairs) != procedure["expected_counts"]:
        raise ValueError("Panel denominators differ from frozen procedure")
    subgroup = procedure["expected_identity_subgroups"]
    if subgroup != {
        "normalized_scaffold_hops": len(pairs),
        "topology_changes": sum(p["topology_change"] for p in pairs),
        "same_topology_scaffold_changes": sum(not p["topology_change"] for p in pairs),
    }:
        raise ValueError("Fixed scaffold/topology subgroup counts differ")
    identity_path = verified_reference(procedure_path, procedure, "identity_manifest")
    historical_path = verified_reference(procedure_path, procedure, "historical_scores")
    verify_inventory(pairs, json.loads(identity_path.read_text()))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    (output / "panel.json").write_bytes(Path(panel_path).read_bytes())
    (output / "procedure.json").write_bytes(Path(procedure_path).read_bytes())
    (output / "identity_manifest.json").write_bytes(identity_path.read_bytes())
    save(output / "pairs.json", pairs)
    protocol = {
        "schema": SPEC["version"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "procedure_sha256": procedure_sha256,
        "panel_sha256": file_hash(panel_path),
        "pairs_sha256": file_hash(output / "pairs.json"),
        "calibration_spec": SPEC,
        "runtime": runtime(),
        "source_sha256": source_hashes(),
        "workers": workers,
        "source_payloads_verified": True,
        "expected_counts": counts(pairs),
        "identity_manifest_sha256": file_hash(identity_path),
        "historical_scores_sha256": file_hash(historical_path),
    }
    save(output / "protocol.json", protocol)
    by_id = {p["pair_id"]: failed_row(p, "Not yet reported", "pending") for p in pairs}
    completed = set()

    def snapshot(status):
        save(output / "scores.json", [by_id[p["pair_id"]] for p in pairs])
        save(
            output / "progress.json",
            {
                "status": status,
                "complete": status == "complete",
                "expected_pairs": len(pairs),
                "reported_pairs": len(completed),
                "expected_trials": len(pairs) * len(SEEDS),
            },
        )

    def retain(pair, row):
        if row.get("pair_id") != pair["pair_id"]:
            row = failed_row(pair, "Worker returned incorrect pair identity", "worker_execution")
        try:
            validate_scored_row(row, pair)
        except Exception as exc:
            row = failed_row(pair, f"Invalid worker payload: {type(exc).__name__}: {exc}", "worker_payload")
        by_id[pair["pair_id"]] = row
        completed.add(pair["pair_id"])
        snapshot("scoring")

    snapshot("scoring")
    try:
        if workers == 1:
            for pair in pairs:
                try:
                    row = score_pair(pair)
                except Exception as exc:
                    row = failed_row(pair, f"{type(exc).__name__}: {exc}", "worker_execution")
                retain(pair, row)
        else:
            try:
                executor = concurrent.futures.ProcessPoolExecutor(
                    max_workers=workers, mp_context=multiprocessing.get_context("spawn")
                )
            except Exception as exc:
                for pair in pairs:
                    retain(pair, failed_row(pair, f"{type(exc).__name__}: {exc}", "worker_initialization"))
                executor = None
            if executor is not None:
                with executor:
                    pending = {}
                    for pair in pairs:
                        try:
                            pending[executor.submit(score_pair, pair)] = pair
                        except Exception as exc:
                            retain(pair, failed_row(pair, str(exc), "worker_submission"))
                    for future in concurrent.futures.as_completed(pending):
                        pair = pending[future]
                        try:
                            row = future.result()
                        except Exception as exc:
                            row = failed_row(pair, f"{type(exc).__name__}: {exc}", "worker_execution")
                        retain(pair, row)
        snapshot("complete")
        save(
            output / "parity.json",
            historical_parity([by_id[p["pair_id"]] for p in pairs], json.loads(historical_path.read_text())),
        )
        save(
            output / "artifacts.json",
            {
                name: file_hash(output / name)
                for name in (
                    "panel.json",
                    "procedure.json",
                    "identity_manifest.json",
                    "pairs.json",
                    "protocol.json",
                    "scores.json",
                    "progress.json",
                    "parity.json",
                )
            },
        )
    except BaseException:
        snapshot("interrupted")
        raise
    return [by_id[p["pair_id"]] for p in pairs]


def inverse_cdf(values_and_weights, quantile):
    grouped = {}
    for value, weight in values_and_weights:
        if not finite(value) or weight <= 0:
            raise ValueError("Inverse CDF requires finite values and positive weights")
        grouped[value] = grouped.get(value, Fraction()) + weight
    total, cumulative = sum(grouped.values(), Fraction()), Fraction()
    if not grouped or not 0 < quantile <= 1:
        raise ValueError("Inverse CDF requires support and quantile in (0,1]")
    for value in sorted(grouped):
        cumulative += grouped[value]
        if cumulative >= quantile * total:
            return value
    raise AssertionError("Unreachable inverse CDF endpoint")


def weighted_grid(rows):
    supported = [r for r in rows if r["label"] == "active" and r.get("features") is not None]
    hierarchy = {}
    for row in supported:
        hierarchy.setdefault(row["target_id"], {}).setdefault(query_key(row), {}).setdefault(
            row["control_scaffold_family"], []
        ).append(row)
    weights = {}
    for queries in hierarchy.values():
        for scaffolds in queries.values():
            for local in scaffolds.values():
                for row in local:
                    weights[row["pair_id"]] = Fraction(
                        1, len(hierarchy) * len(queries) * len(scaffolds) * len(local)
                    )
    grid, level_values = [], {}
    for index, metric in enumerate(METRICS):
        points = [(r["features"][index], weights[r["pair_id"]]) for r in supported]
        levels = {str(q): inverse_cdf(points, q) for q in LEVELS} if points else {}
        level_values[metric] = levels
        grid.append(sorted({v for v in levels.values() if v > 0}))
    total = sum(weights.values(), Fraction())
    squared = sum((w * w for w in weights.values()), Fraction())
    diagnostics = {
        "supported_active_rows": len(supported),
        "excluded_active_rows": sum(r["label"] == "active" for r in rows) - len(supported),
        "supported_targets": sorted(hierarchy),
        "row_weights": {key: str(value) for key, value in sorted(weights.items())},
        "effective_sample_size": str(total * total / squared) if squared else None,
        "maximum_row_weight": str(max(weights.values())) if weights else None,
        "quantile_values": level_values,
        "grid": dict(zip(METRICS, grid)),
        "support_by_query": [
            {
                **dict(zip(QUERY_FIELDS, key)),
                "active_original": sum(r["label"] == "active" for r in local),
                "active_complete": sum(
                    r["label"] == "active" and r.get("features") is not None for r in local
                ),
                "scaffolds": {
                    s: {
                        "original": sum(
                            r["label"] == "active" and r["control_scaffold_family"] == s for r in local
                        ),
                        "complete": sum(
                            r["label"] == "active"
                            and r["control_scaffold_family"] == s
                            and r.get("features") is not None
                            for r in local
                        ),
                    }
                    for s in sorted({r["control_scaffold_family"] for r in local if r["label"] == "active"})
                },
            }
            for key, local in groups(rows, query_key).items()
        ],
    }
    return list(itertools.product(*grid)), diagnostics


def groups(rows, key):
    result = {}
    for row in rows:
        result.setdefault(key(row), []).append(row)
    return dict(sorted(result.items()))


def accepted(row, cutoff):
    return (
        cutoff is not None
        and row.get("features") is not None
        and all(v >= c for v, c in zip(row["features"], cutoff))
    )


def rates(rows, predictions, label, *, topology=False, scaffold=False):
    """Exact full-denominator rates; unsupported original groups remain null."""
    by_target, queries, unsupported_queries = {}, [], 0
    for key, original in groups(rows, query_key).items():
        selected = [r for r in original if r["label"] == label and (not topology or r["topology_change"])]
        if not selected:
            unsupported_queries += 1
            continue
        local = (
            list(groups(selected, lambda r: r["control_scaffold_family"]).values())
            if scaffold
            else [selected]
        )
        score = sum(
            (
                Fraction(
                    sum(
                        predictions[r["pair_id"]] if label == "active" else not predictions[r["pair_id"]]
                        for r in group
                    ),
                    len(group),
                )
                for group in local
            ),
            Fraction(),
        ) / len(local)
        by_target.setdefault(key[0], []).append(score)
        queries.append({"query": list(key), "denominator": len(selected), "rate": str(score)})
    targets = {t: sum(values, Fraction()) / len(values) for t, values in by_target.items()}
    macro = sum(targets.values(), Fraction()) / len(targets) if targets else None
    return {
        "macro": macro,
        "targets": targets,
        "queries": queries,
        "unsupported_queries": unsupported_queries,
        "unsupported_targets": len({r["target_id"] for r in rows} - set(targets)),
    }


def export_rates(report):
    return {
        **report,
        "macro": float(report["macro"]) if report["macro"] is not None else None,
        "macro_rational": str(report["macro"]) if report["macro"] is not None else None,
        "targets": {t: {"rate": float(v), "rational": str(v)} for t, v in report["targets"].items()},
    }


def prediction_report(rows, predictions, *, fitted_targets):
    fitted_targets = set(fitted_targets)
    fit_rows = [r for r in rows if r["target_id"] in fitted_targets]
    complete = [r for r in fit_rows if r.get("features") is not None]
    result = {
        "active_retention": export_rates(rates(rows, predictions, "active")),
        "inactive_operational_rejection": export_rates(rates(rows, predictions, "inactive")),
        "topology_active_retention": export_rates(rates(rows, predictions, "active", topology=True)),
        "scaffold_sensitivity": {
            label: export_rates(rates(rows, predictions, label, scaffold=True))
            for label in ("active", "inactive")
        },
        "fit_coverage": {
            "targets": sorted(fitted_targets),
            "pairs": len(fit_rows),
            "no_fit_pairs": len(rows) - len(fit_rows),
        },
        "computational_failure_pairs": sum(r.get("features") is None for r in rows),
        "geometric_discrimination": {
            "status": "assessed_conditional_on_fit_and_complete_computation" if complete else "not_assessed",
            "pairs": len(complete),
            "rates": {
                label: export_rates(rates(complete, predictions, label)) for label in ("active", "inactive")
            },
        },
        "rejection_counts": {
            "no_fit": sum(r["target_id"] not in fitted_targets for r in rows),
            "computational_failure_in_fitted_targets": sum(r.get("features") is None for r in fit_rows),
            "geometric_threshold_rejection": sum(not predictions[r["pair_id"]] for r in complete),
        },
        "rejection_counts_by_label": {
            label: {
                "no_fit": sum(r["label"] == label and r["target_id"] not in fitted_targets for r in rows),
                "computational_failure_in_fitted_targets": sum(
                    r["label"] == label and r.get("features") is None for r in fit_rows
                ),
                "geometric_threshold_rejection": sum(
                    r["label"] == label and not predictions[r["pair_id"]] for r in complete
                ),
            }
            for label in ("active", "inactive")
        },
    }
    return result


def fit_thresholds(rows):
    candidates, support = weighted_grid(rows)
    target_ids = {r["target_id"] for r in rows}
    evaluated, best, best_key = [], None, None
    for cutoff in candidates:
        predictions = {r["pair_id"]: accepted(r, cutoff) for r in rows}
        active, inactive, topology = (
            rates(rows, predictions, "active"),
            rates(rows, predictions, "inactive"),
            rates(rows, predictions, "active", topology=True),
        )
        reasons = []
        if active["macro"] is None or active["macro"] < Fraction(4, 5):
            reasons.append("active_macro_below_4/5_or_unsupported")
        if set(active["targets"]) != target_ids or any(
            v < Fraction(3, 5) for v in active["targets"].values()
        ):
            reasons.append("training_target_active_below_3/5_or_unsupported")
        if topology["macro"] is None or topology["macro"] < Fraction(4, 5):
            reasons.append("topology_active_macro_below_4/5_or_unsupported")
        if inactive["macro"] is None:
            reasons.append("inactive_objective_unsupported")
        evaluated.append(
            {
                "cutoff": list(cutoff),
                "cutoff_hex": [float(v).hex() for v in cutoff],
                "feasible": not reasons,
                "reasons": reasons,
                "active": export_rates(active),
                "inactive_rejection": export_rates(inactive),
                "topology_active": export_rates(topology),
            }
        )
        if not reasons:
            key = (-inactive["macro"], -active["macro"], *cutoff)
            if best_key is None or key < best_key:
                best, best_key = cutoff, key
    return {
        "fit_status": "fitted" if best is not None else "infeasible",
        "cutoff": list(best) if best is not None else None,
        "support": support,
        "candidates": evaluated,
        "infeasible_reason": None
        if best is not None
        else "empty_positive_grid_or_no_jointly_feasible_triple",
    }


def bootstrap(rows, predictions):
    import numpy as np

    targets = sorted({r["target_id"] for r in rows})
    draws = np.random.Generator(np.random.PCG64(20261016)).integers(
        0, len(targets), size=(2000, len(targets))
    )
    reports = {
        "active_retention": rates(rows, predictions, "active"),
        "inactive_operational_rejection": rates(rows, predictions, "inactive"),
        "topology_active_retention": rates(rows, predictions, "active", topology=True),
    }
    metrics = {}
    for name, report in reports.items():
        samples = []
        for draw in draws:
            values = [report["targets"][targets[i]] for i in draw if targets[i] in report["targets"]]
            if values:
                samples.append(sum(values, Fraction()) / len(values))
        ordered = sorted(samples)
        interval = (
            [float(ordered[math.ceil(q * len(ordered)) - 1]) for q in (Fraction(1, 40), Fraction(39, 40))]
            if ordered
            else None
        )
        metrics[name] = {
            "interval": interval,
            "valid_replicates": len(ordered),
            "unassessed_replicates": len(draws) - len(ordered),
        }
    return {
        "targets": targets,
        "draws": draws.tolist(),
        "draws_sha256": json_hash(draws.tolist()),
        "numpy_version": np.__version__,
        "metrics": metrics,
        "scope": "Conditional descriptive target-cluster intervals of fixed predictions; no refitting or calibrated uncertainty",
    }


def fit_development(rows):
    """Pure development fitting entrypoint; CLI also verifies frozen artifacts."""
    target_ids = sorted({r["target_id"] for r in rows})
    if len({r["pair_id"] for r in rows}) != len(rows) or not target_ids:
        raise ValueError("Fitting requires unique nonempty pair identities")
    for row in rows:
        feature = row.get("features")
        if row["label"] not in {"active", "inactive"} or not isinstance(row["topology_change"], bool):
            raise ValueError("Malformed fixed label/subgroup")
        if feature is not None and (len(feature) != 3 or not all(finite(v) and 0 <= v <= 1 for v in feature)):
            raise ValueError("Malformed calibration features")
    folds, oof, fitted_targets = [], {}, []
    for target in target_ids:
        train = [r for r in rows if r["target_id"] != target]
        test = [r for r in rows if r["target_id"] == target]
        fitted = fit_thresholds(train)
        predictions = {r["pair_id"]: accepted(r, fitted["cutoff"]) for r in test}
        if fitted["cutoff"] is not None:
            fitted_targets.append(target)
        oof.update(predictions)
        folds.append(
            {
                "test_target": target,
                "training_targets": sorted(set(target_ids) - {target}),
                **fitted,
                "test_predictions": predictions,
                "test_report": prediction_report(
                    test, predictions, fitted_targets=[target] if fitted["cutoff"] is not None else []
                ),
            }
        )
    final = fit_thresholds(rows)
    final_predictions = {r["pair_id"]: accepted(r, final["cutoff"]) for r in rows}
    return {
        "scope": SPEC["scope"],
        "folds": folds,
        "oof": {
            "predictions": oof,
            "report": prediction_report(rows, oof, fitted_targets=fitted_targets),
            "bootstrap": bootstrap(rows, oof),
        },
        "final": {
            **final,
            "predictions": final_predictions,
            "report": prediction_report(
                rows, final_predictions, fitted_targets=target_ids if final["cutoff"] is not None else []
            ),
            "bootstrap": bootstrap(rows, final_predictions),
        },
        "minimum_scores": dict(zip(METRICS, final["cutoff"])) if final["cutoff"] is not None else None,
    }


def fit(scores_dir, procedure_path, procedure_sha256, output):
    procedure = load_procedure(procedure_path, procedure_sha256)
    verify_runtime(procedure)
    source = Path(scores_dir)
    artifacts = json.loads((source / "artifacts.json").read_text())
    required = {
        "panel.json",
        "procedure.json",
        "identity_manifest.json",
        "pairs.json",
        "protocol.json",
        "scores.json",
        "progress.json",
        "parity.json",
    }
    if set(artifacts) != required or any(file_hash(source / name) != artifacts[name] for name in required):
        raise ValueError("Scoring artifact integrity mismatch")
    protocol = json.loads((source / "protocol.json").read_text())
    if (
        artifacts["procedure.json"] != procedure_sha256
        or artifacts["panel.json"] != procedure["panel_sha256"]
        or protocol["calibration_spec"] != SPEC
        or protocol["procedure_sha256"] != procedure_sha256
        or protocol["pairs_sha256"] != artifacts["pairs.json"]
    ):
        raise ValueError("Scoring protocol does not match frozen procedure")
    if protocol["source_sha256"] != source_hashes():
        raise ValueError("Scoring/fit source changed; review and freeze a new run")
    current_runtime = runtime()
    if (
        protocol["runtime"]["packages"] != current_runtime["packages"]
        or protocol["runtime"]["pharmacophore_feature_definition"]["sha256"]
        != current_runtime["pharmacophore_feature_definition"]["sha256"]
    ):
        raise ValueError("Scoring/fit scientific runtime changed")
    if not json.loads((source / "progress.json").read_text()).get("complete"):
        raise ValueError("Incomplete scoring run cannot be fitted")
    pairs = panel_pairs(json.loads((source / "panel.json").read_text()))
    identity_path = verified_reference(procedure_path, procedure, "identity_manifest")
    if artifacts["identity_manifest.json"] != file_hash(identity_path):
        raise ValueError("Frozen identity manifest differs from scoring run")
    verify_inventory(pairs, json.loads((source / "identity_manifest.json").read_text()))
    if (
        json.loads((source / "pairs.json").read_text()) != pairs
        or counts(pairs) != procedure["expected_counts"]
    ):
        raise ValueError("Frozen pair identity or denominator mismatch")
    rows = json.loads((source / "scores.json").read_text())
    if len(rows) != len(pairs):
        raise ValueError("Scoring rows omitted original pairs")
    for row, pair in zip(rows, pairs):
        validate_scored_row(row, pair)
    result = fit_development(rows)
    result["provenance"] = {
        "procedure_sha256": procedure_sha256,
        "scoring_artifacts": artifacts,
        "runtime": runtime(),
        "source_sha256": source_hashes(),
        "counts": counts(rows),
    }
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    save(output / "fit.json", result)
    artifacts = {"fit.json": file_hash(output / "fit.json")}
    if result["minimum_scores"] is not None:
        save(
            output / "preset-candidate.json",
            {
                "schema": "qualified-development-preset-candidate-v1",
                "status": "development_fit_requires_independent_review",
                "scope": SPEC["scope"],
                "minimum_scores": result["minimum_scores"],
                "calibration_spec": SPEC,
                "fit_sha256": artifacts["fit.json"],
                "provenance": result["provenance"],
                "final_development_report": result["final"]["report"],
                "oof_development_report": result["oof"]["report"],
                "fold_statuses": [
                    {
                        "test_target": fold["test_target"],
                        "fit_status": fold["fit_status"],
                        "cutoff": fold["cutoff"],
                    }
                    for fold in result["folds"]
                ],
                "limitations": [
                    "Previously inspected development panel; LOTO is not unseen validation",
                    "Observed three-seed minima and spread are not calibrated uncertainty or convergence",
                    "Scalar pose selection may abstain despite another pose meeting component minima",
                    "Real-design overlay availability is an additional production veto not measured by these controls",
                    "No biological activity retention guarantee, synthesis probability, or established 3D advantage",
                    "Independent review and root stage gate remain required before publishing a preset",
                ],
            },
        )
        artifacts["preset-candidate.json"] = file_hash(output / "preset-candidate.json")
    save(output / "artifacts.json", artifacts)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("spec", help="Print mechanical protocol only; no data access or fitting")
    for name in ("rescore", "fit"):
        cmd = sub.add_parser(name)
        cmd.add_argument("--procedure", type=Path, required=True)
        cmd.add_argument("--procedure-sha256", required=True)
        cmd.add_argument("--output", type=Path, required=True)
        if name == "rescore":
            cmd.add_argument("--panel", type=Path, required=True)
            cmd.add_argument("--workers", type=int, default=1)
        else:
            cmd.add_argument("--scores-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "spec":
            print(json.dumps(SPEC, indent=2, sort_keys=True))
        elif args.command == "rescore":
            rescore(args.panel, args.procedure, args.procedure_sha256, args.output, args.workers)
        else:
            fit(args.scores_dir, args.procedure, args.procedure_sha256, args.output)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
