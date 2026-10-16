"""End-to-end synthesis-aware scaffold hopping and hit expansion."""

import csv
import copy
import importlib.metadata
import json
import logging
import math
import subprocess
import time
from collections import deque
from contextlib import ExitStack
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold

from .catalog import Catalog
from .common import FP, canonical, digest, molecule, positive, write_json
from .conformers import (
    ConformerSettings,
    _check_stereo,
    ensemble,
    prepare_reference,
    reference_from_sdf,
    score_candidate,
)
from .continuation import load_carried_routes
from .design_selection import select_designs, selection_policy, sized_count
from .evidence import evaluator_from_config, validate_evidence_config
from .identity import classify_change, family_key, scaffold_family_key
from .pocket import constraints_from_config, validate_pocket_config
from .protocols import apply_synthesis_protocol, protocol_provenance
from .qualification import decide, effective_policy, validate_result
from .reactions import Rules
from .refinement import RefinementBudget, RefinementSettings, refine_candidate
from .route_objective import RootProductObjective, search_policy
from .routes import SearchLimits, expand_route, resolve_proposal, verify_route
from .synllama import SavedProposals, SynLlama

LOGGER = logging.getLogger(__name__)
COLUMNS = [
    "rank",
    "recommendation_status",
    "qualification_policy",
    "smiles",
    "design",
    "score",
    "pareto_front",
    "hit_shape",
    "hit_pharmacophore",
    "design_shape",
    "design_pharmacophore",
    "hit_tanimoto",
    "design_tanimoto",
    "hit_morgan",
    "relative_energy",
    "reaction_steps",
    "route_id",
    "alternative_route_ids",
    "alternative_family_candidates",
    "scoring_stage",
    "pose_seed",
    "seed_score_median",
    "seed_score_min",
    "seed_score_max",
    "seed_successes",
    "seed_count",
    "family_smiles",
    "scaffold_family_smiles",
    "ranking_score",
    "ranking_metrics",
    "seed_score_range",
    "seed_score_std",
    "convergence_status",
    "refinement_conformers",
]


def scaffold(smiles, *, isomeric=True):
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(molecule(smiles)), isomericSmiles=isomeric)


def check_constraints(product, hit, design, mode="hop", core_smiles=None, preserve_design_scaffold=False):
    mol = molecule(product)
    if canonical(product) == canonical(hit) and mode == "hop":
        raise ValueError("Reconstruction returned the original hit")
    if mode == "hop":
        change = classify_change(hit, product)
        if change == "state_only":
            raise ValueError(
                "Reconstructed product changed only chemical state, not the original Murcko scaffold"
            )
        if change != "scaffold":
            raise ValueError("Reconstructed product did not change the original Murcko scaffold")
    elif mode == "expand":
        if classify_change(design, product) in {"unchanged", "state_only"}:
            raise ValueError("Hit expansion returned the unchanged design chemical family")
        if scaffold(product) != scaffold(design):
            raise ValueError("Hit expansion changed the selected design's Murcko scaffold")
    else:
        raise ValueError("mode must be hop or expand")
    if preserve_design_scaffold and scaffold(product) != scaffold(design):
        raise ValueError("Reconstruction did not preserve the selected design scaffold")
    if core_smiles and not mol.HasSubstructMatch(molecule(core_smiles), useChirality=True):
        raise ValueError("Reconstructed product lost the required core/anchor substructure")


def _config(config):
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise ValueError("Synthesis configuration requires schema_version=1")
    config = apply_synthesis_protocol(config)
    permitted = {
        "schema_version",
        "protocol",
        "catalog",
        "rules",
        "provider",
        "search",
        "conformers",
        "mode",
        "core_smiles",
        "preserve_design_scaffold",
        "max_designs",
        "max_design_attempts",
        "max_scored_products",
        "max_results",
        "expansion",
        "reference_sdf",
        "minimum_scores",
        "property_thresholds",
        "refinement",
        "carry_forward",
        "allocation_policy",
        "synthesis_evidence",
        "pocket",
        "recommendation",
        "route_search_policy",
        "route_search_cache_size",
    }
    if set(config) - permitted:
        raise ValueError(f"Unknown synthesis configuration keys: {sorted(set(config) - permitted)}")
    for required in ("catalog", "rules", "provider"):
        if required not in config:
            raise ValueError(f"Missing synthesis configuration: {required}")
    for key in (
        "provider",
        "search",
        "conformers",
        "expansion",
        "minimum_scores",
        "property_thresholds",
        "refinement",
    ):
        if key in config and not isinstance(config[key], dict):
            raise ValueError(f"{key} must be a JSON object")
    if "preserve_design_scaffold" in config and not isinstance(config["preserve_design_scaffold"], bool):
        raise ValueError("preserve_design_scaffold must be true or false")
    if not isinstance(config.get("allocation_policy", "balanced"), str) or config.get(
        "allocation_policy", "balanced"
    ) not in {"balanced", "encounter"}:
        raise ValueError("allocation_policy must be balanced or encounter")
    SearchLimits(**config.get("search", {}))
    ConformerSettings(**config.get("conformers", {}))
    RefinementSettings(**config.get("refinement", {}))
    validate_evidence_config(config.get("synthesis_evidence", {}))
    if "pocket" in config:
        validate_pocket_config(config["pocket"])
        if not config.get("reference_sdf"):
            raise ValueError("pocket requires an explicit reference_sdf")
    carried = config.get("carry_forward", [])
    if not isinstance(carried, list) or not all(isinstance(p, str) and p for p in carried):
        raise ValueError("carry_forward must be a list of run directories")
    if config.get("core_smiles") is not None:
        canonical(config["core_smiles"])
    for key, default in (("max_designs", 20), ("max_scored_products", 256), ("max_results", 50)):
        positive(config.get(key, default), key)
    if config.get("mode", "hop") not in {"hop", "expand"}:
        raise ValueError("mode must be hop or expand")
    minimum = config.get("minimum_scores", {})
    for key, value in minimum.items():
        if (
            key
            not in {
                "hit_shape",
                "hit_pharmacophore",
                "design_shape",
                "design_pharmacophore",
                "hit_tanimoto",
                "design_tanimoto",
            }
            or not isinstance(value, (int, float))
            or isinstance(value, bool)
            or not math.isfinite(value)
            or not 0 <= value <= 1
        ):
            raise ValueError("minimum_scores must contain known score names with values in [0, 1]")
    exp = config.get("expansion", {})
    if set(exp) - {
        "rounds",
        "max_routes",
        "locked_building_block_ids",
        "building_block_top_k",
        "min_tanimoto",
    }:
        raise ValueError("Unknown expansion configuration key")
    if (
        isinstance(exp.get("rounds", 0), bool)
        or not isinstance(exp.get("rounds", 0), int)
        or exp.get("rounds", 0) < 0
    ):
        raise ValueError("expansion rounds must be a nonnegative integer")
    positive(exp.get("max_routes", 64), "expansion.max_routes")
    if "building_block_top_k" in exp:
        positive(exp["building_block_top_k"], "expansion.building_block_top_k")
    minimum_bb = exp.get("min_tanimoto", 0.0)
    if (
        isinstance(minimum_bb, bool)
        or not isinstance(minimum_bb, (int, float))
        or not math.isfinite(minimum_bb)
        or not 0 <= minimum_bb <= 1
    ):
        raise ValueError("expansion.min_tanimoto must be in [0, 1]")
    locks = exp.get("locked_building_block_ids", [])
    if not isinstance(locks, list) or not all(isinstance(value, str) and value for value in locks):
        raise ValueError("locked_building_block_ids must be a list of nonempty catalog IDs")
    # Reuse the existing public property vocabulary, validating before any expensive work.
    from chembounce import _validate_parameters

    _validate_parameters(0.5, 0.3, 42, config.get("property_thresholds", {}))
    effective_policy(config)
    selection_policy(config)
    search_policy(config)
    return config


def load_config(filename):
    path = Path(filename).resolve()
    return _resolve_paths(_config(json.loads(path.read_text())), path.parent)


def _resolve_paths(config, base):
    config = copy.deepcopy(config)
    for key in ("catalog", "rules", "reference_sdf"):
        if config.get(key):
            config[key] = str((base / config[key]).resolve())
    if "carry_forward" in config:
        config["carry_forward"] = [str((base / p).resolve()) for p in config["carry_forward"]]
    if config.get("pocket", {}).get("receptor_pdb"):
        config["pocket"]["receptor_pdb"] = str((base / config["pocket"]["receptor_pdb"]).resolve())
    evidence_predictor = config.get("synthesis_evidence", {}).get("predictor")
    if evidence_predictor is not None:
        evidence_predictor["model_path"] = str((base / evidence_predictor["model_path"]).resolve())
    provider = config["provider"] = dict(config["provider"])
    for key in ("model_path", "prompt_path", "path"):
        if provider.get(key):
            provider[key] = str((base / provider[key]).resolve())
    return config


def _provider(config):
    config = dict(config)
    kind = config.pop("kind", None)
    if kind == "saved":
        if set(config) != {"path"}:
            raise ValueError("saved provider requires only path")
        return SavedProposals(config["path"])
    if kind == "synllama":
        if set(config) - {"model_path", "prompt_path", "device", "samples", "max_new_tokens", "seed"}:
            raise ValueError("Unknown SynLlama provider option")
        return SynLlama(**config)
    raise ValueError("Provider kind must be synllama or saved")


def _evidence_evaluator(config):
    return evaluator_from_config(config.get("synthesis_evidence", {}))


def _source_snapshot():
    return {p.name: digest(p) for p in sorted(Path(__file__).parent.glob("*.py"))}


# Snapshot at import, as with the main workflow, to detect edits in a running process.
_SOURCE_HASHES = _source_snapshot()


def _manifest(config, hit, designs):
    try:
        root = Path(__file__).parent
        revision = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, stderr=subprocess.DEVNULL, text=True
        ).strip()
        dirty = bool(
            subprocess.check_output(
                ["git", "status", "--porcelain"], cwd=root, stderr=subprocess.DEVNULL, text=True
            ).strip()
        )
    except (OSError, subprocess.CalledProcessError):
        revision, dirty = None, None
    return {
        "schema_version": 1,
        "hit": hit,
        "designs": designs,
        "config": config,
        "protocol": protocol_provenance(config.get("protocol", "legacy")),
        "code_sha256": _SOURCE_HASHES,
        "sources_changed_since_import": _source_snapshot() != _SOURCE_HASHES,
        "git_revision": revision,
        "git_dirty": dirty,
        "versions": {name: importlib.metadata.version(name) for name in ("rdkit", "numpy", "scipy")},
        "chemical_state_policy": "Explicit input states; no implicit protonation, tautomer, or stereoisomer enumeration",
        "verification_scope": "Catalog snapshot membership and registered reaction replay; no experimental synthesis or current-stock guarantee",
        "ranking_scope": "Heuristic, bounded search; joint shape/pharmacophore scores are not calibrated activity probabilities",
        "refinement_policy": "Every seed must succeed at the final ensemble size; median-score seed supplies the exported pose. Threshold acceptance defaults to every seed. Optional seed-range convergence and adaptive stages are budgeted globally in requested candidate conformers, excluding reference ensembles. Ranking diagnostics are separate from pose components; failed primary designs permit bounded alternative-design fallback.",
        "fingerprint": {
            "algorithm": "Morgan",
            "radius": 2,
            "bits": 2048,
            "include_chirality": True,
            "similarity": "Tanimoto",
            "hit_morgan_alias": "hit_tanimoto",
        },
        "ranking_metrics": [
            "hit_tanimoto",
            "design_tanimoto",
            "hit_shape",
            "hit_pharmacophore",
            "design_shape",
            "design_pharmacophore",
            "relative_energy",
            "reaction_steps",
        ],
    }


def _pareto(rows):
    if not rows:
        return
    modes = {r.get("scoring_stage") == "2d_screening" for r in rows}
    if len(modes) > 1:
        raise ValueError("Cannot mix 2D and 3D candidates in one ranking")
    if modes == {True}:
        points = np.asarray([(r["hit_tanimoto"], r["design_tanimoto"], -r["reaction_steps"]) for r in rows])
    else:
        points = np.asarray(
            [
                (
                    r.get("ranking_metrics", {}).get("hit_shape", r["hit_shape"]),
                    r.get("ranking_metrics", {}).get("hit_pharmacophore", r["hit_pharmacophore"]),
                    r.get("ranking_metrics", {}).get("design_shape", r["design_shape"]),
                    r.get("ranking_metrics", {}).get("design_pharmacophore", r["design_pharmacophore"]),
                    r["hit_tanimoto"],
                    r["design_tanimoto"],
                    -r.get("ranking_metrics", {}).get("relative_energy", r["relative_energy"]),
                    -r["reaction_steps"],
                )
                for r in rows
            ]
        )
    if not np.isfinite(points).all():
        raise ValueError("Pareto ranking requires finite metrics")
    # Compute dominance once. Rechecking every remaining pair on every front
    # makes a fully ordered pool cubic in the number of candidates.
    dominates = np.asarray(
        [np.all(point >= points, axis=1) & np.any(point > points, axis=1) for point in points]
    )
    counts = dominates.sum(axis=0)
    front = np.flatnonzero(counts == 0)
    level = 0
    while len(front):
        for i in front:
            rows[i]["pareto_front"] = level
        counts[front] = -1
        counts -= dominates[front].sum(axis=0)
        front = np.flatnonzero(counts == 0)
        level += 1


def _rank(rows):
    _pareto(rows)
    rows.sort(
        key=lambda r: (
            r["pareto_front"],
            -r.get("ranking_score", r["seed_score_median"]),
            r["reaction_steps"],
            r["smiles"],
            r["design"],
        )
    )
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
    return rows


def _write_table(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                row
                | {
                    "alternative_route_ids": json.dumps(row["alternative_route_ids"]),
                    "alternative_family_candidates": json.dumps(row.get("alternative_family_candidates", [])),
                    "ranking_metrics": json.dumps(row.get("ranking_metrics", {}), sort_keys=True),
                }
            )


def _balanced_order(items, *, product, design):
    """Round-robin scaffold families and designs, preserving within-group order."""
    groups = {}
    for item in items:
        scaffold_id = scaffold_family_key(product(item))
        groups.setdefault(scaffold_id, {}).setdefault(design(item), deque()).append(item)
    queue = deque(deque(group.values()) for group in groups.values())
    ordered = []
    while queue:
        designs = queue.popleft()
        bucket = designs.popleft()
        ordered.append(bucket.popleft())
        if bucket:
            designs.append(bucket)
        if designs:
            queue.append(designs)
    return ordered


def _candidate_order(candidates, policy):
    # Existing carried evidence has priority, but all eligible designs share
    # each priority tier's compute budget instead of first-design exhaustion.
    ordered = []
    for carried in (True, False):
        tier = [item for item in candidates if any(r.get("carried_from") for r in item["records"]) == carried]
        if policy == "balanced":
            tier = _balanced_order(
                tier, product=lambda item: item["product"], design=lambda item: item["design"]
            )
        ordered.extend(tier)
    return ordered


def _family_groups(rows):
    groups = {}
    for row in rows:
        groups.setdefault(row["family_smiles"], []).append(row)
    for alternatives in groups.values():
        alternatives.sort(
            key=lambda row: (
                row["pareto_front"],
                -row["score"],
                row["reaction_steps"],
                row["design"],
                row["smiles"],
            )
        )
    return groups


def _with_alternative_routes(row, alternatives):
    route_ids = dict.fromkeys(
        route_id
        for alternative in alternatives
        if alternative["smiles"] == row["smiles"]
        for route_id in [alternative["route_id"], *alternative["alternative_route_ids"]]
        if route_id != row["route_id"]
    )
    family_alternatives = [
        {"smiles": alternative["smiles"], "design": alternative["design"], "route_id": route_id}
        for alternative in alternatives
        if alternative["smiles"] != row["smiles"]
        for route_id in [alternative["route_id"], *alternative["alternative_route_ids"]]
    ]
    return {
        **row,
        "alternative_route_ids": list(route_ids),
        "alternative_family_candidates": family_alternatives,
    }


def _design_batches(designs, provider, deadline_monotonic=None, *, carried_only=()):
    # The scoring budget is global: carry every design before starting any new
    # projection. Otherwise early designs can displace later designs' old hits.
    for design in designs:
        yield design, [None]
    for design in designs:
        if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
            return
        if design in carried_only:
            continue
        proposals = (
            provider.propose(design, deadline_monotonic=deadline_monotonic)
            if deadline_monotonic is not None
            else provider.propose(design)
        )
        yield design, proposals


def run(hit, designs, config, output_dir, *, deadline_monotonic=None, scoring_mode="3d"):
    """Write complete evidence for each verified route, failure, and accepted pose.

    A fresh output directory is required. Configuration paths passed directly to
    this API are relative to the current directory; load_config resolves file paths.
    """
    out = Path(output_dir)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Synthesis output directory is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    summary = {
        "schema_version": 1,
        "status": "running",
        "proposals": 0,
        "verified_routes": 0,
        "scored_products": 0,
        "accepted_products": 0,
        "failures": [],
        "search_reports": [],
        "limits_reached": [],
        "deadline_reached": False,
        "scoring_mode": scoring_mode,
        "screening_allocation": [],
    }

    def deadline_reached():
        expired = deadline_monotonic is not None and time.monotonic() >= deadline_monotonic
        if expired:
            summary["deadline_reached"] = True
            if "deadline_monotonic" not in summary["limits_reached"]:
                summary["limits_reached"].append("deadline_monotonic")
        return expired

    write_json(out / "summary.json", summary)
    route_records, score_records, raw_records, refinement_records = [], [], [], []
    rows, poses, similarity_cache, candidates_by_key = {}, {}, {}, {}
    evidence_evaluator = None
    policy = None
    route_policy = None
    root_objectives = {}
    selected = []
    projected_designs = set()
    try:
        if deadline_monotonic is not None and (
            isinstance(deadline_monotonic, bool)
            or not isinstance(deadline_monotonic, (int, float))
            or not math.isfinite(deadline_monotonic)
        ):
            raise ValueError("deadline_monotonic must be a finite monotonic-clock timestamp")
        if scoring_mode not in {"3d", "2d"}:
            raise ValueError("scoring_mode must be 3d or 2d")
        config = _resolve_paths(_config(config), Path.cwd())
        policy = effective_policy(config)
        qualified = policy["mode"] == "qualified"
        route_policy = search_policy(config)
        summary["recommendation_mode"] = policy["mode"]
        if qualified and scoring_mode != "3d":
            raise ValueError("Qualified recommendations require completed 3D scoring")
        if scoring_mode == "2d":
            forbidden = set(config.get("minimum_scores", {})) - {"hit_tanimoto", "design_tanimoto"}
            if "pocket" in config or forbidden:
                raise ValueError(
                    "2D scoring cannot enforce pocket constraints or 3D minimum_scores; remove those explicitly for this ablation"
                )
        requested_hit = hit
        hit = canonical(hit)
        hit_fp = FP.GetFingerprint(molecule(hit))
        design_policy = selection_policy(config)
        bounded_refill = design_policy["version"] == "bounded-refill-v1"
        raw_designs = designs
        if bounded_refill:
            summary["designs_raw_supplied"] = sized_count(raw_designs)
        else:
            requested_designs = list(designs)
            all_designs = list(dict.fromkeys(canonical(s) for s in requested_designs))
            designs = all_designs[: config.get("max_designs", 20)]
            capped_designs = list(designs)
            if len(designs) < len(all_designs):
                summary["limits_reached"].append("max_designs")
            summary["designs_supplied"], summary["designs_used"] = len(all_designs), len(designs)
        if config.get("core_smiles") and not molecule(hit).HasSubstructMatch(
            molecule(config["core_smiles"]), useChirality=True
        ):
            raise ValueError("The original hit does not contain core_smiles")
        limits = SearchLimits(**config.get("search", {}))
        settings = ConformerSettings(**config.get("conformers", {}))
        requested_refinement = RefinementSettings(**config.get("refinement", {}))
        refinement = replace(requested_refinement, top_n=0) if scoring_mode == "2d" else requested_refinement

        def qualify(row, records, result, pose_settings, report=None):
            # Reuse the exact final-product chemistry checks at the export
            # boundary, before family selection or an alternative-design break.
            constraint_error = None
            try:
                check_constraints(
                    row["smiles"],
                    hit,
                    row["design"],
                    "hop",
                    config.get("core_smiles"),
                    config.get("preserve_design_scaffold", False),
                )
            except ValueError as exc:
                constraint_error = str(exc)
            final_fp = FP.GetFingerprint(molecule(row["smiles"]))
            row["hit_tanimoto"] = row["hit_morgan"] = DataStructs.TanimotoSimilarity(hit_fp, final_fp)
            row["design_tanimoto"] = DataStructs.TanimotoSimilarity(
                FP.GetFingerprint(molecule(row["design"])), final_fp
            )
            passed = []
            for record in records:
                decision = decide(row, record, result, policy, pose_settings, refinement, report, hit=hit)
                if constraint_error:
                    decision.update(status="rejected", qualified=False)
                    decision["reasons"].append(constraint_error)
                record["qualification"] = decision
                record["selection_status"] = (
                    "qualified" if decision["qualified"] else "qualification_rejected"
                )
                if decision["qualified"]:
                    passed.append(record)
            if not passed:
                return None
            primary = min(passed, key=lambda r: (len(r["verification"]["steps"]), r["route_id"]))
            return {
                **row,
                "route_id": primary["route_id"],
                "reaction_steps": len(primary["verification"]["steps"]),
                "alternative_route_ids": [r["route_id"] for r in passed if r is not primary],
                "alternative_family_candidates": [],
                "recommendation_status": "qualified",
                "qualification_policy": policy["policy_version"],
            }

        hit_pose = None
        try:
            if bounded_refill and deadline_reached():
                summary["query_preflight"] = {
                    "input": requested_hit,
                    "canonical_smiles": hit,
                    "role": "original_hit",
                    "status": "not_assessed",
                    "reason_code": "deadline",
                }
            elif scoring_mode == "3d":
                if bounded_refill or designs:
                    hit_pose = prepare_reference(hit, settings, config.get("reference_sdf"))
                elif config.get("reference_sdf"):
                    hit_pose = reference_from_sdf(config["reference_sdf"], hit, settings)
        except ValueError as exc:
            summary["query_preflight"] = {
                "input": requested_hit,
                "canonical_smiles": hit,
                "role": "original_hit",
                "status": "ineligible",
                "qualification_status": "not_assessed",
                "reason_code": "unspecified_stereochemistry"
                if "specified stereochemistry" in str(exc)
                else "reference_ineligible",
                "reason": str(exc),
                "interpretation": "Query could not be assessed; this is not evidence that no qualifying chemistry exists",
            }
            raise
        if bounded_refill:
            summary.setdefault(
                "query_preflight",
                {
                    "input": requested_hit,
                    "canonical_smiles": hit,
                    "role": "original_hit",
                    "status": "eligible" if scoring_mode == "3d" else "not_applicable_2d",
                    "stereo_policy": "Exact input state retained; no stereoisomer assignment or enumeration",
                },
            )
        if bounded_refill:
            design_preflight = select_designs(
                raw_designs,
                design_policy,
                lambda design: prepare_reference(design, settings) if scoring_mode == "3d" else None,
                deadline_reached,
            )
            requested_designs = design_preflight["requested_designs"]
            all_designs = capped_designs = design_preflight["canonical_unique_requested_designs"]
            design_checks = design_preflight["checks"]
            designs = design_preflight["eligible_designs"]
            summary["design_selection"] = {
                key: design_preflight[key] for key in ("policy", "counts", "stop_reason", "input_exhausted")
            }
            summary["designs_supplied"] = design_preflight["counts"]["raw_supplied"]
            summary["designs_used"] = len(capped_designs)
            summary["failures"].extend(dict(c) for c in design_checks if c["status"] == "rejected")
            for limit in design_preflight["limits_reached"]:
                if limit not in summary["limits_reached"]:
                    summary["limits_reached"].append(limit)
        else:
            design_checks, eligible_designs = [], []
            for index, design in enumerate(capped_designs):
                check = {"design": design, "capped_design_index": index, "role": "design"}
                if scoring_mode == "3d":
                    try:
                        prepare_reference(design, settings)
                    except ValueError as exc:
                        check.update(status="rejected", stage="design_reference_preflight", reason=str(exc))
                        summary["failures"].append(dict(check))
                        design_checks.append(check)
                        continue
                check["status"] = "eligible" if scoring_mode == "3d" else "not_applicable_2d"
                design_checks.append(check)
                eligible_designs.append(design)
            designs = eligible_designs
            design_preflight = {
                "requested_designs": requested_designs,
                "canonical_unique_requested_designs": all_designs,
                "capped_designs": capped_designs,
                "eligible_designs": designs,
                "checks": design_checks,
                "scope": "Configured generated-reference stereochemistry and alignment prerequisites only; no embedding or stereo assignment. Applied after the original design cap, with no backfill or provider reseeding.",
                "policy": design_policy,
            }
        rejected_designs = {
            check["design"]: check
            for check in design_checks
            if check["status"] == "rejected" and "design" in check
        }
        summary["designs_preflight_rejected"] = len(rejected_designs)
        summary["designs_eligible_for_scoring"] = len(designs)
        # Keep original capped inputs for continuation and preserve carried evidence
        # even when its design reference is unusable under today's scoring policy.
        carried, carry_sources = load_carried_routes(
            config.get("carry_forward", []),
            hit,
            capped_designs,
            **({"include_excluded": True} if bounded_refill else {}),
        )
        summary["carried_routes_loaded"] = sum(len(items) for items in carried.values())
        if bounded_refill:
            # Only already-saved raw prefixes and known carried design identities
            # can extend a later continuation. Never infer an unread input tail.
            supply, known = list(requested_designs), set(all_designs)
            for source in carry_sources:
                prior = source.get("prior_design_preflight") or {}
                saved = prior.get("saved_design_supply", prior.get("requested_designs", []))
                if len(saved) > len(supply) and saved[: len(supply)] == supply:
                    supply = list(saved)
                    known.update(
                        prior.get(
                            "saved_design_identities", prior.get("canonical_unique_requested_designs", [])
                        )
                    )
            for design in carried:
                if design not in known:
                    supply.append(design)
                    known.add(design)
            design_preflight["saved_design_supply"] = supply
            design_preflight["saved_design_identities"] = sorted(known)
        excluded_designs = (
            {
                design: rejected_designs.get(
                    design,
                    {
                        "design": design,
                        "role": "design",
                        "status": "not_assessed",
                        "reason": "Carried design was not selected within this run's raw/eligible/deadline bounds",
                    },
                )
                for design in carried
                if design not in designs
            }
            if bounded_refill
            else rejected_designs
        )
        for design, check in excluded_designs.items():
            for item in carried[design]:
                route_records.append(
                    {
                        "route_id": f"route_{len(route_records):06d}",
                        "proposal_id": None,
                        "design": design,
                        **item,
                        "selection_status": "design_preflight_rejected"
                        if check["status"] == "rejected"
                        else "design_not_assessed",
                        "reference_preflight": dict(check),
                        "reason": check["reason"],
                        "verification_status": "not_evaluated",
                        "eligible_for_evidence": False,
                        "eligible_for_scoring": False,
                    }
                )
        all_designs_rejected = not designs and (bounded_refill or bool(capped_designs))
        catalog = rules = None
        if all_designs_rejected:
            summary["no_candidates_reason"] = (
                "No eligible design was assessed within the configured intake limits; this does not establish that no qualifying chemistry exists"
                if bounded_refill
                else "All capped designs failed the configured 3D reference preflight"
            )
            summary["assets_not_loaded"] = ["catalog", "rules", "provider", "synthesis_evidence"]
            if "pocket" in config:
                summary["assets_not_loaded"].append("pocket")
        else:
            catalog, rules = Catalog(config["catalog"]), Rules(config["rules"])
            if deadline_reached():
                evidence_settings = {
                    key: value
                    for key, value in config.get("synthesis_evidence", {}).items()
                    if key != "predictor"
                }
                evidence_evaluator = evaluator_from_config(evidence_settings)
                summary["evidence_model_not_loaded_deadline"] = True
            else:
                evidence_evaluator = _evidence_evaluator(config)
        manifest = _manifest(config, hit, capped_designs)
        manifest["recommendation_policy"] = policy
        manifest["route_search_policy"] = route_policy
        manifest["design_preflight"] = {
            **design_preflight,
            "reject_unspecified_stereo": settings.reject_unspecified_stereo,
        }
        manifest["assets_not_loaded"] = summary.get("assets_not_loaded", [])
        manifest["scoring_backend"] = {
            "mode": scoring_mode,
            "shared_candidate_path": "After mode-specific reference eligibility, route collection, chemical/property constraints, synthesis evidence, allocation and screening cap use the same implementation",
            "requested_refinement": asdict(requested_refinement),
            "effective_refinement": asdict(refinement),
            "missing_3d_values": "null; no conformers or poses generated" if scoring_mode == "2d" else None,
        }
        if scoring_mode == "2d":
            manifest["ranking_scope"] = (
                "2D ablation: Pareto hit/design Morgan Tanimoto and reaction steps, then geometric mean similarity"
            )
            manifest["ranking_metrics"] = ["hit_tanimoto", "design_tanimoto", "reaction_steps"]
            manifest["refinement_policy"] = "Disabled because the complete 3D evaluation component is removed"
        manifest["deadline_policy"] = {
            "enabled": deadline_monotonic is not None,
            "clock": "time.monotonic",
            "seconds_available_at_start": max(0.0, deadline_monotonic - start)
            if deadline_monotonic is not None
            else None,
            "scope": "Cooperative boundary checks; complete evaluated candidates retained, incomplete seed blocks rejected; single long calls still require an outer hard timeout",
        }
        manifest["synthesis_evidence"] = (
            {
                **evidence_evaluator.provenance,
                "allocation": "One whole route per balanced product/design round before alternate routes; shared input-only prediction cache and explicit global call/token budgets",
            }
            if evidence_evaluator is not None
            else {"status": "not_initialized", "reason": summary["no_candidates_reason"]}
        )
        manifest["catalog_sha256"] = catalog.sha256 if catalog is not None else None
        manifest["rules_sha256"] = rules.sha256 if rules is not None else None
        manifest["catalog_import"] = catalog.import_provenance if catalog is not None else None
        manifest["carry_forward"] = carry_sources
        manifest["effective_search"], manifest["effective_conformers"] = asdict(limits), asdict(settings)
        manifest["effective_workflow"] = {
            "allocation_policy": config.get("allocation_policy", "balanced"),
            "mode": config.get("mode", "hop"),
            "core_smiles": config.get("core_smiles"),
            "preserve_design_scaffold": config.get("preserve_design_scaffold", False),
            "max_designs": config.get("max_designs", 20),
            "design_selection": design_policy,
            "max_scored_products": config.get("max_scored_products", 256),
            "max_results": config.get("max_results", 50),
            "expansion": {
                "rounds": 0,
                "max_routes": 64,
                "building_block_top_k": limits.building_block_top_k,
                "min_tanimoto": 0.0,
                "locked_building_block_ids": [],
                **config.get("expansion", {}),
            },
            "minimum_scores": config.get("minimum_scores", {}),
            "property_thresholds": config.get("property_thresholds", {}),
            "refinement": asdict(refinement),
        }
        pose_constraints = (
            constraints_from_config(config["pocket"], hit_pose, config.get("reference_sdf"))
            if "pocket" in config and not all_designs_rejected
            else None
        )
        scoring_options = {"pose_constraints": pose_constraints} if pose_constraints is not None else {}
        manifest["pocket"] = pose_constraints.provenance if pose_constraints is not None else None
        if config.get("reference_sdf"):
            manifest["reference_sdf_sha256"] = digest(config["reference_sdf"])
        provider = _provider(config["provider"]) if designs and not deadline_reached() else None
        carried_only = {
            design
            for design in designs
            if bounded_refill
            and config.get("carry_forward")
            and isinstance(provider, SavedProposals)
            and design not in provider.proposals
            and carried[design]
        }
        if carried_only:
            summary["carried_only_designs"] = [design for design in designs if design in carried_only]
            manifest["carried_only_policy"] = (
                "No saved proposal for this design in the latest archive; reverify and score available carried graphs without new provider calls or fabricated sampling"
            )
        manifest["provider"] = provider.provenance if provider else None
        write_json(out / "manifest.json", manifest)
        write_json(
            out / "proposals.json", {"schema_version": 1, "provider": manifest["provider"], "proposals": []}
        )
        from utils import _get_molecular_prop_extd_, determ_props

        def check_final_product(product, design):
            check_constraints(
                product,
                hit,
                design,
                config.get("mode", "hop"),
                config.get("core_smiles"),
                config.get("preserve_design_scaffold", False),
            )
            mol = molecule(product)
            if qualified:
                _check_stereo(mol, ConformerSettings())
            props = _get_molecular_prop_extd_(mol)
            accepted, reasons = determ_props(
                props, config.get("property_thresholds", {}), w_all_failed_reason=True
            )
            if not accepted:
                raise ValueError("Product property filter: " + "; ".join(reasons))
            return props

        for design, proposals in _design_batches(
            designs, provider, deadline_monotonic, carried_only=carried_only
        ):
            LOGGER.info("Projecting design %s onto catalog reaction space", design)
            design_fp = FP.GetFingerprint(molecule(design))
            if route_policy["version"] == "original-hit-v1" and design not in root_objectives:
                root_objectives[design] = RootProductObjective(
                    hit,
                    design,
                    config.get("minimum_scores", {}),
                    lambda product, design=design: check_final_product(product, design),
                    route_policy,
                )
            root_objective = root_objectives.get(design)
            search_options = {"root_objective": root_objective} if root_objective is not None else {}
            expansion_cache = {}
            # Carry evidence first, before the fresh search can consume the scoring budget.
            for proposal in proposals:
                if proposal is None:
                    proposal_id, projected = None, []
                    candidates = carried[design]
                else:
                    projected_designs.add(design)
                    proposal_id = len(raw_records)
                    raw_records.append({"design": design, **proposal})
                    write_json(
                        out / "proposals.json",
                        {"schema_version": 1, "provider": provider.provenance, "proposals": raw_records},
                    )
                    summary["proposals"] += 1
                    if deadline_reached():
                        summary["failures"].append(
                            {
                                "stage": "projection",
                                "proposal_id": proposal_id,
                                "status": "not_evaluated_deadline",
                                "reason": "Deadline reached before route projection",
                            }
                        )
                        continue
                    try:
                        projected, report = resolve_proposal(
                            proposal["response"], catalog, rules, limits, **search_options
                        )
                        summary["search_reports"].append(
                            {"proposal_id": proposal_id, "stage": "projection", **asdict(report)}
                        )
                        if not projected:
                            summary["failures"].append(
                                {
                                    "stage": "projection",
                                    "proposal_id": proposal_id,
                                    "reason": "No catalog-supported reaction route survived",
                                }
                            )
                    except (ValueError, KeyError, TypeError) as exc:
                        summary["failures"].append(
                            {"stage": "proposal", "proposal_id": proposal_id, "reason": str(exc)}
                        )
                        continue
                    candidates = list(projected)
                expansion = config.get("expansion", {})
                if expansion.get("rounds", 0) and not deadline_reached():
                    for item in projected:
                        if deadline_reached():
                            break
                        route_key = json.dumps(item["route"], sort_keys=True)
                        cached = route_key in expansion_cache
                        if not cached:
                            expansion_cache[route_key] = expand_route(
                                item["route"],
                                catalog,
                                rules,
                                limits,
                                target_smiles=design,
                                **expansion,
                                **search_options,
                            )
                        expanded, report = expansion_cache[route_key]
                        candidates.extend(expanded)
                        summary["search_reports"].append(
                            {
                                "proposal_id": proposal_id,
                                "stage": "expansion",
                                "cached": cached,
                                **asdict(report),
                            }
                        )
                for item in candidates:
                    route_id = f"route_{len(route_records):06d}"
                    record = {
                        "route_id": route_id,
                        "proposal_id": proposal_id,
                        "design": design,
                        **item,
                        "selection_status": "pending",
                    }
                    route_records.append(record)
                    if deadline_reached():
                        record.update(
                            selection_status="not_verified_deadline",
                            reason="Deadline reached before route verification",
                        )
                        continue
                    # Recheck carried routes against today's catalog and rules as well.
                    try:
                        verification = verify_route(
                            item["route"], catalog, rules, limits.reaction_product_limit
                        )
                        if len(verification["steps"]) > limits.max_steps:
                            raise ValueError("Route exceeds search.max_steps")
                    except (ValueError, KeyError, TypeError) as exc:
                        record.update(selection_status="verification_rejected", reason=str(exc))
                        continue
                    record["verification"] = verification
                    product = verification["product"]
                    summary["verified_routes"] += 1
                    if record.get("search_eligibility", {}).get("eligible") is False:
                        record.update(
                            selection_status="search_seed_only",
                            eligible_for_evidence=False,
                            eligible_for_scoring=False,
                            reason="Root retained only as a bounded expansion seed; failed current cheap final-product requirements",
                        )
                        continue
                    try:
                        props = check_final_product(product, design)
                    except ValueError as exc:
                        record.update(selection_status="constraint_rejected", reason=str(exc))
                        continue
                    key = (product, design)
                    if key not in similarity_cache:
                        product_fp = FP.GetFingerprint(molecule(product))
                        similarity_cache[key] = {
                            "hit_tanimoto": DataStructs.TanimotoSimilarity(hit_fp, product_fp),
                            "design_tanimoto": DataStructs.TanimotoSimilarity(design_fp, product_fp),
                        }
                    tanimoto = similarity_cache[key]
                    record["tanimoto"] = tanimoto
                    if any(
                        tanimoto[k] < v for k, v in config.get("minimum_scores", {}).items() if k in tanimoto
                    ):
                        record.update(selection_status="2d_threshold_rejected", scores=tanimoto)
                        continue
                    record["eligible_for_evidence"] = True
                    entry = candidates_by_key.setdefault(
                        key,
                        {
                            "product": product,
                            "design": design,
                            "properties": props,
                            "tanimoto": tanimoto,
                            "records": [],
                        },
                    )
                    entry["records"].append(record)
        # Collect before allocating expensive 3D work so one early design cannot
        # exhaust the global budget. Routes remain fully auditable when unscored.
        ordered = _candidate_order(
            list(candidates_by_key.values()), config.get("allocation_policy", "balanced")
        )
        # Independent model work is also shared fairly: try one route per
        # product/design before a prolific early candidate's alternative routes.
        evidence_queue = deque(
            (
                entry,
                deque(
                    sorted(
                        entry["records"],
                        key=lambda record: (len(record["verification"]["steps"]), record["route_id"]),
                    )
                ),
            )
            for entry in ordered
        )
        while evidence_queue:
            entry, pending_routes = evidence_queue.popleft()
            record = pending_routes.popleft()
            if deadline_reached():
                record.update(
                    selection_status="not_evaluated_deadline",
                    reason="Deadline reached before independent evidence evaluation",
                )
                record["synthesis_evidence"] = {
                    "schema_version": 1,
                    "status": "not_evaluated",
                    "reason": "Deadline reached",
                    "eligibility": {
                        "passed": False,
                        "reasons": ["Independent evidence was not evaluated before deadline"],
                    },
                }
                if pending_routes:
                    evidence_queue.append((entry, pending_routes))
                continue
            evidence = evidence_evaluator.evaluate(
                record["verification"],
                **({"deadline_monotonic": deadline_monotonic} if deadline_monotonic is not None else {}),
            )
            record["synthesis_evidence"] = evidence
            if evidence["eligibility"]["passed"]:
                record["eligible_for_scoring"] = True
            else:
                record.update(
                    selection_status="evidence_rejected", reasons=evidence["eligibility"]["reasons"]
                )
            if pending_routes:
                evidence_queue.append((entry, pending_routes))
        for entry in ordered:
            product, design = entry["product"], entry["design"]
            records = [record for record in entry["records"] if record.get("eligible_for_scoring")]
            tanimoto = entry["tanimoto"]
            if not records:
                continue
            if deadline_reached():
                for record in records:
                    record.update(
                        selection_status="not_scored_deadline",
                        reason="Deadline reached before candidate screening",
                    )
                continue
            if summary["scored_products"] >= config.get("max_scored_products", 256):
                for record in records:
                    record["selection_status"] = "not_scored_budget"
                if "max_scored_products" not in summary["limits_reached"]:
                    summary["limits_reached"].append("max_scored_products")
                continue
            summary["scored_products"] += 1
            summary["screening_allocation"].append(
                {
                    "product": product,
                    "design": design,
                    "route_ids": [record["route_id"] for record in records],
                }
            )
            if scoring_mode == "3d":
                try:
                    result = score_candidate(product, hit, design, settings, hit_pose, **scoring_options)
                    if qualified:
                        validate_result(result, product, settings, hit=hit, design=design)
                except (ValueError, RuntimeError) as exc:
                    for record in records:
                        record.update(selection_status="3d_failed", reason=str(exc))
                    continue
                info, hit_aligned, design_aligned = result
                scores = {
                    f"{ref}_{metric}": info[ref][metric]
                    for ref in ("hit", "design")
                    for metric in ("shape", "pharmacophore")
                }
            else:
                info = {
                    "score": math.sqrt(tanimoto["hit_tanimoto"] * tanimoto["design_tanimoto"]),
                    "scoring_mode": "2d",
                    "candidate_relative_energy": None,
                }
                hit_aligned = design_aligned = None
                scores = {
                    f"{ref}_{metric}": None
                    for ref in ("hit", "design")
                    for metric in ("shape", "pharmacophore")
                }
            scores.update(tanimoto)
            if not refinement.top_n and any(
                scores[k] < v for k, v in config.get("minimum_scores", {}).items()
            ):
                for record in records:
                    record.update(selection_status="3d_threshold_rejected", scores=scores)
                continue
            for record in records:
                record["selection_status"] = "screened" if refinement.top_n else "accepted"
                score_records.append(
                    {
                        "route_id": record["route_id"],
                        "product": product,
                        "design": design,
                        "properties": entry["properties"],
                        "tanimoto": tanimoto,
                        "scoring_stage": "screening" if scoring_mode == "3d" else "2d_screening",
                        **info,
                        "hit_smiles": hit,
                        "design_smiles": design,
                    }
                )
            primary = min(
                records, key=lambda record: (len(record["verification"]["steps"]), record["route_id"])
            )
            row = {
                "smiles": product,
                "design": design,
                "score": info["score"],
                **scores,
                "hit_morgan": tanimoto["hit_tanimoto"],
                "relative_energy": info["candidate_relative_energy"],
                "reaction_steps": len(primary["verification"]["steps"]),
                "route_id": primary["route_id"],
                "alternative_route_ids": [record["route_id"] for record in records if record is not primary],
                "scoring_stage": "screening" if scoring_mode == "3d" else "2d_screening",
                "pose_seed": settings.seed if scoring_mode == "3d" else None,
                "seed_score_median": info["score"] if scoring_mode == "3d" else None,
                "seed_score_min": info["score"] if scoring_mode == "3d" else None,
                "seed_score_max": info["score"] if scoring_mode == "3d" else None,
                "seed_successes": 1 if scoring_mode == "3d" else 0,
                "seed_count": 1 if scoring_mode == "3d" else 0,
                "family_smiles": family_key(product),
                "scaffold_family_smiles": scaffold_family_key(product),
                "ranking_score": info["score"],
                "seed_score_range": 0.0 if scoring_mode == "3d" else None,
                "seed_score_std": 0.0 if scoring_mode == "3d" else None,
                "convergence_status": "not_assessed" if scoring_mode == "3d" else "not_applicable_2d",
                "refinement_conformers": settings.num_conformers if scoring_mode == "3d" else 0,
                "recommendation_status": "exploratory_verified",
                "qualification_policy": None,
            }
            if qualified and not refinement.top_n:
                row = qualify(row, records, result, settings)
                if row is None:
                    continue
            rows[(product, design)] = row
            if scoring_mode == "3d":
                poses[(product, design)] = (info, hit_aligned, design_aligned)
        all_screened = _rank(list(rows.values()))
        _write_table(out / "screening_alternatives.tsv", all_screened)
        groups = _family_groups(all_screened)
        ranked = _rank(
            [_with_alternative_routes(alternatives[0], alternatives) for alternatives in groups.values()]
        )
        _write_table(out / "screening.tsv", ranked)
        summary["screened_products"] = len(ranked)
        summary["screened_product_design_pairs"] = len(all_screened)
        summary["screened_scaffold_families"] = len({row["scaffold_family_smiles"] for row in ranked})
        if refinement.top_n:
            shortlist = ranked
            if config.get("allocation_policy", "balanced") == "balanced":
                shortlist = _balanced_order(
                    shortlist, product=lambda row: row["smiles"], design=lambda row: row["design"]
                )
            shortlist = shortlist[: refinement.top_n]
            if len(shortlist) < len(ranked):
                summary["limits_reached"].append("refinement.top_n")
            ranked = []
            records_by_id = {record["route_id"]: record for record in route_records}
            budget = RefinementBudget(refinement.max_total_conformers)
            attempted = 0
            base_cost = refinement.num_conformers * len(refinement.seeds)
            requested_shortlist_count = len(shortlist)
            if budget.limit is not None and len(shortlist) * base_cost > budget.limit:
                # Budget permits a balanced prefix only; record every omitted family.
                eligible_count = budget.limit // base_cost
                for omitted in shortlist[eligible_count:]:
                    for coarse in groups[omitted["family_smiles"]]:
                        refinement_records.append(
                            {
                                "route_id": coarse["route_id"],
                                "product": coarse["smiles"],
                                "design": coarse["design"],
                                "family_smiles": coarse["family_smiles"],
                                "status": "budget_exhausted",
                                "convergence_status": "not_evaluated",
                                "reason": "Global budget cannot fund a complete initial seed block",
                                "trials": [],
                                "candidate_conformer_requests": 0,
                            }
                        )
                        for route_id in [coarse["route_id"], *coarse["alternative_route_ids"]]:
                            records_by_id[route_id]["selection_status"] = "not_refined_budget"
                shortlist = shortlist[:eligible_count]
                summary["limits_reached"].append("refinement.max_total_conformers")
            for family_index, primary in enumerate(shortlist):
                # Reserve a complete initial block for every remaining family.
                # Adaptive escalation and fallbacks can only spend unreserved work.
                budget.reserved_for_others = (len(shortlist) - family_index - 1) * base_cost
                alternatives = groups[primary["family_smiles"]]
                # A failed/threshold-rejected design cannot erase another valid
                # route to the same chemical family. Try alternatives deterministically.
                for alternative_index, coarse in enumerate(alternatives):
                    attempted += 1
                    LOGGER.info(
                        "Refining %s against %s (alternative %s)",
                        coarse["smiles"],
                        coarse["design"],
                        alternative_index,
                    )
                    result, report = refine_candidate(
                        coarse["smiles"],
                        hit,
                        coarse["design"],
                        settings,
                        refinement,
                        hit_pose,
                        minimum_scores=config.get("minimum_scores", {}),
                        budget=budget,
                        **(
                            {
                                "result_validator": lambda result, product, settings: validate_result(
                                    result, product, settings, hit=hit, design=coarse["design"]
                                )
                            }
                            if qualified
                            else {}
                        ),
                        **scoring_options,
                        **(
                            {"deadline_monotonic": deadline_monotonic}
                            if deadline_monotonic is not None
                            else {}
                        ),
                    )
                    if report["status"] == "deadline_exhausted":
                        deadline_reached()
                    evidence = {
                        "route_id": coarse["route_id"],
                        "product": coarse["smiles"],
                        "design": coarse["design"],
                        "family_smiles": coarse["family_smiles"],
                        "alternative_index": alternative_index,
                        **report,
                    }
                    refinement_records.append(evidence)
                    related = [
                        records_by_id[route_id]
                        for route_id in [coarse["route_id"], *coarse["alternative_route_ids"]]
                    ]
                    if result is None:
                        status = {
                            "threshold_rejected": "3d_threshold_rejected",
                            "unstable": "refinement_unstable",
                            "budget_exhausted": "not_refined_budget",
                            "deadline_exhausted": "not_refined_deadline",
                        }.get(report["status"], "refinement_failed")
                        for record in related:
                            record.update(selection_status=status, refinement_status=report["status"])
                        if (
                            report["status"] == "budget_exhausted"
                            and "refinement.max_total_conformers" not in summary["limits_reached"]
                        ):
                            summary["limits_reached"].append("refinement.max_total_conformers")
                        continue
                    info, _, _ = result
                    scores = {
                        f"{ref}_{metric}": info[ref][metric]
                        for ref in ("hit", "design")
                        for metric in ("shape", "pharmacophore")
                    }
                    row = {
                        **_with_alternative_routes(coarse, alternatives),
                        **scores,
                        "score": info["score"],
                        "relative_energy": info["candidate_relative_energy"],
                        "scoring_stage": "refined",
                        "refinement_conformers": info["settings"]["num_conformers"],
                        **{
                            key: report[key]
                            for key in (
                                "pose_seed",
                                "seed_score_median",
                                "seed_score_min",
                                "seed_score_max",
                                "seed_successes",
                                "seed_count",
                                "seed_score_range",
                                "seed_score_std",
                                "ranking_score",
                                "convergence_status",
                            )
                        },
                    }
                    if refinement.ranking_statistic == "minimum":
                        row["ranking_metrics"] = {
                            key: values["max" if key == "relative_energy" else "min"]
                            for key, values in report["metric_ranges"].items()
                        }
                    if any(row[k] < v for k, v in config.get("minimum_scores", {}).items()):
                        evidence["status"] = "3d_threshold_rejected"
                        for record in related:
                            record.update(selection_status="3d_threshold_rejected", scores=scores)
                        continue
                    if qualified:
                        row = qualify(row, related, result, ConformerSettings(**info["settings"]), report)
                        if row is None:
                            continue
                    for record in related:
                        if not qualified:
                            record["selection_status"] = "accepted"
                    ranked.append(row)
                    poses[(row["smiles"], row["design"])] = result
                    score_records.append(
                        {
                            "route_id": row["route_id"],
                            "product": row["smiles"],
                            "design": row["design"],
                            "scoring_stage": "refined",
                            **info,
                            "hit_smiles": hit,
                            "design_smiles": row["design"],
                        }
                    )
                    break
            ranked = _rank(ranked)
            summary["refinement_shortlisted_families"] = requested_shortlist_count
            summary["refinement_budgeted_families"] = len(shortlist)
            summary["refinement_attempted"] = attempted
            summary["refinement_accepted"] = len(ranked)
            summary["refinement_budget"] = budget.snapshot()
        if len(ranked) > config.get("max_results", 50):
            summary["limits_reached"].append("max_results")
        selected = ranked[: config.get("max_results", 50)]
        with ExitStack() as stack:
            stream = stack.enter_context((out / "results.tsv").open("w", newline=""))
            if scoring_mode == "3d":
                hit_writer = stack.enter_context(Chem.SDWriter(str(out / "poses_hit.sdf")))
                design_writer = stack.enter_context(Chem.SDWriter(str(out / "poses_design.sdf")))
                hit_ref_writer = stack.enter_context(Chem.SDWriter(str(out / "references_hit.sdf")))
                design_ref_writer = stack.enter_context(Chem.SDWriter(str(out / "references_design.sdf")))
            writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for rank, row in enumerate(selected, 1):
                row["rank"] = rank
                writer.writerow(
                    row
                    | {
                        "alternative_route_ids": json.dumps(row["alternative_route_ids"]),
                        "alternative_family_candidates": json.dumps(
                            row.get("alternative_family_candidates", [])
                        ),
                        "ranking_metrics": json.dumps(row.get("ranking_metrics", {}), sort_keys=True),
                    }
                )
                if scoring_mode == "2d":
                    continue
                info, hit_aligned, design_aligned = poses[(row["smiles"], row["design"])]
                pose_settings = ConformerSettings(**info["settings"])
                for mol, sdf in ((hit_aligned, hit_writer), (design_aligned, design_writer)):
                    mol = Chem.Mol(mol)
                    mol.SetProp("_Name", row["route_id"])
                    mol.SetProp("smiles", row["smiles"])
                    mol.SetProp("score", str(row["score"]))
                    mol.SetProp("scoring_stage", row["scoring_stage"])
                    mol.SetIntProp("pose_seed", row["pose_seed"])
                    for metric in ("seed_score_median", "seed_score_min", "seed_score_max"):
                        mol.SetDoubleProp(metric, row[metric])
                    for metric in (
                        "hit_tanimoto",
                        "design_tanimoto",
                        "hit_shape",
                        "design_shape",
                        "hit_pharmacophore",
                        "design_pharmacophore",
                    ):
                        mol.SetDoubleProp(metric, row[metric])
                    mol.SetProp("verification", "catalog_and_rule_verified; experimental_validation=false")
                    mol.SetProp("recommendation_status", row["recommendation_status"])
                    mol.SetProp("qualification_policy", row["qualification_policy"] or "not_assessed")
                    sdf.write(mol, confId=info["candidate_conformer_id"])
                for role, reference, sdf in (
                    (
                        "hit",
                        hit_pose if hit_pose is not None else ensemble(hit, pose_settings),
                        hit_ref_writer,
                    ),
                    ("design", ensemble(row["design"], pose_settings), design_ref_writer),
                ):
                    mol = Chem.Mol(reference.mol)
                    mol.SetProp("_Name", row["route_id"])
                    mol.SetProp("reference_role", role)
                    cid = info[role]["reference_conformer_id"]
                    mol.SetIntProp("source_conformer_id", cid)
                    mol.SetIntProp("pose_seed", row["pose_seed"])
                    sdf.write(mol, confId=cid)
        deadline_reached()
        summary["designs_not_projected"] = [
            design for design in capped_designs if design not in projected_designs
        ]
        summary.update(
            status="success" if selected else "no_candidates",
            accepted_products=len(ranked),
            returned_products=len(selected),
        )
        summary["qualified_products"] = len(ranked) if qualified else 0
        if qualified and not selected:
            summary.setdefault(
                "no_candidates_reason",
                ("No candidate completed the qualified recommendation contract within configured limits"),
            )
        eligible = {r["verification"]["product"] for r in route_records if r.get("eligible_for_scoring")}
        carried_eligible = {
            r["verification"]["product"]
            for r in route_records
            if r.get("eligible_for_scoring") and r.get("carried_from")
        }
        summary.update(
            eligible_products=len(eligible),
            carried_eligible_products=len(carried_eligible),
            new_eligible_products=len(eligible - carried_eligible),
            carried_screened_products=len(carried_eligible & {key[0] for key in rows}),
            accepted_scaffold_families=len({row["scaffold_family_smiles"] for row in ranked}),
            returned_scaffold_families=len({row["scaffold_family_smiles"] for row in selected}),
        )
        return selected
    except Exception as exc:
        summary.update(status="error", error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        write_json(
            out / "search_objectives.json",
            {
                "schema_version": 1,
                "policy": route_policy,
                "design_objectives": [objective.snapshot() for objective in root_objectives.values()],
            },
        )
        for record in route_records:
            if "qualification" not in record:
                record["qualification"] = {
                    "policy_version": policy["policy_version"] if policy else None,
                    "status": "not_assessed"
                    if not policy or policy["mode"] == "exploratory"
                    else "not_qualified",
                    "qualified": False,
                    "reasons": [record.get("reason", record.get("selection_status", "not_evaluated"))],
                    "route_id": record["route_id"],
                }
        write_json(
            out / "qualification.json",
            {
                "schema_version": 1,
                "policy": policy,
                "decisions": [record["qualification"] for record in route_records],
                "selected_route_ids": [row["route_id"] for row in selected],
            },
        )
        write_json(out / "results.json", {"schema_version": 1, "policy": policy, "results": selected})
        summary["elapsed_seconds"] = time.monotonic() - start
        if evidence_evaluator is not None:
            summary["synthesis_evidence_usage"] = evidence_evaluator.usage
            summary["synthesis_evidence_evaluated_routes"] = sum(
                "synthesis_evidence" in record for record in route_records
            )
        summary["route_rejections"] = {
            status: sum(r["selection_status"] == status for r in route_records)
            for status in sorted({r["selection_status"] for r in route_records})
        }
        write_json(out / "routes.json", {"schema_version": 1, "routes": route_records})
        write_json(out / "scores.json", {"schema_version": 1, "scores": score_records})
        write_json(out / "refinement.json", {"schema_version": 1, "candidates": refinement_records})
        write_json(out / "summary.json", summary)
