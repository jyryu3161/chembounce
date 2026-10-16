"""Versioned final-output contract, separate from search and heuristic ranking."""

import math
from dataclasses import asdict

import numpy as np
from rdkit import Chem

from .common import canonical, molecule
from .conformers import ConformerSettings, _check_stereo
from .evidence import EvidencePolicy, evidence_requirements, matching_prediction_ranks
from .refinement import RefinementSettings
from .profiles import profile_provenance

VERSION = "qualified-v1"
REQUIRED_MINIMA = ("hit_tanimoto", "hit_shape", "hit_pharmacophore")


def effective_policy(config):
    """No empirical cutoff is supplied by this engineering contract."""
    request = config.get("recommendation", {})
    if not isinstance(request, dict) or set(request) - {"mode", "policy_version"}:
        raise ValueError("recommendation requires only mode and policy_version")
    mode = request.get("mode", "exploratory")
    if mode not in ("exploratory", "qualified") or request.get("policy_version", VERSION) != VERSION:
        raise ValueError("Unknown recommendation mode or policy_version")
    minimum = config.get("minimum_scores", {})
    if mode == "qualified":
        if any(not finite(minimum.get(key)) or not 0 < minimum[key] <= 1 for key in REQUIRED_MINIMA):
            raise ValueError(
                "Qualified recommendations require explicit positive finite minimum_scores for "
                + ", ".join(REQUIRED_MINIMA)
            )
        if config.get("mode", "hop") != "hop":
            raise ValueError("Qualified recommendations require mode=hop")
        if not config.get("conformers", {}).get("reject_unspecified_stereo", True):
            raise ValueError("Qualified recommendations require reject_unspecified_stereo=true")
        if config.get("refinement", {}).get("threshold_policy", "all_seeds") != "all_seeds":
            raise ValueError("Qualified recommendations require refinement.threshold_policy=all_seeds")
    return {
        "policy_version": VERSION,
        "mode": mode,
        "calibration": profile_provenance(config),
        "minimum_scores": dict(minimum),
        "evidence_profile": config.get("synthesis_evidence", {}).get("profile", "legacy-conditional-v1"),
        "evidence_policy": asdict(
            EvidencePolicy(
                **config.get("synthesis_evidence", {}).get(
                    "policy",
                    {"require_forward_agreement": True}
                    if config.get("synthesis_evidence", {}).get("profile") == "computational-route-support-v1"
                    else {},
                )
            )
        ),
        "evidence_product_matching": config.get("synthesis_evidence", {}).get("product_matching", "exact"),
        "conformers": asdict(ConformerSettings(**config.get("conformers", {}))),
        "ranking_policy": config.get(
            "ranking_policy", config.get("conformers", {}).get("scoring_objective", "joint-v1")
        ),
        "refinement": asdict(RefinementSettings(**config.get("refinement", {}))),
        "constraints": {
            key: config.get(key, default)
            for key, default in (
                ("mode", "hop"),
                ("core_smiles", None),
                ("preserve_design_scaffold", False),
                ("property_thresholds", {}),
            )
        },
        "reference": "original_hit_and_regenerated_final_product",
        "fingerprint": "binary Morgan radius=2 bits=2048 chirality=true; not sparse-count Morgan",
        "seed_acceptance": "every requested seed at the final ensemble size",
        "scope": "Explicit engineering thresholds or reviewed development fit as recorded in calibration provenance; requested computational evidence, no activity or synthesis probability",
    }


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _validate_score_info(info, settings):
    """Validate score protocol/energy metadata, including persisted seed reports."""
    try:
        if info["settings"] != asdict(settings):
            raise ValueError("Scoring settings do not match the requested seed/ensemble")
        if info.get("scoring_objective") != settings.scoring_objective:
            raise ValueError("Scoring objective does not match the requested protocol")
        cid = info["candidate_conformer_id"]
        if isinstance(cid, bool) or not isinstance(cid, int) or cid < 0:
            raise ValueError("Missing or invalid candidate conformer identity")
        for value in (info["score"], info["candidate_relative_energy"], info["candidate_energy"]):
            if not finite(value):
                raise ValueError("Missing or nonfinite scoring/energy value")
        for role in ("hit", "design"):
            reference_id = info[role]["reference_conformer_id"]
            if isinstance(reference_id, bool) or not isinstance(reference_id, int) or reference_id < 0:
                raise ValueError("Missing or invalid reference conformer identity")
            if not finite(info[role]["reference_relative_energy"]):
                raise ValueError("Missing or nonfinite reference energy")
            for metric in ("shape", "pharmacophore"):
                value = info[role][metric]
                if not finite(value) or not 0 <= value <= 1:
                    raise ValueError("Missing, nonfinite or out-of-range pose score")
        if settings.scoring_objective == "original-hit-v1":
            for count in info["ensemble_sizes"].values():
                if (
                    isinstance(count, bool)
                    or not isinstance(count, int)
                    or not 1 <= count <= settings.num_conformers
                ):
                    raise ValueError("Invalid realized ensemble size")
            if set(info["ensemble_sizes"]) != {"candidate", "hit", "design"}:
                raise ValueError("Incomplete realized ensemble sizes")
            for role in ("hit", "design"):
                overlay = info[role]
                if overlay["reference_relative_energy"] < 0:
                    raise ValueError("Negative reference relative energy")
                joint = math.sqrt(overlay["shape"] * overlay["pharmacophore"])
                adjusted = joint * math.exp(-settings.strain_penalty * overlay["reference_relative_energy"])
                if any(
                    not finite(overlay.get(key))
                    or not math.isclose(overlay[key], expected, rel_tol=1e-10, abs_tol=1e-12)
                    for key, expected in (("joint", joint), ("reference_adjusted_joint", adjusted))
                ):
                    raise ValueError("Overlay objective disagrees with pose metrics/reference energy")
            expected = info["hit"]["reference_adjusted_joint"] * math.exp(
                -settings.strain_penalty * info["candidate_relative_energy"]
            )
            if info["candidate_relative_energy"] < 0 or not math.isclose(
                info["score"], expected, rel_tol=1e-10, abs_tol=1e-12
            ):
                raise ValueError("Original-hit score disagrees with reference/candidate energy objective")
    except (KeyError, TypeError, AttributeError, RuntimeError) as exc:
        raise ValueError(f"Incomplete scoring metadata: {exc}") from exc


def validate_result(result, product, settings, *, hit=None, design=None):
    """Check the scorer's actual export poses and complete score payload."""
    try:
        info, hit_pose, design_pose = result
        for key, expected in (("candidate_smiles", product), ("hit_smiles", hit), ("design_smiles", design)):
            if expected is not None and info.get(key) != expected:
                raise ValueError(f"Scoring reference/product identity mismatch: {key}")
        _check_stereo(molecule(product), ConformerSettings())
        _validate_score_info(info, settings)
        for pose in (hit_pose, design_pose):
            if pose is None or canonical(Chem.MolToSmiles(pose)) != product:
                raise ValueError("Output pose does not match the exact final product")
            conformer = pose.GetConformer(info["candidate_conformer_id"])
            if not conformer.Is3D() or not np.isfinite(conformer.GetPositions()).all():
                raise ValueError("Output pose requires finite 3D coordinates")
    except (KeyError, TypeError, AttributeError, RuntimeError) as exc:
        raise ValueError(f"Incomplete final scoring result: {exc}") from exc


def decide(row, record, result, policy, settings, refinement, report=None, *, hit, evidence_validator=None):
    """Conjunction over final-route identity, requested evidence, scores and seeds.

    Chemical constraints use workflow.check_constraints before this call. Route
    replay and external evidence use the existing verifiers, not new chemistry.
    """
    reasons = []
    verification = record.get("verification", {})
    route = record.get("route", {})
    root = next((node for node in route.get("nodes", []) if node.get("id") == route.get("root")), {})
    try:
        root_product = canonical(root.get("product"))
    except ValueError:
        root_product = None
    if (
        verification.get("status") != "catalog_and_rule_verified"
        or verification.get("product") != row["smiles"]
        or not verification.get("steps")
        or root_product != row["smiles"]
    ):
        reasons.append("Exact final-product route verification is missing or mismatched")
    evidence = record.get("synthesis_evidence", {})
    if (
        evidence.get("product") != row["smiles"]
        or evidence.get("policy") != policy["evidence_policy"]
        or evidence.get("eligibility", {}).get("passed") is not True
        or len(evidence.get("steps", [])) != len(verification.get("steps", []))
    ):
        reasons.append("Requested synthesis evidence is missing, failed or mismatched")
    elif policy.get("evidence_profile") == "computational-route-support-v1":
        try:
            if evidence_validator is None or evidence_validator(record) is not True:
                raise ValueError("Strict evidence requires fresh route/raw prediction validation")
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            reasons.append(str(exc))
    else:
        for step, verified in zip(evidence["steps"], verification["steps"]):
            if (
                step.get("node_id") != verified["node_id"]
                or step.get("product") != verified["product"]
                or step.get("rule_replay") != "matched"
                or step.get("eligibility", {}).get("passed") is not True
            ):
                reasons.append("Route step evidence is incomplete or mismatched")
                break
            context = step.get("context", {}).get("supplied", {})
            requirements = evidence_requirements(
                step.get("forward_prediction", {}).get("status"),
                "reagents" in context and bool(context.get("conditions")),
                step.get("atom_sources", {}).get("status"),
                step.get("observations", {}).get("status"),
            )
            if any(required and not requirements[key] for key, required in policy["evidence_policy"].items()):
                reasons.append("A requested route evidence requirement was not established")
            if policy["evidence_policy"]["require_forward_agreement"]:
                prediction = step.get("forward_prediction", {})
                matches = matching_prediction_ranks(
                    prediction.get("predictions", []),
                    verified["product"],
                    policy["evidence_product_matching"],
                )
                if (
                    not matches
                    or prediction.get("matched_rank") != min(matches)
                    or prediction.get("matching_policy") != policy["evidence_product_matching"]
                ):
                    reasons.append(
                        "Required forward prediction lacks matching product evidence under the requested policy"
                    )
    try:
        validate_result(result, row["smiles"], settings, hit=hit, design=row["design"])
    except ValueError as exc:
        reasons.append(str(exc))
    info = result[0]
    if (
        row.get("pose_seed") != settings.seed
        or row.get("scoring_objective") != settings.scoring_objective
        or row.get("ranking_policy") != policy["ranking_policy"]
        or row.get("recommendation_profile") != (policy["calibration"]["name"] or "none")
        or row.get("calibration_status") != policy["calibration"]["status"]
        or row.get("calibration_artifact_sha256") != (policy["calibration"]["artifact_sha256"] or "none")
        or row.get("score") != info.get("score")
        or row.get("relative_energy") != info.get("candidate_relative_energy")
        or any(
            row.get(f"{role}_{metric}") != info.get(role, {}).get(metric)
            for role in ("hit", "design")
            for metric in ("shape", "pharmacophore")
        )
    ):
        reasons.append("Final row does not describe the exported pose scores and seed")
    minimum = policy["minimum_scores"]
    failed = [key for key, cutoff in minimum.items() if not finite(row.get(key)) or row[key] < cutoff]
    if failed:
        reasons.append("Final-product score minima failed: " + ", ".join(failed))
    if refinement.top_n:
        trials = (report or {}).get("trials", [])
        allowed_sizes = {refinement.num_conformers}
        size = refinement.num_conformers
        maximum = refinement.adaptive_max_conformers or size
        while size < maximum:
            size = min(maximum, size * refinement.conformer_growth_factor)
            allowed_sizes.add(size)
        expected_settings = {
            **policy["conformers"],
            "num_conformers": settings.num_conformers,
            "seed": settings.seed,
        }
        if (
            not report
            or report.get("status") != "success"
            or [trial.get("seed") for trial in trials] != list(refinement.seeds)
            or report.get("seed_successes") != len(refinement.seeds)
            or report.get("seed_count") != len(refinement.seeds)
            or report.get("pose_seed") != settings.seed
            or report.get("threshold_policy") != "all_seeds"
            or report.get("num_conformers") != settings.num_conformers
            or settings.num_conformers not in allowed_sizes
            or asdict(settings) != expected_settings
        ):
            reasons.append("Required final refinement seed block is incomplete")
        for trial in trials:
            try:
                _validate_score_info(
                    trial, ConformerSettings(**{**expected_settings, "seed": trial.get("seed")})
                )
            except (ValueError, TypeError) as exc:
                reasons.append(f"Required seed score metadata: {exc}")
            if (
                trial.get("status") != "success"
                or not finite(trial.get("score"))
                or trial.get("settings", {}).get("seed") != trial.get("seed")
                or trial.get("scoring_objective") != settings.scoring_objective
                or trial.get("candidate_smiles") != row["smiles"]
                or trial.get("hit_smiles") != hit
                or trial.get("design_smiles") != row["design"]
                or trial.get("settings", {}).get("num_conformers") != settings.num_conformers
                or trial.get("settings") != {**expected_settings, "seed": trial.get("seed")}
                or any(
                    not finite(trial.get(role, {}).get(metric))
                    or trial[role][metric] < minimum.get(f"{role}_{metric}", 0)
                    for role in ("hit", "design")
                    for metric in ("shape", "pharmacophore")
                )
            ):
                reasons.append("A required seed failed scoring or a score minimum")
                break
        if refinement.max_score_range is not None or refinement.max_metric_range is not None:
            if (report or {}).get("convergence_status") != "converged":
                reasons.append("Requested seed convergence was not established")
            values = [t.get("score") for t in trials]
            if (
                refinement.max_score_range is not None
                and values
                and all(finite(v) for v in values)
                and max(values) - min(values) > refinement.max_score_range
            ):
                reasons.append("Required seed score range exceeded")
            for role in ("hit", "design"):
                for metric in ("shape", "pharmacophore"):
                    values = [t.get(role, {}).get(metric) for t in trials]
                    if (
                        refinement.max_metric_range is not None
                        and values
                        and all(finite(v) for v in values)
                        and max(values) - min(values) > refinement.max_metric_range
                    ):
                        reasons.append("Required seed metric range exceeded")
        selected = next((t for t in trials if t.get("seed") == settings.seed), None)
        if selected is None or any(
            selected.get(k) != result[0].get(k)
            for k in ("score", "hit", "design", "candidate_conformer_id", "settings")
        ):
            reasons.append("Exported pose scores do not match the accepted representative seed")
    return {
        "policy_version": VERSION,
        "calibration": {key: policy["calibration"][key] for key in ("name", "status", "artifact_sha256")},
        "status": "qualified" if not reasons else "rejected",
        "qualified": not reasons,
        "reasons": reasons,
        "scores": {key: row.get(key) for key in minimum},
        "route_id": record["route_id"],
    }
