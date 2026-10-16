"""Reviewed development profile assets, effective defaults, and scope provenance."""

import copy
import hashlib
import importlib.metadata
import json
import math
import platform
from importlib.resources import files
from pathlib import Path

PROFILE = "qualified-development-v1"
# Reviewed development release; the untouched fit candidate remains separately bound.
PROFILE_SHA256 = {PROFILE: "9d3c60a895912b28c13a29aac0a04f385ddb8a0effc9fcd81d64a50a3d48bdd2"}
KERNEL_FILES = ("conformers.py", "common.py", "identity.py")
METRICS = ("hit_tanimoto", "hit_shape", "hit_pharmacophore")
FIXED_CONFORMERS = {
    "scoring_objective": "original-hit-v1",
    "num_conformers": 8,
    "seed": 42,
    "energy_window": 5.0,
    "strain_penalty": 0.05,
    "max_iterations": 500,
    "pharmacophore_sigma": 1.0,
    "reject_unspecified_stereo": True,
}
FIXED_REFINEMENT = {
    "top_n": 10,
    "num_conformers": 8,
    "seeds": [42, 43, 44],
    "threshold_policy": "all_seeds",
    "ranking_statistic": "minimum",
    "max_score_range": None,
    "max_metric_range": None,
    "adaptive_max_conformers": None,
    "conformer_growth_factor": 2,
    "borderline_margin": 0.0,
    "max_total_conformers": 240,
}


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _resource(name):
    return files("chembounce_synthesis").joinpath("data", name + ".json")


def load_profile(name):
    if not isinstance(name, str) or name not in PROFILE_SHA256:
        raise ValueError(f"Unknown recommendation_profile: {name!r}")
    expected = PROFILE_SHA256[name]
    if expected is None:
        raise ValueError(
            f"Recommendation profile {name} is unavailable: no reviewed feasible calibration artifact is packaged"
        )
    try:
        raw = _resource(name).read_bytes()
    except OSError as exc:
        raise ValueError(f"Missing packaged recommendation profile: {name}") from exc
    if _sha(raw) != expected:
        raise ValueError("Packaged recommendation profile SHA-256 mismatch")
    artifact = json.loads(raw)
    minima = artifact.get("minimum_scores")
    spec = artifact.get("calibration_spec", {})
    if (
        artifact.get("schema") != "qualified-development-reviewed-profile-v1"
        or artifact.get("status") != "reviewed_development_fit"
        or artifact.get("scope") != "previously_inspected_development_only"
        or not isinstance(minima, dict)
        or set(minima) != set(METRICS)
        or any(
            isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or not 0 < v <= 1
            for v in minima.values()
        )
        or spec.get("seeds") != [42, 43, 44]
        or spec.get("conformers") != {k: v for k, v in FIXED_CONFORMERS.items() if k != "seed"}
        or spec.get("fingerprint") != {"kind": "binary_morgan", "radius": 2, "bits": 2048, "chirality": True}
    ):
        raise ValueError(
            "Packaged recommendation profile lacks the reviewed fixed protocol and positive minima"
        )
    review = artifact.get("review", {})
    digests = [artifact.get("fit_sha256"), artifact.get("provenance", {}).get("procedure_sha256")]
    digests.extend(
        review.get(key) for key in ("candidate_sha256", "root_review_sha256", "independent_review_sha256")
    )
    if any(
        not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value)
        for value in digests
    ):
        raise ValueError(
            "Packaged recommendation profile requires candidate/fit/procedure and independent/root review digests"
        )
    return artifact


def profile_defaults(artifact):
    return {
        "recommendation": {"mode": "qualified", "policy_version": "qualified-v1"},
        "mode": "hop",
        "conformers": copy.deepcopy(FIXED_CONFORMERS),
        "refinement": copy.deepcopy(FIXED_REFINEMENT),
        "minimum_scores": dict(artifact["minimum_scores"]),
        "ranking_policy": "original-hit-v1",
        "route_search_policy": "original-hit-v1",
        "max_designs": 20,
        "max_design_attempts": 80,
        "max_scored_products": 256,
        "max_results": 10,
        "allocation_policy": "balanced",
        "synthesis_evidence": {
            "profile": "computational-route-support-v1",
            "predictor": {"kind": "reactiont5"},
            "policy": {"require_forward_agreement": True},
        },
    }


def apply_profile_defaults(config):
    """Apply before generic protocol defaults; explicit user values always survive."""
    result = copy.deepcopy(config)
    if "recommendation_profile" not in config:
        return result
    defaults = profile_defaults(load_profile(config["recommendation_profile"]))
    if not isinstance(result.get("synthesis_evidence", {}), dict):
        raise ValueError("synthesis_evidence must be an object")
    if (
        result.get("synthesis_evidence", {}).get("profile", "computational-route-support-v1")
        != "computational-route-support-v1"
    ):
        raise ValueError(
            "The named recommendation profile requires computational-route-support-v1; use an explicit engineering configuration for legacy evidence"
        )

    def fill(target, source):
        for key, value in source.items():
            if key not in target:
                target[key] = copy.deepcopy(value)
            elif isinstance(value, dict) and isinstance(target[key], dict):
                fill(target[key], value)

    fill(result, defaults)
    return result


def numerical_runtime():
    from rdkit.Chem import RDConfig

    feature_path = Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"
    return {
        "python": platform.python_version(),
        "packages": {name: importlib.metadata.version(name) for name in ("rdkit", "numpy", "scipy")},
        "pharmacophore_feature_definition": {
            "path": str(feature_path),
            "sha256": _sha(feature_path.read_bytes()),
        },
        "source_sha256": {
            "chembounce_synthesis/" + name: _sha(files("chembounce_synthesis").joinpath(name).read_bytes())
            for name in KERNEL_FILES
        },
    }


def _differences(expected, actual, prefix=""):
    differences = {}
    for key in sorted(set(expected) | set(actual)):
        name = prefix + key
        before, after = expected.get(key), actual.get(key)
        if isinstance(before, dict) and isinstance(after, dict):
            differences.update(_differences(before, after, name + "."))
        elif before != after:
            differences[name] = {"profile": before, "effective": after}
    return differences


def profile_provenance(config):
    name = config.get("recommendation_profile")
    if name is None:
        return {"name": None, "status": "not_requested", "artifact_sha256": None}
    artifact = load_profile(name)
    defaults = profile_defaults(artifact)
    # Route support is an additional selection constraint, outside the fitted numerical protocol.
    defaults.pop("synthesis_evidence", None)
    effective = {key: config.get(key) for key in defaults}
    differences = _differences(defaults, effective)
    resource_keys = {
        "max_designs",
        "max_design_attempts",
        "max_scored_products",
        "max_results",
        "allocation_policy",
        "refinement.top_n",
        "refinement.max_total_conformers",
    }
    numerical_changes = {key: value for key, value in differences.items() if key not in resource_keys}
    if config.get("reference_sdf"):
        numerical_changes["reference_sdf"] = {"profile": None, "effective": config["reference_sdf"]}
    if config.get("pocket"):
        numerical_changes["pocket"] = {"profile": None, "effective": config["pocket"]}
    # A profile match requires an actual complete qualified acceptance protocol.
    if config.get("refinement", {}).get("top_n", 0) <= 0:
        numerical_changes["refinement.top_n"] = {
            "profile": "positive",
            "effective": config.get("refinement", {}).get("top_n"),
        }
    budget = config.get("refinement", {}).get("max_total_conformers")
    if budget is None or budget < 24:
        numerical_changes["refinement.max_total_conformers"] = {"profile": "at least 24", "effective": budget}
    current = numerical_runtime()
    recorded = artifact["provenance"]
    mismatches = {}
    for field in ("python", "packages"):
        if current[field] != recorded["runtime"][field]:
            mismatches[field] = {"calibration": recorded["runtime"][field], "current": current[field]}
    feature = "pharmacophore_feature_definition"
    if current[feature]["sha256"] != recorded["runtime"][feature]["sha256"]:
        mismatches[feature] = {"calibration": recorded["runtime"][feature], "current": current[feature]}
    for path, digest in current["source_sha256"].items():
        if digest != recorded["source_sha256"].get(path):
            mismatches[path] = {"calibration": recorded["source_sha256"].get(path), "current": digest}
    return {
        "name": name,
        "status": "custom_unvalidated" if numerical_changes or mismatches else "development_protocol_matched",
        "artifact_sha256": PROFILE_SHA256[name],
        "fit_sha256": artifact["fit_sha256"],
        "review": artifact["review"],
        "scope": artifact["scope"],
        "calibration_spec": artifact["calibration_spec"],
        "calibration_provenance": recorded,
        "current_numerical_runtime": current,
        "material_overrides": numerical_changes,
        "runtime_source_mismatches": mismatches,
        "resource_overrides": {key: value for key, value in differences.items() if key in resource_keys},
        "additional_selection": {
            key: config.get(key)
            for key in (
                "core_smiles",
                "preserve_design_scaffold",
                "property_thresholds",
                "search",
                "expansion",
                "synthesis_evidence",
            )
            if key in config
        },
        "source_match_scope": ["chembounce_synthesis/" + name for name in KERNEL_FILES],
        "limitations": artifact.get("limitations", [])
        + [
            "Profile match describes the development scoring/acceptance protocol, not pipeline retention or biological activity",
            "Real-design diagnostic availability, anchors/properties, finite budgets and route evidence impose additional vetoes",
            "Observed three-seed minima do not establish future-seed robustness or convergence",
        ],
    }


def profile_label(config):
    if not config.get("recommendation_profile"):
        return (
            " [evidence=" + config.get("synthesis_evidence", {}).get("profile", "legacy-conditional-v1") + "]"
        )
    provenance = profile_provenance(config)
    return f" [profile={provenance['name']}; {provenance['status']}; computational-route-support-v1 conditional selection support; permissive development thresholds, not strong-similarity or activity guarantees]"
