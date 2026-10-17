"""Versioned development-fit and declared engineering profiles with provenance."""

import copy
import hashlib
import importlib.metadata
import json
import math
import platform
from importlib.resources import files
from pathlib import Path

PROFILE = "qualified-development-v1"
# Historical benchmark authority: never retarget this name.
SELECTIVE_PROFILE = "qualified-selective-engineering-v2"
GUARDED_SELECTIVE_PROFILE = "qualified-selective-engineering-v3"
DEFAULT_RECOMMENDATION_PROFILE = GUARDED_SELECTIVE_PROFILE
ENGINEERING_PROFILES = (SELECTIVE_PROFILE, GUARDED_SELECTIVE_PROFILE)
# Portable declarations bound by the packaged artifact hash. External release
# review authenticates these prior records; loading never opens historical paths.
GUARDED_SOURCE_REVISION = {
    "guarded_source_commit": "88c420f8bb428f658d174df1cebe12a46870550b",
    "predecessor": {
        "name": SELECTIVE_PROFILE,
        "artifact_sha256": "a8b47a409a21518dd482acbefdce38774414ae0b65f28d6ea7cb52ab63e1ef52",
    },
    "source_reviews": [
        {
            "role": "independent_source_review",
            "record_label": "review-molzip-assembly-source-88c420f/report.json",
            "record_sha256": "46346859a9d762c888156a178a439d0bfe715315877336406f9c2a8ab16d630f",
        },
        {
            "role": "root_source_gate",
            "record_label": "stage14-stereo-develop-merge-05ec327.json",
            "record_sha256": "a16f0137ee4b75c7cf3987155356d50c4dabce5eaf8d43d883bccd69c95d01f6",
        },
    ],
}
ENGINEERING_SCHEMA = "qualified-selective-engineering-profile-v2"
SELECTIVE_KERNEL_FILES = ("selective_scoring.py", "similarity_constraints.py")
SELECTIVE_CONSTRAINTS = {
    "version": "selective-similarity-constraints-v2",
    "polar_recall": 0.75,
    "landmarks": [],
}
# Reviewed development release; the untouched fit candidate remains separately bound.
PROFILE_SHA256 = {
    GUARDED_SELECTIVE_PROFILE: "fd5fb3b65b12307ed1a9b5371e74427eb709d8cf3012bbcba5095022ddb34537",
    SELECTIVE_PROFILE: "a8b47a409a21518dd482acbefdce38774414ae0b65f28d6ea7cb52ab63e1ef52",
    PROFILE: "9d3c60a895912b28c13a29aac0a04f385ddb8a0effc9fcd81d64a50a3d48bdd2",
}
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


def request_recommendations(config, cli_profile=None):
    """Select a fixed policy without concealing either explicit declaration."""
    if not isinstance(config, dict):
        raise ValueError("Recommendation configuration must be an object")
    for origin, value, present in (
        ("configuration", config.get("recommendation_profile"), "recommendation_profile" in config),
        ("CLI", cli_profile, cli_profile is not None),
    ):
        if present and value != DEFAULT_RECOMMENDATION_PROFILE:
            raise ValueError(
                f"--recommend/recommend() conflicts with {origin} recommendation_profile={value!r}; "
                f"requires {DEFAULT_RECOMMENDATION_PROFILE}"
            )
    result = copy.deepcopy(config)
    result["recommendation_profile"] = DEFAULT_RECOMMENDATION_PROFILE
    return result


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
    if name in ENGINEERING_PROFILES:
        _validate_engineering_artifact(artifact, name)
        return artifact
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
    defaults = {
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

    if artifact["schema"] == ENGINEERING_SCHEMA:
        defaults["similarity_constraints"] = copy.deepcopy(SELECTIVE_CONSTRAINTS)
    return defaults


def _validate_engineering_artifact(artifact, name):
    if (
        artifact.get("schema") != ENGINEERING_SCHEMA
        or artifact.get("name") != name
        or artifact.get("status") != "declared_engineering_policy"
        or artifact.get("scope") != "explicit_engineering_preferences_not_fitted"
        or artifact.get("minimum_scores")
        != {"hit_tanimoto": 0.35, "hit_shape": 0.70, "hit_pharmacophore": 0.50}
        or artifact.get("engineering_spec")
        != {
            "conformers": FIXED_CONFORMERS,
            "refinement": FIXED_REFINEMENT,
            "similarity_constraints": SELECTIVE_CONSTRAINTS,
            "pose_selection_policy": "joint-feasible-hit-pose-v2",
            "ranking_policy": "original-hit-v1",
            "route_search_policy": "original-hit-v1",
            "evidence_profile": "computational-route-support-v1",
            "fingerprint": {"kind": "binary_morgan", "radius": 2, "bits": 2048, "chirality": True},
        }
        or any(key in artifact for key in ("fit_sha256", "review", "calibration_spec"))
        or set(artifact.get("provenance", {}).get("source_sha256", {}))
        != {"chembounce_synthesis/" + name for name in KERNEL_FILES + SELECTIVE_KERNEL_FILES}
    ):
        raise ValueError("Packaged engineering profile lacks its declared fixed policy/source scope")
    if (
        name == GUARDED_SELECTIVE_PROFILE
        and artifact["provenance"].get("source_revision") != GUARDED_SOURCE_REVISION
    ):
        raise ValueError(
            "Packaged guarded engineering profile lacks its reviewed source revision declaration"
        )


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


def numerical_runtime(name=None):
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
            for name in (
                KERNEL_FILES + SELECTIVE_KERNEL_FILES if name in ENGINEERING_PROFILES else KERNEL_FILES
            )
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
    # Strict route support remains a mandatory additional selection constraint.
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
    if name not in ENGINEERING_PROFILES and "similarity_constraints" in config:
        numerical_changes["similarity_constraints"] = {
            "profile": None,
            "effective": config["similarity_constraints"],
        }
    # A profile match requires an actual complete qualified acceptance protocol.
    if config.get("refinement", {}).get("top_n", 0) <= 0:
        numerical_changes["refinement.top_n"] = {
            "profile": "positive",
            "effective": config.get("refinement", {}).get("top_n"),
        }
    budget = config.get("refinement", {}).get("max_total_conformers")
    if budget is None or budget < 24:
        numerical_changes["refinement.max_total_conformers"] = {"profile": "at least 24", "effective": budget}
    # Optional declared landmarks are additional vetoes, not replacement floors.
    if name in ENGINEERING_PROFILES:
        numerical_changes.pop("similarity_constraints.landmarks", None)
    current = numerical_runtime(name) if name in ENGINEERING_PROFILES else numerical_runtime()
    recorded = artifact["provenance"]
    mismatches = {}
    provenance_key = "engineering_policy" if name in ENGINEERING_PROFILES else "calibration"
    for field in ("python", "packages"):
        if current[field] != recorded["runtime"][field]:
            mismatches[field] = {provenance_key: recorded["runtime"][field], "current": current[field]}
    feature = "pharmacophore_feature_definition"
    if current[feature]["sha256"] != recorded["runtime"][feature]["sha256"]:
        mismatches[feature] = {provenance_key: recorded["runtime"][feature], "current": current[feature]}
    for path, digest in current["source_sha256"].items():
        if digest != recorded["source_sha256"].get(path):
            mismatches[path] = {provenance_key: recorded["source_sha256"].get(path), "current": digest}
    engineering = name in ENGINEERING_PROFILES
    return {
        "name": name,
        "status": "custom_unvalidated"
        if numerical_changes or mismatches
        else ("engineering_policy_matched" if engineering else "development_protocol_matched"),
        "artifact_sha256": PROFILE_SHA256[name],
        **(
            {"engineering_spec": artifact["engineering_spec"], "engineering_provenance": recorded}
            if engineering
            else {
                "fit_sha256": artifact["fit_sha256"],
                "review": artifact["review"],
                "calibration_spec": artifact["calibration_spec"],
                "calibration_provenance": recorded,
            }
        ),
        "scope": artifact["scope"],
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
                *(("similarity_constraints",) if engineering else ()),
            )
            if key in config
        },
        "source_match_scope": list(current["source_sha256"]),
        "limitations": artifact.get("limitations", [])
        + [
            (
                "Profile match describes declared engineering constraints, not fitted thresholds, activity or experimental synthesis"
                if engineering
                else "Profile match describes the development scoring/acceptance protocol, not pipeline retention or biological activity"
            ),
            "Real-design diagnostic availability, anchors/properties, finite budgets and route evidence impose additional vetoes",
            "Observed three-seed minima do not establish future-seed robustness or convergence",
        ],
    }


def profile_label(config):
    selective = (
        " [selector=joint-feasible-hit-pose-v2; explicit engineering constraints, not a calibrated similarity or activity claim]"
        if "similarity_constraints" in config
        else ""
    )
    if not config.get("recommendation_profile"):
        return (
            " [evidence=" + config.get("synthesis_evidence", {}).get("profile", "legacy-conditional-v1") + "]"
        ) + selective
    provenance = profile_provenance(config)
    if provenance["name"] in ENGINEERING_PROFILES:
        return (
            f" [profile={provenance['name']}; {provenance['status']}; "
            f"hit_tanimoto>={config['minimum_scores']['hit_tanimoto']}; "
            f"hit_shape>={config['minimum_scores']['hit_shape']}; "
            f"hit_pharmacophore>={config['minimum_scores']['hit_pharmacophore']}; "
            f"polar_recall>={config['similarity_constraints']['polar_recall']}; "
            "explicit engineering preferences, not fitted/calibrated; computational-route-support-v1 conditional support, no experimental synthesis or activity claim]"
            + selective
        )
    return (
        f" [profile={provenance['name']}; {provenance['status']}; computational-route-support-v1 conditional selection support; permissive development thresholds, not strong-similarity or activity guarantees]"
        + selective
    )
