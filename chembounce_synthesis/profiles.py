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
DIAGNOSTIC_SELECTIVE_PROFILE = "qualified-selective-engineering-v4"
LANDMARK_SELECTIVE_PROFILE = "qualified-selective-engineering-v5"
DEFAULT_RECOMMENDATION_PROFILE = LANDMARK_SELECTIVE_PROFILE
ENGINEERING_PROFILES = (
    SELECTIVE_PROFILE,
    GUARDED_SELECTIVE_PROFILE,
    DIAGNOSTIC_SELECTIVE_PROFILE,
    LANDMARK_SELECTIVE_PROFILE,
)
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
# No declaration/self hashes: the artifact registry and external review bind these.
DIAGNOSTIC_SOURCE_REVISION = {
    "diagnostic_source_commit": "612596ae8921e663dda5a239c662f50d782704f7",
    "predecessor": {
        "name": GUARDED_SELECTIVE_PROFILE,
        "artifact_sha256": "fd5fb3b65b12307ed1a9b5371e74427eb709d8cf3012bbcba5095022ddb34537",
    },
    "inherited_guard_revision": copy.deepcopy(GUARDED_SOURCE_REVISION),
    "source_reviews": [
        {
            "role": "independent_source_review",
            "record_label": "review-rejected-pose-source-612596a-v1/report.json",
            "record_sha256": "8f535615f9339580fb416bf4d81cbfbb4661402deb7a84bb323bdb7a67b58c92",
        },
        {
            "role": "root_source_gate",
            "record_label": "stage18-rejected-pose-root-source-gate-612596a-v1.json",
            "record_sha256": "ec71e9c72527d12075209ca01e830140970d307e7451f1e478e007c1fdad9fbf",
        },
    ],
}
DIAGNOSTIC_FILES = ("rejected_pose_diagnostics.py", "refinement.py", "workflow.py", "selective_scoring.py")

# Portable prior-source review, not calibration or a self-referential package pin.
LANDMARK_SOURCE_REVISION = {
    "landmark_source_commit": "427fdbd22ba40c3fde2cb4a4a43a3a7b64c85ad8",
    "predecessor": {
        "name": DIAGNOSTIC_SELECTIVE_PROFILE,
        "artifact_sha256": "437b5c11f9b5627122e2ad70ecbc2d4f9f8b2e5294c6fbe8be81058ffb8e48b5",
    },
    "inherited_diagnostic_revision": copy.deepcopy(DIAGNOSTIC_SOURCE_REVISION),
    "source_reviews": [
        {
            "role": "independent_source_review",
            "record_label": "review-landmark-correction-427fdbd-v1/report.json",
            "record_sha256": "f57aa678ba4fdc3977eec3a94961bbeda9d61c6b8fb3fd7959c5c4a74efa3402",
        },
        {
            "role": "root_source_gate",
            "record_label": "stage20-landmark-prototype-root-source-gate-427fdbd-v1.json",
            "record_sha256": "9745470ecfe3e6d7a9a3ba4c7bc57de8e29f7305e9a158ded79ceb0add9b3a33",
        },
    ],
}
# Declarative values: loading the disabled profile must not import or require ScoreMol.
LANDMARK_DECLARATION = {
    "default_enabled": False,
    "version": "landmark-anchored-rigid-v1",
    "limits": {
        "alternatives_per_pair": 2,
        "alternatives_per_call": 8,
        "scoremol_per_call": 16,
        "mapping_edge_visits_per_pair": 128,
        "complete_maps_per_pair": 4,
        "helper_edge_visits_per_pair": 4096,
        "cooperative_seconds_per_call": 1.0,
    },
    "parity_tolerance": 1e-5,
}
LANDMARK_FILE = "chembounce_synthesis/landmark_proposals.py"

ENGINEERING_SCHEMA = "qualified-selective-engineering-profile-v2"
SELECTIVE_KERNEL_FILES = ("selective_scoring.py", "similarity_constraints.py")
SELECTIVE_CONSTRAINTS = {
    "version": "selective-similarity-constraints-v2",
    "polar_recall": 0.75,
    "landmarks": [],
}
# Reviewed development release; the untouched fit candidate remains separately bound.
PROFILE_SHA256 = {
    LANDMARK_SELECTIVE_PROFILE: "1273082948d008dc7dd02daf18ace1bc6aeeea78515ebf074f321a3c28a450ad",
    DIAGNOSTIC_SELECTIVE_PROFILE: "437b5c11f9b5627122e2ad70ecbc2d4f9f8b2e5294c6fbe8be81058ffb8e48b5",
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
    if artifact.get("name") in (DIAGNOSTIC_SELECTIVE_PROFILE, LANDMARK_SELECTIVE_PROFILE):
        defaults["rejected_pose_diagnostics"] = False
    if artifact.get("name") == LANDMARK_SELECTIVE_PROFILE:
        defaults["similarity_constraints"]["landmark_pose_proposals"] = False
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
        != (
            {"chembounce_synthesis/" + path for path in KERNEL_FILES + SELECTIVE_KERNEL_FILES}
            | ({LANDMARK_FILE} if name == LANDMARK_SELECTIVE_PROFILE else set())
        )
    ):
        raise ValueError("Packaged engineering profile lacks its declared fixed policy/source scope")
    if (
        name == GUARDED_SELECTIVE_PROFILE
        and artifact["provenance"].get("source_revision") != GUARDED_SOURCE_REVISION
    ):
        raise ValueError(
            "Packaged guarded engineering profile lacks its reviewed source revision declaration"
        )

    if name in (DIAGNOSTIC_SELECTIVE_PROFILE, LANDMARK_SELECTIVE_PROFILE):
        from .rejected_pose_diagnostics import LIMITS, VERSION

        diagnostic = artifact.get("optional_diagnostics", {}).get("rejected_pose_diagnostics", {})
        if (
            artifact["provenance"].get("source_revision")
            != (
                LANDMARK_SOURCE_REVISION if name == LANDMARK_SELECTIVE_PROFILE else DIAGNOSTIC_SOURCE_REVISION
            )
            or set(artifact.get("optional_diagnostics", {})) != {"rejected_pose_diagnostics"}
            or set(diagnostic) != {"default_enabled", "version", "limits", "source_sha256"}
            or diagnostic["default_enabled"] is not False
            or diagnostic["version"] != VERSION
            or diagnostic["limits"] != LIMITS
            or any(type(diagnostic["limits"][key]) is not type(value) for key, value in LIMITS.items())
            or set(diagnostic["source_sha256"]) != {"chembounce_synthesis/" + n for n in DIAGNOSTIC_FILES}
            or any(
                not isinstance(value, str)
                or len(value) != 64
                or any(c not in "0123456789abcdef" for c in value)
                for value in diagnostic["source_sha256"].values()
            )
        ):
            raise ValueError(
                "Packaged diagnostic engineering profile lacks its reviewed optional source declaration"
            )

    if name == LANDMARK_SELECTIVE_PROFILE:
        declaration = artifact.get("optional_landmark_proposals", {})
        expected = {
            **LANDMARK_DECLARATION,
            "source_sha256": {LANDMARK_FILE: artifact["provenance"]["source_sha256"][LANDMARK_FILE]},
        }
        digest = expected["source_sha256"][LANDMARK_FILE]
        if (
            declaration != expected
            or type(declaration.get("default_enabled")) is not bool
            or any(
                type(declaration.get("limits", {}).get(key)) is not type(value)
                for key, value in LANDMARK_DECLARATION["limits"].items()
            )
            or type(declaration.get("parity_tolerance")) is not float
            or not isinstance(digest, str)
            or len(digest) != 64
            or any(c not in "0123456789abcdef" for c in digest)
        ):
            raise ValueError("Packaged landmark engineering profile lacks its declared opt-in source policy")


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
                KERNEL_FILES
                + SELECTIVE_KERNEL_FILES
                + (("landmark_proposals.py",) if name == LANDMARK_SELECTIVE_PROFILE else ())
                if name in ENGINEERING_PROFILES
                else KERNEL_FILES
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
    # Diagnostic output consumes resources but never changes the numeric specification.
    defaults.pop("rejected_pose_diagnostics", None)
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
    diagnostics = None
    diagnostic_mismatches = {}
    if name in (DIAGNOSTIC_SELECTIVE_PROFILE, LANDMARK_SELECTIVE_PROFILE):
        declaration = copy.deepcopy(artifact["optional_diagnostics"]["rejected_pose_diagnostics"])
        observed = {
            "chembounce_synthesis/" + path: _sha(files("chembounce_synthesis").joinpath(path).read_bytes())
            for path in DIAGNOSTIC_FILES
        }
        diagnostic_mismatches = {
            path: {"declared": expected, "current": observed[path]}
            for path, expected in declaration["source_sha256"].items()
            if observed[path] != expected
        }
        diagnostics = {
            **declaration,
            "enabled": config.get("rejected_pose_diagnostics", False),
            "current_source_sha256": observed,
            "source_mismatches": diagnostic_mismatches,
            "scope": "Optional rejected-pose evidence, not selection, calibration or efficacy. Costs consume the existing deadline and can affect deadline-limited coverage; no extra allowance or identical wall-time outcome claim.",
            "observed_costs": "Profile declaration/source hashing is setup work outside collector statistics and consumes the enclosing workflow deadline when supplied. Run summary.rejected_pose_diagnostics records collector setup/capture/publication time, skips and overruns; final summary-write completion requires its post-return observation.",
        }
    return {
        "name": name,
        "status": "custom_unvalidated"
        if numerical_changes
        or mismatches
        or (diagnostics and diagnostics["enabled"] and diagnostic_mismatches)
        else ("engineering_policy_matched" if engineering else "development_protocol_matched"),
        "artifact_sha256": PROFILE_SHA256[name],
        **(
            {
                "optional_landmark_proposals": {
                    **copy.deepcopy(artifact["optional_landmark_proposals"]),
                    "enabled": config.get("similarity_constraints", {}).get("landmark_pose_proposals", False),
                    "scope": "Explicit landmark proposals are a material search override and remain custom_unvalidated. Default false requires no ScoreMol capability. Source review is not calibration or efficacy evidence.",
                }
            }
            if name == LANDMARK_SELECTIVE_PROFILE
            else {}
        ),
        **({"optional_diagnostics": {"rejected_pose_diagnostics": diagnostics}} if diagnostics else {}),
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
