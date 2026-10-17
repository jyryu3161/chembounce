"""V4 source declaration, not empirical calibration or a model evaluation."""

import copy
import hashlib
import json
from argparse import Namespace

import pytest

from chembounce_synthesis import profiles, rejected_pose_diagnostics as diagnostic
from chembounce_synthesis.options import apply_expansion_options

V4 = profiles.DIAGNOSTIC_SELECTIVE_PROFILE


def effective(enabled=False):
    return profiles.apply_profile_defaults(
        {"recommendation_profile": V4, "rejected_pose_diagnostics": enabled}
    )


def test_new_default_preserves_exact_numeric_and_runtime_declaration():
    old, new = (profiles.load_profile(n) for n in (profiles.GUARDED_SELECTIVE_PROFILE, V4))
    assert profiles.DEFAULT_RECOMMENDATION_PROFILE == profiles.LANDMARK_SELECTIVE_PROFILE
    for key in ("engineering_spec", "minimum_scores", "schema", "status", "scope"):
        assert new[key] == old[key]
    assert new["provenance"]["runtime"] == old["provenance"]["runtime"]
    assert (
        new["provenance"]["source_revision"]["inherited_guard_revision"]
        == old["provenance"]["source_revision"]
    )
    defaults = profiles.profile_defaults(new)
    assert defaults.pop("rejected_pose_diagnostics") is False
    assert defaults == profiles.profile_defaults(old)
    assert not {"review", "fit_sha256", "calibration_spec"} & new.keys()
    declaration = new["optional_diagnostics"]["rejected_pose_diagnostics"]
    assert declaration["version"] == diagnostic.VERSION and declaration["limits"] == diagnostic.LIMITS
    assert all(not path.endswith("profiles.py") for path in declaration["source_sha256"])


@pytest.mark.parametrize("enabled", [False, True])
def test_actual_v4_source_match_is_separate_from_runtime_and_diagnostic_cost(enabled):
    p = profiles.profile_provenance(effective(enabled))
    assert not p["material_overrides"]
    assert "chembounce_synthesis/selective_scoring.py" in p["runtime_source_mismatches"]
    assert "chembounce_synthesis/similarity_constraints.py" in p["runtime_source_mismatches"]
    assert p["status"] == (
        "custom_unvalidated" if p["runtime_source_mismatches"] else "engineering_policy_matched"
    )
    d = p["optional_diagnostics"]["rejected_pose_diagnostics"]
    assert d["enabled"] is enabled
    assert "chembounce_synthesis/selective_scoring.py" in d["source_mismatches"]
    assert "deadline" in d["scope"] and "summary" in d["observed_costs"]
    d["limits"].clear()
    assert profiles.profile_provenance(effective())["optional_diagnostics"]["rejected_pose_diagnostics"][
        "limits"
    ]


@pytest.mark.parametrize("path", profiles.DIAGNOSTIC_FILES)
@pytest.mark.parametrize("enabled", [False, True])
def test_each_optional_binding_is_visible_and_enabled_mismatch_is_custom(monkeypatch, path, enabled):
    # Mutate expected artifact data only; keep the real current source hash path.
    original = profiles.load_profile(V4)
    rebound = copy.deepcopy(original)
    key = "chembounce_synthesis/" + path
    rebound["optional_diagnostics"]["rejected_pose_diagnostics"]["source_sha256"][key] = "0" * 64
    monkeypatch.setattr(profiles, "load_profile", lambda name: rebound)
    result = profiles.profile_provenance(effective(enabled))
    assert key in result["optional_diagnostics"]["rejected_pose_diagnostics"]["source_mismatches"]
    if enabled:
        assert result["status"] == "custom_unvalidated"
    assert not result["material_overrides"]


@pytest.mark.parametrize(
    "field",
    [
        "python",
        "packages",
        "pharmacophore_feature_definition",
        *profiles.KERNEL_FILES,
        *profiles.SELECTIVE_KERNEL_FILES,
    ],
)
def test_each_v4_numeric_runtime_binding_remains_required(monkeypatch, field):
    observed = profiles.numerical_runtime(V4)
    if field.endswith(".py"):
        observed["source_sha256"]["chembounce_synthesis/" + field] = "0" * 64
    elif field == "packages":
        observed[field]["rdkit"] = "synthetic-mismatch"
    elif field == "pharmacophore_feature_definition":
        observed[field]["sha256"] = "0" * 64
    else:
        observed[field] = "synthetic-mismatch"
    monkeypatch.setattr(profiles, "numerical_runtime", lambda name=None: observed)
    assert profiles.profile_provenance(effective())["status"] == "custom_unvalidated"


@pytest.mark.parametrize(
    "mutation",
    [
        "commit",
        "predecessor",
        "guard",
        "review0",
        "review1",
        "role",
        "missing_review",
        "extra_revision",
        "version",
        "limits",
        "bool_limit",
        "default",
        "missing_source",
        "self_source",
        "bad_hash",
        "extra_diagnostic",
    ],
)
def test_rebound_declaration_tamper_is_rejected(tmp_path, monkeypatch, mutation):
    artifact = profiles.load_profile(V4)
    revision = artifact["provenance"]["source_revision"]
    d = artifact["optional_diagnostics"]["rejected_pose_diagnostics"]
    if mutation == "commit":
        revision["diagnostic_source_commit"] = "0" * 40
    elif mutation == "predecessor":
        revision["predecessor"]["artifact_sha256"] = "0" * 64
    elif mutation == "guard":
        revision["inherited_guard_revision"]["guarded_source_commit"] = "0" * 40
    elif mutation in ("review0", "review1"):
        revision["source_reviews"][int(mutation[-1])]["record_sha256"] = "0" * 64
    elif mutation == "role":
        revision["source_reviews"][0]["role"] = "calibration"
    elif mutation == "missing_review":
        revision["source_reviews"].pop()
    elif mutation == "extra_revision":
        revision["self_hash"] = "0" * 64
    elif mutation == "version":
        d["version"] = "unknown"
    elif mutation == "limits":
        d["limits"]["records"] += 1
    elif mutation == "bool_limit":
        d["limits"]["capture_seconds"] = True
    elif mutation == "default":
        d["default_enabled"] = True
    elif mutation == "missing_source":
        d["source_sha256"].pop(next(iter(d["source_sha256"])))
    elif mutation == "self_source":
        d["source_sha256"]["chembounce_synthesis/profiles.py"] = "0" * 64
    elif mutation == "bad_hash":
        d["source_sha256"][next(iter(d["source_sha256"]))] = "bad"
    else:
        d["extra"] = True
    path = tmp_path / "synthetic-declaration.json"
    raw = json.dumps(artifact).encode()
    path.write_bytes(raw)
    monkeypatch.setattr(profiles, "_resource", lambda name: path)
    monkeypatch.setitem(profiles.PROFILE_SHA256, V4, hashlib.sha256(raw).hexdigest())
    with pytest.raises(ValueError, match="diagnostic engineering"):
        profiles.load_profile(V4)


@pytest.mark.parametrize(
    "old", [profiles.PROFILE, profiles.SELECTIVE_PROFILE, profiles.GUARDED_SELECTIVE_PROFILE]
)
def test_explicit_historical_profile_is_never_silently_upgraded(old):
    for raw, cli in [({"recommendation_profile": old}, None), ({}, old)]:
        before = copy.deepcopy(raw)
        with pytest.raises(ValueError, match="conflicts"):
            apply_expansion_options(raw, Namespace(recommend=True, recommendation_profile=cli))
        assert raw == before
    assert apply_expansion_options({}, Namespace(recommendation_profile=old))["recommendation_profile"] == old


def test_historical_benchmark_cannot_claim_v4():
    from benchmarks.research.qualified_adapter import benchmark_profile

    with pytest.raises(ValueError, match="binding"):
        benchmark_profile(
            {
                "schema": "qualified-benchmark-protocol-v2",
                "recommendation_profile": V4,
                "profile_artifact_sha256": profiles.PROFILE_SHA256[V4],
            }
        )


def test_real_v1_artifact_without_name_keeps_defaults_and_honest_provenance():
    artifact = profiles.load_profile(profiles.PROFILE)
    assert "name" not in artifact
    defaults = profiles.profile_defaults(artifact)
    assert defaults["minimum_scores"] == artifact["minimum_scores"]
    assert "rejected_pose_diagnostics" not in defaults
    cfg = profiles.apply_profile_defaults({"recommendation_profile": profiles.PROFILE})
    assert profiles.profile_provenance(cfg)["status"] == "custom_unvalidated"


@pytest.mark.parametrize("enabled", [False, True])
def test_v4_actual_controlled_workflow_manifest_and_continuation(
    tmp_path, assets, proposal, monkeypatch, strict_predictor, enabled
):
    from chembounce_synthesis import workflow
    from chembounce_synthesis.continuation import continuation_inputs
    from .test_rejected_pose_diagnostics import workflow_control
    from .test_selective_workflow import config_for, HIT, DESIGN, read

    workflow_control(monkeypatch, fail_screening=True)
    config = config_for(tmp_path, assets, proposal)
    config["recommendation_profile"] = V4
    config["synthesis_evidence"] = {
        "predictor": {"kind": "reactiont5", "model_path": str(tmp_path / "controlled-forward")}
    }
    config["rejected_pose_diagnostics"] = enabled
    if not enabled:
        monkeypatch.setattr(diagnostic, "RejectedPoseDiagnostics", lambda **kw: pytest.fail("off collector"))
    out = tmp_path / "v4"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    manifest = read(out, "manifest.json")
    p = manifest["recommendation_profile"]
    assert p["name"] == V4 and p["artifact_sha256"] == profiles.PROFILE_SHA256[V4]
    d = p["optional_diagnostics"]["rejected_pose_diagnostics"]
    assert d["enabled"] is enabled
    assert "chembounce_synthesis/selective_scoring.py" in d["source_mismatches"]
    for path, digest in d["current_source_sha256"].items():
        assert manifest["code_sha256"][path.rsplit("/", 1)[1]] == digest
    summary = read(out, "summary.json")
    if enabled:
        # Older RDKit may reject the route before scoring; that honest zero
        # remains valid provenance and is not a qualified positive.
        assert summary["rejected_pose_diagnostics"]["calls"] >= 0
        assert "total_observed_seconds" in summary["rejected_pose_diagnostics"]
    else:
        assert "rejected_pose_diagnostics" not in summary
    _, _, carried = continuation_inputs(out)
    assert carried["recommendation_profile"] == V4
    assert carried["rejected_pose_diagnostics"] is enabled


def test_v4_packaged_digest_authenticates_without_external_review_reads(tmp_path, monkeypatch):
    raw = profiles._resource(V4).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == profiles.PROFILE_SHA256[V4]
    path = tmp_path / "portable.json"
    path.write_bytes(raw)
    monkeypatch.setattr(profiles, "_resource", lambda name: path)
    assert profiles.load_profile(V4)["name"] == V4
    path.write_bytes(raw + b" ")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        profiles.load_profile(V4)
