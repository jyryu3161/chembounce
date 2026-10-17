"""V6 declared method/source binding; no efficacy or calibration inference."""

import argparse
import copy
import hashlib
import json

import pytest

from chembounce_synthesis import profiles
from chembounce_synthesis.options import add_expansion_options, apply_expansion_options

V6 = profiles.ALLOCATION_SELECTIVE_PROFILE
WORKFLOW = "chembounce_synthesis/workflow.py"


def config(**overrides):
    return profiles.apply_profile_defaults({"recommendation_profile": V6, **overrides})


@pytest.fixture
def matching_runtime(monkeypatch):
    """Explicit metadata seam; never authenticate the installed production runtime."""
    artifact = profiles.load_profile(V6)
    runtime = copy.deepcopy(artifact["provenance"]["runtime"])
    runtime["source_sha256"] = copy.deepcopy(artifact["provenance"]["source_sha256"])
    monkeypatch.setattr(profiles, "numerical_runtime", lambda name: copy.deepcopy(runtime))
    return runtime


def test_packaged_v6_adds_only_declared_allocations_and_workflow_scope():
    old = profiles.load_profile(profiles.LANDMARK_SELECTIVE_PROFILE)
    new = profiles.load_profile(V6)
    assert new["engineering_spec"] == old["engineering_spec"] | profiles.ALLOCATION_SPEC
    assert new["minimum_scores"] == old["minimum_scores"]
    assert new["provenance"]["runtime"] == old["provenance"]["runtime"]
    assert new["optional_landmark_proposals"] == old["optional_landmark_proposals"]
    assert (
        new["provenance"]["source_revision"]["inherited_landmark_revision"]
        == old["provenance"]["source_revision"]
    )
    observed = profiles.numerical_runtime(V6)
    assert observed["source_sha256"] == new["provenance"]["source_sha256"]
    assert set(observed["source_sha256"]) == set(
        profiles.numerical_runtime(profiles.LANDMARK_SELECTIVE_PROFILE)["source_sha256"]
    ) | {WORKFLOW}
    assert len(observed["source_sha256"]) == 7
    assert hashlib.sha256(profiles._resource(V6).read_bytes()).hexdigest() == profiles.PROFILE_SHA256[V6]
    result = profiles.profile_provenance(config())
    assert not result["material_overrides"]
    assert not result["optional_diagnostics"]["rejected_pose_diagnostics"]["source_mismatches"]
    assert not result["optional_landmark_proposals"]["enabled"]
    assert result["status"] == (
        "custom_unvalidated" if result["runtime_source_mismatches"] else "engineering_policy_matched"
    )
    # Historical v5 is byte-identical and still the recommendation entrypoint.
    assert (
        hashlib.sha256(profiles._resource(profiles.LANDMARK_SELECTIVE_PROFILE).read_bytes()).hexdigest()
        == "1273082948d008dc7dd02daf18ace1bc6aeeea78515ebf074f321a3c28a450ad"
    )
    assert (
        profiles.request_recommendations({})["recommendation_profile"] == profiles.LANDMARK_SELECTIVE_PROFILE
    )


def test_v6_is_explicit_and_does_not_retarget_recommend():
    parser = argparse.ArgumentParser()
    add_expansion_options(parser)
    args = parser.parse_args(["--recommendation-profile", V6])
    selected = profiles.apply_profile_defaults(apply_expansion_options({}, args))
    assert all(selected[key] == value for key, value in profiles.ALLOCATION_SPEC.items())
    with pytest.raises(ValueError, match="conflicts"):
        apply_expansion_options({}, parser.parse_args(["--recommend", "--recommendation-profile", V6]))
    with pytest.raises(ValueError, match="conflicts"):
        profiles.request_recommendations({"recommendation_profile": V6})


@pytest.mark.parametrize("key", list(profiles.ALLOCATION_SPEC))
@pytest.mark.parametrize("value", ["inherit-v1", "unknown-policy"])
def test_v6_deviation_is_material_even_when_other_policy_matches(matching_runtime, key, value):
    result = profiles.profile_provenance(config(**{key: value}))
    assert result["status"] == "custom_unvalidated"
    assert result["material_overrides"] == {
        key: {"profile": profiles.ALLOCATION_SPEC[key], "effective": value}
    }
    assert key not in result["resource_overrides"]


def test_v6_primary_workflow_pin_cannot_be_hidden_by_disabled_diagnostics(matching_runtime):
    matching_runtime["source_sha256"][WORKFLOW] = "0" * 64
    result = profiles.profile_provenance(config(rejected_pose_diagnostics=False))
    assert result["status"] == "custom_unvalidated"
    assert set(result["runtime_source_mismatches"]) == {WORKFLOW}
    assert not result["material_overrides"]


def test_v6_optional_diagnostics_and_landmarks_keep_their_semantics(matching_runtime):
    result = profiles.profile_provenance(config(rejected_pose_diagnostics=True))
    assert result["status"] == "engineering_policy_matched"
    assert result["optional_diagnostics"]["rejected_pose_diagnostics"]["enabled"]
    result = profiles.profile_provenance(config(similarity_constraints={"landmark_pose_proposals": True}))
    assert result["status"] == "custom_unvalidated"
    assert set(result["material_overrides"]) == {"similarity_constraints.landmark_pose_proposals"}


@pytest.mark.parametrize(
    "attack",
    ["name", "schema", "screening", "evidence", "extra_spec", "missing_workflow", "revision", "landmark"],
)
def test_rehashed_v6_artifact_cannot_change_fixed_contract(tmp_path, monkeypatch, attack):
    artifact = profiles.load_profile(V6)
    if attack in ("name", "schema"):
        artifact[attack] = "unregistered"
    elif attack in ("screening", "evidence"):
        artifact["engineering_spec"][attack + "_allocation_policy"] = "inherit-v1"
    elif attack == "extra_spec":
        artifact["engineering_spec"]["trust_me"] = True
    elif attack == "missing_workflow":
        artifact["provenance"]["source_sha256"].pop(WORKFLOW)
    elif attack == "revision":
        artifact["provenance"]["source_revision"]["allocation_source_commit"] = "0" * 40
    else:
        artifact["optional_landmark_proposals"]["default_enabled"] = True
    path = tmp_path / "tampered.json"
    path.write_text(json.dumps(artifact))
    monkeypatch.setattr(profiles, "_resource", lambda name: path)
    monkeypatch.setitem(profiles.PROFILE_SHA256, V6, hashlib.sha256(path.read_bytes()).hexdigest())
    with pytest.raises(ValueError, match="Packaged"):
        profiles.load_profile(V6)
