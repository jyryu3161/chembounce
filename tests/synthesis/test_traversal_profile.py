"""V7 is a new source binding of v5 policy, never a historical-profile alias."""

import argparse
import copy
import hashlib
import json

import pytest

from chembounce_synthesis import profiles, workflow
from chembounce_synthesis.options import add_expansion_options, apply_expansion_options

V7 = profiles.TRAVERSAL_SELECTIVE_PROFILE
COMMON = "chembounce_synthesis/common.py"
HISTORICAL = {
    profiles.PROFILE: "9d3c60a895912b28c13a29aac0a04f385ddb8a0effc9fcd81d64a50a3d48bdd2",
    profiles.SELECTIVE_PROFILE: "a8b47a409a21518dd482acbefdce38774414ae0b65f28d6ea7cb52ab63e1ef52",
    profiles.GUARDED_SELECTIVE_PROFILE: "fd5fb3b65b12307ed1a9b5371e74427eb709d8cf3012bbcba5095022ddb34537",
    profiles.DIAGNOSTIC_SELECTIVE_PROFILE: "437b5c11f9b5627122e2ad70ecbc2d4f9f8b2e5294c6fbe8be81058ffb8e48b5",
    profiles.LANDMARK_SELECTIVE_PROFILE: "1273082948d008dc7dd02daf18ace1bc6aeeea78515ebf074f321a3c28a450ad",
    profiles.ALLOCATION_SELECTIVE_PROFILE: "c1d257b25286c14dc3f8208036018bf706d289c0a227f903727ada99fdf0c37a",
}


def config(**overrides):
    return profiles.apply_profile_defaults({"recommendation_profile": V7, **overrides})


def test_v7_inherits_v5_policy_and_only_changes_primary_common_binding():
    old = profiles.load_profile(profiles.LANDMARK_SELECTIVE_PROFILE)
    new = profiles.load_profile(V7)
    for key in (
        "engineering_spec",
        "minimum_scores",
        "schema",
        "status",
        "scope",
        "optional_diagnostics",
        "optional_landmark_proposals",
    ):
        assert new[key] == old[key]
    assert new["provenance"]["runtime"] == old["provenance"]["runtime"]
    assert new["provenance"]["source_revision"] == profiles.TRAVERSAL_SOURCE_REVISION
    assert new["provenance"]["source_revision"]["predecessor"] == {
        "name": profiles.LANDMARK_SELECTIVE_PROFILE,
        "artifact_sha256": HISTORICAL[profiles.LANDMARK_SELECTIVE_PROFILE],
    }
    assert (
        new["provenance"]["source_revision"]["inherited_landmark_revision"]
        == old["provenance"]["source_revision"]
    )
    before, after = old["provenance"]["source_sha256"], new["provenance"]["source_sha256"]
    assert set(before) == set(after) and len(after) == 6
    assert {key for key in before if before[key] != after[key]} == {COMMON}
    assert after[COMMON] == "83226cd0a354003a350491ad018546d78a1cfbb86338d855a6216d822552164b"
    assert profiles.profile_defaults(new) == profiles.profile_defaults(old)
    assert not set(profiles.ALLOCATION_SPEC) & new["engineering_spec"].keys()
    assert not {"review", "fit_sha256", "calibration_spec"} & new.keys()


def test_effective_workflow_defaults_are_v5_with_only_profile_identity_changed():
    base = {
        "schema_version": 1,
        "catalog": "unused.tsv",
        "rules": "unused.json",
        "provider": {"kind": "saved", "path": "unused.json"},
        "synthesis_evidence": {"predictor": {"kind": "reactiont5", "model_path": "unused"}},
    }
    old = workflow._config(base | {"recommendation_profile": profiles.LANDMARK_SELECTIVE_PROFILE})
    new = workflow._config(profiles.request_recommendations(base))
    assert old.pop("recommendation_profile") == profiles.LANDMARK_SELECTIVE_PROFILE
    assert new.pop("recommendation_profile") == V7
    assert new == old
    for key in profiles.ALLOCATION_SPEC:
        assert new.get(key, "inherit-v1") == "inherit-v1"
    assert new["conformers"]["scoring_objective"] == "original-hit-v1"
    assert new["ranking_policy"] == new["route_search_policy"] == "original-hit-v1"
    assert new["synthesis_evidence"]["policy"]["require_forward_agreement"] is True


def test_current_primary_sources_match_and_optional_diagnostic_mismatch_is_honest():
    current = profiles.numerical_runtime(V7)
    artifact = profiles.load_profile(V7)
    assert current["source_sha256"] == artifact["provenance"]["source_sha256"]
    result = profiles.profile_provenance(config())
    assert not result["material_overrides"]
    assert not any(key.startswith("chembounce_synthesis/") for key in result["runtime_source_mismatches"])
    # CI may use a different numerical runtime; it must never be relabelled matched.
    assert result["status"] == (
        "custom_unvalidated" if result["runtime_source_mismatches"] else "engineering_policy_matched"
    )
    assert result["screening_allocation_policy"] == result["evidence_allocation_policy"] == "inherit-v1"
    diagnostic = result["optional_diagnostics"]["rejected_pose_diagnostics"]
    assert diagnostic["enabled"] is False
    assert set(diagnostic["source_mismatches"]) == {"chembounce_synthesis/workflow.py"}
    enabled = profiles.profile_provenance(config(rejected_pose_diagnostics=True))
    assert enabled["status"] == "custom_unvalidated"
    assert enabled["runtime_source_mismatches"] == result["runtime_source_mismatches"]
    assert (
        enabled["optional_diagnostics"]["rejected_pose_diagnostics"]["source_mismatches"]
        == diagnostic["source_mismatches"]
    )


@pytest.mark.parametrize("name,digest", HISTORICAL.items())
def test_old_artifact_bytes_and_names_are_preserved_and_source_mismatch_remains(name, digest):
    assert profiles.PROFILE_SHA256[name] == digest
    assert hashlib.sha256(profiles._resource(name).read_bytes()).hexdigest() == digest
    requested = profiles.apply_profile_defaults({"recommendation_profile": name})
    assert requested["recommendation_profile"] == name
    result = profiles.profile_provenance(requested)
    assert result["status"] == "custom_unvalidated"
    assert COMMON in result["runtime_source_mismatches"]
    assert not result["material_overrides"]


@pytest.mark.parametrize("name", HISTORICAL)
@pytest.mark.parametrize("origin", ["config", "cli"])
def test_default_v7_never_conceals_an_explicit_historical_profile(name, origin):
    parser = argparse.ArgumentParser()
    add_expansion_options(parser)
    raw = {"recommendation_profile": name if origin == "config" else V7}
    before = copy.deepcopy(raw)
    cli = V7 if origin == "config" else name
    with pytest.raises(ValueError, match="conflicts.*requires qualified-selective-engineering-v7"):
        apply_expansion_options(raw, parser.parse_args(["--recommend", "--recommendation-profile", cli]))
    assert raw == before
    with pytest.raises(ValueError, match="conflicts"):
        profiles.request_recommendations({"recommendation_profile": name})
    assert apply_expansion_options({}, parser.parse_args(["--recommendation-profile", name])) == {
        "recommendation_profile": name
    }


def test_v7_shortcut_explicit_choice_and_copy_are_consistent():
    parser = argparse.ArgumentParser()
    add_expansion_options(parser)
    raw = {"max_results": 1, "minimum_scores": {"hit_shape": 0.8}}
    before = copy.deepcopy(raw)
    via_cli = apply_expansion_options(raw, parser.parse_args(["--recommend", "--recommendation-profile", V7]))
    via_api = profiles.request_recommendations(raw)
    assert via_api == via_cli == raw | {"recommendation_profile": V7}
    assert profiles.DEFAULT_RECOMMENDATION_PROFILE == V7
    assert raw == before
    result = profiles.profile_provenance(profiles.apply_profile_defaults(via_api))
    assert set(result["material_overrides"]) == {"minimum_scores.hit_shape"}
    assert set(result["resource_overrides"]) == {"max_results"}


@pytest.mark.parametrize("key", tuple(profiles.ALLOCATION_SPEC))
def test_v6_allocations_remain_material_overrides_on_v7(key):
    result = profiles.profile_provenance(config(**{key: profiles.ALLOCATION_SPEC[key]}))
    assert result["status"] == "custom_unvalidated"
    assert result["material_overrides"] == {
        key: {"profile": "inherit-v1", "effective": profiles.ALLOCATION_SPEC[key]}
    }


@pytest.mark.parametrize(
    "path", [*profiles.KERNEL_FILES, *profiles.SELECTIVE_KERNEL_FILES, "landmark_proposals.py"]
)
def test_each_primary_source_drift_still_invalidates_v7(monkeypatch, path):
    current = profiles.numerical_runtime(V7)
    key = "chembounce_synthesis/" + path
    current["source_sha256"][key] = "0" * 64
    monkeypatch.setattr(profiles, "numerical_runtime", lambda name: current)
    result = profiles.profile_provenance(config())
    assert result["status"] == "custom_unvalidated"
    assert key in result["runtime_source_mismatches"]


@pytest.mark.parametrize("change", ["revision", "parent", "allocation", "scope", "self_hash"])
def test_rebound_v7_artifact_cannot_change_its_declaration(tmp_path, monkeypatch, change):
    artifact = profiles.load_profile(V7)
    if change == "revision":
        artifact["provenance"]["source_revision"]["source_reviews"].pop()
    elif change == "parent":
        artifact["provenance"]["source_revision"]["predecessor"]["name"] = (
            profiles.ALLOCATION_SELECTIVE_PROFILE
        )
    elif change == "allocation":
        artifact["engineering_spec"].update(profiles.ALLOCATION_SPEC)
    elif change == "scope":
        artifact["provenance"]["source_sha256"].pop(COMMON)
    else:
        artifact["provenance"]["source_revision"]["self_hash"] = "0" * 64
    raw = json.dumps(artifact).encode()
    path = tmp_path / "SYNTHETIC.json"
    path.write_bytes(raw)
    monkeypatch.setattr(profiles, "_resource", lambda name: path)
    monkeypatch.setitem(profiles.PROFILE_SHA256, V7, hashlib.sha256(raw).hexdigest())
    with pytest.raises(ValueError, match="engineering profile"):
        profiles.load_profile(V7)


def test_v7_resource_is_portable_and_remains_byte_pinned(tmp_path, monkeypatch):
    raw = profiles._resource(V7).read_bytes()
    path = tmp_path / "packaged.json"
    path.write_bytes(raw)
    monkeypatch.setattr(profiles, "_resource", lambda name: path)
    assert profiles.load_profile(V7)["name"] == V7
    path.write_bytes(raw + b" ")
    with pytest.raises(ValueError, match="SHA-256"):
        profiles.load_profile(V7)
