"""Test authority boundaries: synthetic binding is not a released profile match."""

import hashlib
import json

import pytest

from chembounce_synthesis import profiles, workflow
from benchmarks.research.qualified_adapter import common_config

PACKAGE_SHA256 = {
    profiles.PROFILE: "9d3c60a895912b28c13a29aac0a04f385ddb8a0effc9fcd81d64a50a3d48bdd2",
    profiles.SELECTIVE_PROFILE: "a8b47a409a21518dd482acbefdce38774414ae0b65f28d6ea7cb52ab63e1ef52",
}


def config(name):
    return {
        "schema_version": 1,
        "recommendation_profile": name,
        "catalog": "unused.tsv",
        "rules": "unused.json",
        "provider": {"kind": "synllama", "model_path": "unused", "prompt_path": "unused"},
        "synthesis_evidence": {"predictor": {"kind": "reactiont5", "model_path": "unused"}},
    }


@pytest.mark.parametrize("name", list(PACKAGE_SHA256))
def test_real_historical_package_remains_unchanged_and_unmatched(name, tmp_path):
    # No synthetic binding fixture here. A source mismatch cannot certify a fit.
    assert {name: profiles.PROFILE_SHA256[name] for name in PACKAGE_SHA256} == PACKAGE_SHA256
    assert hashlib.sha256(profiles._resource(name).read_bytes()).hexdigest() == PACKAGE_SHA256[name]
    artifact = profiles.load_profile(name)
    assert "synthetic_fixture" not in artifact["provenance"]
    effective = workflow._config(config(name))
    provenance = profiles.profile_provenance(effective)
    assert provenance["status"] == "custom_unvalidated"
    assert "chembounce_synthesis/common.py" in provenance["runtime_source_mismatches"]
    path = tmp_path / "real-profile.json"
    path.write_text(json.dumps(config(name)))
    request = {"assets": {"synthesis_config": str(path)}}
    if name == profiles.SELECTIVE_PROFILE:
        request["benchmark_protocol"] = {
            "schema": "qualified-benchmark-protocol-v2",
            "recommendation_profile": name,
            "profile_artifact_sha256": PACKAGE_SHA256[name],
        }
    with pytest.raises(ValueError, match="must match"):
        common_config(request)


@pytest.mark.parametrize("name", list(PACKAGE_SHA256))
def test_synthetic_binding_uses_real_validation_and_preserves_policy(
    name, synthetic_historical_profile_bindings
):
    fixture = synthetic_historical_profile_bindings
    artifact = profiles.load_profile(name)
    assert "SYNTHETIC" in artifact["provenance"]["synthetic_fixture"]["purpose"]
    assert profiles.PROFILE_SHA256[name] != PACKAGE_SHA256[name]
    # The original artifact is read by file path, independent of the fixture resource lookup.
    from importlib.resources import files

    original = json.loads(files("chembounce_synthesis").joinpath("data", name + ".json").read_bytes())
    assert profiles.profile_defaults(artifact) == profiles.profile_defaults(original)
    assert artifact["minimum_scores"] == original["minimum_scores"]
    assert profiles.profile_provenance(workflow._config(config(name)))["status"] == (
        "engineering_policy_matched" if name == profiles.SELECTIVE_PROFILE else "development_protocol_matched"
    )
    # Corruption is rejected by the actual packaged-profile byte check.
    fixture["paths"][name].write_text("{}")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        profiles.load_profile(name)


@pytest.mark.parametrize("name", list(PACKAGE_SHA256))
def test_rebound_synthetic_wrong_source_is_still_custom(
    name, synthetic_historical_profile_bindings, monkeypatch
):
    path = synthetic_historical_profile_bindings["paths"][name]
    artifact = json.loads(path.read_text())
    artifact["provenance"]["source_sha256"]["chembounce_synthesis/common.py"] = "0" * 64
    path.write_text(json.dumps(artifact))
    with monkeypatch.context() as scoped:
        scoped.setitem(profiles.PROFILE_SHA256, name, hashlib.sha256(path.read_bytes()).hexdigest())
        provenance = profiles.profile_provenance(workflow._config(config(name)))
        assert provenance["status"] == "custom_unvalidated"
        assert "chembounce_synthesis/common.py" in provenance["runtime_source_mismatches"]
