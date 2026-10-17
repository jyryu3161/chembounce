"""Engineering v3 registry/provenance contracts; synthetic mutations are not calibration."""

import copy
import hashlib
import json
from argparse import Namespace

import pytest

from chembounce_synthesis import profiles
from chembounce_synthesis.options import apply_expansion_options

V3 = profiles.GUARDED_SELECTIVE_PROFILE


def effective(name=V3):
    return profiles.apply_profile_defaults({"recommendation_profile": name})


def test_real_packaged_identity_policy_and_historical_authorities():
    expected = {
        profiles.PROFILE: "9d3c60a895912b28c13a29aac0a04f385ddb8a0effc9fcd81d64a50a3d48bdd2",
        profiles.SELECTIVE_PROFILE: "a8b47a409a21518dd482acbefdce38774414ae0b65f28d6ea7cb52ab63e1ef52",
        V3: "fd5fb3b65b12307ed1a9b5371e74427eb709d8cf3012bbcba5095022ddb34537",
    }
    assert profiles.SELECTIVE_PROFILE == "qualified-selective-engineering-v2"
    assert profiles.GUARDED_SELECTIVE_PROFILE == V3
    for name, digest in expected.items():
        assert profiles.PROFILE_SHA256[name] == digest
        assert hashlib.sha256(profiles._resource(name).read_bytes()).hexdigest() == digest
    old, new = (profiles.load_profile(n) for n in (profiles.SELECTIVE_PROFILE, V3))
    assert profiles.profile_defaults(old) == profiles.profile_defaults(new)
    assert old["engineering_spec"] == new["engineering_spec"]
    assert old["minimum_scores"] == new["minimum_scores"]
    assert new["provenance"]["source_revision"] == profiles.GUARDED_SOURCE_REVISION
    assert not {"review", "fit_sha256", "calibration_spec"} & new.keys()
    for name in (profiles.PROFILE, profiles.SELECTIVE_PROFILE):
        provenance = profiles.profile_provenance(effective(name))
        assert provenance["status"] == "custom_unvalidated"
        assert "chembounce_synthesis/common.py" in provenance["runtime_source_mismatches"]


def test_real_v3_runtime_is_honest_on_each_runtime_without_rebinding():
    provenance = profiles.profile_provenance(effective())
    assert not provenance["material_overrides"]
    assert set(provenance["source_match_scope"]) == {
        "chembounce_synthesis/" + n for n in profiles.KERNEL_FILES + profiles.SELECTIVE_KERNEL_FILES
    }
    assert "chembounce_synthesis/selective_scoring.py" in provenance["runtime_source_mismatches"]
    assert provenance["status"] == (
        "custom_unvalidated" if provenance["runtime_source_mismatches"] else "engineering_policy_matched"
    )
    assert "engineering_provenance" in provenance and "calibration_provenance" not in provenance
    assert V3 in profiles.profile_label(effective())


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
def test_individual_runtime_or_kernel_mismatch_is_custom(monkeypatch, field):
    current = profiles.numerical_runtime(V3)
    if field.endswith(".py"):
        current["source_sha256"]["chembounce_synthesis/" + field] = "0" * 64
    elif field == "packages":
        current[field]["rdkit"] = "synthetic-mismatch"
    elif field == "pharmacophore_feature_definition":
        current[field]["sha256"] = "0" * 64
    else:
        current[field] = "synthetic-mismatch"
    monkeypatch.setattr(profiles, "numerical_runtime", lambda name=None: current)
    assert profiles.profile_provenance(effective())["status"] == "custom_unvalidated"


@pytest.mark.parametrize("raw", [profiles.PROFILE, profiles.SELECTIVE_PROFILE, None, "unknown"])
@pytest.mark.parametrize("cli", [None, V3])
def test_original_conflict_cannot_be_concealed(raw, cli):
    config = {"recommendation_profile": raw, "nested": {"unchanged": []}}
    before = copy.deepcopy(config)
    with pytest.raises(ValueError, match="conflicts"):
        apply_expansion_options(config, Namespace(recommend=True, recommendation_profile=cli))
    assert config == before


@pytest.mark.parametrize("cli", [profiles.PROFILE, profiles.SELECTIVE_PROFILE, "unknown"])
def test_cli_conflict_cannot_be_concealed(cli):
    with pytest.raises(ValueError, match="conflicts"):
        profiles.request_recommendations({"recommendation_profile": V3}, cli)


def test_idempotent_copy_custom_and_resource_additional_selection():
    config = {"minimum_scores": {"hit_shape": 0.8}, "max_results": 1, "core_smiles": "CC"}
    requested = profiles.request_recommendations(config)
    assert "recommendation_profile" not in config
    assert profiles.request_recommendations(requested, profiles.DEFAULT_RECOMMENDATION_PROFILE) == requested
    loaded = profiles.apply_profile_defaults(requested)
    result = profiles.profile_provenance(loaded)
    assert result["status"] == "custom_unvalidated"
    assert set(result["material_overrides"]) == {"minimum_scores.hit_shape"}
    assert set(result["resource_overrides"]) == {"max_results"}
    assert result["additional_selection"]["core_smiles"] == "CC"
    assert "hit_shape>=0.8" in profiles.profile_label(loaded)
    requested["minimum_scores"]["hit_shape"] = 0.9
    assert config["minimum_scores"]["hit_shape"] == 0.8


@pytest.mark.parametrize(
    "mutation",
    [
        "name",
        "commit",
        "predecessor_name",
        "predecessor_hash",
        "missing_role",
        "duplicate_role",
        "unknown_role",
        "digest",
        "label",
        "extra",
        "missing",
        "fit",
    ],
)
def test_synthetic_rebound_artifact_rejects_invalid_revision(tmp_path, monkeypatch, mutation):
    # Rebound bytes test schema/value checks only; this is not packaging authority.
    artifact = profiles.load_profile(V3)
    revision = artifact["provenance"]["source_revision"]
    if mutation == "name":
        artifact["name"] = profiles.SELECTIVE_PROFILE
    elif mutation == "commit":
        revision["guarded_source_commit"] = "a" * 40
    elif mutation.startswith("predecessor_"):
        revision["predecessor"]["name" if mutation.endswith("name") else "artifact_sha256"] = "wrong"
    elif mutation == "missing_role":
        revision["source_reviews"].pop()
    elif mutation == "duplicate_role":
        revision["source_reviews"][1] = copy.deepcopy(revision["source_reviews"][0])
    elif mutation == "unknown_role":
        revision["source_reviews"][0]["role"] = "calibration"
    elif mutation in {"digest", "label"}:
        revision["source_reviews"][0]["record_sha256" if mutation == "digest" else "record_label"] = "wrong"
    elif mutation == "extra":
        revision["self_hash"] = "a" * 64
    elif mutation == "missing":
        del artifact["provenance"]["source_revision"]
    else:
        artifact["fit_sha256"] = "a" * 64
    raw = json.dumps(artifact).encode()
    path = tmp_path / "SYNTHETIC.json"
    path.write_bytes(raw)
    monkeypatch.setitem(profiles.PROFILE_SHA256, V3, hashlib.sha256(raw).hexdigest())
    monkeypatch.setattr(profiles, "_resource", lambda name: path)
    with pytest.raises(ValueError):
        profiles.load_profile(V3)


def test_load_has_no_historical_files_dependency_and_rejects_tampering(tmp_path, monkeypatch):
    raw = profiles._resource(V3).read_bytes()
    path = tmp_path / "portable-package.json"
    path.write_bytes(raw)
    calls = []

    def resource(name):
        calls.append(name)
        assert name == V3
        return path

    monkeypatch.setattr(profiles, "_resource", resource)
    assert profiles.load_profile(V3)["name"] == V3
    assert calls == [V3]
    path.write_bytes(raw + b" ")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        profiles.load_profile(V3)
    path.unlink()
    with pytest.raises(ValueError, match="Missing packaged"):
        profiles.load_profile(V3)
    with pytest.raises(ValueError, match="Unknown"):
        profiles.load_profile("qualified-selective-engineering-unknown")


@pytest.mark.parametrize(
    "name,digest_name", [(V3, V3), (V3, profiles.SELECTIVE_PROFILE), (profiles.SELECTIVE_PROFILE, V3)]
)
def test_benchmark_authority_remains_historical(name, digest_name):
    from benchmarks.research.qualified_adapter import benchmark_profile

    old = {
        "schema": "qualified-benchmark-protocol-v2",
        "recommendation_profile": profiles.SELECTIVE_PROFILE,
        "profile_artifact_sha256": profiles.PROFILE_SHA256[profiles.SELECTIVE_PROFILE],
    }
    assert benchmark_profile(old) == profiles.SELECTIVE_PROFILE
    with pytest.raises(ValueError, match="binding"):
        benchmark_profile(
            old
            | {
                "recommendation_profile": name,
                "profile_artifact_sha256": profiles.PROFILE_SHA256[digest_name],
            }
        )


@pytest.mark.parametrize("key,value", [("ranking_policy", "joint-v1"), ("route_search_policy", "design-v1")])
def test_material_ranking_and_search_overrides_do_not_become_resources(key, value):
    config = effective()
    config[key] = value
    provenance = profiles.profile_provenance(config)
    assert provenance["status"] == "custom_unvalidated"
    assert provenance["material_overrides"][key]["effective"] == value
    assert key not in provenance["resource_overrides"]


def test_landmarks_remain_additional_selection_not_material_policy_change():
    config = effective()
    config["similarity_constraints"]["landmarks"] = [
        {"id": "synthetic-oxygen", "family": "Acceptor", "hit_atom_ids": [2], "max_distance": 1.0}
    ]
    provenance = profiles.profile_provenance(config)
    assert not provenance["material_overrides"]
    assert provenance["additional_selection"]["similarity_constraints"] == config["similarity_constraints"]
