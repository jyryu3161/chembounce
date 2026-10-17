"""V5 engineering packaging and controlled integration, not scientific yield."""

import argparse
import copy
import hashlib
import json
import subprocess
import sys
from types import SimpleNamespace

import pytest

from chembounce_synthesis import profiles
from chembounce_synthesis.options import add_expansion_options, apply_expansion_options

V5 = profiles.LANDMARK_SELECTIVE_PROFILE
HISTORICAL = {
    profiles.PROFILE: "9d3c60a895912b28c13a29aac0a04f385ddb8a0effc9fcd81d64a50a3d48bdd2",
    profiles.SELECTIVE_PROFILE: "a8b47a409a21518dd482acbefdce38774414ae0b65f28d6ea7cb52ab63e1ef52",
}
# Explicit historical production bytes, never synthetic runtime rebinding.
HISTORICAL.update(
    {
        profiles.GUARDED_SELECTIVE_PROFILE: "fd5fb3b65b12307ed1a9b5371e74427eb709d8cf3012bbcba5095022ddb34537",
        profiles.DIAGNOSTIC_SELECTIVE_PROFILE: "437b5c11f9b5627122e2ad70ecbc2d4f9f8b2e5294c6fbe8be81058ffb8e48b5",
    }
)


def config(**kwargs):
    return profiles.apply_profile_defaults({"recommendation_profile": V5, **kwargs})


def test_packaged_v5_exact_policy_runtime_and_history():
    old = profiles.load_profile(profiles.DIAGNOSTIC_SELECTIVE_PROFILE)
    new = profiles.load_profile(V5)
    for key in ("engineering_spec", "minimum_scores", "schema", "status", "scope"):
        assert new[key] == old[key]
    assert new["provenance"]["runtime"] == old["provenance"]["runtime"]
    assert (
        new["provenance"]["source_revision"]["inherited_diagnostic_revision"]
        == old["provenance"]["source_revision"]
    )
    assert profiles.SELECTIVE_PROFILE == "qualified-selective-engineering-v2"
    assert profiles.DEFAULT_RECOMMENDATION_PROFILE == V5
    assert not {"review", "calibration_spec", "fit_sha256"} & new.keys()
    for name, digest in HISTORICAL.items():
        assert profiles.PROFILE_SHA256[name] == digest
        assert hashlib.sha256(profiles._resource(name).read_bytes()).hexdigest() == digest
        result = profiles.profile_provenance(
            profiles.apply_profile_defaults({"recommendation_profile": name})
        )
        assert result["status"] == "custom_unvalidated"
    observed = profiles.numerical_runtime(V5)
    assert observed["source_sha256"] == new["provenance"]["source_sha256"]
    assert len(observed["source_sha256"]) == 6
    result = profiles.profile_provenance(config())
    assert not result["material_overrides"]
    assert not any(k.startswith("chembounce_synthesis/") for k in result["runtime_source_mismatches"])
    assert result["status"] == (
        "custom_unvalidated" if result["runtime_source_mismatches"] else "engineering_policy_matched"
    )
    diagnostic = result["optional_diagnostics"]["rejected_pose_diagnostics"]
    # Pure route-computation reuse changes workflow source without rebinding
    # the historical optional diagnostic declaration or numerical policy.
    key = "chembounce_synthesis/workflow.py"
    assert diagnostic["source_mismatches"] == {
        key: {
            "declared": diagnostic["source_sha256"][key],
            "current": diagnostic["current_source_sha256"][key],
        }
    }
    enabled = profiles.profile_provenance(config(rejected_pose_diagnostics=True))
    assert enabled["status"] == "custom_unvalidated"
    assert enabled["runtime_source_mismatches"] == result["runtime_source_mismatches"]
    assert (
        enabled["optional_diagnostics"]["rejected_pose_diagnostics"]["source_mismatches"]
        == diagnostic["source_mismatches"]
    )
    assert result["optional_landmark_proposals"]["enabled"] is False


def test_omitted_false_and_disabled_import_capability_independence():
    code = """
import sys
from chembounce_synthesis import profiles
from chembounce_synthesis import conformers
conformers.rdShapeAlign = None
from chembounce_synthesis.selective_scoring import configured_constraints
cfg=profiles.apply_profile_defaults(profiles.request_recommendations({}))
a=configured_constraints(cfg, "CCN")
explicit=profiles.apply_profile_defaults(profiles.request_recommendations({"similarity_constraints":{"landmark_pose_proposals":False}}))
b=configured_constraints(explicit, "CCN")
assert cfg == explicit and a.policy == b.policy and a.sha256 == b.sha256
assert "landmark_pose_proposals" not in a.policy
profiles.profile_provenance(cfg)
assert "chembounce_synthesis.landmark_proposals" not in sys.modules
"""
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("old", list(HISTORICAL))
@pytest.mark.parametrize("origin", ["config", "cli"])
def test_no_concealed_upgrade(old, origin):
    raw = {"recommendation_profile": old} if origin == "config" else {}
    before = copy.deepcopy(raw)
    args = argparse.Namespace(recommend=True, recommendation_profile=V5 if origin == "config" else old)
    with pytest.raises(ValueError, match="conflicts"):
        apply_expansion_options(raw, args)
    assert raw == before
    assert (
        apply_expansion_options({}, argparse.Namespace(recommendation_profile=old))["recommendation_profile"]
        == old
    )


def test_explicit_options_copy_and_material_search_semantics(monkeypatch):
    from chembounce_synthesis import conformers, selective_scoring as scorer
    from .test_landmark_proposals import FIXTURES

    monkeypatch.setattr(
        conformers, "rdShapeAlign", SimpleNamespace(ScoreMol=lambda *a: None, ShapeInputOptions=lambda: None)
    )
    parser = argparse.ArgumentParser()
    add_expansion_options(parser)
    raw = {"similarity_constraints": {"landmarks": FIXTURES[2]["landmarks"]}}
    before = copy.deepcopy(raw)
    enabled = profiles.apply_profile_defaults(
        apply_expansion_options(raw, parser.parse_args(["--recommend", "--landmark-pose-proposals"]))
    )
    bound = scorer.configured_constraints(enabled, FIXTURES[2]["smiles"])
    result = profiles.profile_provenance(enabled)
    assert result["status"] == "custom_unvalidated"
    assert set(result["material_overrides"]) == {"similarity_constraints.landmark_pose_proposals"}
    assert (
        bound.policy["landmark_pose_proposals"]["source_sha256"]
        == result["optional_landmark_proposals"]["source_sha256"][profiles.LANDMARK_FILE]
    )
    assert bound.policy["landmark_pose_proposals"]["limits"] == profiles.LANDMARK_DECLARATION["limits"]
    assert raw == before
    disabled = apply_expansion_options(enabled, parser.parse_args(["--no-landmark-pose-proposals"]))
    assert disabled["similarity_constraints"]["landmark_pose_proposals"] is False
    assert not profiles.profile_provenance(disabled)["material_overrides"]
    assert apply_expansion_options(enabled, parser.parse_args([])) == enabled


@pytest.mark.parametrize(
    "field",
    [
        "helper_hash",
        "missing_helper",
        "default",
        "limit",
        "bool_limit",
        "parity",
        "revision",
        "predecessor",
        "extra",
        "role",
    ],
)
def test_rebound_schema_tamper(tmp_path, monkeypatch, field):
    a = profiles.load_profile(V5)
    d = a["optional_landmark_proposals"]
    if field == "helper_hash":
        d["source_sha256"][profiles.LANDMARK_FILE] = "0" * 64
    elif field == "missing_helper":
        a["provenance"]["source_sha256"].pop(profiles.LANDMARK_FILE)
    elif field == "default":
        d["default_enabled"] = True
    elif field == "limit":
        d["limits"]["alternatives_per_call"] += 1
    elif field == "bool_limit":
        d["limits"]["cooperative_seconds_per_call"] = True
    elif field == "parity":
        d["parity_tolerance"] = 0.1
    elif field == "revision":
        a["provenance"]["source_revision"]["landmark_source_commit"] = "0" * 40
    elif field == "predecessor":
        a["provenance"]["source_revision"]["predecessor"]["name"] = V5
    elif field == "role":
        a["provenance"]["source_revision"]["source_reviews"][0]["role"] = "calibration"
    else:
        d["self_sha256"] = "0" * 64
    raw = json.dumps(a).encode()
    path = tmp_path / "schema-only.json"
    path.write_bytes(raw)
    monkeypatch.setattr(profiles, "_resource", lambda name: path)
    monkeypatch.setitem(profiles.PROFILE_SHA256, V5, hashlib.sha256(raw).hexdigest())
    with pytest.raises(ValueError, match="engineering profile"):
        profiles.load_profile(V5)


def test_actual_helper_drift_invalidates_even_disabled(monkeypatch):
    observed = profiles.numerical_runtime(V5)
    observed["source_sha256"][profiles.LANDMARK_FILE] = "0" * 64
    monkeypatch.setattr(profiles, "numerical_runtime", lambda name=None: observed)
    p = profiles.profile_provenance(config())
    assert p["status"] == "custom_unvalidated" and profiles.LANDMARK_FILE in p["runtime_source_mismatches"]


def test_portable_resource_byte_binding_and_no_historical_path_reads(tmp_path, monkeypatch):
    raw = profiles._resource(V5).read_bytes()
    path = tmp_path / "packaged.json"
    path.write_bytes(raw)
    monkeypatch.setattr(profiles, "_resource", lambda name: path)
    assert profiles.load_profile(V5)["name"] == V5
    path.write_bytes(raw + b" ")
    with pytest.raises(ValueError, match="SHA-256"):
        profiles.load_profile(V5)


def test_early_missing_api_and_no_landmarks(tmp_path, monkeypatch):
    from chembounce_synthesis import workflow, conformers

    def forbidden(*a, **k):
        pytest.fail("asset/model work before configuration refusal")

    monkeypatch.setattr(workflow, "Catalog", forbidden)
    monkeypatch.setattr(workflow, "Rules", forbidden)
    base = {
        "schema_version": 1,
        "catalog": "unused",
        "rules": "unused",
        "provider": {"kind": "saved", "path": "unused"},
        "recommendation_profile": V5,
        "synthesis_evidence": {"predictor": {"kind": "reactiont5", "model_path": "unused"}},
        "similarity_constraints": {"landmark_pose_proposals": True},
    }
    with pytest.raises(ValueError, match="nonempty landmarks"):
        workflow.run("CCN", [], base, tmp_path / "absent")
    base["similarity_constraints"]["landmarks"] = [
        {"id": "N", "family": "Donor", "hit_atom_ids": [2], "max_distance": 0.25}
    ]
    monkeypatch.setattr(conformers, "rdShapeAlign", None)
    with pytest.raises(ValueError, match="ScoreMol"):
        workflow.run("CCN", [], base, tmp_path / "absent2")
    # Generic run persists its existing zero-work error bundle, unlike public
    # recommendation preflight which rejects before dispatch.
    for directory in ("absent", "absent2"):
        summary = json.loads((tmp_path / directory / "summary.json").read_text())
        assert summary["status"] == "error" and summary["scored_products"] == 0
        assert summary["verified_routes"] == 0 and summary["proposals"] == 0
        assert json.loads((tmp_path / directory / "results.json").read_text())["results"] == []


@pytest.mark.parametrize("enabled", [False, True])
def test_controlled_named_workflow_all_seed_saved_report_and_continuation(
    tmp_path, monkeypatch, assets, proposal, strict_predictor, enabled
):
    """Existing known amide; manual coords/controlled shape and forward beams.

    Full workflow/qualification/export/report executes, not a model/yield study.
    """
    import math
    from rdkit import Chem
    from chembounce_synthesis import conformers, selective_scoring as scorer, workflow, report
    from chembounce_synthesis.continuation import continuation_inputs
    from chembounce_synthesis.conformers import Ensemble
    from .test_selective_scoring import pose
    from .test_selective_workflow import config_for, HIT, DESIGN, read, assert_seed_geometry

    calls = []
    monkeypatch.setattr(workflow, "prepare_reference", lambda *a: None)
    monkeypatch.setattr(scorer, "require_alignment", lambda: None)
    monkeypatch.setattr(
        scorer, "ensemble", lambda smiles, settings: Ensemble(pose(smiles), {0: 0.0}, {0: 0.0}, "MMFF94s", {})
    )

    def overlay(ref, candidate, cid, settings):
        calls.append(settings.seed)
        aligned = Chem.Mol(candidate.mol)
        if enabled and len(calls) % 2:
            for i, point in enumerate(aligned.GetConformer(cid).GetPositions()):
                aligned.GetConformer(cid).SetAtomPosition(i, (point + [0, 2, 0]).tolist())
        pharm, counts = conformers.pharmacophore_overlap(
            ref.mol, aligned, 0, cid, settings.pharmacophore_sigma
        )
        joint = math.sqrt(0.95 * pharm)
        return {
            **counts,
            "shape": 0.95,
            "alignment_color": 0.7,
            "pharmacophore": pharm,
            "joint": joint,
            "reference_adjusted_joint": joint,
            "reference_conformer_id": 0,
            "reference_relative_energy": 0.0,
        }, aligned

    monkeypatch.setattr(scorer, "_overlay", overlay)
    fixed = []

    def score(ref, cand, *args):
        fixed.append(1)
        return (0.95, 0.7) if cand.GetConformer().GetAtomPosition(0).y > 1 else (0.75, 0.7)

    monkeypatch.setattr(
        conformers,
        "rdShapeAlign",
        SimpleNamespace(ScoreMol=score, ShapeInputOptions=lambda: SimpleNamespace()),
    )
    c = config_for(tmp_path, assets, proposal)
    c["recommendation_profile"] = V5
    c["synthesis_evidence"] = {
        "predictor": {"kind": "reactiont5", "model_path": str(tmp_path / "controlled-forward")}
    }
    c["similarity_constraints"].update(
        landmark_pose_proposals=enabled,
        landmarks=[{"id": "amide-O", "family": "Acceptor", "hit_atom_ids": [2], "max_distance": 0.25}],
    )
    out = tmp_path / "out"
    rows = workflow.run(HIT, [DESIGN], c, out)
    manifest = read(out, "manifest.json")
    provenance = manifest["recommendation_profile"]
    assert provenance["name"] == V5
    if enabled:
        assert provenance["status"] == "custom_unvalidated"
    # Minimum runtime has no concrete reaction-origin metadata, independently of ScoreMol.
    if not rows:
        assert not calls
        assert read(out, "summary.json")["status"] == "no_candidates"
        assert any("origin" in json.dumps(r).lower() for r in read(out, "routes.json")["routes"])
        return
    assert len(rows) == 1 and rows[0]["recommendation_status"] == "qualified"
    refinement = assert_seed_geometry(out, c)
    for trial in refinement["trials"]:
        assert ("landmark_pose_proposals" in trial["pose_selection"]) is enabled
        if enabled:
            assert (
                trial["hit"]["landmark_proposal"]["candidate_pose_sha256"]
                == trial["selective_pose_evidence"]["hit_candidate"]["sha256"]
            )
            assert trial["pose_selection"]["pairs"][0]["selected_origin"] == "landmark_proposal"
    assert bool(fixed) is enabled
    view = report.build_report(out)
    assert view["artifact_consistency"] == "consistent_recorded_bundle", view["issues"]
    assert view["policy"]["calibration"]["name"] == V5
    assert view["policy"]["calibration"]["status"] == provenance["status"]
    html = report.render_html(view)
    assert V5 in html and (not enabled or "custom_unvalidated" in html)
    _, designs, carried = continuation_inputs(out)
    assert (
        carried["recommendation_profile"] == V5
        and carried["similarity_constraints"]["landmark_pose_proposals"] is enabled
    )
    before = len(calls)
    second = tmp_path / "continued"
    again = workflow.run(HIT, designs, carried, second)
    assert len(calls) > before and len(again) == 1
    assert (
        read(second, "manifest.json")["recommendation_policy"]["similarity_constraints"]
        == manifest["recommendation_policy"]["similarity_constraints"]
    )
