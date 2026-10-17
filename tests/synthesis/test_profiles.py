"""Profile integration: reviewed release invariants and labelled synthetic overrides."""

import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile
from argparse import Namespace

import pytest
from rdkit import Chem

from chembounce_synthesis import profiles, workflow
from chembounce_synthesis.cli import main
from chembounce_synthesis.continuation import continuation_inputs
from chembounce_synthesis.qualification import effective_policy

ROOT = Path(__file__).resolve().parents[2]
HIT, DESIGN = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"


@pytest.fixture
def synthetic_profile(tmp_path, monkeypatch):
    current = profiles.numerical_runtime()
    artifact = {
        "schema": "qualified-development-reviewed-profile-v1",
        "status": "reviewed_development_fit",
        "review": {
            key: "c" * 64 for key in ("candidate_sha256", "root_review_sha256", "independent_review_sha256")
        },
        "scope": "previously_inspected_development_only",
        "minimum_scores": {key: 0.01 for key in profiles.METRICS},
        "calibration_spec": {
            "version": "SYNTHETIC TEST ONLY; NOT A FIT",
            "seeds": [42, 43, 44],
            "conformers": {k: v for k, v in profiles.FIXED_CONFORMERS.items() if k != "seed"},
            "fingerprint": {"kind": "binary_morgan", "radius": 2, "bits": 2048, "chirality": True},
        },
        "fit_sha256": "a" * 64,
        "provenance": {
            "runtime": {k: v for k, v in current.items() if k != "source_sha256"},
            "source_sha256": current["source_sha256"],
            "procedure_sha256": "b" * 64,
            "counts": {"synthetic_fixture_only": True},
        },
        "limitations": ["SYNTHETIC TEST VALUES; NEVER PRODUCTION OR EMPIRICAL CALIBRATION"],
    }
    path = tmp_path / "synthetic-profile.json"
    path.write_text(json.dumps(artifact))
    monkeypatch.setitem(
        profiles.PROFILE_SHA256, profiles.PROFILE, hashlib.sha256(path.read_bytes()).hexdigest()
    )
    monkeypatch.setattr(profiles, "_resource", lambda name: path)
    return path, artifact


@pytest.fixture(autouse=True)
def controlled_profile_forward(strict_predictor):
    return strict_predictor


def configuration(tmp_path, assets, proposal):
    path = tmp_path / "saved.json"
    path.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": DESIGN, "response": proposal}]})
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "recommendation_profile": profiles.PROFILE,
        "synthesis_evidence": {
            "predictor": {"kind": "reactiont5", "model_path": str(tmp_path / "controlled-model")}
        },
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "max_designs": 1,
        "max_design_attempts": 4,
        "max_scored_products": 1,
        "max_results": 1,
        "refinement": {"top_n": 1, "max_total_conformers": 24},
    }


def read(out, name):
    return json.loads((out / name).read_text())


def test_unavailable_or_corrupt_profile_never_silently_falls_back(tmp_path, monkeypatch):
    monkeypatch.setitem(profiles.PROFILE_SHA256, profiles.PROFILE, None)
    with pytest.raises(ValueError, match="unavailable"):
        profiles.apply_profile_defaults({"recommendation_profile": profiles.PROFILE})
    with pytest.raises(ValueError, match="Unknown"):
        profiles.load_profile("invented")
    path = tmp_path / "bad.json"
    path.write_text("{}")
    monkeypatch.setattr(profiles, "_resource", lambda name: path)
    monkeypatch.setitem(profiles.PROFILE_SHA256, profiles.PROFILE, "0" * 64)
    with pytest.raises(ValueError, match="SHA-256"):
        profiles.load_profile(profiles.PROFILE)
    monkeypatch.setitem(
        profiles.PROFILE_SHA256, profiles.PROFILE, hashlib.sha256(path.read_bytes()).hexdigest()
    )
    with pytest.raises(ValueError, match="positive minima"):
        profiles.load_profile(profiles.PROFILE)


def test_profile_precedes_research_defaults_and_explicit_options_precede_validation(
    tmp_path, assets, proposal, synthetic_profile
):
    config = configuration(tmp_path, assets, proposal)
    config["protocol"] = "research-v1"
    configured = workflow._config(config)
    assert configured["conformers"] == profiles.FIXED_CONFORMERS
    assert configured["refinement"]["num_conformers"] == 8
    assert configured["refinement"]["adaptive_max_conformers"] is None
    assert configured["refinement"]["max_score_range"] is None
    assert configured["refinement"]["max_metric_range"] is None
    assert effective_policy(configured)["calibration"]["status"] == "development_protocol_matched"
    assert workflow._config(configured) == configured
    config.pop("recommendation_profile")
    config["recommendation"] = {"mode": "qualified"}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    # The profile flag is applied before missing-minimum validation or generic defaults.
    loaded = workflow.load_config(path, options=Namespace(recommendation_profile=profiles.PROFILE))
    assert loaded["minimum_scores"] == synthetic_profile[1]["minimum_scores"]
    assert loaded["refinement"]["adaptive_max_conformers"] is None
    assert profiles.profile_provenance(loaded)["status"] == "development_protocol_matched"


@pytest.mark.parametrize(
    "patch",
    [
        {"minimum_scores": {"hit_shape": 0.02}},
        {"conformers": {"num_conformers": 12}},
        {"conformers": {"seed": 99}},
        {"conformers": {"strain_penalty": 0.1}},
        {"conformers": {"scoring_objective": "joint-v1"}, "ranking_policy": "joint-v1"},
        {"refinement": {"seeds": [42, 43]}},
        {"refinement": {"top_n": 0}},
        {"refinement": {"max_total_conformers": 1}},
        {"refinement": {"adaptive_max_conformers": 16}},
        {"refinement": {"max_metric_range": 0.1}},
        {"recommendation": {"mode": "exploratory"}},
        {"ranking_policy": "joint-v1"},
        {"route_search_policy": "design-v1"},
    ],
)
def test_material_overrides_cannot_keep_profile_match(tmp_path, assets, proposal, synthetic_profile, patch):
    config = configuration(tmp_path, assets, proposal)
    config.update(patch)
    provenance = effective_policy(workflow._config(config))["calibration"]
    assert provenance["status"] == "custom_unvalidated"
    assert provenance["material_overrides"]
    assert provenance["artifact_sha256"] == profiles.PROFILE_SHA256[profiles.PROFILE]


@pytest.mark.parametrize("kind", ["packages", "python", "feature", "source"])
def test_runtime_and_numerical_asset_drift_is_unvalidated(
    tmp_path, assets, proposal, synthetic_profile, monkeypatch, kind
):
    config = workflow._config(configuration(tmp_path, assets, proposal))
    current = profiles.numerical_runtime()
    if kind == "packages":
        current["packages"]["rdkit"] = "different"
    elif kind == "python":
        current["python"] = "different"
    elif kind == "feature":
        current["pharmacophore_feature_definition"]["sha256"] = "0" * 64
    else:
        current["source_sha256"]["chembounce_synthesis/conformers.py"] = "0" * 64
    monkeypatch.setattr(profiles, "numerical_runtime", lambda: current)
    provenance = profiles.profile_provenance(config)
    assert provenance["status"] == "custom_unvalidated"
    assert provenance["runtime_source_mismatches"]


def test_bound_pose_and_pocket_override_scope_and_resources_are_explicit(
    tmp_path, assets, proposal, synthetic_profile
):
    config = workflow._config(configuration(tmp_path, assets, proposal))
    assert profiles.profile_provenance(config)["resource_overrides"]
    assert profiles.profile_provenance(config)["status"] == "development_protocol_matched"
    config.update(reference_sdf="declared-bound-pose.sdf", pocket={"declared": "fixture"})
    provenance = profiles.profile_provenance(config)
    assert provenance["status"] == "custom_unvalidated"
    assert (
        "reference_sdf" in provenance["material_overrides"] and "pocket" in provenance["material_overrides"]
    )


def test_profile_api_cli_exports_and_continuation_recheck(
    tmp_path, assets, proposal, synthetic_profile, capsys
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = configuration(tmp_path, assets, proposal)
    out = tmp_path / "api"
    rows = workflow.run(HIT, [DESIGN], config, out)
    assert rows and rows[0]["calibration_status"] == "development_protocol_matched"
    config.pop("recommendation_profile")
    config["recommendation"] = {"mode": "qualified"}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    cli = tmp_path / "cli"
    assert (
        main(
            [
                "--hit",
                HIT,
                "--design",
                DESIGN,
                "--config",
                str(path),
                "--recommendation-profile",
                profiles.PROFILE,
                "--output",
                str(cli),
            ]
        )
        == 0
    )
    assert "development_protocol_matched" in capsys.readouterr().out
    selected = read(cli, "results.json")["results"][0]
    digest = profiles.PROFILE_SHA256[profiles.PROFILE]
    for directory in (out, cli):
        assert read(directory, "manifest.json")["recommendation_profile"]["artifact_sha256"] == digest
        assert (
            read(directory, "qualification.json")["policy"]["calibration"]["status"]
            == "development_protocol_matched"
        )
        with (directory / "results.tsv").open() as stream:
            table = list(csv.DictReader(stream, delimiter="\t"))
        assert table[0]["recommendation_profile"] == profiles.PROFILE
        assert table[0]["calibration_status"] == "development_protocol_matched"
        assert table[0]["calibration_artifact_sha256"] == digest
        for name in ("poses_hit.sdf", "poses_design.sdf", "references_hit.sdf", "references_design.sdf"):
            mol = next(iter(Chem.SDMolSupplier(str(directory / name))))
            assert mol.GetProp("_Name") == table[0]["route_id"]
            assert mol.GetProp("calibration_artifact_sha256") == digest
            assert mol.GetProp("recommendation_profile") == profiles.PROFILE
            assert mol.GetProp("calibration_status") == "development_protocol_matched"
            assert mol.GetProp("evidence_profile") == "computational-route-support-v1"
            assert mol.GetProp("route_support_binding") == table[0]["route_support_binding"]
    hit, designs, carried = continuation_inputs(cli)
    assert carried["recommendation_profile"] == profiles.PROFILE
    carried["minimum_scores"]["hit_shape"] = 0.02
    resumed = tmp_path / "resumed"
    rerun = workflow.run(hit, designs, carried, resumed)
    assert rerun and rerun[0]["smiles"] == selected["smiles"]
    assert rerun[0]["calibration_status"] == "custom_unvalidated"
    assert read(resumed, "manifest.json")["recommendation_profile"]["material_overrides"][
        "minimum_scores.hit_shape"
    ]


@pytest.fixture(params=["synthetic", "reviewed_release", "selective_engineering"])
def integrated_profile(request):
    if request.param == "synthetic":
        request.getfixturevalue("synthetic_profile")
    return request.param


def test_integrated_cli_profile_uses_refill_and_cannot_disable_synthesis(
    tmp_path, assets, proposal, integrated_profile, capsys
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    import chembounce

    profile_name = (
        profiles.SELECTIVE_PROFILE if integrated_profile == "selective_engineering" else profiles.PROFILE
    )
    config = configuration(tmp_path, assets, proposal)
    config.pop("recommendation_profile")
    forward_path = config["synthesis_evidence"]["predictor"].pop("model_path")
    # Hopping can emit any pyridyl regioisomer; provide their saved routes.
    proposals = []
    for amine in ("Nc1ccncc1", "Nc1cccnc1", "Nc1ccccn1"):
        design = "CC(=O)" + amine
        proposals.append(
            {
                "design": design,
                "response": {
                    "reactions": [
                        {"reaction_template": "amide", "reactants": ["CC(=O)O", amine], "product": design}
                    ],
                    "building_blocks": ["CC(=O)O", amine],
                },
            }
        )
    Path(config["provider"]["path"]).write_text(json.dumps({"schema_version": 1, "proposals": proposals}))
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    replacements = tmp_path / "replacements.smi"
    replacements.write_text("c1ccncc1\n")
    out = tmp_path / "integrated"
    args = [
        "-i",
        HIT,
        "-o",
        str(out),
        "--fragments",
        "c1ccccc1",
        "--replace_scaffold_files",
        str(replacements),
        "--cand_max_n__rplc",
        "1",
        "--no-shape",
        "-t",
        "0",
        "--search-threshold",
        "0",
        "--synthesis-config",
        str(path),
        "--recommendation-profile",
        profile_name,
        "--forward-model-path",
        forward_path,
    ]
    assert chembounce.main(args + ["--no-synllama"]) == 2
    assert not out.exists()
    assert chembounce.main(args) == 0
    status = profiles.profile_provenance(
        workflow._config(
            {
                **config,
                "recommendation_profile": profile_name,
                "synthesis_evidence": {"predictor": {"kind": "reactiont5", "model_path": forward_path}},
            }
        )
    )["status"]
    assert status in capsys.readouterr().out
    assert read(out, "run_manifest.json")["synthesis_design_supply"]["effective_overall_max_n"] == 4
    manifest = read(out / "synthesis", "manifest.json")
    assert manifest["recommendation_profile"]["status"] == status
    assert manifest["effective_conformers"] == profiles.FIXED_CONFORMERS
    assert read(out / "synthesis", "results.json")["results"][0]["recommendation_profile"] == profile_name


@pytest.mark.parametrize(
    "profile_name", [profiles.PROFILE, profiles.SELECTIVE_PROFILE, profiles.GUARDED_SELECTIVE_PROFILE]
)
def test_wheel_includes_profile_resource_without_installing(tmp_path, profile_name):
    pytest.importorskip("setuptools")
    pytest.importorskip("wheel")
    build = tmp_path / "source"
    build.mkdir()
    for name in ("pyproject.toml", "README.md"):
        shutil.copyfile(ROOT / name, build / name)
    for path in ROOT.glob("*.py"):
        shutil.copyfile(path, build / path.name)
    shutil.copytree(
        ROOT / "chembounce_synthesis",
        build / "chembounce_synthesis",
        ignore=shutil.ignore_patterns("__pycache__"),
    )
    resource = build / "chembounce_synthesis/data" / (profile_name + ".json")
    assert resource.exists(), "The reviewed real profile must ship in the package"
    expected = hashlib.sha256(resource.read_bytes()).hexdigest()
    assert expected == profiles.PROFILE_SHA256[profile_name]
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run(
        [sys.executable, "-c", "import setuptools.build_meta; setuptools.build_meta.build_wheel('dist')"],
        cwd=build,
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    wheel_path = next((build / "dist").glob("*.whl"))
    with zipfile.ZipFile(wheel_path) as archive:
        assert (
            hashlib.sha256(archive.read("chembounce_synthesis/data/" + resource.name)).hexdigest() == expected
        )
    env["PYTHONPATH"] = str(wheel_path)
    probe = subprocess.run(
        [
            sys.executable,
            "-c",
            f"from importlib.resources import files; import hashlib; print(hashlib.sha256(files('chembounce_synthesis').joinpath('data',{profile_name + '.json'!r}).read_bytes()).hexdigest())",
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )
    assert probe.returncode == 0 and probe.stdout.strip() == expected, probe.stderr
    loaded = subprocess.run(
        [
            sys.executable,
            "-c",
            f"from chembounce_synthesis.profiles import load_profile; print(load_profile({profile_name!r})['status'])",
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )
    assert loaded.returncode == 0 and loaded.stdout.strip() == (
        "declared_engineering_policy"
        if profile_name in profiles.ENGINEERING_PROFILES
        else "reviewed_development_fit"
    ), loaded.stderr


def test_reviewed_release_reconstructs_the_exact_unchanged_fit_candidate():
    artifact = profiles.load_profile(profiles.PROFILE)
    candidate = json.loads(json.dumps(artifact))
    review = candidate.pop("review")
    candidate["schema"] = "qualified-development-preset-candidate-v1"
    candidate["status"] = "development_fit_requires_independent_review"
    candidate["limitations"] = artifact["limitations"][:5] + [
        "Independent review and root stage gate remain required before publishing a preset"
    ]
    raw = (json.dumps(candidate, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    assert (
        hashlib.sha256(raw).hexdigest()
        == review["candidate_sha256"]
        == "56cc18999701745b89eca83a5798d3a5f67d1a2789e80d8b6c87aff3b6e7fad7"
    )
    assert artifact["fit_sha256"] == "2294a4d86c38e2e3c5fede92ab5095f6f563693ee1fd76f04c7a19ee74aae77d"
    assert review["root_review_sha256"] == "264cfaaae02039680f7f12931eef6c9b2e4c8b948c1ac9cc27a5c4a7ad46e39b"
    assert (
        review["independent_review_sha256"]
        == "532faeb3c6e14929633f56063cdd579d7ed9dda7243b3d8835d1ceb0284f7a67"
    )
    assert (
        artifact["provenance"]["procedure_sha256"]
        == "979302be7fe2aa667f70404c3e538335ab893ba3584388459ae350b1c8ed257f"
    )
    assert artifact["provenance"]["counts"] == {
        "active": 111,
        "inactive": 86,
        "pairs": 197,
        "queries": 12,
        "targets": 6,
        "trials": 591,
    }
    assert artifact["final_development_report"]["computational_failure_pairs"] == 3
    assert artifact["oof_development_report"]["computational_failure_pairs"] == 3


def test_real_profile_reports_stereo_ingress_source_mismatch_without_rebinding(tmp_path, assets, proposal):
    config = workflow._config(configuration(tmp_path, assets, proposal))
    provenance = profiles.profile_provenance(config)
    current = profiles.numerical_runtime()
    artifact = profiles.load_profile(profiles.PROFILE)
    mismatched = {
        path
        for path in current["source_sha256"]
        if current["source_sha256"][path] != artifact["provenance"]["source_sha256"][path]
    }
    assert mismatched == {"chembounce_synthesis/common.py", "chembounce_synthesis/conformers.py"}
    assert mismatched <= provenance["runtime_source_mismatches"].keys()
    assert provenance["status"] == "custom_unvalidated"
    assert config["minimum_scores"] == {
        "hit_tanimoto": 0.08695652173913043,
        "hit_shape": 0.5055221213607054,
        "hit_pharmacophore": 0.11206335132501284,
    }
    assert config["conformers"] == profiles.FIXED_CONFORMERS
    assert config["refinement"]["seeds"] == [42, 43, 44]
    assert config["refinement"]["num_conformers"] == 8
    assert config["refinement"]["adaptive_max_conformers"] is None
    assert not provenance["material_overrides"]
    expected = (
        "custom_unvalidated" if provenance["runtime_source_mismatches"] else "development_protocol_matched"
    )
    assert provenance["status"] == expected


def test_real_profile_cli_runs_without_manual_cutoffs_and_exports_review_digest(
    tmp_path, assets, proposal, capsys
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = configuration(tmp_path, assets, proposal)
    config.pop("recommendation_profile")
    forward_path = config["synthesis_evidence"]["predictor"].pop("model_path")
    path = tmp_path / "actual.json"
    path.write_text(json.dumps(config))
    out = tmp_path / "actual"
    assert (
        main(
            [
                "--hit",
                HIT,
                "--design",
                DESIGN,
                "--config",
                str(path),
                "--recommendation-profile",
                profiles.PROFILE,
                "--forward-model-path",
                forward_path,
                "--output",
                str(out),
            ]
        )
        == 0
    )
    output = capsys.readouterr().out
    assert "permissive development thresholds" in output
    assert "not strong-similarity or activity guarantees" in output
    manifest = read(out, "manifest.json")
    provenance = manifest["recommendation_profile"]
    assert provenance["artifact_sha256"] == profiles.PROFILE_SHA256[profiles.PROFILE]
    assert provenance["review"] == profiles.load_profile(profiles.PROFILE)["review"]
    rows = read(out, "results.json")["results"]
    assert len(rows) == 1 and rows[0]["recommendation_profile"] == profiles.PROFILE
    assert rows[0]["calibration_status"] == provenance["status"]
    for name in ("poses_hit.sdf", "poses_design.sdf"):
        pose = next(iter(Chem.SDMolSupplier(str(out / name))))
        assert pose.GetProp("_Name") == rows[0]["route_id"]
        assert pose.GetProp("calibration_artifact_sha256") == provenance["artifact_sha256"]


def test_selective_profile_is_engineering_and_preserves_optional_constraints(tmp_path, assets, proposal):
    from chembounce_synthesis.selective_scoring import configured_constraints

    config = configuration(tmp_path, assets, proposal)
    config.update(recommendation_profile=profiles.SELECTIVE_PROFILE, core_smiles="CC(=O)N")
    landmark = {"id": "amide-oxygen", "family": "Acceptor", "hit_atom_ids": [2], "max_distance": 1.0}
    config["similarity_constraints"] = {"landmarks": [landmark]}
    effective = workflow._config(config)
    assert effective["similarity_constraints"] == profiles.SELECTIVE_CONSTRAINTS | {"landmarks": [landmark]}
    assert effective["core_smiles"] == "CC(=O)N"
    bound = configured_constraints(effective, HIT)
    assert bound.policy["core_smiles"] == "CC(N)=O"
    provenance = profiles.profile_provenance(effective)
    assert provenance["status"] == (
        "custom_unvalidated" if provenance["runtime_source_mismatches"] else "engineering_policy_matched"
    )
    assert not provenance["material_overrides"]
    assert not {"fit_sha256", "review", "calibration_spec", "calibration_provenance"} & provenance.keys()
    assert set(provenance["source_match_scope"]) == {
        "chembounce_synthesis/" + name for name in profiles.KERNEL_FILES + profiles.SELECTIVE_KERNEL_FILES
    }
    assert set(profiles.numerical_runtime()["source_sha256"]) == {
        "chembounce_synthesis/" + name for name in profiles.KERNEL_FILES
    }
    assert "permissive" not in profiles.profile_label(effective)
    assert "not fitted/calibrated" in profiles.profile_label(effective)


@pytest.mark.parametrize(
    "change,field",
    [
        ({"minimum_scores": {"hit_shape": 0.8}}, "minimum_scores.hit_shape"),
        ({"similarity_constraints": {"polar_recall": 0.8}}, "similarity_constraints.polar_recall"),
        ({"route_search_policy": "design-v1"}, "route_search_policy"),
        ({"ranking_policy": "joint-v1"}, "ranking_policy"),
    ],
)
def test_selective_profile_material_overrides_survive_and_are_custom(
    tmp_path, assets, proposal, change, field
):
    config = configuration(tmp_path, assets, proposal)
    config.update(recommendation_profile=profiles.SELECTIVE_PROFILE, **change)
    effective = workflow._config(config)
    provenance = profiles.profile_provenance(effective)
    assert provenance["status"] == "custom_unvalidated"
    assert field in provenance["material_overrides"]
    for section, value in change.items():
        assert (
            effective[section] == value
            if not isinstance(value, dict)
            else all(effective[section][k] == v for k, v in value.items())
        )


@pytest.mark.parametrize(
    "change",
    [
        {"conformers": {"num_conformers": 9}},
        {"conformers": {"reject_unspecified_stereo": False}},
        {"refinement": {"seeds": [42, 43]}},
        {"similarity_constraints": None},
        {"synthesis_evidence": {"profile": "legacy-conditional-v1"}},
    ],
)
def test_selective_profile_unsupported_overrides_reject(tmp_path, assets, proposal, change):
    config = configuration(tmp_path, assets, proposal)
    config.update(recommendation_profile=profiles.SELECTIVE_PROFILE, **change)
    with pytest.raises(ValueError):
        workflow._config(config)


@pytest.mark.parametrize("kind", ["runtime", "selective_scoring.py", "similarity_constraints.py"])
def test_selective_profile_binds_new_sources_and_runtime(tmp_path, assets, proposal, monkeypatch, kind):
    config = configuration(tmp_path, assets, proposal)
    config["recommendation_profile"] = profiles.SELECTIVE_PROFILE
    effective = workflow._config(config)
    current = profiles.numerical_runtime(profiles.SELECTIVE_PROFILE)
    if kind == "runtime":
        current["packages"]["rdkit"] = "changed"
    else:
        current["source_sha256"]["chembounce_synthesis/" + kind] = "0" * 64
    monkeypatch.setattr(profiles, "numerical_runtime", lambda name=None: current)
    assert profiles.profile_provenance(effective)["status"] == "custom_unvalidated"


def test_named_selective_api_and_synthesis_cli_use_same_final_contract(tmp_path, assets, proposal, capsys):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = configuration(tmp_path, assets, proposal)
    config["recommendation_profile"] = profiles.SELECTIVE_PROFILE
    expected = profiles.profile_provenance(workflow._config(config))["status"]
    api = workflow.run(HIT, [DESIGN], config, tmp_path / "api-selective")
    assert len(api) == 1 and api[0]["calibration_status"] == expected
    config.pop("recommendation_profile")
    path = tmp_path / "cli-selective.json"
    path.write_text(json.dumps(config))
    out = tmp_path / "cli-selective"
    assert (
        main(
            [
                "--hit",
                HIT,
                "--design",
                DESIGN,
                "--config",
                str(path),
                "--recommendation-profile",
                profiles.SELECTIVE_PROFILE,
                "--output",
                str(out),
            ]
        )
        == 0
    )
    label = capsys.readouterr().out
    assert expected in label and "permissive" not in label and "not fitted/calibrated" in label
    row = read(out, "results.json")["results"][0]
    assert row["similarity_witness_sha256"] == api[0]["similarity_witness_sha256"]
    assert row["recommendation_profile"] == profiles.SELECTIVE_PROFILE
    assert [t["seed"] for t in read(out, "refinement.json")["candidates"][0]["trials"]] == [42, 43, 44]
