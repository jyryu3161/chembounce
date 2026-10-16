"""Profile integration using explicitly synthetic temporary fit artifacts only."""

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


def test_integrated_cli_profile_uses_refill_and_cannot_disable_synthesis(
    tmp_path, assets, proposal, synthetic_profile, capsys
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    import chembounce

    config = configuration(tmp_path, assets, proposal)
    config.pop("recommendation_profile")
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
        profiles.PROFILE,
    ]
    assert chembounce.main(args + ["--no-synllama"]) == 2
    assert not out.exists()
    assert chembounce.main(args) == 0
    assert "development_protocol_matched" in capsys.readouterr().out
    assert read(out, "run_manifest.json")["synthesis_design_supply"]["effective_overall_max_n"] == 4
    manifest = read(out / "synthesis", "manifest.json")
    assert manifest["recommendation_profile"]["status"] == "development_protocol_matched"
    assert manifest["effective_conformers"] == profiles.FIXED_CONFORMERS
    assert read(out / "synthesis", "results.json")["results"][0]["recommendation_profile"] == profiles.PROFILE


def test_wheel_includes_profile_resource_without_installing(tmp_path):
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
    resource = build / "chembounce_synthesis/data" / (profiles.PROFILE + ".json")
    resource.parent.mkdir(exist_ok=True)
    if not resource.exists():
        # Only a temporary packaging fixture while the empirical asset is unavailable.
        resource.write_text('{"synthetic_packaging_fixture_only":true}\n')
    expected = hashlib.sha256(resource.read_bytes()).hexdigest()
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
            "from importlib.resources import files; import hashlib; print(hashlib.sha256(files('chembounce_synthesis').joinpath('data','qualified-development-v1.json').read_bytes()).hexdigest())",
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )
    assert probe.returncode == 0 and probe.stdout.strip() == expected, probe.stderr
