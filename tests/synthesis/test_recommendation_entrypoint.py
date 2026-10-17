"""Explicit recommendation intent; controlled fixtures, no scientific campaign."""

import copy
import json
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import pytest

from chembounce_synthesis import recommend, conformers, profiles, workflow
from chembounce_synthesis.cli import main as synthesis_main
from chembounce_synthesis.options import apply_expansion_options

HIT, DESIGN = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"


@pytest.fixture
def config(tmp_path, assets, proposal):
    saved = tmp_path / "proposals.json"
    saved.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": DESIGN, "response": proposal}]})
    )
    model = tmp_path / "forward-model"
    model.mkdir()
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(saved)},
        "synthesis_evidence": {"predictor": {"kind": "reactiont5", "model_path": str(model)}},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "max_designs": 1,
        "max_design_attempts": 4,
        "max_scored_products": 1,
        "max_results": 1,
        "refinement": {"top_n": 1, "max_total_conformers": 24},
    }


@pytest.fixture
def no_work(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("Expensive constructor/generator invoked during request preflight")

    import chembounce

    for target, name in [
        (chembounce, "chembounce"),
        (workflow, "Catalog"),
        (workflow, "Rules"),
        (workflow, "SynLlama"),
        (workflow, "_evidence_evaluator"),
        (workflow, "prepare_reference"),
        (workflow, "ensemble"),
    ]:
        monkeypatch.setattr(target, name, forbidden)
    # Pure schema/path tests do not claim alignment support on the minimum runtime.
    monkeypatch.setattr(conformers, "require_alignment", lambda: None)


def untouched_designs():
    pytest.fail("Design iterator consumed before request validation")
    yield DESIGN


@pytest.mark.parametrize(
    "raw,cli",
    [
        (profiles.PROFILE, None),
        (profiles.GUARDED_SELECTIVE_PROFILE, None),
        (None, profiles.GUARDED_SELECTIVE_PROFILE),
        (profiles.PROFILE, profiles.SELECTIVE_PROFILE),
        (profiles.SELECTIVE_PROFILE, profiles.PROFILE),
        (None, profiles.PROFILE),
    ],
)
def test_both_original_profile_declarations_checked_before_precedence(config, no_work, raw, cli):
    if raw is not None:
        config["recommendation_profile"] = raw
    before = copy.deepcopy(config)
    with pytest.raises(ValueError, match="conflicts"):
        apply_expansion_options(config, Namespace(recommend=True, recommendation_profile=cli))
    assert config == before
    # Existing ordinary CLI precedence is unaffected.
    if cli:
        assert (
            apply_expansion_options(config, Namespace(recommendation_profile=cli))["recommendation_profile"]
            == cli
        )


def test_api_and_cli_resolve_same_config_without_consumption_mutation_or_deadline_change(
    config, no_work, monkeypatch, tmp_path
):
    calls = []

    def run(hit, designs, effective, output_dir, **kwargs):
        calls.append((hit, designs, effective, output_dir, kwargs))
        return ["controlled sentinel"]

    monkeypatch.setattr(workflow, "run", run)
    supplied = untouched_designs()
    before = copy.deepcopy(config)
    assert recommend(HIT, supplied, config, tmp_path / "api", deadline_monotonic=123.45) == [
        "controlled sentinel"
    ]
    assert len(calls) == 1 and calls[0][1] is supplied
    assert calls[0][4] == {"scoring_mode": "3d", "deadline_monotonic": 123.45}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    loaded = workflow.load_config(path, options=Namespace(recommend=True))
    assert calls[0][2] == workflow.preflight_recommendation(HIT, loaded)
    assert calls[0][2]["recommendation_profile"] == profiles.DEFAULT_RECOMMENDATION_PROFILE
    assert config == before
    assert profiles.request_recommendations(loaded) == loaded


@pytest.mark.parametrize(
    "change",
    [
        {"recommendation": {"mode": "exploratory"}},
        {"mode": "expand"},
        {"refinement": {"top_n": 0}},
        {"refinement": {"seeds": [42, 43]}},
        {"conformers": {"num_conformers": 12}},
        {"conformers": {"reject_unspecified_stereo": False}},
        {"conformers": {"scoring_objective": "joint-v1"}},
        {"refinement": {"adaptive_max_conformers": 16}},
        {"reference_sdf": "missing.sdf"},
        {"pocket": {}},
        {"synthesis_evidence": {"profile": "legacy-conditional-v1"}},
        {"max_scored_products": 0},
        {"minimum_scores": {"hit_shape": None}},
        {"minimum_scores": {"hit_shape": 0}},
        {"minimum_scores": {"hit_shape": float("nan")}},
        {"minimum_scores": {"hit_shape": float("inf")}},
        {"minimum_scores": {"hit_shape": True}},
        {"similarity_constraints": {"polar_recall": 0}},
    ],
)
def test_invalid_effective_policy_fails_before_iterator_or_assets(config, no_work, tmp_path, change):
    config.update(change)
    with pytest.raises(ValueError):
        recommend(HIT, untouched_designs(), config, tmp_path / "out")


@pytest.mark.parametrize(
    "change,field",
    [
        ({"minimum_scores": {"hit_shape": 0.8}}, "minimum_scores.hit_shape"),
        ({"similarity_constraints": {"polar_recall": 0.8}}, "similarity_constraints.polar_recall"),
        ({"route_search_policy": "design-v1"}, "route_search_policy"),
        ({"ranking_policy": "joint-v1"}, "ranking_policy"),
    ],
)
def test_explicit_material_overrides_remain_visible_and_custom(config, no_work, change, field):
    config.update(change)
    effective = workflow.preflight_recommendation(HIT, config)
    provenance = profiles.profile_provenance(effective)
    assert field in provenance["material_overrides"]
    assert provenance["status"] == "custom_unvalidated"
    label = profiles.profile_label(effective)
    for metric, value in effective["minimum_scores"].items():
        assert f"{metric}>={value}" in label
    assert f"polar_recall>={effective['similarity_constraints']['polar_recall']}" in label
    assert "custom_unvalidated" in label


def test_resource_and_additional_constraints_do_not_become_material(config, no_work):
    config.update(core_smiles="CC(=O)N", property_thresholds={"mw": [0, 400]}, expansion={"rounds": 0})
    config["similarity_constraints"] = {
        "landmarks": [{"id": "oxygen", "family": "Acceptor", "hit_atom_ids": [2]}]
    }
    effective = workflow.preflight_recommendation(HIT, config)
    p = profiles.profile_provenance(effective)
    assert not p["material_overrides"]
    assert p["resource_overrides"]
    assert p["additional_selection"]["core_smiles"] == config["core_smiles"]
    assert "search" in p["additional_selection"] and "expansion" in p["additional_selection"]


@pytest.mark.parametrize(
    "change",
    [
        {"catalog": "missing"},
        {"rules": "missing"},
        {"provider": {"kind": "unknown"}},
        {"provider": {"kind": "saved"}},
        {"provider": {"kind": "saved", "path": "missing"}},
        {"provider": {"kind": "saved", "path": True}},
        {"provider": {"kind": "synllama"}},
        {"synthesis_evidence": {"predictor": {"model_path": "missing"}}},
    ],
)
def test_invalid_provider_and_paths_fail_before_work(config, no_work, change, tmp_path):
    config.update(change)
    with pytest.raises(ValueError):
        recommend(HIT, untouched_designs(), config, tmp_path / "out")


def test_file_directory_types_and_path_authority(config, no_work, monkeypatch, tmp_path):
    parent = tmp_path / "config-dir"
    parent.mkdir()
    # Absolute '..' paths exercise JSON-parent authority, separate from caller cwd.
    raw = copy.deepcopy(config)
    raw["catalog"] = "../catalog.tsv"
    raw["rules"] = "../rules.json"
    raw["provider"]["path"] = "../proposals.json"
    raw["synthesis_evidence"]["predictor"]["model_path"] = "../forward-model"
    path = parent / "config.json"
    path.write_text(json.dumps(raw))
    monkeypatch.chdir(parent)
    cli = workflow.load_config(path, options=Namespace(recommend=True, forward_model_path="../forward-model"))
    assert workflow.preflight_recommendation(HIT, cli)["catalog"] == config["catalog"]
    assert workflow.preflight_recommendation(HIT, raw)["provider"]["path"] == config["provider"]["path"]
    for key in ("catalog", "rules"):
        wrong = copy.deepcopy(config)
        wrong[key] = str(parent)
        with pytest.raises(ValueError, match="file"):
            workflow.preflight_recommendation(HIT, wrong)
    wrong = copy.deepcopy(config)
    wrong["synthesis_evidence"]["predictor"]["model_path"] = config["catalog"]
    with pytest.raises(ValueError, match="directory"):
        workflow.preflight_recommendation(HIT, wrong)


def test_alignment_capability_checked_without_prepare_reference(config, no_work, monkeypatch, tmp_path):
    def unavailable():
        raise ValueError("controlled missing rdShapeAlign")

    monkeypatch.setattr(conformers, "require_alignment", unavailable)
    with pytest.raises(ValueError, match="rdShapeAlign"):
        recommend(HIT, untouched_designs(), config, tmp_path / "out")


@pytest.mark.parametrize(
    "hit,extra",
    [
        ("NC(C)C(=O)O", {}),
        (HIT, {"core_smiles": "NCCN"}),
        (
            HIT,
            {
                "similarity_constraints": {
                    "landmarks": [{"id": "bad", "family": "Donor", "hit_atom_ids": [100]}]
                }
            },
        ),
    ],
)
def test_invalid_hit_core_landmarks_before_designs(config, no_work, tmp_path, hit, extra):
    config.update(extra)
    with pytest.raises(ValueError):
        recommend(hit, untouched_designs(), config, tmp_path / "out")


@pytest.mark.parametrize("args", [["--no-synllama"], ["--no-synllama", "--synthesis-config", "missing"], []])
def test_integrated_requires_enabled_config_even_when_unreadable(tmp_path, no_work, args):
    import chembounce

    assert chembounce.main(["-i", HIT, "-o", str(tmp_path / "out"), "--recommend", *args]) == 2


def test_integrated_intersection_and_alignment_fail_before_raw_generation(
    config, no_work, tmp_path, monkeypatch
):
    import chembounce

    config["property_thresholds"] = {"mw": [200, 300]}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    args = [
        "-i",
        HIT,
        "-o",
        str(tmp_path / "out"),
        "--recommend",
        "--synthesis-config",
        str(path),
        "--no-shape",
    ]
    assert chembounce.main(args + ["--mw_max", "100"]) == 2
    assert chembounce.main(args + ["--mw_min", "nan"]) == 2
    assert chembounce.main(args + ["--mw_min", "inf"]) == 2

    def unavailable():
        raise ValueError("controlled missing rdShapeAlign")

    monkeypatch.setattr(conformers, "require_alignment", unavailable)
    assert chembounce.main(args) == 2


def test_integrated_runs_exact_preflight_intersection(config, no_work, tmp_path, monkeypatch):
    import chembounce
    import pandas as pd

    config["property_thresholds"] = {"mw": [150, 300]}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    out = tmp_path / "out"
    effective = []
    preflight = workflow.preflight_recommendation

    def check(hit, cfg):
        result = preflight(hit, cfg)
        effective.append(copy.deepcopy(result))
        return result

    def hopping(*args, **kwargs):
        assert len(effective) == 1
        out.mkdir()
        (out / "run_manifest.json").write_text("{}")
        return pd.DataFrame({"Final structure": [DESIGN]}), {}

    def run(hit, designs, cfg, output):
        assert cfg == effective[0]
        assert cfg["recommendation_profile"] == profiles.DEFAULT_RECOMMENDATION_PROFILE
        assert cfg["property_thresholds"]["mw"] == [160.0, 250.0]
        assert cfg["core_smiles"] == "CC(=O)N"
        return []

    monkeypatch.setattr(workflow, "preflight_recommendation", check)
    monkeypatch.setattr(chembounce, "chembounce", hopping)
    monkeypatch.setattr(workflow, "run", run)
    assert (
        chembounce.main(
            [
                "-i",
                HIT,
                "-o",
                str(out),
                "--recommend",
                "--synthesis-config",
                str(path),
                "--mw_min",
                "160",
                "--mw_max",
                "250",
                "--core_smiles",
                "CC(=O)N",
            ]
        )
        == 3
    )


def test_lazy_public_import_and_parser_do_not_import_workflow_or_models():
    code = """
import sys
from chembounce_synthesis import recommend
from chembounce_synthesis.options import add_expansion_options
import argparse
p = argparse.ArgumentParser(); add_expansion_options(p)
assert p.parse_args(['--recommend']).recommend
assert callable(recommend)
assert not any(n in sys.modules for n in ['chembounce_synthesis.workflow','torch','transformers','rdkit'])
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_actual_api_cli_and_continuation_keep_final_three_seed_contract(
    config, strict_predictor, tmp_path, capsys
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    api_out = tmp_path / "api"
    rows = recommend(HIT, [DESIGN], config, api_out)
    assert len(rows) == 1 and rows[0]["recommendation_status"] == "qualified"
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    cli_out = tmp_path / "cli"
    assert (
        synthesis_main(
            ["--hit", HIT, "--design", DESIGN, "--recommend", "--config", str(path), "--output", str(cli_out)]
        )
        == 0
    )
    label = capsys.readouterr().out
    for value in ("hit_tanimoto>=0.35", "hit_shape>=0.7", "hit_pharmacophore>=0.5", "polar_recall>=0.75"):
        assert value in label
    cli_rows = json.loads((cli_out / "results.json").read_text())["results"]
    assert cli_rows[0]["similarity_witness_sha256"] == rows[0]["similarity_witness_sha256"]
    report = json.loads((cli_out / "refinement.json").read_text())["candidates"][0]
    assert [t["seed"] for t in report["trials"]] == [42, 43, 44]
    assert all(t["status"] == "success" and t["selective_pose_evidence"] for t in report["trials"])
    assert (cli_out / "poses_hit.sdf").read_text().strip()
    from rdkit import Chem
    import csv
    from chembounce_synthesis.report import build_report, render_html

    expected_name = profiles.DEFAULT_RECOMMENDATION_PROFILE
    provenance = json.loads((cli_out / "manifest.json").read_text())["recommendation_profile"]
    assert provenance["name"] == expected_name
    assert provenance["optional_diagnostics"]["rejected_pose_diagnostics"]["enabled"] is False
    assert provenance["artifact_sha256"] == profiles.PROFILE_SHA256[expected_name]
    assert cli_rows[0]["recommendation_profile"] == rows[0]["recommendation_profile"] == expected_name
    assert json.loads((cli_out / "summary.json").read_text())["recommendation_profile"] == {
        key: provenance[key] for key in ("name", "status", "artifact_sha256")
    }
    with (cli_out / "results.tsv").open() as handle:
        tsv = list(csv.DictReader(handle, delimiter="\t"))
    assert tsv[0]["recommendation_profile"] == expected_name
    for filename in ("poses_hit.sdf", "poses_design.sdf"):
        for mol in Chem.SDMolSupplier(str(cli_out / filename), removeHs=False):
            assert mol.GetProp("recommendation_profile") == expected_name
            assert mol.GetProp("calibration_status") == provenance["status"]
    view = build_report(cli_out)
    assert view["artifact_consistency"] == "consistent_recorded_bundle", view
    assert expected_name in render_html(view)
    continued = tmp_path / "continued"
    assert synthesis_main(["--continue-from", str(cli_out), "--recommend", "--output", str(continued)]) == 0
    assert json.loads((continued / "routes.json").read_text())["routes"][0]["carried_from"]
    # No overwrite/fallback on a populated output.
    assert synthesis_main(["--continue-from", str(cli_out), "--recommend", "--output", str(continued)]) == 2


def test_explicit_null_profile_and_missing_provider_fields_do_not_default(config, no_work):
    config["recommendation_profile"] = None
    with pytest.raises(ValueError, match="conflicts"):
        workflow.preflight_recommendation(HIT, config)
    config.pop("recommendation_profile")
    config["provider"] = {
        "kind": "synllama",
        "model_path": config["synthesis_evidence"]["predictor"]["model_path"],
    }
    with pytest.raises(ValueError, match="prompt_path"):
        workflow.preflight_recommendation(HIT, config)


@pytest.mark.parametrize(
    "historical", [profiles.PROFILE, profiles.SELECTIVE_PROFILE, profiles.GUARDED_SELECTIVE_PROFILE]
)
def test_continuation_conflict_and_replacement_are_distinct(
    config, no_work, tmp_path, monkeypatch, historical
):
    from chembounce_synthesis import cli

    inherited = copy.deepcopy(config)
    inherited["recommendation_profile"] = historical
    monkeypatch.setattr(cli, "continuation_inputs", lambda directory: (HIT, [DESIGN], inherited))
    monkeypatch.setattr(cli, "run", lambda *a, **k: pytest.fail("Inherited conflicting profile ran"))
    old = tmp_path / "old"
    args = ["--continue-from", str(old), "--recommend", "--output", str(tmp_path / "new")]
    assert synthesis_main(args + ["--recommendation-profile", profiles.DEFAULT_RECOMMENDATION_PROFILE]) == 2
    path = tmp_path / "replacement.json"
    path.write_text(json.dumps(config))

    def replaced(hit, designs, effective, out):
        assert effective["recommendation_profile"] == profiles.DEFAULT_RECOMMENDATION_PROFILE
        assert effective["carry_forward"] == [str(old)]
        return []

    monkeypatch.setattr(cli, "run", replaced)
    assert synthesis_main(args + ["--config", str(path)]) == 3
    assert inherited["recommendation_profile"] == historical


def test_standalone_preflight_precedes_design_tsv_read(config, no_work, tmp_path, monkeypatch):
    def unavailable():
        raise ValueError("alignment unavailable before TSV")

    monkeypatch.setattr(conformers, "require_alignment", unavailable)
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    assert (
        synthesis_main(
            [
                "--hit",
                HIT,
                "--designs-tsv",
                str(tmp_path / "never-read.tsv"),
                "--config",
                str(path),
                "--recommend",
                "--output",
                str(tmp_path / "out"),
            ]
        )
        == 2
    )


@pytest.mark.parametrize("criterion", ["similarity", "scaffold", "route_support"])
def test_no_fallback_when_one_required_criterion_fails(config, strict_predictor, tmp_path, criterion):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    if criterion == "similarity":
        config["minimum_scores"] = {"hit_tanimoto": 0.99}
    elif criterion == "scaffold":
        saved = Path(config["provider"]["path"])
        data = json.loads(saved.read_text())
        data["proposals"][0]["response"] = {
            "reactions": [
                {"reaction_template": "amide", "reactants": ["CC(=O)O", "Nc1ccccc1"], "product": HIT}
            ],
            "building_blocks": ["CC(=O)O", "Nc1ccccc1"],
        }
        saved.write_text(json.dumps(data))
    else:
        strict_predictor.products = ["CC", "CCC", "CO"]
    out = tmp_path / criterion
    assert recommend(HIT, [DESIGN], config, out) == []
    assert json.loads((out / "summary.json").read_text())["status"] == "no_candidates"
    assert json.loads((out / "results.json").read_text())["results"] == []


def test_integrated_shortcut_runs_actual_final_gate(config, strict_predictor, tmp_path, monkeypatch):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    import chembounce
    import pandas as pd

    out = tmp_path / "integrated"

    def known_design(*args, **kwargs):
        out.mkdir()
        (out / "run_manifest.json").write_text("{}")
        return pd.DataFrame({"Final structure": [DESIGN]}), {}

    monkeypatch.setattr(chembounce, "chembounce", known_design)
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    assert (
        chembounce.main(
            ["-i", HIT, "-o", str(out), "--recommend", "--synthesis-config", str(path), "--no-shape"]
        )
        == 0
    )
    rows = json.loads((out / "synthesis/results.json").read_text())["results"]
    assert len(rows) == 1 and rows[0]["recommendation_status"] == "qualified"
    assert rows[0]["evidence_profile"] == "computational-route-support-v1"
    assert len(strict_predictor.calls) == 1
    assert rows[0]["smiles"] == DESIGN


def test_carry_only_saved_archive_still_replays_under_new_entrypoint(config, strict_predictor, tmp_path):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    first = tmp_path / "first"
    assert recommend(HIT, [DESIGN], config, first)
    empty = tmp_path / "empty.json"
    empty.write_text(json.dumps({"schema_version": 1, "proposals": []}))
    config["provider"] = {"kind": "saved", "path": str(empty)}
    path = tmp_path / "replacement.json"
    path.write_text(json.dumps(config))
    out = tmp_path / "carried"
    assert (
        synthesis_main(
            ["--continue-from", str(first), "--config", str(path), "--recommend", "--output", str(out)]
        )
        == 0
    )
    manifest = json.loads((out / "manifest.json").read_text())
    assert "without new provider calls" in manifest["carried_only_policy"]
    assert json.loads((out / "summary.json").read_text())["scored_products"] == 1


def test_conflicting_cli_profile_without_config_returns_error_not_traceback(tmp_path, no_work):
    import chembounce

    assert (
        chembounce.main(
            [
                "-i",
                HIT,
                "-o",
                str(tmp_path / "out"),
                "--recommend",
                "--recommendation-profile",
                profiles.PROFILE,
            ]
        )
        == 2
    )


def test_minimum_runtime_fails_closed_without_consuming_designs(config, tmp_path):
    if conformers.rdShapeAlign is not None:
        pytest.skip("Minimum-runtime capability assertion")
    with pytest.raises(ValueError, match="rdShapeAlign"):
        recommend(HIT, untouched_designs(), config, tmp_path / "out")


def test_packaged_asset_template_uses_same_policy_after_explicit_local_paths(config, no_work):
    template = Path(__file__).resolve().parents[2] / "examples/synthesis/strict-recommendation-config.json"
    raw = json.loads(template.read_text())
    for key in ("catalog", "rules", "provider", "synthesis_evidence"):
        raw[key] = config[key]
    effective = workflow.preflight_recommendation(HIT, raw)
    assert effective["recommendation_profile"] == profiles.DEFAULT_RECOMMENDATION_PROFILE
    assert effective["minimum_scores"] == {"hit_tanimoto": 0.35, "hit_shape": 0.7, "hit_pharmacophore": 0.5}
    assert effective["refinement"]["seeds"] == [42, 43, 44]


@pytest.mark.parametrize("explicit_v3", [False, True])
def test_historical_continuation_without_shortcut_preserves_or_explicitly_migrates(
    config, no_work, monkeypatch, tmp_path, explicit_v3
):
    from chembounce_synthesis import cli

    inherited = workflow._config(config | {"recommendation_profile": profiles.SELECTIVE_PROFILE})
    inherited["carry_forward"] = [str(tmp_path / "old")]
    before = copy.deepcopy(inherited)
    monkeypatch.setattr(cli, "continuation_inputs", lambda directory: (HIT, [DESIGN], inherited))
    received = []

    def run(hit, designs, effective, output):
        received.append(effective)
        return []

    monkeypatch.setattr(cli, "run", run)
    args = ["--continue-from", str(tmp_path / "old"), "--output", str(tmp_path / "new")]
    if explicit_v3:
        args += ["--recommendation-profile", profiles.DEFAULT_RECOMMENDATION_PROFILE]
    assert synthesis_main(args) == 3
    assert inherited == before
    expected = profiles.DEFAULT_RECOMMENDATION_PROFILE if explicit_v3 else profiles.SELECTIVE_PROFILE
    assert received[0]["recommendation_profile"] == expected
    for key in ("carry_forward", "max_design_attempts", "max_scored_products", "refinement"):
        assert received[0][key] == before[key]
