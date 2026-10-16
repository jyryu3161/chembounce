"""Actual workflow boundary with small known chemistry and controlled corruptions.

No model inference or empirical calibration. Legacy route replay is explicitly
used in engineering fixtures; the strict evidence fixture uses controlled beams.
"""

import csv
import json

import pytest
from rdkit import Chem

from chembounce_synthesis import selective_scoring, workflow
from chembounce_synthesis.cli import main
from chembounce_synthesis.common import canonical
from chembounce_synthesis.conformers import ConformerSettings
from chembounce_synthesis.continuation import continuation_inputs
from chembounce_synthesis.similarity_constraints import BoundSimilarity, MINIMUM_SCORES, VERSION

HIT, DESIGN = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"
ALTERNATE = "CCC(=O)Nc1ccncc1"


@pytest.fixture
def alignment():
    pytest.importorskip("rdkit.Chem.rdShapeAlign")


def config_for(tmp_path, assets, proposal, designs=(DESIGN,)):
    proposals = tmp_path / "proposals.json"
    proposals.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": s, "response": proposal} for s in designs]})
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(proposals)},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "conformers": {"num_conformers": 8, "scoring_objective": "original-hit-v1", "seed": 42},
        "refinement": {
            "top_n": 2,
            "num_conformers": 8,
            "seeds": [42, 43, 44],
            "max_total_conformers": 48,
            "ranking_statistic": "minimum",
        },
        "recommendation": {"mode": "qualified"},
        "minimum_scores": dict(MINIMUM_SCORES),
        "similarity_constraints": {"version": VERSION},
        "max_scored_products": 2,
        "max_results": 2,
    }


def read(out, name):
    return json.loads((out / name).read_text())


def assert_empty(out):
    assert read(out, "summary.json")["status"] == "no_candidates"
    assert read(out, "results.json")["results"] == []
    assert len((out / "results.tsv").read_text().splitlines()) == 1
    for role in ("hit", "design"):
        assert not (out / f"poses_{role}.sdf").read_text().strip()
        assert not (out / f"references_{role}.sdf").read_text().strip()
    assert not any(d["qualified"] for d in read(out, "qualification.json")["decisions"])


def assert_seed_geometry(out, config):
    reports = read(out, "refinement.json")["candidates"]
    report = next(r for r in reports if r["status"] == "success")
    bound = selective_scoring.configured_constraints(workflow._config(config), HIT)
    for trial in report["trials"]:
        assert trial["status"] == "success"
        restored = selective_scoring.restore_result(trial)
        selective_scoring.validate_result(
            restored,
            DESIGN,
            ConformerSettings(**trial["settings"]),
            hit=HIT,
            design=report["design"],
            similarity=bound,
        )
        assert trial["similarity_constraints"]["witness"]["reference_pose"]["coordinates"]
    assert [t["seed"] for t in report["trials"]] == [42, 43, 44]
    return report


def test_actual_cli_selective_contract_and_all_exported_seed_bindings(
    tmp_path, assets, proposal, alignment, capsys
):
    config = config_for(tmp_path, assets, proposal)
    config["core_smiles"] = "CC(=O)N"
    config["similarity_constraints"]["landmarks"] = [
        {"id": "amide-oxygen", "family": "Acceptor", "hit_atom_ids": [2], "max_distance": 1.0}
    ]
    path, out = tmp_path / "config.json", tmp_path / "out"
    path.write_text(json.dumps(config))
    assert main(["--hit", HIT, "--design", DESIGN, "--config", str(path), "--output", str(out)]) == 0
    assert "explicit engineering constraints" in capsys.readouterr().out
    row = read(out, "results.json")["results"][0]
    assert row["recommendation_status"] == "qualified"
    report = assert_seed_geometry(out, config)
    assert row["pose_seed"] == report["pose_seed"]
    assert (
        row["similarity_constraint_sha256"]
        == read(out, "manifest.json")["recommendation_policy"]["similarity_constraints"]["sha256"]
    )
    with (out / "results.tsv").open() as stream:
        table = list(csv.DictReader(stream, delimiter="\t"))[0]
    keys = (
        "pose_selection_policy",
        "similarity_constraint_policy",
        "similarity_constraint_sha256",
        "similarity_witness_sha256",
        "recommendation_profile",
        "calibration_status",
    )
    for key in keys:
        assert table[key] == row[key]
    for file in ("poses_hit.sdf", "poses_design.sdf", "references_hit.sdf", "references_design.sdf"):
        mol = list(Chem.SDMolSupplier(str(out / file)))[0]
        assert mol.GetProp("_Name") == row["route_id"]
        assert mol.GetIntProp("pose_seed") == row["pose_seed"]
        for key in keys:
            assert mol.GetProp(key) == row[key]
    selected_trial = next(t for t in report["trials"] if t["seed"] == row["pose_seed"])
    assert (
        row["similarity_witness_sha256"]
        == selective_scoring.result_fields(selected_trial)["similarity_witness_sha256"]
    )


def test_integrated_cli_uses_selective_final_workflow(tmp_path, assets, proposal, alignment):
    import chembounce

    config = config_for(tmp_path, assets, proposal, (DESIGN, "CC(=O)Nc1cccnc1", "CC(=O)Nc1ccccn1"))
    path, replacements = tmp_path / "integrated.json", tmp_path / "replacements.smi"
    path.write_text(json.dumps(config))
    replacements.write_text("c1ccncc1\n")
    out = tmp_path / "integrated"
    assert (
        chembounce.main(
            [
                "-i",
                HIT,
                "-o",
                str(out),
                "--fragments",
                "c1ccccc1",
                "--replace_scaffold_files",
                str(replacements),
                "--overall_max_n",
                "1",
                "--cand_max_n__rplc",
                "1",
                "--no-shape",
                "-t",
                "0.1",
                "--search-threshold",
                "0.1",
                "--synthesis-config",
                str(path),
            ]
        )
        == 0
    )
    rows = read(out / "synthesis", "results.json")["results"]
    assert len(rows) == 1
    assert rows[0]["pose_selection_policy"] == selective_scoring.VERSION
    assert rows[0]["recommendation_status"] == "qualified"
    assert_seed_geometry(out / "synthesis", config)


@pytest.mark.parametrize("mutation", ["missing_geometry", "coordinate", "witness", "selector"])
def test_final_boundary_rejects_nonrepresentative_seed_corruption(
    tmp_path, assets, proposal, alignment, monkeypatch, mutation
):
    config = config_for(tmp_path, assets, proposal)
    original = workflow.refine_candidate

    def corrupted(*args, **kwargs):
        result, report = original(*args, **kwargs)
        assert result is not None
        trial = next(t for t in report["trials"] if t["seed"] != report["pose_seed"])
        if mutation == "missing_geometry":
            del trial["selective_pose_evidence"]
        elif mutation == "coordinate":
            trial["selective_pose_evidence"]["hit_candidate"]["coordinates"][0][0] += 0.5
        elif mutation == "witness":
            trial["similarity_constraints"]["witness"]["polar"]["recall"] = 1.0
        else:
            trial["pose_selection_policy"] = "legacy"
        return result, report

    monkeypatch.setattr(workflow, "refine_candidate", corrupted)
    out = tmp_path / "rejected"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    assert_empty(out)
    assert any(
        "Selective final contract" in reason
        for d in read(out, "qualification.json")["decisions"]
        for reason in d["reasons"]
    )


def test_failed_required_seed_does_not_hide_behind_good_seeds_and_audit_survives(
    tmp_path, assets, proposal, alignment, monkeypatch
):
    original = selective_scoring.score_candidate

    def failure(smiles, hit, design, settings, **kwargs):
        result = original(smiles, hit, design, settings, **kwargs)
        if settings.seed == 43:
            raise selective_scoring.SimilarityAbstention(
                "controlled_required_seed_failure", result[0]["pose_selection"]
            )
        return result

    monkeypatch.setattr(selective_scoring, "score_candidate", failure)
    out = tmp_path / "failed"
    assert workflow.run(HIT, [DESIGN], config_for(tmp_path, assets, proposal), out) == []
    assert_empty(out)
    report = read(out, "refinement.json")["candidates"][0]
    assert [t["status"] for t in report["trials"]] == ["success", "failed", "success"]
    failure = report["trials"][1]["similarity_failure"]
    assert failure["ensemble_request_status"] == "requested"
    assert failure["scorer_audit"]["hit_pairs_attempted"] > 0


@pytest.mark.parametrize("limit", [24, 48])
def test_failed_primary_design_refills_from_valid_alternate_with_same_caps(
    tmp_path, assets, proposal, alignment, monkeypatch, limit
):
    original = selective_scoring.score_candidate

    def failure(smiles, hit, design, settings, **kwargs):
        result = original(smiles, hit, design, settings, **kwargs)
        if design == DESIGN and settings.seed == 43:
            raise selective_scoring.SimilarityAbstention(
                "controlled_primary_failure", result[0]["pose_selection"]
            )
        return result

    monkeypatch.setattr(selective_scoring, "score_candidate", failure)
    config = config_for(tmp_path, assets, proposal, (DESIGN, ALTERNATE))
    config["refinement"]["max_total_conformers"] = limit
    out = tmp_path / "fallback"
    rows = workflow.run(HIT, [DESIGN, ALTERNATE], config, out)
    if limit == 48:
        assert len(rows) == 1 and rows[0]["design"] == ALTERNATE
        assert_seed_geometry(out, config)
    else:
        assert rows == []
        assert_empty(out)
    assert read(out, "summary.json")["refinement_budget"]["candidate_conformer_requests"] == limit
    assert read(out, "summary.json")["scored_products"] == 2


def test_continuation_recomputes_changed_constraints_and_ignores_old_scores(
    tmp_path, assets, proposal, alignment
):
    config = config_for(tmp_path, assets, proposal)
    first = tmp_path / "first"
    assert workflow.run(HIT, [DESIGN], config, first)
    report = assert_seed_geometry(first, config)
    assert all(t["similarity_constraints"]["witness"]["polar"]["recall"] < 1 for t in report["trials"])
    # These files are not imported as scorer evidence when continuing.
    (first / "scores.json").write_text('{"pretend_all_constraints_passed":true}')
    (first / "refinement.json").write_text('{"pretend_all_seeds_passed":true}')
    hit, designs, resumed = continuation_inputs(first)
    resumed["similarity_constraints"]["polar_recall"] = 1.0
    out = tmp_path / "continued"
    assert workflow.run(hit, designs, resumed, out) == []
    assert_empty(out)
    assert (
        read(out, "manifest.json")["recommendation_policy"]["similarity_constraints"]["sha256"]
        != read(first, "manifest.json")["recommendation_policy"]["similarity_constraints"]["sha256"]
    )
    assert any(
        r.get("similarity_failure", {}).get("scorer_audit", {}).get("hit_pairs_attempted", 0) > 0
        for r in read(out, "routes.json")["routes"]
    )


@pytest.mark.parametrize("phase", ["preflight", "validation"])
def test_screening_failure_accounting_distinguishes_preflight_from_post_score(
    tmp_path, assets, proposal, alignment, monkeypatch, phase
):
    if phase == "preflight":
        original = BoundSimilarity.preflight_candidate

        def fail(self, smiles):
            if canonical(smiles) == DESIGN:
                raise ValueError("controlled candidate preflight rejection")
            return original(self, smiles)

        monkeypatch.setattr(BoundSimilarity, "preflight_candidate", fail)
        monkeypatch.setattr(
            selective_scoring,
            "score_candidate",
            lambda *a, **k: pytest.fail("scorer called after failed preflight"),
        )
    else:
        original = selective_scoring.validate_result

        def fail(result, *args, **kwargs):
            original(result, *args, **kwargs)
            raise ValueError("controlled post-score validation rejection")

        monkeypatch.setattr(selective_scoring, "validate_result", fail)
    out = tmp_path / phase
    assert workflow.run(HIT, [DESIGN], config_for(tmp_path, assets, proposal), out) == []
    failure = read(out, "routes.json")["routes"][0]["similarity_failure"]
    assert failure["phase"] == phase
    assert failure["ensemble_request_status"] == ("none_requested" if phase == "preflight" else "requested")
    assert (failure["scorer_audit"] is None) == (phase == "preflight")


@pytest.mark.parametrize(
    "change",
    [
        {"reference_sdf": "never-read.sdf"},
        {"pocket": {}},
        {"conformers": {"num_conformers": 12}},
        {"conformers": {"reject_unspecified_stereo": False}},
        {"refinement": {"adaptive_max_conformers": 16}},
        {"refinement": {"top_n": 0}},
        {"refinement": {"threshold_policy": "representative"}},
        {"refinement": {"seeds": [42, 43]}},
        {"similarity_constraints": {"version": "unknown"}},
    ],
)
def test_unsupported_modes_fail_before_any_assets(tmp_path, assets, proposal, monkeypatch, change):
    config = config_for(tmp_path, assets, proposal)
    for key, value in change.items():
        if isinstance(value, dict):
            config.setdefault(key, {}).update(value)
        else:
            config[key] = value
    for name in ("Catalog", "Rules", "_provider", "_evidence_evaluator"):
        monkeypatch.setattr(workflow, name, lambda *a, **k: pytest.fail("assets loaded"))
    with pytest.raises(ValueError):
        workflow.run(HIT, [DESIGN], config, tmp_path / "bad")


@pytest.mark.parametrize(
    "change",
    [
        {"core_smiles": "S"},
        {
            "similarity_constraints": {
                "version": VERSION,
                "landmarks": [{"id": "missing", "family": "Acceptor", "hit_atom_ids": [999]}],
            }
        },
    ],
)
def test_query_selectors_core_preflight_in_api_and_integrated_cli(
    tmp_path, assets, proposal, monkeypatch, change
):
    import chembounce

    config = config_for(tmp_path, assets, proposal)
    config.update(change)
    for name in ("Catalog", "Rules", "_provider", "_evidence_evaluator"):
        monkeypatch.setattr(workflow, name, lambda *a, **k: pytest.fail("assets loaded"))
    monkeypatch.setattr(chembounce, "chembounce", lambda *a, **k: pytest.fail("upstream generation started"))
    with pytest.raises(ValueError):
        workflow.run(HIT, [DESIGN], config, tmp_path / "bad-api")
    assert (
        read(tmp_path / "bad-api", "summary.json")["similarity_preflight"]["ensemble_request_status"]
        == "none_requested"
    )
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    assert chembounce.main(["-i", HIT, "-o", str(tmp_path / "bad-cli"), "--synthesis-config", str(path)]) == 2


def test_v1_named_profile_with_changed_selector_never_claims_calibration_match(tmp_path, assets, proposal):
    from chembounce_synthesis.profiles import profile_provenance, load_profile

    config = config_for(tmp_path, assets, proposal)
    config["recommendation_profile"] = "qualified-development-v1"
    config["minimum_scores"] = load_profile("qualified-development-v1")["minimum_scores"]
    config["synthesis_evidence"] = {"predictor": {"model_path": "/nonexistent/no-load"}}
    resolved = workflow._config(config)
    provenance = profile_provenance(resolved)
    assert provenance["status"] == "custom_unvalidated"
    assert "similarity_constraints" in provenance["material_overrides"]


def test_design_feature_preflight_failure_does_not_load_assets(
    tmp_path, assets, proposal, alignment, monkeypatch
):
    original = selective_scoring.feature_records

    def bounded(mol, *args):
        if canonical(Chem.MolToSmiles(mol)) == DESIGN:
            raise ValueError("controlled design feature bound rejection")
        return original(mol, *args)

    monkeypatch.setattr(selective_scoring, "feature_records", bounded)
    for name in ("Catalog", "Rules", "_provider", "_evidence_evaluator"):
        monkeypatch.setattr(
            workflow, name, lambda *a, **k: pytest.fail("assets loaded after all designs failed")
        )
    out = tmp_path / "no-eligible-designs"
    assert workflow.run(HIT, [DESIGN], config_for(tmp_path, assets, proposal), out) == []
    summary = read(out, "summary.json")
    assert summary["scored_products"] == 0
    assert summary["assets_not_loaded"] == ["catalog", "rules", "provider", "synthesis_evidence"]
    assert summary["designs_preflight_rejected"] == 1


@pytest.mark.parametrize("matches", [True, False])
def test_strict_route_evidence_still_required_no_alternate_borrowing(
    tmp_path, assets, proposal, alignment, strict_predictor, matches
):
    config = config_for(tmp_path, assets, proposal)
    config["synthesis_evidence"] = {
        "profile": "computational-route-support-v1",
        "predictor": {"kind": "reactiont5", "model_path": "/controlled-no-model"},
        "max_predictor_calls": 1,
    }
    strict_predictor.products = None if matches else ["CCO", "CC", "CO"]
    out = tmp_path / "strict"
    rows = workflow.run(HIT, [DESIGN], config, out)
    if matches:
        assert len(rows) == 1
        assert rows[0]["evidence_profile"] == "computational-route-support-v1"
        assert rows[0]["route_support_binding"]
        assert_seed_geometry(out, config)
    else:
        assert rows == []
        assert_empty(out)
        assert read(out, "summary.json")["scored_products"] == 0
    assert len(strict_predictor.calls) == 1
