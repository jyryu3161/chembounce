"""Exercise the final contract through real route replay and workflow exports.

Tiny positive thresholds here are engineering fixtures, not empirical cutoffs.
Scoring/evidence perturbations deliberately simulate incomplete upstream work.
"""

import csv
import json

import pytest
from rdkit import Chem

pytest.importorskip("rdkit.Chem.rdShapeAlign")

from chembounce_synthesis import refinement as refinement_module, workflow
from chembounce_synthesis.cli import main
from chembounce_synthesis.continuation import continuation_inputs
from chembounce_synthesis.evidence import SynthesisEvidenceEvaluator

HIT, DESIGN = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"


def configuration(tmp_path, assets, proposal, designs=(DESIGN,)):
    path = tmp_path / "proposals.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "proposals": [{"design": design, "response": proposal} for design in designs],
            }
        )
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "conformers": {"num_conformers": 2},
        "recommendation": {"mode": "qualified", "policy_version": "qualified-v1"},
        "minimum_scores": {"hit_tanimoto": 0.01, "hit_shape": 0.01, "hit_pharmacophore": 0.01},
    }


def read(out, name):
    return json.loads((out / name).read_text())


def assert_empty(out):
    assert read(out, "summary.json")["status"] == "no_candidates"
    assert read(out, "results.json")["results"] == []
    assert len((out / "results.tsv").read_text().splitlines()) == 1
    assert not (out / "poses_hit.sdf").read_text().strip()
    assert all(not d["qualified"] for d in read(out, "qualification.json")["decisions"])


def test_real_qualified_route_has_same_selection_and_policy_in_every_artifact(tmp_path, assets, proposal):
    config = configuration(tmp_path, assets, proposal)
    config["core_smiles"] = "CC(=O)N"
    out = tmp_path / "qualified"
    rows = workflow.run(HIT, [DESIGN], config, out)
    assert len(rows) == 1 and rows[0]["recommendation_status"] == "qualified"
    assert read(out, "results.json")["results"] == rows
    assert read(out, "manifest.json")["recommendation_policy"] == read(out, "qualification.json")["policy"]
    with (out / "results.tsv").open() as stream:
        table = list(csv.DictReader(stream, delimiter="\t"))
    assert table[0]["route_id"] == rows[0]["route_id"]
    assert table[0]["recommendation_status"] == "qualified"
    for role in ("hit", "design"):
        poses = list(Chem.SDMolSupplier(str(out / f"poses_{role}.sdf")))
        assert len(poses) == 1 and poses[0].GetProp("_Name") == rows[0]["route_id"]
        assert poses[0].GetProp("recommendation_status") == "qualified"
        assert Chem.MolToSmiles(poses[0]) == rows[0]["smiles"]
    decision = read(out, "qualification.json")["decisions"][0]
    assert decision["qualified"] and decision["reasons"] == []
    assert read(out, "summary.json")["qualified_products"] == 1
    assert rows[0]["convergence_status"] == "not_assessed"  # no stability claim


@pytest.mark.parametrize(
    "minimum",
    [
        {},
        {"hit_tanimoto": 0},
        {"hit_shape": float("nan")},
        {"hit_pharmacophore": float("inf")},
        {"hit_shape": True},
    ],
)
def test_missing_zero_nonfinite_minima_are_configuration_errors(tmp_path, assets, proposal, minimum):
    config = configuration(tmp_path, assets, proposal)
    if minimum:
        config["minimum_scores"].update(minimum)
    else:
        config["minimum_scores"] = {}
    with pytest.raises(ValueError):
        workflow.run(HIT, [DESIGN], config, tmp_path / "bad")
    assert read(tmp_path / "bad", "summary.json")["status"] == "error"


@pytest.mark.parametrize(
    "patch",
    [
        {"mode": "expand"},
        {"conformers": {"reject_unspecified_stereo": False}},
        {"refinement": {"top_n": 1, "threshold_policy": "representative"}},
        {"recommendation": {"mode": "qualified", "policy_version": "invented"}},
    ],
)
def test_incompatible_qualified_configuration_fails_early(tmp_path, assets, proposal, patch):
    config = configuration(tmp_path, assets, proposal)
    config.update(patch)
    with pytest.raises(ValueError):
        workflow.run(HIT, [DESIGN], config, tmp_path / "bad")


def test_2d_cannot_be_qualified(tmp_path, assets, proposal):
    with pytest.raises(ValueError, match="completed 3D"):
        workflow.run(
            HIT, [DESIGN], configuration(tmp_path, assets, proposal), tmp_path / "2d", scoring_mode="2d"
        )


def test_similar_design_without_route_cannot_recommend(tmp_path, assets, proposal):
    config = configuration(tmp_path, assets, {"reactions": [], "building_blocks": []})
    out = tmp_path / "no-route"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    assert_empty(out)
    assert read(out, "summary.json")["failures"]


def test_design_identity_cannot_compensate_for_low_original_hit_similarity(
    tmp_path, assets, proposal, monkeypatch
):
    config = configuration(tmp_path, assets, proposal)
    config["minimum_scores"]["hit_tanimoto"] = 0.99
    config["route_search_policy"] = (
        "design-v1"  # Exercise the unchanged final gate, not early search rejection.
    )
    monkeypatch.setattr(workflow, "score_candidate", lambda *a, **k: pytest.fail("2D cutoff must precede 3D"))
    out = tmp_path / "distant"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    route = read(out, "routes.json")["routes"][0]
    assert route["tanimoto"]["design_tanimoto"] == 1
    assert route["tanimoto"]["hit_tanimoto"] < 0.99
    assert_empty(out)


def test_state_only_change_is_not_a_qualified_hop(tmp_path, assets, proposal):
    out = tmp_path / "state"
    config = configuration(tmp_path, assets, proposal)
    config["route_search_policy"] = "design-v1"
    assert workflow.run("CC(=O)Nc1cc[nH+]cc1", [DESIGN], config, out) == []
    assert read(out, "routes.json")["routes"][0]["selection_status"] == "constraint_rejected"
    assert_empty(out)


@pytest.mark.parametrize(
    "defect",
    [
        "nan",
        "missing_metric",
        "missing_reference",
        "missing_pose",
        "wrong_pose",
        "wrong_hit",
        "low_shape",
        "low_pharmacophore",
    ],
)
def test_invalid_final_scoring_never_exports(tmp_path, assets, proposal, monkeypatch, defect):
    real = workflow.score_candidate

    def score(*args, **kwargs):
        info, hit_pose, design_pose = real(*args, **kwargs)
        if defect == "nan":
            info["hit"]["shape"] = float("nan")
        elif defect == "missing_metric":
            del info["hit"]["pharmacophore"]
        elif defect == "missing_reference":
            del info["hit"]["reference_conformer_id"]
        elif defect == "missing_pose":
            hit_pose = None
        elif defect == "wrong_pose":
            hit_pose = Chem.MolFromSmiles(HIT)
        elif defect == "wrong_hit":
            info["hit_smiles"] = DESIGN
        else:
            info["hit"][defect.removeprefix("low_")] = 0.0
        return info, hit_pose, design_pose

    monkeypatch.setattr(workflow, "score_candidate", score)
    out = tmp_path / defect
    assert workflow.run(HIT, [DESIGN], configuration(tmp_path, assets, proposal), out) == []
    assert_empty(out)


@pytest.mark.parametrize(
    "defect",
    ["missing_seed", "wrong_seed", "failed_seed", "low_seed", "wrong_ensemble", "unassessed_convergence"],
)
def test_final_boundary_rechecks_required_seed_evidence(tmp_path, assets, proposal, monkeypatch, defect):
    real = workflow.refine_candidate

    def corrupt(*args, **kwargs):
        result, report = real(*args, **kwargs)
        assert result is not None
        if defect == "missing_seed":
            report["trials"].pop()
        elif defect == "wrong_seed":
            report["trials"][0]["settings"]["seed"] = 1234
        elif defect == "failed_seed":
            report["trials"][0]["status"] = "failed"
        elif defect == "low_seed":
            report["trials"][0]["hit"]["shape"] = 0.0
        elif defect == "wrong_ensemble":
            report["num_conformers"] = 12345
        else:
            report["convergence_status"] = "not_assessed"
        return result, report

    config = configuration(tmp_path, assets, proposal)
    config["refinement"] = {"top_n": 1, "num_conformers": 2, "seeds": [42, 43, 44], "max_score_range": 1.0}
    monkeypatch.setattr(workflow, "refine_candidate", corrupt)
    out = tmp_path / defect
    assert workflow.run(HIT, [DESIGN], config, out) == []
    assert_empty(out)
    assert read(out, "qualification.json")["decisions"][0]["status"] == "rejected"


def test_successful_best_seed_cannot_hide_failed_required_seed(tmp_path, assets, proposal, monkeypatch):
    real = refinement_module.score_candidate

    def score(smiles, hit, design, settings, *args, **kwargs):
        if settings.seed == 43:
            raise ValueError("controlled failed required seed")
        return real(smiles, hit, design, settings, *args, **kwargs)

    monkeypatch.setattr(refinement_module, "score_candidate", score)
    config = configuration(tmp_path, assets, proposal)
    config["refinement"] = {"top_n": 1, "num_conformers": 2, "seeds": [42, 43, 44]}
    out = tmp_path / "failed-seed"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    assert_empty(out)


def test_complete_refined_protocol_qualifies_and_budget_shortfall_does_not(tmp_path, assets, proposal):
    config = configuration(tmp_path, assets, proposal)
    config["refinement"] = {"top_n": 1, "num_conformers": 2, "seeds": [42, 43, 44], "max_score_range": 1.0}
    out = tmp_path / "complete"
    rows = workflow.run(HIT, [DESIGN], config, out)
    assert len(rows) == 1 and rows[0]["scoring_stage"] == "refined"
    assert rows[0]["seed_successes"] == rows[0]["seed_count"] == 3
    assert rows[0]["convergence_status"] == "converged"
    config["refinement"]["max_total_conformers"] = 5
    out = tmp_path / "budget"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    assert_empty(out)
    assert "refinement.max_total_conformers" in read(out, "summary.json")["limits_reached"]


def test_required_anchor_cannot_disappear_from_final_product(tmp_path, assets, proposal):
    config = configuration(tmp_path, assets, proposal)
    config["core_smiles"] = "c1ccccc1"
    config["route_search_policy"] = "design-v1"
    out = tmp_path / "anchor"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    assert "core/anchor" in read(out, "routes.json")["routes"][0]["reason"]
    assert_empty(out)


@pytest.mark.parametrize("defect", ["product", "policy", "step", "claimed_pass"])
def test_claimed_evidence_pass_cannot_bypass_exact_requested_evidence(
    tmp_path, assets, proposal, monkeypatch, defect
):
    config = configuration(tmp_path, assets, proposal)
    if defect == "claimed_pass":
        config["synthesis_evidence"] = {"policy": {"require_forward_agreement": True}}
    evaluator = SynthesisEvidenceEvaluator(policy=config.get("synthesis_evidence", {}).get("policy"))
    real = evaluator.evaluate

    def evaluate(*args, **kwargs):
        evidence = real(*args, **kwargs)
        if defect == "product":
            evidence["product"] = HIT
        elif defect == "policy":
            evidence["policy"] = {}
        elif defect == "step":
            evidence["steps"][0]["product"] = HIT
        else:
            evidence["eligibility"]["passed"] = True
            evidence["steps"][0]["eligibility"]["passed"] = True
        return evidence

    monkeypatch.setattr(evaluator, "evaluate", evaluate)
    monkeypatch.setattr(workflow, "_evidence_evaluator", lambda config: evaluator)
    out = tmp_path / defect
    assert workflow.run(HIT, [DESIGN], config, out) == []
    assert_empty(out)


def test_alternate_route_survives_a_failed_final_evidence_gate(tmp_path, assets, proposal, monkeypatch):
    config = configuration(tmp_path, assets, proposal)
    path = tmp_path / "proposals.json"
    data = json.loads(path.read_text())
    data["proposals"] *= 2
    path.write_text(json.dumps(data))
    evaluator = SynthesisEvidenceEvaluator()
    real = evaluator.evaluate
    calls = []

    def evaluate(*args, **kwargs):
        evidence = real(*args, **kwargs)
        calls.append(1)
        if len(calls) == 1:
            evidence["product"] = HIT
        return evidence

    monkeypatch.setattr(evaluator, "evaluate", evaluate)
    monkeypatch.setattr(workflow, "_evidence_evaluator", lambda config: evaluator)
    out = tmp_path / "alternate-route"
    rows = workflow.run(HIT, [DESIGN], config, out)
    assert rows[0]["route_id"] == "route_000001"
    assert rows[0]["alternative_route_ids"] == []
    assert [d["qualified"] for d in read(out, "qualification.json")["decisions"]] == [False, True]


@pytest.mark.parametrize("defect", [None, "missing_predictions", "wrong_matching_policy"])
def test_required_forward_evidence_is_bound_to_actual_prediction(
    tmp_path, assets, proposal, monkeypatch, defect
):
    class Predictor:
        provenance = {"provider": "controlled_test_predictor", "uses_conditions": False}

        def token_budget(self, reactants, reagents, conditions):
            return {"input_tokens": 8, "max_output_tokens": 24}

        def predict(self, reactants, reagents, conditions):
            return {"products": [DESIGN], "complete": [True], "input_tokens": 8, "output_tokens": 12}

    config = configuration(tmp_path, assets, proposal)
    config["synthesis_evidence"] = {"policy": {"require_forward_agreement": True}}
    evaluator = SynthesisEvidenceEvaluator(Predictor(), config["synthesis_evidence"]["policy"])
    real = evaluator.evaluate

    def evaluate(*args, **kwargs):
        evidence = real(*args, **kwargs)
        prediction = evidence["steps"][0]["forward_prediction"]
        if defect == "missing_predictions":
            prediction.pop("predictions")
        elif defect == "wrong_matching_policy":
            prediction["matching_policy"] = "component"
        return evidence

    monkeypatch.setattr(evaluator, "evaluate", evaluate)
    monkeypatch.setattr(workflow, "_evidence_evaluator", lambda config: evaluator)
    out = tmp_path / "forward"
    rows = workflow.run(HIT, [DESIGN], config, out)
    if defect:
        assert rows == []
        assert_empty(out)
    else:
        assert len(rows) == 1 and rows[0]["recommendation_status"] == "qualified"


def test_alternate_design_survives_refined_qualification_before_family_selection(
    tmp_path, assets, proposal, monkeypatch
):
    other = "CCC(=O)Nc1ccncc1"
    config = configuration(tmp_path, assets, proposal, designs=(DESIGN, other))
    config["refinement"] = {"top_n": 1, "num_conformers": 2, "seeds": [42, 43]}
    real = workflow.refine_candidate
    seen = []

    def refine(*args, **kwargs):
        result, report = real(*args, **kwargs)
        seen.append(args[2])
        if len(seen) == 1:
            report["seed_successes"] = 1
        return result, report

    monkeypatch.setattr(workflow, "refine_candidate", refine)
    out = tmp_path / "alternate-design"
    rows = workflow.run(HIT, [DESIGN, other], config, out)
    assert len(seen) == 2 and len(rows) == 1 and rows[0]["design"] == seen[1]
    assert rows[0]["alternative_route_ids"] == []
    assert sum(d["qualified"] for d in read(out, "qualification.json")["decisions"]) == 1


def test_continuation_requalifies_and_cli_distinguishes_exploration(tmp_path, assets, proposal, capsys):
    config = configuration(tmp_path, assets, proposal)
    config.pop("recommendation")
    config.pop("minimum_scores")
    old = tmp_path / "old"
    rows = workflow.run(HIT, [DESIGN], config, old)
    assert rows[0]["recommendation_status"] == "exploratory_verified"
    assert not read(old, "qualification.json")["decisions"][0]["qualified"]
    assert (
        main(
            [
                "--continue-from",
                str(old),
                "--output",
                str(tmp_path / "bad"),
                "--recommendation-mode",
                "qualified",
            ]
        )
        == 2
    )
    hit, designs, resumed = continuation_inputs(old)
    resumed.update(
        recommendation={"mode": "qualified"},
        minimum_scores={"hit_tanimoto": 0.99, "hit_shape": 0.01, "hit_pharmacophore": 0.01},
    )
    out = tmp_path / "requalified"
    assert workflow.run(hit, designs, resumed, out) == []
    assert_empty(out)
    assert (
        main(
            [
                "--continue-from",
                str(old),
                "--output",
                str(tmp_path / "exploratory"),
                "--recommendation-mode",
                "exploratory",
            ]
        )
        == 0
    )
    assert "not qualified recommendations" in capsys.readouterr().out


def test_standalone_and_integrated_cli_reach_qualified_contract(tmp_path, assets, proposal, capsys):
    import chembounce

    config = configuration(tmp_path, assets, proposal)
    # Real hopping can choose any positional aminopyridine.
    path = tmp_path / "proposals.json"
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
    path.write_text(json.dumps({"schema_version": 1, "proposals": proposals}))
    config.pop("recommendation")  # exercise shared CLI override
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    assert (
        main(
            [
                "--hit",
                HIT,
                "--design",
                DESIGN,
                "--config",
                str(config_path),
                "--output",
                str(tmp_path / "standalone"),
                "--recommendation-mode",
                "qualified",
            ]
        )
        == 0
    )
    replacements = tmp_path / "replacements.smi"
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
                str(config_path),
                "--recommendation-mode",
                "qualified",
            ]
        )
        == 0
    )
    assert read(out / "synthesis", "results.json")["results"][0]["recommendation_status"] == "qualified"
    assert capsys.readouterr().out.count("Qualified recommendations [evidence=legacy-conditional-v1]:") == 2
    assert (
        chembounce.main(
            [
                "-i",
                HIT,
                "-o",
                str(tmp_path / "disabled"),
                "--synthesis-config",
                str(config_path),
                "--recommendation-mode",
                "qualified",
                "--no-synllama",
            ]
        )
        == 2
    )
