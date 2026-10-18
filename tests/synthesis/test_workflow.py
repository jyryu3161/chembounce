import json
from pathlib import Path

import pytest
from rdkit import Chem, DataStructs

pytest.importorskip(
    "rdkit.Chem.rdShapeAlign", reason="Optional synthesis alignment requires a newer RDKit wheel"
)

from chembounce_synthesis import proposal_json
from chembounce_synthesis.cli import main
from chembounce_synthesis.common import FP, canonical, digest, molecule
from chembounce_synthesis.proposal_json import POLICY, bind_row
from chembounce_synthesis.routes import verify_route
from chembounce_synthesis.synllama import load_prompt
from chembounce_synthesis.workflow import check_constraints, run

HIT = "CC(=O)Nc1ccccc1"
DESIGN = "CC(=O)Nc1ccncc1"


def config_for(tmp_path, assets, proposal):
    path = tmp_path / "proposals.json"
    path.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": DESIGN, "response": proposal}]})
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "search": {"building_block_top_k": 3, "beam_width": 4},
        "conformers": {"num_conformers": 2},
        "expansion": {"rounds": 1, "max_routes": 6},
        "max_scored_products": 10,
    }


def test_pipeline_outputs_replayable_routes_and_actual_poses(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, proposal)
    out = tmp_path / "out"
    rows = run(HIT, [DESIGN], config, out)
    assert rows and len({r["smiles"] for r in rows}) == len(rows)
    assert canonical(HIT) not in {r["smiles"] for r in rows}
    summary = json.loads((out / "summary.json").read_text())
    assert summary["status"] == "success"
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["catalog_sha256"] == assets[0].sha256
    assert manifest["provider"]["provider"] == "saved_proposals"
    routes = json.loads((out / "routes.json").read_text())["routes"]
    lookup = {r["route_id"]: r for r in routes}
    for row in rows:
        assert verify_route(lookup[row["route_id"]]["route"], *assets)["product"] == row["smiles"]
        product_fp = FP.GetFingerprint(molecule(row["smiles"]))
        for ref, reference in (("hit", HIT), ("design", DESIGN)):
            expected = DataStructs.TanimotoSimilarity(product_fp, FP.GetFingerprint(molecule(reference)))
            assert row[f"{ref}_tanimoto"] == pytest.approx(expected)
        assert row["hit_morgan"] == row["hit_tanimoto"]
    poses = list(Chem.SDMolSupplier(str(out / "poses_hit.sdf")))
    assert len(poses) == len(rows)
    assert all(m is not None and m.GetConformer().Is3D() for m in poses)
    for pose, row in zip(poses, rows):
        assert pose.GetDoubleProp("design_tanimoto") == pytest.approx(row["design_tanimoto"])
    references = list(Chem.SDMolSupplier(str(out / "references_hit.sdf")))
    assert len(references) == len(poses)
    from chembounce_synthesis.conformers import pharmacophore_overlap

    for row, pose, ref in zip(rows, poses, references):
        overlap, _ = pharmacophore_overlap(ref, pose, 0, 0)
        assert overlap == pytest.approx(row["hit_pharmacophore"], abs=1e-4)
    with pytest.raises(ValueError, match="not empty"):
        run(HIT, [DESIGN], config, out)


def test_no_result_has_header_and_explicit_reasons(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, proposal)
    config["property_thresholds"] = {"mw": [None, 10]}
    out = tmp_path / "none"
    assert run(HIT, [DESIGN], config, out) == []
    summary = json.loads((out / "summary.json").read_text())
    assert summary["status"] == "no_candidates"
    assert summary["route_rejections"]["constraint_rejected"] > 0
    assert len((out / "results.tsv").read_text().splitlines()) == 1


def test_invalid_model_response_is_visible(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, {"reactions": [], "building_blocks": []})
    out = tmp_path / "badmodel"
    assert run(HIT, [DESIGN], config, out) == []
    assert json.loads((out / "summary.json").read_text())["failures"][0]["stage"] == "proposal"
    assert (out / "proposals.json").exists()


def test_opt_in_saved_proposal_reaches_route_gate_with_raw_audit(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, proposal)
    changed = json.loads(json.dumps(proposal))
    changed["reactions"][0]["product"] = "F/C=C\\Cl"
    raw = json.dumps(changed).replace("C\\\\Cl", "C\\Cl")
    provider = {
        "response_json_dialect": POLICY,
        "prompt_sha256": "a" * 64,
        "response_decoder_source_sha256": digest(Path(proposal_json.__file__)),
    }
    row = bind_row(raw, {"seed": 1, "sample_index": 0}, "a" * 64, provider["response_decoder_source_sha256"])
    assert row["compatibility"]["status"] == "compatibility_decoded"
    source = tmp_path / "proposals.json"
    source.write_text(
        json.dumps({"schema_version": 1, "provider": provider, "proposals": [{"design": DESIGN, **row}]})
    )
    out = tmp_path / "compatibility"
    rows = run(HIT, [DESIGN], config, out)
    assert rows
    recorded = json.loads((out / "proposals.json").read_text())["proposals"][0]
    assert recorded["raw_response"] == raw
    assert recorded["compatibility"] == row["compatibility"]
    summary = json.loads((out / "summary.json").read_text())
    assert any(item["stage"] == "projection" for item in summary["search_reports"])
    assert all(item["stage"] != "proposal" for item in summary["failures"])


def test_bounded_workflow_and_cli(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, proposal)
    config.update(max_scored_products=1, max_results=1)
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    out = tmp_path / "cli"
    assert main(["--hit", HIT, "--design", DESIGN, "--config", str(config_path), "--output", str(out)]) == 0
    summary = json.loads((out / "summary.json").read_text())
    assert summary["scored_products"] == 1
    assert "max_scored_products" in summary["limits_reached"]


def test_final_structure_constraints():
    check_constraints(DESIGN, HIT, DESIGN, core_smiles="CC(=O)N")
    with pytest.raises(ValueError, match="original hit"):
        check_constraints(HIT, HIT, DESIGN)
    with pytest.raises(ValueError, match="Murcko"):
        check_constraints("CCC(=O)Nc1ccccc1", HIT, DESIGN)
    check_constraints("CCC(=O)Nc1ccncc1", HIT, DESIGN, mode="expand")
    with pytest.raises(ValueError, match="Murcko"):
        check_constraints("CCC(=O)Nc1ccccc1", HIT, DESIGN, mode="expand")
    with pytest.raises(ValueError, match="anchor"):
        check_constraints(DESIGN, HIT, DESIGN, core_smiles="S")


@pytest.mark.parametrize(
    "hit,product",
    [
        (HIT, "CC(=O)Nc1cc[13cH]cc1"),
        ("C1CC[C@H]2CCCC[C@@H]2C1", "C1CC[C@@H]2CCCC[C@@H]2C1"),
    ],
)
def test_isotope_or_stereo_change_alone_is_not_scaffold_hopping(hit, product):
    with pytest.raises(ValueError, match="Murcko"):
        check_constraints(product, hit, product)


def test_prompt_loader_does_not_execute_upstream_source(tmp_path):
    path = tmp_path / "vars.py"
    path.write_text(
        'raise RuntimeError("must never execute")\nTEMPLATE = {"instruction":"test", "input":"SMILES_STRING"}\n'
    )
    assert load_prompt(path)["instruction"] == "test"
    path.write_text('TEMPLATE = __import__("os").getcwd()')
    with pytest.raises(ValueError):
        load_prompt(path)


def test_bad_configuration_records_error(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, proposal)
    config["invented_option"] = True
    out = tmp_path / "badconfig"
    with pytest.raises(ValueError, match="Unknown synthesis"):
        run(HIT, [DESIGN], config, out)
    assert json.loads((out / "summary.json").read_text())["status"] == "error"


def test_main_cli_connects_real_hopping_and_reapplies_filters(tmp_path, assets, proposal):
    import chembounce

    config = config_for(tmp_path, assets, proposal)
    # Every positional aminopyridine that hopping may choose has a saved proposal.
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
    proposal_path = tmp_path / "all-proposals.json"
    proposal_path.write_text(json.dumps({"schema_version": 1, "proposals": proposals}))
    config["provider"]["path"] = str(proposal_path)
    config["max_scored_products"] = 2
    config_path = tmp_path / "integrated.json"
    config_path.write_text(json.dumps(config))
    replacements = tmp_path / "replacements.smi"
    replacements.write_text("c1ccncc1\n")
    out = tmp_path / "integrated"
    status = chembounce.main(
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
            "--mw_max",
            "140",
            "--core_smiles",
            "CC(=O)N",
        ]
    )
    assert status == 0
    assert (out / "overall_result.txt").is_file()
    manifest = json.loads((out / "synthesis/manifest.json").read_text())
    assert manifest["effective_workflow"]["property_thresholds"]["mw"][1] == 140
    assert manifest["effective_workflow"]["core_smiles"] == "CC(=O)N"
    routes = json.loads((out / "synthesis/routes.json").read_text())["routes"]
    assert any(
        r["selection_status"] == "constraint_rejected" and "property filter" in r["reason"] for r in routes
    )


def test_tanimoto_filter_is_opt_in_and_applied_before_3d(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, proposal)
    config["minimum_scores"] = {"design_tanimoto": 1.0}
    out = tmp_path / "tanimoto-filter"
    rows = run(HIT, [DESIGN], config, out)
    assert rows and all(row["design_tanimoto"] == 1.0 for row in rows)
    summary = json.loads((out / "summary.json").read_text())
    assert summary["route_rejections"]["2d_threshold_rejected"] > 0
    assert summary["scored_products"] == 1


def test_pareto_uses_tanimoto_alongside_3d():
    from chembounce_synthesis.workflow import _pareto

    reference = dict(
        hit_shape=0.8,
        hit_pharmacophore=0.6,
        design_shape=0.8,
        design_pharmacophore=0.6,
        relative_energy=0,
        reaction_steps=1,
        hit_tanimoto=0.5,
        design_tanimoto=0.5,
    )
    rows = [reference.copy(), reference | {"design_tanimoto": 0.8}]
    _pareto(rows)
    assert [r["pareto_front"] for r in rows] == [1, 0]


def test_standalone_cli_expansion_overrides(tmp_path, assets, proposal):
    config = config_for(tmp_path, assets, proposal)
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    out = tmp_path / "expanded-cli"
    assert (
        main(
            [
                "--hit",
                HIT,
                "--design",
                DESIGN,
                "--config",
                str(path),
                "--output",
                str(out),
                "--expansion-top-k",
                "6",
                "--expansion-rounds",
                "2",
                "--expansion-max-routes",
                "20",
            ]
        )
        == 0
    )
    effective = json.loads((out / "manifest.json").read_text())["effective_workflow"]["expansion"]
    assert effective["building_block_top_k"] == 6
    assert effective["rounds"] == 2
    assert effective["max_routes"] == 20


def test_refined_outputs_are_actual_median_seed_poses(tmp_path, assets, proposal):
    from chembounce_synthesis.conformers import pharmacophore_overlap

    config = config_for(tmp_path, assets, proposal)
    config["refinement"] = {"top_n": 2, "num_conformers": 4, "seeds": [43, 44, 45]}
    out = tmp_path / "refined"
    rows = run(HIT, [DESIGN], config, out)
    assert len(rows) == 2
    reports = json.loads((out / "refinement.json").read_text())["candidates"]
    by_id = {r["route_id"]: r for r in reports}
    for role in ("hit", "design"):
        poses = list(Chem.SDMolSupplier(str(out / f"poses_{role}.sdf")))
        references = list(Chem.SDMolSupplier(str(out / f"references_{role}.sdf")))
        for row, pose, reference in zip(rows, poses, references):
            assert row["scoring_stage"] == "refined"
            assert row["seed_count"] == row["seed_successes"] == 3
            assert row["pose_seed"] in (43, 44, 45)
            trials = by_id[row["route_id"]]["trials"]
            assert row["score"] == row["seed_score_median"] == sorted(r["score"] for r in trials)[1]
            assert pose.GetIntProp("pose_seed") == row["pose_seed"]
            overlap, _ = pharmacophore_overlap(reference, pose, 0, 0)
            assert overlap == pytest.approx(row[f"{role}_pharmacophore"], abs=1e-4)
    assert len((out / "screening.tsv").read_text().splitlines()) > len(rows) + 1


def test_refinement_defers_3d_thresholds_until_reevaluation(tmp_path, assets, proposal, monkeypatch):
    import chembounce_synthesis.workflow as workflow

    real_score = workflow.score_candidate

    def poor_screen(*args, **kwargs):
        info, *poses = real_score(*args, **kwargs)
        info["hit"]["shape"] = 0.0
        return info, *poses

    monkeypatch.setattr(workflow, "score_candidate", poor_screen)
    config = config_for(tmp_path, assets, proposal)
    config["minimum_scores"] = {"hit_shape": 0.1}
    config["refinement"] = {"top_n": 1, "num_conformers": 3, "seeds": [42, 43, 44]}
    rows = run(HIT, [DESIGN], config, tmp_path / "rescue")
    assert rows and rows[0]["hit_shape"] > 0.1


def test_continuation_retains_candidates_without_model_and_rechecks_core(
    tmp_path, assets, proposal, monkeypatch
):
    from chembounce_synthesis.common import digest
    from chembounce_synthesis.continuation import continuation_inputs

    config = config_for(tmp_path, assets, proposal)
    first = tmp_path / "first"
    old_rows = run(HIT, [DESIGN], config, first)
    before = {p.name: digest(p) for p in first.iterdir()}
    # Continuing must use saved responses even when the recorded run used a model.
    manifest_path = first / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["config"]["provider"] = {"kind": "synllama", "model_path": "not-installed"}
    manifest_path.write_text(json.dumps(manifest))
    before["manifest.json"] = digest(manifest_path)
    monkeypatch.setattr(
        "chembounce_synthesis.workflow.SynLlama", lambda **kwargs: pytest.fail("model loaded")
    )
    second = tmp_path / "second"
    assert main(["--continue-from", str(first), "--output", str(second), "--expansion-rounds", "0"]) == 0
    import csv

    with (second / "screening.tsv").open() as stream:
        retained = {r["smiles"] for r in csv.DictReader(stream, delimiter="\t")}
    assert {r["smiles"] for r in old_rows} <= retained
    assert {p.name: digest(p) for p in first.iterdir()} == before
    assert json.loads((second / "manifest.json").read_text())["provider"]["provider"] == "saved_proposals"
    hit, designs, resumed = continuation_inputs(second)
    resumed["core_smiles"] = "c1ccccc1"  # Present in hit; only some expanded products retain a phenyl ring.
    third = tmp_path / "third"
    rows = run(hit, designs, resumed, third)
    assert len(rows) < len(retained)
    assert all(molecule(row["smiles"]).HasSubstructMatch(molecule(resumed["core_smiles"])) for row in rows)
    records = json.loads((third / "routes.json").read_text())["routes"]
    assert any(r.get("carried_from") and r["selection_status"] == "constraint_rejected" for r in records)


def test_carried_routes_are_not_trusted_and_hit_must_match(tmp_path, assets, proposal):
    from chembounce_synthesis.continuation import continuation_inputs

    first = tmp_path / "first"
    run(HIT, [DESIGN], config_for(tmp_path, assets, proposal), first)
    path = first / "routes.json"
    data = json.loads(path.read_text())
    # A saved verification claim cannot authorize an unregistered reaction.
    for record in data["routes"]:
        for node in record["route"]["nodes"]:
            if "rule_id" in node:
                node["rule_id"] = "invented"
    path.write_text(json.dumps(data))
    hit, designs, config = continuation_inputs(first)
    second = tmp_path / "second"
    assert run(hit, designs, config, second)
    records = json.loads((second / "routes.json").read_text())["routes"]
    assert all(r["selection_status"] == "verification_rejected" for r in records if r.get("carried_from"))
    with pytest.raises(ValueError, match="different original hit"):
        run("CCC(=O)Nc1ccccc1", designs, config, tmp_path / "wrong-hit")


def test_continuation_can_score_previously_unscored_routes_and_override_refinement(
    tmp_path, assets, proposal
):
    config = config_for(tmp_path, assets, proposal)
    config["max_scored_products"] = 1
    first = tmp_path / "limited"
    assert len(run(HIT, [DESIGN], config, first)) == 1
    second = tmp_path / "larger"
    assert (
        main(
            [
                "--continue-from",
                str(first),
                "--output",
                str(second),
                "--expansion-rounds",
                "0",
                "--synthesis-max-scored",
                "10",
                "--refine-top-n",
                "4",
                "--refine-conformers",
                "3",
                "--refine-seeds",
                "43",
                "44",
            ]
        )
        == 0
    )
    manifest = json.loads((second / "manifest.json").read_text())
    effective_refinement = manifest["effective_workflow"]["refinement"]
    assert {key: effective_refinement[key] for key in ("top_n", "num_conformers", "seeds")} == {
        "top_n": 4,
        "num_conformers": 3,
        "seeds": [43, 44],
    }
    assert effective_refinement["threshold_policy"] == "all_seeds"
    summary = json.loads((second / "summary.json").read_text())
    assert summary["returned_products"] == 4
    assert summary["carried_screened_products"] > 1


def test_refinement_failure_keeps_screening_pool_without_exporting_coarse_pose(
    tmp_path, assets, proposal, monkeypatch
):
    monkeypatch.setattr(
        "chembounce_synthesis.workflow.refine_candidate",
        lambda *a, **kw: (None, {"status": "failed", "trials": []}),
    )
    config = config_for(tmp_path, assets, proposal)
    config["refinement"] = {"top_n": 1}
    out = tmp_path / "refinement-failed"
    assert run(HIT, [DESIGN], config, out) == []
    summary = json.loads((out / "summary.json").read_text())
    assert summary["screened_products"] > 0
    assert summary["status"] == "no_candidates"
    assert summary["refinement_attempted"] == 1
    assert summary["route_rejections"]["refinement_failed"] >= 1
    assert len((out / "results.tsv").read_text().splitlines()) == 1
    assert len((out / "screening.tsv").read_text().splitlines()) > 1


def test_continuation_scores_all_designs_carried_routes_before_new_candidates(tmp_path, assets, proposal):
    from chembounce_synthesis.continuation import continuation_inputs

    second_design = "CC(=O)Nc1ccccn1"
    second_proposal = {
        "reactions": [
            {"reaction_template": "amide", "reactants": ["CC(=O)O", "Nc1ccccn1"], "product": second_design}
        ],
        "building_blocks": ["CC(=O)O", "Nc1ccccn1"],
    }
    config = config_for(tmp_path, assets, proposal)
    config["search"] = {"building_block_top_k": 1, "beam_width": 1}
    config["expansion"] = {"rounds": 0}
    config["max_scored_products"] = 2
    path = tmp_path / "proposals.json"
    data = json.loads(path.read_text())
    data["proposals"].append({"design": second_design, "response": second_proposal})
    path.write_text(json.dumps(data))
    first = tmp_path / "first"
    old = run(HIT, [DESIGN, second_design], config, first)
    assert len(old) == 2
    hit, designs, resumed = continuation_inputs(first)
    resumed["search"] = {"building_block_top_k": 3, "beam_width": 4}
    continued = run(hit, designs, resumed, tmp_path / "continued")
    assert {r["smiles"] for r in old} <= {r["smiles"] for r in continued}


def test_cli_reports_model_runtime_failure_without_traceback(tmp_path, assets, proposal, monkeypatch, caplog):
    def fail(config):
        raise RuntimeError("Install the synllama optional dependency group")

    monkeypatch.setattr("chembounce_synthesis.workflow._provider", fail)
    config = config_for(tmp_path, assets, proposal)
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    out = tmp_path / "runtime-error"
    assert main(["--hit", HIT, "--design", DESIGN, "--config", str(path), "--output", str(out)]) == 2
    assert "Install the synllama" in caplog.text
    assert json.loads((out / "summary.json").read_text())["status"] == "error"
