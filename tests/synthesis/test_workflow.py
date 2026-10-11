import json

import pytest
from rdkit import Chem

pytest.importorskip(
    "rdkit.Chem.rdShapeAlign", reason="Optional synthesis alignment requires a newer RDKit wheel"
)

from chembounce_synthesis.cli import main
from chembounce_synthesis.common import canonical
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
    poses = list(Chem.SDMolSupplier(str(out / "poses_hit.sdf")))
    assert len(poses) == len(rows)
    assert all(m is not None and m.GetConformer().Is3D() for m in poses)
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
