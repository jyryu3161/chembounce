"""Reject only unusable design references before spending synthesis/scoring work."""

import json
from unittest.mock import Mock

import numpy as np
import pytest
from rdkit import Chem

from chembounce_synthesis import conformers, workflow
from chembounce_synthesis.common import canonical
from chembounce_synthesis.continuation import continuation_inputs

HIT = "CC(=O)Nc1ccccc1"
DESIGN = "CC(=O)Nc1ccncc1"
OTHER = "CCC(=O)Nc1ccncc1"
UNSPECIFIED = "CC(O)c1ccncc1"
CAGE = "NC12CC3CC(CC(C3)C1)C2"
REAL_HIT = "COc1ccc2c(c1)CCCN2c1nc(Cl)nc2ccccc12"
REAL_DESIGN = "COC1CC1CCNS(=O)(=O)CC1CCOCC1Cl"
REAL_PRODUCT = "CO[C@H]1C[C@@H](CCNS(=O)(=O)CC2CCOCC2)C1"


def read(output, name):
    return json.loads((output / name).read_text())


def config_for(tmp_path, assets, proposal):
    saved = tmp_path / "saved.json"
    saved.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "proposals": [
                    {"design": design, "response": proposal, "sampling": {"seed": 17, "sample": 0}}
                    for design in [DESIGN, OTHER, UNSPECIFIED]
                ],
            }
        )
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(saved)},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "conformers": {"num_conformers": 2, "max_iterations": 1000, "seed": 17},
        "max_scored_products": 1,
    }


def forbid_assets(monkeypatch):
    calls = {}
    for name in (
        "Catalog",
        "Rules",
        "_provider",
        "_evidence_evaluator",
        "resolve_proposal",
        "verify_route",
        "score_candidate",
        "constraints_from_config",
    ):
        calls[name] = Mock(side_effect=AssertionError(f"Unexpected expensive work: {name}"))
        monkeypatch.setattr(workflow, name, calls[name])
    calls["embedding"] = Mock(side_effect=AssertionError("Preflight cannot embed or assign stereo"))
    monkeypatch.setattr(conformers.AllChem, "EmbedMultipleConfs", calls["embedding"])
    return calls


def carry_source(tmp_path, config, route, designs, hit=HIT):
    source = tmp_path / "previous"
    source.mkdir()
    (source / "manifest.json").write_text(
        json.dumps({"schema_version": 1, "hit": hit, "designs": designs, "config": config})
    )
    (source / "routes.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "routes": [
                    {
                        "route_id": f"old_{index}",
                        "design": design,
                        "route": route,
                        # Neither this stale status nor fabricated evidence may be trusted.
                        "selection_status": "accepted",
                        "verification": {"product": "CC"},
                    }
                    for index, design in enumerate(designs)
                ],
            }
        )
    )
    return source


@pytest.mark.parametrize("designs", [[UNSPECIFIED], [UNSPECIFIED, CAGE], [REAL_DESIGN]])
def test_all_invalid_designs_skip_assets_and_export_explicit_empty_outcome(
    tmp_path, assets, proposal, monkeypatch, designs
):
    monkeypatch.setattr(conformers, "rdShapeAlign", object())
    config = config_for(tmp_path, assets, proposal)
    # Even missing asset files must not be opened when no capped design can score.
    config.update(catalog="absent-catalog.tsv", rules="absent-rules.json")
    config["provider"]["path"] = "absent-proposals.json"
    calls = forbid_assets(monkeypatch)
    output = tmp_path / "all-invalid"
    hit = REAL_HIT if designs == [REAL_DESIGN] else HIT
    assert workflow.run(hit, designs, config, output) == []
    assert not any(call.called for call in calls.values())
    summary, manifest = read(output, "summary.json"), read(output, "manifest.json")
    assert summary["status"] == "no_candidates"
    assert summary["scored_products"] == summary["verified_routes"] == summary["proposals"] == 0
    assert summary["designs_used"] == summary["designs_preflight_rejected"] == len(designs)
    assert summary["designs_eligible_for_scoring"] == 0
    assert "All capped designs" in summary["no_candidates_reason"]
    assert summary["designs_not_projected"] == [canonical(design) for design in designs]
    assert manifest["designs"] == [canonical(design) for design in designs]
    assert manifest["design_preflight"]["requested_designs"] == designs
    assert manifest["design_preflight"]["eligible_designs"] == []
    assert manifest["design_preflight"]["checks"] == summary["failures"]
    for index, failure in enumerate(summary["failures"]):
        assert failure["design"] == canonical(designs[index])
        assert failure["role"] == "design" and failure["capped_design_index"] == index
        assert failure["stage"] == "design_reference_preflight" and failure["status"] == "rejected"
        assert "specified stereochemistry" in failure["reason"]
    assert manifest["provider"] is None
    assert manifest["catalog_sha256"] is manifest["rules_sha256"] is manifest["catalog_import"] is None
    assert manifest["synthesis_evidence"]["status"] == "not_initialized"
    assert manifest["assets_not_loaded"] == ["catalog", "rules", "provider", "synthesis_evidence"]
    assert read(output, "proposals.json")["proposals"] == []
    assert read(output, "routes.json")["routes"] == []
    assert read(output, "scores.json")["scores"] == []
    for table in ["results.tsv", "screening.tsv", "screening_alternatives.tsv"]:
        assert len((output / table).read_text().splitlines()) == 1
    for sdf in ["poses_hit.sdf", "poses_design.sdf", "references_hit.sdf", "references_design.sdf"]:
        assert (output / sdf).read_text() == ""


def test_cap_is_applied_before_rejection_without_backfill_and_preserves_requested_inputs(
    tmp_path, assets, proposal, monkeypatch
):
    monkeypatch.setattr(conformers, "rdShapeAlign", object())
    config = config_for(tmp_path, assets, proposal)
    config["max_designs"] = 1
    designs = [UNSPECIFIED, "OC(C)c1ccncc1", DESIGN]
    calls = forbid_assets(monkeypatch)
    output = tmp_path / "capped"
    assert workflow.run(HIT, iter(designs), config, output) == []
    assert not any(call.called for call in calls.values())
    manifest, summary = read(output, "manifest.json"), read(output, "summary.json")
    assert manifest["designs"] == [canonical(UNSPECIFIED)]
    preflight = manifest["design_preflight"]
    assert preflight["requested_designs"] == designs
    assert preflight["canonical_unique_requested_designs"] == [canonical(UNSPECIFIED), DESIGN]
    assert len(preflight["checks"]) == 1 and preflight["eligible_designs"] == []
    assert summary["designs_supplied"] == 2 and summary["designs_used"] == 1
    assert "max_designs" in summary["limits_reached"]


def test_mixed_designs_preserve_valid_results_poses_sampling_and_budget(
    tmp_path, assets, proposal, monkeypatch
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)
    provider_factory, evidence_factory = workflow._provider, workflow._evidence_evaluator
    provider_calls, evidence_calls = [], []

    def provider_factory_spy(config):
        provider = provider_factory(config)
        propose = provider.propose

        def record(design):
            provider_calls.append(design)
            return propose(design)

        provider.propose = record
        return provider

    def evidence_factory_spy(config):
        evaluator = evidence_factory(config)
        evaluate = evaluator.evaluate

        def record(verification):
            evidence_calls.append(verification)
            return evaluate(verification)

        evaluator.evaluate = record
        return evaluator

    monkeypatch.setattr(workflow, "_provider", provider_factory_spy)
    monkeypatch.setattr(workflow, "_evidence_evaluator", evidence_factory_spy)
    baseline, mixed = tmp_path / "baseline", tmp_path / "mixed"
    conformers.ensemble.cache_clear()
    expected = workflow.run(HIT, [DESIGN, OTHER], config, baseline)
    baseline_evidence = list(evidence_calls)
    provider_calls.clear()
    evidence_calls.clear()
    conformers.ensemble.cache_clear()
    actual = workflow.run(HIT, [UNSPECIFIED, DESIGN, CAGE, OTHER], config, mixed)
    assert actual == expected and len(actual) == 1
    assert provider_calls == [DESIGN, OTHER]
    assert evidence_calls == baseline_evidence
    assert read(mixed, "proposals.json")["proposals"] == read(baseline, "proposals.json")["proposals"]
    assert read(mixed, "scores.json") == read(baseline, "scores.json")
    assert read(mixed, "routes.json") == read(baseline, "routes.json")
    for sdf in ["poses_hit.sdf", "poses_design.sdf", "references_hit.sdf", "references_design.sdf"]:
        before = next(iter(Chem.SDMolSupplier(str(baseline / sdf))))
        after = next(iter(Chem.SDMolSupplier(str(mixed / sdf))))
        assert np.array_equal(before.GetConformer().GetPositions(), after.GetConformer().GetPositions())
    summary = read(mixed, "summary.json")
    assert summary["scored_products"] == 1 and summary["designs_preflight_rejected"] == 2
    assert summary["designs_eligible_for_scoring"] == 2
    assert [entry["design"] for entry in summary["screening_allocation"]] == [DESIGN]


@pytest.mark.parametrize("mode,allow", [("3d", True), ("2d", False), ("2d", True)])
def test_permissive_and_2d_designs_still_produce_results(
    tmp_path, assets, proposal, monkeypatch, mode, allow
):
    if mode == "3d":
        pytest.importorskip("rdkit.Chem.rdShapeAlign")
    else:
        monkeypatch.setattr(conformers, "rdShapeAlign", None)
        monkeypatch.setattr(workflow, "prepare_reference", Mock(side_effect=AssertionError("2D preflight")))
        monkeypatch.setattr(workflow, "score_candidate", Mock(side_effect=AssertionError("2D scoring")))
    config = config_for(tmp_path, assets, proposal)
    config["conformers"]["reject_unspecified_stereo"] = not allow
    output = tmp_path / "permissive"
    rows = workflow.run(HIT, [UNSPECIFIED], config, output, scoring_mode=mode)
    assert rows and rows[0]["design"] == canonical(UNSPECIFIED)
    summary, manifest = read(output, "summary.json"), read(output, "manifest.json")
    assert summary["status"] == "success" and summary["designs_preflight_rejected"] == 0
    assert summary["scored_products"] == 1 and summary["failures"] == []
    assert manifest["design_preflight"]["checks"][0]["status"] == (
        "eligible" if mode == "3d" else "not_applicable_2d"
    )
    assert read(output, "proposals.json")["proposals"][0]["sampling"]["seed"] == 17
    assert (output / "poses_hit.sdf").exists() == (mode == "3d")


def test_mixed_carried_routes_reject_only_invalid_reference_and_preserve_old_graph(
    tmp_path, assets, proposal, route
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)
    source = carry_source(tmp_path, config, route, [UNSPECIFIED, DESIGN])
    config["carry_forward"] = [str(source)]
    output = tmp_path / "carried"
    rows = workflow.run(HIT, [UNSPECIFIED, DESIGN], config, output)
    assert len(rows) == 1 and rows[0]["design"] == DESIGN
    records = read(output, "routes.json")["routes"]
    rejected = [record for record in records if record["design"] == canonical(UNSPECIFIED)]
    assert len(rejected) == 1
    rejected = rejected[0]
    assert rejected["route"] == route and rejected["carried_from"]["route_id"] == "old_0"
    assert rejected["selection_status"] == "design_preflight_rejected"
    assert rejected["reference_preflight"]["role"] == "design"
    assert rejected["verification_status"] == "not_evaluated"
    assert not rejected["eligible_for_evidence"] and not rejected["eligible_for_scoring"]
    assert "verification" not in rejected and "synthesis_evidence" not in rejected
    primary = next(record for record in records if record["route_id"] == rows[0]["route_id"])
    assert primary["carried_from"]["route_id"] == "old_1"
    assert primary["verification"]["product"] == DESIGN
    summary = read(output, "summary.json")
    assert summary["carried_routes_loaded"] == 2 and summary["carried_screened_products"] == 1
    assert summary["route_rejections"]["design_preflight_rejected"] == 1
    assert continuation_inputs(output)[1] == [canonical(UNSPECIFIED), DESIGN]


def test_all_invalid_carried_routes_remain_available_for_explicit_2d_continuation(
    tmp_path, assets, proposal, route, monkeypatch
):
    config = config_for(tmp_path, assets, proposal)
    source = carry_source(tmp_path, config, route, [UNSPECIFIED])
    config["carry_forward"] = [str(source)]
    output = tmp_path / "rejected"
    with monkeypatch.context() as patch:
        patch.setattr(conformers, "rdShapeAlign", object())
        calls = forbid_assets(patch)
        assert workflow.run(HIT, [UNSPECIFIED], config, output) == []
        assert not any(call.called for call in calls.values())
    records = read(output, "routes.json")["routes"]
    assert len(records) == 1 and records[0]["route"] == route
    assert records[0]["verification_status"] == "not_evaluated"
    manifest = read(output, "manifest.json")
    assert manifest["carry_forward"][0]["loaded_routes"] == 1
    hit, designs, continued = continuation_inputs(output)
    # Preflight did not invent a saved proposal; explicitly supply the old source.
    continued["provider"] = config["provider"]
    next_output = tmp_path / "continued"
    rows = workflow.run(hit, designs, continued, next_output, scoring_mode="2d")
    assert len(rows) == 1 and rows[0]["smiles"] == DESIGN
    primary = next(
        record
        for record in read(next_output, "routes.json")["routes"]
        if record["route_id"] == rows[0]["route_id"]
    )
    assert primary["carried_from"]["run"] == str(output)
    assert primary["verification"]["product"] == DESIGN


def test_all_invalid_does_not_bypass_carried_source_hit_validation(
    tmp_path, assets, proposal, route, monkeypatch
):
    monkeypatch.setattr(conformers, "rdShapeAlign", object())
    config = config_for(tmp_path, assets, proposal)
    source = carry_source(tmp_path, config, route, [UNSPECIFIED], hit=DESIGN)
    config["carry_forward"] = [str(source)]
    calls = forbid_assets(monkeypatch)
    output = tmp_path / "wrong-source"
    with pytest.raises(ValueError, match="different original hit"):
        workflow.run(HIT, [UNSPECIFIED], config, output)
    assert not any(call.called for call in calls.values())
    assert read(output, "summary.json")["status"] == "error"


def test_real_diagnostic_is_a_design_reference_failure_not_a_hit_or_product_failure(monkeypatch):
    # This observed case uses the supported 3D runtime's stereo perception.
    # RDKit 2023 (without rdShapeAlign) perceives an extra potential center here.
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    monkeypatch.setattr(conformers, "rdShapeAlign", object())
    monkeypatch.setattr(
        conformers.AllChem, "EmbedMultipleConfs", Mock(side_effect=AssertionError("embedding"))
    )
    assert conformers.prepare_reference(REAL_HIT) is None
    assert conformers.prepare_reference(REAL_PRODUCT) is None
    with pytest.raises(ValueError, match="specified stereochemistry"):
        conformers.prepare_reference(REAL_DESIGN)


@pytest.mark.parametrize("with_pocket", [False, True])
def test_bound_hit_never_supplies_a_design_reference_pose(
    tmp_path, assets, proposal, monkeypatch, with_pocket
):
    monkeypatch.setattr(conformers, "rdShapeAlign", object())
    config = config_for(tmp_path, assets, proposal)
    mol = Chem.MolFromSmiles(HIT)
    coordinates = Chem.Conformer(mol.GetNumAtoms())
    coordinates.Set3D(True)
    for index in range(mol.GetNumAtoms()):
        coordinates.SetAtomPosition(index, (float(index), 0.0, 0.0))
    mol.AddConformer(coordinates)
    sdf = tmp_path / "bound.sdf"
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(mol)
    config["reference_sdf"] = str(sdf)
    if with_pocket:
        config["pocket"] = {
            "receptor_pdb": "absent-receptor.pdb",
            "reference_ligand": {"chain": "A", "residue_number": 1, "residue_name": "LIG"},
        }
    calls = forbid_assets(monkeypatch)
    preflight = Mock(wraps=workflow.prepare_reference)
    monkeypatch.setattr(workflow, "prepare_reference", preflight)
    output = tmp_path / "bound-hit"
    assert workflow.run(HIT, [UNSPECIFIED], config, output) == []
    assert not any(call.called for call in calls.values())
    assert preflight.call_args_list[0].args == (
        HIT,
        conformers.ConformerSettings(**config["conformers"]),
        str(sdf),
    )
    assert preflight.call_args_list[1].args == (canonical(UNSPECIFIED), preflight.call_args_list[0].args[1])
    manifest = read(output, "manifest.json")
    assert manifest["reference_sdf_sha256"]
    assert ("pocket" in manifest["assets_not_loaded"]) == with_pocket
    assert manifest["pocket"] is None
