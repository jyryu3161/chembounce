"""Cooperative time limits preserve complete evidence and reject partial seeds."""

import json
import time
from types import SimpleNamespace

import pytest

from chembounce_synthesis import refinement, workflow
from chembounce_synthesis.conformers import ConformerSettings
from chembounce_synthesis.refinement import RefinementBudget, RefinementSettings


def test_deadline_interrupting_seed_block_never_returns_success(monkeypatch):
    clock = {"now": 0.0}
    monkeypatch.setattr(refinement, "time", SimpleNamespace(monotonic=lambda: clock["now"]))

    def score(smiles, hit, design, settings, pose):
        clock["now"] = 10.0
        return {"score": 0.9}, None, None

    monkeypatch.setattr(refinement, "score_candidate", score)
    budget = RefinementBudget(100)
    result, report = refinement.refine_candidate(
        "CC",
        "CC",
        "CC",
        ConformerSettings(),
        RefinementSettings(num_conformers=2),
        budget=budget,
        deadline_monotonic=5.0,
    )
    assert result is None and report["status"] == "deadline_exhausted"
    assert report["convergence_status"] == "incomplete_seed_block"
    assert report["seed_count"] == 3 and report["seed_successes"] == 1
    assert [trial["status"] for trial in report["trials"]] == ["success", "not_evaluated", "not_evaluated"]
    assert budget.requested == 2 and budget.score_calls == 1


def test_expired_refinement_does_not_consume_budget_or_run_score(monkeypatch):
    monkeypatch.setattr(refinement, "score_candidate", lambda *a, **k: pytest.fail("expired score called"))
    budget = RefinementBudget(100)
    result, report = refinement.refine_candidate(
        "CC",
        "CC",
        "CC",
        ConformerSettings(),
        RefinementSettings(),
        budget=budget,
        deadline_monotonic=time.monotonic() - 1,
    )
    assert result is None and report["status"] == "deadline_exhausted"
    assert budget.requested == 0 and budget.score_calls == 0


def config_for(tmp_path, assets, proposal):
    path = tmp_path / "proposals.json"
    path.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": "CC(=O)Nc1ccncc1", "response": proposal}]})
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "search": {"building_block_top_k": 3, "beam_width": 4},
        "conformers": {"num_conformers": 2},
    }


def test_expired_workflow_records_unassessed_design_without_loading_provider(
    tmp_path, assets, proposal, monkeypatch
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    monkeypatch.setattr(workflow, "_provider", lambda *a, **k: pytest.fail("expired provider loaded"))
    monkeypatch.setattr(
        workflow, "_evidence_evaluator", lambda *a, **k: pytest.fail("expired evidence model loaded")
    )
    out = tmp_path / "expired"
    rows = workflow.run(
        "CC(=O)Nc1ccccc1",
        ["CC(=O)Nc1ccncc1"],
        config_for(tmp_path, assets, proposal),
        out,
        deadline_monotonic=time.monotonic() - 1,
    )
    assert rows == []
    report = json.loads((out / "summary.json").read_text())
    assert report["status"] == "no_candidates" and report["deadline_reached"]
    assert report["designs_not_projected"] == ["CC(=O)Nc1ccncc1"]
    assert (out / "results.tsv").read_text().count("\n") == 1
    assert json.loads((out / "manifest.json").read_text())["deadline_policy"]["enabled"]


def test_complete_screened_candidate_survives_deadline_before_next_candidate(
    tmp_path, assets, proposal, monkeypatch
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    start = time.monotonic()
    clock = {"now": start}
    monkeypatch.setattr(workflow, "time", SimpleNamespace(monotonic=lambda: clock["now"]))
    original = workflow.score_candidate

    def score(*args, **kwargs):
        result = original(*args, **kwargs)
        clock["now"] = start + 200
        return result

    monkeypatch.setattr(workflow, "score_candidate", score)
    out = tmp_path / "partial"
    rows = workflow.run(
        "CC(=O)Nc1ccccc1",
        ["CC(=O)Nc1ccncc1"],
        config_for(tmp_path, assets, proposal),
        out,
        deadline_monotonic=start + 100,
    )
    assert len(rows) == 1
    report = json.loads((out / "summary.json").read_text())
    assert report["status"] == "success" and report["deadline_reached"]
    assert report["scored_products"] == 1
    assert report["route_rejections"]["not_scored_deadline"] >= 1
    assert len((out / "results.tsv").read_text().splitlines()) == 2
    from rdkit import Chem

    assert len(list(Chem.SDMolSupplier(str(out / "poses_hit.sdf")))) == 1
