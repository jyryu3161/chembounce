"""Evidence scheduling through real toy route replay; controlled predictors only."""

import copy
import json
import time
from types import SimpleNamespace

import pytest

from chembounce_synthesis import profiles, workflow
from chembounce_synthesis.common import canonical
from chembounce_synthesis.evidence import SynthesisEvidenceEvaluator

OPTION = profiles.DESIGN_EVIDENCE_ALLOCATION
HIT = "CC(=O)Nc1ccccc1"
FIRST = "CC(=O)Nc1ccncc1"
LATER = "CCC(=O)Nc1cccnc1"


def read(out, name):
    return json.loads((out / name).read_text())


def configuration(tmp_path, assets, *, duplicate=False):
    with assets[0].path.open("a") as stream:
        stream.write("furyl\tNCc1ccoc1\ttest fixture\t2026-10-17\n")
        stream.write("pyrimidyl\tNc1ncccn1\ttest fixture\t2026-10-17\n")
    definitions = [
        (FIRST, "CC(=O)O", "Nc1ccncc1"),
        (FIRST, "CC(=O)O", "NCc1ccoc1"),
        (FIRST, "CC(=O)O", "Nc1ncccn1"),
        (LATER, "CCC(=O)O", "Nc1cccnc1"),
    ]
    if duplicate:
        definitions.append((LATER, "CC(=O)O", "Nc1ccncc1"))
    proposals = []
    for design, acid, amine in definitions:
        product = acid.removesuffix("O") + amine
        proposals.append(
            dict(
                design=design,
                response={
                    "reactions": [
                        {"reaction_template": "amide", "reactants": [acid, amine], "product": product}
                    ],
                    "building_blocks": [acid, amine],
                },
            )
        )
    path = tmp_path / "proposals.json"
    path.write_text(json.dumps({"schema_version": 1, "proposals": proposals}))
    return dict(
        schema_version=1,
        catalog=str(assets[0].path),
        rules=str(assets[1].path),
        provider={"kind": "saved", "path": str(path)},
        search={"building_block_top_k": 1, "beam_width": 1},
        conformers={"num_conformers": 2},
        max_scored_products=5,
        synthesis_evidence={"policy": {"require_forward_agreement": True}, "max_predictor_calls": 2},
    )


def engine(monkeypatch):
    calls = []

    class Predictor:
        provenance = {"provider": "controlled_evidence_allocation", "uses_conditions": False}

        def token_budget(self, *args):
            return {"input_tokens": 8, "max_output_tokens": 24}

        def predict(self, reactants, *args):
            calls.append(tuple(reactants))
            acid = next(r for r in reactants if r.endswith("(=O)O"))
            amine = next(r for r in reactants if r.startswith("N"))
            return {
                "products": [canonical(acid.removesuffix("O") + amine)],
                "complete": [True],
                "input_tokens": 8,
                "output_tokens": 12,
            }

    evaluator = SynthesisEvidenceEvaluator(
        Predictor(), {"require_forward_agreement": True}, max_predictor_calls=2
    )
    monkeypatch.setattr(workflow, "_evidence_evaluator", lambda config: evaluator)
    return evaluator, calls


@pytest.mark.parametrize("mode", ["2d", "3d"])
@pytest.mark.parametrize("allocation", ["balanced", "encounter"])
def test_actual_workflow_coverage_default_parity_and_independent_screen_order(
    tmp_path, assets, monkeypatch, mode, allocation
):
    if mode == "3d":
        pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = configuration(tmp_path, assets)
    snapshots = []
    for name, policy in [("omitted", None), ("inherit", "inherit-v1"), ("rotate", OPTION)]:
        current = copy.deepcopy(config)
        current["allocation_policy"] = allocation
        if policy is not None:
            current["evidence_allocation_policy"] = policy
        evaluator, calls = engine(monkeypatch)
        out = tmp_path / name
        rows = workflow.run(HIT, [FIRST, LATER], current, out, scoring_mode=mode)
        records = read(out, "routes.json")["routes"]
        summary = read(out, "summary.json")
        snapshots.append((rows, records, summary, calls))
        assert evaluator.usage["predictor_calls"] == 2
        assert len(records) == 4
        assert evaluator.usage["cache_hits"] == 0
        assert summary["evidence_allocation_policy"] == (policy or "inherit-v1")
        manifest = read(out, "manifest.json")
        assert manifest["effective_workflow"]["evidence_allocation_policy"] == (policy or "inherit-v1")
        assert manifest["effective_workflow"]["screening_allocation_policy"] == "inherit-v1"
    assert snapshots[0][3] == snapshots[1][3]
    assert snapshots[0][2]["screening_allocation"] == snapshots[1][2]["screening_allocation"]
    assert [r["smiles"] for r in snapshots[0][0]] == [r["smiles"] for r in snapshots[1][0]]
    assert snapshots[0][1][3]["synthesis_evidence"]["steps"][0]["forward_prediction"]["status"] == "unknown"
    assert snapshots[2][1][3]["eligible_for_scoring"]
    assert "CCC(=O)O" in snapshots[2][3][1]
    assert snapshots[2][1][1]["synthesis_evidence"]["steps"][0]["forward_prediction"]["status"] == "unknown"
    # Evidence rotation is FIRST,LATER; original screening order is retained.
    screen = [(r["product"], r["design"]) for r in snapshots[2][2]["screening_allocation"]]
    expected = [(canonical(FIRST), canonical(FIRST)), (canonical(LATER), canonical(LATER))]
    assert screen == expected
    assert any(r["smiles"] == canonical(LATER) for r in snapshots[2][0])


def test_duplicate_design_associations_share_input_cache_without_claiming_equal_coverage(
    tmp_path, assets, monkeypatch
):
    config = dict(configuration(tmp_path, assets, duplicate=True), evidence_allocation_policy=OPTION)
    evaluator, calls = engine(monkeypatch)
    out = tmp_path / "duplicates"
    workflow.run(HIT, [FIRST, LATER], config, out, scoring_mode="2d")
    records = read(out, "routes.json")["routes"]
    assert len(records) == 5 and len(calls) == 2 and evaluator.usage["cache_hits"] == 1
    assert records[0]["eligible_for_scoring"] and records[4]["eligible_for_scoring"]
    assert records[0]["design"] != records[4]["design"]
    assert records[1]["synthesis_evidence"]["steps"][0]["forward_prediction"]["status"] == "unknown"
    assert (
        "Configured evaluator cache"
        in read(out, "manifest.json")["evidence_allocation_policy"]["cache_scope"]
    )


def test_one_design_per_carried_tier_is_exact_noop_and_lists_are_not_mutated(monkeypatch):
    monkeypatch.setattr(workflow, "scaffold_family_key", lambda product: product[0])

    def row(product, design, carried=False):
        return dict(product=product, design=design, records=[{"carried_from": "old"} if carried else {}])

    original = [
        row("A1", "old", True),
        row("A2", "old", True),
        row("B1", "old", True),
        row("C1", "new"),
        row("C2", "new"),
        row("D1", "new"),
    ]
    snapshot = copy.deepcopy(original)
    assert workflow._evidence_order(original, "inherit-v1") is original
    result = workflow._evidence_order(original, OPTION)
    assert result == original == snapshot and result is not original
    # Distinct design rotation then family rotation, with exact group references.
    original.extend([row("E1", "other"), row("E2", "other")])
    result = workflow._evidence_order(original, OPTION)
    assert [r["product"] for r in result] == ["A1", "A2", "B1", "C1", "E1", "D1", "E2", "C2"]
    assert all(any(r is o for o in original) for r in result)
    assert [r["product"] for r in original][-5:] == ["C1", "C2", "D1", "E1", "E2"]


def test_actual_workflow_deadline_keeps_rotated_remaining_routes_unknown(tmp_path, assets, monkeypatch):
    config = dict(configuration(tmp_path, assets), evidence_allocation_policy=OPTION)
    evaluator, calls = engine(monkeypatch)
    clock = [time.monotonic()]
    deadline = clock[0] + 50
    original = evaluator.evaluate

    def evaluate(*args, **kwargs):
        result = original(*args, **kwargs)
        clock[0] = deadline + 1
        return result

    monkeypatch.setattr(evaluator, "evaluate", evaluate)
    monkeypatch.setattr(workflow, "time", SimpleNamespace(monotonic=lambda: clock[0]))
    out = tmp_path / "deadline"
    assert (
        workflow.run(HIT, [FIRST, LATER], config, out, scoring_mode="2d", deadline_monotonic=deadline) == []
    )
    records = read(out, "routes.json")["routes"]
    assert len(calls) == 1 and records[0]["selection_status"] == "not_scored_deadline"
    assert all(r["selection_status"] == "not_evaluated_deadline" for r in records[1:])


@pytest.mark.parametrize("value", [None, True, {}, "balanced-v1"])
def test_invalid_evidence_policy_fails_before_assets(value):
    with pytest.raises(ValueError, match="evidence_allocation_policy"):
        workflow._config(
            {
                "schema_version": 1,
                "catalog": "missing",
                "rules": "missing",
                "provider": {"kind": "saved", "path": "missing"},
                "evidence_allocation_policy": value,
            }
        )


@pytest.mark.parametrize(
    "name",
    [
        profiles.SELECTIVE_PROFILE,
        profiles.GUARDED_SELECTIVE_PROFILE,
        profiles.DIAGNOSTIC_SELECTIVE_PROFILE,
        profiles.LANDMARK_SELECTIVE_PROFILE,
    ],
)
def test_explicit_evidence_option_is_material_with_labelled_synthetic_runtime(monkeypatch, name):
    artifact = profiles.load_profile(name)
    current = copy.deepcopy(artifact["provenance"]["runtime"])
    current["source_sha256"] = copy.deepcopy(artifact["provenance"]["source_sha256"])
    monkeypatch.setattr(profiles, "numerical_runtime", lambda *args: current)
    config = dict(
        profiles.profile_defaults(artifact),
        recommendation_profile=name,
        evidence_allocation_policy="inherit-v1",
    )
    config["refinement"]["max_total_conformers"] = 24
    assert profiles.profile_provenance(config)["status"] == "engineering_policy_matched"
    config["evidence_allocation_policy"] = OPTION
    result = profiles.profile_provenance(config)
    assert result["status"] == "custom_unvalidated"
    assert result["material_overrides"]["evidence_allocation_policy"] == {
        "profile": "inherit-v1",
        "effective": OPTION,
    }
    assert result["evidence_allocation_policy"] == OPTION
