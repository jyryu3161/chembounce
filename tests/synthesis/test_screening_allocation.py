"""Optional allocation: actual toy workflow/evidence, no model or production assets."""

import copy
import json
import time
from types import SimpleNamespace

import pytest

from chembounce_synthesis import profiles, workflow
from chembounce_synthesis.common import canonical
from chembounce_synthesis.evidence import SynthesisEvidenceEvaluator

HIT = "CC(=O)Nc1ccccc1"
HIGH = "CC(=O)Nc1ccncc1"
MID = "CCC(=O)Nc1ccncc1"
LOW = "CCCC(=O)Nc1ccncc1"
OPTION = profiles.HIT_SCREENING_ALLOCATION


@pytest.fixture
def alignment_api():
    pytest.importorskip("rdkit.Chem.rdShapeAlign", reason="Actual3D workflow requires optional alignment API")


def configuration(tmp_path, assets, *, calls=3, cap=1):
    proposals = []
    for design, acid, product in [
        (LOW, "CCCC(=O)O", LOW),
        (LOW, "CC(=O)O", HIGH),
        (HIGH, "CCCC(=O)O", LOW),
        (HIGH, "CCC(=O)O", MID),
    ]:
        proposals.append(
            dict(
                design=design,
                response={
                    "reactions": [
                        {"reaction_template": "amide", "reactants": [acid, "Nc1ccncc1"], "product": product}
                    ],
                    "building_blocks": [acid, "Nc1ccncc1"],
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
        max_scored_products=cap,
        synthesis_evidence={"policy": {"require_forward_agreement": True}, "max_predictor_calls": calls},
    )


def engine(monkeypatch, calls=3):
    predictions = []

    class Predictor:
        provenance = {"provider": "controlled_allocation_predictor", "uses_conditions": False}

        def token_budget(self, *args):
            return {"input_tokens": 8, "max_output_tokens": 24}

        def predict(self, reactants, *args):
            predictions.append(tuple(reactants))
            acid = next(r for r in reactants if r.endswith("(=O)O"))
            product = {"CC(=O)O": HIGH, "CCC(=O)O": MID, "CCCC(=O)O": LOW}[acid]
            return {"products": [product], "complete": [True], "input_tokens": 8, "output_tokens": 12}

    value = SynthesisEvidenceEvaluator(
        Predictor(), {"require_forward_agreement": True}, max_predictor_calls=calls
    )
    monkeypatch.setattr(workflow, "_evidence_evaluator", lambda config: value)
    return value, predictions


def read(path, name):
    return json.loads((path / name).read_text())


@pytest.mark.parametrize("legacy", ["balanced", "encounter"])
def test_inherit_parity_and_hit_rank_changes_only_screening(
    tmp_path, assets, monkeypatch, legacy, alignment_api
):
    config = configuration(tmp_path, assets, calls=2)
    config["allocation_policy"] = legacy
    snapshots = []
    for name, policy in [("omitted", None), ("inherit", "inherit-v1"), ("ranked", OPTION)]:
        current = copy.deepcopy(config)
        if policy is not None:
            current["screening_allocation_policy"] = policy
        evaluator, predictions = engine(monkeypatch, 2)
        out = tmp_path / name
        rows = workflow.run(HIT, [LOW, HIGH], current, out)
        records = read(out, "routes.json")["routes"]
        summary = read(out, "summary.json")
        snapshots.append((rows, records, summary, predictions, evaluator.usage))
    assert snapshots[0][2]["screening_allocation"] == snapshots[1][2]["screening_allocation"]
    assert (
        [r["smiles"] for r in snapshots[0][0]] == [r["smiles"] for r in snapshots[1][0]] == [canonical(LOW)]
    )
    assert [r["smiles"] for r in snapshots[2][0]] == [canonical(HIGH)]
    assert snapshots[0][3] == snapshots[1][3] == snapshots[2][3]
    assert [s[4]["predictor_calls"] for s in snapshots] == [2, 2, 2]
    for _, records, summary, _, _ in snapshots:
        assert len(records) == 4
        assert [r.get("eligible_for_scoring", False) for r in records] == [True, True, True, False]
        assert records[3]["selection_status"] == "evidence_rejected"
        assert "budget" in records[3]["synthesis_evidence"]["steps"][0]["forward_prediction"]["reason"]
        assert summary["scored_products"] == summary["screening_allocation_counts"]["scorer_invocations"] == 1
    ranked = snapshots[2][2]["screening_allocation_counts"]
    assert ranked["eligible_distinct_products"] == 2 and ranked["deferred_design_groups"] == 1
    assert ranked["deferred_design_attempts"] == 0 and ranked["distinct_products_attempted"] == 1
    assert (
        snapshots[2][1][0]["selection_status"]
        == snapshots[2][1][2]["selection_status"]
        == "not_scored_budget"
    )
    manifest = read(tmp_path / "ranked", "manifest.json")
    assert manifest["effective_workflow"]["allocation_policy"] == legacy
    assert manifest["effective_workflow"]["screening_allocation_policy"] == OPTION


def test_deferred_design_is_independently_scored_after_distinct_products(
    tmp_path, assets, monkeypatch, alignment_api
):
    config = configuration(tmp_path, assets, cap=4)
    config["screening_allocation_policy"] = OPTION
    engine(monkeypatch)
    original = workflow.score_candidate
    calls = []

    def score(product, hit, design, *args, **kwargs):
        calls.append((product, design))
        if product == canonical(LOW) and design == canonical(LOW):
            raise ValueError("controlled design-specific failure")
        return original(product, hit, design, *args, **kwargs)

    monkeypatch.setattr(workflow, "score_candidate", score)
    out = tmp_path / "fallback"
    rows = workflow.run(HIT, [LOW, HIGH], config, out)
    assert calls == [
        (canonical(HIGH), canonical(LOW)),
        (canonical(MID), canonical(HIGH)),
        (canonical(LOW), canonical(LOW)),
        (canonical(LOW), canonical(HIGH)),
    ]
    records = read(out, "routes.json")["routes"]
    assert records[0]["selection_status"] == "3d_failed" and records[2]["selection_status"] == "accepted"
    selected = next(r for r in rows if r["smiles"] == canonical(LOW))
    assert selected["design"] == canonical(HIGH) and selected["route_id"] == records[2]["route_id"]
    assert records[0]["route_id"] not in selected["alternative_route_ids"]
    summary = read(out, "summary.json")
    assert summary["screening_allocation_counts"]["distinct_products_attempted"] == 3
    assert summary["screening_allocation_counts"]["deferred_design_attempts"] == 1
    assert summary["screening_allocation_counts"]["scorer_invocations"] == 4
    assert summary["screening_allocation"][-1]["allocation_pass"] == "deferred_design"


def test_ties_follow_first_evidence_eligible_route_not_balanced_or_rejected_order():
    records = [{"route_id": f"r{i}", "eligible_for_scoring": i != 0} for i in range(4)]

    def entry(product, design, items):
        return dict(product=product, design=design, records=items, tanimoto={"hit_tanimoto": 0.5})

    a = entry("A", "first", [records[0], records[2]])
    b = entry("B", "first", [records[1]])
    alt = entry("B", "second", [records[3]])
    original = [alt, a, b]
    inherited, deferred = workflow._screening_order(original, "inherit-v1", records)
    assert inherited is original and not deferred
    ordered, deferred = workflow._screening_order(original, OPTION, records)
    assert ordered == [b, a, alt] and deferred == {("B", "second")}
    assert records[0]["eligible_for_scoring"] is False


@pytest.mark.parametrize("policy", ["inherit-v1", OPTION])
def test_deadline_does_not_promote_unevaluated_or_unscored_routes(
    tmp_path, assets, monkeypatch, policy, alignment_api
):
    config = configuration(tmp_path, assets, cap=4)
    config["screening_allocation_policy"] = policy
    evaluator, _ = engine(monkeypatch)
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
    assert workflow.run(HIT, [LOW, HIGH], config, out, deadline_monotonic=deadline) == []
    records = read(out, "routes.json")["routes"]
    assert records[0]["selection_status"] == "not_scored_deadline"
    assert all(r["selection_status"] == "not_evaluated_deadline" for r in records[1:])
    summary = read(out, "summary.json")
    assert summary["deadline_reached"] and summary["scored_products"] == 0
    assert summary["screening_allocation_counts"]["scorer_invocations"] == 0


def test_carried_priority_inherited_and_ranked_routes_stay_auditable(
    tmp_path, assets, monkeypatch, alignment_api
):
    config = configuration(tmp_path, assets)
    engine(monkeypatch)
    first = tmp_path / "first"
    assert workflow.run(HIT, [LOW, HIGH], config, first)[0]["smiles"] == canonical(LOW)
    config["carry_forward"] = [str(first)]
    for name, policy in [("inherited", "inherit-v1"), ("ranked", OPTION)]:
        engine(monkeypatch)
        current = dict(config, screening_allocation_policy=policy)
        out = tmp_path / name
        rows = workflow.run(HIT, [LOW, HIGH], current, out)
        assert rows[0]["smiles"] == canonical(LOW if name == "inherited" else HIGH)
        carried = [r for r in read(out, "routes.json")["routes"] if r.get("carried_from")]
        assert carried and all(r["synthesis_evidence"]["eligibility"]["passed"] for r in carried)
        if name == "ranked":
            carried_low = [r for r in carried if r["verification"]["product"] == canonical(LOW)]
            assert carried_low and all(r["selection_status"] == "not_scored_budget" for r in carried_low)
            assert not any(r.get("qualification", {}).get("qualified") for r in carried_low)


@pytest.mark.parametrize(
    "name",
    [
        profiles.SELECTIVE_PROFILE,
        profiles.GUARDED_SELECTIVE_PROFILE,
        profiles.DIAGNOSTIC_SELECTIVE_PROFILE,
        profiles.LANDMARK_SELECTIVE_PROFILE,
    ],
)
def test_explicit_option_is_material_even_with_synthetic_matching_runtime(monkeypatch, name):
    artifact = profiles.load_profile(name)
    current = {
        **copy.deepcopy(artifact["provenance"]["runtime"]),
        "source_sha256": copy.deepcopy(artifact["provenance"]["source_sha256"]),
    }
    monkeypatch.setattr(profiles, "numerical_runtime", lambda *args: current)
    config = dict(profiles.profile_defaults(artifact), recommendation_profile=name)
    config["refinement"]["max_total_conformers"] = 24
    assert profiles.profile_provenance(config)["status"] == "engineering_policy_matched"
    config["screening_allocation_policy"] = "inherit-v1"
    assert profiles.profile_provenance(config)["status"] == "engineering_policy_matched"
    config["screening_allocation_policy"] = OPTION
    result = profiles.profile_provenance(config)
    assert result["status"] == "custom_unvalidated"
    assert result["material_overrides"]["screening_allocation_policy"] == {
        "profile": "inherit-v1",
        "effective": OPTION,
    }
    assert result["screening_allocation_policy"] == OPTION


@pytest.mark.parametrize("value", [None, True, {}, "balanced-v1", "invented"])
def test_invalid_policy_rejected_before_workflow_assets(tmp_path, assets, value):
    config = configuration(tmp_path, assets)
    config["screening_allocation_policy"] = value
    with pytest.raises(ValueError, match="screening_allocation_policy"):
        workflow._config(config)


def test_selective_preflight_failure_counts_attempt_not_scorer_call(
    tmp_path, assets, monkeypatch, alignment_api
):
    from chembounce_synthesis.similarity_constraints import BoundSimilarity, MINIMUM_SCORES, VERSION

    config = configuration(tmp_path, assets)
    config.update(
        screening_allocation_policy=OPTION,
        recommendation={"mode": "qualified"},
        conformers={"num_conformers": 8, "seed": 42, "scoring_objective": "original-hit-v1"},
        refinement={"top_n": 1, "num_conformers": 8, "seeds": [42, 43, 44], "max_total_conformers": 24},
        minimum_scores=dict(MINIMUM_SCORES),
        similarity_constraints={"version": VERSION},
    )
    engine(monkeypatch)
    screening_started = []
    original_order = workflow._screening_order
    original_preflight = BoundSimilarity.preflight_candidate

    def order(*args):
        screening_started.append(True)
        return original_order(*args)

    def reject(self, product):
        # Earlier route-capacity admission also invokes this method. Activate
        # only after the actual evidence queue reaches its screening boundary.
        if screening_started:
            raise ValueError("controlled screening preflight failure")
        return original_preflight(self, product)

    monkeypatch.setattr(workflow, "_screening_order", order)
    monkeypatch.setattr(BoundSimilarity, "preflight_candidate", reject)
    monkeypatch.setattr(
        workflow.selective_scoring, "score_candidate", lambda *a, **k: pytest.fail("no scorer call")
    )
    out = tmp_path / "preflight"
    assert workflow.run(HIT, [LOW, HIGH], config, out) == []
    summary = read(out, "summary.json")
    assert summary["scored_products"] == summary["screening_allocation_counts"]["screening_attempts"] == 1
    assert summary["screening_allocation_counts"]["scorer_invocations"] == 0
    failed = [r for r in read(out, "routes.json")["routes"] if r["selection_status"] == "3d_failed"]
    assert failed and all(
        r["similarity_failure"]["ensemble_request_status"] == "none_requested" for r in failed
    )


def test_actual_packaged_v5_preserves_historical_source_mismatch():
    # Actual installed research runtime/profile; no synthetic source/runtime patch.
    artifact = profiles.load_profile(profiles.LANDMARK_SELECTIVE_PROFILE)
    current = profiles.numerical_runtime(profiles.LANDMARK_SELECTIVE_PROFILE)
    assert {
        path
        for path, digest in current["source_sha256"].items()
        if digest != artifact["provenance"]["source_sha256"][path]
    } == {"chembounce_synthesis/common.py"}
    config = dict(
        profiles.profile_defaults(artifact), recommendation_profile=profiles.LANDMARK_SELECTIVE_PROFILE
    )
    config["refinement"]["max_total_conformers"] = 24
    provenance = profiles.profile_provenance(config)
    diagnostic = provenance["optional_diagnostics"]["rejected_pose_diagnostics"]
    assert "chembounce_synthesis/workflow.py" in diagnostic["source_mismatches"]
    assert not diagnostic["enabled"]
    assert provenance["status"] == "custom_unvalidated"
    assert "chembounce_synthesis/common.py" in provenance["runtime_source_mismatches"]
    # The historical source mismatch never conceals a material allocation override.
    assert "screening_allocation_policy" not in provenance["material_overrides"]
    config["screening_allocation_policy"] = OPTION
    assert (
        profiles.profile_provenance(config)["material_overrides"]["screening_allocation_policy"]["effective"]
        == OPTION
    )


def test_experimental_policy_refuses_2d_before_loading_assets(tmp_path, assets, monkeypatch):
    config = configuration(tmp_path, assets)
    config["screening_allocation_policy"] = OPTION
    monkeypatch.setattr(workflow, "Catalog", lambda *a, **k: pytest.fail("must refuse before catalog load"))
    with pytest.raises(ValueError, match="requires3d"):
        workflow.run(HIT, [LOW, HIGH], config, tmp_path / "two-dimensional", scoring_mode="2d")
