"""Pure recorded-data diagnostics; no chemistry or model calls."""

import copy

import pytest

from chembounce_synthesis import evidence_outcomes as e


def step(status, code=None, node="r0", reason=None):
    prediction = {"status": status}
    if code is not None:
        prediction["reason_code"] = code
    if reason is not None:
        prediction["reason"] = reason
    return {"node_id": node, "forward_prediction": prediction}


def test_mixed_route_and_record_counts_keep_denominators_and_unknowns():
    supported = step("matched", "exact_full_product_match")
    unknown = step("unknown", "input_token_budget_exhausted", "r1", "Forward prediction budget exhausted")
    failed = step("failed", "prediction_execution_failed", reason="controlled failure")
    partial = step("partially_evaluated", "incomplete_or_invalid_hypotheses")
    nonmatch = step("mismatched", "completed_top_three_nonmatch")
    records = [
        {"synthesis_evidence": {"steps": [supported, unknown]}},
        {"synthesis_evidence": {"steps": [supported]}},
        {"synthesis_evidence": {"steps": [failed]}},
        {"synthesis_evidence": {"steps": [partial]}},
        {"synthesis_evidence": {"steps": [nonmatch]}},
        {},
    ]
    value = e.summarize(records)
    assert value["route_records_with_evidence"] == sum(value["route_record_outcomes"].values()) == 5
    assert value["route_records_without_evidence"] == 1
    assert value["step_records"] == sum(value["step_record_outcomes"].values()) == 6
    assert sum(value["step_record_causes"].values()) == 6
    assert value["route_record_outcomes"]["mixed"] == 1
    assert value["step_record_outcomes"]["supported"] == 2
    assert value["step_record_outcomes"]["unassessed"] == 1
    assert "not unique products or predictor calls" in value["counting_scope"]
    e.validate_summary({"forward_evidence_outcomes": value}, records)
    value["step_records"] = 6.0
    with pytest.raises(ValueError, match="counts differ"):
        e.validate_summary({"forward_evidence_outcomes": value}, records)


def test_legacy_cause_absence_is_not_inferred_from_text():
    legacy = step("unknown", reason="Forward prediction budget exhausted")
    outcome = e.prediction_outcome(legacy["forward_prediction"])
    assert outcome["status"] == "unassessed" and outcome["cause_code"] == "not_recorded"
    e.validate_recorded({"steps": [legacy]})
    e.validate_summary({}, [{"synthesis_evidence": {"steps": [legacy]}}])
    assert e.route_outcome([])["status"] == "unrecorded"


@pytest.mark.parametrize(
    "status,code",
    [
        ("matched", "predictor_call_budget_exhausted"),
        ("unknown", "exact_full_product_match"),
        ("mismatched", "incomplete_or_invalid_hypotheses"),
        ("failed", "prediction_deadline_exhausted"),
        ("matched", "future-code"),
        ("matched", False),
        ("matched", 0),
        ("matched", []),
    ],
)
def test_present_contradictory_or_malformed_cause_refused_even_without_projection(status, code):
    item = step(status)
    item["forward_prediction"]["reason_code"] = code
    with pytest.raises(ValueError, match="reason_code"):
        e.validate_recorded({"steps": [item]})


def test_projection_cannot_relabel_incomplete_as_supported():
    item = step("partially_evaluated", "incomplete_or_invalid_hypotheses")
    item["forward_outcome"] = e.prediction_outcome(item["forward_prediction"])
    evidence = {"steps": [item], "forward_outcome": e.route_outcome([item])}
    e.validate_recorded(evidence)
    for level in ("step", "route"):
        forged = copy.deepcopy(evidence)
        target = forged["steps"][0] if level == "step" else forged
        target["forward_outcome"]["status"] = "supported"
        with pytest.raises(ValueError, match="outcome differs"):
            e.validate_recorded(forged)


def test_fallback_preserves_plural_reasons_and_singular_context_without_duplicates():
    assert e.fallback_reasons(
        {"reason": "error", "reasons": ["error", "budget"], "selection_status": "evidence_rejected"}
    ) == ["error", "budget"]
    assert e.fallback_reasons({"selection_status": "not_scored_budget"}) == ["not_scored_budget"]


def test_legacy_nonmatch_does_not_imply_complete_strict_top_three():
    prediction = {
        "status": "mismatched",
        "reason": "No completed exact full-product hypothesis in bounded top three",
    }
    value = e.prediction_outcome(prediction)
    assert value["status"] == "nonmatch_recorded"
    assert value["cause_code"] == "not_recorded"
    assert "completion and operational cause were not recorded" in e.unsupported_reason(prediction)
    records = [{"synthesis_evidence": {"steps": [{"forward_prediction": prediction}]}}]
    summary = e.summarize(records)
    assert summary["step_record_outcomes"] == {"nonmatch_recorded": 1}
    assert summary["step_record_causes"] == {"not_recorded": 1}
    prediction["reason_code"] = "completed_top_three_nonmatch"
    assert e.prediction_outcome(prediction)["status"] == "completed_nonmatch"


@pytest.mark.parametrize("location", ["route", "decision", "both"])
def test_qualification_projection_rederived_from_predictions(location):
    item = step("unknown", "predictor_call_budget_exhausted")
    record = {"route_id": "unselected", "synthesis_evidence": {"steps": [item]}, "qualification": {}}
    decision = {"route_id": "unselected"}
    e.validate_records([record], [decision])  # Legacy absence is valid.
    expected = e.route_outcome([item])
    record["qualification"]["forward_outcome"] = copy.deepcopy(expected)
    decision["forward_outcome"] = copy.deepcopy(expected)
    e.validate_records([record], [decision])
    if location in ("route", "both"):
        record["qualification"]["forward_outcome"]["status"] = "supported"
    if location in ("decision", "both"):
        decision["forward_outcome"]["status"] = "supported"
    with pytest.raises(ValueError, match="qualification forward outcome"):
        e.validate_records([record], [decision])


@pytest.mark.parametrize("route_present", [False, True])
def test_qualification_projection_without_source_evidence_refused(route_present):
    record = {"route_id": "r0"}
    decision = {"route_id": "r0", "forward_outcome": e.route_outcome([])}
    with pytest.raises(ValueError, match="no route evidence"):
        e.validate_records([record] if route_present else [], [decision])
