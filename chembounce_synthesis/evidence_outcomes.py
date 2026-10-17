"""Recorded forward-support diagnostics, never an acceptance or chemistry gate.

Pure standard library so saved reports need no model or chemistry dependencies.
Missing historical cause codes remain missing; reason text is never reverse parsed.
"""

from collections import Counter
import json

VERSION = "forward-evidence-outcomes-v1"
CAUSES = {
    "matched": {"exact_full_product_match"},
    "mismatched": {"completed_top_three_nonmatch"},
    "partially_evaluated": {"incomplete_or_invalid_hypotheses"},
    "unknown": {
        "predictor_unavailable",
        "prediction_deadline_exhausted",
        "predictor_preparation_failed",
        "predictor_call_budget_exhausted",
        "input_token_budget_exhausted",
        "output_token_budget_exhausted",
    },
    "failed": {"prediction_execution_failed", "prediction_validation_failed"},
}
SCOPE = "Bounded computational support only; neither physical synthesis failure nor experimental validation"


def prediction_outcome(prediction):
    """Describe existing forward status without manufacturing a historical cause."""
    status = prediction.get("status")
    code = prediction.get("reason_code", "not_recorded")
    if "reason_code" in prediction and (not isinstance(code, str) or code not in CAUSES.get(status, set())):
        raise ValueError("Recorded forward reason_code is malformed or contradicts prediction status")
    outcome = {
        "matched": "supported",
        "mismatched": "completed_nonmatch",
        "partially_evaluated": "incomplete_prediction",
        "unknown": "unassessed",
        "failed": "execution_error",
    }.get(status, "unrecorded")
    if status == "mismatched" and code == "not_recorded":
        outcome = "nonmatch_recorded"
    if status == "unknown" and code == "predictor_preparation_failed":
        outcome = "execution_error"
    return {
        "status": outcome,
        "prediction_status": status if status is not None else "not_recorded",
        "cause_code": code,
        "reason": prediction.get("reason"),
    }


def route_outcome(steps):
    per_step = [
        {"node_id": step.get("node_id"), **prediction_outcome(step.get("forward_prediction", {}))}
        for step in steps
    ]
    statuses = {item["status"] for item in per_step}
    return {
        "version": VERSION,
        "status": next(iter(statuses)) if len(statuses) == 1 else "mixed" if statuses else "unrecorded",
        "step_outcomes": per_step,
        "scope": SCOPE,
    }


def unsupported_reason(prediction):
    """Truthful wording for an already-failed exact forward-support requirement."""
    status = prediction.get("status")
    if status == "mismatched":
        if prediction.get("reason_code") == "completed_top_three_nonmatch":
            return "Completed bounded top-three predictions do not match the exact full product"
        return "A forward nonmatch was recorded; completion and operational cause were not recorded"
    if status == "partially_evaluated":
        return "Forward prediction is incomplete; exact full-product support was not established"
    outcome = prediction_outcome(prediction)
    label = "error" if outcome["status"] == "execution_error" else "unassessed"
    detail = prediction.get("reason") or "No operational reason recorded"
    return f"Forward prediction {label} ({outcome['cause_code']}): {detail}"


def fallback_reasons(record):
    """Preserve singular and plural causes; use selection status only as fallback."""
    reasons = []
    singular = record.get("reason")
    if isinstance(singular, str) and singular:
        reasons.append(singular)
    plural = record.get("reasons", [])
    if isinstance(plural, list):
        for reason in plural:
            if isinstance(reason, str) and reason and reason not in reasons:
                reasons.append(reason)
    return reasons or [record.get("selection_status", "not_evaluated")]


def summarize(records):
    """Exclusive record counts at each level; route and step counts are separate."""
    route_counts, step_counts, causes = Counter(), Counter(), Counter()
    with_evidence = 0
    for record in records:
        if "synthesis_evidence" not in record:
            continue
        with_evidence += 1
        outcome = route_outcome(record["synthesis_evidence"].get("steps", []))
        route_counts[outcome["status"]] += 1
        for step in outcome["step_outcomes"]:
            step_counts[step["status"]] += 1
            causes[step["cause_code"]] += 1
    return {
        "version": VERSION,
        "route_records_with_evidence": with_evidence,
        "route_records_without_evidence": len(records) - with_evidence,
        "route_record_outcomes": dict(sorted(route_counts.items())),
        "step_records": sum(step_counts.values()),
        "step_record_outcomes": dict(sorted(step_counts.items())),
        "step_record_causes": dict(sorted(causes.items())),
        "counting_scope": "Repeated route/step records, not unique products or predictor calls; mixed routes retain every step outcome",
        "scope": SCOPE,
    }


def _same(left, right):
    return json.dumps(left, sort_keys=True, allow_nan=False) == json.dumps(
        right, sort_keys=True, allow_nan=False
    )


def validate_recorded(evidence):
    """Check optional saved diagnostic projections; absence is legacy, not proof."""
    steps = evidence.get("steps", [])
    for step in steps:
        expected = prediction_outcome(step.get("forward_prediction", {}))
        if "forward_outcome" in step and not _same(step["forward_outcome"], expected):
            raise ValueError("Recorded step forward outcome differs from its prediction status/cause")
    if "forward_outcome" in evidence and not _same(evidence["forward_outcome"], route_outcome(steps)):
        raise ValueError("Recorded route forward outcome differs from its steps")


def validate_summary(summary, records):
    if "forward_evidence_outcomes" in summary and not _same(
        summary["forward_evidence_outcomes"], summarize(records)
    ):
        raise ValueError("Recorded forward outcome counts differ from route/step records")


def validate_records(records, decisions=()):
    """Check all retained projections, including duplicated qualification decisions."""
    by_id = {record.get("route_id"): record for record in records}

    def qualification_projection(decision, record):
        if "forward_outcome" not in decision:
            return
        if record is None or "synthesis_evidence" not in record:
            raise ValueError("Recorded qualification forward outcome has no route evidence")
        expected = route_outcome(record["synthesis_evidence"].get("steps", []))
        if not _same(decision["forward_outcome"], expected):
            raise ValueError("Recorded qualification forward outcome differs from route predictions")

    for record in records:
        if "synthesis_evidence" in record:
            validate_recorded(record["synthesis_evidence"])
        qualification_projection(record.get("qualification", {}), record)
    for decision in decisions:
        qualification_projection(decision, by_id.get(decision.get("route_id")))
