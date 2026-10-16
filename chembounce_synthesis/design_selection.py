"""Bounded raw-design intake and deterministic reference eligibility/refill."""

from collections.abc import Sized

from .common import canonical, positive


def selection_policy(config):
    quota = config.get("max_designs", 20)
    positive(quota, "max_designs")
    bounded = config.get("recommendation", {}).get("mode") == "qualified" or "max_design_attempts" in config
    attempts = config.get("max_design_attempts", 4 * quota) if bounded else None
    if bounded:
        positive(attempts, "max_design_attempts")
    return {
        "version": "bounded-refill-v1" if bounded else "legacy-capped-v1",
        "max_eligible_designs": quota,
        "max_raw_attempts": attempts,
        "default_attempt_budget": "4 * max_designs; engineering resource bound" if bounded else None,
        "stereo_policy": "Preserve explicit states; never assign or enumerate unspecified stereochemistry",
        "provider_sampling": "Original per-design, per-sample identity/order; rejected entries cause no provider calls",
    }


def sized_count(designs):
    return len(designs) if isinstance(designs, Sized) else None


def select_designs(designs, policy, prepare, deadline_reached):
    """Consume at most the raw bound, with no exhaustion lookahead or retries.

    ``prepare`` performs only reference preflight, never embedding or scoring.
    Unsized/unconsumed tails remain explicitly unknown and are not reconstructible.
    """
    total = sized_count(designs)
    iterator = iter(designs)
    checks, raw, unique, eligible = [], [], [], []
    seen = {}
    stop = "raw_attempt_limit"
    exhausted = None
    for index in range(policy["max_raw_attempts"]):
        if deadline_reached():
            stop = "deadline"
            break
        if len(eligible) >= policy["max_eligible_designs"]:
            stop = "eligible_limit"
            break
        try:
            value = next(iterator)
        except StopIteration:
            stop, exhausted = "input_exhausted", True
            break
        # Invalid non-string values retain their type, without coercing them to
        # valid chemistry or requiring arbitrary objects to be JSON serializable.
        saved = value if isinstance(value, str) else None
        raw.append(saved)
        check = {"raw_index": index, "input": saved, "input_type": type(value).__name__, "role": "design"}
        checks.append(check)
        try:
            design = canonical(value)
        except ValueError as exc:
            check.update(
                status="rejected", stage="design_parse", reason_code="invalid_design", reason=str(exc)
            )
            continue
        check["design"] = design
        if design in seen:
            check.update(
                status="duplicate", duplicate_of_raw_index=seen[design], reason_code="exact_state_duplicate"
            )
            continue
        seen[design] = index
        unique.append(design)
        try:
            prepare(design)
        except ValueError as exc:
            check.update(
                status="rejected",
                stage="design_reference_preflight",
                reason_code="unspecified_stereochemistry"
                if "specified stereochemistry" in str(exc)
                else "reference_ineligible",
                reason=str(exc),
            )
            continue
        check["status"] = "eligible"
        eligible.append(design)
    if total is not None:
        exhausted = len(raw) == total
        if exhausted:
            stop = "input_exhausted"
    elif len(eligible) >= policy["max_eligible_designs"] and stop == "raw_attempt_limit":
        stop = "eligible_limit"
    counts = {
        "raw_supplied": total,
        "raw_consumed": len(raw),
        "parse_attempted": len(raw),
        "preflight_attempted": len(unique),
        "invalid": sum(c.get("stage") == "design_parse" for c in checks),
        "duplicates": sum(c["status"] == "duplicate" for c in checks),
        "preflight_rejected": sum(c.get("stage") == "design_reference_preflight" for c in checks),
        "eligible": len(eligible),
        "unattempted": total - len(raw) if total is not None else None,
    }
    limits = []
    if exhausted is not True:
        if len(raw) == policy["max_raw_attempts"]:
            limits.append("max_design_attempts")
        if len(eligible) == policy["max_eligible_designs"]:
            limits.append("max_designs")
        if stop == "deadline":
            limits.append("deadline_monotonic")
    return {
        "policy": policy,
        "requested_designs": raw,
        "canonical_unique_requested_designs": unique,
        "capped_designs": unique,
        "eligible_designs": eligible,
        "checks": checks,
        "counts": counts,
        "stop_reason": stop,
        "input_exhausted": exhausted,
        "limits_reached": limits,
        "unconsumed_tail": "Not read or canonicalized; cannot be reconstructed by continuation",
        "scope": "Bounded raw attempts including invalid inputs and exact-state duplicates; stop at eligible quota; no stereo assignment, embedding or provider retries",
    }
