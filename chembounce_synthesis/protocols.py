"""Versioned, explicit computational protocols with no assumed activity cutoff."""

import copy

PROTOCOLS = ("legacy", "research-v1")
CORE_DEFAULTS = {"search_policy": "legacy", "exploration_fraction": 0.25, "search_candidate_budget": 4096}
RESEARCH_REFINEMENT = {
    "top_n": 10,
    "num_conformers": 8,
    "seeds": [42, 43, 44],
    "threshold_policy": "all_seeds",
    "ranking_statistic": "minimum",
    "max_score_range": 0.10,
    "max_metric_range": 0.15,
    "adaptive_max_conformers": 32,
    "conformer_growth_factor": 2,
    "borderline_margin": 0.02,
    "max_total_conformers": 1680,
}


def validate_protocol(name):
    if not isinstance(name, str) or name not in PROTOCOLS:
        raise ValueError(f"protocol must be one of {PROTOCOLS}")


def core_protocol(name="legacy"):
    validate_protocol(name)
    result = dict(CORE_DEFAULTS)
    if name == "research-v1":
        result["search_policy"] = "diverse"
    return result


def apply_synthesis_protocol(config):
    """Explicit config fields override preset defaults; no cutoffs are invented."""
    name = config.get("protocol", "legacy")
    validate_protocol(name)
    result = copy.deepcopy(config)
    if name == "research-v1":
        result.setdefault("allocation_policy", "balanced")
        refinement = result.setdefault("refinement", {})
        if not isinstance(refinement, dict):
            raise ValueError("refinement must be a JSON object")
        for key, value in RESEARCH_REFINEMENT.items():
            refinement.setdefault(key, copy.deepcopy(value))
        result.setdefault("minimum_scores", {})
    return result


def protocol_provenance(name="legacy"):
    validate_protocol(name)
    return {
        "name": name,
        "core_defaults": core_protocol(name),
        "synthesis_refinement_defaults": copy.deepcopy(RESEARCH_REFINEMENT) if name == "research-v1" else {},
        "override_policy": "Explicit API/config/CLI values override preset defaults; effective settings are recorded separately",
        "tolerance_interpretation": "Research spread tolerances and borderline margin are engineering sampling criteria, not calibrated activity, affinity or synthesis-success thresholds",
        "similarity_acceptance": "No hard similarity minimum is assumed; minimum_scores must be specified explicitly when scientifically justified",
        "budget_scope": (
            "1680 refinement candidate-conformer requests permit ten complete 8+16+32 conformer stages over three seeds; reference ensembles and actual CPU time are reported separately, not equated to this budget"
            if name == "research-v1"
            else "Explicit refinement settings determine candidate-conformer request bounds"
        ),
    }
