"""Typed work accounting from synthetic artifacts; no chemistry or model runs."""

import pytest

from benchmarks.research.protocol import write_json
from benchmarks.research.qualified_analysis import COSTS, _cost_report, _validate_costs
from benchmarks.research.qualified_runner import unit_counts


def observed_counts(tmp_path, **artifacts):
    out = tmp_path / "synthesis"
    out.mkdir(exist_ok=True)
    for name, value in artifacts.items():
        write_json(out / (name + ".json"), value)
    counts = unit_counts(tmp_path)
    _validate_costs(
        {
            "wall_seconds": None,
            "cpu_seconds": None,
            "gpu_seconds": None,
            "peak_rss_mb": None,
            "peak_gpu_mb": None,
            "native_units": counts,
            "common_units": {},
        }
    )
    return counts


def test_missing_artifacts_do_not_impute_work(tmp_path):
    assert unit_counts(tmp_path) == {}
    assert observed_counts(tmp_path, summary={"scored_products": 2}) == {"scored_products": 2}


@pytest.mark.parametrize("artifact", [{}, None, [], {"candidates": None}, {"candidates": {}}])
def test_refinement_requires_an_observed_collection(tmp_path, artifact):
    counts = observed_counts(tmp_path, summary={"scored_products": 2}, refinement=artifact)
    assert counts == {"scored_products": 2}


@pytest.mark.parametrize(
    "bad_row",
    [
        {},
        None,
        {"candidate_conformer_requests": True},
        {"candidate_conformer_requests": -1},
        {"candidate_conformer_requests": 1.5},
    ],
)
def test_partial_refinement_does_not_publish_a_partial_sum(tmp_path, bad_row):
    counts = observed_counts(
        tmp_path,
        summary={"refinement_budget": {"limit": 144}},
        refinement={"candidates": [{"candidate_conformer_requests": 24}, bad_row]},
    )
    assert counts == {"refinement_candidate_conformer_limit": 144}


@pytest.mark.parametrize("values", [[], [0], [24, 48]])
def test_observed_refinement_and_proposal_totals_without_summary(tmp_path, values):
    counts = observed_counts(
        tmp_path,
        refinement={"candidates": [{"candidate_conformer_requests": n} for n in values]},
        proposals={"proposals": [{"sampling": {"generated_tokens": n}} for n in values]},
    )
    assert counts == {
        "refinement_candidate_conformer_requests": sum(values),
        "returned_proposal_records": len(values),
        "returned_proposals_with_token_counts": len(values),
        "returned_proposal_generated_tokens": sum(values),
    }


@pytest.mark.parametrize("value", [True, False, -1, 1.5, "3", None, {}, []])
def test_invalid_scalar_counts_are_omitted_without_losing_known_counts(tmp_path, value):
    counts = observed_counts(
        tmp_path,
        summary={
            "scored_products": value,
            "proposals": value,
            "verified_routes": value,
            "designs_used": value,
            "designs_supplied": value,
            "design_selection": {"counts": {"unknown": value, "observed": 2}},
            "synthesis_evidence_usage": {"unknown": value, "observed": 0},
            "refinement_budget": {"limit": value},
        },
    )
    assert counts == {"design_observed": 2, "forward_observed": 0}


@pytest.mark.parametrize("value", [0, 3])
def test_observed_scalar_counts_are_retained(tmp_path, value):
    summary = dict.fromkeys(
        ("scored_products", "proposals", "verified_routes", "designs_used", "designs_supplied"), value
    )
    counts = observed_counts(
        tmp_path,
        summary={
            **summary,
            "design_selection": {"counts": {"observed": value}},
            "synthesis_evidence_usage": {"observed": value},
            "refinement_budget": {"limit": value},
        },
    )
    assert counts == {
        **summary,
        "design_observed": value,
        "forward_observed": value,
        "refinement_candidate_conformer_limit": value,
    }


@pytest.mark.parametrize(
    "artifact", [{}, None, [], {"proposals": None}, {"proposals": {}}, {"proposals": [None]}]
)
def test_proposals_require_an_observed_record_collection(tmp_path, artifact):
    assert observed_counts(tmp_path, proposals=artifact) == {}


@pytest.mark.parametrize(
    "row",
    [
        {},
        {"sampling": None},
        {"sampling": {}},
        {"sampling": {"generated_tokens": True}},
        {"sampling": {"generated_tokens": -1}},
    ],
)
def test_partial_proposal_tokens_keep_known_record_counts(tmp_path, row):
    counts = observed_counts(tmp_path, proposals={"proposals": [{"sampling": {"generated_tokens": 8}}, row]})
    assert counts == {"returned_proposal_records": 2, "returned_proposals_with_token_counts": 1}


@pytest.mark.parametrize("value", [True, -1, None])
def test_search_totals_require_valid_counts_in_every_uncached_report(tmp_path, value):
    counts = observed_counts(
        tmp_path,
        summary={
            "search_reports": [
                {"stage": "expansion", "combinations": 8, "root_state_visits": 3},
                {"stage": "expansion", "combinations": value, "root_state_visits": 2},
                {"stage": "expansion", "cached": True, "combinations": 99, "root_state_visits": 99},
            ]
        },
    )
    assert counts["expansion_search_reports"] == 2
    assert counts["expansion_root_state_visits"] == 5
    assert "expansion_combinations" not in counts
    assert "expansion_beam_truncations" not in counts


def test_missing_search_row_count_is_unknown(tmp_path):
    counts = observed_counts(
        tmp_path,
        summary={"search_reports": [{"stage": "projection", "combinations": 8}, {"stage": "projection"}]},
    )
    assert counts["projection_search_reports"] == 2
    assert "projection_combinations" not in counts


def test_observed_empty_search_reports_are_zero(tmp_path):
    counts = observed_counts(tmp_path, summary={"search_reports": []})
    assert counts["projection_search_reports"] == counts["expansion_search_reports"] == 0
    assert counts["projection_combinations"] == counts["expansion_combinations"] == 0
    assert set(counts.values()) == {0}


@pytest.mark.parametrize("reports", [None, {}, [None], [{}], [{"stage": "unrecognized"}]])
def test_malformed_search_evidence_does_not_impute_zero(tmp_path, reports):
    assert observed_counts(tmp_path, summary={"scored_products": 2, "search_reports": reports}) == {
        "scored_products": 2
    }


@pytest.mark.parametrize("cached", [1, 0, "true", None])
def test_malformed_cache_flag_does_not_discard_or_count_work(tmp_path, cached):
    counts = observed_counts(
        tmp_path,
        summary={
            "search_reports": [
                {"stage": "projection", "combinations": 8},
                {"stage": "expansion", "cached": cached, "combinations": 5},
            ]
        },
    )
    assert counts["projection_combinations"] == 8
    assert not any(key.startswith("expansion_") for key in counts)


@pytest.mark.parametrize("value", [None, [], 3])
def test_malformed_summary_objects_preserve_independent_evidence(tmp_path, value):
    counts = observed_counts(
        tmp_path, summary=value, refinement={"candidates": [{"candidate_conformer_requests": 24}]}
    )
    assert counts == {"refinement_candidate_conformer_requests": 24}
    counts = observed_counts(
        tmp_path,
        summary={
            "scored_products": 2,
            "design_selection": value,
            "synthesis_evidence_usage": value,
            "refinement_budget": value,
        },
    )
    assert counts == {"scored_products": 2, "refinement_candidate_conformer_requests": 24}


def test_analysis_coverage_distinguishes_unknown_from_observed_zero(tmp_path):
    rows = []
    for name, artifacts in (
        ("unknown", {}),
        ("empty", {"refinement": {"candidates": []}}),
        ("positive", {"refinement": {"candidates": [{"candidate_conformer_requests": 24}]}}),
    ):
        directory = tmp_path / name
        directory.mkdir()
        counts = observed_counts(directory, **artifacts)
        rows.append({"costs": {**dict.fromkeys(COSTS), "native_units": counts, "common_units": {}}})
    assert _cost_report(rows)["native_units"]["refinement_candidate_conformer_requests"] == {
        "observed_total": 24,
        "observed_attempts": 2,
        "planned_attempts": 3,
    }
