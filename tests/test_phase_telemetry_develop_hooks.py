"""Tiny inert native/workflow hook controls; model/catalog/3D work is forbidden."""

import json

import pandas as pd
import pytest
from rdkit import Chem

import chembounce as cb
from benchmarks.research.qualified_adapter import native_designs
from chembounce_synthesis import conformers, phase_telemetry as pt, profiles, workflow


def native_fixture(monkeypatch, path, telemetry, *, fail=None):
    path.mkdir()
    calls = []
    monkeypatch.setattr(cb, "_read_replacements", lambda *a, **kw: pd.Series([1.0], index=["c1ccncc1"]))
    monkeypatch.setattr(cb.chemistry, "prepare_replacement", lambda *a: None)

    def generate(*args, **kwargs):
        calls.append(args[3])
        if fail is not None:
            raise fail
        return [(0.0, s) for s in ("CC(F)c1ccncc1", "COC1CCCC1", "NC1CCCCC1")]

    monkeypatch.setattr(cb, "replace_molecule", generate)
    report = {}
    cb.make_scaffold_hopping(
        "CCc1ccccc1",
        Chem.MolFromSmiles("CCc1ccccc1"),
        [Chem.MolFromSmiles("c1ccccc1")],
        tanimoto_threshold=0,
        overall_max_n=2,
        frag_max_n=2,
        scaffold_top_n=1,
        cand_max_n__rplc=3,
        _merge_structure_top_n_=8,
        output_dir=path,
        candidate_thresholds={},
        calculate_shape=False,
        diagnostics=report,
        search_policy="diverse",
        exploration_fraction=0.5,
        search_candidate_budget=3,
        **({"telemetry": telemetry} if telemetry is not None else {}),
    )
    return calls, report


def test_native_hook_parity_exact_outputs_generation_and_no_added_checks(tmp_path, monkeypatch):
    for name in ("ensemble", "score_candidate"):
        monkeypatch.setattr(conformers, name, lambda *a, **k: pytest.fail("3D work forbidden"))
    first, a = native_fixture(monkeypatch, tmp_path / "off", None)
    timing = pt.PhaseTelemetry("native")
    with timing.span("native.worker"):
        with timing.span("native.api"):
            second, b = native_fixture(monkeypatch, tmp_path / "on", timing)
    assert first == second == [3]
    for name in ("overall_result.txt", "fragment_result_0001.txt", "search_decisions.json"):
        assert (tmp_path / "off" / name).read_bytes() == (tmp_path / "on" / name).read_bytes()
    assert a["candidate_selection"] == b["candidate_selection"]
    result = timing.snapshot()
    assert result["status"] == "closed" and not result["interrupted"]
    for name in ("setup", "generation", "allocation", "output"):
        assert result["buckets"][f"native.{name}"]["completed"] == 1
    assert result["buckets"]["native.reference_preflight"]["status"] == "not_reached"
    assert result["buckets"]["native.reference_preflight"]["inclusive_wall_ns"] is None
    assert not any("timing" in p.name or "telemetry" in p.name for p in (tmp_path / "on").iterdir())


def test_native_error_is_interrupted_and_retains_exception(tmp_path, monkeypatch):
    failure = TimeoutError("synthetic owner boundary")
    timing = pt.PhaseTelemetry("native")
    with pytest.raises(TimeoutError) as caught:
        with timing.span("native.api"):
            native_fixture(monkeypatch, tmp_path / "out", timing, fail=failure)
    assert caught.value is failure
    result = timing.snapshot()["buckets"]
    assert result["native.generation"]["status"] == "interrupted"
    assert result["native.allocation"]["status"] == "not_reached"
    assert result["native.output"]["status"] == "not_reached"


def test_native_disabled_no_new_clock_reads(tmp_path, monkeypatch):
    monkeypatch.setattr(pt.PhaseTelemetry, "__init__", lambda *a, **k: pytest.fail("new context/proc reads"))
    monkeypatch.setattr(pt.time, "monotonic_ns", lambda: pytest.fail("new clock"))
    monkeypatch.setattr(pt.time, "process_time_ns", lambda: pytest.fail("new clock"))
    native_fixture(monkeypatch, tmp_path / "off", None)


def test_adapter_explicit_context_forwarding_and_no_json_callback(tmp_path, monkeypatch):
    calls = []

    def native(hit, **kwargs):
        calls.append(kwargs)
        return pd.DataFrame({"Final structure": []}), {}

    monkeypatch.setattr(cb, "chembounce", native)
    request = {
        "query": {"smiles": "CCO"},
        "seed": 42,
        "output_dir": str(tmp_path),
        "assets": {},
        "settings": {"generator_options": {}},
    }
    native_designs(request)
    assert "telemetry" not in calls[-1]
    timing = pt.PhaseTelemetry("native")
    native_designs(request, telemetry=timing)
    assert calls[-1]["telemetry"] is timing
    assert request["settings"]["generator_options"] == {}
    request["settings"]["generator_options"]["telemetry"] = {}
    with pytest.raises(ValueError, match="process-local"):
        native_designs(request)


def inert_workflow(monkeypatch, path, telemetry=None, *, late=False, fail=None):
    # Tiny explicit strings only; no catalog, provider, embedding or scoring.
    for name in ("Catalog", "Rules", "_provider", "_evidence_evaluator", "score_candidate", "ensemble"):
        monkeypatch.setattr(workflow, name, lambda *a, **k: pytest.fail("scientific work forbidden"))

    def prepare(smiles, *args):
        if smiles == "CCN":
            if fail is not None:
                raise fail
            raise ValueError("controlled reference rejection")
        return None

    monkeypatch.setattr(workflow, "prepare_reference", prepare)
    config = {
        "schema_version": 1,
        "catalog": "unused",
        "rules": "unused",
        "provider": {"kind": "saved", "path": "unused"},
    }
    if late:
        config["max_design_attempts"] = 1
    rows = workflow.run(
        "CCO", ["CCN"], config, path, telemetry=telemetry, **({"deadline_monotonic": -1.0} if late else {})
    )
    return rows, json.loads((path / "summary.json").read_text())


def test_workflow_empty_admission_parity_assets_unreached_and_finalization(tmp_path, monkeypatch):
    a, first = inert_workflow(monkeypatch, tmp_path / "off")
    timing = pt.PhaseTelemetry("common")
    with timing.span("common.worker"):
        with timing.span("common.workflow"):
            b, second = inert_workflow(monkeypatch, tmp_path / "on", timing)
    assert a == b == []
    first.pop("elapsed_seconds")
    second.pop("elapsed_seconds")
    # Existing per-run measurement is intentionally not deterministic.
    for summary in (first, second):
        summary["identity_comparison_export"].pop("selected_comparison_serialization_seconds")
    assert first == second
    result = timing.snapshot()
    assert result["status"] == "closed" and not result["interrupted"]
    for name in ("reference_setup", "reference_admission", "workflow_finalization"):
        assert result["buckets"][f"common.{name}"]["status"] == "complete"
    for name in ("catalog", "rules", "evidence_prepare"):
        assert result["buckets"][f"common.{name}"]["status"] == "not_reached"
    assert result["buckets"]["common.refinement"]["status"] == "skipped_by_policy"
    assert not any("telemetry" in p.name for p in (tmp_path / "on").iterdir())


def test_workflow_exception_closes_current_scope_and_finalization_still_runs(tmp_path, monkeypatch):
    failure = RuntimeError("controlled reference dependency failure")
    timing = pt.PhaseTelemetry("common")
    with pytest.raises(RuntimeError) as caught:
        with timing.span("common.workflow"):
            inert_workflow(monkeypatch, tmp_path / "out", timing, fail=failure)
    assert caught.value is failure
    result = timing.snapshot()["buckets"]
    assert result["common.reference_admission"]["status"] == "interrupted"
    assert result["common.catalog"]["status"] == "not_reached"
    assert result["common.workflow_finalization"]["status"] == "complete"


def test_workflow_disabled_new_clock_reads_forbidden(tmp_path, monkeypatch):
    monkeypatch.setattr(pt.PhaseTelemetry, "__init__", lambda *a, **k: pytest.fail("new context/proc reads"))
    monkeypatch.setattr(pt.time, "monotonic_ns", lambda: pytest.fail("new clock"))
    monkeypatch.setattr(pt.time, "process_time_ns", lambda: pytest.fail("new clock"))
    inert_workflow(monkeypatch, tmp_path / "out")


def test_primary_v7_source_binding_and_diagnostic_mismatch_remain_honest():
    config = profiles.apply_profile_defaults({"recommendation_profile": profiles.TRAVERSAL_SELECTIVE_PROFILE})
    current = profiles.numerical_runtime(profiles.TRAVERSAL_SELECTIVE_PROFILE)
    artifact = profiles.load_profile(profiles.TRAVERSAL_SELECTIVE_PROFILE)
    assert current["source_sha256"] == artifact["provenance"]["source_sha256"]
    ordinary = profiles.profile_provenance(config)
    assert not ordinary["material_overrides"]
    assert ordinary["status"] == (
        "custom_unvalidated" if ordinary["runtime_source_mismatches"] else "engineering_policy_matched"
    )
    diagnostic = profiles.profile_provenance(config | {"rejected_pose_diagnostics": True})
    assert diagnostic["status"] == "custom_unvalidated"
    mismatches = diagnostic["optional_diagnostics"]["rejected_pose_diagnostics"]["source_mismatches"]
    assert "chembounce_synthesis/workflow.py" in mismatches


def test_existing_workflow_deadline_is_interrupted_not_completed_scientific_work(tmp_path, monkeypatch):
    timing = pt.PhaseTelemetry("common")
    with timing.span("common.workflow"):
        rows, summary = inert_workflow(monkeypatch, tmp_path / "deadline", timing, late=True)
    assert rows == [] and summary["deadline_reached"]
    result = timing.snapshot()["buckets"]
    assert result["common.reference_admission"]["status"] == "interrupted"
    assert result["common.route_search"]["status"] == "interrupted"
    assert result["common.catalog"]["status"] == "not_reached"
    # Core publication can finish despite an incomplete scientific phase.
    assert result["common.workflow_finalization"]["status"] == "complete"


@pytest.mark.parametrize("policy", ["legacy", "diverse"])
def test_public_native_propagation_and_sequential_repeated_scopes(tmp_path, monkeypatch, policy):
    def fragments(**kwargs):
        kwargs["diagnostics"]["database"] = {"path": None}
        return [], 1, 1, 1, 1, 1

    monkeypatch.setattr(cb, "get_frags_cands", fragments)
    monkeypatch.setattr(cb, "replace_molecule", lambda *a, **k: pytest.fail("generation forbidden"))
    monkeypatch.setattr(cb.utils, "calc_electron_shape", lambda *a, **k: pytest.fail("3D forbidden"))
    timing = pt.PhaseTelemetry("native")
    with timing.span("native.worker"):
        with timing.span("native.api"):
            frame, _ = cb.chembounce(
                "CCO",
                output_dir=tmp_path,
                candidate_thresholds={},
                calculate_shape=False,
                search_policy=policy,
                search_candidate_budget=1,
                telemetry=timing,
            )
    assert frame.empty
    result = timing.snapshot()["buckets"]
    assert result["native.output"]["status"] == "complete"
    assert result["native.output"]["calls"] == (2 if policy == "diverse" else 1)
    assert result["native.setup"]["calls"] == (2 if policy == "diverse" else 1)
    assert result["native.generation"]["status"] == ("complete" if policy == "diverse" else "unavailable")
    assert result["native.reference_preflight"]["status"] == "not_reached"
