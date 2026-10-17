"""Synthetic source-only five-arm orchestration; no model or tokenizer loads."""

import json
from types import SimpleNamespace

import pytest

# Test-only artifacts retain real loaders, policy checks and collector validation.
pytestmark = pytest.mark.usefixtures("synthetic_benchmark_profile_bindings")

from benchmarks.research import qualified_adapter as adapter, qualified_runner as runner
from benchmarks.research.protocol import write_json


def request(tmp_path, arm="full", seed=97):
    return {
        "query": {"id": "q", "target_id": "t", "smiles": "CC(=O)Nc1ccccc1"},
        "seed": seed,
        "settings": {"arm": arm},
        "assets": {},
        "output_dir": str(tmp_path),
        "deadline_monotonic": 999999999999,
    }


def config_file(tmp_path, profile="qualified-development-v1"):
    config = {
        "schema_version": 1,
        "recommendation_profile": profile,
        "catalog": "catalog.tsv",
        "rules": "rules.json",
        "provider": {
            "kind": "synllama",
            "model_path": "unused-model",
            "prompt_path": "unused-prompt",
            "samples": 2,
            "max_new_tokens": 1200,
            "seed": 42,
        },
        "max_design_attempts": 12,
        "max_designs": 3,
        "max_results": 3,
        "max_scored_products": 12,
        "refinement": {"top_n": 6, "max_total_conformers": 144},
        "search": {
            "building_block_top_k": 3,
            "beam_width": 6,
            "max_combinations": 256,
            "max_steps": 12,
            "reaction_product_limit": 64,
        },
        "expansion": {"rounds": 1, "max_routes": 12, "building_block_top_k": 6, "min_tanimoto": 0.0},
        "synthesis_evidence": {
            "profile": "computational-route-support-v1",
            "predictor": {"kind": "reactiont5", "model_path": "unused-forward"},
            "max_predictor_calls": 24,
            "max_input_tokens": 12288,
            "max_output_tokens": 18432,
        },
    }
    path = tmp_path / "config.json"
    write_json(path, config)
    from chembounce_synthesis.workflow import load_config
    from chembounce_synthesis.profiles import profile_provenance

    assert profile_provenance(load_config(path))["status"] == (
        "engineering_policy_matched"
        if adapter.selective_benchmark_profile(profile)
        else "development_protocol_matched"
    )
    return path


@pytest.mark.parametrize("arm,samples", [("full", 2), ("crem", 2), ("mol2mol", 2), ("synllama", 6)])
def test_generation_seed_never_overwrites_fixed_final_protocol(tmp_path, arm, samples):
    r = request(tmp_path, arm)
    r["assets"]["synthesis_config"] = str(config_file(tmp_path))
    cfg = adapter.common_config(r)
    assert cfg["provider"]["seed"] == 97 and cfg["provider"]["samples"] == samples
    assert cfg["conformers"]["seed"] == 42 and cfg["conformers"]["num_conformers"] == 8
    assert cfg["refinement"]["seeds"] == [42, 43, 44] and cfg["refinement"]["max_total_conformers"] == 144


def test_ablation_only_changes_search_or_raw_limit(tmp_path):
    r = request(tmp_path)
    r["assets"]["synthesis_config"] = str(config_file(tmp_path))
    original = adapter.common_config(r)
    for ablation, key, value in [
        ("no_extra_refill_capacity", "max_design_attempts", 3),
        ("design_only_route_search", "route_search_policy", "design-v1"),
    ]:
        r["settings"]["ablation"] = ablation
        changed = adapter.common_config(r)
        assert changed[key] == value
        changed[key] = original[key]
        assert changed == original


def test_baseline_bridge_excludes_unselected_routes_and_models(tmp_path):
    native, output = tmp_path / "native", tmp_path / "common"
    (native / "synthesis").mkdir(parents=True)
    output.mkdir()
    route = {"root": "r", "nodes": [{"id": "r", "product": "CC(=O)Nc1ccncc1"}]}
    row = {"route_id": "selected", "design": "CC(=O)Nc1ccncc1", "smiles": "CC(=O)Nc1ccncc1"}
    (native / "synthesis/results.tsv").write_text("\t".join(row) + "\n" + "\t".join(row.values()) + "\n")
    records = [
        {**row, "route": route, "verification": {"product": row["smiles"]}},
        {**row, "route_id": "unselected", "route": {"do_not_import": True}},
    ]
    write_json(native / "synthesis/routes.json", {"routes": records})
    response = {"candidates": [{"smiles": row["smiles"], "design": row["design"], "route": route}]}
    r = request(output, "baseline")
    designs, bridge, counts = adapter.baseline_bridge(native, response, r, output)
    assert designs == [row["design"]] and counts["raw_count"] == 1
    assert [r["route_id"] for r in adapter.read(bridge / "routes.json")["routes"]] == ["selected"]
    r["assets"]["synthesis_config"] = str(config_file(tmp_path))
    cfg = adapter.common_config(r, bridge)
    assert cfg["provider"]["kind"] == "saved" and cfg["expansion"]["rounds"] == 0
    assert adapter.read(cfg["provider"]["path"])["proposals"] == []
    response["candidates"][0]["route"] = {"different": True}
    with pytest.raises(ValueError, match="primary graph"):
        adapter.baseline_bridge(native, response, r, output)


def test_crem_complete_raw_pool_reselects_beyond_native_three(tmp_path):
    native = tmp_path / "native"
    native.mkdir()
    r = request(tmp_path, "crem")
    write_json(native / "request.json", {"query": r["query"], "seed": 97})
    raw = ["invalid", "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1", "CC(=O)Nc1cccnc1"]
    (native / "response.raw.jsonl").write_text(
        "".join(json.dumps({"index": i, "raw_smiles": s}) + "\n" for i, s in enumerate(raw))
    )
    write_json(
        native / "response.json",
        {"candidates": [], "evaluations": 4, "metadata": {"completed_within_budget": True}},
    )
    r["native"] = {"directory": str(native), "files_sha256": runner.file_manifest(native)}
    selected, bridge, counts = adapter.native_supply(r, tmp_path)
    assert len(selected) == 2 and bridge is None and counts["raw_count"] == 4
    (native / "response.raw.jsonl").write_text("")
    with pytest.raises(ValueError, match="changed"):
        adapter.native_supply(r, tmp_path)


def protocol():
    return {
        "schema": "qualified-benchmark-protocol-v1",
        "collector_catalog_policy": runner.CATALOG_POLICY,
        "run_id": "synthetic",
        "budget": runner.BUDGET,
        "source_sha256": {},
        "queries": [{"id": f"q{i}", "target_id": f"t{i}", "smiles": "CC"} for i in range(2)],
        "seeds": [97],
        "methods": [
            {
                "id": arm,
                "arm": arm,
                "resolved_assets": {},
                "native_method": {
                    "kind": "external",
                    "command": ["unused"],
                    "asset_names": {},
                    "settings": {},
                },
            }
            for arm in ["full", "baseline", "crem", "mol2mol", "synllama"]
        ],
    }


def test_plan_blocks_complete_and_rotates():
    p = runner.make_plan(protocol(), "a" * 64)
    assert len(p["slot_order"]) == 10
    assert [s["method_id"] for s in p["slot_order"][:5]] == p["methods"]
    assert [s["method_id"] for s in p["slot_order"][5:]] == p["methods"][1:] + p["methods"][:1]


def test_two_phase_budget_failure_does_not_rescue_partial(tmp_path, monkeypatch):
    p = protocol()
    plan = runner.make_plan(p, "a" * 64)
    calls = []

    def bounded(method, req, out, **kwargs):
        calls.append(req)
        out.mkdir()
        write_json(out / "request.json", req)
        return {
            "status": "timeout",
            "response": {},
            "resources": {"sampled_cpu_seconds_lower_bound": 7, "sampled_peak_rss_mb": 200},
        }

    monkeypatch.setattr(runner, "bounded_method", bounded)
    result = runner.attempt(
        p["methods"][0], p["queries"][0], 97, tmp_path / "attempt", plan, SimpleNamespace(signal_number=None)
    )
    assert len(calls) == 1 and calls[0]["budget"]["wall_seconds"] == 360
    assert result["status"] == "timeout" and result["qualified"] == [] and result["costs"]["cpu_seconds"] == 7


def test_common_budget_no_unused_native_transfer(tmp_path, monkeypatch):
    p = protocol()
    plan = runner.make_plan(p, "a" * 64)
    calls = []

    def bounded(method, req, out, **kwargs):
        calls.append(req)
        out.mkdir()
        write_json(out / "request.json", req)
        write_json(out / "response.json", {"candidates": [], "evaluations": 0})
        return {
            "status": "no_candidates",
            "response": {"candidates": [], "evaluations": 0},
            "resources": {"sampled_cpu_seconds_lower_bound": 2, "sampled_peak_rss_mb": 100},
        }

    monkeypatch.setattr(runner, "bounded_method", bounded)
    monkeypatch.setattr(runner.qualified_collector, "checked_products", lambda *args, **kwargs: [])
    result = runner.attempt(
        p["methods"][0], p["queries"][0], 97, tmp_path / "attempt", plan, SimpleNamespace(signal_number=None)
    )
    assert len(calls) == 2 and calls[1]["budget"]["wall_seconds"] <= 600
    assert result["status"] == "no_candidates" and result["costs"]["cpu_seconds"] == 4


def test_source_drift_fails_before_execution(tmp_path):
    p = protocol()
    with pytest.raises(ValueError, match="source pins"):
        runner.resolve_protocol(p, tmp_path)


def test_baseline_request_kind_is_native_pipeline(tmp_path):
    method = protocol()["methods"][1]
    method["native_method"]["request_kind"] = "pipeline"
    r = runner.request_for(method, protocol()["queries"][0], 97, tmp_path, 360, "native", {})
    assert r["kind"] == "pipeline" and r["seed"] == 97


@pytest.mark.parametrize("artifact", ["summary", "refinement", "proposals"])
@pytest.mark.parametrize("response", [{}, {"candidates": [], "evaluations": 2}])
def test_partial_cost_json_preserves_failed_attempt_and_measured_resources(
    tmp_path, monkeypatch, artifact, response
):
    p = protocol()
    plan = runner.make_plan(p, "a" * 64)

    def bounded(method, req, out, **kwargs):
        (out / "synthesis").mkdir(parents=True)
        write_json(out / "synthesis/summary.json", {"scored_products": 1})
        (out / "synthesis" / (artifact + ".json")).write_text('{"truncated":')
        return {
            "status": "timeout",
            "response": response,
            "resources": {"sampled_cpu_seconds_lower_bound": 7, "sampled_peak_rss_mb": 200},
        }

    monkeypatch.setattr(runner, "bounded_method", bounded)
    result = runner.attempt(
        p["methods"][0], p["queries"][0], 97, tmp_path / "attempt", plan, SimpleNamespace(signal_number=None)
    )
    assert result["status"] == "timeout" and result["costs"]["cpu_seconds"] == 7
    assert result["costs"]["peak_rss_mb"] == 200
    assert result["costs"]["native_units"] == ({"returned_evaluations": 2} if response else {})
    assert result["failures"][-1]["stage"] == "native_cost_accounting"
    assert adapter.read(tmp_path / "attempt/collection-attempt.json") == result


def test_actual_native_cli_obeys_supply_quota_and_seed(tmp_path):
    import subprocess
    import sys

    req = request(tmp_path)
    req["assets"] = {"scaffold_db": "fixture-only"}
    write_json(tmp_path / "request.json", req)
    code = """
import sys, types
from pathlib import Path
module = types.ModuleType("chembounce")
def generate(hit, **options):
    assert options["overall_max_n"] == 12 and options["seed"] == 97
    assert options["scaffold_db_file"] == "fixture-only"
    return {"Final structure": ["CC", "CO"]}, {"fixture_calls": 1}
module.chembounce = generate
sys.modules["chembounce"] = module
from benchmarks.research.qualified_adapter import main
main()
"""
    subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            "--phase",
            "native",
            "--request",
            str(tmp_path / "request.json"),
            "--response",
            str(tmp_path / "response.json"),
        ],
        check=True,
    )
    assert len(adapter.read(tmp_path / "response.json")["candidates"]) == 2


def test_runner_executes_full_rotated_grid_and_checkpoints(tmp_path, monkeypatch):
    p = protocol()
    path = tmp_path / "protocol.json"
    write_json(path, p)
    monkeypatch.setattr(runner, "resolve_protocol", lambda raw, base: p)
    monkeypatch.setattr(runner, "preload_catalogs", lambda *args: {m["id"]: None for m in p["methods"]})
    called = []

    def attempted(method, query, seed, output, plan, termination, **kwargs):
        called.append((method["id"], query["id"], seed))
        return {
            "method_id": method["id"],
            "target_id": query["target_id"],
            "query_id": query["id"],
            "query_sha256": runner.canonical_sha256(query),
            "seed": seed,
            "status": "no_candidates",
            "qualified": [],
            "failures": [{"stage": "native", "reason": "synthetic empty"}],
            "costs": {
                "wall_seconds": 1,
                "cpu_seconds": None,
                "gpu_seconds": None,
                "peak_rss_mb": None,
                "peak_gpu_mb": None,
                "native_units": {},
                "common_units": {},
            },
        }

    monkeypatch.setattr(runner, "attempt", attempted)
    result = runner.run(path, tmp_path / "run")
    plan = adapter.read(tmp_path / "run/plan.json")
    assert called == [(s["method_id"], s["query_id"], s["seed"]) for s in plan["slot_order"]]
    assert len(result["attempts"]) == 10
    assert adapter.read(tmp_path / "run/analysis.json")["planned_attempts"] == 10
    assert result["measurement_policy"]["gpu_seconds"] == "not measured; null"


@pytest.mark.parametrize(
    "key,value",
    [
        ("minimum_scores", {"hit_tanimoto": 0.9, "hit_shape": 0.9, "hit_pharmacophore": 0.9}),
        ("ranking_policy", "joint-v1"),
        ("route_search_policy", "design-v1"),
    ],
)
def test_undeclared_scientific_overrides_fail_before_work(tmp_path, key, value):
    path = config_file(tmp_path)
    config = adapter.read(path)
    config[key] = value
    write_json(path, config)
    r = request(tmp_path)
    r["assets"]["synthesis_config"] = str(path)
    with pytest.raises(ValueError, match="frozen calibration"):
        adapter.common_config(r)


def test_mol2mol_native_selection_preserves_nll_order_and_failed_raw(tmp_path):
    from chembounce_synthesis.common import canonical
    from benchmarks.research.protocol import digest

    native = tmp_path / "native"
    native.mkdir()
    r = request(tmp_path, "mol2mol")
    write_json(native / "request.json", {"query": r["query"], "seed": 97})
    raw = {
        "input_smiles": canonical(r["query"]["smiles"]),
        "seed": 97,
        "rows": [
            {"smiles": s, "native_nll": n, "truncated": truncated, "terminated": not truncated}
            for s, n, truncated in [
                ("CC(=O)Nc1ccncc1", 2.0, False),
                ("invalid", -10.0, True),
                ("CC(=O)Nc1cccnc1", 1.0, False),
            ]
        ],
    }
    write_json(native / "mol2mol-raw.json", raw)
    write_json(
        native / "response.json",
        {"candidates": [], "evaluations": 3, "metadata": {"raw_sha256": digest(native / "mol2mol-raw.json")}},
    )
    r["native"] = {"directory": str(native), "files_sha256": runner.file_manifest(native)}
    designs, bridge, counts = adapter.native_supply(r, tmp_path)
    assert designs == ["CC(=O)Nc1cccnc1", "CC(=O)Nc1ccncc1"] and counts["raw_count"] == 3
    decisions = adapter.read(tmp_path / "native-selection.json")["decisions"]
    assert decisions[1]["selection_status"] == "unterminated_at_model_limit"


def test_typed_search_work_does_not_double_count_cached_expansions(tmp_path):
    out = tmp_path / "synthesis"
    out.mkdir()
    write_json(
        out / "summary.json",
        {
            "search_reports": [
                {"stage": "projection", "combinations": 4, "root_state_visits": 2},
                {"stage": "expansion", "combinations": 5, "root_state_visits": 3},
                {"stage": "expansion", "cached": True, "combinations": 5, "root_state_visits": 3},
            ],
            "refinement_budget": {"limit": 504},
        },
    )
    write_json(out / "refinement.json", {"candidates": [{"candidate_conformer_requests": 168}]})
    write_json(out / "proposals.json", {"proposals": [{"sampling": {"generated_tokens": 99}}]})
    counts = runner.unit_counts(tmp_path)
    assert counts["projection_combinations"] == 4 and counts["expansion_combinations"] == 5
    assert counts["expansion_root_state_visits"] == 3 and "projection_beam_truncations" not in counts
    assert (
        counts["refinement_candidate_conformer_requests"] == 168
        and counts["refinement_candidate_conformer_limit"] == 504
    )
    assert counts["returned_proposal_generated_tokens"] == 99
    write_json(out / "proposals.json", {"proposals": [{"sampling": {}}]})
    assert "returned_proposal_generated_tokens" not in runner.unit_counts(tmp_path)


def test_common_cleanup_overrun_preserves_costs_and_continues_next_slot(tmp_path, monkeypatch):
    p = protocol()
    p["queries"] = p["queries"][:1]
    p["methods"] = [{"id": name, "arm": "synllama", "resolved_assets": {}} for name in ("first", "next")]
    path = tmp_path / "protocol.json"
    write_json(path, p)
    monkeypatch.setattr(runner, "resolve_protocol", lambda raw, base: p)
    monkeypatch.setattr(runner, "preload_catalogs", lambda *args: {m["id"]: None for m in p["methods"]})
    clock, requests, collected = SimpleNamespace(now=1000.0), [], []
    monkeypatch.setattr(runner.time, "monotonic", lambda: clock.now)

    def bounded(method, req, out, **kwargs):
        requests.append(req)
        out.mkdir()
        write_json(out / "request.json", req)
        duration = req["budget"]["wall_seconds"] + 1 if len(requests) == 1 else 10
        clock.now += duration
        return {
            "status": "no_candidates",
            "response": {"candidates": [], "evaluations": 0},
            "resources": {
                "wall_seconds": duration,
                "sampled_cpu_seconds_lower_bound": 5,
                "sampled_peak_rss_mb": 100,
            },
        }

    def collect(*args, **kwargs):
        collected.append(args[0])
        return []

    monkeypatch.setattr(runner, "bounded_method", bounded)
    monkeypatch.setattr(runner.qualified_collector, "checked_products", collect)
    result = runner.run(path, tmp_path / "run")
    first, second = result["attempts"]
    assert len(requests) == 2 and len(collected) == 1
    assert first["status"] == "timeout" and first["qualified"] == []
    assert first["failures"][0]["stage"] == "common" and "including cleanup" in first["failures"][0]["reason"]
    assert first["costs"]["wall_seconds"] == 601 and first["costs"]["cpu_seconds"] == 5
    assert second["status"] == "no_candidates"
    assert adapter.read(tmp_path / "run/attempt-0000/phases.json")["phases"]["common"]["allocation_exceeded"]


@pytest.fixture(params=["v2", "v3"])
def profile_binding(request):
    return selective_binding(request.param)


def selective_binding(version="v2"):
    from chembounce_synthesis.profiles import SELECTIVE_PROFILE, LANDMARK_SELECTIVE_PROFILE, PROFILE_SHA256

    name = {"v2": SELECTIVE_PROFILE, "v3": LANDMARK_SELECTIVE_PROFILE}[version]

    return {
        "schema": "qualified-benchmark-protocol-" + version,
        "recommendation_profile": name,
        "profile_artifact_sha256": PROFILE_SHA256[name],
    }


def test_selective_protocol_binding_reaches_common_and_v3_full_native(
    tmp_path, profile_binding, benchmark_metadata_guard
):
    p = protocol() | profile_binding
    plan = runner.make_plan(p, "a" * 64)
    assert plan["benchmark_protocol"] == profile_binding
    assert len(plan["slot_order"]) == 10
    for phase in ("native", "common"):
        req = runner.request_for(p["methods"][0], p["queries"][0], 97, tmp_path, 360, phase, plan)
        assert req["seed"] == 97
        bound = phase == "common" or profile_binding["schema"].endswith("v3")
        assert ("benchmark_protocol" in req) == bound
        if bound:
            assert req["benchmark_protocol"] == profile_binding


@pytest.mark.parametrize("mutation", ["missing_name", "wrong_digest", "v1_upgrade", "wrong_profile"])
def test_selective_protocol_rejects_missing_or_mismatched_binding(
    tmp_path, mutation, profile_binding, benchmark_metadata_guard
):
    p = protocol() | profile_binding
    if mutation == "missing_name":
        p.pop("recommendation_profile")
    elif mutation == "wrong_digest":
        p["profile_artifact_sha256"] = "0" * 64
    elif mutation == "v1_upgrade":
        p["schema"] = "qualified-benchmark-protocol-v1"
    else:
        p["recommendation_profile"] = "qualified-development-v1"
    with pytest.raises(ValueError, match="binding|protocol v2"):
        runner.resolve_protocol(p, tmp_path)


def test_selective_common_config_requires_exact_profile_and_preserves_anchors(tmp_path, profile_binding):
    r = request(tmp_path)
    r["assets"]["synthesis_config"] = str(config_file(tmp_path, profile_binding["recommendation_profile"]))
    with pytest.raises(ValueError, match="protocol-bound"):
        adapter.common_config(r)
    r["benchmark_protocol"] = profile_binding
    cfg = adapter.common_config(r)
    assert cfg["provider"]["seed"] == 97
    assert cfg["conformers"]["seed"] == 42 and cfg["refinement"]["seeds"] == [42, 43, 44]
    assert cfg["minimum_scores"] == {"hit_tanimoto": 0.35, "hit_shape": 0.70, "hit_pharmacophore": 0.50}
    assert cfg["similarity_constraints"]["polar_recall"] == 0.75
    config = adapter.read(r["assets"]["synthesis_config"])
    config["core_smiles"] = "CC(=O)N"
    landmark = {"id": "oxygen", "family": "Acceptor", "hit_atom_ids": [2], "max_distance": 1.0}
    config["similarity_constraints"] = {"landmarks": [landmark]}
    write_json(r["assets"]["synthesis_config"], config)
    anchored = adapter.common_config(r)
    assert anchored["core_smiles"] == "CC(=O)N"
    assert anchored["similarity_constraints"]["landmarks"] == [landmark]
    config["similarity_constraints"]["polar_recall"] = 0.5
    write_json(r["assets"]["synthesis_config"], config)
    with pytest.raises(ValueError, match="frozen engineering profile"):
        adapter.common_config(r)


def test_selective_runner_failure_keeps_denominator_and_fixed_budgets(
    tmp_path, monkeypatch, profile_binding, benchmark_metadata_guard
):
    from benchmarks.research.qualified_analysis import analyze

    p = protocol() | profile_binding
    plan = runner.make_plan(p, "a" * 64)
    calls = []

    def failed(method, req, out, **kwargs):
        calls.append(req)
        out.mkdir()
        write_json(out / "request.json", req)
        return {
            "status": "error",
            "response": {},
            "resources": {"sampled_cpu_seconds_lower_bound": 0, "sampled_peak_rss_mb": 0},
        }

    monkeypatch.setattr(runner, "bounded_method", failed)
    result = runner.attempt(
        p["methods"][-1], p["queries"][0], 97, tmp_path / "failed", plan, SimpleNamespace(signal_number=None)
    )
    assert len(calls) == 1 and calls[0]["benchmark_protocol"] == profile_binding
    assert calls[0]["budget"]["wall_seconds"] <= 600
    assert result["status"] == "error" and result["qualified"] == []
    report = analyze(
        plan,
        {
            "schema": "qualified-benchmark-collection-v1",
            "plan_sha256": runner.canonical_sha256(plan),
            "collector_source_sha256": plan["collector_source_sha256"],
            "attempts": [result],
        },
    )
    assert report["planned_attempts"] == 10


@pytest.fixture
def benchmark_metadata_guard(monkeypatch):
    """Native modules are imported, but these tests must perform no chemistry or inference."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from chembounce_synthesis import conformers, selective_scoring, workflow
    from chembounce_synthesis.catalog import Catalog
    from chembounce_synthesis.reactions import Rules
    from chembounce_synthesis.forward_support import PinnedReactionT5Forward

    def forbidden(*args, **kwargs):
        pytest.fail("Metadata-only test attempted chemistry, model/catalog loading or a workflow run")

    # Supported minimum RDKit may omit alignment; metadata checks still run.
    if conformers.rdShapeAlign is not None:
        monkeypatch.setattr(conformers.rdShapeAlign, "AlignMol", forbidden)
    for target, names in (
        (Chem, ("MolFromSmiles", "SDMolSupplier")),
        (AllChem, ("EmbedMolecule", "EmbedMultipleConfs")),
        (conformers, ("ensemble", "score_candidate", "_overlay")),
        (selective_scoring, ("ensemble", "score_candidate", "_overlay", "restore_result", "validate_result")),
        (workflow, ("run", "score_candidate")),
        (Catalog, ("__init__",)),
        (Rules, ("__init__",)),
        (PinnedReactionT5Forward, ("__init__",)),
        (adapter, ("native_designs", "native_supply")),
        (runner, ("bounded_method",)),
    ):
        for name in names:
            monkeypatch.setattr(target, name, forbidden)


@pytest.mark.parametrize("attack", ["hash", "name", "schema", "extra", "missing", "cross_version"])
def test_benchmark_metadata_exact_binding(profile_binding, attack, benchmark_metadata_guard):
    binding = dict(profile_binding)
    assert adapter.benchmark_profile(binding) == binding["recommendation_profile"]
    if attack == "hash":
        binding["profile_artifact_sha256"] = "0" * 64
    elif attack == "name":
        binding["recommendation_profile"] = "qualified-selective-engineering-v4"
    elif attack == "schema":
        binding["schema"] = "qualified-benchmark-protocol-v4"
    elif attack == "extra":
        binding["trust_me"] = True
    elif attack == "missing":
        binding.pop("profile_artifact_sha256")
    else:
        binding["schema"] = "qualified-benchmark-protocol-" + (
            "v3" if binding["schema"].endswith("v2") else "v2"
        )
    with pytest.raises(ValueError, match="binding"):
        adapter.benchmark_profile(binding)


def test_benchmark_metadata_default_cannot_upgrade_legacy(monkeypatch, benchmark_metadata_guard):
    from chembounce_synthesis import profiles

    monkeypatch.setattr(profiles, "DEFAULT_RECOMMENDATION_PROFILE", "future-profile")
    assert adapter.benchmark_profile() == profiles.PROFILE
    assert adapter.selective_benchmark_profile(profiles.PROFILE) is False
    for name in (profiles.SELECTIVE_PROFILE, profiles.LANDMARK_SELECTIVE_PROFILE):
        assert adapter.selective_benchmark_profile(name) is True
    for name in (profiles.GUARDED_SELECTIVE_PROFILE, profiles.DIAGNOSTIC_SELECTIVE_PROFILE, "future-profile"):
        with pytest.raises(ValueError, match="Unsupported"):
            adapter.selective_benchmark_profile(name)
    assert "benchmark_protocol" not in runner.make_plan(protocol(), "a" * 64)


@pytest.mark.parametrize("arm", ["full", "synllama"])
def test_benchmark_metadata_resolve_exact_config(tmp_path, profile_binding, benchmark_metadata_guard, arm):
    """Use tiny inert asset bytes: no catalog or model is constructed."""
    from pathlib import Path

    config = config_file(tmp_path, profile_binding["recommendation_profile"])
    assets = {"synthesis_config": {"artifact": str(config), "sha256": runner.digest(config)}}
    for name in ("catalog.tsv", "rules.json", "unused-model", "unused-prompt", "unused-forward"):
        path = tmp_path / name
        path.write_text("inert binding fixture, not a usable scientific asset\n")
        alias = {"catalog.tsv": "catalog", "rules.json": "rules"}.get(name, name)
        assets[alias] = {"artifact": str(path), "sha256": runner.digest(path)}
    if arm == "full":
        scaffold = tmp_path / "inert-scaffolds.txt"
        scaffold.write_text("inert binding fixture; never loaded")
        assets["scaffold_db"] = {"artifact": str(scaffold), "sha256": runner.digest(scaffold)}
    root = Path(runner.__file__).resolve().parents[2]
    sources = (
        list((root / "chembounce_synthesis").glob("*.py"))
        + [
            root / "benchmarks/research" / name
            for name in (
                "qualified_runner.py",
                "qualified_collector.py",
                "qualified_adapter.py",
                "qualified_analysis.py",
                "runner.py",
            )
        ]
        + [root / "chembounce.py"]
    )
    p = protocol() | profile_binding
    p["source_sha256"] = {str(path): runner.digest(path) for path in sources}
    p["methods"] = [{"id": arm, "arm": arm, "assets": assets}]
    resolved = runner.resolve_protocol(p, tmp_path)
    assert resolved["methods"][0]["resolved_assets"]["synthesis_config"] == str(config)
    assert runner.make_plan(resolved, "a" * 64)["benchmark_protocol"] == profile_binding
    helper = str(root / "chembounce_synthesis/landmark_proposals.py")
    expected = p["source_sha256"].pop(helper)
    with pytest.raises(ValueError, match="required current source pins"):
        runner.resolve_protocol(p, tmp_path)
    p["source_sha256"][helper] = "0" * 64
    with pytest.raises(ValueError, match="Protocol source changed"):
        runner.resolve_protocol(p, tmp_path)
    p["source_sha256"][helper] = expected
    # A separately rehashed config still cannot switch the bound profile.
    config_file(tmp_path)
    p["methods"][0]["assets"]["synthesis_config"]["sha256"] = runner.digest(config)
    with pytest.raises(ValueError, match="Method config profile"):
        runner.resolve_protocol(p, tmp_path)


@pytest.mark.parametrize("tamper", [False, True])
def test_benchmark_metadata_runner_to_collector_binding(
    tmp_path, monkeypatch, profile_binding, benchmark_metadata_guard, tamper
):
    """Real request/envelope/config path, fake execution and no scientific artifact validation."""
    p = protocol() | profile_binding
    config = config_file(tmp_path, profile_binding["recommendation_profile"])
    p["methods"][0]["resolved_assets"] = {
        "synthesis_config": str(config),
        "scaffold_db": str(tmp_path / "inert-scaffolds.txt"),
    }
    plan, requests, collected = runner.make_plan(p, "a" * 64), [], []

    def bounded(method, req, out, **kwargs):
        requests.append(req)
        out.mkdir()
        if tamper and len(requests) == 2 and "benchmark_protocol" in req:
            req["benchmark_protocol"] = {**req["benchmark_protocol"], "profile_artifact_sha256": "0" * 64}
        write_json(out / "request.json", req)
        return {
            "status": "no_candidates",
            "response": {"candidates": []},
            "resources": {"sampled_cpu_seconds_lower_bound": 0, "sampled_peak_rss_mb": 0},
        }

    def artifacts(directory, response, query, expected_source, **kwargs):
        collected.append(kwargs)
        return []

    monkeypatch.setattr(runner, "bounded_method", bounded)
    monkeypatch.setattr(runner.qualified_collector, "checked_workflow_artifacts", artifacts)
    result = runner.attempt(
        p["methods"][0], p["queries"][0], 97, tmp_path / "attempt", plan, SimpleNamespace(signal_number=None)
    )
    assert len(requests) == 2 and all(req["seed"] == 97 for req in requests)
    assert ("benchmark_protocol" in requests[0]) == profile_binding["schema"].endswith("v3")
    assert plan["benchmark_protocol"] == profile_binding  # request replacement cannot alter authority
    if tamper:
        assert not collected and result["status"] == "error"
        assert "frozen runner authority" in result["failures"][0]["reason"]
    else:
        assert result["status"] == "no_candidates" and len(collected) == 1
        assert collected[0]["profile_name"] == profile_binding["recommendation_profile"]
        cfg = collected[0]["expected_config"]
        assert cfg["provider"]["seed"] == 97 and cfg["conformers"]["seed"] == 42
        assert cfg["refinement"]["seeds"] == [42, 43, 44]
        assert requests[1]["benchmark_protocol"] == profile_binding


@pytest.mark.parametrize("attack", [None, "policy", "outcome", "missing_exports", "header", "nonempty_sdf"])
def test_benchmark_metadata_selective_empty_never_uses_legacy_return(
    tmp_path, monkeypatch, profile_binding, benchmark_metadata_guard, attack
):
    """Exercise real empty-artifact checks; mock only chemistry-derived similarity identity.

    Full saved pose/seed/route/raw evidence tests for both versions live in
    synthesis/test_qualified_collector.py and require scientific toy fixtures.
    """
    import copy
    from chembounce_synthesis import selective_scoring, workflow
    from chembounce_synthesis.qualification import effective_policy

    collector = runner.qualified_collector
    config = config_file(tmp_path, profile_binding["recommendation_profile"])
    req = request(tmp_path)
    req["assets"]["synthesis_config"] = str(config)
    req["benchmark_protocol"] = profile_binding
    cfg = workflow._config(adapter.common_config(req))
    similarity = SimpleNamespace(policy={"synthetic_metadata_only": True}, sha256="a" * 64)
    monkeypatch.setattr(selective_scoring, "configured_constraints", lambda *args: similarity)
    monkeypatch.setattr(collector, "canonical", lambda value: value)
    policy = effective_policy(cfg)
    policy["similarity_constraints"] = {"policy": similarity.policy, "sha256": similarity.sha256}
    manifest = {
        "sources_changed_since_import": False,
        "code_sha256": workflow._SOURCE_HASHES,
        "hit": req["query"]["smiles"],
        "config": cfg,
        "recommendation_policy": copy.deepcopy(policy),
        "recommendation_profile": policy["calibration"],
    }
    summary = {
        "status": "no_candidates",
        "returned_products": 0,
        "accepted_products": 0,
        "qualified_products": 0,
        "refinement_accepted": 0,
        "refinement_attempted": 0,
        "route_rejections": {},
        "refinement_budget": {"candidate_conformer_requests": 0},
    }
    if attack == "policy":
        policy["unbound_extra_policy"] = True
    if attack == "outcome":
        summary["qualified_products"] = 1
    out = tmp_path / "synthesis"
    out.mkdir()
    for name, value in {
        "manifest.json": manifest,
        "results.json": {"results": [], "policy": policy},
        "qualification.json": {"policy": policy, "selected_route_ids": [], "decisions": []},
        "routes.json": {"routes": []},
        "refinement.json": {"candidates": []},
        "summary.json": summary,
    }.items():
        write_json(out / name, value)
    if attack != "missing_exports":
        (out / "results.tsv").write_text(
            "smiles\n"
            if attack == "header"
            else "\t".join(workflow.COLUMNS + workflow.SELECTIVE_COLUMNS) + "\n"
        )
        for name in ("poses_hit.sdf", "poses_design.sdf", "references_hit.sdf", "references_design.sdf"):
            (out / name).write_text("forged\n" if attack == "nonempty_sdf" else "")
    response = {"candidates": [], "metadata": {"files_sha256": runner.file_manifest(tmp_path)}}

    def collect():
        return collector.checked_workflow_artifacts(
            tmp_path,
            response,
            req["query"],
            workflow._SOURCE_HASHES,
            expected_config=cfg,
            profile_name=profile_binding["recommendation_profile"],
        )

    if attack is None:
        assert collect() == []
    else:
        with pytest.raises(ValueError, match="Selective|bound exports|Empty selective"):
            collect()


@pytest.mark.parametrize("attack", ["enabled_proposals", "missing_helper", "changed_helper"])
def test_benchmark_metadata_v5_rejects_unmatched_configuration(
    tmp_path, monkeypatch, benchmark_metadata_guard, attack
):
    from chembounce_synthesis import profiles, landmark_proposals

    req = request(tmp_path)
    req["benchmark_protocol"] = selective_binding("v3")
    path = config_file(tmp_path, profiles.LANDMARK_SELECTIVE_PROFILE)
    req["assets"]["synthesis_config"] = str(path)
    if attack == "enabled_proposals":
        config = adapter.read(path)
        config["similarity_constraints"] = {
            "landmarks": [{"id": "oxygen", "family": "Acceptor", "hit_atom_ids": [2], "max_distance": 1.0}],
            "landmark_pose_proposals": True,
        }
        write_json(path, config)
        # Isolate profile rejection from optional native capability availability.
        monkeypatch.setattr(landmark_proposals, "require_capability", lambda: None)
    elif attack == "missing_helper":
        from pathlib import Path

        original_read = Path.read_bytes

        def missing(path):
            if path == Path(profiles.__file__).with_name("landmark_proposals.py"):
                raise FileNotFoundError("controlled missing landmark helper source")
            return original_read(path)

        # Real numerical_runtime always reads the fixed helper list; a missing
        # file raises instead of returning a fabricated incomplete source map.
        monkeypatch.setattr(Path, "read_bytes", missing)
        with pytest.raises(FileNotFoundError, match="missing landmark helper"):
            adapter.common_config(req)
        return
    else:
        original = profiles.numerical_runtime

        def changed(name):
            current = original(name)
            current["source_sha256"][profiles.LANDMARK_FILE] = "0" * 64
            return current

        monkeypatch.setattr(profiles, "numerical_runtime", changed)
    with pytest.raises(ValueError, match="frozen engineering profile"):
        adapter.common_config(req)
