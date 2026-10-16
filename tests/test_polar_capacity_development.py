"""Controlled source tests. No scientific assets, chemistry, model or tokenizer calls."""

import copy
import json
import os
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
from functools import lru_cache

import pytest

from benchmarks.research import polar_capacity_development as study
from benchmarks.research import polar_capacity_observer as observation
from benchmarks.research import polar_capacity_supervisor as supervisor


def config(archive="/archive.json"):
    return {
        "schema_version": 1,
        "recommendation_profile": study.PROFILE,
        "mode": "hop",
        "catalog": "/catalog.tsv",
        "rules": "/rules.json",
        "provider": {"kind": "saved", "path": archive},
        "minimum_scores": {"hit_tanimoto": 0.35, "hit_shape": 0.7, "hit_pharmacophore": 0.5},
        "similarity_constraints": {
            "version": "selective-similarity-constraints-v2",
            "polar_recall": 0.75,
            "landmarks": [],
        },
        "route_search_policy": "original-hit-v1",
        "ranking_policy": "original-hit-v1",
        "route_search_cache_size": 4096,
        "allocation_policy": "balanced",
        "max_design_attempts": 12,
        "max_designs": 3,
        "max_scored_products": 12,
        "max_results": 3,
        "conformers": {
            "scoring_objective": "original-hit-v1",
            "num_conformers": 8,
            "seed": 42,
            "energy_window": 5.0,
            "strain_penalty": 0.05,
            "max_iterations": 500,
            "pharmacophore_sigma": 1.0,
            "reject_unspecified_stereo": True,
        },
        "refinement": {
            "top_n": 6,
            "num_conformers": 8,
            "seeds": [42, 43, 44],
            "threshold_policy": "all_seeds",
            "ranking_statistic": "minimum",
            "max_total_conformers": 144,
            "max_score_range": None,
            "max_metric_range": None,
            "adaptive_max_conformers": None,
            "conformer_growth_factor": 2,
            "borderline_margin": 0.0,
        },
        "search": {
            "beam_width": 6,
            "building_block_top_k": 3,
            "max_combinations": 256,
            "max_steps": 12,
            "reaction_product_limit": 64,
        },
        "expansion": {"building_block_top_k": 6, "max_routes": 12, "min_tanimoto": 0.0, "rounds": 1},
        "synthesis_evidence": {
            "profile": "computational-route-support-v1",
            "policy": {"require_forward_agreement": True},
            "max_predictor_calls": 24,
            "max_input_tokens": 12288,
            "max_output_tokens": 18432,
            "predictor": {
                "kind": "reactiont5",
                "device": "cpu",
                "top_k": 3,
                "max_new_tokens": 256,
                "max_input_tokens": 512,
                "model_revision": "933114058cb2604dc1bf536dbebdfcefbe83d4fc",
                "model_path": "/model",
            },
        },
    }


@pytest.mark.parametrize(
    "field,value",
    [
        ("reference_pose", None),
        ("reference_sdf", "/ref.sdf"),
        ("pocket", {}),
        ("core_smiles", "C"),
        ("carry_forward", ["old"]),
        ("max_design_attempts", 13),
        ("max_scored_products", 13),
        ("allocation_policy", "encounter"),
        ("recommendation_profile", "qualified-development-v1"),
    ],
)
def test_fixed_config_conflicts(field, value):
    data = config()
    study.validate_config(data, "/archive.json")
    data[field] = value
    with pytest.raises(ValueError):
        study.validate_config(data, "/archive.json")


@pytest.mark.parametrize("field,value", [("polar_recall", 0.74), ("landmarks", [{}])])
def test_similarity_not_tuned(field, value):
    data = config()
    data["similarity_constraints"][field] = value
    with pytest.raises(ValueError):
        study.validate_config(data, "/archive.json")


def test_null_actual_reference_and_exact_seed_config():
    data = config()
    data["reference_sdf"] = None
    study.validate_config(data, "/archive.json")
    data["refinement"]["seeds"] = [42, 43]
    with pytest.raises(ValueError):
        study.validate_config(data, "/archive.json")


def archive(tmp_path):
    pool = ["literal-" + str(i) for i in range(12)]
    rows = [
        {"design": pool[i], "response": "malformed-but-preserved", "sampling": {"seed": seed}}
        for i, seed in zip([1, 1, 2, 2, 3, 3], [42, 43, 42, 43, 42, 43])
    ]
    query = {"id": "known", "target_id": "known", "smiles": "literal"}
    for name, value in [
        ("pool", {"designs": pool}),
        ("proposals", {"proposals": rows}),
        ("request", {"query": query}),
    ]:
        supervisor.write(tmp_path / (name + ".json"), value)
    return {
        "design_pool": str(tmp_path / "pool.json"),
        "proposals": str(tmp_path / "proposals.json"),
        "request": str(tmp_path / "request.json"),
        "query": query,
    }


@pytest.mark.parametrize("mutation", ["reorder", "seed", "drop", "query", "pool"])
def test_archive_retains_all_coverage_and_rejects_mutation(tmp_path, mutation):
    q = archive(tmp_path)
    assert len(study.validate_archive(q)) == 12
    path = Path(q["proposals"])
    data = study.read(path)
    if mutation == "reorder":
        data["proposals"].reverse()
    elif mutation == "seed":
        data["proposals"][0]["sampling"]["seed"] = 44
    elif mutation == "drop":
        data["proposals"].pop()
    elif mutation == "query":
        q["query"] = dict(q["query"], id="changed")
    else:
        supervisor.write(q["design_pool"], {"designs": ["only-covered"] * 3})
    supervisor.write(path, data)
    with pytest.raises(ValueError):
        study.validate_archive(q)


def fake_plan(tmp_path):
    q = archive(tmp_path)
    q["config"] = config(q["proposals"])
    return {
        "queries": {"attempt-0000": q, "attempt-0009": copy.deepcopy(q)},
        "sources": {
            arm: {"commit": commit, "synthesis_source": {"file.py": "hash"}}
            for arm, commit in study.ARMS.items()
        },
        "pins": {q[k]: study.digest(q[k]) for k in ["design_pool", "proposals", "request"]},
        "observer": {"baseline": [{"kind": "fake"}], "capacity": [{"kind": "fake"}]},
    }


@pytest.mark.parametrize(
    "key",
    ["index", "arm", "historical_slot", "config", "source_commit", "archive", "observer_sha256", "query"],
)
def test_development_request_authenticates_before_dispatch(tmp_path, key):
    plan = fake_plan(tmp_path)
    request = study.request_for(plan, "plan-hash", 0)
    assert study.authenticate_request(request, plan, "plan-hash") == 0
    request[key] = "changed"
    with pytest.raises(ValueError):
        study.authenticate_request(request, plan, "plan-hash")


@pytest.mark.parametrize(
    "allocation,elapsed,status",
    [
        (600, 600, "completed"),
        (600, 601, "allocation_exceeded"),
        (300, 300, "completed"),
        (300, 301, "allocation_exceeded"),
    ],
)
def test_exact_success_allocation(allocation, elapsed, status):
    assert supervisor.terminal_status(None, 0, elapsed, allocation, True, True) == status


@pytest.mark.parametrize(
    "reason,code,clean,durable,expected",
    [
        (None, 1, True, True, "child_error"),
        (None, 0, False, True, "unsafe_cleanup"),
        (None, 0, True, False, "missing_terminal"),
        ("memory_limit", 0, True, True, "memory_limit"),
        ("resource_verification_failed", 0, True, True, "resource_verification_failed"),
    ],
)
def test_terminal_failures_never_success(reason, code, clean, durable, expected):
    assert supervisor.terminal_status(reason, code, 1, 600, clean, durable) == expected


def launch(tmp_path, body, **kwargs):
    out = tmp_path / "child"
    code = "import os,sys,time,signal,json\nfrom pathlib import Path\n" + body
    return supervisor.supervise(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=os.environ.copy(),
        output=out,
        allocation=kwargs.pop("allocation", 2),
        watchdog=kwargs.pop("watchdog", 3),
        grace=0.15,
        interval=0.02,
        **kwargs,
    )


def test_real_fake_child_completion_and_receipt_hash(tmp_path):
    r = launch(tmp_path, "Path('child/terminal.json').write_text(json.dumps({'status':'fake'}))")
    assert r["status"] == "completed" and r["cleanup"]["complete"]
    assert r["terminal_sha256"] == study.digest(tmp_path / "child/terminal.json")
    assert r["resources"]["lifetime_ru_maxrss_kib"] > 0
    assert r["resources"]["sampled_tree_peak_rss_bytes"] is not None


def test_real_fake_timeout_ignoring_term(tmp_path):
    r = launch(
        tmp_path,
        "signal.signal(signal.SIGTERM,signal.SIG_IGN)\ntime.sleep(10)",
        allocation=0.05,
        watchdog=0.1,
    )
    assert r["status"] == "watchdog_timeout"
    assert r["cleanup"]["kill_monotonic"] is not None
    assert r["wall_seconds"] < 1


def test_reuse_output_refused(tmp_path):
    launch(tmp_path, "Path('child/terminal.json').write_text('{}')")
    with pytest.raises(FileExistsError):
        launch(tmp_path, "pass")


def test_workload_blocks_launch(tmp_path):
    r = launch(tmp_path, "raise AssertionError('must not launch')", check_workload=lambda b: [dict(pid=123)])
    assert r["status"] == "workload_overlap" and r["identity"] is None


def test_rss_gap_fails_not_zero(tmp_path):
    r = launch(
        tmp_path, "time.sleep(10)", sample=lambda p: {"rss_bytes": None, "gaps": [{"error": "AccessDenied"}]}
    )
    assert r["status"] == "resource_verification_failed"
    assert r["resources"]["sampled_tree_peak_rss_bytes"] is None
    assert r["resources"]["gaps"]


def test_rss_exceedance_stops_only_owned(tmp_path):
    r = launch(tmp_path, "time.sleep(10)", memory_bytes=1)
    assert r["status"] == "memory_limit" and r["cleanup"]["complete"]


def test_pidfd_rechecks_identity_before_signal(monkeypatch):
    entry = {"pid": 123, "start_ticks": 99}
    closed = []
    monkeypatch.setattr(supervisor, "pidfd_open", lambda pid: 7)
    monkeypatch.setattr(supervisor, "identity", lambda pid: {"pid": 123, "start_ticks": 100})
    monkeypatch.setattr(supervisor.os, "close", closed.append)
    monkeypatch.setattr(supervisor, "pidfd_send_signal", lambda *a: pytest.fail("foreign PID signaled"))
    assert supervisor.signal_owned(entry, 15) is False and closed == [7]


def test_fake_surviving_descendant_is_owned_and_cleaned(tmp_path):
    r = launch(
        tmp_path,
        "pid=os.fork()\nif pid==0:\n os.setsid()\n time.sleep(10)\nelse:\n Path('child/terminal.json').write_text('{}')\n time.sleep(.05)",
    )
    assert r["status"] == "owned_descendants_after_exit"
    assert r["cleanup"]["complete"]


def test_preexisting_child_not_signaled(tmp_path):
    import subprocess

    p = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(10)"])
    try:
        r = launch(tmp_path, "pass")
        assert r["status"] == "preexisting_children" and p.poll() is None
    finally:
        p.terminate()
        p.wait()


def test_observer_alias_lru_results_and_exception_transparency(tmp_path, monkeypatch):
    module = ModuleType("observer_fixture")
    module.__file__ = __file__
    calls = []

    @lru_cache(maxsize=2)
    def original(n):
        calls.append(n)
        if n == 9:
            raise ValueError("preserved")
        return {"n": n}

    module.fn = original
    alias = ModuleType("observer_fixture.alias")
    alias.__file__ = __file__
    alias.fn = original
    monkeypatch.setitem(sys.modules, module.__name__, module)
    monkeypatch.setitem(sys.modules, alias.__name__, alias)
    observer = observation.Observer(tmp_path / "events.jsonl", module_prefix="observer_fixture")
    entry = observer.install(module, "fn", kind="ensemble")
    observer.expected = [copy.deepcopy(entry)]
    observer.seal()
    assert module.fn is alias.fn
    assert module.fn(1) is alias.fn(1)
    module.fn(2)
    module.fn(3)
    module.fn(1)
    assert original.cache_info().hits == 1 and original.cache_info().misses == 4
    with pytest.raises(ValueError, match="preserved"):
        module.fn(9)
    observer.close()
    assert module.fn is alias.fn is original
    events = [json.loads(s) for s in (tmp_path / "events.jsonl").read_text().splitlines()]
    assert study.event_summary(events, expected_inventory=[entry])["calls"]["ensemble"] == {
        "started": 6,
        "returned": 5,
        "failed": 1,
        "nonreturned": 0,
    }
    assert calls == [1, 2, 3, 1, 9]


def test_observer_inventory_mismatch_is_before_replacement(tmp_path):
    module = ModuleType("fake")
    module.__file__ = __file__
    module.fn = lambda x: x
    original = module.fn
    observer = observation.Observer(tmp_path / "events", expected_inventory=[{}])
    with pytest.raises(ValueError):
        observer.install(module, "fn", kind="fake")
    assert module.fn is original
    observer.close()


def test_observer_unknown_nonreturn_retained():
    started = {"event": "started", "id": 1, "kind": "generation", "call_key": "key", "physical": True}
    assert study.event_summary([started])["calls"]["generation"]["nonreturned"] == 1


def product(name):
    return {
        "smiles": name,
        "family": name,
        "scaffold": name,
        "generic_scaffold": "generic",
        "route_id": "route-" + name,
        "route_support_binding": "binding-" + name,
    }


def valid_receipt(allocation, status="completed"):
    return {
        "status": status,
        "returncode": 0,
        "allocation_seconds": allocation,
        "wall_seconds": 1.0,
        "parent_observed_phase_seconds": 1.1,
        "cleanup": {"complete": True, "remaining": []},
        "identity": {"pid": 123, "start_ticks": 1},
        "terminal_sha256": "a" * 64,
        "resources": {"samples": 1, "gaps": [], "sampled_tree_peak_rss_bytes": 1024},
    }


def outcome(index, products=(), status="completed"):
    return {
        "index": index,
        "status": status,
        "worker": valid_receipt(600, status),
        "assessor": valid_receipt(300, status),
        "attestation": {
            "status": "validated",
            "index": index,
            "coverage_failures": [],
            "qualified": list(products),
        },
    }


def test_four_row_oracle_failures_and_diversity():
    rows = [
        outcome(0, [product("a")]),
        outcome(1, [product("a"), product("b"), product("c")]),
        outcome(2, [product("must-not-count")], status="allocation_exceeded"),
        {"index": 3, "status": "not_started"},
    ]
    r = study.summarize(rows)
    assert (
        r["arms"]["baseline"]["primary"]
        == r["arms"]["capacity"]["primary"]
        == {"numerator": 1, "denominator": 2}
    )
    assert r["arms"]["baseline"]["family_yield"] == {"numerator": 1, "denominator": 6}
    assert r["arms"]["capacity"]["family_yield"] == {"numerator": 3, "denominator": 6}
    assert r["rows"][0]["diversity"] is None
    assert r["rows"][1]["diversity"] == {"normalized_scaffolds": 3, "generic_scaffolds": 1}
    assert r["pairs"][0]["new"] == ["b", "c"]
    assert r["rows"][2]["primary"] == 0


@pytest.mark.parametrize("bad", [[], [{"index": 0}] * 4, [{"index": i} for i in range(5)]])
def test_ledger_duplicate_missing_extra_rejected(bad):
    with pytest.raises(ValueError):
        study.summarize(bad)


def test_call_multiset_counts_not_rejections():
    def e(k, physical=True):
        return {"event": "started", "kind": "generation", "call_key": k, "physical": physical}

    r = study.call_difference([e("a"), e("a"), e("b")], [e("a"), e("c"), e("rejection", False)])
    assert r["removed_started"] == [
        {"kind": "generation", "call_key": "a", "count": 1},
        {"kind": "generation", "call_key": "b", "count": 1},
    ]
    assert r["new_started"][0]["call_key"] == "c"


def test_default_collector_authenticates_before_common_delegate(tmp_path, monkeypatch):
    from benchmarks.research import qualified_collector as collector
    from chembounce_synthesis import workflow

    calls = []
    monkeypatch.setattr(
        collector, "checked_workflow_artifacts", lambda *a, **k: calls.append(k) or ["checked"]
    )
    monkeypatch.setattr(collector, "common_config", lambda *a: {"fixed": "config"})
    monkeypatch.setattr(workflow, "_config", lambda value: value)
    supervisor.write(tmp_path / "request.json", {"query": {"id": "q"}, "settings": {"arm": "full"}})
    assert collector.checked_products(tmp_path, {}, {"id": "q"}, {}) == ["checked"]
    assert calls[0]["expected_config"] == {"fixed": "config"}
    with pytest.raises(ValueError):
        collector.checked_products(tmp_path, {}, {"id": "other"}, {})
    assert len(calls) == 1


def test_shared_checker_rejects_partial_before_assets(tmp_path, monkeypatch):
    from benchmarks.research import qualified_collector as collector

    monkeypatch.setattr(collector, "Catalog", lambda *a: pytest.fail("asset constructor"))
    with pytest.raises(ValueError, match="Incomplete"):
        collector.checked_workflow_artifacts(
            tmp_path,
            {"metadata": {"files_sha256": {}}},
            {},
            {},
            expected_config={},
            profile_name=study.PROFILE,
        )


def test_development_checker_rejects_request_before_import_or_assets(tmp_path, monkeypatch):
    plan = fake_plan(tmp_path)
    request = study.request_for(plan, "hash", 0)
    request["query"] = {"changed": True}
    monkeypatch.setattr(study, "forbid_assessor_compute", lambda: pytest.fail("too early"))
    with pytest.raises(ValueError):
        study.assess(plan, "hash", request, tmp_path, {}, 0)


def complete_plan(tmp_path):
    plan = fake_plan(tmp_path)
    files = {}
    kernels = (
        "common.py",
        "identity.py",
        "conformers.py",
        "selective_scoring.py",
        "similarity_constraints.py",
    )
    for arm, source in [*plan["sources"].items(), ("checker", {})]:
        root = tmp_path / arm
        source.update(
            path=str(root), runtime={"fake": arm}, files={}, synthesis_source=dict.fromkeys(kernels, "kernel")
        )
        source.setdefault("commit", "future-reviewed-checker")
        for name, value in source["synthesis_source"].items():
            source["files"][str(root / "chembounce_synthesis" / name)] = value
        source["files"][str(root / "chembounce_synthesis/data" / f"{study.PROFILE}.json")] = study.PROFILE_SHA
        source["files"][str(root / "chembounce_synthesis/data/qualified-development-v1.json")] = (
            "9d3c60a895912b28c13a29aac0a04f385ddb8a0effc9fcd81d64a50a3d48bdd2"
        )
        files.update(source["files"])
        if arm == "checker":
            plan["checker"] = source
    for name in (
        "polar_capacity_development.py",
        "polar_capacity_observer.py",
        "polar_capacity_supervisor.py",
        "qualified_collector.py",
        "qualified_adapter.py",
    ):
        filename = str(Path(study.__file__).with_name(name))
        plan["checker"]["files"][filename] = "source"
        files[filename] = "source"
    plan.update(
        schema="polar-capacity-development-plan-v2",
        design={"path": "/design", "sha256": study.DESIGN_SHA},
        gate={"path": "/gate", "sha256": study.GATE_SHA},
        python=sys.executable,
        output=str(tmp_path / "run"),
        workloads={"identities": [{"pid": 0, "start_ticks": 0}], "command_markers": ["known-science"]},
        administrative_runtime={"fake": "parent"},
        feature_definition={"path": "/features", "sha256": "feature"},
    )
    plan["pins"].update(files)
    plan["pins"].update(
        {
            "/design": study.DESIGN_SHA,
            "/gate": study.GATE_SHA,
            "/features": "feature",
            str(Path(sys.executable).resolve()): "python",
            "/catalog.tsv": "catalog",
            "/catalog.tsv.manifest.json": "sidecar",
            "/rules.json": "rules",
        }
    )
    for arm in study.ARMS:
        plan["observer"][arm] = [{"source_path": "/hook", "source_sha256": "hook"}]
    plan["pins"]["/hook"] = "hook"
    return plan


def test_plan_structural_authority_and_kernels(tmp_path):
    plan = complete_plan(tmp_path)
    study.validate_plan(plan, check_files=False)
    plan["sources"]["capacity"]["synthesis_source"]["common.py"] = "mutated"
    with pytest.raises(ValueError):
        study.validate_plan(plan, check_files=False)


@pytest.mark.parametrize(
    "mutation", ["extra", "arm", "gate", "pin", "observer", "feature", "sidecar", "source", "config"]
)
def test_plan_rejects_unresolved_or_conflicting_authority(tmp_path, mutation):
    plan = complete_plan(tmp_path)
    if mutation == "extra":
        plan["unknown"] = True
    elif mutation == "arm":
        plan["sources"]["other"] = plan["sources"]["baseline"]
    elif mutation == "gate":
        plan["gate"]["sha256"] = "bad"
    elif mutation == "pin":
        plan["pins"].pop("/hook")
    elif mutation == "observer":
        plan["observer"]["capacity"] = []
    elif mutation == "feature":
        plan["pins"].pop("/features")
    elif mutation == "sidecar":
        plan["pins"].pop("/catalog.tsv.manifest.json")
    elif mutation == "source":
        plan["sources"]["capacity"]["commit"] = "different"
    else:
        plan["queries"]["attempt-0009"]["config"]["max_results"] = 4
    with pytest.raises(ValueError):
        study.validate_plan(plan, check_files=False)


def test_four_slot_runner_fake_lifecycle_and_denominators(tmp_path, monkeypatch):
    plan = complete_plan(tmp_path)
    plan_path = tmp_path / "plan.json"
    supervisor.write(plan_path, plan)
    sha = study.digest(plan_path)
    auth = tmp_path / "authorization.json"
    supervisor.write(
        auth,
        {
            "schema": "polar-capacity-development-launch-v1",
            "plan_sha256": sha,
            "source_review_passed": True,
            "resource_clearance": True,
        },
    )
    monkeypatch.setattr(study, "runtime", lambda: plan["administrative_runtime"])
    calls = []

    def fake_supervise(argv, **kw):
        command = argv[2]
        out = kw["output"]
        out.mkdir(parents=True)
        ledger = study.read(Path(plan["output"]) / "ledger.json")
        assert len(ledger["rows"]) == 4
        assert [r["index"] for r in ledger["rows"]] == [0, 1, 2, 3]
        calls.append(command)
        terminal = {"status": "verified"}
        status = "completed"
        if command == "worker":
            request = study.read(argv[argv.index("--request") + 1])
            index = request["index"]
            assert ledger["rows"][index]["status"] == "worker_running"
            if index == 0:
                status = "allocation_exceeded"
            terminal = {"status": "returned"}
        if command == "assess":
            request = study.read(argv[argv.index("--request") + 1])
            index = request["index"]
            terminal = {
                "status": "validated",
                "index": index,
                "coverage_failures": [],
                "qualified": [product(str(index))] if index == 1 else [],
            }
        supervisor.write(out / "terminal.json", terminal)
        result = valid_receipt(kw["allocation"], status)
        result["wall_seconds"] = 0.0
        return result

    monkeypatch.setattr(study, "supervise", fake_supervise)
    args = SimpleNamespace(
        authorization=str(auth),
        authorization_sha256=study.digest(auth),
        plan=str(plan_path),
        plan_sha256=sha,
        output=plan["output"],
    )
    assert study.run(args) == 0
    ledger = study.read(Path(plan["output"]) / "ledger.json")
    assert calls == ["admin", "worker", "worker", "assess", "worker", "assess", "worker", "assess", "admin"]
    assert ledger["summary"]["arms"]["capacity"]["primary"] == {"numerator": 1, "denominator": 2}
    assert ledger["summary"]["arms"]["baseline"]["primary"] == {"numerator": 0, "denominator": 2}
    assert ledger["rows"][0]["status"] == "allocation_exceeded"


def test_shared_selective_empty_exports_still_checked(tmp_path, monkeypatch):
    from benchmarks.research import qualified_collector as collector
    from chembounce_synthesis import workflow, qualification, selective_scoring

    cfg = {"catalog": "never-opened", "rules": "never-opened"}
    policy = {
        "mode": "qualified",
        "evidence_profile": collector.PROFILE,
        "conformers": collector.FIXED_CONFORMERS,
        "calibration": {
            "artifact_sha256": study.PROFILE_SHA,
            "current_numerical_runtime": {"source_sha256": {}},
        },
        "minimum_scores": {"a": 0.5},
        "similarity_constraints": {"policy": {}, "sha256": "policy"},
    }
    monkeypatch.setattr(collector, "canonical", lambda x: x)
    monkeypatch.setattr(collector, "Catalog", lambda *a: pytest.fail("asset construction"))
    monkeypatch.setattr(
        selective_scoring, "configured_constraints", lambda *a: SimpleNamespace(policy={}, sha256="policy")
    )
    monkeypatch.setattr(qualification, "effective_policy", lambda *a: copy.deepcopy(policy))
    monkeypatch.setattr(collector, "load_profile", lambda *a: {"minimum_scores": {"a": 0.5}})
    out = tmp_path / "synthesis"
    out.mkdir()
    documents = {
        "manifest.json": {
            "sources_changed_since_import": False,
            "code_sha256": {},
            "hit": "literal",
            "config": cfg,
            "recommendation_policy": policy,
            "recommendation_profile": policy["calibration"],
        },
        "results.json": {"results": [], "policy": policy},
        "qualification.json": {"policy": policy, "selected_route_ids": [], "decisions": []},
        "routes.json": {"routes": []},
        "refinement.json": {"candidates": []},
        "summary.json": {
            "status": "no_candidates",
            "route_rejections": {},
            "returned_products": 0,
            "accepted_products": 0,
            "qualified_products": 0,
            "refinement_accepted": 0,
            "refinement_attempted": 0,
            "refinement_budget": {"candidate_conformer_requests": 0},
        },
    }
    for name, value in documents.items():
        supervisor.write(out / name, value)
    (out / "results.tsv").write_text("\t".join(workflow.COLUMNS + workflow.SELECTIVE_COLUMNS) + "\n")
    for name in ["poses_hit.sdf", "poses_design.sdf", "references_hit.sdf", "references_design.sdf"]:
        (out / name).write_text("")

    def check():
        response = {"metadata": {"files_sha256": study.inventory(tmp_path)}, "candidates": []}
        return collector.checked_workflow_artifacts(
            tmp_path, response, {"smiles": "literal"}, {}, expected_config=cfg, profile_name=study.PROFILE
        )

    assert check() == []
    for name in ["poses_hit.sdf", "poses_design.sdf", "references_hit.sdf", "references_design.sdf"]:
        (out / name).write_text("unexpected")
        with pytest.raises(ValueError, match="nonempty"):
            check()
        (out / name).write_text("")
    documents["summary.json"]["status"] = "error"
    supervisor.write(out / "summary.json", documents["summary.json"])
    with pytest.raises(ValueError, match="outcome"):
        check()


def test_stale_session_number_does_not_adopt_foreign_process(monkeypatch):
    leader = {"pid": 40, "start_ticks": 1}
    foreign = {"pid": 40, "start_ticks": 2, "session": 40, "ppid": 1, "state": "S"}
    monkeypatch.setattr(supervisor.os, "listdir", lambda p: ["40"])
    monkeypatch.setattr(supervisor, "identity", lambda pid: foreign)
    assert supervisor.owned_snapshot(leader, {40: leader}) == {}


def test_owned_pidfd_stays_on_kernel_object_after_recheck(monkeypatch):
    entry = {"pid": 12, "start_ticks": 34}
    calls = []
    monkeypatch.setattr(supervisor, "pidfd_open", lambda p: 99)
    monkeypatch.setattr(supervisor, "identity", lambda p: entry)
    monkeypatch.setattr(supervisor, "pidfd_send_signal", lambda fd, n: calls.append((fd, n)))
    monkeypatch.setattr(supervisor.os, "close", lambda fd: calls.append(("close", fd)))
    assert supervisor.signal_owned(entry, 15)
    assert calls == [(99, 15), ("close", 99)]


def test_observer_python_implementation_detects_wraps_disguise():
    import functools

    def original(x):
        return x

    @functools.wraps(original)
    def replacement(x):
        return x + 1

    assert original.__name__ == replacement.__name__
    assert observation.implementation_identity(original) != observation.implementation_identity(replacement)


def test_embedding_keys_have_exact_enclosing_smiles_and_seed(tmp_path, monkeypatch):
    module = ModuleType("context_fixture")
    module.__file__ = __file__
    monkeypatch.setitem(sys.modules, module.__name__, module)
    module.embed = lambda same_size, params=None: [0]

    def ensemble(smiles, settings):
        return module.embed(6, params=object())

    module.ensemble = ensemble
    observer = observation.Observer(tmp_path / "events", module_prefix="context_fixture")
    a = observer.install(module, "ensemble", kind="ensemble")
    b = observer.install(module, "embed", kind="embedding", physical=True)
    observer.expected = [a, b]
    observer.seal()
    for smiles, seed in [("same-size-A", 42), ("same-size-A", 43), ("same-size-B", 42)]:
        module.ensemble(smiles, {"seed": seed, "num_conformers": 8})
    observer.close()
    events = [json.loads(line) for line in (tmp_path / "events").read_text().splitlines()]
    starts = [e for e in events if e["event"] == "started" and e["kind"] == "embedding"]
    assert len({e["call_key"] for e in starts}) == 3
    assert all(e["scientific_context"][0]["kind"] == "ensemble" for e in starts)


@pytest.mark.parametrize(
    "mutation", ["late", "cleanup", "coverage", "label_only", "assessor_late", "no_parent_time"]
)
def test_reducer_does_not_trust_success_labels(mutation):
    rows = [outcome(i, [product(str(i))]) for i in range(4)]
    row = rows[0]
    if mutation == "late":
        row["worker"]["wall_seconds"] = 601
    elif mutation == "cleanup":
        row["worker"]["cleanup"]["complete"] = False
    elif mutation == "coverage":
        row["attestation"]["coverage_failures"] = ["missing design"]
    elif mutation == "label_only":
        row["worker"] = {"status": "completed"}
    elif mutation == "assessor_late":
        row["assessor"]["parent_observed_phase_seconds"] = 301
    else:
        row["worker"].pop("parent_observed_phase_seconds")
    assert study.summarize(rows)["rows"][0]["primary"] == 0


def test_scientific_child_dispatch_requires_launch_authority(tmp_path):
    args = SimpleNamespace(
        authorization=None,
        authorization_sha256=None,
        output=str(tmp_path),
        deadline=999,
        plan_sha256="x",
        command="worker",
    )
    with pytest.raises(ValueError, match="authorization"):
        study.child_authority(args, {"workloads": {}})


def test_direct_child_cannot_spoof_parent_identity(tmp_path, monkeypatch):
    monkeypatch.setattr(study, "launch_authority", lambda args: {})
    supervisor.write(
        tmp_path / "receipt.json",
        {
            "status": "running",
            "supervisor_identity": {"pid": 1, "start_ticks": 1},
            "identity": {"pid": 2, "start_ticks": 2},
        },
    )
    monkeypatch.setattr(study, "identity", lambda pid: {"pid": 3, "start_ticks": 3})
    args = SimpleNamespace(output=str(tmp_path), deadline=100)
    with pytest.raises(ValueError, match="supervisor"):
        study.child_authority(args, {"workloads": {}})


def test_supervisor_late_durable_receipt_is_not_success(tmp_path, monkeypatch):
    import time

    offset = [0.0]
    original_write = supervisor.write

    def delayed_receipt(path, value):
        original_write(path, value)
        if value.get("status") == "completed":
            offset[0] += 3

    monkeypatch.setattr(supervisor, "write", delayed_receipt)
    result = launch(
        tmp_path, "Path('child/terminal.json').write_text('{}')", clock=lambda: time.monotonic() + offset[0]
    )
    assert result["status"] == "allocation_exceeded"
    assert result["wall_seconds"] > 2


def test_actual_missing_archive_exception_is_classified_without_science(tmp_path, monkeypatch):
    from chembounce_synthesis import workflow
    from rdkit import RDConfig

    plan = fake_plan(tmp_path)
    plan["checker"] = {}
    feature = tmp_path / "BaseFeatures.fdef"
    feature.write_text("synthetic-feature-binding")
    plan["feature_definition"] = {"path": str(feature), "sha256": study.digest(feature)}
    monkeypatch.setattr(RDConfig, "RDDataDir", str(tmp_path))
    monkeypatch.setattr(study, "context_check", lambda *a: None)
    monkeypatch.setattr(workflow, "_config", lambda x: x)
    monkeypatch.setattr(
        workflow, "effective_policy", lambda x: {"calibration": {"status": "engineering_policy_matched"}}
    )
    monkeypatch.setattr(workflow, "_resolve_paths", lambda x, p: x)
    monkeypatch.setattr(observation, "install_scientific", lambda o: None)

    def absent(*a, **k):
        raise ValueError("No saved proposals for design uncovered")

    monkeypatch.setattr(workflow, "run", absent)
    request = study.request_for(plan, "hash", 0)
    with pytest.raises(study.ArchiveCoverageError):
        study.worker(plan, "hash", request, tmp_path, 100)


def test_reused_terminal_is_not_overwritten(tmp_path):
    terminal = tmp_path / "terminal.json"
    terminal.write_text("preserve")
    with pytest.raises(FileExistsError):
        study.main(["worker", "--plan", "unused", "--plan-sha256", "unused", "--output", str(tmp_path)])
    assert terminal.read_text() == "preserve"


def test_complete_tracked_inventory_rejects_dirty_or_omitted_root(tmp_path, monkeypatch):
    source = {"path": str(tmp_path), "commit": "commit", "tree": "tree", "files": {}}

    def answer(argv, **kwargs):
        args = argv[1:]
        if args == ["rev-parse", "HEAD"]:
            return b"commit\n"
        if args == ["rev-parse", "HEAD^{tree}"]:
            return b"tree\n"
        if args[0] == "status":
            return b""
        if args == ["ls-files", "-z"]:
            return b"chembounce.py\0utils.py\0"
        raise AssertionError(args)

    monkeypatch.setattr(study.subprocess, "check_output", answer)
    with pytest.raises(study.IntegrityError, match="Complete tracked"):
        study.tracked_inventory(source)
    monkeypatch.setattr(
        study.subprocess,
        "check_output",
        lambda argv, **kw: b" M utils.py\n" if argv[1] == "status" else answer(argv, **kw),
    )
    with pytest.raises(study.IntegrityError, match="not clean"):
        study.tracked_inventory(source)


def test_direct_root_import_path_is_checked(tmp_path, monkeypatch):
    module = ModuleType("utils")
    module.__file__ = str(tmp_path / "outside.py")
    monkeypatch.setitem(sys.modules, "utils", module)
    current = {
        "python": "3.10.12",
        "packages": {
            "rdkit": "2025.9.6",
            "numpy": "1.26.4",
            "torch": "2.8.0",
            "transformers": "4.51.3",
            "tokenizers": "0.21.4",
            "safetensors": "0.6.2",
        },
    }
    monkeypatch.setattr(study, "runtime", lambda: current)
    with pytest.raises(study.IntegrityError, match="escaped"):
        study.context_check({"path": str(tmp_path / "source"), "runtime": current, "files": {}})


def test_outer_timer_is_literal_and_restored():
    import time

    with pytest.raises(TimeoutError, match="outer deadline"):
        with study.outer_watchdog(time.monotonic() + 0.02):
            time.sleep(0.2)
    assert study.signal.getitimer(study.signal.ITIMER_REAL) == (0.0, 0.0)


def test_assessor_guards_direct_entrypoints_with_only_fake_libraries(monkeypatch):
    # Import no torch/transformers/tokenizer library: this tests the actual guard
    # installer against a complete controlled module tree.
    names = [
        "torch",
        "transformers",
        "transformers.generation",
        "transformers.generation.utils",
        "rdkit",
        "rdkit.Chem",
        "rdkit.Chem.AllChem",
        "rdkit.Chem.rdMolAlign",
        "rdkit.Chem.rdDistGeom",
        "rdkit.Chem.rdShapeAlign",
        "safetensors",
        "safetensors.torch",
        "chembounce_synthesis",
        "chembounce_synthesis.conformers",
        "chembounce_synthesis.selective_scoring",
        "chembounce_synthesis.workflow",
        "chembounce_synthesis.forward_support",
    ]
    modules = {name: ModuleType(name) for name in names}
    for name, module in modules.items():
        module.__path__ = []
        monkeypatch.setitem(sys.modules, name, module)
        if "." in name:
            parent, attr = name.rsplit(".", 1)
            setattr(modules[parent], attr, module)
    modules["torch"].nn = SimpleNamespace(Module=type("FakeModule", (), {}))
    modules["transformers.generation.utils"].GenerationMixin = type("FakeGeneration", (), {})
    modules["chembounce_synthesis.forward_support"].PinnedReactionT5Forward = type("FakeForward", (), {})
    study.forbid_assessor_compute()
    targets = [
        modules["torch"].nn.Module.__init__,
        modules["torch"].load,
        modules["safetensors"].safe_open,
        modules["safetensors.torch"].load_file,
        modules["transformers.generation.utils"].GenerationMixin.generate,
        modules["rdkit.Chem.rdShapeAlign"].AlignMol,
        modules["rdkit.Chem.rdMolAlign"].GetO3A,
        modules["rdkit.Chem.rdDistGeom"].EmbedMultipleConfs,
        modules["rdkit.Chem.AllChem"].EmbedMolecule,
        modules["chembounce_synthesis.conformers"].ensemble,
        modules["chembounce_synthesis.workflow"].run,
        modules["chembounce_synthesis.forward_support"].PinnedReactionT5Forward.predict_raw,
    ]
    for target in targets:
        with pytest.raises(RuntimeError, match="forbids"):
            target()


@pytest.mark.parametrize(
    "scenario,expected_calls",
    [
        ("unsafe_setup", ["admin"]),
        ("failed_setup", ["admin", "admin"]),
        ("unsafe_worker", ["admin", "worker"]),
        ("drift", ["admin", "worker", "admin"]),
        ("interrupted", ["admin", "worker", "admin"]),
    ],
)
def test_systemic_stop_preserves_unstarted_and_suppresses_unsafe_finalization(
    tmp_path, monkeypatch, scenario, expected_calls
):
    plan = complete_plan(tmp_path)
    path = tmp_path / "plan.json"
    supervisor.write(path, plan)
    sha = study.digest(path)
    auth = tmp_path / "auth.json"
    supervisor.write(
        auth,
        {
            "schema": "polar-capacity-development-launch-v1",
            "plan_sha256": sha,
            "source_review_passed": True,
            "resource_clearance": True,
        },
    )
    monkeypatch.setattr(study, "runtime", lambda: plan["administrative_runtime"])
    calls = []

    def fake(argv, **kw):
        command = argv[2]
        calls.append(command)
        out = kw["output"]
        out.mkdir()
        result = valid_receipt(kw["allocation"])
        result["wall_seconds"] = 0
        terminal = {"status": "verified"}
        if len(calls) == 1 and scenario.endswith("setup"):
            result["status"] = "unsafe_cleanup" if scenario == "unsafe_setup" else "child_error"
            result["cleanup"]["complete"] = scenario != "unsafe_setup"
        elif command == "worker":
            result["status"] = (
                "unsafe_cleanup"
                if scenario == "unsafe_worker"
                else "interrupted"
                if scenario == "interrupted"
                else "child_error"
            )
            result["cleanup"]["complete"] = scenario != "unsafe_worker"
            terminal = {
                "status": "failed",
                "failure_category": "integrity" if scenario == "drift" else "execution",
            }
        supervisor.write(out / "terminal.json", terminal)
        return result

    monkeypatch.setattr(study, "supervise", fake)
    args = SimpleNamespace(
        authorization=str(auth),
        authorization_sha256=study.digest(auth),
        plan=str(path),
        plan_sha256=sha,
        output=plan["output"],
    )
    assert study.run(args) == 1
    ledger = study.read(Path(plan["output"]) / "ledger.json")
    assert calls == expected_calls
    assert all(r["status"] == "not_started" for r in ledger["rows"][1:])
    assert all(v["primary"]["numerator"] == 0 for v in ledger["summary"]["arms"].values())


def test_subreaper_restored_after_success_and_preexisting_rejection(tmp_path):
    import subprocess

    before = supervisor.subreaper()
    parent = tmp_path / "success"
    parent.mkdir()
    result = launch(parent, "Path('child/terminal.json').write_text('{}')")
    assert result["status"] == "completed"
    assert supervisor.subreaper() == before
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(10)"])
    try:
        other = tmp_path / "rejected"
        other.mkdir()
        assert launch(other, "pass")["status"] == "preexisting_children"
        assert supervisor.subreaper() == before
    finally:
        child.terminate()
        child.wait()
