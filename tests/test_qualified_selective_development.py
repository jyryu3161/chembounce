"""Synthetic lifecycle/denominator tests and one known tiny chemistry fixture."""

import copy
import importlib
import os
import signal
import time
from types import SimpleNamespace

import pytest

from benchmarks.research import qualified_selective_development as subject


def fixture_rows():
    pairs = []
    for index, (target, query, label) in enumerate(
        (
            ("A", "q1", "active"),
            ("A", "q1", "active"),
            ("A", "q2", "active"),
            ("B", "q3", "active"),
            ("B", "q3", "inactive"),
        )
    ):
        pairs.append(
            dict(
                pair_id=str(index),
                target_id=target,
                assay_id="assay",
                context_id="context",
                query_id=query,
                query_smiles="CCO",
                control_smiles="CCO",
                label=label,
                topology_change=True,
            )
        )
    old = [
        {
            **p,
            "hit_tanimoto": 1.0,
            "trials": [
                {"seed": s, "status": "success", "info": {"hit": {"shape": 1.0, "pharmacophore": 1.0}}}
                for s in subject.SEEDS
            ],
        }
        for p in pairs
    ]
    return {"pairs": pairs, "slots": subject.ordered_slots(pairs)}, old


class Fake:
    def __init__(self, output, mode=None):
        self.output, self.mode = output, mode
        self.references, self.trials, self.finished = [], [], []
        self.closed = False

    def start(self):
        ledger = subject.read(self.output / "ledger.json")
        assert len(ledger["trials"]) == 15
        assert all(r["status"] == "not_started" for r in ledger["trials"])

    def reference(self, pair, seed, deadline):
        self.references.append((pair["query_id"], seed))
        if self.mode == "reference_failure":
            raise ValueError("unavailable reference")
        return seed, {"seed": seed}

    def finish_reference(self, ref):
        self.finished.append(ref)

    def trial(self, pair, seed, old, old_pair, ref, deadline):
        self.trials.append((pair["pair_id"], seed))
        if self.mode == "signal":
            os.kill(os.getpid(), signal.SIGTERM)
        if self.mode == "exception":
            raise ValueError("known failure")
        if self.mode == "rebind":
            return {"status": "success", "slot": {"index": -1}}
        return {
            "status": "abstained" if pair["pair_id"] == "0" else "success",
            "old_comparison": {"status": "assessed", "passed": False},
        }

    def costs(self):
        return {"calls": len(self.trials)}

    def close(self):
        self.closed = True


def execute(tmp_path, mode=None, clock=time.monotonic, wall=30):
    plan, old = fixture_rows()
    output = tmp_path / "out"
    backend = Fake(output, mode)
    result = subject._execute(plan, old, output, backend, lambda: None, wall, clock)
    return plan, result, backend, subject.read(output / "summary.json")


def test_frozen_slot_order_and_duplicate_rejection():
    plan, _ = fixture_rows()
    assert [s["seed"] for s in plan["slots"][:6]] == [42, 42, 43, 43, 44, 44]
    with pytest.raises(ValueError, match="Duplicate"):
        subject.ordered_slots(plan["pairs"] * 2)
    plan["pairs"][1]["query_smiles"] = "CO"
    with pytest.raises(ValueError, match="conflicting"):
        subject.ordered_slots(plan["pairs"])


def test_lifecycle_complete_and_macro_full_denominators(tmp_path):
    _, _, backend, summary = execute(tmp_path)
    assert backend.closed and len(backend.references) == 9 and len(backend.finished) == 9
    assert len(backend.trials) == 15
    assert summary["run_status"] == "complete"
    rate = summary["rates"]["new_joint"]["active_acceptance"]
    assert rate["macro_rational"] == "7/8"
    assert rate["targets"]["A"]["rational"] == "3/4"
    assert rate["denominator"] == 4
    ledger = subject.read(tmp_path / "out/ledger.json")
    assert all(r["reference_artifact"] for r in ledger["trials"])
    assert all(
        subject.base.file_hash(r["artifact"]["path"]) == r["artifact"]["sha256"] for r in ledger["trials"]
    )


@pytest.mark.parametrize("mode", ["reference_failure", "exception", "rebind"])
def test_all_failures_retained_no_retries(tmp_path, mode):
    _, _, backend, summary = execute(tmp_path, mode)
    assert summary["trial_status_counts"] == {"failed": 15}
    assert summary["rates"]["new_joint"]["active_acceptance"]["denominator"] == 4
    assert summary["rates"]["new_joint"]["active_acceptance"]["macro"] == 0
    assert summary["old_joint_unassessed_pairs"] == 5
    assert len(backend.references) == 9
    assert len(backend.trials) == (0 if mode == "reference_failure" else 15)


def test_signal_retains_remaining_slots(tmp_path):
    _, _, backend, summary = execute(tmp_path, "signal")
    assert summary["run_status"] == "interrupted"
    assert sum(summary["trial_status_counts"].values()) == 15
    assert summary["trial_status_counts"]["not_started"] == 14
    assert backend.closed


def test_deadline_before_embedding_keeps_zeros(tmp_path):
    ticks = iter(range(1000))
    _, _, backend, summary = execute(tmp_path, clock=lambda: next(ticks), wall=1)
    assert summary["run_status"] == "deadline"
    assert summary["trial_status_counts"] == {"not_started": 15}
    assert not backend.references and not backend.trials


def test_source_failure_invalidates_interpretation(tmp_path):
    plan, old = fixture_rows()
    backend = Fake(tmp_path / "out")
    checks = []

    def verify():
        checks.append(1)
        if len(checks) > 1:
            raise ValueError("source changed")

    with pytest.raises(ValueError, match="source changed"):
        subject._execute(plan, old, tmp_path / "out", backend, verify, 30)
    report = subject.read(tmp_path / "out/summary.json")
    assert report["run_status"] == "integrity_or_cleanup_failed"
    assert report["valid_for_interpretation"] is False
    assert backend.closed


def test_output_never_reused(tmp_path):
    plan, old = fixture_rows()
    with pytest.raises(FileExistsError):
        subject._execute(plan, old, tmp_path, Fake(tmp_path), lambda: None, 30)


@pytest.mark.parametrize("change", ["missing", "extra", "rebound"])
def test_summary_slot_identity_fail_closed(change):
    plan, old = fixture_rows()
    rows = {s["index"]: {"slot": copy.deepcopy(s), "status": "not_started"} for s in plan["slots"]}
    if change == "missing":
        rows.pop(0)
    elif change == "extra":
        rows[100] = rows[0]
    else:
        rows[0]["slot"]["pair_id"] = "elsewhere"
    with pytest.raises(ValueError):
        subject.summarize(plan["pairs"], plan["slots"], rows, old)


def test_model_import_and_swallowed_network_fail_closed():
    with pytest.raises(RuntimeError, match="Forbidden operation"):
        with subject.cpu_only():
            with pytest.raises(RuntimeError, match="Model import"):
                importlib.import_module("transformers")
    import socket

    with pytest.raises(RuntimeError, match="Forbidden operation"):
        with subject.cpu_only():
            with socket.socket() as sock:
                with pytest.raises(RuntimeError):
                    sock.bind(("127.0.0.1", 0))


def test_measured_embedding_calls_are_not_cache_misses():
    backend = subject.ScientificBackend()
    backend.original_embed = lambda *args, **kw: [0, 1]
    backend.role = "reference_prewarm"
    assert backend._embed(None, numConfs=8, params=SimpleNamespace(numThreads=1)) == [0, 1]
    assert backend.costs()["actual_embedding_calls"] == 1
    assert backend.embeddings[0]["returned_conformers"] == 2
    assert backend.requests == []
    with pytest.raises(ValueError, match="bound"):
        backend._embed(None, numConfs=9, params=SimpleNamespace(numThreads=1))


def test_measured_embedding_exception_retains_unknown_returned_count():
    backend = subject.ScientificBackend()

    def fail(*args, **kw):
        raise ValueError("embed failed")

    backend.original_embed = fail
    with pytest.raises(ValueError):
        backend._embed(None, numConfs=8, params=SimpleNamespace(numThreads=1))
    assert backend.embeddings[0]["status"] == "failed"
    assert backend.embeddings[0]["returned_conformers"] is None


def test_known_tiny_fixture_reference_cache_and_saved_old_pose():
    from chembounce_synthesis import conformers

    if conformers.rdShapeAlign is None:
        pytest.skip("RDKit alignment capability unavailable")
    backend = subject.ScientificBackend()
    pair = {"query_smiles": "CCO", "control_smiles": "CCO"}
    backend.start()
    try:
        ref, _ = backend.reference(pair, 42, time.monotonic() + 30)
        old_result = conformers.score_candidate("CCO", "CCO", "CCO", ref["settings"])
        info, hit, design = old_result
        cid = info["candidate_conformer_id"]
        old = {
            "status": "success",
            "info": info,
            "poses": {
                "hit": subject.base.pose_record(hit, cid),
                "design": subject.base.pose_record(design, cid),
            },
        }
        result = backend.trial(pair, 42, old, {"hit_tanimoto": 1.0}, ref, time.monotonic() + 30)
        assert result["status"] == "success", result
        assert result["old_comparison"]["status"] == "assessed"
        assert result["old_comparison"]["passed"]
        assert "regenerated_reference" in result["old_comparison"]["scope"]
        assert backend.costs()["actual_embedding_calls"] == 1
        assert sum(r["cache_misses"] for r in backend.requests) == 1
        backend.finish_reference(ref)
        restored = (
            result["info"],
            subject.base.restore_pose(result["poses"]["hit"]),
            subject.base.restore_pose(result["poses"]["design"]),
        )
        backend.scorer.validate_result(
            restored, "CCO", ref["settings"], hit="CCO", design="CCO", similarity=ref["policy"]
        )
        count = backend.costs()["actual_embedding_calls"]
        failure = backend.trial(
            {**pair, "control_smiles": "CC(O)C(=O)O"},
            42,
            {"status": "failed"},
            {},
            ref,
            time.monotonic() + 30,
        )
        assert failure["status"] == "failed"
        assert backend.costs()["actual_embedding_calls"] == count
        abstained = backend.trial(
            {**pair, "control_smiles": "c1ccccc1"}, 42, {"status": "failed"}, {}, ref, time.monotonic() + 30
        )
        assert abstained["status"] == "abstained"
        assert abstained["audit"]["status"] == "abstained"
        assert abstained["audit"]["ensemble_requests"] == []
        assert backend.costs()["actual_embedding_calls"] == count
        bad = copy.deepcopy(old)
        bad["info"]["ensemble_diagnostics"]["hit"]["embedded"] += 1
        unavailable = backend.old_comparison(pair, bad, {"hit_tanimoto": 1.0}, ref)
        assert unavailable["status"] == "unavailable"
        ref["ensemble"].relative_energies[0] += 1
        with pytest.raises(ValueError, match="mutated"):
            backend.finish_reference(ref)
    finally:
        backend.close()
    assert conformers.ensemble.cache_info().currsize == 0


def declaration_fixture(tmp_path, monkeypatch):
    from benchmarks.research import assay_matched

    plan, old = fixture_rows()
    paths = {}
    payloads = {key: {} for key in subject.INPUTS}
    runtime = {"python": "3.10.12", "packages": {"rdkit": "2025.9.6", "numpy": "1.26.4", "scipy": "1.13.1"}}
    sources = {
        "chembounce_synthesis/" + name: "synthetic_hash"
        for name in ("conformers.py", "common.py", "identity.py")
    }
    payloads.update(
        pairs=plan["pairs"],
        old_scores=old,
        old_protocol={
            "runtime": runtime,
            "source_sha256": sources,
            "calibration_spec": {"conformers": subject.base.SETTINGS, "seeds": list(subject.SEEDS)},
        },
    )
    for key, value in payloads.items():
        path = tmp_path / (key + ".json")
        subject.base.save(path, value)
        paths[key] = {"path": str(path), "sha256": subject.base.file_hash(path)}
    monkeypatch.setattr(subject, "INPUTS", {k: v["sha256"] for k, v in paths.items()})
    monkeypatch.setattr(subject, "COUNTS", subject.base.counts(plan["pairs"]))
    monkeypatch.setattr(subject, "sources", lambda: sources)
    monkeypatch.setattr(subject.base, "runtime", lambda: runtime)
    monkeypatch.setattr(subject.base, "verify_inventory", lambda pairs, manifest: None)
    monkeypatch.setattr(assay_matched, "control_pairs", lambda panel: plan["pairs"])
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        monkeypatch.setenv(key, "1")
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setenv("PYTHONDONTWRITEBYTECODE", "1")
    declaration = {
        "schema": "qualified-selective-development-declaration-v2",
        "spec": subject.SPEC,
        "inputs": paths,
        "source_sha256": sources,
        "runtime": runtime,
    }
    path = tmp_path / "declaration.json"
    subject.base.save(path, declaration)
    return path, declaration


def test_prepare_run_exact_plan_and_before_after_bindings(tmp_path, monkeypatch):
    path, _ = declaration_fixture(tmp_path, monkeypatch)
    plan = subject.prepare(path, subject.base.file_hash(path), tmp_path / "prepared")
    assert len(plan["slots"]) == 15
    plan_path = tmp_path / "prepared/plan.json"
    backend = Fake(tmp_path / "execution")
    monkeypatch.setattr(subject, "ScientificBackend", lambda: backend)
    result = subject.run(plan_path, subject.base.file_hash(plan_path), tmp_path / "execution")
    assert len(result) == 15
    integrity = subject.read(tmp_path / "execution/integrity.json")
    assert len(integrity["checks"]) == 4
    assert integrity["forbidden_attempts"] == []


@pytest.mark.parametrize("field", ["input", "source", "runtime", "spec", "cpu"])
def test_declaration_drift_rejected(tmp_path, monkeypatch, field):
    _, declaration = declaration_fixture(tmp_path, monkeypatch)
    declaration = copy.deepcopy(declaration)
    if field == "input":
        open(declaration["inputs"]["panel"]["path"], "a").write(" ")
    elif field == "source":
        declaration["source_sha256"]["extra"] = "unbound"
    elif field == "runtime":
        declaration["runtime"]["packages"]["rdkit"] = "other"
    elif field == "spec":
        declaration["spec"]["minimum_polar_reference_recall"] = 0.5
    else:
        monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    with pytest.raises(ValueError):
        subject.verify_declaration(declaration)


@pytest.mark.parametrize("field", ["slots", "pairs", "declaration_bytes", "external_digest"])
def test_plan_tampering_fails_before_backend(tmp_path, monkeypatch, field):
    path, _ = declaration_fixture(tmp_path, monkeypatch)
    plan = subject.prepare(path, subject.base.file_hash(path), tmp_path / "prepared")
    plan_path = tmp_path / "prepared/plan.json"
    if field == "slots":
        plan["slots"].pop()
    elif field == "pairs":
        plan["pairs"][0]["label"] = "inactive"
    elif field == "declaration_bytes":
        with open(path, "a") as stream:
            stream.write(" ")
    subject.base.save(plan_path, plan)
    expected = "invalid" if field == "external_digest" else subject.base.file_hash(plan_path)
    with pytest.raises(ValueError):
        subject.run(plan_path, expected, tmp_path / "execution")
    assert not (tmp_path / "execution").exists()
