"""Portable execution retains assay-local denominators, seed errors and inputs."""

import concurrent.futures
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from benchmarks.research import assay_control_scoring as scoring
from benchmarks.research import assay_matched as matched
from benchmarks.research.protocol import digest, write_json
from chembounce_synthesis import conformers
from test_assay_matched import write_sources


@pytest.fixture
def inputs(tmp_path):
    # Reuse synthetic source provenance only; these labels are not evidence.
    pool = matched.curate_sources(write_sources(tmp_path / "sources"))
    panel = matched.select_panel(
        pool,
        matched.Selection(
            min_active_controls=1,
            min_inactive_controls=1,
            min_active_scaffolds=1,
            max_controls_per_label=1,
            max_queries_per_context=1,
            max_contexts_per_target=1,
        ),
    )
    panel_path, settings_path = tmp_path / "panel.json", tmp_path / "conformers.json"
    write_json(panel_path, panel)
    write_json(settings_path, {"num_conformers": 2, "max_iterations": 500})
    return panel_path, settings_path, panel


def read(output, name):
    return json.loads((output / name).read_text())


def test_actual_common_scorer_retains_2d_after_seed_failure_and_full_denominators(
    inputs, tmp_path, monkeypatch, capsys
):
    panel_path, settings_path, panel = inputs
    pairs = list(matched.control_pairs(panel))
    failed_product = pairs[0]["control_smiles"]

    def candidate(product, hit, design, settings, hit_pose):
        if product == failed_product and settings.seed == 43:
            raise ValueError("Deliberate conformer failure")
        return {"score": 0.5, "hit": {"shape": 0.4, "pharmacophore": 0.6}}, None, None

    monkeypatch.setattr(conformers, "score_candidate", candidate)
    output = tmp_path / "scores"
    rows = scoring.run_scoring(panel_path, output, settings_path, [42, 43], 1, resamples=20)
    assert [r["pair_id"] for r in rows] == [p["pair_id"] for p in pairs]
    assert read(output, "scores.json") == rows
    assert all(row["status"] == "success" for row in rows)
    failed = [row for row in rows if not row["all_3d_seeds_succeeded"]]
    assert failed
    for row in failed:
        assert set(row["scores"]) == {"bit_tanimoto", "count_tanimoto"}
        assert [t["seed"] for t in row["trials"]] == [42, 43]
        assert [t["status"] for t in row["trials"]] == ["success", "failed"]
        assert row["trials"][1]["reason"] == "Deliberate conformer failure"
    full, paired = read(output, "framework-summary.json"), read(output, "framework-paired-summary.json")
    assert full["bit_tanimoto"]["missing_or_failed_pairs"] == 0
    assert full["score"]["missing_or_failed_pairs"] == len(failed)
    assert paired["bit_tanimoto"]["missing_or_failed_pairs"] == len(failed)
    for report in [*full.values(), *paired.values()]:
        assert report["pair_count"] == len(pairs) and report["unreturned_pairs"] == 0
        assert len(report["targets"]) == 2
    protocol = read(output, "protocol.json")
    assert protocol["panel_sha256"] == digest(panel_path) == digest(output / "panel.json")
    assert protocol["conformers_sha256"] == digest(settings_path) == digest(output / "conformers.json")
    assert protocol["pairs_sha256"] == digest(output / "control_pairs.json")
    assert protocol["settings_sha256"] == scoring._json_hash(protocol["effective_settings"])
    assert "seed" not in protocol["effective_settings"]["conformers"]
    assert protocol["source_payloads_verified"] is True
    assert protocol["source_sha256"] == {s["id"]: s["sha256"] for s in panel["sources"]}
    root = Path(__file__).resolve().parents[1]
    assert all(digest(root / name) == value for name, value in protocol["scoring_code_sha256"].items())
    updates = [json.loads(line) for line in capsys.readouterr().err.splitlines() if line.startswith("{")]
    assert updates[0]["scored_pairs"] == 0 and not updates[0]["complete"]
    assert {u["scored_pairs"] for u in updates} == set(range(len(pairs) + 1))
    assert updates[-1]["complete"] and updates[-1]["failed_3d_pairs"] == len(failed)


def test_worker_execution_errors_do_not_disappear(inputs, tmp_path, monkeypatch):
    class BrokenExecutor:
        def __init__(self, **kwargs):
            assert kwargs["mp_context"].get_start_method() == "spawn"

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def submit(self, function, payload):
            future = concurrent.futures.Future()
            future.set_exception(RuntimeError("Worker crashed"))
            return future

    monkeypatch.setattr(scoring.concurrent.futures, "ProcessPoolExecutor", BrokenExecutor)
    panel_path, settings_path, panel = inputs
    output = tmp_path / "scores"
    rows = scoring.run_scoring(panel_path, output, settings_path, [42, 43], 3, resamples=10)
    assert len(rows) == len(list(matched.control_pairs(panel)))
    assert all(r["status"] == "error" and r["scores"] == {} for r in rows)
    assert all(r["error_stage"] == "worker_execution" for r in rows)
    assert all([t["seed"] for t in r["trials"]] == [42, 43] for r in rows)
    assert all(t["status"] == "unreported" for r in rows for t in r["trials"])
    progress = read(output, "progress.json")
    assert progress["complete"] and progress["unexpected_errors"] == len(rows)
    for report in read(output, "framework-summary.json").values():
        assert report["missing_or_failed_pairs"] == len(rows)
        assert report["target_macro"]["roc_auc"] is None
        assert report["target_macro"]["roc_auc_bounds"] == [0.0, 1.0]


def test_worker_initialization_and_submission_failures_preserve_pairs(monkeypatch):
    pairs = [{"pair_id": "a"}, {"pair_id": "b"}]

    def broken_init(**kwargs):
        raise OSError("Cannot start workers")

    monkeypatch.setattr(scoring.concurrent.futures, "ProcessPoolExecutor", broken_init)
    rows = list(scoring._pair_results(pairs, {}, [42], 3))
    assert [r["pair_id"] for r in rows] == ["a", "b"]
    assert all(r["error_stage"] == "worker_initialization" for r in rows)

    class BrokenSubmission:
        def __init__(self, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def submit(self, *args):
            raise RuntimeError("Cannot submit")

    monkeypatch.setattr(scoring.concurrent.futures, "ProcessPoolExecutor", BrokenSubmission)
    rows = list(scoring._pair_results(pairs, {}, [42], 3))
    assert [r["pair_id"] for r in rows] == ["a", "b"]
    assert all(r["error_stage"] == "worker_submission" for r in rows)


def test_uncaught_scorer_error_records_unknown_seed_outcomes(monkeypatch):
    def fail(*args):
        raise KeyError("unexpected scoring error")

    monkeypatch.setattr(scoring, "similarity_scores", fail)
    row = scoring.score_pair(({"query_smiles": "C", "control_smiles": "CC"}, {}, [42, 43]))
    assert row["status"] == "error" and not row["all_3d_seeds_succeeded"]
    assert row["scores"] == {} and "KeyError" in row["error"]
    assert [t["seed"] for t in row["trials"]] == [42, 43]
    assert all(t["status"] == "unreported" for t in row["trials"])


def test_missing_seed_never_exports_partial_3d_aggregate(monkeypatch):
    monkeypatch.setattr(
        scoring,
        "similarity_scores",
        lambda *args: (
            {"bit_tanimoto": 0.4, "score": 0.8, "score_stability": {}},
            [{"seed": 42, "status": "success"}],
        ),
    )
    row = scoring.score_pair(({"query_smiles": "C", "control_smiles": "CC"}, {}, [42, 43]))
    assert row["status"] == "success" and not row["all_3d_seeds_succeeded"]
    assert row["scores"] == {"bit_tanimoto": 0.4}
    assert row["trials"][1]["seed"] == 43 and row["trials"][1]["status"] == "unreported"


@pytest.mark.parametrize(
    "options",
    [
        {"seeds": []},
        {"seeds": [42, 42]},
        {"seeds": [True]},
        {"seeds": [-1]},
        {"seeds": [2**31]},
        {"workers": 0},
        {"workers": True},
        {"resamples": 0},
        {"bootstrap_seed": -1},
        {"bootstrap_seed": True},
    ],
)
def test_invalid_runtime_settings_rejected_before_creating_output(inputs, tmp_path, options):
    panel_path, settings_path, _ = inputs
    output = tmp_path / "scores"
    with pytest.raises(ValueError):
        scoring.run_scoring(panel_path, output, settings_path, **options)
    assert not output.exists()


@pytest.mark.parametrize(
    "settings", [[], {"num_conformers": 0}, {"unknown": 1}, {"energy_window": float("nan")}]
)
def test_invalid_conformer_settings_rejected_before_scoring(inputs, tmp_path, settings):
    panel_path, settings_path, _ = inputs
    settings_path.write_text(json.dumps(settings))
    output = tmp_path / "scores"
    with pytest.raises((TypeError, ValueError)):
        scoring.run_scoring(panel_path, output, settings_path)
    assert not output.exists()


def test_source_hash_verification_is_mandatory(inputs, tmp_path):
    panel_path, settings_path, panel = inputs
    Path(panel["sources"][0]["artifact"]).write_text("Changed source")
    output = tmp_path / "scores"
    with pytest.raises(ValueError, match="Source hash mismatch"):
        scoring.run_scoring(panel_path, output, settings_path)
    assert not output.exists()


def test_duplicate_pair_ids_rejected_before_output(inputs, tmp_path, monkeypatch):
    panel_path, settings_path, panel = inputs
    pair = next(matched.control_pairs(panel))
    monkeypatch.setattr(scoring, "control_pairs", lambda panel: iter([pair, pair]))
    output = tmp_path / "scores"
    with pytest.raises(ValueError, match="duplicate pair IDs"):
        scoring.run_scoring(panel_path, output, settings_path)
    assert not output.exists()


@pytest.mark.parametrize("kind", ["file", "directory", "dangling_symlink"])
def test_existing_output_is_never_overwritten(tmp_path, kind):
    output = tmp_path / "scores"
    if kind == "file":
        output.write_text("Retain original")
    elif kind == "directory":
        output.mkdir()
        (output / "record").write_text("Retain original")
    else:
        output.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError, match="Output already exists"):
        scoring.run_scoring(tmp_path / "missing-panel", output, tmp_path / "missing-settings")
    if kind == "file":
        assert output.read_text() == "Retain original"
    elif kind == "directory":
        assert (output / "record").read_text() == "Retain original"
    else:
        assert output.is_symlink()


def test_portable_cli_actual_small_scorer_with_spawn_workers(inputs, tmp_path):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    panel_path, settings_path, panel = inputs
    output = tmp_path / "spawn-smoke"
    root = Path(__file__).resolve().parents[1]
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "benchmarks.research.assay_control_scoring",
            "--panel",
            str(panel_path),
            "--output",
            str(output),
            "--conformers",
            str(settings_path),
            "--seeds",
            "42",
            "--workers",
            "2",
            "--resamples",
            "20",
        ],
        cwd=root,
        env={**os.environ, "PYTHONPATH": str(root)},
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stderr
    rows = read(output, "scores.json")
    assert [r["pair_id"] for r in rows] == [p["pair_id"] for p in matched.control_pairs(panel)]
    assert all(row["status"] == "success" and row["all_3d_seeds_succeeded"] for row in rows)
    assert all(set(scoring.METRICS).issubset(row["scores"]) for row in rows)
    assert read(output, "progress.json")["complete"]
    assert read(output, "framework-summary.json")["score"]["missing_or_failed_pairs"] == 0
