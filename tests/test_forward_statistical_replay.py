"""Model-free synthetic fixtures; trust anchors are fixture-only, never study evidence."""

import copy
import builtins
import importlib.util
from pathlib import Path
import shutil
import sys

import pytest

TOOL = Path(__file__).resolve().parents[1] / "benchmarks/diagnostics/forward_replay"
SPEC = importlib.util.spec_from_file_location("forward_replay", TOOL / "replay.py")
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)
io = tool.io


def fixture_data():
    blobs = {n: io.read_bytes(p) for n, (p, _) in tool.fixed_files().items()}
    plan = io.parse(blobs["forward-analysis-plan.json"])
    protocol = {"expected_attempts": tool.SPLITS, "model": {"scope": "SYNTHETIC metadata; no model exists"}}
    blobs["audit-protocol.json"] = io.encoded(protocol)
    plan["audit_protocol_sha256"] = io.sha(blobs["audit-protocol.json"])
    plan["source_sha256"] = {n: io.sha(("synthetic source " + n).encode()) for n in plan["source_sha256"]}
    blobs["forward-analysis-plan.json"] = io.encoded(plan)
    attempts, occurrences = [], []
    for split, settings in plan["expected_grid"].items():
        methods = [m for spec in settings["groups"].values() for m in spec["methods"]]
        for method_index, method in enumerate(methods):
            for query_index, (query, target) in enumerate(settings["queries"].items()):
                for seed in [42, 43, 44]:
                    case = (query_index * 3 + seed - 42 + method_index) % 8
                    status = (
                        "timeout"
                        if case == 4
                        else "infrastructure_interruption"
                        if case == 5
                        else "no_candidates"
                        if case == 6
                        else "success"
                    )
                    attempt = {
                        "split": split,
                        "method_id": method,
                        "query_id": query,
                        "target_id": target,
                        "seed": seed,
                        "status": status,
                    }
                    attempts.append(attempt)
                    categories = (
                        ["no-route", "matched"]
                        if case == 0
                        else ["mismatched"]
                        if case == 1
                        else ["unknown"]
                        if case == 2
                        else ["failed"]
                        if case == 3
                        else ["matched"]
                        if case == 7
                        else []
                    )
                    for rank, category in enumerate(categories, 1):
                        row = {
                            "split": split,
                            "method_id": method,
                            "query_id": query,
                            "target_id": target,
                            "generation_seed": seed,
                            "rank": rank,
                            "smiles": "synthetic-product-" + category,
                            "status": "no_route_supplied"
                            if category == "no-route"
                            else "audit_failed"
                            if category == "failed"
                            else "assessed",
                            "forward_all_steps_matched": None,
                            "matched": False,
                            "evaluated": False,
                            "input_keys": [],
                            "route_key": None,
                        }
                        if category != "no-route":
                            row["route_key"] = io.sha(("synthetic-route-" + category).encode())
                            row["route_cache_key"] = row["route_key"]
                            row["all_steps_evaluated"] = None
                            row["error"] = "synthetic failure" if category == "failed" else None
                        if row["status"] == "assessed":
                            row["matched"] = row["forward_all_steps_matched"] = category == "matched"
                            row["evaluated"] = row["all_steps_evaluated"] = category != "unknown"
                            row["input_keys"] = [io.sha(b"synthetic shared context")]
                            if category == "matched":
                                row["input_keys"].append(io.sha(b"synthetic second shared context"))
                        occurrences.append(row)
    blobs["attempts.json"] = io.encoded(attempts)
    blobs["normalized-occurrences.json"] = io.encoded(occurrences)
    run_hashes = {
        f"comparison-{g}" if s == "scaffold" else f"comparison-{s}-{g}": io.sha((s + g).encode())
        for s, settings in plan["expected_grid"].items()
        for g in settings["groups"]
    }
    progress = {
        "complete": True,
        "observed_attempts": tool.SPLITS,
        "candidate_occurrences": len(occurrences),
        "status_counts": dict(tool.Counter(r["status"] for r in occurrences)),
        "input_run_sha256": run_hashes,
        "unique_cached_routes": len({r["route_key"] for r in occurrences if r["route_key"]}),
        "inference_usage_this_process": {
            "predictor_calls": 0,
            "cache_hits": 0,
            "input_tokens": 0,
            "output_tokens": 0,
        },
        "updated_utc": "synthetic fixture",
    }
    blobs["audit-progress.json"] = io.encoded(progress)
    source = io.parse(blobs[io.PROVENANCE])["code_sha256"]
    hashes = {"/synthetic/study/" + n: v for n, v in plan["source_sha256"].items()}
    hashes.update({"/synthetic/study/" + n + "/runs.json": v for n, v in run_hashes.items()})
    hashes.update({"/synthetic/source/" + n: v for n, v in source.items()})
    hashes.update(
        {
            "/synthetic/study/common-forward-audit-cpu/" + n + ".json": io.sha(blobs["audit-" + n + ".json"])
            for n in ("protocol", "progress")
        }
    )
    hashes["/synthetic/study/common-forward-audit-cpu/occurrences.json"] = io.sha(
        b"synthetic original occurrences"
    )
    validated = {
        "passed": True,
        "total_attempts": 630,
        "analysis_helper_sha256": tool.HELPER_SHA,
        "analysis_plan_sha256": io.sha(blobs["forward-analysis-plan.json"]),
        "scientific_comparison_emitted": True,
        "analysis_runtime": {
            "forward_inference_performed": False,
            "python": "synthetic",
            "numpy": "synthetic",
            "rdkit": "synthetic metadata only",
        },
        "audit_model_provenance": protocol["model"],
        "source_sha256": hashes,
        "resolved_inputs": {
            "source_root": "/synthetic/source",
            "frozen_code_file_count": 39,
            "frozen_code_inventory_sha256": tool.content_key(source),
            "model_directory": "/synthetic/nonexistent",
            "original_model_path_provenance": "/synthetic/nonexistent",
            "replay_assets": {},
        },
        "grid": {
            s: {"expected": n, "observed": n, "complete": True, "duplicates": [], "missing": []}
            for s, n in tool.SPLITS.items()
        },
    }
    blobs["recorded-validation.json"] = io.encoded(validated)
    integrity = {
        "passed": True,
        "attempts": 630,
        "candidate_occurrences": len(occurrences),
        "unique_routes": progress["unique_cached_routes"],
        "files_sha256_verified_before_after": dict(hashes),
        "outputs_sha256": {n: io.sha(blobs[n]) for n in ("attempts.json", "normalized-occurrences.json")},
        "scope": "SYNTHETIC fixture; no raw chemistry or scientific study was evaluated",
    }
    blobs["integrity-validation.json"] = io.encoded(integrity)
    stat = tool.module_bytes(
        blobs["statistics_extracted.py"], "synthetic_forward_statistics", TOOL / "statistics_extracted.py"
    )
    rows, outputs, targets, caps = tool.validate_normalized(attempts, occurrences, plan)
    result = {
        "analysis_plan": plan,
        "validation": validated,
        **{k: plan[k] for k in ("interpretation", "dependence", "conditional_metric", "infrastructure")},
        "splits": {
            s: stat.aggregate_split(rows[s], outputs[s], targets[s], caps[s], plan["bootstrap"])
            for s in tool.SPLITS
        },
    }
    blobs["recorded-analysis.json"] = io.encoded(result)
    integrity["files_sha256_verified_before_after"].update(
        {
            f"/synthetic/study/forward-cpu-runtime/final-analysis-20261016/{role}.json": io.sha(
                blobs[f"recorded-{role}.json"]
            )
            for role in ("validation", "analysis")
        }
    )
    blobs["integrity-validation.json"] = io.encoded(integrity)
    evidence = io.encoded({"scope": "SYNTHETIC test assessment only", "passed": True})
    blobs["parent-review-evidence/statistical-audit.json"] = evidence
    blobs["parent-review.json"] = io.encoded(
        {
            "schema_version": 1,
            "status": "passed",
            "total_attempts": 630,
            "expected_attempts": tool.SPLITS,
            "analysis_helper_sha256": tool.HELPER_SHA,
            "analysis_plan_sha256": io.sha(blobs["forward-analysis-plan.json"]),
            "extracted_statistics_sha256": tool.EXTRACTED_SHA,
            "original_source_commit": tool.COMMIT,
            "normalization_source_validation_sha256": io.sha(blobs["recorded-validation.json"]),
            "input_sha256": {n: io.sha(blobs[n]) for n in tool.DATA},
            "assessment_evidence": [{"artifact": "statistical-audit.json", "sha256": io.sha(evidence)}],
        }
    )
    return blobs


@pytest.fixture(scope="module")
def fixture_blobs():
    return fixture_data()


@pytest.fixture
def synthetic(tmp_path, monkeypatch, fixture_blobs):
    """Change only fixture plan trust anchor/path; never patch validation/statistics."""
    blobs = dict(fixture_blobs)
    plan_file = tmp_path / "synthetic-plan.json"
    plan_file.write_bytes(blobs["forward-analysis-plan.json"])
    monkeypatch.setattr(tool, "PLAN_SHA", io.sha(blobs["forward-analysis-plan.json"]))
    original = tool.fixed_files

    def fixture_files():
        files = original()
        files["forward-analysis-plan.json"] = (plan_file, tool.PLAN_SHA)
        return files

    monkeypatch.setattr(tool, "fixed_files", fixture_files)
    return blobs


def refresh(blobs):
    integrity = io.parse(blobs["integrity-validation.json"])
    integrity["outputs_sha256"] = {
        n: io.sha(blobs[n]) for n in ("attempts.json", "normalized-occurrences.json")
    }
    for role in ("validation", "analysis"):
        integrity["files_sha256_verified_before_after"][
            f"/synthetic/study/forward-cpu-runtime/final-analysis-20261016/{role}.json"
        ] = io.sha(blobs[f"recorded-{role}.json"])
    blobs["integrity-validation.json"] = io.encoded(integrity)
    review = io.parse(blobs["parent-review.json"])
    review["input_sha256"] = {n: io.sha(blobs[n]) for n in tool.DATA}
    blobs["parent-review.json"] = io.encoded(review)


def alter(blobs, name, callback):
    value = io.parse(blobs[name])
    callback(value)
    blobs[name] = io.encoded(value)
    refresh(blobs)


def write_source(root, blobs):
    for folder in ("inputs", "analysis", "audit", "review"):
        (root / folder).mkdir(parents=True)
    for name in ("attempts.json", "normalized-occurrences.json", "integrity-validation.json"):
        (root / "inputs" / name).write_bytes(blobs[name])
    for role in ("analysis", "validation"):
        (root / "analysis" / (role + ".json")).write_bytes(blobs["recorded-" + role + ".json"])
    for role in ("protocol", "progress"):
        (root / "audit" / (role + ".json")).write_bytes(blobs["audit-" + role + ".json"])
    (root / "review/parent.json").write_bytes(blobs["parent-review.json"])
    (root / "review/statistical-audit.json").write_bytes(
        blobs["parent-review-evidence/statistical-audit.json"]
    )


def test_pack_relocation_and_complete_analysis_replay(tmp_path, synthetic, monkeypatch):
    original_import = builtins.__import__

    def no_scientific_models(name, *args, **kwargs):
        assert name.split(".")[0] not in {"torch", "rdkit", "transformers", "chembounce_synthesis"}
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_scientific_models)
    source = tmp_path / "source"
    write_source(source, synthetic)
    packed = tool.pack(
        source / "inputs",
        source / "analysis",
        source / "audit",
        source / "review/parent.json",
        io.sha(synthetic["parent-review.json"]),
        tmp_path / "bundle",
    )
    relocated = tmp_path / "relocated"
    (tmp_path / "bundle").rename(relocated)
    shutil.rmtree(source)
    assert (
        tool.check(relocated, packed["manifest_sha256"], tmp_path / "verify.json", False)[
            "comparison_recomputed"
        ]
        is False
    )
    result = tool.check(relocated, packed["manifest_sha256"], tmp_path / "replayed", True)
    assert result["comparison_recomputed"] is True and result["attempts"] == 630
    assert result["integrity_recorded_analysis_sha256"] == {
        n: io.sha(synthetic[n]) for n in ("recorded-validation.json", "recorded-analysis.json")
    }
    for name, raw in synthetic.items():
        assert (relocated / name).read_bytes() == raw
    io.same(
        io.parse((tmp_path / "replayed/analysis.json").read_bytes()),
        io.parse(synthetic["recorded-analysis.json"]),
        tolerance=1e-12,
    )
    recorded = io.parse(synthetic["recorded-analysis.json"])["splits"]["scaffold"]
    coverage = recorded["methods"]["pipeline"]["coverage_counts"]
    assert coverage["no_route_supplied"] > 0 and coverage["route_audit_failed"] > 0
    assert recorded["dependence_diagnostics"]["input_keys_shared_across_targets"] > 0
    assert (
        recorded["complete_case_secondary_selection_biased"]["methods"]["pipeline"][
            "excluded_infrastructure_attempts"
        ]
        > 0
    )
    with pytest.raises(ValueError, match="Output exists"):
        tool.check(relocated, packed["manifest_sha256"], tmp_path / "replayed", True)
    with pytest.raises(ValueError, match="inside bundle"):
        tool.check(relocated, packed["manifest_sha256"], relocated / "extra", True)
    with pytest.raises(ValueError, match="external manifest SHA"):
        tool.check(relocated, "0" * 64, tmp_path / "wrong.json", False)


def test_extracted_statistics_equal_original_on_synthetic_630(fixture_blobs):
    before = list(sys.path)
    try:
        original = tool.module_bytes(
            fixture_blobs["analyze_forward_historical.py"],
            "historical_equivalence_only",
            TOOL / "analyze_forward_historical.py",
        )
    finally:
        sys.path[:] = before
    extracted = tool.module_bytes(
        fixture_blobs["statistics_extracted.py"], "extracted_equivalence", TOOL / "statistics_extracted.py"
    )
    plan = io.parse(fixture_blobs["forward-analysis-plan.json"])
    rows, outputs, targets, caps = tool.validate_normalized(
        io.parse(fixture_blobs["attempts.json"]), io.parse(fixture_blobs["normalized-occurrences.json"]), plan
    )
    for split in tool.SPLITS:
        args = rows[split], outputs[split], targets[split], caps[split], plan["bootstrap"]
        io.same(extracted.aggregate_split(*args), original.aggregate_split(*args), tolerance=1e-12)


def test_exact_query_seed_pairs_do_not_subtract_different_conditional_populations(fixture_blobs):
    stat = tool.module_bytes(
        fixture_blobs["statistics_extracted.py"], "pairing_test", TOOL / "statistics_extracted.py"
    )

    def row(seed, conditional):
        metrics = dict.fromkeys(stat.METRICS, 0)
        metrics[stat.CONDITIONAL] = conditional
        return {"query_id": "q", "target_id": "t", "seed": seed, "status": "success", "metrics": metrics}

    left, right = [row(42, None), row(43, 1.0)], [row(42, 1.0), row(43, None)]
    assert all(r["metrics"][stat.CONDITIONAL] is None for r in stat.paired_rows(left, right))


@pytest.mark.parametrize(
    "change",
    [
        "missing",
        "duplicate",
        "target",
        "seed-bool",
        "unknown-status",
        "noncontiguous",
        "duplicate-occurrence",
        "failed-attempt-output",
        "no-route-promoted",
        "matched-not-evaluated",
        "raw-flag",
        "bad-input-key",
        "bad-route-key",
        "inconsistent-reuse",
        "cap",
        "incomplete",
        "coverage",
        "run-hash",
        "integrity",
        "parent630",
        "parent-evidence",
        "analysis-value",
        "analysis-type",
    ],
)
def test_fail_closed_while_updated_parent_snapshot_hashes_do_not_bypass_rules(synthetic, change):
    blobs = synthetic
    if change in {"missing", "duplicate", "target", "seed-bool", "unknown-status", "failed-attempt-output"}:

        def modify(rows):
            if change == "missing":
                rows.pop()
            elif change == "duplicate":
                rows.append(rows[0])
            elif change == "target":
                rows[0]["target_id"] = "wrong"
            elif change == "seed-bool":
                rows[0]["seed"] = True
            elif change == "unknown-status":
                rows[0]["status"] = "running"
            else:
                rows[0]["status"] = "timeout"

        alter(blobs, "attempts.json", modify)
    elif change in {
        "noncontiguous",
        "duplicate-occurrence",
        "no-route-promoted",
        "matched-not-evaluated",
        "raw-flag",
        "bad-input-key",
        "bad-route-key",
        "inconsistent-reuse",
        "cap",
    }:

        def modify(rows):
            good = next(r for r in rows if r["matched"])
            if change == "noncontiguous":
                rows[0]["rank"] = 3
            elif change == "duplicate-occurrence":
                rows.append(rows[0])
            elif change == "no-route-promoted":
                rows[0]["matched"] = rows[0]["evaluated"] = True
            elif change == "matched-not-evaluated":
                good["evaluated"] = False
            elif change == "raw-flag":
                good["forward_all_steps_matched"] = False
            elif change == "bad-input-key":
                good["input_keys"] = ["bad"]
            elif change == "bad-route-key":
                good["route_key"] = "0" * 64
            elif change == "inconsistent-reuse":
                good["input_keys"] = [io.sha(b"changed-input")]
            else:
                extra = copy.deepcopy(rows[0])
                extra["rank"] = 3
                rows.append(extra)
                extra = copy.deepcopy(extra)
                extra["rank"] = 4
                rows.append(extra)

        alter(blobs, "normalized-occurrences.json", modify)
    elif change in {"incomplete", "coverage", "run-hash"}:

        def modify(value):
            if change == "incomplete":
                value["complete"] = False
            elif change == "coverage":
                value["candidate_occurrences"] += 1
            else:
                value["input_run_sha256"]["comparison-gpu0"] = "0" * 64

        alter(blobs, "audit-progress.json", modify)
    elif change == "integrity":
        alter(blobs, "integrity-validation.json", lambda v: v.update(passed=False))
    elif change in {"parent630", "parent-evidence"}:
        alter(
            blobs,
            "parent-review.json",
            lambda v: v.update(
                {"total_attempts": 576} if change == "parent630" else {"assessment_evidence": []}
            ),
        )
    else:

        def modify(value):
            result = value["splits"]["scaffold"]["methods"]["pipeline"]
            if change == "analysis-type":
                result["attempts"] = 36.0
            else:
                result["target_macro"]["route_supplied"]["mean"] += 0.01

        alter(blobs, "recorded-analysis.json", modify)
    with pytest.raises((ValueError, KeyError)):
        tool.validate(blobs, calculate=True)


def test_real_plan_and_archive_pins_cannot_be_replaced(fixture_blobs):
    with pytest.raises(ValueError, match="forward-analysis-plan"):
        tool.check_fixed(fixture_blobs)
    raw = {n: io.read_bytes(p) for n, (p, _) in tool.fixed_files().items()}
    tool.check_fixed(raw)
    for name in (
        "analyze_forward_historical.py",
        "statistics_extracted.py",
        "extraction.json",
        "transport_helpers.py",
    ):
        changed = dict(raw)
        changed[name] += b"\n"
        with pytest.raises(ValueError):
            tool.check_fixed(changed)


@pytest.mark.parametrize("role", ["validation", "analysis"])
def test_integrity_must_bind_each_original_recorded_file(synthetic, role):
    integrity = io.parse(synthetic["integrity-validation.json"])
    integrity["files_sha256_verified_before_after"][
        f"/synthetic/study/forward-cpu-runtime/final-analysis-20261016/{role}.json"
    ] = "0" * 64
    synthetic["integrity-validation.json"] = io.encoded(integrity)
    review = io.parse(synthetic["parent-review.json"])
    review["input_sha256"] = {n: io.sha(synthetic[n]) for n in tool.DATA}
    synthetic["parent-review.json"] = io.encoded(review)
    with pytest.raises(ValueError, match="normalization/original recorded " + role + " hash"):
        tool.validate(synthetic, calculate=False)


@pytest.mark.parametrize("raw", [b'{"x":1,"x":1}', b'{"x":NaN}', b'{"x":1e309}', b'{"x":Infinity}'])
def test_strict_json(raw):
    with pytest.raises(ValueError):
        io.parse(raw)


def test_no_statistics_start_before_full_grid_validation(synthetic, monkeypatch):
    alter(synthetic, "attempts.json", lambda rows: rows.pop())

    def forbidden(*args, **kwargs):
        pytest.fail("Statistical code executed before completion gate")

    monkeypatch.setattr(tool, "module_bytes", forbidden)
    with pytest.raises(ValueError, match="Incomplete"):
        tool.validate(synthetic, calculate=True)


def test_changed_input_during_pack_and_output_preservation(tmp_path, synthetic, monkeypatch):
    source = tmp_path / "source"
    write_source(source, synthetic)
    original = tool.validate

    def changed(blobs, calculate=False):
        result = original(blobs, calculate=calculate)
        with (source / "inputs/attempts.json").open("ab") as stream:
            stream.write(b" ")
        return result

    monkeypatch.setattr(tool, "validate", changed)
    args = (
        source / "inputs",
        source / "analysis",
        source / "audit",
        source / "review/parent.json",
        io.sha(synthetic["parent-review.json"]),
    )
    with pytest.raises(ValueError, match="input changed"):
        tool.pack(*args, tmp_path / "bundle")
    assert not (tmp_path / "bundle").exists()
    (tmp_path / "existing").mkdir()
    (tmp_path / "existing/sentinel").write_text("untouched")
    with pytest.raises(ValueError, match="already exists"):
        tool.pack(*args, tmp_path / "existing")
    assert (tmp_path / "existing/sentinel").read_text() == "untouched"


def write_bundle(path, blobs):
    path.mkdir()
    for name, raw in blobs.items():
        p = path / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(raw)
    raw = io.encoded(tool.manifest_for(blobs))
    (path / "manifest.json").write_bytes(raw)
    return io.sha(raw)


@pytest.mark.parametrize("change", ["symlink", "escape", "tamper", "extra", "duplicate-manifest-key"])
def test_bundle_storage_boundaries(tmp_path, synthetic, change):
    path = tmp_path / "bundle"
    digest = write_bundle(path, synthetic)
    if change == "symlink":
        (path / "attempts.json").rename(tmp_path / "outside.json")
        (path / "attempts.json").symlink_to(tmp_path / "outside.json")
    elif change == "escape":
        manifest = io.parse((path / "manifest.json").read_bytes())
        manifest["files"]["../outside.json"] = manifest["files"].pop("attempts.json")
        raw = io.encoded(manifest)
        (path / "manifest.json").write_bytes(raw)
        digest = io.sha(raw)
    elif change == "tamper":
        (path / "attempts.json").write_bytes(b"[]")
    elif change == "extra":
        (path / "unexpected.json").write_text("{}")
    else:
        raw = (
            (path / "manifest.json")
            .read_bytes()
            .replace(b'"schema_version": 1,', b'"schema_version": 1, "schema_version": 1,')
        )
        (path / "manifest.json").write_bytes(raw)
        digest = io.sha(raw)
    with pytest.raises(ValueError):
        tool.check(path, digest, tmp_path / "verified.json", False)
