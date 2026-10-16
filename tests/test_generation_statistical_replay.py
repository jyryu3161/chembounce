"""Synthetic 630-slot records exercise storage/statistics; no chemistry is executed."""

import copy
import importlib.util
import json
from pathlib import Path
import shutil
import statistics

import pytest

TOOL = Path(__file__).resolve().parents[1] / "benchmarks/diagnostics/generation_replay"
SPEC = importlib.util.spec_from_file_location("generation_replay", TOOL / "replay.py")
replay = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(replay)


def make_fixture():
    blobs = {n: (TOOL / n).read_bytes() for n in (replay.HELPER, replay.ATTRIBUTION, replay.PROVENANCE)}
    helper = replay.historical(blobs[replay.HELPER])
    code = replay.parse(blobs[replay.PROVENANCE])["code_sha256"]
    sources = replay.parse(blobs[replay.ATTRIBUTION])["sources"]
    source_inventory = [
        {
            "id": source["name"],
            "sha256": source["sha256"],
            "url": source["url"],
            "artifact": f"/synthetic/original-study/sources/chembl/{source['name']}.json",
        }
        for source in sources
    ]
    for split, (_, query_count) in replay.SPLITS.items():
        names = replay.paths(split)
        queries = [{"id": f"q{i}", "compound_id": f"c{i}"} for i in range(query_count)]
        panel = {
            "split": {"kind": split},
            "queries": queries,
            "compounds": [
                {"id": f"c{i}", "target_id": f"synthetic-target-{i // 2}", "smiles": f"synthetic-query-{i}"}
                for i in range(query_count)
            ],
            "sources": source_inventory,
        }
        blobs[names["panel"]] = replay.encoded(panel)
        all_runs, input_hashes = [], {}
        for group in replay.GROUPS:
            mids = replay.METHODS[group]
            if split != "scaffold" and group != "baselines":
                mids = mids[:1]
            protocol = {
                "seeds": replay.SEEDS,
                "methods": [{"id": m} for m in mids],
                "budget": {
                    "max_candidates": 3,
                    "max_evaluations": 1000000,
                    "max_memory_mb": 24000,
                    "wall_seconds": 360.0,
                },
            }
            blobs[names[f"protocol_{group}"]] = replay.encoded(protocol)
            normalized = copy.deepcopy(protocol)
            for method in normalized["methods"]:
                method["resolved_assets"] = {"unused_recorded_path": "/private/metadata/only"}
            manifest = {
                "code_sha256": code,
                "protocol": normalized,
                "protocol_sha256": replay.sha(blobs[names[f"protocol_{group}"]]),
                "dataset": {
                    "dataset_sha256": replay.sha(blobs[names["panel"]]),
                    "source_artifacts": {s["id"]: s["artifact"] for s in source_inventory},
                },
                "versions": dict.fromkeys(["numpy", "rdkit", "scipy", "psutil"], "synthetic-fixture"),
                "python": "synthetic-fixture",
                "platform": "synthetic-fixture",
            }
            blobs[names[f"manifest_{group}"]] = replay.encoded(manifest)
            runs = []
            for method in mids:
                for q in range(query_count):
                    for seed in replay.SEEDS:
                        case = (q * 3 + seed - 42) % 6
                        status = [
                            "success",
                            "timeout",
                            "infrastructure_interruption",
                            "success",
                            "success",
                            "no_candidates",
                        ][case]
                        count = 2 if case == 0 else 1 if case == 3 else 0
                        candidates = [
                            {
                                "smiles": f"synthetic-candidate-{q}-{j}",
                                "family": f"family-{q}-{j}",
                                # Raw rank three is valid-output ordinal two: test this exact distinction.
                                "rank": 1 if j == 0 else 3,
                            }
                            for j in range(count)
                        ]
                        metrics = {m: 0 for m in helper.EMPTY_OUTPUT_METRICS}
                        metrics.update(
                            run_success=int(bool(count)),
                            candidate_count=3 if count == 2 else count,
                            valid_candidate_count=count,
                            valid_fraction=count / 3 if count == 2 else float(bool(count)),
                            unique_families=count,
                            hop_families=count,
                            verified_rule_routes=int(bool(count)),
                            reference_active_scaffolds=2,
                            active_scaffold_recall_proxy=0.0,
                        )
                        resources = {
                            "wall_seconds": None
                            if case == 2
                            else 1.5 + q / 10 + (0.1 if method != "pipeline" else 0)
                        }
                        if case == 2:
                            resources["measurement_status"] = "lost_with_runner"
                        runs.append(
                            {
                                "method_id": method,
                                "query_id": f"q{q}",
                                "target_id": f"synthetic-target-{q // 2}",
                                "seed": seed,
                                "status": status,
                                "resources": resources,
                                "evaluation": {"candidates": candidates, "metrics": metrics},
                            }
                        )
            blobs[names[f"runs_{group}"]] = replay.encoded(runs)
            input_hashes[group] = replay.sha(blobs[names[f"runs_{group}"]])
            all_runs.extend(runs)
        common = {
            "dataset_sha256": replay.sha(blobs[names["panel"]]),
            "seeds": replay.SEEDS,
            "expected_runs": len(all_runs),
            "conformers": {
                "num_conformers": 8,
                "energy_window": 5.0,
                "max_iterations": 500,
                "pharmacophore_sigma": 1.0,
                "reject_unspecified_stereo": True,
                "strain_penalty": 0.05,
            },
            "code_sha256": {
                k: code[k] for k in ("benchmarks/research/controls.py", "chembounce_synthesis/conformers.py")
            },
        }
        blobs[names["common_protocol"]] = replay.encoded(common)
        scored = []
        for run in all_runs:
            for ordinal, candidate in enumerate(run["evaluation"]["candidates"], 1):
                success = ordinal == 1
                scores = {"bit_tanimoto": 0.2, "count_tanimoto": 0.25}
                if success:
                    for metric in helper.SCORES:
                        scores[metric] = 0.5
                        scores[metric + "_stability"] = {
                            "n": 3,
                            "mean": statistics.mean([0.4, 0.5, 0.7]),
                            "median": 0.5,
                            "min": 0.4,
                            "max": 0.7,
                            "range": 0.7 - 0.4,
                            "sample_sd": statistics.stdev([0.4, 0.5, 0.7]),
                        }
                trials = []
                for s, value in zip(replay.SEEDS, [0.4, 0.5, 0.7]):
                    trial = {
                        "seed": s,
                        "status": "success" if success or s != 43 else "failed",
                        "seconds": 0.1,
                    }
                    if trial["status"] == "success":
                        trial["metrics"] = dict.fromkeys(helper.SCORES, value)
                    else:
                        trial["reason"] = "synthetic failure"
                    trials.append(trial)
                q = int(run["query_id"][1:])
                cache = replay.sha(
                    json.dumps(
                        [
                            replay.sha(blobs[names["common_protocol"]]),
                            panel["compounds"][q]["smiles"],
                            candidate["smiles"],
                        ],
                        separators=(",", ":"),
                    ).encode()
                )
                scored.append(
                    {
                        "method_id": run["method_id"],
                        "query_id": run["query_id"],
                        "target_id": run["target_id"],
                        "generation_seed": run["seed"],
                        "rank": ordinal,
                        "smiles": candidate["smiles"],
                        "cache_key": cache,
                        "all_3d_seeds_succeeded": success,
                        "trials": trials,
                        "scores": scores,
                        "error": None,
                    }
                )
        blobs[names["common_generated_scores"]] = replay.encoded(scored)
        blobs[names["common_progress"]] = replay.encoded(
            {
                "protocol_sha256": replay.sha(blobs[names["common_protocol"]]),
                "input_runs_sha256": input_hashes,
                "runs_observed": len(all_runs),
                "expected_runs": len(all_runs),
                "candidate_occurrences": len(scored),
                "unique_query_candidate_pairs": len({s["cache_key"] for s in scored}),
                "assessed_occurrences": len(scored),
                "failed_3d_occurrences": sum(not s["all_3d_seeds_succeeded"] for s in scored),
                "queued_unique_pairs": 0,
                "complete": True,
            }
        )
        blobs[names["comparison"]] = b"{}"
        runs, scores, targets, caps, metadata = replay.validate_split(split, blobs, code)
        blobs[names["comparison"]] = replay.encoded(
            {
                **metadata,
                **helper.analyze_runs(runs, scores, targets, caps),
                **replay.report_text(blobs[replay.HELPER]),
            }
        )
    evidence = replay.encoded({"scope": "SYNTHETIC fixture only; not real study evidence", "attempts": 630})
    blobs["parent-review-evidence/audit.json"] = evidence
    blobs["parent-review.json"] = replay.encoded(
        {
            "schema_version": 1,
            "status": "passed",
            "total_runs": 630,
            "split_counts": {s: n for s, (n, _) in replay.SPLITS.items()},
            "historical_helper_sha256": replay.HELPER_SHA,
            "original_source_commit": replay.COMMIT,
            "source39_sha256": code,
            "input_sha256": {n: replay.sha(blobs[n]) for n in replay.input_names()},
            "assessment_evidence": [{"artifact": "audit.json", "sha256": replay.sha(evidence)}],
        }
    )
    return blobs


@pytest.fixture(scope="module")
def fixture_blobs():
    return make_fixture()


def write_source(root, blobs):
    root.mkdir()
    for name in replay.input_names() - {replay.ATTRIBUTION}:
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(blobs[name])
    (root / "parent-review.json").write_bytes(blobs["parent-review.json"])
    (root / "audit.json").write_bytes(blobs["parent-review-evidence/audit.json"])
    return root / "parent-review.json"


def refresh_review(blobs):
    review = replay.parse(blobs["parent-review.json"])
    review["input_sha256"] = {n: replay.sha(blobs[n]) for n in replay.input_names()}
    blobs["parent-review.json"] = replay.encoded(review)


def mutate(blobs, name, change):
    value = replay.parse(blobs[name])
    change(value)
    blobs[name] = replay.encoded(value)
    refresh_review(blobs)


def test_full_630_pack_relocation_replay_preserves_bytes_and_failures(tmp_path, fixture_blobs):
    source = tmp_path / "source"
    review = write_source(source, fixture_blobs)
    packed = replay.pack(source, review, replay.sha(review.read_bytes()), tmp_path / "bundle")
    relocated = tmp_path / "relocated"
    shutil.move(tmp_path / "bundle", relocated)
    shutil.rmtree(source)  # Neither recorded private paths nor original source paths can be dependencies.
    verified = replay.verify(relocated, packed["manifest_sha256"])
    assert verified["attempts"] == 630 and verified["comparison_recomputed"] is False
    result = replay.replay(relocated, packed["manifest_sha256"], tmp_path / "recomputed")
    assert result["status"] == "passed"
    for name, raw in fixture_blobs.items():
        assert (relocated / name).read_bytes() == raw
    comparison = replay.parse((tmp_path / "recomputed/scaffold-comparison.json").read_bytes())
    pipeline = comparison["methods"]["pipeline"]
    assert pipeline["attempts"] == 36
    assert pipeline["status_counts"]["timeout"] == 6
    assert pipeline["runtime_missing_attempts"] == 6
    assert pipeline["target_macro"]["common_3d_failed_candidates"]["mean"] > 0
    wall_interval = pipeline["target_macro"]["wall_seconds"]["ci95"]
    assert wall_interval[0] < wall_interval[1]
    paired_wall = comparison["paired_pipeline_minus_comparator"]["nearest_catalog"]["wall_seconds"]
    assert paired_wall["mean"] == pytest.approx(-0.1, abs=1e-12)
    assert "infrastructure_missing_outcome_sensitivity" in comparison
    with pytest.raises(ValueError, match="already exists"):
        replay.replay(relocated, packed["manifest_sha256"], tmp_path / "recomputed")
    with pytest.raises(ValueError, match="already exists"):
        replay.pack(relocated, relocated / "parent-review.json", packed["manifest_sha256"], relocated)
    with pytest.raises(ValueError, match="manifest SHA"):
        replay.verify(relocated, "0" * 64)
    path = relocated / replay.paths("target")["runs_gpu1"]
    path.write_bytes(path.read_bytes() + b" ")
    with pytest.raises(ValueError, match="artifact"):
        replay.verify(relocated, packed["manifest_sha256"])


@pytest.mark.parametrize(
    "change",
    [
        "missing-run",
        "duplicate-run",
        "seed",
        "target",
        "method",
        "cap",
        "incomplete",
        "ordinal",
        "score-extra",
        "progress-sha",
        "source39",
        "null-time",
        "failed-accepted",
        "comparison-field",
        "comparison-number",
        "comparison-type",
        "review576",
        "review-evidence",
        "helper",
    ],
)
def test_fail_closed_semantics_even_after_input_hashes_updated(fixture_blobs, change):
    blobs = dict(fixture_blobs)
    names = replay.paths("target")
    if change in {"missing-run", "duplicate-run", "seed", "target", "failed-accepted", "null-time"}:

        def modify(rows):
            if change == "missing-run":
                rows.pop()
            elif change == "duplicate-run":
                rows.append(rows[0])
            elif change == "seed":
                rows[0]["seed"] = 45
            elif change == "target":
                rows[0]["target_id"] = "wrong-target"
            elif change == "failed-accepted":
                rows[0]["status"] = "timeout"
            else:
                rows[0]["resources"]["wall_seconds"] = None

        mutate(blobs, names["runs_gpu1"], modify)
    elif change in {"method", "cap"}:
        mutate(
            blobs,
            names["protocol_gpu1"],
            lambda v: v["methods"][0].update(
                {"id": "invented"} if change == "method" else {"budget_overrides": {"max_candidates": 4}}
            ),
        )
    elif change in {"incomplete", "progress-sha"}:
        mutate(
            blobs,
            names["common_progress"],
            lambda v: v.update(
                {"complete": False} if change == "incomplete" else {"protocol_sha256": "0" * 64}
            ),
        )
    elif change in {"ordinal", "score-extra"}:
        mutate(
            blobs,
            names["common_generated_scores"],
            lambda v: v[1].update(rank=3) if change == "ordinal" else v.append(v[0]),
        )
    elif change == "source39":
        mutate(blobs, names["manifest_gpu1"], lambda v: v["code_sha256"].update({"cli.py": "0" * 64}))
    elif change.startswith("comparison"):

        def modify(v):
            if change == "comparison-field":
                v.pop("complete_case_secondary")
            elif change == "comparison-type":
                v["methods"]["pipeline"]["attempts"] = 18.0
            else:
                v["methods"]["pipeline"]["target_macro"]["wall_seconds"]["mean"] += 0.01

        mutate(blobs, names["comparison"], modify)
    elif change in {"review576", "review-evidence"}:
        review = replay.parse(blobs["parent-review.json"])
        review.update({"total_runs": 576} if change == "review576" else {"assessment_evidence": []})
        blobs["parent-review.json"] = replay.encoded(review)
    else:
        blobs[replay.HELPER] += b"\n"
    with pytest.raises((ValueError, KeyError)):
        replay.validate(blobs, calculate=True)


@pytest.mark.parametrize(
    "raw", [b'{"x":1,"x":2}', b'{"x":NaN}', b'{"x":Infinity}', b'{"x":-Infinity}', b'{"x":1e309}']
)
def test_json_duplicates_and_nonfinite_rejected(raw):
    with pytest.raises(ValueError):
        replay.parse(raw)


@pytest.mark.parametrize("name", ["../escape", "/absolute", "a/../../b", "a//b", "./x", "a\\b", "C:/x"])
def test_artifact_path_escape_rejected(name):
    with pytest.raises(ValueError):
        replay.safe_relative(name)


def test_symlink_and_parent_evidence_escape_rejected(tmp_path, fixture_blobs):
    (tmp_path / "original").write_text("{}")
    (tmp_path / "linked").symlink_to(tmp_path / "original")
    with pytest.raises(ValueError, match="Symlink"):
        replay.read_bytes(tmp_path / "linked")
    blobs = dict(fixture_blobs)
    review = replay.parse(blobs["parent-review.json"])
    review["assessment_evidence"][0]["artifact"] = "../outside.json"
    blobs["parent-review.json"] = replay.encoded(review)
    source = tmp_path / "source"
    rp = write_source(source, blobs)
    with pytest.raises(ValueError, match="Unsafe"):
        replay.pack(source, rp, replay.sha(rp.read_bytes()), tmp_path / "bundle")
    assert not (tmp_path / "bundle").exists()


def test_mutation_during_operation_fails_before_manifest(tmp_path, fixture_blobs, monkeypatch):
    source = tmp_path / "source"
    rp = write_source(source, fixture_blobs)
    real_validate = replay.validate

    def during(blobs, calculate=False):
        real_validate(blobs, calculate=calculate)
        p = source / replay.paths("target")["runs_gpu0"]
        p.write_bytes(p.read_bytes() + b" ")

    monkeypatch.setattr(replay, "validate", during)
    with pytest.raises(ValueError, match="input changed"):
        replay.pack(source, rp, replay.sha(rp.read_bytes()), tmp_path / "bundle")
    assert not (tmp_path / "bundle").exists()


def test_tolerance_types_and_fields():
    replay.same({"a": 0.0}, {"a": 1e-12}, tolerance=1e-12)
    for actual in (
        {"a": 1.01e-12},
        {"a": 0},
        {"a": False},
        {"b": 0.0},
        {"a": float("nan")},
        {"a": float("inf")},
    ):
        with pytest.raises(ValueError):
            replay.same(actual, {"a": 0.0}, tolerance=1e-12)


def test_explicit_runtime_fields_cannot_be_silently_projected_out():
    raw = {"methods": [{"id": "pipeline", "resolved_assets": {"x": "explicit"}}]}
    changed = copy.deepcopy(raw)
    changed["methods"][0]["resolved_assets"] = {"x": "changed"}
    with pytest.raises(ValueError):
        replay.normalized_protocol(raw, changed)


def write_bundle(path, blobs):
    path.mkdir()
    for name, raw in blobs.items():
        file = path / name
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_bytes(raw)
    manifest = {
        "schema_version": 1,
        "kind": "original-630-statistical-replay",
        "scope": replay.SCOPE,
        "historical_helper_sha256": replay.HELPER_SHA,
        "original_source_commit": replay.COMMIT,
        "files": {n: {"sha256": replay.sha(v), "bytes": len(v)} for n, v in blobs.items()},
    }
    (path / "manifest.json").write_bytes(replay.encoded(manifest))
    return manifest


@pytest.mark.parametrize("failure", ["escape", "symlink", "extra", "duplicate-json", "missing-comparison"])
def test_relocated_bundle_boundary_checks(tmp_path, fixture_blobs, failure):
    bundle = tmp_path / "bundle"
    manifest = write_bundle(bundle, fixture_blobs)
    file = replay.paths("target")["comparison"]
    if failure == "escape":
        manifest["files"]["../outside.json"] = manifest["files"].pop(file)
    elif failure == "symlink":
        outside = tmp_path / "outside.json"
        (bundle / file).rename(outside)
        (bundle / file).symlink_to(outside)
    elif failure == "extra":
        (bundle / "extra.json").write_text("{}")
    elif failure == "missing-comparison":
        (bundle / file).unlink()
        manifest["files"].pop(file)
    raw = replay.encoded(manifest)
    if failure == "duplicate-json":
        raw = raw.replace(b'"schema_version": 1,', b'"schema_version": 1, "schema_version": 1,')
    (bundle / "manifest.json").write_bytes(raw)
    with pytest.raises((ValueError, KeyError)):
        replay.verify(bundle, replay.sha(raw))


def test_incomplete_target_blocks_all_computation_and_output(tmp_path, fixture_blobs, monkeypatch):
    blobs = dict(fixture_blobs)
    mutate(blobs, replay.paths("target")["common_progress"], lambda v: v.update(complete=False))
    bundle = tmp_path / "bundle"
    raw = replay.encoded(write_bundle(bundle, blobs))

    def forbidden(_):
        pytest.fail("Historical statistics must not start before ALL three splits validate")

    monkeypatch.setattr(replay, "historical", forbidden)
    with pytest.raises(ValueError, match="common completion"):
        replay.replay(bundle, replay.sha(raw), tmp_path / "output")
    assert not (tmp_path / "output").exists()
    with pytest.raises(ValueError, match="inside the bundle"):
        replay.replay(bundle, replay.sha(raw), bundle / "output")


def test_external_review_hash_and_existing_output_preserved(tmp_path, fixture_blobs):
    review = write_source(tmp_path / "source", fixture_blobs)
    with pytest.raises(ValueError, match="external parent review SHA"):
        replay.pack(tmp_path / "source", review, "0" * 64, tmp_path / "out")
    existing = tmp_path / "out"
    existing.mkdir()
    sentinel = existing / "unchanged"
    sentinel.write_bytes(b"preserve me")
    with pytest.raises(ValueError, match="already exists"):
        replay.pack(tmp_path / "source", review, replay.sha(review.read_bytes()), existing)
    assert sentinel.read_bytes() == b"preserve me"


@pytest.mark.parametrize("field", ["median", "range", "n", "failed-with-summary"])
def test_stored_trial_arithmetic_and_failure_summary(fixture_blobs, field):
    blobs = dict(fixture_blobs)
    name = replay.paths("target")["common_generated_scores"]

    def change(rows):
        if field == "failed-with-summary":
            rows[1]["scores"]["hit_shape"] = 0.5
        elif field == "median":
            rows[0]["scores"]["hit_shape"] = 0.6
        else:
            rows[0]["scores"]["hit_shape_stability"][field] = 1 if field == "n" else 0.4

    mutate(blobs, name, change)
    with pytest.raises(ValueError, match="3D trial|Failed 3D"):
        replay.validate(blobs)


@pytest.mark.parametrize("change", ["digest-instead-of-path", "changed-path", "missing-source-id"])
def test_source_path_provenance_matches_actual_manifest_schema(fixture_blobs, change):
    blobs = dict(fixture_blobs)

    def modify(manifest):
        paths = manifest["dataset"]["source_artifacts"]
        if change == "missing-source-id":
            paths.pop("status")
        else:
            paths["status"] = (
                "0" * 64 if change == "digest-instead-of-path" else "/different-origin/status.json"
            )

    mutate(blobs, replay.paths("target")["manifest_gpu0"], modify)
    with pytest.raises(ValueError, match="source path"):
        replay.validate(blobs)


@pytest.mark.parametrize("field", ["sha256", "url"])
def test_panel_source_hash_and_url_remain_separately_bound(fixture_blobs, field):
    blobs = dict(fixture_blobs)
    mutate(
        blobs,
        replay.paths("target")["panel"],
        lambda p: p["sources"][0].update(
            {field: "0" * 64 if field == "sha256" else "https://invalid.example/source"}
        ),
    )
    with pytest.raises(ValueError, match="ChEMBL source"):
        replay.validate(blobs)


def test_recorded_path_construction_is_lexical_and_never_dereferenced(monkeypatch):
    def no_filesystem(*args, **kwargs):
        pytest.fail("Recorded provenance paths must not be read or resolved")

    monkeypatch.setattr(Path, "open", no_filesystem)
    monkeypatch.setattr(Path, "resolve", no_filesystem)
    replay.recorded_source_paths(
        [{"id": "a", "artifact": "/synthetic/a/../b.json"}], {"a": "/synthetic/b.json"}
    )
    with pytest.raises(ValueError, match="relative provenance"):
        replay.recorded_source_paths([{"id": "a", "artifact": "../../b.json"}], {"a": "/synthetic/b.json"})
