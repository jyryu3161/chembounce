"""Synthetic stored records only; no study outcomes or chemistry are generated."""

import builtins
import copy
import importlib.util
from pathlib import Path
import shutil

import pytest

HERE = Path(__file__).resolve().parents[1] / "benchmarks/diagnostics/supplement_replay"
spec = importlib.util.spec_from_file_location("supplement_replay_fixture", HERE / "replay.py")
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)


def put(blobs, name, value):
    blobs[name] = r.io.encoded(value)


def refresh(blobs):
    for prefix in ("common", "analysis"):
        m = r.io.parse(blobs[prefix + "/manifest.json"])
        m["files_sha256"] = {
            n[len(prefix) + 1 :]: r.io.sha(v)
            for n, v in blobs.items()
            if n.startswith(prefix + "/") and n != prefix + "/manifest.json"
        }
        if prefix == "analysis":
            m["common_scores_manifest_sha256"] = r.io.sha(blobs["common/manifest.json"])
        put(blobs, prefix + "/manifest.json", m)
    inputs = {n: r.io.sha(v) for n, v in blobs.items() if n.startswith(("analysis/", "common/"))}
    integrity = r.io.parse(blobs["integrity-validation.json"])
    integrity["source_sha256"] = inputs
    integrity["output_sha256"] = {
        "normalized-occurrences.json": r.io.sha(blobs["normalized-occurrences.json"])
    }
    put(blobs, "integrity-validation.json", integrity)
    review = r.io.parse(blobs["parent-review.json"])
    review["input_sha256"] = inputs | {
        n: r.io.sha(blobs[n]) for n in ("integrity-validation.json", "normalized-occurrences.json")
    }
    put(blobs, "parent-review.json", review)


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    root = tmp_path_factory.mktemp("supplement-synthetic")
    fixed = r.fixed_files()
    blobs = {n: p.read_bytes() for n, (p, _) in fixed.items()}
    grid = []
    for split, count in r.SPLITS.items():
        for method in r.METHODS:
            for q in range(count // 9):
                for seed in r.SEEDS:
                    grid.append(
                        {
                            "split": split,
                            "method_id": method,
                            "query_id": f"synthetic-q{q}",
                            "target_id": f"synthetic-T{q // 2}",
                            "seed": seed,
                            "query_index": q,
                            "output_dir": f"/historical/nonexistent/{split}/{method}/{q}/{seed}",
                        }
                    )
    plan = r.io.parse(blobs["prepared/analysis-plan.json"])
    plan["expected_grid"] = grid
    put(blobs, "prepared/grid.json", grid)
    put(blobs, "prepared/analysis-plan.json", plan)
    common = r.io.parse(blobs["prepared/common-3d-plan.json"])
    common["source"]["root"] = "/historical/nonexistent/scientific-source"
    put(blobs, "prepared/common-3d-plan.json", common)
    execution = {
        "source": common["source"],
        "commands": [{"argv": ["/historical/nonexistent/python"]}],
        "runtimes": {"/historical/nonexistent/python": {"packages": {"numpy": "1.26.4"}}},
    }
    put(blobs, "prepared/execution-plan.json", execution)
    pm = r.io.parse(blobs["prepared/manifest.json"])
    for name in ("grid.json", "analysis-plan.json", "common-3d-plan.json", "execution-plan.json"):
        pm["files_sha256"][name] = r.io.sha(blobs["prepared/" + name])
    put(blobs, "prepared/manifest.json", pm)
    prepared_sha = r.io.sha(blobs["prepared/manifest.json"])
    put(blobs, "generation/generation.json", {"expected_attempts": 270, "campaigns": [], "synthetic": True})
    gm = r.io.parse(blobs["generation/manifest.json"])
    gm["study_manifest_sha256"] = prepared_sha
    gm["files_sha256"] = {"generation.json": r.io.sha(blobs["generation/generation.json"])}
    put(blobs, "generation/manifest.json", gm)
    generation_sha = r.io.sha(blobs["generation/manifest.json"])
    gate = r.io.encoded({"passed": True, "synthetic": True, "total_attempts": 270})
    gate_sha = r.io.sha(gate)
    # Only frozen metadata/parent trust anchors are replaced in fixtures. Original
    # scientific helpers and statistical extraction bytes are unchanged.
    for name, (path, pin) in list(fixed.items()):
        if name.startswith(("prepared/", "generation/")):
            p = root / "fixed" / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(blobs[name])
            fixed[name] = (p, r.io.sha(blobs[name]))
    stats = r.module_bytes(blobs["statistics_extracted.py"], "synthetic_statistics")
    stats.preparation = r.SimpleNamespace(SPLITS=tuple(r.SPLITS), METHODS=r.METHODS)
    attempts = []
    occurrences = []
    joined = []
    scores = {}
    retained = []
    evaluator_metrics = plan["evaluator_metrics"]
    for g in grid:
        q = g["query_index"]
        seed = g["seed"]
        method = g["method_id"]
        status = "no_candidates"
        if q == 0 and method == "pipeline_current":
            status = {42: "success", 43: "worker_error", 44: "infrastructure_interruption"}[seed]
        if q == 0 and method == "crem_train_partial_exo" and seed in (42, 43):
            status = "success"
        candidates = []
        selected = (
            [("A", 1), ("B", 3)]
            if status == "success" and method == "pipeline_current"
            else ([("C" + str(seed), 1)] if status == "success" else [])
        )
        for label, rank in selected:
            candidates.append(
                {
                    "rank": rank,
                    "smiles": "synthetic-" + label,
                    "family": "synthetic-family-" + label,
                    "hit_bit_tanimoto": 0.7,
                    "hit_count_tanimoto": 0.6,
                    "route_evidence": {"level": "unassessed"},
                    "declared_external_synthesis_evidence": None,
                }
            )
        metrics = {n: 0 for n in evaluator_metrics}
        metrics.update(
            run_success=int(status in {"success", "no_candidates"}),
            candidate_count=3 if method == "pipeline_current" and candidates else len(candidates),
            valid_candidate_count=len(candidates),
            valid_fraction=len(candidates) / 3 if candidates else 0.0,
            unique_families=len(candidates),
            hop_families=len(candidates),
            reference_active_scaffolds=0,
            active_scaffold_recall_proxy=None,
        )
        a = {k: g[k] for k in r.FIELDS} | {
            "status": status,
            "evaluation": {"candidates": candidates, "metrics": metrics},
            "resources": {"wall_seconds": 1.0},
            "raw_accounting": {"raw_count_scope": "synthetic"},
            "response_sha256": "a" * 64,
        }
        attempts.append(a)
        for c in candidates:
            retention = {
                n: None if c["smiles"].endswith("B") else n != "query_scaffold_not_contained"
                for n in stats.RETENTION_FIELDS
            }
            retention |= {"query_scaffold_heavy_atoms": 6, "candidate_scaffold_heavy_atoms": 7}
            o = {k: g[k] for k in r.FIELDS} | {
                "rank": c["rank"],
                "smiles": c["smiles"],
                "query_exact": "synthetic-query",
                "candidate": c,
                "route_supplied": False,
                "serialized_route_sha256": None,
                "scaffold_retention": retention,
            }
            occurrences.append(o)
            payload = {
                "query_exact": o["query_exact"],
                "candidate_exact": o["smiles"],
                "source": common["source"],
                "conformers": common["conformers"],
                "seeds": r.SEEDS,
                "identity_policy": "protonation-tautomer-family-v1",
                "rdkit": "2025.09.6",
            }
            key = r.content_key(payload)
            trials = []
            middle = 0.4 if c["smiles"].endswith("A") else 0.2 if c["smiles"].endswith("42") else 0.9
            for s in r.SEEDS:
                if c["smiles"].endswith("B") and s == 43:
                    trials.append(
                        {"seed": s, "status": "failed", "reason": "synthetic 3D failure", "seconds": 0.1}
                    )
                else:
                    trials.append(
                        {
                            "seed": s,
                            "status": "success",
                            "metrics": {"hit_shape": middle, "hit_pharmacophore": 0.5, "score": 0.3},
                            "seconds": 0.1,
                        }
                    )
            success = all(t["status"] == "success" for t in trials)
            values = {"bit_tanimoto": 0.7, "count_tanimoto": 0.6}
            if success:
                for name, value in [("hit_shape", middle), ("hit_pharmacophore", 0.5), ("score", 0.3)]:
                    values[name] = value
                    values[name + "_stability"] = stats.score_stability([value] * 3)
            score = {
                "payload": payload,
                "cache_key": key,
                "scores": values,
                "trials": trials,
                "all_3d_seeds_succeeded": success,
                "pair_error": None,
                "wall_seconds": 0.3,
            }
            put(blobs, "common/cache/" + key + ".json", score)
            scores[stats.occurrence_key(o)] = score
            joined.append(
                {k: o[k] for k in (*r.FIELDS, "rank", "smiles", "query_exact")} | {"cache_key": key}
            )
            retained.append(
                {k: o[k] for k in (*r.FIELDS, "rank", "query_exact", "smiles", "scaffold_retention")}
                | {"family": c["family"]}
            )
    rows = stats.attempt_metrics(attempts, occurrences, scores)
    put(blobs, "analysis/attempts.json", rows)
    put(blobs, "analysis/recomputed-attempts.json", attempts)
    blobs["common/recomputed-attempts.json"] = blobs["analysis/recomputed-attempts.json"]
    put(blobs, "analysis/comparison.json", stats.comparisons(rows, grid, stats.cluster_interval))
    for name in ("raw-accounting-summary.json", "unique-identities.json", "scaffold-retention-summary.json"):
        put(blobs, "analysis/" + name, {"parent_preserved_synthetic": True})
    blobs["analysis/analysis-plan.json"] = blobs["prepared/analysis-plan.json"]
    put(blobs, "analysis/scaffold-retention-occurrences.json", retained)
    put(blobs, "normalized-occurrences.json", occurrences)
    put(blobs, "common/occurrences.json", joined)
    put(blobs, "common/context.json", {"complete": False, "completion_requires_final_manifest": True})
    put(
        blobs,
        "common/coverage.json",
        {
            "attempts": 270,
            "valid_candidate_occurrences": len(occurrences),
            "unique_exact_pairs": len({j["cache_key"] for j in joined}),
            "successful_3d_occurrences": sum(s["all_3d_seeds_succeeded"] for s in scores.values()),
        },
    )
    cm = {
        "schema_version": 1,
        "stage": "common_3d",
        "complete": True,
        "assessment_code_sha256": r.HELPER_SHA,
        "study_manifest_sha256": prepared_sha,
        "generation_manifest_sha256": generation_sha,
        "common_plan_sha256": r.io.sha(blobs["prepared/common-3d-plan.json"]),
        "source": common["source"],
        "runtime_python": "/historical/nonexistent/python",
        "runtime_packages": {"numpy": "1.26.4"},
        "execution": "CPU only, sequential, generation budgets excluded",
    }
    put(blobs, "common/manifest.json", cm)
    put(
        blobs,
        "analysis/manifest.json",
        {
            "schema_version": 1,
            "stage": "final_comparison",
            "complete_requires_parent_review": True,
            "expected_attempts": 270,
            "assessment_code_sha256": r.HELPER_SHA,
            "study_manifest_sha256": prepared_sha,
            "generation_manifest_sha256": generation_sha,
        },
    )
    put(
        blobs,
        "integrity-validation.json",
        {
            "passed": True,
            "total_attempts": 270,
            "candidate_occurrences": len(occurrences),
            "raw_candidate_joins_verified": True,
            "route_flags_verified": True,
            "scaffold_annotations_verified": True,
        },
    )
    blobs["parent-review-evidence/generation.json"] = gate
    put(
        blobs,
        "parent-review-evidence/final.json",
        {"passed": True, "synthetic": True, "no_real_study_claim": True},
    )
    put(
        blobs,
        "parent-review.json",
        {
            "schema_version": 1,
            "status": "passed",
            "complete": True,
            "stage": "final_supplementary_270_completion",
            "total_attempts": 270,
            "split_attempts": r.SPLITS,
            "analysis_helper_sha256": r.HELPER_SHA,
            "prepared_manifest_sha256": prepared_sha,
            "generation_manifest_sha256": generation_sha,
            "generation_parent_gate_sha256": gate_sha,
            "source_commit": r.SOURCE_COMMIT,
            "assessment_evidence": [
                {"role": "generation_parent_gate", "artifact": "generation.json", "sha256": gate_sha},
                {
                    "role": "final_integrity_and_statistics",
                    "artifact": "final.json",
                    "sha256": r.io.sha(blobs["parent-review-evidence/final.json"]),
                },
            ],
        },
    )
    refresh(blobs)
    return root, blobs, fixed, prepared_sha, generation_sha, gate_sha


@pytest.fixture
def fixture(synthetic, monkeypatch):
    root, blobs, fixed, prepared, generation, gate = synthetic
    monkeypatch.setattr(r, "fixed_files", lambda: fixed)
    monkeypatch.setattr(r, "PREPARED_SHA", prepared)
    monkeypatch.setattr(r, "GENERATION_SHA", generation)
    monkeypatch.setattr(r, "GENERATION_GATE_SHA", gate)
    return dict(blobs)


def test_full_comparison_hand_contrasts_missingness_and_raw_rank_gaps(fixture):
    comparison = r.validate(fixture)
    scaffold = comparison["scaffold"]
    pipeline = scaffold["methods"]["pipeline_current"]
    assert pipeline["attempts"] == 36
    assert pipeline["metrics"]["attempt_with_hop"]["target_macro"]["mean"] == pytest.approx(1 / 36)
    assert pipeline["metrics"]["common_3d_output_slot_coverage"]["target_macro"]["mean"] == pytest.approx(
        1 / 108
    )
    assert pipeline["metrics"]["common_hit_shape_mean_success_conditional"]["target_macro"]["mean"] == 0.4
    assert (
        pipeline["metrics"]["common_hit_shape_mean_success_conditional"]["availability"]["observed_slots"]
        == 1
    )
    assert pipeline["metrics"]["common_hit_shape_mean_success_conditional"]["target_macro"]["ci95"] is None
    assert (
        pipeline["metrics"]["common_hit_shape_mean_success_conditional"]["target_macro"][
            "independent_targets"
        ]
        == 1
    )
    assert (
        scaffold["methods"]["mol2mol_pretrained"]["metrics"]["common_hit_shape_mean_success_conditional"][
            "target_macro"
        ]["mean"]
        is None
    )
    paired = scaffold["paired_pipeline_current_minus_comparator"]["crem_train_partial_exo"]
    assert paired["metrics"]["common_hit_shape_mean_success_conditional"]["target_macro"][
        "mean"
    ] == pytest.approx(0.2)
    assert pipeline["infrastructure_counterfactual_bounds"]["attempt_with_hop"] == pytest.approx(
        [1 / 36, 2 / 36]
    )
    occurrences = r.io.parse(fixture["normalized-occurrences.json"])
    assert [o["rank"] for o in occurrences[:2]] == [1, 3]
    assert occurrences[1]["scaffold_retention"]["query_scaffold_subgraph_retained"] is None


def test_portable_roundtrip_no_chemistry_and_preserved_original_bytes(fixture, tmp_path, monkeypatch):
    root = tmp_path / "input"
    root.mkdir()
    for name, raw in fixture.items():
        if name.startswith(("analysis/", "common/")):
            p = root / name
        elif name.startswith("parent-review-evidence/"):
            p = root / "reviews" / name.split("/", 1)[1]
        elif name == "parent-review.json":
            p = root / "reviews" / "completion.json"
        elif name in ("normalized-occurrences.json", "integrity-validation.json"):
            p = root / "normalized" / name
        else:
            continue
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(raw)
    importer = builtins.__import__

    def guarded(name, *args, **kwargs):
        assert not name.startswith(("rdkit", "torch", "chembounce_synthesis", "benchmarks.research"))
        return importer(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded)
    bundle = tmp_path / "bundle"
    packed = r.pack(
        root / "normalized",
        root / "analysis",
        root / "common",
        root / "reviews/completion.json",
        r.io.sha(fixture["parent-review.json"]),
        bundle,
    )
    relocated = tmp_path / "relocated"
    bundle.rename(relocated)
    shutil.rmtree(root)
    report = r.check(relocated, packed["manifest_sha256"], tmp_path / "replay")
    assert report["comparison_recomputed"] is True
    r.io.same(
        r.io.parse((tmp_path / "replay/comparison.json").read_bytes()),
        r.io.parse(fixture["analysis/comparison.json"]),
    )
    assert (relocated / "analysis/raw-accounting-summary.json").read_bytes() == fixture[
        "analysis/raw-accounting-summary.json"
    ]
    assert (relocated / "assess_historical.py").read_bytes() == fixture["assess_historical.py"]
    with pytest.raises(ValueError, match="Output exists"):
        r.check(relocated, packed["manifest_sha256"], tmp_path / "replay")
    (relocated / "unexpected.json").write_text("{}")
    with pytest.raises(ValueError, match="bundle membership"):
        r.check(relocated, packed["manifest_sha256"], tmp_path / "check.json", False)


@pytest.mark.parametrize(
    "change",
    [
        "missing",
        "duplicate",
        "order",
        "target",
        "rank",
        "candidate",
        "cache_seed",
        "cache_config",
        "cache_median",
        "conditional",
        "comparison",
        "incomplete",
        "review",
        "raw_preserved",
    ],
)
def test_adversarial_even_with_refreshed_container_hashes(fixture, change):
    if change in {"missing", "duplicate", "order", "target"}:
        x = r.io.parse(fixture["analysis/recomputed-attempts.json"])
        if change == "missing":
            x.pop()
        elif change == "duplicate":
            x[-1] = copy.deepcopy(x[0])
        elif change == "order":
            x[0], x[1] = x[1], x[0]
        else:
            x[0]["target_id"] = "wrong-target"
        put(fixture, "analysis/recomputed-attempts.json", x)
        fixture["common/recomputed-attempts.json"] = fixture["analysis/recomputed-attempts.json"]
    elif change in {"rank", "candidate"}:
        x = r.io.parse(fixture["normalized-occurrences.json"])
        if change == "rank":
            x[1]["rank"] = 2
        else:
            x[0]["candidate"]["family"] = "different-family"
        put(fixture, "normalized-occurrences.json", x)
    elif change.startswith("cache_"):
        name = next(
            n
            for n in fixture
            if n.startswith("common/cache/") and r.io.parse(fixture[n])["all_3d_seeds_succeeded"]
        )
        x = r.io.parse(fixture[name])
        if change == "cache_seed":
            x["trials"].pop()
        elif change == "cache_config":
            x["payload"]["conformers"]["num_conformers"] = 9
        else:
            x["scores"]["hit_shape"] += 0.1
        put(fixture, name, x)
    elif change == "conditional":
        x = r.io.parse(fixture["analysis/attempts.json"])
        x[1]["metrics"]["common_hit_shape_mean_success_conditional"] = 0.0
        put(fixture, "analysis/attempts.json", x)
    elif change == "comparison":
        x = r.io.parse(fixture["analysis/comparison.json"])
        x["scaffold"]["methods"]["pipeline_current"]["attempts"] = 35
        put(fixture, "analysis/comparison.json", x)
    elif change == "incomplete":
        x = r.io.parse(fixture["common/manifest.json"])
        x["complete"] = False
        put(fixture, "common/manifest.json", x)
    elif change == "review":
        x = r.io.parse(fixture["parent-review.json"])
        x["complete"] = False
        put(fixture, "parent-review.json", x)
    else:
        fixture["analysis/raw-accounting-summary.json"] += b" "
        # Do not refresh hashes: preserved chemistry is protected by parent hashes.

    if change != "raw_preserved":
        refresh(fixture)
    with pytest.raises((ValueError, KeyError)):
        r.validate(fixture)


def test_extract_only_allowlisted_pure_nodes_and_exact_frozen_bytes(fixture):
    raw, records = r.extraction(fixture)
    assert raw == fixture["statistics_extracted.py"]
    assert all(
        n["node"]
        not in {"recompute", "scientific_api", "score_pair", "score_generation", "analyze_generation", "main"}
        for n in records
    )
    fixture["assess_historical.py"] += b'\nraise RuntimeError("must never execute")\n'
    assert r.extraction(fixture)[0] == raw
    with pytest.raises(ValueError, match="fixed artifact"):
        r.validate(fixture)


def test_parent_completion_gate_precedes_live_archive_reads(tmp_path):
    parent = tmp_path / "review" / "parent.json"
    parent.parent.mkdir()
    parent.write_bytes(r.io.encoded({"status": "pending", "complete": False}))
    with pytest.raises(ValueError, match="Final parent completion"):
        r.pack(
            tmp_path / "absent-inputs",
            tmp_path / "absent-analysis",
            tmp_path / "live-scores-not-read",
            parent,
            r.io.sha(parent.read_bytes()),
            tmp_path / "output",
        )
    assert not (tmp_path / "output").exists()
    with pytest.raises(ValueError, match="external final parent review SHA"):
        r.pack(
            tmp_path / "absent-inputs",
            tmp_path / "absent-analysis",
            tmp_path / "live-scores-not-read",
            parent,
            "0" * 64,
            tmp_path / "output",
        )


def test_exact_hierarchy_and_absolute_tolerance(fixture):
    stats, _, _, _ = r.check_fixed(fixture)
    slots = [
        {"query_id": q, "target_id": target, "seed": seed}
        for q, target in [("q1", "A"), ("q2", "A"), ("q3", "B")]
        for seed in [42, 43, 44]
    ]
    values = {
        ("q1", 42): 1.0,
        ("q1", 43): 0.0,
        ("q1", 44): None,
        ("q2", 42): 0.5,
        ("q2", 43): 0.5,
        ("q2", 44): 0.5,
        ("q3", 42): 0.0,
        ("q3", 43): 0.0,
        ("q3", 44): 0.0,
    }
    result = stats.hierarchical(values, slots, stats.cluster_interval)
    assert result["per_target"] == {"A": 0.5, "B": 0.0}
    assert result["target_macro"]["mean"] == 0.25
    assert result["target_macro"]["ci95"] == [0.0, 0.5]
    assert result["availability"]["missing_slots"] == 1
    r.io.same(0.1 + 5e-13, 0.1, tolerance=1e-12)
    with pytest.raises(ValueError):
        r.io.same(1e12 + 0.1, 1e12, tolerance=1e-12)


@pytest.mark.parametrize("raw", [b'{"x":1,"x":2}', b'{"x":NaN}', b'{"x":1e999}'])
def test_strict_json(raw):
    with pytest.raises(ValueError):
        r.io.parse(raw)


@pytest.mark.parametrize("name", ["../escape", "/absolute", "a//b", "a\\b", "a/./b"])
def test_unsafe_paths(name):
    with pytest.raises(ValueError):
        r.io.safe_relative(name)


def test_symlinks_and_external_anchor_fail_before_replay(fixture, tmp_path):
    root = tmp_path / "bundle"
    root.mkdir()
    for name, raw in fixture.items():
        p = root / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(raw)
    manifest = r.io.encoded(r.manifest_for(fixture))
    (root / "manifest.json").write_bytes(manifest)
    with pytest.raises(ValueError, match="external bundle manifest"):
        r.load_bundle(root, "0" * 64)
    original = root / "normalized-occurrences.json"
    saved = tmp_path / "saved"
    original.rename(saved)
    original.symlink_to(saved)
    with pytest.raises(ValueError, match="Symlink"):
        r.load_bundle(root, r.io.sha(manifest))
