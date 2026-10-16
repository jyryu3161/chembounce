"""Hand-derived contrasts and adversarial complete-study fixtures; no models/DBs."""

import copy
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest
from rdkit import Chem

from benchmarks.research.controls import similarity_scores
from benchmarks.research.evaluation import cluster_interval, evaluate_candidates, score_stability
from benchmarks.research.protocol import load_dataset
from benchmarks.research.runner import evidence_verifier, public_request
from chembounce_synthesis.identity import chemical_identity

MODULE = Path(__file__).resolve().parents[1] / "benchmarks/research/assess_comparator_followup.py"
spec = importlib.util.spec_from_file_location("followup_assessment_fixture", MODULE)
assess = importlib.util.module_from_spec(spec)
spec.loader.exec_module(assess)


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(assess.encode(value))


@pytest.fixture
def api():
    calls = []

    def synthetic_3d(hit, product, seeds, conformers):
        # Actual RDKit fingerprints; deliberately hand-set conformer outcomes.
        calls.append((hit, product))
        scores, _ = similarity_scores(hit, product, seeds)
        trials = []
        for index, seed in enumerate(seeds):
            if "s" in product and index == 1:
                trials.append({"seed": seed, "status": "failed", "reason": "synthetic conformer failure"})
            else:
                trials.append(
                    {
                        "seed": seed,
                        "status": "success",
                        "metrics": {
                            "hit_shape": [0.2, 0.4, 0.6][index],
                            "hit_pharmacophore": [0.8, 0.7, 0.9][index],
                            "score": [0.3, 0.5, 0.7][index],
                        },
                    }
                )
        if all(t["status"] == "success" for t in trials):
            for metric, median in (("hit_shape", 0.4), ("hit_pharmacophore", 0.8), ("score", 0.5)):
                scores[metric] = median
                scores[metric + "_stability"] = score_stability([t["metrics"][metric] for t in trials])
        return scores, trials

    return SimpleNamespace(
        load_dataset=load_dataset,
        public_request=public_request,
        evidence_verifier=evidence_verifier,
        evaluate=evaluate_candidates,
        identity=chemical_identity,
        similarity=synthetic_3d,
        stability=score_stability,
        interval=cluster_interval,
        calls=calls,
        mol_from_smiles=Chem.MolFromSmiles,
    )


@pytest.fixture
def completed(tmp_path, monkeypatch, api):
    root = tmp_path / "prepared"
    root.mkdir()
    source = tmp_path / "source"
    source.mkdir()
    public_source = tmp_path / "source.txt"
    public_source.write_text("Synthetic independent-assessment fixture; no experimental claim.\n")
    decision = tmp_path / "decision.json"
    amendment = tmp_path / "amendment.json"
    save(decision, {"source_root": str(tmp_path / "original")})
    save(
        amendment,
        {
            "execution_policy": {"protected_resources": [str(tmp_path / "protected")]},
            "current_pipeline": {"source_worktree": str(source)},
        },
    )
    packages = assess.runtime_packages()
    plan = {
        "parent_decision": {"artifact": str(decision)},
        "amendment": {"artifact": str(amendment)},
        "source": {"root": str(source), "commit": "fixture", "files_sha256": {"fixture.py": "a" * 64}},
        "runtimes": {sys.executable: {"version": sys.version, "packages": packages}},
        "commands": [],
    }
    grid = []
    for split in assess.preparation.SPLITS:
        nq = assess.preparation.COUNTS[split]
        compounds = []
        queries = []
        for i in range(nq):
            compounds.append(
                {
                    "id": f"input-{i}",
                    "target_id": f"T{i // 2}",
                    "smiles": "c1ccccc1",
                    "label": "active",
                    "split": "train",
                    "source_id": "fixture",
                    "assay_id": "fixture",
                    "label_basis": "synthetic",
                }
            )
            queries.append({"id": f"query-{i}", "compound_id": f"input-{i}"})
        for i in range(nq // 2):
            for label, smiles in (("active", "c1ccncc1"), ("inactive", "c1ccsc1")):
                compounds.append(
                    {
                        "id": f"{label}-{i}",
                        "target_id": f"T{i}",
                        "smiles": smiles,
                        "label": label,
                        "split": "test",
                        "source_id": "fixture",
                        "assay_id": "fixture",
                        "label_basis": "synthetic",
                    }
                )
        panel = tmp_path / f"{split}-panel.json"
        save(
            panel,
            {
                "schema_version": 1,
                "dataset_id": split + "-fixture",
                "scope": "synthetic_fixture",
                "sources": [
                    {
                        "id": "fixture",
                        "url": "https://example.invalid/fixture",
                        "artifact": str(public_source),
                        "sha256": assess.digest(public_source),
                    }
                ],
                "compounds": compounds,
                "queries": queries,
                "split": {"kind": "scaffold"},
            },
        )
        dataset, _ = load_dataset(panel)
        save(
            root / f"public-queries/{split}.json",
            [{"id": q["id"], "smiles": "c1ccccc1", "target_id": f"T{i // 2}"} for i, q in enumerate(queries)],
        )
        for method in assess.preparation.METHODS:
            output = tmp_path / "runs" / split / method
            method_object = {
                "id": method,
                "kind": "external",
                "assets": {},
                "settings": {},
                "budget_overrides": {"max_evaluations": 12 if method == "pipeline_current" else 256},
            }
            protocol = {
                "schema_version": 1,
                "methods": [method_object],
                "budget": assess.preparation.BUDGET,
                "seeds": [42, 43, 44],
            }
            protocol_path = root / f"protocols/{split}-{method}.json"
            save(protocol_path, protocol)
            command = {
                "split": split,
                "method_id": method,
                "argv": [
                    sys.executable,
                    "-m",
                    "benchmarks.research.runner",
                    "--dataset",
                    str(panel),
                    "--protocol",
                    str(protocol_path),
                    "--output",
                    str(output),
                ],
            }
            plan["commands"].append(command)
            runs = []
            for index, query in enumerate(queries):
                for seed in (42, 43, 44):
                    directory = output / method / f"query-{index:04d}-seed-{seed}"
                    directory.mkdir(parents=True)
                    chosen = (
                        method == "pipeline_current"
                        and split == "scaffold"
                        and index == 0
                        and seed in (42, 43)
                    )
                    candidates = (
                        [{"smiles": "c1ccncc1"}, {"smiles": "definitely-invalid"}, {"smiles": "c1ccsc1"}]
                        if chosen
                        else []
                    )
                    response = {
                        "candidates": candidates,
                        "evaluations": len(candidates),
                        "evaluation_unit": sorted(assess.UNITS[method])[0],
                        "metadata": {},
                    }
                    if method == "crem_train_partial_exo":
                        response["metadata"] = {"raw_products": 0, "decisions": [], "postfilter_counts": {}}
                        (directory / "response.raw.jsonl").write_text("")
                    elif method == "mol2mol_pretrained":
                        save(
                            directory / "mol2mol-raw.json",
                            {"rows": [], "seed": seed, "input_smiles": "c1ccccc1"},
                        )
                        save(directory / "mol2mol-selection.json", [])
                        response["metadata"] = {
                            "raw_sha256": assess.digest(directory / "mol2mol-raw.json"),
                            "selection_sha256": assess.digest(directory / "mol2mol-selection.json"),
                            "selection_counts": {},
                        }
                    save(directory / "response.json", response)
                    resolved = {**method_object, "resolved_assets": {}}
                    request = public_request(dataset, query, resolved, seed, protocol["budget"], directory)
                    request.update(deadline_monotonic=1000000.0, finalization_margin_seconds=2.0)
                    save(directory / "request.json", request)
                    status = "success" if candidates else "no_candidates"
                    evaluation = evaluate_candidates(
                        dataset, query, candidates, status, evidence_verifier({})
                    )
                    run = {
                        "method_id": method,
                        "query_id": query["id"],
                        "target_id": f"T{index // 2}",
                        "seed": seed,
                        "status": status,
                        "response": response,
                        "response_sha256": assess.digest(directory / "response.json"),
                        "evaluation": evaluation,
                        "resources": {
                            "wall_seconds": 1.0,
                            "sampled_peak_rss_mb": 100.0,
                            "sampled_cpu_seconds_lower_bound": 0.5,
                            "includes_import_and_model_loading": True,
                        },
                    }
                    save(directory / "evaluation.json", run)
                    runs.append(run)
                    grid.append(
                        {
                            "split": split,
                            "method_id": method,
                            "query_id": query["id"],
                            "target_id": f"T{index // 2}",
                            "seed": seed,
                            "query_index": index,
                            "output_dir": str(directory),
                        }
                    )
            save(output / "runs.json", runs)
            save(output / "campaign.json", {"status": "completed", "attempted_runs": len(runs)})
            save(
                output / "manifest.json",
                {
                    "protocol_sha256": assess.digest(protocol_path),
                    "dataset": {"dataset_sha256": assess.digest(panel)},
                    "code_sha256": plan["source"]["files_sha256"],
                    "python": sys.version,
                    "versions": {p: packages[p] for p in ("rdkit", "numpy", "scipy", "psutil")},
                },
            )
    common = {"conformers": assess.preparation.CONFORMERS, "seeds": [42, 43, 44]}
    addendum_path = tmp_path / "retention-addendum.json"
    addendum = {
        "candidate_annotations": {name: "synthetic fixture definition" for name in assess.RETENTION_FIELDS},
        "eligibility_ranking_budgets_primary_unchanged": True,
    }
    save(addendum_path, addendum)
    addendum_sha = assess.digest(addendum_path)
    monkeypatch.setattr(assess.preparation, "RETENTION_ADDENDUM_SHA", addendum_sha)
    plan["input_files_sha256"] = {str(addendum_path): addendum_sha}
    analysis = {
        "design": "synthetic fixture",
        "secondary_scaffold_retention": addendum,
        "parent_scaffold_retention_addendum": {"artifact": str(addendum_path), "sha256": addendum_sha},
    }
    save(root / "common-3d-plan.json", common)
    save(root / "execution-plan.json", plan)
    save(root / "grid.json", grid)
    save(root / "analysis-plan.json", analysis)
    monkeypatch.setattr(assess, "verify_science", lambda study: None)
    return {
        "root": root,
        "anchor": "parent-fixture-study-anchor",
        "plan": plan,
        "grid": grid,
        "common": common,
        "analysis": analysis,
        "tmp": tmp_path,
    }


def test_complete_270_recompute_raw_ranks_and_cache_coverage(completed, api):
    records, campaigns = assess.collect_generation(completed)
    assert len(records) == 270 and len(campaigns) == 9
    rebuilt, occurrences = assess.recompute(completed, records, api)
    assert len(occurrences) == 4
    assert [o["rank"] for o in occurrences] == [1, 3, 1, 3]
    first = rebuilt[0]["evaluation"]["metrics"]
    assert first["candidate_count"] == 3 and first["valid_candidate_count"] == 2
    assert first["active_scaffold_recall_proxy"] == 1
    assert first["active_scaffold_precision_proxy"] == 0.5
    assert first["exact_known_active_matches"] == first["known_inactive_family_matches"] == 1
    generation = completed["tmp"] / "frozen"
    generation_sha = assess.freeze_generation(completed, generation)
    scores = completed["tmp"] / "common"
    scores_sha = assess.score_generation(completed, generation, generation_sha, scores, api)
    assert len(api.calls) == 2  # four occurrences share two exact pairs
    observed = assess.load_scores(completed, scores, scores_sha, generation_sha, occurrences, api)
    assert len(observed) == 4 and sum(s["all_3d_seeds_succeeded"] for s in observed.values()) == 2
    metrics = assess.attempt_metrics(rebuilt, occurrences, observed)[0]["metrics"]
    assert metrics["attempt_with_hop"] == 1
    assert metrics["common_3d_output_slot_coverage"] == pytest.approx(1 / 3)
    assert metrics["common_3d_failed_candidates"] == 1
    assert metrics["common_hit_shape_mean_success_conditional"] == 0.4
    assert metrics["no_route_supplied_occurrences"] == 2
    assert metrics["route_rejected_occurrences"] == metrics["verified_rule_routes"] == 0
    failed = assess.attempt_metrics(rebuilt, occurrences, observed)[2]["metrics"]
    assert failed["attempt_with_hop"] == 0
    assert failed["hit_bit_tanimoto_mean_conditional"] is None


@pytest.mark.parametrize(
    "alter", ["missing", "duplicate", "old_method", "response_hash", "evaluation", "source", "runtime"]
)
def test_generation_gate_rejects_tampering(completed, alter):
    command = completed["plan"]["commands"][0]
    directory = Path(command["argv"][-1])
    runs = assess.read_json(directory / "runs.json")
    if alter == "missing":
        runs.pop()
    elif alter == "duplicate":
        runs[-1] = runs[0]
    elif alter == "old_method":
        runs[0]["method_id"] = "pipeline"
    elif alter == "response_hash":
        attempt = Path(completed["grid"][0]["output_dir"])
        (attempt / "response.json").write_text("{}")
    elif alter == "evaluation":
        runs[0]["evaluation"]["metrics"]["hop_families"] = 3
    else:
        manifest = assess.read_json(directory / "manifest.json")
        if alter == "source":
            manifest["code_sha256"] = {"old_source.py": "b" * 64}
        else:
            manifest["versions"]["rdkit"] = "2026.9.1"
        save(directory / "manifest.json", manifest)
    save(directory / "runs.json", runs)
    with pytest.raises(ValueError):
        assess.collect_generation(completed)


def test_recomputed_metrics_catch_self_consistent_report_tampering(completed, api):
    directory = Path(completed["plan"]["commands"][0]["argv"][-1])
    runs = assess.read_json(directory / "runs.json")
    runs[0]["evaluation"]["metrics"]["active_scaffold_recall_proxy"] = 0.123
    save(directory / "runs.json", runs)
    save(Path(completed["grid"][0]["output_dir"]) / "evaluation.json", runs[0])
    records, _ = assess.collect_generation(completed)
    with pytest.raises(ValueError, match="Recomputed chemical"):
        assess.recompute(completed, records, api)


def test_private_request_fields_rejected_after_snapshot(completed, api):
    path = Path(completed["grid"][0]["output_dir"]) / "request.json"
    request = assess.read_json(path)
    request["query"]["label"] = "active"
    save(path, request)
    records, _ = assess.collect_generation(completed)
    with pytest.raises(ValueError, match="private-field"):
        assess.recompute(completed, records, api)


def test_failed_partial_response_never_becomes_accepted(completed, api):
    campaign = Path(completed["plan"]["commands"][0]["argv"][-1])
    runs = assess.read_json(campaign / "runs.json")
    first = runs[0]
    dataset, _ = load_dataset(completed["tmp"] / "scaffold-panel.json")
    first.update(
        status="timeout",
        response={},
        evaluation=evaluate_candidates(dataset, dataset["queries"][0], [], "timeout", evidence_verifier({})),
    )
    save(campaign / "runs.json", runs)
    save(Path(completed["grid"][0]["output_dir"]) / "evaluation.json", first)
    records, _ = assess.collect_generation(completed)
    rebuilt, occurrences = assess.recompute(completed, records, api)
    assert rebuilt[0]["evaluation"]["metrics"]["hop_families"] == 0
    assert rebuilt[0]["raw_accounting"]["raw_returned_candidates"] == 3
    assert not rebuilt[0]["raw_accounting"]["accepted_response"]
    assert len(occurrences) == 2


def test_native_raw_selection_counts_are_not_trusted(completed, api):
    slot = next(g for g in completed["grid"] if g["method_id"] == "mol2mol_pretrained")
    directory = Path(slot["output_dir"])
    response = assess.read_json(directory / "response.json")
    response["metadata"]["selection_counts"] = {"selected": 5}
    save(directory / "response.json", response)
    campaign = next(
        Path(c["argv"][-1])
        for c in completed["plan"]["commands"]
        if c["split"] == slot["split"] and c["method_id"] == slot["method_id"]
    )
    runs = assess.read_json(campaign / "runs.json")
    runs[0].update(response=response, response_sha256=assess.digest(directory / "response.json"))
    save(directory / "evaluation.json", runs[0])
    save(campaign / "runs.json", runs)
    records, _ = assess.collect_generation(completed)
    with pytest.raises(ValueError, match="selection counts differ"):
        assess.recompute(completed, records, api)


def test_parent_snapshot_hash_catches_generation_change(completed):
    frozen = completed["tmp"] / "frozen"
    sha = assess.freeze_generation(completed, frozen)
    (Path(completed["grid"][0]["output_dir"]) / "new-debug.log").write_text("added after parent review")
    with pytest.raises(ValueError, match="changed after parent-pinned"):
        assess.load_generation(completed, frozen, sha)


def test_hierarchy_uses_query_then_target_means_and_preserves_null():
    slots = [
        {"query_id": q, "target_id": t, "seed": s}
        for q, t in (("q1", "A"), ("q2", "A"), ("q3", "B"))
        for s in (42, 43, 44)
    ]
    values = {
        ("q1", 42): 1,
        ("q1", 43): 1,
        ("q1", 44): 1,
        ("q2", 42): 0,
        ("q3", 42): 0,
        ("q3", 43): 0,
        ("q3", 44): 0,
    }
    report = assess.hierarchical(values, slots, cluster_interval)
    # q1=1, q2=0 -> target A=.5; target B=0 -> target macro=.25.
    assert report["target_macro"]["mean"] == 0.25
    assert report["target_macro"]["ci95"] == [0, 0.5]
    assert report["availability"]["missing_slots"] == 2
    assert report["availability"]["observed_seeds_per_query"]["q2"] == 1
    absent = assess.hierarchical({}, slots, cluster_interval)
    assert absent["target_macro"]["mean"] is None and absent["target_macro"]["ci95"] is None
    singleton = assess.hierarchical({("q1", 42): 0.7}, slots[:3], cluster_interval)
    assert singleton["target_macro"]["mean"] == 0.7 and singleton["target_macro"]["ci95"] is None


def test_paired_conditional_missing_and_infrastructure_bounds_are_hand_calculated():
    rows, grid = [], []
    for split in assess.preparation.SPLITS:
        for method in assess.preparation.METHODS:
            for seed in (42, 43, 44):
                grid.append(
                    {"split": split, "method_id": method, "query_id": "q", "target_id": "T", "seed": seed}
                )
                primary = method == "pipeline_current"
                value = (1 if seed == 42 else 0) if primary else 0
                conditional = {42: 10, 43: None, 44: 30}[seed] if primary else {42: 1, 43: 5, 44: None}[seed]
                status = assess.INFRA if primary and seed == 43 else ("success" if seed == 42 else "timeout")
                rows.append(
                    {
                        **grid[-1],
                        "status": status,
                        "metrics": {"attempt_with_hop": value, "conditional": conditional},
                    }
                )
    result = assess.comparisons(rows, grid, cluster_interval)["scaffold"]
    primary = result["methods"]["pipeline_current"]
    assert primary["metrics"]["attempt_with_hop"]["target_macro"]["mean"] == pytest.approx(1 / 3)
    assert primary["infrastructure_counterfactual_bounds"]["attempt_with_hop"] == pytest.approx(
        [1 / 3, 2 / 3]
    )
    assert (
        primary["complete_case_excluding_infrastructure_only"]["attempt_with_hop"]["target_macro"]["mean"]
        == 0.5
    )
    pair = result["paired_pipeline_current_minus_comparator"]["crem_train_partial_exo"]
    # Paired observed seed42 gives 10-1=9; subtracting marginal means would give 17.
    assert pair["metrics"]["conditional"]["target_macro"]["mean"] == 9
    assert pair["metrics"]["conditional"]["availability"]["observed_slots"] == 1
    assert pair["infrastructure_counterfactual_bounds"]["attempt_with_hop"] == pytest.approx([1 / 3, 2 / 3])
    assert (
        pair["complete_case_excluding_either_infrastructure"]["attempt_with_hop"]["target_macro"]["mean"]
        == 0.5
    )
    assert not result["counterfactual_bounds_are_confidence_intervals"]


@pytest.mark.parametrize("change", ["missing_seed", "median", "config", "pair_error", "success_flag"])
def test_common_score_seed_provenance_and_failure_flags(completed, api, change):
    records, _ = assess.collect_generation(completed)
    _, occurrences = assess.recompute(completed, records, api)
    occurrence = occurrences[0]
    payload = assess.cache_payload(completed, occurrence)
    score = assess.score_pair(payload, api)
    if change == "missing_seed":
        score["trials"].pop()
    elif change == "median":
        score["scores"]["hit_shape"] = 0.9
    elif change == "config":
        score["payload"] = copy.deepcopy(payload)
        score["payload"]["conformers"]["num_conformers"] = 2
    elif change == "pair_error":
        score["pair_error"] = "hidden failure"
    else:
        score["all_3d_seeds_succeeded"] = False
    with pytest.raises(ValueError):
        assess.validate_score(score, payload, occurrence["candidate"], api)


def test_common_score_join_rejects_duplicate_and_missing_rank(completed, api):
    frozen, scores = completed["tmp"] / "frozen", completed["tmp"] / "scores"
    generation_sha = assess.freeze_generation(completed, frozen)
    scores_sha = assess.score_generation(completed, frozen, generation_sha, scores, api)
    records, _ = assess.collect_generation(completed)
    _, occurrences = assess.recompute(completed, records, api)
    rows = assess.read_json(scores / "occurrences.json")
    rows[-1] = rows[0]
    save(scores / "occurrences.json", rows)
    manifest = assess.read_json(scores / "manifest.json")
    manifest["files_sha256"]["occurrences.json"] = assess.digest(scores / "occurrences.json")
    save(scores / "manifest.json", manifest)
    with pytest.raises(ValueError, match="Missing, duplicate"):
        assess.load_scores(
            completed, scores, assess.digest(scores / "manifest.json"), generation_sha, occurrences, api
        )
    with pytest.raises(ValueError, match="parent-pinned"):
        assess.verify_archive(scores, scores_sha)


def test_assessment_cannot_write_into_campaign(completed):
    with pytest.raises(ValueError, match="cannot modify"):
        assess.new_output(Path(completed["plan"]["commands"][0]["argv"][-1]) / "new-analysis", completed)


def test_original_630_cannot_substitute_for_supplementary_grid(completed, monkeypatch):
    monkeypatch.setattr(assess.preparation, "verify_prepared", lambda *args: None)
    assert len(assess.load_study(completed["root"], completed["anchor"])["grid"]) == 270
    grid = assess.read_json(completed["root"] / "grid.json")
    for row in grid:
        if row["method_id"] == "pipeline_current":
            row["method_id"] = "pipeline"
    save(completed["root"] / "grid.json", grid)
    with pytest.raises(ValueError, match="exact 270-slot"):
        assess.load_study(completed["root"], completed["anchor"])


@pytest.mark.parametrize(
    "query,candidate,retained,generic,expansion,not_contained,q_atoms,c_atoms",
    [
        ("c1ccccc1", "c1ccc(-c2ccccc2)cc1", True, True, True, False, 6, 12),
        ("c1ccccc1", "c1ccc2ccccc2c1", True, True, True, False, 6, 10),
        ("c1ccccc1", "c1ccncc1", False, True, False, True, 6, 6),
        ("C1CCCCC1", "C1CCCC1", False, False, False, True, 6, 5),
        ("C1CCCCC1", "C1CCCCCC1", False, False, False, True, 6, 7),
        ("c1ccccc1", "C1CCCCC1", False, True, False, True, 6, 6),
        ("c1ccc(-c2ccccc2)cc1", "c1ccccc1", False, False, False, True, 12, 6),
        ("Cc1ccccc1", "Nc1ccccc1", True, True, False, False, 6, 6),
    ],
)
def test_secondary_scaffold_annotations_use_actual_graphs(
    api, query, candidate, retained, generic, expansion, not_contained, q_atoms, c_atoms
):
    observed = assess.scaffold_retention(query, candidate, api)
    assert observed["query_scaffold_subgraph_retained"] is retained
    assert observed["query_generic_ring_skeleton_subgraph_retained"] is generic
    assert observed["core_retaining_scaffold_expansion"] is expansion
    assert observed["query_scaffold_not_contained"] is not_contained
    assert observed["query_scaffold_heavy_atoms"] == q_atoms
    assert observed["candidate_scaffold_heavy_atoms"] == c_atoms
    assert observed["atom_bond_status"] == observed["generic_status"] == "assessed"


@pytest.mark.parametrize("query,candidate", [("CC", "c1ccccc1"), ("c1ccccc1", "CC"), ("invalid", "c1ccccc1")])
def test_missing_scaffolds_are_unassessed_not_false(api, query, candidate):
    observed = assess.scaffold_retention(query, candidate, api)
    assert all(observed[name] is None for name in assess.RETENTION_FIELDS)
    assert observed["atom_bond_status"] == observed["generic_status"] == "unassessed"
    assert observed["reason"]


def test_scaffold_annotation_does_not_redefine_hop_or_primary(completed, api):
    records, _ = assess.collect_generation(completed)
    rebuilt, occurrences = assess.recompute(completed, records, api)
    scores = {}
    for occurrence in occurrences:
        payload = assess.cache_payload(completed, occurrence)
        scores[assess.occurrence_key(occurrence)] = assess.score_pair(payload, api)
    metrics = assess.attempt_metrics(rebuilt, occurrences, scores)[0]["metrics"]
    assert metrics["hop_families"] == 2 and metrics["attempt_with_hop"] == 1
    assert metrics["query_scaffold_not_contained_occurrences"] == 2
    assert metrics["query_scaffold_not_contained_unique_families"] == 2
    assert metrics["query_generic_ring_skeleton_subgraph_retained_occurrences"] == 1
    assert metrics["core_retaining_scaffold_expansion_occurrences"] == 0
    assert metrics["attempt_with_query_scaffold_not_contained"] == 1


def test_full_analysis_keeps_all_270_denominators_and_secondary_annotations(completed, api):
    frozen, scores, output = (completed["tmp"] / name for name in ("frozen", "scores", "analysis"))
    generation_sha = assess.freeze_generation(completed, frozen)
    scores_sha = assess.score_generation(completed, frozen, generation_sha, scores, api)
    analysis_sha = assess.analyze_generation(
        completed, frozen, generation_sha, scores, scores_sha, output, api
    )
    manifest = assess.verify_archive(output, analysis_sha)
    assert manifest["complete_requires_parent_review"] is True
    assert len(assess.read_json(output / "attempts.json")) == 270
    comparison = assess.read_json(output / "comparison.json")
    primary = comparison["scaffold"]["methods"]["pipeline_current"]
    assert primary["attempts"] == 36
    assert primary["metrics"]["attempt_with_hop"]["target_macro"]["mean"] == pytest.approx(2 / 36)
    assert primary["metrics"]["common_3d_output_slot_coverage"]["target_macro"]["mean"] == pytest.approx(
        2 / 108
    )
    assert primary["metrics"]["common_hit_shape_mean_success_conditional"]["target_macro"]["mean"] == 0.4
    assert (
        primary["metrics"]["common_hit_shape_mean_success_conditional"]["availability"]["observed_slots"] == 2
    )
    summary = assess.read_json(output / "scaffold-retention-summary.json")["scaffold"]["pipeline_current"]
    assert summary["attempts"] == 36 and summary["valid_candidate_occurrences"] == 4
    annotation = summary["annotations"]["query_scaffold_not_contained"]
    assert annotation == {
        "qualifying_occurrences": 4,
        "unique_qualifying_families": 2,
        "assessed_occurrences": 4,
        "unassessed_occurrences": 0,
        "attempts_with_qualifying_output": 2,
    }
    with pytest.raises(ValueError, match="new"):
        assess.analyze_generation(completed, frozen, generation_sha, scores, scores_sha, output, api)


def test_pair_error_preserves_all_unassessed_seed_records(completed, api):
    records, _ = assess.collect_generation(completed)
    _, occurrences = assess.recompute(completed, records, api)
    payload = assess.cache_payload(completed, occurrences[0])

    def fail(*args):
        raise RuntimeError("synthetic pair-level evaluator failure")

    api.similarity = fail
    scored = assess.score_pair(payload, api)
    assess.validate_score(scored, payload, occurrences[0]["candidate"], api)
    assert [t["seed"] for t in scored["trials"]] == [42, 43, 44]
    assert all(t["status"] == "unassessed_due_to_pair_error" for t in scored["trials"])
    assert scored["scores"] == {} and scored["all_3d_seeds_succeeded"] is False


def test_interrupted_common_pass_preserves_only_completed_pairs(completed, api):
    original = api.similarity
    calls = 0

    def interrupted(*args):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise KeyboardInterrupt("synthetic interruption after first pair")
        return original(*args)

    api.similarity = interrupted
    generation = completed["tmp"] / "frozen"
    anchor = assess.freeze_generation(completed, generation)
    output = completed["tmp"] / "interrupted-common"
    with pytest.raises(KeyboardInterrupt):
        assess.score_generation(completed, generation, anchor, output, api)
    assert len(list((output / "cache").glob("*.json"))) == 1
    cached = assess.read_json(next((output / "cache").glob("*.json")))
    assert cached["all_3d_seeds_succeeded"] is True
    assert len(cached["trials"]) == 3
    assert assess.read_json(output / "context.json")["complete"] is False
    assert not (output / "manifest.json").exists()
    assert not (output / "coverage.json").exists()
    with pytest.raises(ValueError, match="no overwrite or resume"):
        assess.score_generation(completed, generation, anchor, output, api)


def test_top_level_helper_imports_must_come_from_pinned_source(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    helper = root / "utils.py"
    helper.write_text("PINNED = True\n")
    source = {"root": str(root), "files_sha256": {"utils.py": assess.digest(helper)}}
    assess.verify_import_origins(source, {"utils": SimpleNamespace(__file__=str(helper))})
    elsewhere = tmp_path / "utils.py"
    elsewhere.write_bytes(helper.read_bytes())
    with pytest.raises(ValueError, match="outside pinned scientific source"):
        assess.verify_import_origins(source, {"utils": SimpleNamespace(__file__=str(elsewhere))})
    helper.write_text("PINNED = False\n")
    with pytest.raises(ValueError, match="outside pinned scientific source"):
        assess.verify_import_origins(source, {"utils": SimpleNamespace(__file__=str(helper))})
