"""Adversarial protocol tests plus actual catalog/reaction adapter execution."""

import copy
import json
from pathlib import Path
import sys

import pytest

from benchmarks.research.adapters import nearest_catalog, reaction_enumeration, synllama_only
from benchmarks.research.evaluation import aggregate_runs, binary_ranking_metrics, evaluate_candidates
from benchmarks.research.import_chembl import activity_label
from benchmarks.research.protocol import audit_inventory, digest, load_dataset, write_json
from benchmarks.research.runner import (
    bounded_method,
    evidence_verifier,
    load_protocol,
    public_request,
    run_benchmark,
)


@pytest.fixture
def panel(tmp_path):
    source = tmp_path / "source.txt"
    source.write_text("Synthetic test fixture; no experimental claim.\n")
    data = {
        "schema_version": 1,
        "dataset_id": "synthetic",
        "scope": "synthetic_fixture",
        "sources": [
            {
                "id": "fixture",
                "url": "https://example.invalid/fixture",
                "artifact": str(source),
                "sha256": digest(source),
            }
        ],
        "compounds": [
            {
                "id": "q",
                "smiles": "CC(=O)Nc1ccccc1",
                "target_id": "T",
                "label": "active",
                "split": "train",
                "source_id": "fixture",
                "assay_id": "fixture",
                "label_basis": "synthetic fixture",
            },
            {
                "id": "a",
                "smiles": "CC(=O)Nc1ccncc1",
                "target_id": "T",
                "label": "active",
                "split": "test",
                "source_id": "fixture",
                "assay_id": "fixture",
                "label_basis": "synthetic fixture",
            },
            {
                "id": "n",
                "smiles": "CC(=O)Nc1ccsc1",
                "target_id": "T",
                "label": "inactive",
                "split": "test",
                "source_id": "fixture",
                "assay_id": "fixture",
                "label_basis": "synthetic fixture",
            },
            {
                "id": "d",
                "smiles": "CC(=O)NC1CCCC1",
                "target_id": "T",
                "label": "decoy",
                "split": "test",
                "source_id": "fixture",
            },
        ],
        "queries": [{"id": "query", "compound_id": "q"}],
        "split": {"kind": "scaffold"},
    }
    path = tmp_path / "panel.json"
    write_json(path, data)
    return data, path


def test_query_is_public_input_not_required_to_be_test(panel):
    data, path = panel
    _, report = load_dataset(path)
    assert report["query_coverage"][0]["novel_active_references"] == 1
    assert not report["scientific_evidence"]
    method = {"id": "baseline", "kind": "nearest_catalog", "resolved_assets": {}}
    request = public_request(data, data["queries"][0], method, 42, {}, path.parent)
    serialized = json.dumps(request)
    assert data["compounds"][1]["smiles"] not in serialized
    assert "assay_id" not in serialized and "label" not in serialized


def test_source_tampering_rejected(panel):
    data, path = panel
    Path(data["sources"][0]["artifact"]).write_text("changed")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        load_dataset(path)


def test_chemical_family_leakage_includes_protonation(panel):
    data, path = panel
    data["compounds"][0]["smiles"] = "CC(=O)Nc1cc[nH+]cc1"
    write_json(path, data)
    with pytest.raises(ValueError, match="Cross-partition family leakage"):
        load_dataset(path)


def test_scaffold_leakage_and_target_partition_are_explicit(panel):
    data, path = panel
    data["compounds"][1]["smiles"] = "CCC(=O)Nc1ccccc1"
    write_json(path, data)
    with pytest.raises(ValueError, match="split_group leakage"):
        load_dataset(path)
    data["compounds"][1]["smiles"] = "CC(=O)Nc1ccncc1"
    data["split"]["kind"] = "target"
    write_json(path, data)
    with pytest.raises(ValueError, match="split_group leakage"):
        load_dataset(path)
    data["compounds"][0]["split"] = "query"
    write_json(path, data)
    load_dataset(path)


def test_publication_year_not_fabricated_as_exact_date(panel):
    data, path = panel
    data["split"] = {"kind": "time", "cutoff_year": 2020, "date_kind": "publication_year"}
    for row in data["compounds"]:
        row["publication_year"] = 2010 if row["split"] == "train" else 2022
    write_json(path, data)
    load_dataset(path)
    data["compounds"][1]["publication_year"] = 2019
    write_json(path, data)
    with pytest.raises(ValueError, match="on/before"):
        load_dataset(path)


def test_unavailable_or_partial_inventory_never_proves_clean_holdout(panel, tmp_path):
    data, _ = panel
    unknown = audit_inventory(
        {"kind": "training", "status": "unknown", "reason": "Upstream inventory not published"},
        data,
        tmp_path,
    )
    assert unknown["status"] == "unknown"
    inv = tmp_path / "inventory.json"
    write_json(inv, ["C1CCC1"])
    spec = {
        "kind": "training",
        "status": "available",
        "coverage": "partial",
        "artifact": str(inv),
        "sha256": digest(inv),
    }
    audit = audit_inventory(spec, data, tmp_path)
    assert audit["status"] == "no_overlap_in_declared_inventory"
    assert not audit["complete_inventory_no_detected_overlap"]
    write_json(inv, ["CC(=O)Nc1cc[nH+]cc1"])
    spec["sha256"] = digest(inv)
    audit = audit_inventory(spec, data, tmp_path)
    assert audit["status"] == "overlap_detected"
    assert any(r["compound_id"] == "a" and r["family_overlap"] for r in audit["matches"])


def test_inventory_unassessable_records_are_partial_with_original_indices(panel, tmp_path):
    import hashlib

    data, _ = panel
    inv = tmp_path / "inventory.json"
    values = ["C1CCC1", "[1*+]c1cccc2c1NCC21CCCC1", "not_smiles", {}]
    write_json(inv, values)
    report = audit_inventory(
        {
            "kind": "training",
            "status": "available",
            "coverage": "complete",
            "artifact": str(inv),
            "sha256": digest(inv),
        },
        data,
        tmp_path,
    )
    assert report["records"] == 1 and report["selected_records"] == 4
    assert report["declared_coverage"] == "complete" and report["coverage"] == "partial"
    assert report["unassessed_records"] == 3
    assert [row["record_index"] for row in report["unassessed"]] == [1, 2, 3]
    assert (
        report["unassessed"][0]["input_sha256"] == hashlib.sha256(json.dumps(values[1]).encode()).hexdigest()
    )
    assert not report["complete_inventory_no_detected_overlap"]
    assert "generalization_claim_supported" not in report
    assert report["status"] == "no_overlap_in_assessed_subset"


def test_invalid_duplicates_and_decoys_cannot_inflate_hop_recovery(panel):
    data, _ = panel
    candidates = [
        {"smiles": data["compounds"][0]["smiles"]},
        {"smiles": data["compounds"][1]["smiles"], "score_samples": [0.7, 0.99, 0.98]},
        {"smiles": "CC(=O)Nc1cc[nH+]cc1"},
        {"smiles": data["compounds"][2]["smiles"]},
        {"smiles": data["compounds"][3]["smiles"]},
        {"smiles": "invalid"},
    ]
    report = evaluate_candidates(data, data["queries"][0], candidates)
    metric = report["metrics"]
    assert metric["candidate_count"] == 6 and metric["valid_candidate_count"] == 5
    assert metric["duplicate_chemical_states"] == 1
    assert metric["unchanged_query_families"] == 1
    assert metric["recovered_active_scaffolds"] == 1
    assert metric["known_inactive_family_matches"] == 1
    assert metric["decoy_family_matches"] == 1
    assert report["candidates"][1]["score_stability"]["range"] == pytest.approx(0.29)
    failed = evaluate_candidates(data, data["queries"][0], candidates, "timeout")
    assert failed["metrics"]["active_scaffold_recall_proxy"] == 0
    assert failed["metrics"]["candidate_count"] == 0


def test_target_resampling_not_candidate_or_seed_pseudoreplication(panel):
    data, _ = panel
    positive = evaluate_candidates(data, data["queries"][0], [{"smiles": data["compounds"][1]["smiles"]}])
    failure = evaluate_candidates(data, data["queries"][0], [], "timeout")
    runs = [
        {
            "method_id": "m",
            "query_id": "q1",
            "target_id": "T1",
            "seed": i,
            "status": "success",
            "evaluation": positive,
            "resources": {"wall_seconds": 1},
        }
        for i in range(20)
    ]
    runs.append(
        {
            "method_id": "m",
            "query_id": "q2",
            "target_id": "T2",
            "seed": 1,
            "status": "timeout",
            "evaluation": failure,
            "resources": {"wall_seconds": 1},
        }
    )
    result = aggregate_runs(runs)["m"]
    recall = result["target_macro_metrics"]["active_scaffold_recall_proxy"]
    assert result["attempts"] == 21 and result["status_counts"]["timeout"] == 1
    assert recall["mean"] == 0.5 and recall["independent_targets"] == 2


def test_controls_separate_discrimination_from_probability_calibration():
    result = binary_ranking_metrics([("active", 0.9), ("inactive", 0.1), ("decoy", 1.0)])
    assert result["roc_auc"] == 1 and result["average_precision"] == 1
    assert result["inactive_count"] == 1
    assert result["probability_calibration"] == "not_applicable_to_similarity_scores"
    assert binary_ranking_metrics([("active", 1), ("decoy", 0)])["roc_auc"] is None
    ties = binary_ranking_metrics([("active", 1), ("inactive", 1)])
    assert ties["roc_auc"] == ties["average_precision"] == 0.5


@pytest.mark.parametrize(
    "relation,value,expected",
    [
        ("<", 100, "active"),
        ("<=", 1000, "active"),
        (">", 10000, "inactive"),
        (">=", 20000, "inactive"),
        (">", 100, None),
        ("<", 20000, None),
        ("=", 3000, None),
        ("=", 0, None),
        ("~", 20, None),
    ],
)
def test_censored_assay_labels(relation, value, expected):
    row = {
        "standard_units": "nM",
        "standard_type": "IC50",
        "standard_flag": 1,
        "assay_type": "B",
        "target_organism": "Homo sapiens",
        "standard_value": value,
        "standard_relation": relation,
    }
    assert activity_label(row)[0] == expected


@pytest.fixture
def assets(tmp_path):
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text(
        "id\tsmiles\tsource\tsnapshot\na\tCC(=O)O\tfixture\t2026-01-01\np\tNc1ccncc1\tfixture\t2026-01-01\nb\tCCC(=O)O\tfixture\t2026-01-01\n"
    )
    rules = tmp_path / "rules.json"
    write_json(
        rules,
        {
            "schema_version": 1,
            "rules": [
                {
                    "id": "amide",
                    "smarts": "[C:1](=[O:2])[O;H1:3].[N;H1,H2:4]>>[C:1](=[O:2])[N:4]",
                    "source": "test fixture",
                }
            ],
        },
    )
    return {"catalog": str(catalog), "rules": str(rules)}


def request(assets, out):
    return {
        "query": {"id": "q", "smiles": "CC(=O)Nc1ccccc1", "target_id": "T"},
        "assets": assets,
        "seed": 42,
        "budget": {"wall_seconds": 20, "max_candidates": 10, "max_evaluations": 20},
        "settings": {},
        "output_dir": str(out),
    }


def test_actual_catalog_and_reaction_adapters_verify_routes(assets, tmp_path):
    req = request(assets, tmp_path)
    nearest = nearest_catalog(req)
    assert nearest["evaluations"] == 3 and len(nearest["candidates"]) == 3
    result = reaction_enumeration(req)
    assert result["metadata"]["reaction_combinations"] == 2
    assert {r["smiles"] for r in result["candidates"]} == {"CC(=O)Nc1ccncc1", "CCC(=O)Nc1ccncc1"}
    verify = evidence_verifier(assets)
    assert all(verify(row)["level"] == "catalog_and_rule_verified" for row in result["candidates"])
    altered = copy.deepcopy(result["candidates"][0])
    altered["smiles"] = "CCC"
    assert verify(altered)["level"] == "route_rejected"


def test_reaction_rules_prepare_lazily_and_keep_completed_partial_results(assets, tmp_path, monkeypatch):
    import time
    from benchmarks.research import adapters

    rule_path = Path(assets["rules"])
    registry = json.loads(rule_path.read_text())
    template = registry["rules"][0]
    for number in (5, 6):
        registry["rules"].append(
            {**template, "id": f"amide-{number}", "smarts": template["smarts"].replace(":4]", f":{number}]")}
        )
    write_json(rule_path, registry)
    req = request(assets, tmp_path)
    req["budget"]["max_evaluations"] = 1
    result = reaction_enumeration(req)
    assert result["candidates"] and result["evaluations"] == 1
    meta = result["metadata"]
    assert len(meta["initialized_rules"]) == 1 and len(meta["unvisited_rules"]) == 2
    assert len(meta["reaction_combinations_by_rule"]) == 1
    assert len(meta["rule_pool_preparation"]) == 2

    # An actual completed chemical tuple is retained when time expires before
    # another rule can be initialized. No mocked candidate or route is returned.
    real_apply = adapters.Rules.apply
    req["budget"]["max_evaluations"] = 20
    req["deadline_monotonic"] = time.monotonic() + 30

    def expire_after_reaction(self, *args, **kwargs):
        products = real_apply(self, *args, **kwargs)
        req["deadline_monotonic"] = time.monotonic() - 1
        return products

    monkeypatch.setattr(adapters.Rules, "apply", expire_after_reaction)
    result = reaction_enumeration(req)
    assert result["candidates"] and result["evaluations"] == 1
    assert result["metadata"]["anytime_stop"] == "wall_margin"
    assert len(result["metadata"]["unvisited_rules"]) == 2
    assert all(
        evidence_verifier(assets)(row)["level"] == "catalog_and_rule_verified" for row in result["candidates"]
    )


def test_synllama_adapter_reconstructs_and_ranks_actual_routes(assets, tmp_path, monkeypatch):
    # Only expensive stochastic generation is substituted; graph resolution,
    # catalog lookup, reaction execution, verification and ranking are real.
    from chembounce_synthesis import synllama

    class Provider:
        provenance = {"provider": "test fixture; not a model benchmark"}

        def __init__(self, **kwargs):
            pass

        def propose(self, smiles, deadline_monotonic=None):
            return [
                {
                    "response": {
                        "reactions": [
                            {
                                "reaction_template": "amide",
                                "reactants": ["CC(=O)O", "Nc1ccncc1"],
                                "product": "CC(=O)Nc1ccncc1",
                            }
                        ],
                        "building_blocks": ["CC(=O)O", "Nc1ccncc1"],
                    }
                }
            ]

    monkeypatch.setattr(synllama, "SynLlama", Provider)
    req = request(assets | {"model_directory": "fixture", "prompt": "fixture"}, tmp_path)
    result = synllama_only(req)
    assert result["candidates"] and result["candidates"][0]["smiles"] == "CC(=O)Nc1ccncc1"
    assert result["candidates"][0]["route"]["schema_version"] == 1


def test_subprocess_timeout_and_over_budget_are_terminal_failures(tmp_path):
    req = request({}, tmp_path / "timeout")
    req["budget"]["wall_seconds"] = 0.15
    method = {
        "kind": "external",
        "command": [sys.executable, "-c", "import time; time.sleep(10)", "{request}", "{response}"],
    }
    result = bounded_method(method, req, tmp_path / "timeout")
    assert result["status"] == "timeout" and result["resources"]["wall_seconds"] < 3
    req = request({}, tmp_path / "over")
    code = "import json,sys; json.dump({'candidates':[], 'evaluations':21, 'evaluation_unit':'test'},open(sys.argv[2],'w'))"
    method["command"] = [sys.executable, "-c", code, "{request}", "{response}"]
    assert bounded_method(method, req, tmp_path / "over")["status"] == "budget_exceeded"


def test_real_subprocess_benchmark_all_runs_and_sources_recorded(panel, assets, tmp_path):
    _, dataset_path = panel
    protocol = {
        "schema_version": 1,
        "seeds": [42, 43],
        "budget": {"wall_seconds": 30, "max_candidates": 10, "max_evaluations": 20},
        "methods": [
            {
                "id": "enumeration",
                "kind": "reaction_enumeration",
                "assets": {name: {"artifact": path, "sha256": digest(path)} for name, path in assets.items()},
            }
        ],
        "planned_unexecuted_methods": [{"id": "rush", "reason": "ROCS license unavailable"}],
    }
    protocol_path = tmp_path / "protocol.json"
    write_json(protocol_path, protocol)
    load_protocol(protocol_path, panel[0])
    runs = run_benchmark(dataset_path, protocol_path, tmp_path / "runs")
    assert len(runs) == 2 and all(r["status"] == "success" for r in runs)
    assert all(r["evaluation"]["metrics"]["active_scaffold_recall_proxy"] == 1 for r in runs)
    summary = json.loads((tmp_path / "runs/summary.json").read_text())
    assert list(summary) == ["enumeration"]
    assert summary["enumeration"]["target_macro_metrics"]["run_success"]["ci95"] is None


def test_failed_control_bounds_keep_all_labels():
    from benchmarks.research.evaluation import binary_ranking_with_failures

    report = binary_ranking_with_failures(
        [("active", 0.9), ("inactive", 0.1), ("active", None), ("inactive", None)]
    )
    assert report["roc_auc"] == 1
    assert report["roc_auc_bounds"] == [0.25, 1.0]
    assert report["average_precision_bounds"][0] < report["average_precision"]
    assert report["total_active_controls"] == report["total_inactive_controls"] == 2


def test_common_3d_evaluator_uses_real_score_candidate_tuple():
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from benchmarks.research.controls import similarity_scores

    scores, trials = similarity_scores(
        "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1", [42, 43], {"num_conformers": 2, "max_iterations": 1000}
    )
    assert all(t["status"] == "success" for t in trials)
    assert 0 <= scores["hit_shape"] <= 1 and 0 <= scores["hit_pharmacophore"] <= 1
    assert scores["score_stability"]["n"] == 2


def test_training_library_contains_no_labels_or_test_families(panel, tmp_path):
    from benchmarks.research.training_library import export_library

    _, path = panel
    output = tmp_path / "train.tsv"
    report = export_library(path, output)
    assert report["records"] == 1
    assert "active" not in output.read_text() and "inactive" not in output.read_text()
    assert "CC(=O)Nc1ccncc1" not in output.read_text()


def test_full_pipeline_and_no3d_ablation_use_real_generation_and_saved_routes(assets, tmp_path, monkeypatch):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from chembounce import chembounce
    from benchmarks.research.adapters import pipeline

    hit = "CC(=O)Nc1ccccc1"
    library = tmp_path / "scaffolds.smi"
    library.write_text("c1ccncc1\n")
    options = {
        "tanimoto_threshold": 0.0,
        "_search_scf_thr_": 0.0,
        "calculate_shape": False,
        "candidate_thresholds": {},
        "frag_max_n": 3,
        "scaffold_top_n": 3,
    }
    designs, _ = chembounce(
        hit, output_dir=tmp_path / "preflight", scaffold_db_file=library, overall_max_n=3, **options
    )
    assert len(designs)
    proposal = {
        "reactions": [
            {
                "reaction_template": "amide",
                "reactants": ["CC(=O)O", "Nc1ccncc1"],
                "product": "CC(=O)Nc1ccncc1",
            }
        ],
        "building_blocks": ["CC(=O)O", "Nc1ccncc1"],
    }
    proposals = tmp_path / "proposals.json"
    write_json(
        proposals,
        {
            "schema_version": 1,
            "provider": {"provider": "synthetic fixture"},
            "proposals": [{"design": smiles, "response": proposal} for smiles in designs["Final structure"]],
        },
    )
    config_path = tmp_path / "synthesis.json"
    write_json(
        config_path,
        {
            "schema_version": 1,
            "catalog": assets["catalog"],
            "rules": assets["rules"],
            "provider": {"kind": "saved", "path": str(proposals)},
            "search": {"building_block_top_k": 1, "beam_width": 2},
            "conformers": {"num_conformers": 2, "max_iterations": 1000},
            "refinement": {"top_n": 0},
        },
    )
    req = request(
        assets | {"scaffold_db": str(library), "synthesis_config": str(config_path)}, tmp_path / "full"
    )
    req["budget"]["max_candidates"] = 3
    req["settings"] = {"generator_options": options}
    result = pipeline(req)
    assert result["candidates"]
    assert all(
        evidence_verifier(assets)(row)["level"] == "catalog_and_rule_verified" for row in result["candidates"]
    )
    from chembounce_synthesis import conformers

    monkeypatch.setattr(
        conformers, "score_candidate", lambda *a, **kw: pytest.fail("no-3D ablation called 3D scoring")
    )
    req["settings"]["disable_3d"] = True
    req["output_dir"] = str(tmp_path / "no-3d")
    result = pipeline(req)
    assert result["candidates"] and result["metadata"]["ablation"] == "3d_scoring_disabled"


def test_curator_preserves_years_exclusions_and_checks_size_before_normalization(tmp_path, monkeypatch):
    from benchmarks.research import import_chembl

    common = {
        "standard_units": "nM",
        "standard_type": "IC50",
        "standard_flag": 1,
        "assay_type": "B",
        "target_organism": "Homo sapiens",
        "standard_relation": "=",
        "target_chembl_id": "T",
        "assay_chembl_id": "A",
        "molecule_chembl_id": "M",
    }
    entries = [
        common
        | {
            "activity_id": 1,
            "canonical_smiles": "CC(=O)Nc1ccccc1",
            "document_year": 2010,
            "standard_value": "100",
        },
        common
        | {
            "activity_id": 2,
            "canonical_smiles": "CC(=O)Nc1ccncc1",
            "document_year": 2022,
            "standard_value": "50",
        },
        common
        | {
            "activity_id": 3,
            "canonical_smiles": "CC(=O)Nc1ccsc1",
            "document_year": 2022,
            "standard_value": "20000",
        },
        common
        | {"activity_id": 4, "canonical_smiles": "C" * 80, "document_year": 2022, "standard_value": "30"},
        common
        | {
            "activity_id": 5,
            "canonical_smiles": "CC(=O)NC1CCOC1",
            "document_year": 2022,
            "standard_value": "50",
        },
        common
        | {
            "activity_id": 6,
            "canonical_smiles": "CC(=O)NC1CCOC1",
            "document_year": 2022,
            "standard_value": "50000",
        },
    ]
    entries.append(
        common
        | {
            "activity_id": 7,
            "canonical_smiles": "CC(=O)Nc1ccccc1",
            "document_year": 2022,
            "standard_value": "50000",
        }
    )
    source = tmp_path / "raw.json"
    write_json(
        source,
        {
            "activities": entries,
            "assays": [
                {
                    "assay_chembl_id": "A",
                    "confidence_score": 9,
                    "relationship_type": "D",
                    "assay_organism": "Homo sapiens",
                    "target_chembl_id": "T",
                    "assay_type": "B",
                    "variant_sequence": None,
                }
            ],
        },
    )
    manifest = tmp_path / "manifest.json"
    write_json(
        manifest,
        {
            "scope": "synthetic_fixture",
            "sampling": "unit-test synthetic rows",
            "sources": [
                {
                    "name": "fixture",
                    "url": "https://example.invalid",
                    "path": str(source),
                    "sha256": digest(source),
                }
            ],
        },
    )
    actual = import_chembl.chemical_identity

    def check(smiles):
        assert smiles != "C" * 80, "Out-of-scope molecule reached expensive standardization"
        return actual(smiles)

    monkeypatch.setattr(import_chembl, "chemical_identity", check)
    result, report = import_chembl.curate(
        manifest, tmp_path / "curated", split_kind="time", queries_per_target=1
    )
    assert len(result["queries"]) == 1
    training = [row for row in result["compounds"] if row["split"] == "train"]
    assert len(training) == 1 and training[0]["label"] == "active"
    assert training[0]["supporting_activity_ids"] == [1]
    assert report["excluded_counts"]["family_seen_in_development"] == 1
    assert result["scope"] == "synthetic_fixture"
    assert report["excluded_counts"]["declared_heavy_atom_scope"] == 1
    assert report["excluded_counts"]["conflicting_family_labels_or_endpoints"] == 2
    assert all("date" not in row and isinstance(row["publication_year"], int) for row in result["compounds"])


def test_launch_error_and_nonfinite_response_do_not_abort_panel(tmp_path):
    req = request({}, tmp_path / "missing")
    method = {"kind": "external", "command": ["/nonexistent/benchmark/worker", "{request}", "{response}"]}
    assert bounded_method(method, req, tmp_path / "missing")["status"] == "launch_error"
    code = 'import sys; open(sys.argv[2],\'w\').write(\'{"candidates":[],"evaluations":1,"evaluation_unit":"test","bad":NaN}\')'
    method["command"] = [sys.executable, "-c", code, "{request}", "{response}"]
    assert (
        bounded_method(method, request({}, tmp_path / "nan"), tmp_path / "nan")["status"]
        == "invalid_response"
    )


def test_fixed_public_bound_pose_is_hashed_and_only_query_pose_forwarded(panel, tmp_path):
    data, path = panel
    pose = tmp_path / "query.sdf"
    pose.write_text("hash-only fixture; chemistry scorer validates the actual SDF separately")
    data["queries"][0]["reference_sdf"] = {"artifact": str(pose), "sha256": digest(pose)}
    write_json(path, data)
    validated, _ = load_dataset(path)
    req = public_request(
        validated,
        validated["queries"][0],
        {"id": "m", "kind": "external", "resolved_assets": {}},
        42,
        {},
        tmp_path,
    )
    assert req["query"]["reference_sdf_sha256"] == digest(pose)
    assert str(pose) == req["query"]["reference_sdf"]
    assert data["compounds"][1]["smiles"] not in json.dumps(req)
    assert "identities" not in json.dumps(validated)


def test_method_work_caps_differ_but_wall_and_output_caps_cannot(panel, assets, tmp_path):
    data, _ = panel
    base = {
        "schema_version": 1,
        "seeds": [42],
        "budget": {"wall_seconds": 10, "max_candidates": 3, "max_evaluations": 1000},
        "methods": [
            {
                "id": "m",
                "kind": "nearest_catalog",
                "budget_overrides": {"max_evaluations": 17},
                "assets": {"catalog": {"artifact": assets["catalog"], "sha256": digest(assets["catalog"])}},
            }
        ],
    }
    path = tmp_path / "caps.json"
    write_json(path, base)
    protocol = load_protocol(path, data)
    req = public_request(data, data["queries"][0], protocol["methods"][0], 42, protocol["budget"], tmp_path)
    assert req["budget"] == {"wall_seconds": 10, "max_candidates": 3, "max_evaluations": 17}
    base["methods"][0]["budget_overrides"]["wall_seconds"] = 20
    write_json(path, base)
    with pytest.raises(ValueError, match="Only max_evaluations"):
        load_protocol(path, data)
    base["methods"][0]["budget_overrides"] = {"max_evaluations": 1001}
    write_json(path, base)
    with pytest.raises(ValueError, match="exceeds common"):
        load_protocol(path, data)


def test_bad_generated_identity_does_not_erase_other_routes(assets, tmp_path, monkeypatch):
    from benchmarks.research import adapters

    real = adapters.chemical_identity

    def reject(smiles):
        if smiles == "CCC(=O)Nc1ccncc1":
            raise ValueError("injected bounded identity failure")
        return real(smiles)

    monkeypatch.setattr(adapters, "chemical_identity", reject)
    result = adapters.reaction_enumeration(request(assets, tmp_path))
    assert [r["smiles"] for r in result["candidates"]] == ["CC(=O)Nc1ccncc1"]
    assert result["evaluations"] == 2 and len(result["metadata"]["failures"]) == 1


def test_directory_asset_hashes_and_nested_model_dependencies_pinned(panel, assets, tmp_path):
    from benchmarks.research.build_protocol import artifact

    model = tmp_path / "model"
    model.mkdir()
    (model / "weights.safetensors").write_bytes(b"fixture weights")
    prompt = tmp_path / "prompt.json"
    write_json(prompt, {"instruction": "fixture", "input": "SMILES_STRING"})
    config = tmp_path / "synthesis.json"
    write_json(
        config,
        {
            "schema_version": 1,
            "catalog": assets["catalog"],
            "rules": assets["rules"],
            "provider": {"kind": "synllama", "model_path": str(model), "prompt_path": str(prompt)},
        },
    )
    scaffold = tmp_path / "scaffolds.smi"
    scaffold.write_text("c1ccncc1\n")
    declared = assets | {
        "synthesis_config": str(config),
        "scaffold_db": str(scaffold),
        "pipeline_provider_model_path": str(model),
        "pipeline_provider_prompt_path": str(prompt),
    }
    protocol = {
        "schema_version": 1,
        "seeds": [42],
        "budget": {"wall_seconds": 10, "max_candidates": 3, "max_evaluations": 1000},
        "methods": [
            {
                "id": "pipeline",
                "kind": "pipeline",
                "assets": {name: artifact(path, tmp_path) for name, path in declared.items()},
            }
        ],
    }
    path = tmp_path / "pinned.json"
    write_json(path, protocol)
    load_protocol(path, panel[0])
    (model / "weights.safetensors").write_bytes(b"changed model")
    with pytest.raises(ValueError, match="directory contents"):
        load_protocol(path, panel[0])
    del protocol["methods"][0]["assets"]["pipeline_provider_model_path"]
    write_json(path, protocol)
    with pytest.raises(ValueError, match="nested asset"):
        load_protocol(path, panel[0])
