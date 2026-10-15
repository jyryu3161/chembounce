"""Assay-local source/label isolation, score-blind selection and clustered summaries."""

import copy
import json

import pytest

from benchmarks.research import assay_matched as matched
from benchmarks.research.protocol import digest, write_json

STRUCTURES = [
    ("c1ccccc1", 100),
    ("c1ccncc1", 200),
    ("c1ccsc1", 300),
    ("C1CCCCC1", 20000),
    ("c1ccc2ccccc2c1", 30000),
]


def write_sources(directory, *, changes=None):
    directory.mkdir(parents=True, exist_ok=True)
    assays, activities = [], []
    for index, (assay_id, target_id) in enumerate((("A1", "T1"), ("A2", "T1"), ("A3", "T2"))):
        assays.append(
            {
                "assay_chembl_id": assay_id,
                "target_chembl_id": target_id,
                "confidence_score": 9,
                "relationship_type": "D",
                "assay_organism": "Homo sapiens",
                "assay_tax_id": 9606,
                "assay_type": "B",
                "variant_sequence": None,
                "description": "Inhibition in preparation of unknown origin; context unresolved",
                "bao_format": "BAO:0000219",
                "bao_label": "cell-based format",
            }
        )
        for j, (smiles, value) in enumerate(STRUCTURES):
            activities.append(
                {
                    "activity_id": 100 * index + j + 1,
                    "assay_chembl_id": assay_id,
                    "target_chembl_id": target_id,
                    "target_organism": "Homo sapiens",
                    "assay_type": "B",
                    "standard_units": "nM",
                    "standard_type": "IC50",
                    "standard_flag": 1,
                    "standard_value": value,
                    "standard_relation": "=",
                    "canonical_smiles": smiles,
                    "document_year": 2020,
                    "molecule_chembl_id": f"M{j}",
                    "bao_label": "single protein format",
                }
            )
    if changes:
        changes(assays, activities)
    assay_path, activity_path = directory / "assays.json", directory / "activities.json"
    write_json(assay_path, {"assays": assays})
    write_json(activity_path, {"activities": activities})
    manifest = directory / "manifest.json"
    write_json(
        manifest,
        {
            "sampling": "Synthetic test fixture only",
            "license": {"name": "fixture"},
            "sources": [
                {
                    "name": name,
                    "path": path.name,
                    "url": f"https://example.invalid/{name}",
                    "sha256": digest(path),
                }
                for name, path in (("assays", assay_path), ("activities", activity_path))
            ],
        },
    )
    return manifest


@pytest.fixture
def pool(tmp_path):
    return matched.curate_sources(write_sources(tmp_path / "sources"))


def test_same_assay_source_labels_and_unresolved_metadata_are_preserved(pool):
    assert len(pool["records"]) == 15
    assert all(row["both_labels"] for row in pool["availability"]["contexts"])
    panel = matched.select_panel(pool)
    matched.validate_panel(panel, verify_sources=True)
    assert len(panel["queries"]) == 6
    for context in panel["contexts"]:
        raw = context["assay_metadata"]["record"]
        assert raw["assay_organism"] == "Homo sapiens" and "unknown origin" in raw["description"]
        assert raw["variant_sequence"] is None and raw["bao_label"] == "cell-based format"
    assert "not experimentally established wild-type" in " ".join(panel["interpretation"])
    for pair in matched.control_pairs(panel):
        assert pair["query_family"] != pair["control_family"]
        assert pair["query_scaffold_family"] != pair["control_scaffold_family"]
    assert all(set(q) == {"id", "smiles", "target_id", "context_id"} for q in matched.public_queries(panel))


def test_label_disagreement_between_assays_is_kept_but_within_assay_family_rejected(tmp_path):
    def changes(assays, activities):
        activities[5]["standard_value"] = 20000
        activities.append(
            {**activities[1], "activity_id": 99, "canonical_smiles": "c1cc[nH+]cc1", "standard_value": 20000}
        )

    pool = matched.curate_sources(write_sources(tmp_path / "sources", changes=changes))
    phenyl = [r for r in pool["records"] if r["family"] == "c1ccccc1" and r["target_id"] == "T1"]
    assert {(r["assay_id"], r["label"]) for r in phenyl} == {("A1", "active"), ("A2", "inactive")}
    assert not [r for r in pool["records"] if r["assay_id"] == "A1" and r["family"] == "c1ccncc1"]
    assert pool["availability"]["exclusion_counts"]["conflicting_labels_within_assay_family"] == 1
    assert pool["availability"]["repeated_target_families_with_different_context_labels"] == 1


def test_endpoint_isolation_censoring_and_direct_assignment_reused(tmp_path):
    def changes(assays, activities):
        activities.extend(
            [
                {**activities[0], "activity_id": 80, "standard_type": "Ki", "standard_value": 20000},
                {**activities[0], "activity_id": 81, "standard_value": 100, "standard_relation": ">"},
                {**activities[0], "activity_id": 82, "standard_value": 20000, "standard_relation": "<"},
                {**activities[0], "activity_id": 83, "target_chembl_id": "OTHER"},
                {**activities[0], "activity_id": 84, "standard_units": "pM"},
            ]
        )
        assays[2]["confidence_score"] = 8

    pool = matched.curate_sources(write_sources(tmp_path / "sources", changes=changes))
    a1 = [r for r in pool["records"] if r["assay_id"] == "A1" and r["family"] == "c1ccccc1"]
    assert {(r["endpoint"], r["label"]) for r in a1} == {("IC50", "active"), ("Ki", "inactive")}
    counts = pool["availability"]["exclusion_counts"]
    assert counts["ambiguous_or_intermediate_activity"] == 2
    assert counts["assay_target_mismatch"] == 1 and counts["assay_not_confidence9"] == 5
    assert counts["unsupported_endpoint_or_units"] == 1


def test_source_tampering_fails_and_duplicate_payload_conflict_fails(tmp_path):
    manifest = write_sources(tmp_path / "sources")
    (manifest.parent / "assays.json").write_text("{}")
    with pytest.raises(ValueError, match="hash mismatch"):
        matched.curate_sources(manifest)

    def changes(assays, activities):
        activities.append({**activities[0], "standard_value": 50000})

    with pytest.raises(ValueError, match="Conflicting source payload"):
        matched.curate_sources(write_sources(tmp_path / "second", changes=changes))


def test_selection_is_order_independent_and_meets_scaffold_minimum_after_cap(pool):
    policy = matched.Selection(max_controls_per_label=2, max_queries_per_context=1, max_contexts_per_target=1)
    first = matched.select_panel(pool, policy)
    shuffled = copy.deepcopy(pool)
    shuffled["records"].reverse()
    shuffled["contexts"].reverse()
    second = matched.select_panel(shuffled, policy)
    assert first["queries"] == second["queries"]
    assert len(first["queries"]) == 2
    for query in first["queries"]:
        pairs = [pair for pair in matched.control_pairs(first) if pair["query_id"] == query["id"]]
        assert len(pairs) == 4
        assert len({p["control_scaffold_family"] for p in pairs if p["label"] == "active"}) == 2
    assert first["selection_audit"]["exclusion_counts"]["target_context_cap"] == 1


def test_insufficient_targets_are_explicit_and_not_filled_by_lowering_minima(pool):
    panel = matched.select_panel(
        pool,
        matched.Selection(
            min_active_controls=5,
            min_inactive_controls=5,
            min_active_scaffolds=1,
            max_contexts_per_target=2,
            max_queries_per_context=1,
        ),
    )
    assert panel["queries"] == []
    assert panel["selection_audit"]["targets_available"] == ["T1", "T2"]
    assert panel["selection_audit"]["targets_selected"] == []
    assert panel["selection_audit"]["exclusion_counts"]["no_eligible_assay_context"] == 2
    assert matched.aggregate_scores(panel, [])["target_macro"]["roc_auc"] is None


def test_cross_context_reference_or_tampered_label_is_rejected(pool):
    panel = matched.select_panel(pool)
    query = panel["queries"][0]
    other = next(
        r for r in panel["records"] if r["context_id"] != query["context_id"] and r["label"] == "inactive"
    )
    query["control_record_ids"][0] = other["id"]
    with pytest.raises(ValueError, match="Controls cross context"):
        matched.validate_panel(panel)
    panel = matched.select_panel(pool)
    panel["records"][0]["label"] = "inactive" if panel["records"][0]["label"] == "active" else "active"
    with pytest.raises(ValueError, match="does not support"):
        matched.validate_panel(panel, verify_sources=True)


def test_activity_payload_cannot_be_changed_while_attaching_valid_source_hash(pool):
    panel = matched.select_panel(pool)
    panel["records"][0]["measurements"][0]["activity"]["standard_value"] *= 1.01
    with pytest.raises(ValueError, match="differs from retained source"):
        matched.validate_panel(panel, verify_sources=True)


def test_target_assay_hierarchy_and_missing_pairs_preserve_denominators(pool):
    panel = matched.select_panel(pool, matched.Selection(max_queries_per_context=1))
    pairs = list(matched.control_pairs(panel))
    # Synthetic unit-test scores: target macro must be 0.5, not assay-pooled 2/3.
    rows = [
        {
            "pair_id": p["pair_id"],
            "scores": {"score": int((p["label"] == "active") == (p["target_id"] == "T1"))},
        }
        for p in pairs
    ]
    result = matched.aggregate_scores(panel, rows, resamples=100)
    assert len(result["assays"]) == 3 and len(result["targets"]) == 2
    assert result["target_macro"]["roc_auc"] == 0.5
    assert result["success_conditional_target_bootstrap"]["roc_auc"]["target_groups"] == 2
    missing = matched.aggregate_scores(panel, [], resamples=100)
    assert (
        missing["pair_count"]
        == missing["missing_or_failed_pairs"]
        == missing["unreturned_pairs"]
        == len(pairs)
    )
    assert missing["target_macro"]["roc_auc"] is None
    assert missing["target_macro"]["roc_auc_bounds"] == [0, 1]
    rows[0]["status"] = "failed"
    assert matched.aggregate_scores(panel, rows, resamples=100)["missing_or_failed_pairs"] == 1
    with pytest.raises(ValueError, match="Unknown or duplicate"):
        matched.aggregate_scores(panel, rows + rows[:1])


def test_cli_pins_pool_and_emits_evaluator_pairs_separately_from_public_queries(pool, tmp_path):
    path = tmp_path / "pool.json"
    write_json(path, pool)
    output = tmp_path / "selected"
    with pytest.raises(SystemExit):
        matched.main(["select", "--pool", str(path), "--pool-sha256", "0" * 64, "--output", str(output)])
    matched.main(["select", "--pool", str(path), "--pool-sha256", digest(path), "--output", str(output)])
    assert (output / "control_pairs.json").is_file()
    public = json.loads((output / "public_queries.json").read_text())
    assert public and all("label" not in row and "control_record_ids" not in row for row in public)
    assert json.loads((output / "artifacts.json").read_text())["panel_sha256"] == digest(
        output / "panel.json"
    )


def test_two_endpoints_in_one_assay_do_not_double_its_target_weight(tmp_path):
    def changes(assays, activities):
        activities.extend(
            {**row, "activity_id": row["activity_id"] + 1000, "standard_type": "Ki"}
            for row in list(activities)
            if row["assay_chembl_id"] == "A1"
        )

    pool = matched.curate_sources(write_sources(tmp_path / "sources", changes=changes))
    panel = matched.select_panel(pool, matched.Selection(max_queries_per_context=1))
    pairs = list(matched.control_pairs(panel))
    scores = []
    for pair in pairs:
        perfect = pair["target_id"] == "T1" and not (pair["assay_id"] == "A1" and pair["endpoint"] == "Ki")
        scores.append(
            {"pair_id": pair["pair_id"], "scores": {"score": int((pair["label"] == "active") == perfect)}}
        )
    report = matched.aggregate_scores(panel, scores, resamples=10)
    assert len(report["contexts"]) == 4 and len(report["assays"]) == 3
    assert next(r for r in report["assays"] if r["assay_id"] == "A1")["metrics"]["roc_auc"] == 0.5
    assert next(r for r in report["targets"] if r["target_id"] == "T1")["metrics"]["roc_auc"] == 0.75
    assert report["target_macro"]["roc_auc"] == 0.375


def test_changing_exact_protonation_state_cannot_bypass_source_measurement(pool):
    panel = matched.select_panel(pool)
    row = next(r for r in panel["records"] if r["smiles"] == "c1ccncc1")
    row["smiles"] = "c1cc[nH+]cc1"
    with pytest.raises(ValueError, match="Exact query/control state differs"):
        matched.validate_panel(panel, verify_sources=True)


def test_audit_cli_records_hashes_and_all_context_counts(tmp_path):
    manifest = write_sources(tmp_path / "sources")
    output = tmp_path / "audit"
    matched.main(["audit", "--manifest", str(manifest), "--output", str(output)])
    audit = json.loads((output / "availability.json").read_text())
    assert audit["raw_unique_measurements"] == 15
    assert len(audit["contexts"]) == 3
    hashes = json.loads((output / "artifacts.json").read_text())
    assert hashes["pool_sha256"] == digest(output / "pool.json")
    with pytest.raises(SystemExit):
        matched.main(["audit", "--manifest", str(manifest), "--output", str(output)])
