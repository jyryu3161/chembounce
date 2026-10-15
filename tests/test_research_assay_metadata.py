import copy
import json

import pytest

from benchmarks.research.assay_metadata import assess_direct_human_assay, classify_assay, load_assay_metadata
from benchmarks.research.protocol import digest


@pytest.fixture
def assay_files(tmp_path):
    record = {
        "assay_chembl_id": "CHEMBL_ASSAY",
        "confidence_score": 9,
        "relationship_type": "D",
        "assay_organism": "Homo sapiens",
        "assay_tax_id": 9606,
        "target_chembl_id": "CHEMBL_TARGET",
        "assay_type": "B",
        "variant_sequence": None,
    }
    activity = {
        "assay_chembl_id": "CHEMBL_ASSAY",
        "target_chembl_id": "CHEMBL_TARGET",
        "target_organism": "Homo sapiens",
        "assay_type": "B",
    }
    assay = tmp_path / "assays.json"
    assay.write_text(json.dumps({"assays": [record]}))
    activities = tmp_path / "activities.json"
    activities.write_text(json.dumps({"activities": [activity]}))
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "sources": [
                    {
                        "name": "assays",
                        "url": "https://example.org/assays",
                        "path": str(assay),
                        "sha256": digest(assay),
                    },
                    {
                        "name": "activities",
                        "url": "https://example.org/activities",
                        "path": str(activities),
                        "sha256": digest(activities),
                    },
                ]
            }
        )
    )
    return manifest, assay, activity


def test_verified_assay_join_retains_source_provenance(assay_files):
    manifest, assay, activity = assay_files
    mapping, provenance = load_assay_metadata(manifest)
    accepted, reason, metadata = assess_direct_human_assay(activity, mapping)
    assert accepted and reason is None
    assert metadata["source_id"] == "assays"
    assert metadata["source_sha256"] == digest(assay)
    assert metadata["record"]["assay_organism"] == "Homo sapiens"
    assert provenance["assay_count"] == provenance["activity_assay_count"] == 1
    assert provenance["missing_activity_assay_ids"] == []
    assert classify_assay(activity, mapping) is None


@pytest.mark.parametrize(
    "change,reason",
    [
        ({"confidence_score": 8}, "assay_not_confidence9"),
        ({"confidence_score": None}, "assay_not_confidence9"),
        ({"relationship_type": "H"}, "assay_not_direct_target"),
        ({"assay_organism": "Mus musculus"}, "assay_not_explicitly_human"),
        ({"assay_organism": None}, "assay_not_explicitly_human"),
        ({"assay_tax_id": 10090}, "assay_taxonomy_conflicts_with_human"),
        ({"target_chembl_id": "OTHER_TARGET"}, "assay_target_mismatch"),
        ({"assay_type": "F"}, "assay_not_binding"),
        ({"variant_sequence": {"mutation": "T790M"}}, "assay_variant_sequence"),
    ],
)
def test_human_target_annotation_cannot_override_assay_metadata(assay_files, change, reason):
    manifest, _, activity = assay_files
    mapping, _ = load_assay_metadata(manifest)
    mapping[activity["assay_chembl_id"]]["record"].update(change)
    assert activity["target_organism"] == "Homo sapiens"
    accepted, actual, metadata = assess_direct_human_assay(activity, mapping)
    assert not accepted and actual == reason
    assert metadata["source_sha256"]


def test_missing_assay_and_missing_target_join_fail_closed(assay_files):
    manifest, _, activity = assay_files
    mapping, _ = load_assay_metadata(manifest)
    assert classify_assay({**activity, "assay_chembl_id": "UNKNOWN"}, mapping) == "missing_assay_metadata"
    assert classify_assay({**activity, "assay_chembl_id": None}, mapping) == "missing_assay_identifier"
    assert classify_assay({**activity, "target_chembl_id": None}, mapping) == "assay_target_mismatch"
    mapping[activity["assay_chembl_id"]]["record"].pop("variant_sequence")
    assert classify_assay(activity, mapping) == "missing_assay_variant_annotation"


def test_tampered_or_absent_metadata_inventory_fails_before_curation(assay_files):
    manifest, assay, _ = assay_files
    original = assay.read_text()
    assay.write_text("{}")
    with pytest.raises(ValueError, match="hash mismatch"):
        load_assay_metadata(manifest)
    assay.write_text(original)
    data = json.loads(manifest.read_text())
    data["sources"] = data["sources"][1:]
    manifest.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="No retained assay metadata"):
        load_assay_metadata(manifest)


def test_duplicate_conflicting_assay_records_are_rejected(assay_files):
    manifest, assay, _ = assay_files
    data = json.loads(assay.read_text())
    duplicate = copy.deepcopy(data["assays"][0])
    duplicate["confidence_score"] = 8
    data["assays"].append(duplicate)
    assay.write_text(json.dumps(data))
    metadata = json.loads(manifest.read_text())
    metadata["sources"][0]["sha256"] = digest(assay)
    manifest.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="Conflicting retained"):
        load_assay_metadata(manifest)


def test_portable_relative_manifest_and_explicit_coverage_gaps(assay_files):
    manifest, assay, activity = assay_files
    data = json.loads(manifest.read_text())
    data["sources"][0]["path"] = assay.name
    activities = manifest.parent / "activities.json"
    activities.write_text(json.dumps({"activities": [activity, {**activity, "assay_chembl_id": "MISSING"}]}))
    data["sources"][1]["sha256"] = digest(activities)
    manifest.write_text(json.dumps(data))
    _, provenance = load_assay_metadata(manifest)
    assert provenance["missing_activity_assay_ids"] == ["MISSING"]
