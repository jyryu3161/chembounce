"""Join independently retained assay annotations before claiming human-target labels.

An activity row's target annotation is insufficient: ChEMBL can map homologous
assays to that target. Direct human assignment requires the separate, hashed
assay record and explicit confidence/relationship/organism checks.
"""

from collections import Counter
import copy
import json
from pathlib import Path

from benchmarks.research.protocol import digest


DIRECT_HUMAN_POLICY = {
    "confidence_score": 9,
    "relationship_type": "D",
    "assay_organism": "Homo sapiens",
    "assay_type": "B",
    "target_match": "exact_activity_target_chembl_id",
    "variant_sequence": "must_be_null",
    "missing_metadata": "reject",
}


def _source_path(manifest, source):
    declared = source.get("path")
    expected = source.get("sha256")
    if not isinstance(declared, str) or not declared.strip():
        raise ValueError("Assay source manifests require a path for every retained source")
    if not isinstance(expected, str) or len(expected) != 64:
        raise ValueError("Assay source manifests require a SHA-256 for every retained source")
    raw = Path(declared)
    for path in (raw, manifest.parent / raw, manifest.parent / raw.name):
        if path.is_file():
            if digest(path) != expected:
                raise ValueError(f"Assay metadata source hash mismatch: {path}")
            return path.resolve()
    raise ValueError(f"Retained assay source is missing: {declared}")


def load_assay_metadata(manifest_path):
    """Return ``(assay_id -> sourced record, provenance)`` from verified artifacts.

    Both assay-list API responses and individual assay JSON records are accepted.
    Every listed source is hash-checked, including activity sources used to audit
    assay-ID coverage. Missing individual joins remain explicit rejections;
    absent assay metadata altogether is an early configuration error.
    """
    manifest = Path(manifest_path).resolve()
    data = json.loads(manifest.read_text())
    sources = data.get("sources") if isinstance(data, dict) else None
    if not isinstance(sources, list) or not sources:
        raise ValueError("Assay metadata manifest requires a nonempty sources list")
    assays, inventory, requested, source_ids = {}, [], set(), set()
    for index, source in enumerate(sources):
        if not isinstance(source, dict):
            raise ValueError("Source metadata must be an object")
        source_id = source.get("name", f"source-{index}")
        if not isinstance(source_id, str) or not source_id or source_id in source_ids:
            raise ValueError("Source metadata requires unique nonempty source IDs")
        source_ids.add(source_id)
        if not isinstance(source.get("url"), str) or not source["url"].strip():
            raise ValueError("Assay source provenance requires a source URL")
        path = _source_path(manifest, source)
        content = json.loads(path.read_text())
        if not isinstance(content, dict):
            continue
        if isinstance(content.get("activities"), list):
            requested.update(
                record["assay_chembl_id"]
                for record in content["activities"]
                if isinstance(record, dict) and isinstance(record.get("assay_chembl_id"), str)
            )
        if "assays" in content:
            records = content["assays"]
            if not isinstance(records, list):
                raise ValueError("Retained assays field must be a list")
        elif "assay_chembl_id" in content and "confidence_score" in content:
            records = [content]
        else:
            continue
        provenance = {
            "source_id": source_id,
            "source_sha256": source["sha256"],
            "source_url": source["url"],
            "source_artifact": str(path),
        }
        inventory.append({**provenance, "records": len(records)})
        for record in records:
            if (
                not isinstance(record, dict)
                or not isinstance(record.get("assay_chembl_id"), str)
                or not record["assay_chembl_id"]
            ):
                raise ValueError("Assay records require a nonempty assay_chembl_id")
            ident = record["assay_chembl_id"]
            if ident in assays and assays[ident]["record"] != record:
                raise ValueError(f"Conflicting retained assay metadata for {ident}")
            assays.setdefault(ident, {"record": copy.deepcopy(record), **provenance})
    if not assays:
        raise ValueError("No retained assay metadata; cannot establish direct human target assignment")
    return assays, {
        "manifest_path": str(manifest),
        "manifest_sha256": digest(manifest),
        "assay_sources": inventory,
        "assay_count": len(assays),
        "activity_assay_count": len(requested),
        "missing_activity_assay_ids": sorted(requested - set(assays)),
        "confidence_counts": dict(
            Counter(str(value["record"].get("confidence_score")) for value in assays.values())
        ),
        "policy": DIRECT_HUMAN_POLICY.copy(),
    }


def assess_direct_human_assay(activity, by_assay_id):
    """Return ``(accepted, rejection_reason, sourced_metadata)``; missing fails closed."""
    if not isinstance(activity, dict):
        return False, "invalid_activity_record", None
    ident = activity.get("assay_chembl_id")
    if not isinstance(ident, str) or not ident:
        return False, "missing_assay_identifier", None
    entry = by_assay_id.get(ident)
    if entry is None:
        return False, "missing_assay_metadata", None
    assay = entry["record"]
    metadata = {**copy.deepcopy(entry), "assignment_policy": DIRECT_HUMAN_POLICY.copy()}
    checks = (
        (assay.get("assay_chembl_id") == ident, "assay_identifier_mismatch"),
        (assay.get("confidence_score") == 9, "assay_not_confidence9"),
        (assay.get("relationship_type") == "D", "assay_not_direct_target"),
        (assay.get("assay_organism") == "Homo sapiens", "assay_not_explicitly_human"),
        (assay.get("assay_tax_id") in (None, 9606), "assay_taxonomy_conflicts_with_human"),
        (
            isinstance(activity.get("target_chembl_id"), str)
            and bool(activity["target_chembl_id"])
            and assay.get("target_chembl_id") == activity["target_chembl_id"],
            "assay_target_mismatch",
        ),
        (assay.get("assay_type") == "B", "assay_not_binding"),
        ("variant_sequence" in assay, "missing_assay_variant_annotation"),
        (assay.get("variant_sequence") is None, "assay_variant_sequence"),
    )
    for passed, reason in checks:
        if not passed:
            return False, reason, metadata
    return True, None, metadata


def classify_assay(activity, by_assay_id):
    """Convenience rejection-reason API: ``None`` only for an accepted assay."""
    return assess_direct_human_assay(activity, by_assay_id)[1]
