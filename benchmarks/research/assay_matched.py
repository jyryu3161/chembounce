"""Construct score-blind, same-assay similarity controls from pinned observations.

This schema deliberately does not use the target-pooled generation benchmark.
An assay context is (assay ID, exact target, endpoint, units); conflicting labels
are resolved only within that context and chemical family. No chemical scores
are calculated here.
"""

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

from rdkit import Chem, rdBase

from benchmarks.research.assay_metadata import (
    _source_path,
    assess_direct_human_assay,
    load_assay_metadata,
)
from benchmarks.research.import_chembl import activity_label
from benchmarks.research.protocol import digest, write_json
from chembounce_synthesis.identity import chemical_identity, IDENTITY_POLICY


SCHEMA = "assay-matched-controls-v1"
CONTEXT_FIELDS = ("assay_chembl_id", "target_chembl_id", "standard_type", "standard_units")
INTERPRETATION = [
    "Same assay ID narrows the experimental context; it does not prove identical conditions for every measurement.",
    "Structured human assignment is retained; assay prose is preserved without using regex matches as biological truth.",
    "Null variant_sequence means no variant annotation, not experimentally established wild-type sequence.",
    "BAO labels may be incomplete or inconsistent with prose; they are audit metadata, not condition equivalence.",
    "Repeated chemical families across assays and targets remain correlated observations, not independent replicates.",
    "This retrospective same-assay control panel does not replace a generation benchmark train/test partition.",
]
AGGREGATION = {
    "hierarchy": ["query", "context", "assay", "target"],
    "weights": "Equal child weights at each level; contexts within an assay cannot give that assay extra target weight",
    "uncertainty_unit": "Resample target groups only, never query-control pairs, assays or conformer seeds as independent targets",
    "missing_scores": "Keep every selected pair; missing/failed scores enter pessimistic and optimistic ranking bounds",
    "conditional_estimates": "Report success-conditional estimates and unestimable groups explicitly; do not replace missing estimates with zero",
    "dependence": "Target clustering contains within-target assay/family repeats; shared chemistry and related targets can still correlate across clusters",
    "interpretation": "Descriptive uncertainty for sampled targets, not proof of broad generalization or calibrated activity probability",
}


def stable_hash(*values):
    return hashlib.sha256(json.dumps(values, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def context_key(row):
    values = tuple(row.get(field) for field in CONTEXT_FIELDS)
    return values if all(isinstance(v, str) and v.strip() for v in values) else None


def context_record(key, metadata):
    assay, target, endpoint, units = key
    return {
        "id": "context-" + stable_hash(*key)[:24],
        "assay_id": assay,
        "target_id": target,
        "endpoint": endpoint,
        "units": units,
        "assay_metadata": metadata,
        "annotation_interpretation": INTERPRETATION[:4],
    }


def eligible_controls(query, records):
    """Both labels must differ from the query's chemical family AND scaffold."""
    return [
        row
        for row in records
        if row["context_id"] == query["context_id"]
        and row["family"] != query["family"]
        and row["scaffold_family"] != query["scaffold_family"]
    ]


@dataclass(frozen=True)
class Selection:
    seed: int = 42
    min_active_controls: int = 2
    min_inactive_controls: int = 2
    min_active_scaffolds: int = 2
    max_controls_per_label: int = 10
    max_queries_per_context: int = 2
    max_contexts_per_assay: int = 2
    max_assays_per_target: int = 4
    max_contexts_per_target: int = 4

    def __post_init__(self):
        for key, value in asdict(self).items():
            if isinstance(value, bool) or not isinstance(value, int) or (key != "seed" and value <= 0):
                raise ValueError(f"{key} must be {'an integer' if key == 'seed' else 'a positive integer'}")
        if self.max_controls_per_label < max(
            self.min_active_controls, self.min_inactive_controls, self.min_active_scaffolds
        ):
            raise ValueError("Control cap cannot be smaller than minimum controls or active scaffolds")


def query_availability(query, records):
    controls = eligible_controls(query, records)
    return {
        "query_record_id": query["id"],
        "active_controls": sum(r["label"] == "active" for r in controls),
        "inactive_controls": sum(r["label"] == "inactive" for r in controls),
        "active_scaffolds": len({r["scaffold_family"] for r in controls if r["label"] == "active"}),
    }


def inadequacies(available, policy):
    return [
        name
        for name, value, minimum in (
            ("insufficient_active_controls", available["active_controls"], policy.min_active_controls),
            ("insufficient_inactive_controls", available["inactive_controls"], policy.min_inactive_controls),
            (
                "insufficient_active_reference_scaffolds",
                available["active_scaffolds"],
                policy.min_active_scaffolds,
            ),
        )
        if value < minimum
    ]


def curate_sources(manifest_path, *, active_nm=1000.0, inactive_nm=10000.0, max_heavy_atoms=60):
    """Verify every retained source, then collapse labels inside assay contexts."""
    if not 0 < active_nm < inactive_nm:
        raise ValueError("Require 0 < active_nm < inactive_nm")
    if isinstance(max_heavy_atoms, bool) or not isinstance(max_heavy_atoms, int) or max_heavy_atoms < 1:
        raise ValueError("max_heavy_atoms must be a positive integer")
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    metadata, provenance = load_assay_metadata(manifest_path)
    sources, seen, excluded, contexts = [], {}, [], {}
    raw_context_counts = Counter()
    for index, source in enumerate(manifest["sources"], 1):
        path = _source_path(manifest_path, source)
        source_id = source["name"]
        sources.append(
            {"id": source_id, "url": source["url"], "artifact": str(path), "sha256": source["sha256"]}
        )
        payload = json.loads(path.read_text())
        activities = payload.get("activities", []) if isinstance(payload, dict) else []
        if not isinstance(activities, list):
            continue
        print(
            f"Assay-matched source {index}/{len(manifest['sources'])}: {source_id}, {len(activities)} activities",
            file=sys.stderr,
            flush=True,
        )
        for row in activities:
            if not isinstance(row, dict) or not isinstance(row.get("activity_id"), (int, str)):
                excluded.append({"source_id": source_id, "reason": "missing_activity_identifier"})
                continue
            activity_id = str(row["activity_id"])
            if activity_id in seen:
                if seen[activity_id]["row"] != row:
                    # A duplicate ID with different contents is provenance
                    # corruption, not an assay-label disagreement to average.
                    raise ValueError(f"Conflicting source payload for activity ID {activity_id}")
                seen[activity_id]["source_ids"].add(source_id)
            else:
                seen[activity_id] = {"row": row, "source_ids": {source_id}}
    grouped = defaultdict(list)
    for index, (activity_id, item) in enumerate(sorted(seen.items()), 1):
        if index % 500 == 0:
            print(f"Assay-matched chemical curation {index}/{len(seen)}", file=sys.stderr, flush=True)
        row = item["row"]
        key = context_key(row)
        if key is not None:
            contexts.setdefault(key, context_record(key, metadata.get(key[0])))
            raw_context_counts[key] += 1
        label, reason = activity_label(row, active_nm, inactive_nm)
        if not reason:
            _, reason, _ = assess_direct_human_assay(row, metadata)
        if not reason and key is None:
            reason = "missing_assay_context_identifier"
        smiles = row.get("canonical_smiles")
        if not reason:
            mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
            if mol is None:
                reason = "invalid_smiles"
            elif len(Chem.GetMolFrags(mol)) != 1:
                reason = "disconnected_chemical_scope"
            elif mol.GetNumHeavyAtoms() > max_heavy_atoms:
                reason = "declared_heavy_atom_scope"
            elif not isinstance(row.get("document_year"), int) or isinstance(row.get("document_year"), bool):
                reason = "publication_year_missing"
        detail = {}
        if not reason:
            try:
                identity = chemical_identity(smiles)
                if not identity.has_scaffold:
                    reason = "acyclic_not_scaffold_control"
            except (ValueError, RuntimeError) as exc:
                reason, detail = "chemical_identity_failed", {"error": str(exc)}
        if reason:
            excluded.append(
                {
                    "activity_id": activity_id,
                    "source_ids": sorted(item["source_ids"]),
                    "context_id": contexts[key]["id"] if key else None,
                    "reason": reason,
                    **detail,
                }
            )
            continue
        grouped[(key, identity.family_smiles)].append({**item, "label": label, "identity": identity})
    records = []
    for (key, family), items in sorted(grouped.items()):
        if len({item["label"] for item in items}) != 1:
            excluded.append(
                {
                    "context_id": contexts[key]["id"],
                    "family": family,
                    "activity_ids": sorted(str(item["row"]["activity_id"]) for item in items),
                    "reason": "conflicting_labels_within_assay_family",
                }
            )
            continue
        chosen = min(
            items,
            key=lambda item: (
                item["row"]["document_year"],
                str(item["row"]["activity_id"]),
                item["identity"].exact_smiles,
            ),
        )
        identity = chosen["identity"]
        records.append(
            {
                "id": "record-" + stable_hash(key, family)[:24],
                "context_id": contexts[key]["id"],
                "assay_id": key[0],
                "target_id": key[1],
                "endpoint": key[2],
                "units": key[3],
                "smiles": identity.exact_smiles,
                "family": family,
                "scaffold_family": identity.scaffold_family_smiles,
                "generic_skeleton": identity.ring_skeleton_smiles,
                "label": chosen["label"],
                "source_ids": sorted(set().union(*(item["source_ids"] for item in items))),
                "measurements": [
                    {"activity": item["row"], "source_ids": sorted(item["source_ids"])}
                    for item in sorted(items, key=lambda item: str(item["row"]["activity_id"]))
                ],
            }
        )
    by_context = defaultdict(list)
    for row in records:
        by_context[row["context_id"]].append(row)
    availability = []
    for key, context in sorted(contexts.items()):
        local = by_context[context["id"]]
        available = [query_availability(row, local) for row in local if row["label"] == "active"]
        counts = Counter(row["label"] for row in local)
        availability.append(
            {
                "context_id": context["id"],
                "assay_id": key[0],
                "target_id": key[1],
                "endpoint": key[2],
                "units": key[3],
                "raw_measurements": raw_context_counts[key],
                "family_records": len(local),
                "active_families": counts["active"],
                "inactive_families": counts["inactive"],
                "active_scaffolds": len(
                    {row["scaffold_family"] for row in local if row["label"] == "active"}
                ),
                "both_labels": bool(counts["active"] and counts["inactive"]),
                "query_availability": available,
            }
        )
    per_target = {}
    for target in sorted({key[1] for key in contexts}):
        local = [row for row in availability if row["target_id"] == target]
        per_target[target] = {
            "contexts": len(local),
            "assays": len({row["assay_id"] for row in local}),
            "family_context_records": sum(row["family_records"] for row in local),
            "both_label_contexts": sum(row["both_labels"] for row in local),
            "both_label_assays": len({row["assay_id"] for row in local if row["both_labels"]}),
            "query_feasibility": {
                f"min_{minimum}_per_label_and_{scaffolds}_active_scaffolds": sum(
                    not inadequacies(
                        query,
                        Selection(
                            min_active_controls=minimum,
                            min_inactive_controls=minimum,
                            min_active_scaffolds=scaffolds,
                        ),
                    )
                    for row in local
                    for query in row["query_availability"]
                )
                for minimum, scaffolds in ((1, 1), (2, 2), (5, 1), (5, 2))
            },
        }
    family_contexts = defaultdict(list)
    for row in records:
        family_contexts[(row["target_id"], row["family"])].append(row)
    repeated = [items for items in family_contexts.values() if len(items) > 1]
    availability_report = {
        "scope": "Availability in this retained bounded sample, not all ChEMBL; no similarity scores computed",
        "raw_unique_measurements": len(seen),
        "curated_family_context_records": len(records),
        "per_target": per_target,
        "contexts": availability,
        "exclusion_counts": dict(Counter(row["reason"] for row in excluded)),
        "exclusion_count_unit": "One row per excluded measurement, except one row per conflicting context-family",
        "repeated_target_families_across_contexts": len(repeated),
        "repeated_target_families_with_different_context_labels": sum(
            len({r["label"] for r in items}) > 1 for items in repeated
        ),
    }
    return {
        "schema": SCHEMA,
        "stage": "curated_pool",
        "sources": sources,
        "source_manifest_sha256": digest(manifest_path),
        "source_sampling": manifest.get("sampling"),
        "source_license": manifest.get("license"),
        "assay_source_provenance": provenance,
        "curation": {
            "active_nm": active_nm,
            "inactive_nm": inactive_nm,
            "max_heavy_atoms": max_heavy_atoms,
            "chemical_identity_policy": IDENTITY_POLICY,
            "rdkit_version": rdBase.rdkitVersion,
            "conflict_scope": list(CONTEXT_FIELDS) + ["chemical_family"],
            "representative_state": "Earliest publication year, then activity ID; no score or potency selection",
        },
        "contexts": list(sorted(contexts.values(), key=lambda row: row["id"])),
        "records": records,
        "availability": availability_report,
        "exclusions": excluded,
        "interpretation": INTERPRETATION,
    }


def _ranked(rows, policy, namespace, key=lambda row: row["id"]):
    return sorted(rows, key=lambda row: (stable_hash(policy.seed, namespace, key(row)), key(row)))


def _select_controls(query, records, policy):
    pools = defaultdict(list)
    for row in eligible_controls(query, records):
        pools[row["label"]].append(row)
    active = _ranked(pools["active"], policy, [query["id"], "active_controls"])
    # Reserve distinct active scaffolds before filling the remaining cap. A
    # naive family sample could erase the diversity required for eligibility.
    selected, scaffolds = [], set()
    for row in active:
        if row["scaffold_family"] not in scaffolds:
            selected.append(row)
            scaffolds.add(row["scaffold_family"])
        if len(scaffolds) == policy.min_active_scaffolds:
            break
    selected_ids = {row["id"] for row in selected}
    selected.extend(row for row in active if row["id"] not in selected_ids)
    inactive = _ranked(pools["inactive"], policy, [query["id"], "inactive_controls"])
    return selected[: policy.max_controls_per_label] + inactive[: policy.max_controls_per_label]


def select_panel(pool, policy=None):
    """Select after availability review, before any chemical score is observed."""
    policy = policy or Selection()
    if pool.get("schema") != SCHEMA or pool.get("stage") != "curated_pool":
        raise ValueError("Expected an assay-matched curated pool")
    records = {row["id"]: row for row in pool["records"]}
    if len(records) != len(pool["records"]):
        raise ValueError("Duplicate curated record IDs")
    by_context = defaultdict(list)
    for row in records.values():
        by_context[row["context_id"]].append(row)
    eligible, exclusions, context_counts = {}, [], []
    for context in sorted(pool["contexts"], key=lambda row: row["id"]):
        local = by_context[context["id"]]
        options = []
        for query in sorted(local, key=lambda row: row["id"]):
            if query["label"] != "active":
                continue
            available = query_availability(query, local)
            reasons = inadequacies(available, policy)
            if reasons:
                exclusions.append(
                    {"level": "query", "context_id": context["id"], **available, "reasons": reasons}
                )
            else:
                options.append(query)
        context_counts.append(
            {
                "context_id": context["id"],
                "target_id": context["target_id"],
                "assay_id": context["assay_id"],
                "eligible_queries": len(options),
            }
        )
        if options:
            eligible[context["id"]] = (context, options)
        else:
            exclusions.append(
                {
                    "level": "context",
                    "context_id": context["id"],
                    "target_id": context["target_id"],
                    "assay_id": context["assay_id"],
                    "reason": "no_query_meets_reference_requirements",
                    "family_records": len(local),
                }
            )
    by_target = defaultdict(lambda: defaultdict(list))
    for context, options in eligible.values():
        by_target[context["target_id"]][context["assay_id"]].append((context, options))
    selected_contexts, queries, chosen_records = [], [], {}
    for target in sorted({context["target_id"] for context in pool["contexts"]}):
        target_contexts = 0
        assay_ids = sorted(by_target[target])
        ranked_assays = _ranked(assay_ids, policy, [target, "assays"], key=lambda value: value)
        if not ranked_assays:
            exclusions.append({"level": "target", "target_id": target, "reason": "no_eligible_assay_context"})
        for assay_id in ranked_assays[policy.max_assays_per_target :]:
            exclusions.append(
                {"level": "assay", "target_id": target, "assay_id": assay_id, "reason": "assay_cap"}
            )
        for assay_id in ranked_assays[: policy.max_assays_per_target]:
            contexts = _ranked(
                by_target[target][assay_id],
                policy,
                [target, assay_id, "contexts"],
                key=lambda value: value[0]["id"],
            )
            for context, _ in contexts[policy.max_contexts_per_assay :]:
                exclusions.append({"level": "context", "context_id": context["id"], "reason": "context_cap"})
            for context, options in contexts[: policy.max_contexts_per_assay]:
                if target_contexts >= policy.max_contexts_per_target:
                    exclusions.append(
                        {
                            "level": "context",
                            "context_id": context["id"],
                            "target_id": target,
                            "assay_id": assay_id,
                            "reason": "target_context_cap",
                        }
                    )
                    continue
                target_contexts += 1
                selected_contexts.append(context)
                seen_scaffolds, chosen_queries = set(), []
                for row in _ranked(options, policy, [context["id"], "queries"]):
                    reason = None
                    if row["scaffold_family"] in seen_scaffolds:
                        reason = "query_scaffold_already_selected"
                    elif len(chosen_queries) >= policy.max_queries_per_context:
                        reason = "query_cap"
                    if reason:
                        exclusions.append(
                            {
                                "level": "query",
                                "context_id": context["id"],
                                "query_record_id": row["id"],
                                "reason": reason,
                            }
                        )
                        continue
                    chosen_queries.append(row)
                    seen_scaffolds.add(row["scaffold_family"])
                    controls = _select_controls(row, by_context[context["id"]], policy)
                    query_id = "query-" + stable_hash(context["id"], row["id"])[:24]
                    queries.append(
                        {
                            "id": query_id,
                            "record_id": row["id"],
                            "context_id": context["id"],
                            "assay_id": context["assay_id"],
                            "target_id": target,
                            "available": query_availability(row, by_context[context["id"]]),
                            "control_record_ids": [control["id"] for control in controls],
                        }
                    )
                    chosen_records.update({value["id"]: value for value in [row, *controls]})
    panel = {
        "schema": SCHEMA,
        "stage": "selected_panel",
        "scope": "same_assay_similarity_controls",
        "sources": pool["sources"],
        "source_manifest_sha256": pool["source_manifest_sha256"],
        "source_sampling": pool["source_sampling"],
        "source_license": pool.get("source_license"),
        "curation": pool["curation"],
        "selection": asdict(policy),
        "selection_rule": "SHA-256 ranking of [seed, namespace, stable identifier], no similarity/3D/potency-score ordering; active scaffold reservation before control cap",
        "contexts": selected_contexts,
        "records": sorted(chosen_records.values(), key=lambda row: row["id"]),
        "queries": queries,
        "aggregation": AGGREGATION,
        "interpretation": INTERPRETATION,
        "selection_audit": {
            "contexts": context_counts,
            "exclusions": exclusions,
            "exclusion_counts": dict(
                Counter(
                    row.get("reason", reason)
                    for row in exclusions
                    for reason in row.get("reasons", [row.get("reason")])
                )
            ),
            "targets_available": sorted({c["target_id"] for c in pool["contexts"]}),
            "targets_selected": sorted({q["target_id"] for q in queries}),
        },
    }
    validate_panel(panel)
    return panel


def validate_panel(panel, *, verify_sources=False):
    """Validate same-context labels/identities and query/control separation.

    Optional source verification checks exact raw activity payloads and assay
    metadata as well as hashes; merely attaching a valid source hash is not
    sufficient evidence for a modified label or molecule.
    """
    if panel.get("schema") != SCHEMA or panel.get("stage") != "selected_panel":
        raise ValueError("Expected a selected assay-matched panel")
    policy = Selection(**panel["selection"])
    contexts = {c["id"]: c for c in panel["contexts"]}
    records = {r["id"]: r for r in panel["records"]}
    if len(contexts) != len(panel["contexts"]) or len(records) != len(panel["records"]):
        raise ValueError("Duplicate context or record ID")
    sources, observed_activities, observed_assays = {}, {}, {}
    for source in panel["sources"]:
        if source["id"] in sources:
            raise ValueError("Duplicate source ID")
        sources[source["id"]] = source
        if verify_sources:
            path = Path(source["artifact"])
            if digest(path) != source["sha256"]:
                raise ValueError(f"Source hash mismatch: {path}")
            payload = json.loads(path.read_text())
            if isinstance(payload, dict) and "assay_chembl_id" in payload and "confidence_score" in payload:
                observed_assays[(source["id"], payload["assay_chembl_id"])] = payload
            for field, lookup, ident in (
                ("activities", observed_activities, "activity_id"),
                ("assays", observed_assays, "assay_chembl_id"),
            ):
                values = payload.get(field, []) if isinstance(payload, dict) else []
                if isinstance(values, list):
                    for row in values:
                        lookup[(source["id"], str(row[ident]))] = row
    unique_families = set()
    for row in records.values():
        context = contexts.get(row["context_id"])
        if context is None or any(
            row[field] != context[field] for field in ("assay_id", "target_id", "endpoint", "units")
        ):
            raise ValueError("Record crosses assay context")
        identity = chemical_identity(row["smiles"])
        if not identity.has_scaffold or (row["family"], row["scaffold_family"], row["generic_skeleton"]) != (
            identity.family_smiles,
            identity.scaffold_family_smiles,
            identity.ring_skeleton_smiles,
        ):
            raise ValueError("Record chemical identity mismatch")
        if (row["context_id"], row["family"]) in unique_families:
            raise ValueError("Duplicate family within assay context")
        unique_families.add((row["context_id"], row["family"]))
        metadata = context.get("assay_metadata")
        if not metadata or metadata.get("source_id") not in sources:
            raise ValueError("Missing sourced assay metadata")
        if metadata["source_sha256"] != sources[metadata["source_id"]]["sha256"]:
            raise ValueError("Assay metadata source hash differs")
        if (
            verify_sources
            and observed_assays.get((metadata["source_id"], context["assay_id"])) != metadata["record"]
        ):
            raise ValueError("Assay metadata differs from retained source")
        if not row.get("measurements"):
            raise ValueError("Record has no supporting measurements")
        representatives = []
        for measurement in row["measurements"]:
            activity = measurement["activity"]
            if context_key(activity) != (
                context["assay_id"],
                context["target_id"],
                context["endpoint"],
                context["units"],
            ):
                raise ValueError("Supporting measurement crosses assay context")
            label, reason = activity_label(
                activity, panel["curation"]["active_nm"], panel["curation"]["inactive_nm"]
            )
            accepted, _, _ = assess_direct_human_assay(activity, {context["assay_id"]: metadata})
            if reason or label != row["label"] or not accepted:
                raise ValueError("Supporting measurement does not support the context label")
            observed_identity = chemical_identity(activity["canonical_smiles"])
            if observed_identity.family_smiles != row["family"]:
                raise ValueError("Supporting measurement belongs to another chemical family")
            year = activity.get("document_year")
            if isinstance(year, bool) or not isinstance(year, int):
                raise ValueError("Supporting measurement lacks a publication year")
            if (
                Chem.MolFromSmiles(observed_identity.exact_smiles).GetNumHeavyAtoms()
                > panel["curation"]["max_heavy_atoms"]
            ):
                raise ValueError("Supporting measurement exceeds chemical size scope")
            representatives.append((year, str(activity["activity_id"]), observed_identity.exact_smiles))
            if not measurement["source_ids"] or any(
                source_id not in sources for source_id in measurement["source_ids"]
            ):
                raise ValueError("Unknown measurement source")
            if verify_sources and any(
                observed_activities.get((source_id, str(activity["activity_id"]))) != activity
                for source_id in measurement["source_ids"]
            ):
                raise ValueError("Measurement differs from retained source")
        if row["smiles"] != min(representatives)[2]:
            raise ValueError("Exact query/control state differs from the declared representative measurement")
    query_ids = set()
    query_counts, query_scaffolds = Counter(), defaultdict(set)
    for query in panel["queries"]:
        if query["id"] in query_ids:
            raise ValueError("Duplicate query ID")
        query_ids.add(query["id"])
        row = records[query["record_id"]]
        if row["label"] != "active" or any(
            query[field] != row[field] for field in ("context_id", "target_id", "assay_id")
        ):
            raise ValueError("Query must be an active record in its own context")
        if query["id"] != "query-" + stable_hash(query["context_id"], row["id"])[:24]:
            raise ValueError("Query ID does not match its fixed context and record")
        query_counts[query["context_id"]] += 1
        if row["scaffold_family"] in query_scaffolds[query["context_id"]]:
            raise ValueError("Repeated query scaffold within context")
        query_scaffolds[query["context_id"]].add(row["scaffold_family"])
        control_ids = query["control_record_ids"]
        if len(control_ids) != len(set(control_ids)):
            raise ValueError("Duplicate control record for query")
        controls = [records[ident] for ident in control_ids]
        if len(eligible_controls(row, controls)) != len(controls):
            raise ValueError("Controls cross context or share the query family/scaffold")
        if inadequacies(query_availability(row, controls), policy):
            raise ValueError("Selected controls do not meet declared minima")
        if any(
            count > policy.max_controls_per_label for count in Counter(r["label"] for r in controls).values()
        ):
            raise ValueError("Selected controls exceed label cap")
    if any(n > policy.max_queries_per_context for n in query_counts.values()):
        raise ValueError("Selected queries exceed context cap")
    assay_contexts, target_assays = defaultdict(set), defaultdict(set)
    for context in contexts.values():
        assay_contexts[(context["target_id"], context["assay_id"])].add(context["id"])
        target_assays[context["target_id"]].add(context["assay_id"])
    if any(len(values) > policy.max_contexts_per_assay for values in assay_contexts.values()):
        raise ValueError("Selected contexts exceed assay cap")
    if any(len(values) > policy.max_assays_per_target for values in target_assays.values()):
        raise ValueError("Selected assays exceed target cap")
    if any(
        count > policy.max_contexts_per_target
        for count in Counter(c["target_id"] for c in contexts.values()).values()
    ):
        raise ValueError("Selected contexts exceed target cap")
    return panel


def control_pairs(panel):
    """Evaluator-only rows; never provide their labels/structures to a generator."""
    records = {row["id"]: row for row in panel["records"]}
    contexts = {row["id"]: row for row in panel["contexts"]}
    for query in panel["queries"]:
        hit = records[query["record_id"]]
        context = contexts[query["context_id"]]
        for ident in query["control_record_ids"]:
            control = records[ident]
            yield {
                "pair_id": "pair-" + stable_hash(query["id"], ident)[:24],
                "query_id": query["id"],
                "query_record_id": hit["id"],
                "query_smiles": hit["smiles"],
                "query_family": hit["family"],
                "query_scaffold_family": hit["scaffold_family"],
                "control_record_id": ident,
                "control_smiles": control["smiles"],
                "label": control["label"],
                "control_family": control["family"],
                "control_scaffold_family": control["scaffold_family"],
                "context_id": query["context_id"],
                "assay_id": query["assay_id"],
                "target_id": query["target_id"],
                "endpoint": context["endpoint"],
                "units": context["units"],
                "query_source_ids": hit["source_ids"],
                "control_source_ids": control["source_ids"],
                "assay_source_id": context["assay_metadata"]["source_id"],
            }


def public_queries(panel):
    records = {row["id"]: row for row in panel["records"]}
    return [
        {
            "id": query["id"],
            "smiles": records[query["record_id"]]["smiles"],
            "target_id": query["target_id"],
            "context_id": query["context_id"],
        }
        for query in panel["queries"]
    ]


def aggregate_scores(panel, scored_pairs, metric="score", bootstrap_seed=20261016, resamples=10000):
    """Aggregate external scores without computing or silently dropping failures."""
    from benchmarks.research.evaluation import binary_ranking_with_failures, cluster_interval

    if isinstance(resamples, bool) or not isinstance(resamples, int) or resamples < 1:
        raise ValueError("resamples must be a positive integer")
    pairs = list(control_pairs(panel))
    pair_ids = {row["pair_id"] for row in pairs}
    lookup = {}
    for row in scored_pairs:
        ident = row.get("pair_id")
        if ident not in pair_ids or ident in lookup:
            raise ValueError("Unknown or duplicate scored pair ID")
        value = row.get("scores", {}).get(metric)
        if row.get("status", "success") not in {"success", "failed", "timeout", "error", "not_evaluated"}:
            raise ValueError("Unknown scored-pair status")
        if value is not None and (
            isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
        ):
            raise ValueError("Scores must be finite numeric values or null")
        # A failed pair cannot contribute a stale, partially completed score.
        lookup[ident] = (
            None if row.get("status") in {"failed", "timeout", "error", "not_evaluated"} else value
        )
    by_query = defaultdict(list)
    for row in pairs:
        by_query[row["query_id"]].append((row["label"], lookup.get(row["pair_id"])))
    queries = []
    for query in panel["queries"]:
        report = binary_ranking_with_failures(by_query[query["id"]])
        queries.append(
            {
                "query_id": query["id"],
                "context_id": query["context_id"],
                "assay_id": query["assay_id"],
                "target_id": query["target_id"],
                "metrics": report,
            }
        )

    def average(rows):
        result = {}
        for name in ("roc_auc", "average_precision"):
            present = [row["metrics"][name] for row in rows if row["metrics"][name] is not None]
            result[name] = statistics.mean(present) if present else None
            result[name + "_estimable_children"] = len(present)
            bounds = [row["metrics"][name + "_bounds"] for row in rows]
            result[name + "_bounds"] = (
                [statistics.mean(value[i] for value in bounds) for i in (0, 1)] if bounds else None
            )
        result["children"] = len(rows)
        return result

    def group(rows, fields):
        buckets = defaultdict(list)
        for row in rows:
            buckets[tuple(row[field] for field in fields)].append(row)
        return [
            {**dict(zip(fields, key)), "metrics": average(values)} for key, values in sorted(buckets.items())
        ]

    contexts = group(queries, ("target_id", "assay_id", "context_id"))
    assays = group(contexts, ("target_id", "assay_id"))
    targets = group(assays, ("target_id",))
    macro = average(targets)
    uncertainty = {}
    for name in ("roc_auc", "average_precision"):
        result = cluster_interval(
            [row["metrics"][name] for row in targets if row["metrics"][name] is not None],
            seed=bootstrap_seed,
            resamples=resamples,
        )
        result["target_groups"] = result.pop("independent_targets")
        uncertainty[name] = result
    return {
        "schema": SCHEMA,
        "metric": metric,
        "pair_count": len(pairs),
        "missing_or_failed_pairs": sum(lookup.get(row["pair_id"]) is None for row in pairs),
        "unreturned_pairs": len(pair_ids - set(lookup)),
        "queries": queries,
        "contexts": contexts,
        "assays": assays,
        "targets": targets,
        "target_macro": macro,
        "success_conditional_target_bootstrap": uncertainty,
        "bootstrap": {"seed": bootstrap_seed, "resamples": resamples},
        "aggregation_contract": AGGREGATION,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    audit = commands.add_parser("audit", help="Verify/curate all retained activities; no scores or sampling")
    audit.add_argument("--manifest", required=True, type=Path)
    audit.add_argument("--output", required=True, type=Path)
    audit.add_argument("--active-nm", type=float, default=1000.0)
    audit.add_argument("--inactive-nm", type=float, default=10000.0)
    audit.add_argument("--max-heavy-atoms", type=int, default=60)
    select = commands.add_parser(
        "select", help="Freeze score-blind query/control choices after availability review"
    )
    select.add_argument("--pool", required=True, type=Path)
    select.add_argument("--pool-sha256", required=True, help="Pin the reviewed curated pool")
    select.add_argument("--output", required=True, type=Path)
    for name, default in asdict(Selection()).items():
        select.add_argument("--" + name.replace("_", "-"), type=int, default=default)
    summarize = commands.add_parser(
        "summarize", help="Aggregate externally computed scores; no chemical scoring"
    )
    summarize.add_argument("--panel", required=True, type=Path)
    summarize.add_argument("--scores", required=True, type=Path)
    summarize.add_argument("--metric", default="score")
    summarize.add_argument("--output", required=True, type=Path)
    summarize.add_argument("--bootstrap-seed", type=int, default=20261016)
    summarize.add_argument("--resamples", type=int, default=10000)
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error("Output already exists; retain previous audit artifacts")
    if args.command == "audit":
        pool = curate_sources(
            args.manifest,
            active_nm=args.active_nm,
            inactive_nm=args.inactive_nm,
            max_heavy_atoms=args.max_heavy_atoms,
        )
        args.output.mkdir(parents=True)
        write_json(args.output / "pool.json", pool)
        write_json(args.output / "availability.json", pool["availability"])
        write_json(args.output / "exclusions.json", pool["exclusions"])
        write_json(
            args.output / "artifacts.json",
            {"pool_sha256": digest(args.output / "pool.json"), "module_sha256": digest(__file__)},
        )
        print(json.dumps(pool["availability"]["per_target"], indent=2))
    elif args.command == "select":
        if digest(args.pool) != args.pool_sha256:
            parser.error("Curated pool SHA-256 differs from reviewed artifact")
        policy = Selection(**{key: getattr(args, key) for key in asdict(Selection())})
        panel = select_panel(json.loads(args.pool.read_text()), policy)
        validate_panel(panel, verify_sources=True)
        panel["curated_pool_sha256"] = args.pool_sha256
        args.output.mkdir(parents=True)
        write_json(args.output / "panel.json", panel)
        write_json(args.output / "control_pairs.json", list(control_pairs(panel)))
        write_json(args.output / "public_queries.json", public_queries(panel))
        write_json(
            args.output / "artifacts.json",
            {
                "panel_sha256": digest(args.output / "panel.json"),
                "module_sha256": digest(__file__),
                "selection_sha256": stable_hash(asdict(policy)),
            },
        )
        print(
            json.dumps(
                {
                    "queries": len(panel["queries"]),
                    "pairs": sum(len(q["control_record_ids"]) for q in panel["queries"]),
                    "targets": panel["selection_audit"]["targets_selected"],
                },
                indent=2,
            )
        )
    else:
        panel = json.loads(args.panel.read_text())
        validate_panel(panel, verify_sources=True)
        report = aggregate_scores(
            panel, json.loads(args.scores.read_text()), args.metric, args.bootstrap_seed, args.resamples
        )
        report.update(panel_sha256=digest(args.panel), scores_sha256=digest(args.scores))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        write_json(args.output, report)


if __name__ == "__main__":
    main()
