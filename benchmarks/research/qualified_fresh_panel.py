"""Prospective bounded ChEMBL acquisition and deterministic offline query selection.

No outcome, model, conformer or synthesis result is an input to this module.
Acquisition is an explicit command; importing or offline selection never fetches.
"""

import argparse
from collections import defaultdict
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import time
from urllib.error import HTTPError
from urllib.parse import urlencode, urljoin, urlparse
from urllib.request import Request, urlopen

from rdkit import Chem, rdBase

from benchmarks.research.assay_metadata import assess_direct_human_assay
from benchmarks.research.import_chembl import activity_label
from benchmarks.research.protocol import digest, write_json
from benchmarks.research.sources import _safe_artifact
from chembounce_synthesis.identity import chemical_identity

BASE = "https://www.ebi.ac.uk/chembl/api/data/"
RELEASE = "ChEMBL_37"
KNOWN_TARGETS = {"CHEMBL203", "CHEMBL205", "CHEMBL2147", "CHEMBL220", "CHEMBL251", "CHEMBL2637"}
POLICY = {
    "version": "qualified-fresh-panel-v1",
    "seed": 20261016,
    "targets": 6,
    "queries_per_target": 2,
    "target_cap": 60,
    "activity_cap": 1000,
    "target_catalog_cap": 25000,
    "assay_batch_size": 50,
    "release": RELEASE,
    "global_scaffold_separation": True,
    "active_nm": 1000,
    "inactive_nm": 10000,
    "max_heavy_atoms": 60,
}
SCOPE = "Absent from named resolved inspected-source exact/family/atom-bond scaffolds and explicitly target-labeled studies; not exhaustive inspection, training, library or target-family independence"


def sha(value):
    return hashlib.sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True, allow_nan=False).encode()
    ).hexdigest()


def read(path):
    def invalid(value):
        raise ValueError(f"Nonfinite JSON: {value}")

    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON field: {key}")
            result[key] = value
        return result

    return json.loads(Path(path).read_text(), parse_constant=invalid, object_pairs_hook=unique)


def bound_file(spec, base=Path(".")):
    path = (base / spec["path"]).resolve()
    if digest(path) != spec["sha256"]:
        raise ValueError(f"Source hash mismatch: {path}")
    return path


def target_identity(row):
    if (
        row.get("organism") != "Homo sapiens"
        or row.get("target_type") != "SINGLE PROTEIN"
        or row.get("species_group_flag") is not False
        or row.get("tax_id") not in (None, 9606)
    ):
        raise ValueError("not_explicit_human_single_protein")
    components = row.get("target_components")
    if not isinstance(components, list) or len(components) != 1:
        raise ValueError("missing_or_ambiguous_protein_component")
    component = components[0]
    if (
        component.get("component_type") != "PROTEIN"
        or component.get("relationship") != "SINGLE PROTEIN"
        or not isinstance(component.get("accession"), str)
        or not component["accession"].strip()
        or type(component.get("component_id")) is not int
        or component["component_id"] <= 0
    ):
        raise ValueError("missing_or_conflicting_accession_component_mapping")
    ident = row.get("target_chembl_id")
    if not isinstance(ident, str) or not ident.startswith("CHEMBL") or not ident[6:].isdigit():
        raise ValueError("invalid_target_identifier")
    aliases = {row.get("pref_name"), component.get("component_description")}
    aliases.update(s.get("component_synonym") for s in component.get("target_component_synonyms", []))
    return {
        "target_id": ident,
        "accession": component["accession"],
        "component_id": component["component_id"],
        "aliases": sorted({s.strip().casefold() for s in aliases if isinstance(s, str) and s.strip()}),
        "mapping_scope": "sourced accession/component identity; sequence bytes not supplied by target endpoint",
    }


def load_exclusions(spec):
    """Require sourced resolution for every known ID/alias before any selection."""
    manifest_path = bound_file(spec["preparation_manifest"])
    manifest = read(manifest_path)
    union = manifest["final_union"]
    ids = set(union["excluded_target_ids"])
    if not KNOWN_TARGETS <= ids:
        raise ValueError("Exclusion manifest omits known inspected target IDs")
    aliases = {s.strip().casefold() for s in union["target_aliases"]}
    artifacts = {}
    for name in (union["identities_artifact"], "normalization-failures.private.json"):
        path = manifest_path.parent / name
        if digest(path) != manifest["artifacts"][name]["sha256"]:
            raise ValueError("Exclusion identity/failure artifact changed")
        artifacts[name] = read(path)
    identities = artifacts[union["identities_artifact"]]
    expected = manifest["reviewed_identity_kernel"]
    root = Path(__file__).resolve().parents[2]
    if set(expected) != {"chembounce_synthesis/common.py", "chembounce_synthesis/identity.py"} or any(
        digest(root / name) != value for name, value in expected.items()
    ):
        raise ValueError("Exclusion identity kernel changed")
    if rdBase.rdkitVersion != manifest["runtime"]["rdkit"]:
        raise ValueError("Exclusion identity runtime changed")
    resolved_ids, resolved_aliases, accessions, components, mappings, seen = (
        set(),
        set(),
        set(),
        set(),
        [],
        {},
    )
    for source in spec["target_sources"]:
        path = bound_file(source)
        if not source.get("url", "").startswith(BASE):
            raise ValueError("Target alias mapping requires retained official ChEMBL provenance")
        payload = read(path)
        rows = payload.get("targets", [payload])
        for row in rows:
            ident = target_identity(row)
            key = ident["target_id"]
            if key in seen and seen[key] != row:
                raise ValueError("Conflicting excluded target metadata")
            seen[key] = row
            matching = aliases.intersection(ident["aliases"])
            if key in ids or matching:
                resolved_ids.add(key)
                resolved_aliases.update(matching)
                accessions.add(ident["accession"])
                components.add(ident["component_id"])
                mappings.append({**ident, "source_sha256": source["sha256"], "source_url": source["url"]})
    alias_identities = defaultdict(set)
    for mapping in mappings:
        for alias in aliases.intersection(mapping["aliases"]):
            alias_identities[alias].add((mapping["accession"], mapping["component_id"]))
    if any(len(values) != 1 for values in alias_identities.values()):
        raise ValueError("Conflicting required target alias/accession mapping")
    if ids - resolved_ids or aliases - resolved_aliases:
        raise ValueError(
            "Unresolved required target exclusion mappings: "
            + repr(sorted(ids - resolved_ids))
            + "; aliases="
            + repr(sorted(aliases - resolved_aliases))
        )
    return {
        "exact": {r["exact_smiles"] for r in identities},
        "families": {r["family_smiles"] for r in identities},
        "scaffolds": {r["scaffold_family_smiles"] for r in identities if r["scaffold_family_smiles"]},
        "generic": {r["ring_skeleton_smiles"] for r in identities if r["ring_skeleton_smiles"]},
        "unresolved_inputs": {r["input_smiles"] for r in artifacts["normalization-failures.private.json"]},
        "target_ids": ids | resolved_ids,
        "accessions": accessions,
        "components": components,
        "provenance": {
            "spec_sha256": sha(spec),
            "manifest_sha256": digest(manifest_path),
            "target_mappings": mappings,
            "unresolved_normalization_counts": manifest["failed_unique_inputs_by_reason"],
            "scope": SCOPE,
            "remaining_gaps": spec["remaining_gaps"],
        },
    }


def rank(kind, *parts):
    return sha([POLICY["seed"], kind, *parts])


def ordered_targets(rows, exclusions):
    if len(rows) > POLICY["target_catalog_cap"]:
        raise ValueError("Target catalog exceeds frozen cap")
    seen, eligible, ledger = {}, [], []
    for row in rows:
        key = row.get("target_chembl_id")
        if key in seen:
            raise ValueError("Duplicate/conflicting target catalog record")
        seen[key] = row
        entry = {"target_id": key, "record_sha256": sha(row)}
        try:
            ident = target_identity(row)
            entry.update(ident)
            if (
                key in exclusions["target_ids"]
                or ident["accession"] in exclusions["accessions"]
                or ident["component_id"] in exclusions["components"]
            ):
                raise ValueError("inspected_target_or_same_protein_alias")
            entry.update(status="eligible", rank_sha256=rank("fresh-target-v1", key))
            eligible.append(entry)
        except (ValueError, KeyError, TypeError) as exc:
            entry.update(status="rejected", reason=str(exc))
        ledger.append(entry)
    by_accession, by_component = defaultdict(set), defaultdict(set)
    for item in ledger:
        if "accession" in item:
            by_accession[item["accession"]].add(item["component_id"])
            by_component[item["component_id"]].add(item["accession"])
    for item in eligible:
        if len(by_accession[item["accession"]]) != 1 or len(by_component[item["component_id"]]) != 1:
            item.update(status="rejected", reason="conflicting_catalog_accession_component_mapping")
    eligible = [row for row in eligible if row["status"] == "eligible"]
    eligible.sort(key=lambda row: (row["rank_sha256"], row["target_id"]))
    return eligible, ledger


def curate_target(target, activities, assay_records, exclusions):
    """All bounded input records are assessed before family conflict decisions."""
    if len(activities) > POLICY["activity_cap"]:
        raise ValueError("Activity sample exceeds frozen cap")
    assay_map = {}
    for entry in assay_records:
        ident = entry["record"].get("assay_chembl_id")
        if not isinstance(ident, str) or not ident:
            raise ValueError("Invalid assay join identifier")
        if ident in assay_map and assay_map[ident]["record"] != entry["record"]:
            raise ValueError("Conflicting assay join records")
        assay_map.setdefault(ident, entry)
    groups, ledger, seen = defaultdict(list), [], {}
    previous = -1
    for row in activities:
        activity_id = row.get("activity_id")
        if type(activity_id) is not int or activity_id < 0 or activity_id < previous:
            raise ValueError("Activities must be ordered by nonnegative integer activity_id")
        previous = activity_id
        entry = {"activity_id": activity_id, "record_sha256": sha(row), "status": "rejected"}
        ledger.append(entry)
        if activity_id in seen:
            if seen[activity_id] != row:
                raise ValueError("Conflicting duplicate activity ID")
            entry["reason"] = "identical_duplicate_activity"
            continue
        seen[activity_id] = row
        if row.get("target_chembl_id") != target:
            raise ValueError("Activity target differs from requested sample")
        label, reason = activity_label(row, POLICY["active_nm"], POLICY["inactive_nm"])
        if not reason:
            accepted, reason, assignment = assess_direct_human_assay(row, assay_map)
        if reason:
            entry["reason"] = reason
            continue
        smiles = row.get("canonical_smiles")
        try:
            if smiles in exclusions["unresolved_inputs"]:
                raise ValueError("unresolved_inspected_identity")
            mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
            if (
                mol is None
                or len(Chem.GetMolFrags(mol)) != 1
                or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())
            ):
                raise ValueError("invalid_disconnected_or_dummy_structure")
            if mol.GetNumHeavyAtoms() > POLICY["max_heavy_atoms"]:
                raise ValueError("heavy_atom_cap")
            if type(row.get("document_year")) is not int:
                raise ValueError("publication_year_missing_or_invalid")
            if any(str(s.specified) != "Specified" for s in Chem.FindPotentialStereo(mol)):
                raise ValueError("unspecified_or_unknown_stereochemistry")
            identity = chemical_identity(smiles)
            if not identity.has_scaffold:
                raise ValueError("no_normalized_ring_scaffold")
            entry.update(identity=asdict(identity), label=label, assay_assignment=assignment)
            # Group before historical exclusions so another state cannot hide a conflicting family label.
            groups[identity.family_smiles].append((row, entry, identity))
            entry["status"] = "curated"
        except (ValueError, RuntimeError, TypeError) as exc:
            entry["reason"] = str(exc)
    representatives = []
    for family, items in sorted(groups.items()):
        reason = None
        if len({e[1]["label"] for e in items}) != 1 or len({e[0]["standard_type"] for e in items}) != 1:
            reason = "conflicting_family_labels_or_endpoints"
        elif items[0][1]["label"] != "active":
            reason = "inactive_family_not_query"
        items.sort(key=lambda item: (item[0]["document_year"], item[0]["activity_id"]))
        row, entry, identity = items[0]
        if reason is None:
            if any(item[2].exact_smiles in exclusions["exact"] for item in items):
                reason = "inspected_exact_state"
            elif family in exclusions["families"]:
                reason = "inspected_chemical_family"
            elif identity.scaffold_family_smiles in exclusions["scaffolds"]:
                reason = "inspected_normalized_scaffold"
        for _, item, _ in items:
            item.update(
                status="rejected" if reason else "family_support", reason=reason or "same_family_support"
            )
        if reason:
            continue
        entry.update(status="eligible_representative", reason=None)
        representatives.append(
            {
                "target_id": target,
                "smiles": identity.exact_smiles,
                "family": family,
                "scaffold": identity.scaffold_family_smiles,
                "generic_scaffold": identity.ring_skeleton_smiles,
                "generic_overlap_with_inspected": identity.ring_skeleton_smiles in exclusions["generic"],
                "activity_id": row["activity_id"],
                "assay_id": row["assay_chembl_id"],
                "rank_sha256": rank("fresh-query-v1", target, family),
            }
        )
    representatives.sort(key=lambda row: (row["rank_sha256"], row["family"], row["smiles"]))
    return representatives, ledger


def select_panel(targets, get_target_sample, exclusions):
    ordered, catalog_ledger = ordered_targets(targets, exclusions)
    chosen, families, scaffolds, proteins, target_ledger = [], set(), set(), set(), []
    for target in ordered[: POLICY["target_cap"]]:
        target_id = target["target_id"]
        entry = {"target_id": target_id, "rank_sha256": target["rank_sha256"], "status": "rejected"}
        target_ledger.append(entry)
        if target["accession"] in proteins:
            entry["reason"] = "same_protein_as_selected_target"
            continue
        activities, assays, provenance = get_target_sample(target_id)
        candidates, activities_ledger = curate_target(target_id, activities, assays, exclusions)
        entry.update(provenance=provenance, activity_decisions=activities_ledger, query_decisions=[])
        local, local_scaffolds = [], set()
        for candidate in candidates:
            reason = (
                "selected_family_overlap"
                if candidate["family"] in families
                else "global_scaffold_overlap"
                if candidate["scaffold"] in scaffolds | local_scaffolds
                else "outside_two_query_quota"
                if len(local) == POLICY["queries_per_target"]
                else None
            )
            entry["query_decisions"].append(
                {**candidate, "status": "not_selected" if reason else "tentative", "reason": reason}
            )
            if reason is None:
                local.append(candidate)
                local_scaffolds.add(candidate["scaffold"])
        for decision in entry["query_decisions"]:
            if decision["status"] == "tentative":
                decision["status"] = (
                    "selected" if len(local) == POLICY["queries_per_target"] else "not_selected"
                )
                decision["reason"] = (
                    None
                    if len(local) == POLICY["queries_per_target"]
                    else "target_has_insufficient_global_scaffolds"
                )
        if len(local) != POLICY["queries_per_target"]:
            entry["reason"] = "fewer_than_two_available_global_scaffolds"
            continue
        entry["status"] = "selected"
        for index, row in enumerate(local, 1):
            chosen.append(
                {"id": f"{target_id}-query-{index}", "target_id": target_id, "smiles": row["smiles"]}
            )
            families.add(row["family"])
            scaffolds.add(row["scaffold"])
        proteins.add(target["accession"])
        if len(proteins) == POLICY["targets"]:
            break
    complete = len(chosen) == POLICY["targets"] * POLICY["queries_per_target"]
    return {
        "schema": "qualified-fresh-panel-selection-v1",
        "status": "complete" if complete else "infeasible",
        "policy": POLICY,
        "scope": SCOPE,
        "queries": chosen if complete else [],
        "partial_queries_not_authorized": [] if complete else chosen,
        "catalog_decisions": catalog_ledger,
        "target_decisions": target_ledger,
        "eligible_catalog_targets": len(ordered),
        "examined_targets": len(target_ledger),
        "eligible_targets_not_examined": [r["target_id"] for r in ordered[len(target_ledger) :]],
        "exclusions": exclusions["provenance"],
        "no_replacement_policy": "No method outcomes are inputs; no post-generation query replacement",
    }


def api_url(resource, **query):
    return BASE + resource + ".json" + ("?" + urlencode(sorted(query.items())) if query else "")


def official_transport(url, timeout, max_bytes):
    """One explicit official request; no retry and no response rewriting."""
    request = Request(
        url, headers={"User-Agent": "ChemBounce-qualified-fresh-panel/1", "Accept": "application/json"}
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            final_url = response.geturl()
            if not final_url.startswith(BASE):
                raise ValueError("ChEMBL response redirected outside official API")
            data = response.read(max_bytes + 1)
            return data, final_url, response.status, dict(response.headers)
    except HTTPError as exc:
        # Preserve the bounded actual error body and status without retrying.
        return exc.read(max_bytes + 1), exc.geturl(), exc.code, dict(exc.headers)


def validate_protocol(protocol):
    import platform

    if protocol.get("schema") != "qualified-fresh-panel-protocol-v1" or protocol.get("policy") != POLICY:
        raise ValueError("Fresh-panel protocol differs from frozen selection policy")
    root = Path(__file__).resolve().parents[2]
    required = {
        "benchmarks/research/qualified_fresh_panel.py",
        "benchmarks/research/import_chembl.py",
        "benchmarks/research/assay_metadata.py",
        "chembounce_synthesis/identity.py",
        "chembounce_synthesis/common.py",
    }
    sources = protocol["source_sha256"]
    if not required <= sources.keys() or any(digest(root / p) != value for p, value in sources.items()):
        raise ValueError("Fresh-panel source pins missing or changed")
    if protocol["runtime"] != {"python": platform.python_version(), "rdkit": rdBase.rdkitVersion}:
        raise ValueError("Fresh-panel numerical runtime changed")
    bounds = protocol["acquisition_limits"]
    for key, maximum in {
        "max_requests": 1287,
        "max_response_bytes": 100_000_000,
        "max_total_bytes": 10_000_000_000,
        "timeout_seconds": 60,
        "wall_seconds": 86400,
    }.items():
        value = bounds.get(key)
        if type(value) is not int or not 0 < value <= maximum:
            raise ValueError(f"Missing or invalid acquisition bound: {key}")
    return load_exclusions(protocol["exclusions"])


class Acquisition:
    def __init__(self, output, protocol, transport):
        self.output, self.protocol, self.transport = output, protocol, transport
        self.started, self.total_bytes = time.monotonic(), 0
        self.records = []

    def fetch(self, resource, role, **query):
        limits = self.protocol["acquisition_limits"]
        remaining = limits["wall_seconds"] - (time.monotonic() - self.started)
        if len(self.records) >= limits["max_requests"] or remaining <= 0:
            raise ValueError("Acquisition request/wall bound exhausted; no retry")
        url = api_url(resource, **query)
        record = {
            "index": len(self.records),
            "role": role,
            "requested_url": url,
            "parameters": query,
            "resource": resource,
            "status": "started",
        }
        self.records.append(record)
        write_json(self.output / "requests.json", self.records)
        try:
            raw, final, status, headers = self.transport(
                url, min(limits["timeout_seconds"], remaining), limits["max_response_bytes"]
            )
            artifact = f"raw-{record['index']:04d}.json"
            (self.output / artifact).write_bytes(raw)
            record.update(
                artifact=artifact,
                sha256=digest(self.output / artifact),
                bytes=len(raw),
                final_url=final,
                http_status=status,
                headers=headers,
            )
            self.total_bytes += len(raw)
            if (
                not final.startswith(BASE)
                or status != 200
                or len(raw) > limits["max_response_bytes"]
                or self.total_bytes > limits["max_total_bytes"]
                or time.monotonic() - self.started > limits["wall_seconds"]
            ):
                raise ValueError("Official response status/location/byte/wall bound failed")
            payload = read(self.output / artifact)
            record["status"] = "complete"
            return payload, record
        except Exception as exc:
            record.update(status="failed", reason=f"{type(exc).__name__}: {exc}")
            raise
        finally:
            write_json(self.output / "requests.json", self.records)


def verify_page(data, key, offset, limit, expected_url, *, total=None):
    from urllib.parse import parse_qs

    meta, rows = data["page_meta"], data[key]
    if (
        not isinstance(rows, list)
        or type(meta.get("total_count")) is not int
        or meta["total_count"] < 0
        or type(meta.get("offset")) is not int
        or type(meta.get("limit")) is not int
        or meta.get("offset") != offset
        or meta.get("limit") != limit
        or len(rows) != min(limit, max(0, meta["total_count"] - offset))
        or (total is not None and meta["total_count"] != total)
    ):
        raise ValueError("Incomplete or inconsistent official pagination")
    next_url = meta.get("next")
    more = offset + len(rows) < meta["total_count"]
    if bool(next_url) != more:
        raise ValueError("Missing or unexpected next-page link")
    if more:
        next_url = urljoin(BASE, next_url)
        actual, expected = urlparse(next_url), urlparse(expected_url)
        parameters = parse_qs(expected.query)
        parameters["offset"] = [str(offset + limit)]
        if (actual.scheme, actual.netloc, actual.path) != (
            expected.scheme,
            expected.netloc,
            expected.path,
        ) or parse_qs(actual.query) != parameters:
            raise ValueError("Pagination changed official endpoint or frozen query")
    return rows, meta["total_count"]


def assay_entries(payload, source, requested):
    rows, _ = verify_page(payload, "assays", 0, 50, source["requested_url"])
    if payload["page_meta"]["total_count"] > len(requested) or any(
        row.get("assay_chembl_id") not in requested for row in rows
    ):
        raise ValueError("Assay response contains an unrequested ID")
    return [
        {
            "record": row,
            "source_sha256": source["sha256"],
            "source_url": source["requested_url"],
            "source_artifact": source["artifact"],
        }
        for row in rows
    ]


def acquire(protocol, output, *, transport=official_transport):
    """Explicit bounded acquisition; the shared pure selector determines only its frozen early stop."""
    exclusions = validate_protocol(protocol)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    acquisition = Acquisition(output, protocol, transport)
    snapshot = {
        "schema": "qualified-fresh-acquisition-v1",
        "status": "incomplete",
        "protocol_sha256": sha(protocol),
        "policy": POLICY,
        "target_pages": [],
        "target_samples": {},
        "sources": acquisition.records,
    }
    write_json(output / "protocol.json", protocol)
    try:
        before, source = acquisition.fetch("status", "release_before")
        snapshot["release_before"] = source["index"]
        if before.get("chembl_db_version") != RELEASE or before.get("status") != "UP":
            raise ValueError("ChEMBL release/status differs before acquisition")
        targets, offset, total = [], 0, None
        while total is None or offset < total:
            payload, source = acquisition.fetch(
                "target",
                "target_catalog",
                organism="Homo sapiens",
                target_type="SINGLE PROTEIN",
                limit=1000,
                offset=offset,
                order_by="target_chembl_id",
            )
            rows, total = verify_page(payload, "targets", offset, 1000, source["requested_url"], total=total)
            if total > POLICY["target_catalog_cap"]:
                raise ValueError("Complete target catalog exceeds declared cap")
            targets.extend(rows)
            snapshot["target_pages"].append(source["index"])
            offset += len(rows)
            if not rows:
                break
        if [row["target_chembl_id"] for row in targets] != sorted(row["target_chembl_id"] for row in targets):
            raise ValueError("Target catalog does not honor frozen order")

        def sample(target):
            payload, source = acquisition.fetch(
                "activity",
                "target_activities",
                target_chembl_id=target,
                assay_type="B",
                standard_flag=1,
                standard_type__in="IC50,Ki,Kd,EC50",
                standard_units="nM",
                limit=1000,
                offset=0,
                order_by="activity_id",
            )
            rows, count = verify_page(payload, "activities", 0, 1000, source["requested_url"])
            indices, assays = [], []
            assay_ids = sorted(
                {
                    r["assay_chembl_id"]
                    for r in rows
                    if isinstance(r.get("assay_chembl_id"), str) and r["assay_chembl_id"]
                }
            )
            for start in range(0, len(assay_ids), 50):
                batch = assay_ids[start : start + 50]
                data, record = acquisition.fetch(
                    "assay",
                    "assay_metadata",
                    assay_chembl_id__in=",".join(batch),
                    limit=50,
                    offset=0,
                    order_by="assay_chembl_id",
                )
                assays.extend(assay_entries(data, record, set(batch)))
                indices.append(record["index"])
            provenance = {
                "activity_source": source["index"],
                "assay_sources": indices,
                "total_reported_activities": count,
                "retained_first": len(rows),
                "unconsumed_tail": max(0, count - len(rows)),
            }
            snapshot["target_samples"][target] = provenance
            return rows, assays, provenance

        draft = select_panel(targets, sample, exclusions)
        after, source = acquisition.fetch("status", "release_after")
        snapshot["release_after"] = source["index"]
        if after.get("chembl_db_version") != RELEASE or after.get("status") != "UP":
            raise ValueError("ChEMBL release/status changed after acquisition")
        snapshot.update(
            status="complete",
            selection_status=draft["status"],
            total_bytes=acquisition.total_bytes,
            elapsed_seconds=time.monotonic() - acquisition.started,
        )
        write_json(output / "acquisition-selection.private.json", draft)
        snapshot["acquisition_selection_sha256"] = digest(output / "acquisition-selection.private.json")
    except Exception as exc:
        snapshot["failure"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        write_json(output / "snapshot.json", snapshot)
    return snapshot


def select_offline(protocol, snapshot_dir, output):
    exclusions = validate_protocol(protocol)
    root = Path(snapshot_dir).resolve()
    snapshot = read(root / "snapshot.json")
    if (
        snapshot.get("schema") != "qualified-fresh-acquisition-v1"
        or snapshot.get("status") != "complete"
        or snapshot["protocol_sha256"] != sha(protocol)
        or snapshot["policy"] != POLICY
    ):
        raise ValueError("Incomplete or mismatched acquisition snapshot")
    records = snapshot["sources"]
    limits = protocol["acquisition_limits"]
    if (
        snapshot["total_bytes"] != sum(r["bytes"] for r in records)
        or len(records) > limits["max_requests"]
        or sum(r["bytes"] for r in records) > limits["max_total_bytes"]
        or any(r["bytes"] > limits["max_response_bytes"] for r in records)
        or snapshot["release_before"] != 0
        or snapshot["release_after"] != len(records) - 1
    ):
        raise ValueError("Retained acquisition bounds/release ordering failed")
    used = []

    def expect(index, resource, role, parameters):
        source = records[index]
        if source["resource"] != resource or source["role"] != role or source["parameters"] != parameters:
            raise ValueError("Retained source differs from frozen acquisition request")
        used.append(index)

    payloads = []
    for index, source in enumerate(records):
        path = _safe_artifact(root, source["artifact"])
        if (
            source["index"] != index
            or source["status"] != "complete"
            or digest(path) != source["sha256"]
            or path.stat().st_size != source["bytes"]
            or source["http_status"] != 200
        ):
            raise ValueError("Raw acquisition bytes/provenance changed")
        if source["requested_url"] != api_url(source["resource"], **source["parameters"]) or not source[
            "final_url"
        ].startswith(BASE):
            raise ValueError("Raw acquisition URL binding changed")
        payloads.append(read(path))
    for key in ("release_before", "release_after"):
        expect(snapshot[key], "status", key, {})
        payload = payloads[snapshot[key]]
        if payload.get("chembl_db_version") != RELEASE or payload.get("status") != "UP":
            raise ValueError("Retained release guard failed")
    targets, total = [], None
    for index in snapshot["target_pages"]:
        expect(
            index,
            "target",
            "target_catalog",
            {
                "organism": "Homo sapiens",
                "target_type": "SINGLE PROTEIN",
                "limit": 1000,
                "offset": len(targets),
                "order_by": "target_chembl_id",
            },
        )
        rows, total = verify_page(
            payloads[index], "targets", len(targets), 1000, records[index]["requested_url"], total=total
        )
        targets.extend(rows)
    if (
        total is None
        or len(targets) != total
        or [r["target_chembl_id"] for r in targets] != sorted(r["target_chembl_id"] for r in targets)
    ):
        raise ValueError("Incomplete target catalog")
    consumed = []

    def sample(target):
        if target not in snapshot["target_samples"]:
            raise ValueError("Missing required ranked target sample; no replacement")
        provenance = snapshot["target_samples"][target]
        index = provenance["activity_source"]
        source = records[index]
        expect(
            index,
            "activity",
            "target_activities",
            {
                "target_chembl_id": target,
                "assay_type": "B",
                "standard_flag": 1,
                "standard_type__in": "IC50,Ki,Kd,EC50",
                "standard_units": "nM",
                "limit": 1000,
                "offset": 0,
                "order_by": "activity_id",
            },
        )
        if source["parameters"].get("target_chembl_id") != target:
            raise ValueError("Retained activity source target mismatch")
        rows, count = verify_page(payloads[index], "activities", 0, 1000, source["requested_url"])
        expected_provenance = {
            "activity_source": index,
            "assay_sources": provenance["assay_sources"],
            "total_reported_activities": count,
            "retained_first": len(rows),
            "unconsumed_tail": max(0, count - len(rows)),
        }
        if provenance != expected_provenance:
            raise ValueError("Retained activity sample counts disagree with raw pagination")
        ids = sorted(
            {
                r["assay_chembl_id"]
                for r in rows
                if isinstance(r.get("assay_chembl_id"), str) and r["assay_chembl_id"]
            }
        )
        batches = [ids[i : i + 50] for i in range(0, len(ids), 50)]
        if len(batches) != len(provenance["assay_sources"]):
            raise ValueError("Retained assay batches incomplete")
        assays = []
        for batch, index in zip(batches, provenance["assay_sources"]):
            source = records[index]
            expect(
                index,
                "assay",
                "assay_metadata",
                {
                    "assay_chembl_id__in": ",".join(batch),
                    "limit": 50,
                    "offset": 0,
                    "order_by": "assay_chembl_id",
                },
            )
            if source["parameters"].get("assay_chembl_id__in") != ",".join(batch):
                raise ValueError("Assay batch differs from frozen lexicographic join")
            assays.extend(assay_entries(payloads[index], source, set(batch)))
        consumed.append(target)
        return rows, assays, provenance

    selection = select_panel(targets, sample, exclusions)
    if set(consumed) != set(snapshot["target_samples"]):
        raise ValueError(
            "Snapshot includes additional outcome-independent target samples outside frozen stop"
        )
    expected_order = [used[0], *used[2:], used[1]]
    if expected_order != list(range(len(records))):
        raise ValueError("Missing, duplicate, reordered or extra acquisition source calls")
    draft = root / "acquisition-selection.private.json"
    if digest(draft) != snapshot["acquisition_selection_sha256"] or read(draft) != selection:
        raise ValueError("Offline selection differs from frozen acquisition decision ledger")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "selection.private.json", selection)
    write_json(
        output / "public-queries.json",
        {
            "schema": "qualified-public-queries-v1",
            "queries": selection["queries"],
            "selection_sha256": digest(output / "selection.private.json"),
        },
    )
    write_json(
        output / "provenance.json",
        {
            "protocol_sha256": sha(protocol),
            "snapshot_sha256": digest(root / "snapshot.json"),
            "status": selection["status"],
            "scope": SCOPE,
        },
    )
    return selection


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("acquire", "select"):
        item = sub.add_parser(name)
        item.add_argument("--protocol", required=True)
        item.add_argument("--output", required=True)
        if name == "select":
            item.add_argument("--snapshot", required=True)
    args = parser.parse_args(argv)
    protocol = read(args.protocol)
    if args.command == "acquire":
        acquire(protocol, args.output)
    else:
        select_offline(protocol, args.snapshot, args.output)


if __name__ == "__main__":
    main()
