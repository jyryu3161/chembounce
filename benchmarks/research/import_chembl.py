"""Curate pinned ChEMBL activity samples into leakage-audited retrospective panels.

Only explicitly sourced assay measurements produce active/inactive labels. The
retained publication year is not represented as an invented experimental date.
"""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import sys

from rdkit import Chem

from benchmarks.research.protocol import digest, load_dataset, write_json
from chembounce_synthesis.identity import chemical_identity


def activity_label(record, active_nm=1000.0, inactive_nm=10000.0):
    """Respect inequality direction; an upper bound cannot prove inactivity."""
    if not 0 < active_nm < inactive_nm:
        raise ValueError("Require 0 < active_nm < inactive_nm")
    if record.get("standard_units") != "nM" or record.get("standard_type") not in {
        "IC50",
        "Ki",
        "Kd",
        "EC50",
    }:
        return None, "unsupported_endpoint_or_units"
    if record.get("standard_flag") != 1 or record.get("assay_type") != "B":
        return None, "nonstandard_or_nonbinding_assay"
    if record.get("data_validity_comment") or record.get("potential_duplicate"):
        return None, "validity_flag_or_potential_duplicate"
    if record.get("assay_variant_mutation") or record.get("assay_variant_accession"):
        return None, "variant_target"
    if record.get("target_organism") != "Homo sapiens":
        return None, "nonhuman_or_unknown_organism"
    if record.get("standard_upper_value") is not None:
        return None, "range_measurement_not_supported"
    try:
        value = float(record["standard_value"])
    except (TypeError, KeyError, ValueError):
        return None, "missing_numeric_value"
    if not math.isfinite(value) or value <= 0:
        return None, "invalid_numeric_value"
    relation = record.get("standard_relation")
    if relation in {"=", "<", "<="} and value <= active_nm:
        return "active", None
    if relation in {"=", ">", ">="} and value >= inactive_nm:
        return "inactive", None
    return None, "ambiguous_or_intermediate_activity"


def curate(
    manifest_path,
    output,
    split_kind="scaffold",
    seed=42,
    queries_per_target=2,
    active_nm=1000.0,
    inactive_nm=10000.0,
    cutoff_year=2020,
    max_heavy_atoms=60,
    require_direct_human_assay=True,
):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    assay_metadata, assay_provenance = {}, {"status": "not_checked_diagnostic_only"}
    if require_direct_human_assay:
        from benchmarks.research.assay_metadata import load_assay_metadata, assess_direct_human_assay

        assay_metadata, assay_provenance = load_assay_metadata(manifest_path)
        if not assay_metadata:
            raise ValueError(
                "Strict panel requires retained assay metadata; use diagnostic override only for an explicitly unverified panel"
            )
    sources, raw, excluded = [], [], []
    for source_index, source in enumerate(manifest["sources"], 1):
        print(
            f"Source {source_index}/{len(manifest['sources'])}: {source['name']}", file=sys.stderr, flush=True
        )
        path = Path(source["path"])
        if not path.is_file():
            path = manifest_path.parent / Path(source["path"]).name
        path = path.resolve()
        if digest(path) != source["sha256"]:
            raise ValueError(f"Raw source hash mismatch: {path}")
        data = json.loads(path.read_text())
        source_id = source["name"]
        sources.append(
            {"id": source_id, "url": source["url"], "artifact": str(path), "sha256": source["sha256"]}
        )
        activities = data.get("activities", [])
        if not isinstance(activities, list):
            continue
        for row in activities:
            label, reason = activity_label(row, active_nm, inactive_nm)
            key = str(row.get("activity_id"))
            if reason:
                excluded.append({"activity_id": key, "source_id": source_id, "reason": reason})
                continue
            assignment = {"status": "not_checked_diagnostic_only"}
            if require_direct_human_assay:
                passed, reason, assignment = assess_direct_human_assay(row, assay_metadata)
                if not passed:
                    excluded.append({"activity_id": key, "source_id": source_id, "reason": reason})
                    continue
            smiles = row.get("canonical_smiles")
            mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
            cheap_reason = None
            if mol is None:
                cheap_reason = "invalid_smiles"
            elif len(Chem.GetMolFrags(mol)) != 1:
                cheap_reason = "disconnected_chemical_scope"
            elif mol.GetNumHeavyAtoms() > max_heavy_atoms:
                cheap_reason = "declared_heavy_atom_scope"
            elif not isinstance(row.get("document_year"), int):
                cheap_reason = "publication_year_missing"
            if cheap_reason:
                excluded.append({"activity_id": key, "source_id": source_id, "reason": cheap_reason})
                continue
            try:
                identity = chemical_identity(smiles)
            except (KeyError, ValueError, RuntimeError) as exc:
                excluded.append(
                    {
                        "activity_id": key,
                        "source_id": source_id,
                        "reason": "chemical_identity_failed",
                        "error": str(exc),
                    }
                )
                continue
            mol = Chem.MolFromSmiles(identity.exact_smiles)
            reason = None
            if not identity.has_scaffold:
                reason = "acyclic_not_scaffold_hopping_panel"
            elif mol.GetNumHeavyAtoms() > max_heavy_atoms:
                reason = "declared_heavy_atom_scope"
            elif not isinstance(row.get("document_year"), int):
                reason = "publication_year_missing"
            if reason:
                excluded.append({"activity_id": key, "source_id": source_id, "reason": reason})
                continue
            raw.append(
                {
                    "row": row,
                    "identity": identity,
                    "label": label,
                    "source_id": source_id,
                    "assay_assignment": assignment,
                }
            )
    # Collapse chemical families only after inspecting every retained activity.
    # Conflicting labels and mixed endpoints within one family-target are excluded.
    print(
        f"Grouping {len(raw)} qualifying measurements; {len(excluded)} exclusions recorded",
        file=sys.stderr,
        flush=True,
    )
    groups = defaultdict(list)
    for item in raw:
        group = (item["row"]["target_chembl_id"], item["identity"].family_smiles)
        if split_kind == "time":
            # A later contradictory assay must not retroactively remove/change
            # a historical development record. Curate each era independently.
            group += (item["row"]["document_year"] > cutoff_year,)
        groups[group].append(item)
    curated = []
    for key, items in sorted(groups.items()):
        if len({i["label"] for i in items}) != 1 or len({i["row"]["standard_type"] for i in items}) != 1:
            excluded.extend(
                {
                    "activity_id": str(i["row"]["activity_id"]),
                    "source_id": i["source_id"],
                    "reason": "conflicting_family_labels_or_endpoints",
                }
                for i in items
            )
            continue
        items.sort(key=lambda i: (i["row"]["document_year"], i["row"]["activity_id"]))
        chosen = items[0]
        row, identity = chosen["row"], chosen["identity"]
        family_id = hashlib.sha256(identity.family_smiles.encode()).hexdigest()[:16]
        curated.append(
            {
                "id": f"{row['target_chembl_id']}-{family_id}",
                "smiles": identity.exact_smiles,
                "target_id": row["target_chembl_id"],
                "label": chosen["label"],
                "source_id": chosen["source_id"],
                "assay_id": row["assay_chembl_id"],
                "assay_assignment": chosen["assay_assignment"],
                "label_basis": f"{'Direct human assay' if require_direct_human_assay else 'Unverified assay assignment, human target annotation only'} binding {row['standard_type']} in nM; active <= {active_nm:g}, inactive >= {inactive_nm:g}; censored relations respected",
                "publication_year": row["document_year"],
                "molecule_chembl_id": row["molecule_chembl_id"],
                "activity": {
                    "type": row["standard_type"],
                    "relation": row["standard_relation"],
                    "value": float(row["standard_value"]),
                    "units": "nM",
                    "organism": row["target_organism"],
                },
                "supporting_activity_ids": sorted({i["row"]["activity_id"] for i in items}),
                "supporting_assay_ids": sorted({i["row"]["assay_chembl_id"] for i in items}),
                "supporting_source_ids": sorted({i["source_id"] for i in items}),
                "_family": identity.family_smiles,
                "_scaffold": identity.scaffold_family_smiles,
            }
        )
    targets = sorted({r["target_id"] for r in curated})
    rng = random.Random(seed)
    shuffled_targets = targets.copy()
    rng.shuffle(shuffled_targets)
    test_targets = set(shuffled_targets[: max(1, len(targets) // 2)])
    for row in curated:
        if split_kind == "scaffold":
            value = hashlib.sha256(f"{seed}:{row['_scaffold']}".encode()).digest()[0] / 256
            row["split"] = "test" if value < 0.25 else "train"
        elif split_kind == "target":
            row["split"] = "test" if row["target_id"] in test_targets else "train"
        elif split_kind == "time":
            row["split"] = "train" if row["publication_year"] <= cutoff_year else "test"
        else:
            raise ValueError("Unknown split_kind")
    # Families observed in development are removed from later held-out partitions
    # across targets too. Do not silently move observations across their dates.
    train_families = {r["_family"] for r in curated if r["split"] == "train"}
    retained = []
    for row in curated:
        if row["split"] == "test" and row["_family"] in train_families:
            excluded.append({"compound_id": row["id"], "reason": "family_seen_in_development"})
        else:
            retained.append(row)
    curated = retained
    queries = []
    for target in targets:
        refs = [
            r for r in curated if r["target_id"] == target and r["split"] == "test" and r["label"] == "active"
        ]
        options = [
            r
            for r in curated
            if r["target_id"] == target
            and r["label"] == "active"
            and r["split"] == ("test" if split_kind == "target" else "train")
        ]
        rng.shuffle(options)
        chosen_scaffolds = set()
        for row in options:
            if row["_scaffold"] in chosen_scaffolds or not any(
                r["_scaffold"] != row["_scaffold"] and r["split"] == "test" for r in refs
            ):
                continue
            if split_kind == "target":
                remaining = [r for r in refs if r["id"] != row["id"] and r["split"] == "test"]
                if any(
                    not any(r["_scaffold"] != scaffold for r in remaining) for scaffold in chosen_scaffolds
                ):
                    continue
                row["split"] = "query"
            queries.append({"id": f"{target}-query-{len(chosen_scaffolds) + 1}", "compound_id": row["id"]})
            chosen_scaffolds.add(row["_scaffold"])
            if len(chosen_scaffolds) >= queries_per_target:
                break
    for row in curated:
        del row["_family"], row["_scaffold"]
    split = {"kind": split_kind, "assignment_seed": seed}
    if split_kind == "time":
        split.update(cutoff_year=cutoff_year, date_kind="publication_year")
    panel = {
        "schema_version": 1,
        "dataset_id": f"chembl-human-binding-{split_kind}-{seed}",
        "scope": manifest.get("scope", "retrospective"),
        "sources": sources,
        "compounds": curated,
        "queries": queries,
        "split": split,
        "curation": {
            "require_direct_human_assay": require_direct_human_assay,
            "assay_assignment_provenance": assay_provenance,
            "active_nm": active_nm,
            "inactive_nm": inactive_nm,
            "max_heavy_atoms": max_heavy_atoms,
            "source_sampling": manifest.get("sampling", "See retained source manifest"),
            "source_manifest_sha256": digest(manifest_path),
            "date_semantics": "Document publication year, not assay date or proof of pretraining holdout",
            "assay_caveat": "Binding endpoints pooled by target; assay heterogeneity remains and is not removed by family-level concordance",
        },
    }
    write_json(output / "panel.json", panel)
    report = {
        "raw_qualifying_measurements": len(raw),
        "curated_families": len(curated),
        "queries": len(queries),
        "excluded_counts": dict(Counter(r["reason"] for r in excluded)),
        "excluded": excluded,
        "per_target": {
            target: dict(
                Counter((r["split"] + ":" + r["label"]) for r in curated if r["target_id"] == target)
            )
            for target in targets
        },
    }
    write_json(output / "curation.json", report)
    print(
        f"Validating {len(curated)} curated compounds and {len(queries)} queries", file=sys.stderr, flush=True
    )
    _, audit = load_dataset(output / "panel.json")
    write_json(output / "audit.json", audit)
    return panel, report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--split", choices=("scaffold", "target", "time"), default="scaffold")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--queries-per-target", type=int, default=2)
    parser.add_argument("--active-nm", type=float, default=1000)
    parser.add_argument("--inactive-nm", type=float, default=10000)
    parser.add_argument("--cutoff-year", type=int, default=2020)
    parser.add_argument("--max-heavy-atoms", type=int, default=60)
    parser.add_argument(
        "--allow-unverified-assay-targets",
        action="store_true",
        help="Diagnostic only: do not claim direct human assay assignment without assay metadata",
    )
    args = parser.parse_args(argv)
    if args.queries_per_target <= 0 or args.max_heavy_atoms <= 0:
        parser.error("Query and size caps must be positive")
    _, report = curate(
        args.manifest,
        args.output,
        args.split,
        args.seed,
        args.queries_per_target,
        args.active_nm,
        args.inactive_nm,
        args.cutoff_year,
        args.max_heavy_atoms,
        not args.allow_unverified_assay_targets,
    )
    print(json.dumps({k: v for k, v in report.items() if k != "excluded"}, indent=2))


if __name__ == "__main__":
    main()
