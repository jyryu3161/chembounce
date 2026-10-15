"""Dataset provenance, chemical partitions, and model-overlap audits.

Benchmark development partitions and a pretrained model's training inventory are
separate facts. An unavailable/partial model inventory can never prove no leakage.
"""

import csv
from datetime import date
import hashlib
import json
import random
from pathlib import Path

from chembounce_synthesis.identity import chemical_identity


LABELS = {"active", "inactive", "decoy", "unlabelled"}
PARTITIONS = {"train", "validation", "test", "query"}


class ValidatedDataset(dict):
    """Plain JSON-compatible data plus per-load caches, never serialized labels in requests."""

    def __init__(self, data, identities):
        super().__init__(data)
        self.identities = identities
        self.compounds_by_id = {r["id"]: r for r in self["compounds"]}
        self.reference_cache = {}


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def require_text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def verified_artifact(record, base):
    path = (Path(base) / require_text(record.get("artifact"), "artifact")).resolve()
    expected = record.get("sha256")
    if not isinstance(expected, str) or len(expected) != 64:
        raise ValueError(f"Missing SHA-256 for {path}")
    actual = digest(path)
    if actual != expected:
        raise ValueError(f"SHA-256 mismatch for {path}")
    return path


def load_dataset(path):
    """Validate sources, labels and declared development partitions before running.

    Query compounds are explicit public inputs. A query may be drawn from a
    development partition; its scaffold is excluded from recovery ground truth.
    The separate ``query`` partition supports target-held-out public inputs.
    Activity is only asserted for records with an assay identifier and source;
    decoys are a distinct label and are never treated as experimentally inactive.
    """
    path = Path(path).resolve()
    data = json.loads(path.read_text())
    if data.get("schema_version") != 1:
        raise ValueError("Research dataset requires schema_version=1")
    require_text(data.get("dataset_id"), "dataset_id")
    if data.get("scope") not in {"retrospective", "synthetic_fixture"}:
        raise ValueError("scope must be retrospective or synthetic_fixture")
    sources = {}
    for source in data.get("sources", []):
        ident = require_text(source.get("id"), "source id")
        if ident in sources:
            raise ValueError("Duplicate source id")
        require_text(source.get("url"), "source URL")
        sources[ident] = str(verified_artifact(source, path.parent))
    if not sources:
        raise ValueError("At least one retained, hashed source artifact is required")
    records, identities = {}, {}
    for record in data.get("compounds", []):
        ident = require_text(record.get("id"), "compound id")
        if ident in records:
            raise ValueError("Duplicate compound id")
        if record.get("source_id") not in sources:
            raise ValueError(f"Unknown source for {ident}")
        require_text(record.get("target_id"), "target_id")
        if record.get("label") not in LABELS or record.get("split") not in PARTITIONS:
            raise ValueError(f"Invalid label or partition for {ident}")
        if record["label"] in {"active", "inactive"}:
            require_text(record.get("assay_id"), f"assay_id for {ident}")
            require_text(record.get("label_basis"), f"label_basis for {ident}")
        if record.get("publication_year") is not None:
            year = record["publication_year"]
            if isinstance(year, bool) or not isinstance(year, int) or not 1800 <= year <= 2200:
                raise ValueError("publication_year must be an integer year")
        if record.get("date") is not None:
            date.fromisoformat(record["date"])
            require_text(record.get("date_kind"), f"date_kind for {ident}")
        identities[ident] = chemical_identity(record["smiles"])
        records[ident] = record
    if not records:
        raise ValueError("No compounds")
    queries, query_ids = [], set()
    for query in data.get("queries", []):
        ident = require_text(query.get("id"), "query id")
        if ident in query_ids or query.get("compound_id") not in records:
            raise ValueError("Duplicate query id or missing query compound")
        query_ids.add(ident)
        if query.get("reference_sdf") is not None:
            query["resolved_reference_sdf"] = str(verified_artifact(query["reference_sdf"], path.parent))
        queries.append(query)
    if not queries:
        raise ValueError("No queries")
    # Conflicting assay labels must be represented as distinct endpoints/panels,
    # not silently collapsed into a single target-level ground truth.
    labels = {}
    for ident, record in records.items():
        key = (record["target_id"], identities[ident].family_smiles)
        if record["label"] in {"active", "inactive"}:
            if key in labels and labels[key] != record["label"]:
                raise ValueError("Conflicting active/inactive labels within a target chemical family")
            labels[key] = record["label"]
    split = data.get("split", {})
    if split.get("kind") not in {"scaffold", "target", "time"}:
        raise ValueError("Declare split.kind: scaffold, target, or time")
    grouped = [r for r in records.values() if r["split"] != "query"]
    for key_name in ("family", "split_group"):
        groups = {}
        for record in grouped:
            identity = identities[record["id"]]
            if key_name == "family":
                key = identity.family_smiles
            elif split["kind"] == "scaffold":
                key = identity.scaffold_family_smiles or "acyclic:" + identity.family_smiles
            elif split["kind"] == "target":
                key = record["target_id"]
            else:
                continue
            if key in groups and groups[key] != record["split"]:
                raise ValueError(f"Cross-partition {key_name} leakage: {key}")
            groups[key] = record["split"]
    if split["kind"] == "time":
        year_precision = "cutoff_year" in split
        cutoff = split["cutoff_year"] if year_precision else date.fromisoformat(split["cutoff"])
        kind = require_text(split.get("date_kind"), "split date_kind")
        if year_precision and (
            kind != "publication_year" or not isinstance(cutoff, int) or isinstance(cutoff, bool)
        ):
            raise ValueError("cutoff_year requires date_kind=publication_year and integer cutoff")
        for record in grouped:
            if year_precision:
                if record.get("publication_year") is None:
                    raise ValueError("Time split needs sourced publication years")
                value = record["publication_year"]
            else:
                if record.get("date") is None or record.get("date_kind") != kind:
                    raise ValueError("Time split needs sourced dates of one explicit date_kind")
                value = date.fromisoformat(record["date"])
            if record["split"] == "train" and value > cutoff:
                raise ValueError("Training record occurs after time cutoff")
            if record["split"] in {"validation", "test"} and value <= cutoff:
                raise ValueError("Held-out record occurs on/before time cutoff")
    coverage = []
    for query in queries:
        record = records[query["compound_id"]]
        identity = identities[record["id"]]
        references = [
            r
            for r in records.values()
            if r["target_id"] == record["target_id"]
            and r["split"] == "test"
            and r["label"] == "active"
            and identities[r["id"]].scaffold_family_smiles != identity.scaffold_family_smiles
            and identities[r["id"]].has_scaffold
        ]
        coverage.append({"query_id": query["id"], "novel_active_references": len(references)})
    return ValidatedDataset(data, identities), {
        "dataset_sha256": digest(path),
        "source_artifacts": sources,
        "chemical_identity_policy": next(iter(identities.values())).policy_version,
        "split": split,
        "query_coverage": coverage,
        "scientific_evidence": data["scope"] != "synthetic_fixture",
    }


def read_inventory(path):
    """Explicit JSON list, TSV smiles column, or one SMILES per text line."""
    path = Path(path)
    if path.suffix == ".json":
        values = json.loads(path.read_text())
        if not isinstance(values, list):
            raise ValueError("Inventory JSON must be a list of SMILES or objects with smiles")
        return [v.get("smiles") if isinstance(v, dict) else v for v in values]
    if path.suffix == ".tsv":
        with path.open(newline="") as stream:
            return [r["smiles"] for r in csv.DictReader(stream, delimiter="\t")]
    return [line.strip().split()[0] for line in path.read_text().splitlines() if line.strip()]


def audit_inventory(spec, dataset, base):
    """Unknown is an explicit outcome, not a successful overlap check."""
    kind = spec.get("kind")
    if kind not in {"training", "prompt"}:
        raise ValueError("Inventory kind must be training or prompt")
    status = spec.get("status", "unknown")
    if status == "unknown":
        return {"kind": kind, "status": "unknown", "reason": require_text(spec.get("reason"), "reason")}
    if status != "available" or spec.get("coverage") not in {"complete", "partial"}:
        raise ValueError("Available inventory requires complete/partial coverage")
    path = verified_artifact(spec, base)
    entries = list(enumerate(read_inventory(path)))
    full_size = len(entries)
    limit = spec.get("max_records")
    coverage = spec["coverage"]
    if limit is not None:
        if isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0:
            raise ValueError("Inventory max_records must be a positive integer")
        if len(entries) > limit:
            entries = random.Random(42).sample(entries, limit)
            coverage = "partial"
    identity, unassessed = [], []
    for index, smiles in entries:
        try:
            identity.append(chemical_identity(smiles))
        except (ValueError, TypeError) as exc:
            # Search libraries can include attachment-point fragments and valid
            # molecules outside the bounded normalization policy. Neither these
            # nor malformed inventory records establish a clean holdout.
            unassessed.append(
                {
                    "record_index": index,
                    "input_sha256": hashlib.sha256(
                        json.dumps(smiles, ensure_ascii=False, sort_keys=True).encode()
                    ).hexdigest(),
                    "input": smiles,
                    "reason": str(exc),
                }
            )
    if unassessed:
        coverage = "partial"
    families = {i.family_smiles for i in identity}
    scaffolds = {i.scaffold_family_smiles for i in identity if i.has_scaffold}
    matches = []
    for record in dataset["compounds"]:
        ident = (
            dataset.identities[record["id"]]
            if isinstance(dataset, ValidatedDataset)
            else chemical_identity(record["smiles"])
        )
        family = ident.family_smiles in families
        scaffold = ident.has_scaffold and ident.scaffold_family_smiles in scaffolds
        if family or scaffold:
            matches.append(
                {"compound_id": record["id"], "family_overlap": family, "scaffold_overlap": scaffold}
            )
    return {
        "kind": kind,
        "status": "overlap_detected"
        if matches
        else ("no_overlap_in_assessed_subset" if unassessed else "no_overlap_in_declared_inventory"),
        "coverage": coverage,
        "declared_coverage": spec["coverage"],
        "analysis_coverage": "partial" if unassessed else "complete_for_selected_records",
        "declared_inventory_records": full_size,
        "sampling": "seed-42 uniform sample without replacement"
        if len(entries) < full_size
        else "complete declared file",
        "sha256": digest(path),
        "selected_records": len(entries),
        "records": len(identity),
        "unassessed_records": len(unassessed),
        "unassessed": unassessed,
        "record_index_convention": "zero-based index in parsed inventory; input_sha256 hashes JSON encoding of input",
        "matches": matches,
        "complete_inventory_no_detected_overlap": not matches and coverage == "complete",
        "interpretation": "Overlap check for this declared inventory only; absence of overlap does not establish model generalization or unseen pretraining",
    }


def synthesis_asset_paths(config_path):
    """External files referenced indirectly by a synthesis configuration."""
    from chembounce_synthesis.workflow import load_config

    config = load_config(config_path)
    result = {}
    for name in ("model_path", "prompt_path", "path"):
        if config["provider"].get(name):
            result["pipeline_provider_" + name] = str(Path(config["provider"][name]).resolve())
    predictor = config.get("synthesis_evidence", {}).get("predictor")
    if predictor:
        result["pipeline_evidence_model"] = str(Path(predictor["model_path"]).resolve())
    if config.get("reference_sdf"):
        result["pipeline_reference_sdf"] = str(Path(config["reference_sdf"]).resolve())
    for index, carry in enumerate(config.get("carry_forward", [])):
        result[f"pipeline_carry_{index}"] = str(Path(carry).resolve())
    return result
