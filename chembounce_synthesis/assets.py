"""Convert external SynLlama assets to declarative ChemBounce input files.

Pickle index conversion substitutes inert placeholders for EVERY global. No
classes or functions from the pickle are imported or executed. Only the upstream
Molecule's stored SMILES string is retained; embeddings are discarded.
"""

import argparse
import csv
import hashlib
import pickle
import pickletools
from pathlib import Path

from rdkit import Chem

from .common import absolute_smiles, canonical, digest, write_json
from .reactions import Rules


class _Placeholder:
    module = ""
    name = ""

    def __new__(cls, *args, **kwargs):
        return object.__new__(cls)

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __setstate__(self, state):
        self.state = state


class _InertUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        # Never delegate to pickle.find_class or import a name requested by data.
        return type("Inert", (_Placeholder,), {"module": module, "name": name})

    def persistent_load(self, pid):
        raise ValueError("Persistent pickle references are not supported")


def read_index_smiles(path):
    """Extract only exact known structural fields from an inert SynLlama index."""
    with Path(path).open("rb") as stream:
        # EXT opcodes can bypass find_class through Python's extension cache.
        # Reject them before allowing the C unpickler to read the index.
        for op, _, _ in pickletools.genops(stream):
            if op.name in {"EXT1", "EXT2", "EXT4", "PERSID", "BINPERSID"}:
                raise ValueError(f"Unsupported pickle opcode: {op.name}")
        stream.seek(0)
        data = _InertUnpickler(stream).load()
    if not isinstance(data, _Placeholder) or (data.module, data.name) not in {
        ("synllama.chem.fpindex", "FingerprintIndex"),
        ("synllama.chem.matrix", "ReactantReactionMatrix"),
    }:
        raise ValueError("Expected a SynLlama FingerprintIndex or ReactantReactionMatrix")
    state = getattr(data, "state", None)
    if not isinstance(state, dict):
        raise ValueError("Unsupported index state")
    molecules = state.get("_molecules", state.get("_reactants"))
    if not isinstance(molecules, (list, tuple)):
        raise ValueError("Index contains no molecule sequence")
    for item in molecules:
        if (
            not isinstance(item, _Placeholder)
            or (item.module, item.name) != ("synllama.chem.mol", "Molecule")
            or not isinstance(getattr(item, "state", None), str)
        ):
            raise ValueError("Unsupported molecule state in index")
        yield item.state


def convert_catalog(inputs, output, source, snapshot, smiles_column="smiles", id_column="id"):
    """Convert SDF/CSV/TSV or inert SynLlama indices, recording skipped records.

    Supplied identifiers are stripped consistently with Catalog. Missing or blank
    IDs receive a structure-derived chemical identity, never a vendor SKU.
    """
    output = Path(output)
    if output.exists() or output.with_suffix(output.suffix + ".manifest.json").exists():
        raise ValueError("Refusing to overwrite converted catalog or manifest")
    if any(
        not isinstance(value, str) or not value.strip()
        for value in (source, snapshot, smiles_column, id_column)
    ):
        raise ValueError("Catalog source, snapshot identifier and column names must be nonempty strings")
    records, origins, failures, inputs_evidence = {}, {}, [], []
    for input_index, filename in enumerate(inputs):
        path = Path(filename)
        file_hash = digest(path)
        inputs_evidence.append({"path": str(path.resolve()), "sha256": file_hash})
        if path.suffix in {".pkl", ".pickle"}:
            entries = ((s, None) for s in read_index_smiles(path))
        elif path.suffix.lower() in {".smi", ".smiles", ".txt"}:
            entries = []
            for line in path.read_text().splitlines():
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                fields = line.split(maxsplit=1)
                if len(fields) == 2 and fields[1].startswith("|"):
                    end = fields[1].find("|", 1)
                    if end < 0:
                        entries.append((line.strip(), None))
                        continue
                    # Keep the complete annotation until the molecular-state check;
                    # it is chemistry metadata, never a supplied identifier.
                    suffix = fields[1][end + 1 :].strip()
                    if "|" in suffix:
                        # A second annotation must not hide inside the ID field.
                        entries.append((line.strip(), None))
                    else:
                        entries.append((fields[0] + " " + fields[1][: end + 1], suffix or None))
                else:
                    entries.append((fields[0], fields[1] if len(fields) == 2 else None))
        elif path.suffix.lower() == ".sdf":

            def sdf_entries():
                for mol in Chem.SDMolSupplier(str(path), removeHs=True):
                    yield (
                        mol if mol is not None else "",
                        mol.GetProp(id_column) if mol is not None and mol.HasProp(id_column) else None,
                    )

            entries = sdf_entries()
        else:
            # Materialize while the file is open; no implicit string splitting of salts.
            with path.open(newline="", encoding="utf-8-sig") as stream:
                reader = csv.DictReader(stream, delimiter="\t" if path.suffix == ".tsv" else ",")
                if smiles_column not in (reader.fieldnames or []):
                    raise ValueError(f"Missing SMILES column {smiles_column}: {path}")
                entries = [(r[smiles_column], r.get(id_column)) for r in reader]
        for index, (smiles, ident) in enumerate(entries):
            try:
                smiles = canonical(absolute_smiles(smiles) if isinstance(smiles, Chem.Mol) else smiles)
            except ValueError as exc:
                failures.append(
                    {"input_index": input_index, "file": str(path), "record": index, "reason": str(exc)}
                )
                continue
            ident = ident.strip() if ident else ""
            identifier_origin = "input_provided" if ident else "generated_structure_sha256"
            if not ident:
                ident = "structure-sha256:" + hashlib.sha256(smiles.encode()).hexdigest()
            if ident in records and records[ident]["smiles"] != smiles:
                raise ValueError(f"Conflicting structures for catalog ID {ident}")
            records[ident] = {
                "id": ident,
                "smiles": smiles,
                "source": source,
                "snapshot": snapshot,
                "availability": "catalogued",
            }
            origins.setdefault(ident, []).append(
                {
                    "input_index": input_index,
                    "record_index": index,
                    "identifier_origin": identifier_origin,
                }
            )
        if digest(path) != file_hash:
            raise ValueError(f"Catalog input changed during conversion: {path}")
    if not records:
        raise ValueError("No valid connected building blocks in inputs")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["id", "smiles", "source", "snapshot", "availability"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(records[k] for k in sorted(records))
    manifest = {
        "schema_version": 2,
        "inputs": inputs_evidence,
        "output_sha256": digest(output),
        "records": len(records),
        "rejected_records": failures,
        "source": source,
        "snapshot": snapshot,
        "smiles_column": smiles_column,
        "id_column": id_column,
        "record_index_policy": "zero_based_logical_records",
        "record_origins": {ident: origins[ident] for ident in sorted(origins)},
        "id_policy": "Strip input ID when nonblank, otherwise structure-sha256 (not a verified supplier SKU)",
        "availability_policy": "Snapshot membership only; live stock is not checked",
    }
    write_json(output.with_suffix(output.suffix + ".manifest.json"), manifest)
    return manifest


def convert_rules(filename, output, source):
    path, output = Path(filename), Path(output)
    if output.exists():
        raise ValueError("Refusing to overwrite reaction rules")
    if not source:
        raise ValueError("Reaction source is required")
    rows, seen, duplicates = [], {}, []
    for line, text in enumerate(path.read_text().splitlines(), 1):
        if text.strip() and not text.lstrip().startswith("#"):
            if text.strip() in seen:
                duplicates.append({"source_line": line, "same_as_rule": seen[text.strip()]})
                continue
            seen[text.strip()] = f"synllama-rule-{line:03d}"
            rows.append({"id": f"synllama-rule-{line:03d}", "smarts": text.strip(), "source": source})
    data = {
        "schema_version": 1,
        "source_sha256": digest(path),
        "source_path": str(path.resolve()),
        "duplicate_source_templates": duplicates,
        "rules": rows,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    # Validate before publishing a usable registry.
    temp = output.with_suffix(output.suffix + ".tmp")
    try:
        write_json(temp, data)
        Rules(temp)
        temp.replace(output)
    finally:
        temp.unlink(missing_ok=True)
    return data


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    catalog = sub.add_parser("catalog")
    catalog.add_argument("--input", action="append", required=True)
    catalog.add_argument("--output", required=True)
    catalog.add_argument("--source", required=True)
    catalog.add_argument("--snapshot", required=True)
    catalog.add_argument("--smiles-column", default="smiles")
    catalog.add_argument("--id-column", default="id")
    rules = sub.add_parser("rules")
    rules.add_argument("--input", required=True)
    rules.add_argument("--output", required=True)
    rules.add_argument("--source", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "catalog":
            result = convert_catalog(
                args.input, args.output, args.source, args.snapshot, args.smiles_column, args.id_column
            )
            print(f"Catalog entries: {result['records']}; rejected: {len(result['rejected_records'])}")
        else:
            result = convert_rules(args.input, args.output, args.source)
            print(f"Registered reaction rules: {len(result['rules'])}")
    except (ValueError, OSError, pickle.UnpicklingError, TypeError) as exc:
        parser.exit(2, f"Asset conversion failed: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
