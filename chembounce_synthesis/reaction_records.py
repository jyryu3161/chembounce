"""Optional offline ORD records. Report-only: never synthesis acceptance evidence.

Hashes bind retained bytes to a caller's receipt, not to experimental truth.
Every new match replays the dataset coordinates and derives chemistry from raw
protobuf. PyArrow/protobuf are lazy optional dependencies, never default imports.
"""

import base64
import copy
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import stat
import sys
import subprocess
import tempfile

VERSION = "public-reaction-records-v1"
SCHEMA_REVISION = "51b3a4968835f0090268cab9c0be034b48ac62d4"
DESCRIPTOR_SHA256 = "386431d9cb8ffeaabeaeea0b0ab2369221394466eeb81d2c0b374df95c37aeb9"
RESOURCE = Path(__file__).parent / "data" / "reaction_records"
MAX_FILE = 16 * 1024 * 1024
MAX_DATASET = 2 * 1024 * 1024
MAX_DECOMPRESSED = 16 * 1024 * 1024
MAX_ROWS = 4096
MAX_RECORDS = 64
MAX_RECORD = 256 * 1024
MAX_DEPTH = 64
LIMITS = {
    "dataset_bytes": MAX_DATASET,
    "decompressed_bytes": MAX_DECOMPRESSED,
    "dataset_rows": MAX_ROWS,
    "selected_records": MAX_RECORDS,
    "record_bytes": MAX_RECORD,
    "parser_extra_address_space_bytes": 256 * 1024 * 1024,
    "parser_cpu_seconds": 25,
    "parser_wall_seconds": 35,
    "parser_output_bytes": MAX_FILE,
}
INTERPRETATION = (
    "Reported public source records only. Byte integrity is relative to declared/external receipts; "
    "it does not authenticate acquisition, experimenter, structural characterization or laboratory truth. "
    "No synthesis success/failure, isolated-product guarantee, training independence or coverage claim. "
    "No match means no_match_in_retained_snapshot. Computational acceptance is unchanged."
)


class RecordError(ValueError):
    def __init__(self, code, message):
        self.code = code
        super().__init__(f"{code}: {message}")


def require(condition, code, message):
    if not condition:
        raise RecordError(code, message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def binding(value):
    return sha(encoded(value))


def _object(pairs):
    value = {}
    for key, child in pairs:
        require(key not in value, "duplicate_json_key", key)
        value[key] = child
    return value


def json_read(raw):
    try:
        value = json.loads(
            raw, object_pairs_hook=_object, parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x))
        )
        encoded(value)  # finite JSON only
        return value
    except (ValueError, TypeError, RecursionError, UnicodeError) as exc:
        raise RecordError("invalid_json", str(exc)) from exc


def read_local(root, name, limit=MAX_FILE):
    """A single regular contained file, no path components or symlink following."""
    require(
        isinstance(name, str)
        and name not in {"", ".", ".."}
        and Path(name).name == name
        and "/" not in name
        and "\\" not in name,
        "unsafe_path",
        "Expected a local basename",
    )
    root = Path(root).resolve(strict=True)
    path = root / name
    try:
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        with os.fdopen(fd, "rb") as stream:
            before = os.fstat(stream.fileno())
            require(
                stat.S_ISREG(before.st_mode) and before.st_size <= limit,
                "file_bound",
                f"{name} is not a regular file within {limit} bytes",
            )
            raw = stream.read(limit + 1)
            after = os.fstat(stream.fileno())
        stamp = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
        require(
            len(raw) <= limit and stamp(before) == stamp(after) == stamp(path.lstat()), "file_changed", name
        )
        return raw
    except OSError as exc:
        raise RecordError("unsafe_or_missing_file", str(exc)) from exc


def _dependencies():
    try:
        import pyarrow
        import pyarrow.parquet as pq
        import google.protobuf
        from google.protobuf import descriptor_pool, json_format, message_factory
    except ImportError as exc:
        raise RecordError(
            "optional_dependency_missing",
            "ORD import/replay/matching requires the optional reaction-records extra "
            "(PyArrow and protobuf) in a separate prepared environment; no auto-install",
        ) from exc
    require(
        hasattr(message_factory, "GetMessageClass"), "unsupported_runtime", "protobuf GetMessageClass needed"
    )
    return pyarrow, pq, google.protobuf, descriptor_pool, json_format, message_factory


def _runtime():
    arrow, _, proto, *_ = _dependencies()
    from rdkit import rdBase

    return {
        "executable": sys.executable,
        "executable_realpath": os.path.realpath(sys.executable),
        "python": sys.version,
        "pyarrow": arrow.__version__,
        "pyarrow_module": arrow.__file__,
        "protobuf": proto.__version__,
        "protobuf_module": proto.__file__,
        "rdkit": rdBase.rdkitVersion,
        "extractor": VERSION,
        "extractor_source_sha256": sha(Path(__file__).read_bytes()),
        "worker_source_sha256": sha(Path(__file__).with_name("reaction_record_worker.py").read_bytes()),
        "parser_platform": "Linux RLIMIT_AS/CPU/FSIZE + wall timeout; limits apply only in child",
        "schema_revision": SCHEMA_REVISION,
        "validation": "bounded descriptor parsing; not full ORD semantic validation or experimental authentication",
    }


def _integer(value):
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _check_snapshot(snapshot):
    require(
        isinstance(snapshot, dict)
        and set(snapshot)
        == {
            "schema",
            "provider",
            "data_revision",
            "schema_revision",
            "dataset",
            "descriptor",
            "selectors",
            "source_metadata",
        },
        "snapshot_schema",
        "Unexpected snapshot fields",
    )
    require(
        snapshot["schema"] == "ord-retained-snapshot-v1"
        and snapshot["provider"] == "open-reaction-database"
        and snapshot["schema_revision"] == SCHEMA_REVISION,
        "unsupported_source_schema",
        "Unsupported source/schema version",
    )
    require(
        isinstance(snapshot["data_revision"], str)
        and re.fullmatch(r"[0-9a-f]{40}", snapshot["data_revision"]),
        "source_revision",
        "Immutable data revision required",
    )
    require(
        isinstance(snapshot["dataset"], dict)
        and set(snapshot["dataset"]) == {"path", "sha256", "bytes", "dataset_id", "url"},
        "snapshot_schema",
        "Dataset fields differ",
    )
    require(
        isinstance(snapshot["descriptor"], dict) and set(snapshot["descriptor"]) == {"path", "sha256"},
        "snapshot_schema",
        "Descriptor fields differ",
    )
    require(
        snapshot["descriptor"]["sha256"] == DESCRIPTOR_SHA256,
        "unsupported_schema",
        "Only the audited descriptor is supported",
    )
    require(
        isinstance(snapshot["source_metadata"], dict)
        and snapshot["source_metadata"].get("license") == "CC-BY-SA-4.0",
        "source_license",
        "ORD data license/attribution metadata required",
    )
    selectors = snapshot["selectors"]
    require(
        isinstance(selectors, list) and 0 < len(selectors) <= MAX_RECORDS,
        "selector_bound",
        "Explicit bounded selectors required",
    )
    seen = set()
    for selected in selectors:
        require(
            isinstance(selected, dict) and set(selected) == {"row_group", "row", "reaction_id"},
            "selector",
            "Expected row_group/row/reaction_id",
        )
        require(
            _integer(selected["row_group"])
            and _integer(selected["row"])
            and isinstance(selected["reaction_id"], str)
            and re.fullmatch(r"ord-[0-9a-f]{32}", selected["reaction_id"]),
            "selector",
            "Invalid selector",
        )
        key = (selected["row_group"], selected["row"])
        require(key not in seen, "duplicate_selector", str(key))
        seen.add(key)


def _wire_metadata(message):
    """Retain presence semantics and actual float32/64 values, not JSON precision guesses."""
    presence, floats = {}, []

    def visit(msg, path, depth):
        require(depth <= MAX_DEPTH, "protobuf_depth", "Message nesting limit")
        populated = {field.name for field, _ in msg.ListFields()}
        for field in msg.DESCRIPTOR.fields:
            here = f"{path}/{field.name}"
            repeated = getattr(field, "is_repeated", None)
            if repeated is None:  # protobuf 6 descriptor API
                repeated = field.label == field.LABEL_REPEATED
            if repeated:
                presence[here] = {
                    "semantics": "repeated_or_map",
                    "state": "nonempty" if field.name in populated else "empty_or_omitted_indistinguishable",
                }
            elif field.has_presence:
                presence[here] = {
                    "semantics": "explicit",
                    "state": "present" if msg.HasField(field.name) else "absent",
                }
            else:
                presence[here] = {
                    "semantics": "implicit_scalar",
                    "state": "nondefault"
                    if field.name in populated
                    else "default_or_omitted_indistinguishable",
                }
            if field.name not in populated:
                continue
            value = getattr(msg, field.name)
            if field.enum_type:
                enum_values = value if repeated else [value]
                require(all(v in field.enum_type.values_by_number for v in enum_values), "unknown_enum", here)
            if field.message_type and field.message_type.GetOptions().map_entry:
                value_field = field.message_type.fields_by_name["value"]
                for key in sorted(value):
                    if value_field.message_type:
                        visit(value[key], f"{here}/{key}", depth + 1)
                continue
            values = list(enumerate(value)) if repeated else [(None, value)]
            for index, item in values:
                child = here if index is None else f"{here}/{index}"
                if field.message_type:
                    visit(item, child, depth + 1)
                elif field.type in {field.TYPE_FLOAT, field.TYPE_DOUBLE}:
                    require(math.isfinite(item), "nonfinite_record", child)
                    floats.append(
                        {
                            "path": child,
                            "wire_type": "float32" if field.type == field.TYPE_FLOAT else "float64",
                            "decoded_value": item,
                            "interpretation": "wire representation, not measurement precision",
                        }
                    )

    visit(message, "", 0)
    return presence, floats


def identity(identifiers):
    """Exact full-material identity; missing/ambiguous stereo and conflicting IDs abstain."""
    from rdkit import Chem
    from .common import parse_absolute_smiles, require_absolute_stereo, absolute_smiles

    keys, issues, unverified = [], [], []
    for item in identifiers:
        kind = item.get("type")
        descriptive = {
            "IUPAC_NAME",
            "NAME",
            "CAS_NUMBER",
            "PUBCHEM_CID",
            "CHEMSPIDER_ID",
            "UNIPROT_ID",
            "PDB_ID",
            "MDL",
        }
        if kind in descriptive or (kind == "CUSTOM" and item.get("details") == "Molecular formula"):
            unverified.append(kind)
            continue
        if kind not in {"SMILES", "CXSMILES", "MOLBLOCK"}:
            issues.append(f"unsupported_identity_bearing_identifier:{kind}")
            continue
        try:
            value = item["value"]
            require(
                isinstance(value, str) and len(value) <= MAX_RECORD, "structure_bound", "Identifier too large"
            )
            mol = (
                parse_absolute_smiles(value)
                if kind in {"SMILES", "CXSMILES"}
                else Chem.MolFromMolBlock(value, sanitize=True, removeHs=True, strictParsing=True)
            )
            require(mol is not None and mol.GetNumAtoms() > 0, "invalid_structure", kind)
            require_absolute_stereo(mol)
            require(not any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()), "dummy_atoms", kind)
            potential = Chem.FindPotentialStereo(mol)
            require(not any(str(s.specified) != "Specified" for s in potential), "unspecified_stereo", kind)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
            keys.append(absolute_smiles(mol, canonical=True, isomericSmiles=True))
        except (ValueError, RuntimeError, KeyError, TypeError) as exc:
            issues.append(str(exc))
    if not keys:
        issues.append("no_assessable_structural_identifier")
    if len(set(keys)) > 1:
        issues.append("conflicting_structural_identifiers")
    return {
        "status": "unassessable" if issues else "exact_absolute_state",
        "key": None if issues else keys[0],
        "reasons": issues,
        "comparison_scope": "All supported structural identifiers; descriptive names/accessions are retained but unverified",
        "unverified_descriptive_identifier_types": unverified,
    }


def smiles_identity(smiles):
    return identity([{"type": "SMILES", "value": smiles}])


def _quantitative_amount(amount):
    values = [amount[k] for k in ("mass", "moles", "volume") if k in amount]
    return (
        len(values) == 1
        and "value" in values[0]
        and values[0]["value"] >= 0
        and bool(values[0].get("units"))
        and values[0].get("precision", 0) >= 0
    )


def _measurement_reasons(measurement):
    kind = measurement.get("type", "UNSPECIFIED")
    reasons = []
    expected = {
        "YIELD": "percentage",
        "PURITY": "percentage",
        "AREA": "float_value",
        "COUNTS": "float_value",
        "INTENSITY": "float_value",
        "AMOUNT": "amount",
    }
    if kind in {"CUSTOM", "UNSPECIFIED"}:
        reasons.append("unsupported_measurement_type")
    if kind != "IDENTITY" and not any(
        k in measurement for k in ("percentage", "float_value", "string_value", "amount")
    ):
        reasons.append("measurement_value_missing")
    if "amount" in measurement and not _quantitative_amount(measurement["amount"]):
        reasons.append("invalid_or_unmeasured_product_amount")
    if "retention_time" in measurement:
        time = measurement["retention_time"]
        if "value" not in time or time["value"] < 0 or not time.get("units") or time.get("precision", 0) < 0:
            reasons.append("invalid_retention_time")
    value_key = expected.get(kind)
    if value_key and value_key not in measurement:
        reasons.append("measurement_value_type_mismatch_or_missing")
    for key in ("percentage", "float_value"):
        if key in measurement:
            value = measurement[key]
            if "value" not in value or value["value"] < 0:
                reasons.append("measurement_value_missing_or_negative")
            if key == "percentage" and value.get("value", 0) > 100:
                reasons.append("percentage_out_of_range")
            if value.get("precision", 0) < 0:
                reasons.append("negative_measurement_precision")
    return reasons


def _derive(raw, reaction_cls, json_format):
    require(0 < len(raw) <= MAX_RECORD, "record_bound", "Reaction bytes exceed limit")
    try:
        message = reaction_cls.FromString(raw)
        original = message.SerializeToString(deterministic=True)
        known = reaction_cls.FromString(raw)
        known.DiscardUnknownFields()
        require(
            original == known.SerializeToString(deterministic=True),
            "unknown_protobuf_fields",
            "Cannot safely interpret unknown schema fields",
        )
        projection = json_format.MessageToDict(message, preserving_proto_field_name=True)
        presence, floats = _wire_metadata(message)
    except RecordError:
        raise
    except Exception as exc:
        raise RecordError("protobuf_parse", str(exc)) from exc
    inputs, products, rejected = [], [], []
    for group, data in sorted(projection.get("inputs", {}).items()):
        for index, component in enumerate(data.get("components", [])):
            inputs.append(
                {
                    "path": f"/inputs/{group}/components/{index}",
                    "role": component.get("reaction_role", "UNSPECIFIED"),
                    "identity": identity(component.get("identifiers", [])),
                    "component": component,
                }
            )
    for oi, outcome in enumerate(projection.get("outcomes", [])):
        for pi, product in enumerate(outcome.get("products", [])):
            for measurement in product.get("measurements", []):
                for reason in _measurement_reasons(measurement):
                    rejected.append({"code": reason, "path": f"/outcomes/{oi}/products/{pi}"})
                key = measurement.get("analysis_key")
                if not key or key not in outcome.get("analyses", {}):
                    rejected.append(
                        {
                            "code": "unresolved_analysis_reference",
                            "path": f"/outcomes/{oi}/products/{pi}",
                            "analysis_key": key,
                        }
                    )
            products.append(
                {
                    "path": f"/outcomes/{oi}/products/{pi}",
                    "outcome_index": oi,
                    "product_index": pi,
                    "role": product.get("reaction_role", "UNSPECIFIED"),
                    "identity": identity(product.get("identifiers", [])),
                    "product": product,
                }
            )
    return {
        "reaction_id": projection.get("reaction_id"),
        "projection": projection,
        "presence": presence,
        "wire_floats": floats,
        "inputs": inputs,
        "products": products,
        "rejections": rejected,
        "projection_format": "protobuf JSON (float32 short decimals); see wire_floats and raw protobuf",
        "presence_limit": "Implicit scalars and empty repeated/maps cannot distinguish omitted from explicit defaults",
    }


def _replay_local(snapshot, dataset_raw, descriptor_raw):
    _check_snapshot(snapshot)
    dataset = snapshot["dataset"]
    require(
        len(dataset_raw) <= MAX_DATASET
        and len(dataset_raw) == dataset["bytes"]
        and sha(dataset_raw) == dataset["sha256"],
        "dataset_hash_or_size",
        "Dataset bytes differ",
    )
    require(sha(descriptor_raw) == DESCRIPTOR_SHA256, "schema_hash", "Descriptor differs")
    arrow, pq, _, descriptor_pool, json_format, message_factory = _dependencies()
    pool = descriptor_pool.DescriptorPool()
    pool.AddSerializedFile(descriptor_raw)
    reaction_cls = message_factory.GetMessageClass(pool.FindMessageTypeByName("ord.Reaction"))
    try:
        file = pq.ParquetFile(
            io.BytesIO(dataset_raw),
            thrift_string_size_limit=MAX_FILE,
            thrift_container_size_limit=MAX_ROWS * 16,
        )
        require(
            file.metadata.num_rows <= MAX_ROWS and file.num_row_groups <= MAX_ROWS,
            "row_bound",
            "Dataset rows/groups exceed cap",
        )
        require(
            file.schema_arrow.names == ["reaction_id", "reaction"]
            and file.schema_arrow.field("reaction_id").type == arrow.string()
            and file.schema_arrow.field("reaction").type == arrow.binary(),
            "parquet_schema",
            "Expected ORD reaction_id/string and reaction/binary columns",
        )
        metadata = {
            k.decode(): v.decode()
            for k, v in (file.schema_arrow.metadata or {}).items()
            if k.startswith(b"ord.")
        }
        require(
            metadata.get("ord.schema_version") == "1"
            and metadata.get("ord.dataset_id") == dataset["dataset_id"],
            "dataset_metadata",
            "Dataset ID or storage version differs",
        )
        decompressed = sum(
            file.metadata.row_group(g).column(c).total_uncompressed_size
            for g in range(file.num_row_groups)
            for c in range(2)
        )
        require(0 <= decompressed <= MAX_DECOMPRESSED, "decompressed_bound", "Parquet decompression cap")
        require(
            all(
                0 <= file.metadata.row_group(g).total_byte_size <= MAX_DECOMPRESSED
                for g in range(file.num_row_groups)
            ),
            "decompressed_bound",
            "Row-group decompression cap",
        )
        groups = {}
        for selected in snapshot["selectors"]:
            group, row = selected["row_group"], selected["row"]
            require(
                group < file.num_row_groups and row < file.metadata.row_group(group).num_rows,
                "selector_range",
                "Row coordinates out of bounds",
            )
            groups.setdefault(group, {})[row] = selected
        records, raw_records, seen, actual_bytes = [], [], {}, 0
        for group, selection in sorted(groups.items()):
            # Only selected row groups decode; small batches bound returned Arrow buffers.
            index = 0
            for batch in file.iter_batches(batch_size=1, row_groups=[group], use_threads=False):
                actual_bytes += batch.nbytes
                require(
                    actual_bytes <= MAX_DECOMPRESSED, "decompressed_bound", "Decoded Arrow bytes exceed cap"
                )
                if index in selection:
                    row = batch.to_pylist()[0]
                    raw = row["reaction"]
                    selected = selection[index]
                    require(
                        isinstance(raw, bytes) and row["reaction_id"] == selected["reaction_id"],
                        "record_coordinate_id",
                        "Row ID does not match selector",
                    )
                    digest = sha(raw)
                    require(
                        row["reaction_id"] not in seen or seen[row["reaction_id"]] == digest,
                        "conflicting_record_id",
                        "Same ID has different raw bytes",
                    )
                    seen[row["reaction_id"]] = digest
                    derived = _derive(raw, reaction_cls, json_format)
                    require(
                        derived["reaction_id"] == selected["reaction_id"],
                        "record_payload_id",
                        "Protobuf ID differs from row/selector",
                    )
                    records.append({"selector": selected, "raw_sha256": digest, "derived": derived})
                    raw_records.append(raw)
                index += 1
                if index > max(selection):
                    break
        require(len(records) == len(snapshot["selectors"]), "selector_coverage", "Missing selected records")
        return records, raw_records, metadata
    except RecordError:
        raise
    except Exception as exc:
        raise RecordError("parquet_replay", str(exc)) from exc


def _replay(snapshot, dataset_raw, descriptor_raw):
    """Native parsing in a bounded child; never change caller resource limits."""
    _dependencies()  # clear optional-extra error before spawning
    request = encoded(
        {
            "snapshot": snapshot,
            "dataset": base64.b64encode(dataset_raw).decode(),
            "descriptor": base64.b64encode(descriptor_raw).decode(),
        }
    )
    require(len(request) <= MAX_FILE, "request_bound", "Worker request exceeds cap")
    env = dict(
        os.environ,
        PYTHONDONTWRITEBYTECODE="1",
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        CUDA_VISIBLE_DEVICES="",
    )
    # A file path avoids cwd/PYTHONPATH importing a different package. Worker imports this checkout.
    try:
        with tempfile.TemporaryFile() as output, tempfile.TemporaryFile() as errors:
            completed = subprocess.run(
                [sys.executable, str(Path(__file__).with_name("reaction_record_worker.py"))],
                input=request,
                stdout=output,
                stderr=errors,
                env=env,
                timeout=35,
            )
            output.seek(0)
            raw = output.read(MAX_FILE + 1)
            require(
                completed.returncode == 0,
                "parser_process_limit_or_failure",
                f"Bounded parser exited {completed.returncode}; no records accepted",
            )
            require(len(raw) <= MAX_FILE, "worker_output_bound", "Parser output exceeds cap")
        result = json_read(raw)
        if "error" in result:
            raise RecordError(result["error"]["code"], result["error"]["message"])
        return (
            result["records"],
            [base64.b64decode(v, validate=True) for v in result["raw_records"]],
            result["metadata"],
            result["parser_observation"],
        )
    except subprocess.TimeoutExpired as exc:
        raise RecordError("parser_timeout", "Bounded parser exceeded 35 seconds") from exc


def _receipt(raw, expected):
    if expected is not None:
        require(
            isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected) and sha(raw) == expected,
            "external_receipt_mismatch",
            "Snapshot receipt differs",
        )
    return {
        "snapshot_sha256": sha(raw),
        "external_receipt_supplied_and_matched": expected is not None,
        "meaning": "Bytes match declared hashes and, if supplied, an external snapshot digest; "
        "the caller is responsible for trusting that receipt. No independent acquisition authentication.",
    }


def import_snapshot(source_dir, output_dir, *, expected_snapshot_sha256=None):
    """Import explicitly selected records. Network acquisition is deliberately absent."""
    snapshot_raw = read_local(source_dir, "snapshot.json")
    snapshot = json_read(snapshot_raw)
    _check_snapshot(snapshot)
    receipt = _receipt(snapshot_raw, expected_snapshot_sha256)
    dataset = read_local(source_dir, snapshot["dataset"]["path"], MAX_DATASET)
    descriptor = read_local(source_dir, snapshot["descriptor"]["path"], MAX_RECORD)
    records, raws, metadata, parser_observation = _replay(snapshot, dataset, descriptor)
    target = Path(output_dir).absolute()
    require(not target.exists() and not target.is_symlink(), "output_exists", "Use a new bundle directory")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".ord-bundle-", dir=target.parent))
    try:
        # Preserve original source receipt bytes and paths; reserved-name collisions reject.
        names = [snapshot["dataset"]["path"], snapshot["descriptor"]["path"]]
        reserved = {"snapshot.json", "bundle.json", "NOTICE.txt", "LICENSE-DATA", "LICENSE-SCHEMA"}
        require(
            len(set(names)) == 2 and not any(n in reserved or n.startswith("record-") for n in names),
            "asset_collision",
            "Source asset name collides with bundle files",
        )
        (temporary / "snapshot.json").write_bytes(snapshot_raw)
        for name, raw in zip(names, (dataset, descriptor)):
            (temporary / name).write_bytes(raw)
        for index, (record, raw) in enumerate(zip(records, raws)):
            record["raw_path"] = f"record-{index}.pb"
            (temporary / record["raw_path"]).write_bytes(raw)
        for name in ("NOTICE.txt", "LICENSE-DATA", "LICENSE-SCHEMA"):
            shutil.copyfile(RESOURCE / name, temporary / name)
        bundle = {
            "schema": VERSION,
            "receipt": receipt,
            "import_runtime": {**_runtime(), "parser_observation": parser_observation},
            "limits": LIMITS,
            "dataset_metadata": metadata,
            "records": records,
            "interpretation": INTERPRETATION,
        }
        (temporary / "bundle.json").write_bytes(encoded(bundle) + b"\n")
        temporary.rename(target)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return {
        "bundle": str(target),
        "selected_records": len(records),
        "rejected_records": sum(bool(r["derived"]["rejections"]) for r in records),
        "receipt": receipt,
    }


def load_bundle(path, *, expected_snapshot_sha256=None):
    """Always recheck bytes, coordinates and all derived values. No unsafe skip flag."""
    raw_bundle = read_local(path, "bundle.json")
    bundle = json_read(raw_bundle)
    require(
        isinstance(bundle, dict)
        and set(bundle)
        == {"schema", "receipt", "import_runtime", "limits", "dataset_metadata", "records", "interpretation"}
        and bundle["schema"] == VERSION,
        "bundle_schema",
        "Unsupported bundle",
    )
    snapshot_raw = read_local(path, "snapshot.json")
    snapshot = json_read(snapshot_raw)
    _check_snapshot(snapshot)
    receipt = _receipt(snapshot_raw, expected_snapshot_sha256)
    require(
        isinstance(bundle["receipt"], dict) and bundle["receipt"].get("snapshot_sha256") == sha(snapshot_raw),
        "snapshot_binding",
        "Snapshot receipt differs from imported bundle",
    )
    dataset = read_local(path, snapshot["dataset"]["path"], MAX_DATASET)
    descriptor = read_local(path, snapshot["descriptor"]["path"], MAX_RECORD)
    records, raws, metadata, parser_observation = _replay(snapshot, dataset, descriptor)
    require(
        isinstance(bundle["records"], list) and len(bundle["records"]) == len(records),
        "record_count",
        "Bundle record count differs",
    )
    for index, (record, raw, saved) in enumerate(zip(records, raws, bundle["records"])):
        record["raw_path"] = f"record-{index}.pb"
        require(
            read_local(path, record["raw_path"], MAX_RECORD) == raw,
            "raw_record_binding",
            "Retained protobuf differs from selected dataset row",
        )
        require(saved == record, "derived_record_mismatch", "Saved fields differ from raw rederivation")
    require(
        bundle["dataset_metadata"] == metadata
        and bundle["limits"] == LIMITS
        and bundle["interpretation"] == INTERPRETATION,
        "derived_bundle_mismatch",
        "Saved metadata differs from replay",
    )
    # Never inherit importer/caller assertions of trust. The current invocation supplies the receipt.
    return {
        "schema": VERSION,
        "bundle_sha256": sha(raw_bundle),
        "receipt": receipt,
        "source": snapshot,
        "metadata": metadata,
        "records": records,
        "runtime": {**_runtime(), "parser_observation": parser_observation},
        "import_runtime_as_reported": bundle["import_runtime"],
        "interpretation": INTERPRETATION,
    }


def record_context(record, outcome_index):
    """Lossless reported protocol view; equality is deliberately stricter than unit conversion."""
    projection = record["derived"]["projection"]
    outcome = projection["outcomes"][outcome_index]
    return {
        **{
            key: copy.deepcopy(projection[key])
            for key in ("inputs", "setup", "conditions", "notes", "workups")
            if key in projection
        },
        "outcome": {
            key: copy.deepcopy(outcome[key]) for key in ("reaction_time", "analyses") if key in outcome
        },
    }


def _context_reasons(record, context):
    reasons = []
    for component in record["derived"]["inputs"]:
        if component["role"] in {"UNSPECIFIED", "CUSTOM"}:
            reasons.append(component["path"] + ":unknown_role")
        if component["identity"]["status"] != "exact_absolute_state":
            reasons.append(component["path"] + ":unresolved_identity")
        if not _quantitative_amount(component["component"].get("amount", {})):
            reasons.append(component["path"] + ":amount_missing_invalid_or_unmeasured")
    for key in ("inputs", "setup", "conditions", "workups"):
        if not context.get(key):
            reasons.append(key + ":missing_or_empty_source_cannot_establish_complete_context")
    time = context["outcome"].get("reaction_time", {})
    if "value" not in time or not time.get("units"):
        reasons.append("outcome:time_missing_or_unspecified_units")
    temperature = context.get("conditions", {}).get("temperature", {}).get("setpoint", {})
    if "value" not in temperature or not temperature.get("units"):
        reasons.append("temperature:missing_or_unspecified_units")
    protocol_inputs = list(context.get("inputs", {}).values()) + [
        w["input"] for w in context.get("workups", []) if "input" in w
    ]
    for index, protocol_input in enumerate(protocol_inputs):
        if protocol_input.get("crude_components"):
            reasons.append(f"protocol_input/{index}:crude_composition_unresolved")
        for component in protocol_input.get("components", []):
            if not _quantitative_amount(component.get("amount", {})):
                reasons.append(f"protocol_input/{index}:amount_missing_invalid_or_unmeasured")
            if identity(component.get("identifiers", []))["key"] is None:
                reasons.append(f"protocol_input/{index}:identity_unresolved")
            if component.get("reaction_role", "UNSPECIFIED") in {"UNSPECIFIED", "CUSTOM"}:
                reasons.append(f"protocol_input/{index}:role_unresolved")
    if not any(x["role"] == "SOLVENT" for x in record["derived"]["inputs"]):
        reasons.append("solvent_inventory:not_explicitly_established")
    return reasons


def _match(loaded, product, reactants=None, context=None):
    target = smiles_identity(product)
    require(
        reactants is None or isinstance(reactants, list) and 0 < len(reactants) <= 256,
        "query_reactants",
        "Reactants must be a nonempty explicit instance list or absent",
    )
    require(
        context is None or isinstance(context, dict), "query_context", "Context must be an object or absent"
    )
    query_inputs = None if reactants is None else [smiles_identity(s) for s in reactants]
    matches = []
    if target["key"] is not None:
        for record in loaded["records"]:
            derived = record["derived"]
            if derived["rejections"]:
                continue
            inventory = derived["inputs"]
            reported = [x["identity"] for x in inventory if x["role"] == "REACTANT"]
            unknown = any(x["role"] in {"UNSPECIFIED", "CUSTOM"} for x in inventory)
            for item in derived["products"]:
                if item["role"] != "PRODUCT" or item["identity"]["key"] != target["key"]:
                    continue
                reasons = []
                if query_inputs is None:
                    reactant_status = "not_compared"
                elif not reported or any(x["key"] is None for x in reported + query_inputs):
                    reactant_status = "unassessable"
                else:
                    reactant_status = (
                        "reported_multiset_matched"
                        if sorted(x["key"] for x in reported) == sorted(x["key"] for x in query_inputs)
                        else "different"
                    )
                if unknown:
                    reasons.append("unknown_input_roles_prevent_complete_reactant_inventory_claim")
                recorded_context = record_context(record, item["outcome_index"])
                reasons.extend(_context_reasons(record, recorded_context))
                context_status = (
                    "not_compared"
                    if context is None
                    else ("same_recorded_context" if context == recorded_context else "different_or_partial")
                )
                if context is None:
                    reasons.append("query_context_absent")
                elif context_status != "same_recorded_context":
                    reasons.append("query_context_differs_or_is_partial")
                scoped = reactant_status == "reported_multiset_matched" and not reasons
                outcome = derived["projection"]["outcomes"][item["outcome_index"]]
                measurements = copy.deepcopy(item["product"].get("measurements", []))
                matches.append(
                    {
                        "reaction_id": derived["reaction_id"],
                        "record_sha256": record["raw_sha256"],
                        "selector": record["selector"],
                        "product_path": item["path"],
                        "product_identity": item["identity"],
                        "reactants": reactant_status,
                        "context": context_status,
                        "unresolved_reasons": reasons,
                        "claim_scope": "reported_reaction_with_matching_recorded_context"
                        if scoped
                        else "exact_product_recorded_precedent",
                        "scope_limit": "Reported compound entry; not isolated material, whole-output equality, "
                        "route validation or laboratory authentication",
                        "measurements": measurements,
                        "analyses": outcome.get("analyses", {}),
                        "wire_floats": [
                            x for x in derived["wire_floats"] if x["path"].startswith(item["path"] + "/")
                        ],
                        "presence": {
                            k: v for k, v in derived["presence"].items() if k.startswith(item["path"] + "/")
                        },
                        "recorded_context": recorded_context,
                        "provenance_as_reported": derived["projection"].get("provenance", {}),
                    }
                )
    unassessable = [
        {
            "reaction_id": row["derived"]["reaction_id"],
            "path": item["path"],
            "reasons": item["identity"]["reasons"],
        }
        for row in loaded["records"]
        for item in row["derived"]["products"]
        if item["role"] in {"PRODUCT", "UNSPECIFIED", "CUSTOM"}
        and (item["identity"]["key"] is None or item["role"] != "PRODUCT")
    ]
    rejected = [
        {"reaction_id": row["derived"]["reaction_id"], "reasons": row["derived"]["rejections"]}
        for row in loaded["records"]
        if row["derived"]["rejections"]
    ]
    return {
        "query": {"product": product, "product_identity": target, "reactants": reactants, "context": context},
        "status": "query_identity_unassessable"
        if target["key"] is None
        else "matched_recorded_product"
        if matches
        else "search_incomplete_no_assessable_match"
        if rejected or unassessable
        else "no_match_in_retained_snapshot",
        "matches": matches,
        "rejected_records": rejected,
        "unassessable_products": unassessable,
        "search_completeness": "incomplete"
        if rejected or unassessable
        else "assessable_selected_entries_only",
    }


def match_product(bundle_path, product, *, reactants=None, context=None, expected_snapshot_sha256=None):
    """Public exact-product lookup; cannot consume a caller-edited normalized bundle object."""
    loaded = load_bundle(bundle_path, expected_snapshot_sha256=expected_snapshot_sha256)
    return {
        "schema": VERSION,
        "interpretation": INTERPRETATION,
        "integrity": {k: loaded[k] for k in ("bundle_sha256", "receipt", "runtime")},
        "source_as_declared": loaded["source"],
        "result": _match(loaded, product, reactants, context),
    }
