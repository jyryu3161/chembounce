"""Source-record traceability fixtures, not supplier or synthesis validation."""

import copy
import csv
import hashlib
import json
from pathlib import Path

import pytest

from chembounce_synthesis import assets as assets_module
from chembounce_synthesis.assets import convert_catalog
from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.common import digest
from chembounce_synthesis.route_audit import audit_route
from chembounce_synthesis.routes import verify_route


def converted(tmp_path, rows=(("a", "CCO"), ("b", "CCN"), ("bad", "invalid"))):
    source = tmp_path / "source.csv"
    with source.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("sku", "structure"))
        writer.writerows(rows)
    output = tmp_path / "converted.tsv"
    manifest = convert_catalog([source], output, "fixture", "snapshot", "structure", "sku")
    return source, output, output.with_suffix(".tsv.manifest.json"), manifest


def test_mixed_records_and_multiline_id_preserve_exact_source_coordinates(tmp_path):
    source, output, sidecar, manifest = converted(
        tmp_path, ((" supplied ", "CCO"), ("bad", "invalid"), ("multi\nline", "CCN"), ("  ", "CCC"))
    )
    catalog = Catalog(output)
    generated = "structure-sha256:" + hashlib.sha256(b"CCC").hexdigest()
    assert manifest["schema_version"] == 2
    assert manifest["records"] == 3
    assert manifest["smiles_column"] == "structure" and manifest["id_column"] == "sku"
    assert manifest["rejected_records"][0]["input_index"] == 0
    assert manifest["rejected_records"][0]["record"] == 1
    for ident, index, origin in (
        ("supplied", 0, "input_provided"),
        ("multi\nline", 2, "input_provided"),
        (generated, 3, "generated_structure_sha256"),
    ):
        evidence = catalog.evidence(ident)
        assert evidence["identifier_provenance"] == {
            "status": "recorded",
            "manifest_sha256": digest(sidecar),
            "source_records": [
                {
                    "input_index": 0,
                    "input_path": str(source.resolve()),
                    "input_sha256": digest(source),
                    "record_index": index,
                    "identifier_origin": origin,
                }
            ],
        }
        assert evidence["supplier_sku_verified"] is False
        assert evidence["current_stock_verified"] is False
    assert "record_origins" not in catalog.import_provenance


def test_all_duplicate_origins_including_generated_supplied_collision(tmp_path):
    generated = "structure-sha256:" + hashlib.sha256(b"CCO").hexdigest()
    first = tmp_path / "first.csv"
    first.write_text(f"id,smiles\n,CCO\n{generated},OCC\n")
    second = tmp_path / "second.csv"
    second.write_text(f"id,smiles\n {generated} ,CCO\nstructure-sha256:supplied,CCN\n")
    output = tmp_path / "catalog.tsv"
    manifest = convert_catalog([first, second, first], output, "fixture", "snapshot")
    assert manifest["records"] == 2
    records = Catalog(output).evidence(generated)["identifier_provenance"]["source_records"]
    assert [(r["input_index"], r["record_index"], r["identifier_origin"]) for r in records] == [
        (0, 0, "generated_structure_sha256"),
        (0, 1, "input_provided"),
        (1, 0, "input_provided"),
        (2, 0, "generated_structure_sha256"),
        (2, 1, "input_provided"),
    ]
    prefix = Catalog(output).evidence("structure-sha256:supplied")
    assert prefix["identifier_provenance"]["source_records"][0]["identifier_origin"] == "input_provided"


@pytest.mark.parametrize("legacy", ["absent", "v1", "unversioned"])
def test_legacy_never_infers_identifier_origin_and_stock_is_unverified(tmp_path, legacy):
    output = tmp_path / "catalog.tsv"
    output.write_text(
        "id\tsmiles\tsource\tsnapshot\tavailability\n"
        "structure-sha256:looks-generated\tCCO\tfixture\tsnapshot\tin_stock\n"
    )
    sidecar = output.with_suffix(".tsv.manifest.json")
    if legacy != "absent":
        manifest = {"output_sha256": digest(output)}
        if legacy == "v1":
            manifest["schema_version"] = 1
        sidecar.write_text(json.dumps(manifest))
    catalog = Catalog(output)
    evidence = catalog.evidence("structure-sha256:looks-generated")
    assert evidence["identifier_provenance"] == {
        "status": "unknown_legacy",
        "source_records": [],
        "manifest_sha256": digest(sidecar) if legacy != "absent" else None,
    }
    assert evidence["availability"] == "in_stock"
    assert evidence["current_stock_verified"] is False
    assert evidence["supplier_sku_verified"] is False


@pytest.mark.parametrize(
    "change",
    [
        "bool_count",
        "wrong_count",
        "bool_input",
        "bool_record",
        "negative_record",
        "input_out_of_range",
        "duplicate_coordinate",
        "duplicate_across_ids",
        "accepted_rejected_overlap",
        "rejected_duplicate",
        "rejected_bool_index",
        "missing_id",
        "extra_id",
        "empty_origins",
        "unknown_origin",
        "bad_hash",
        "bad_path",
        "bad_generated_identity",
        "wrong_source",
        "wrong_snapshot",
        "wrong_index_policy",
        "missing_column_setting",
        "missing_origins",
        "missing_rejections",
    ],
)
def test_malformed_v2_fails_closed(tmp_path, change):
    _, output, sidecar, manifest = converted(tmp_path)
    origins = manifest["record_origins"]
    first = origins["a"][0]
    if change == "bool_count":
        manifest["records"] = True
    elif change == "wrong_count":
        manifest["records"] = 3
    elif change == "bool_input":
        first["input_index"] = False
    elif change == "bool_record":
        first["record_index"] = False
    elif change == "negative_record":
        first["record_index"] = -1
    elif change == "input_out_of_range":
        first["input_index"] = 1
    elif change == "duplicate_coordinate":
        origins["a"].append(copy.deepcopy(first))
    elif change == "duplicate_across_ids":
        origins["b"] = [copy.deepcopy(first)]
    elif change == "accepted_rejected_overlap":
        manifest["rejected_records"][0]["record"] = 0
    elif change == "rejected_duplicate":
        manifest["rejected_records"] *= 2
    elif change == "rejected_bool_index":
        manifest["rejected_records"][0]["input_index"] = False
    elif change == "missing_id":
        del origins["b"]
    elif change == "extra_id":
        origins["unknown"] = [copy.deepcopy(first)]
    elif change == "empty_origins":
        origins["a"] = []
    elif change == "unknown_origin":
        first["identifier_origin"] = "verified_supplier"
    elif change == "bad_hash":
        manifest["inputs"][0]["sha256"] = "x" * 64
    elif change == "bad_path":
        manifest["inputs"][0]["path"] = ""
    elif change == "bad_generated_identity":
        first["identifier_origin"] = "generated_structure_sha256"
    elif change == "wrong_source":
        manifest["source"] = "other"
    elif change == "wrong_snapshot":
        manifest["snapshot"] = "other"
    elif change == "wrong_index_policy":
        manifest["record_index_policy"] = "physical_lines"
    elif change == "missing_column_setting":
        del manifest["id_column"]
    elif change == "missing_origins":
        del manifest["record_origins"]
    elif change == "missing_rejections":
        del manifest["rejected_records"]
    sidecar.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Invalid catalog v2 record provenance"):
        Catalog(output)


@pytest.mark.parametrize("version", [3, "2", True, None])
def test_unknown_explicit_schema_rejected(tmp_path, version):
    _, output, sidecar, manifest = converted(tmp_path)
    manifest["schema_version"] = version
    sidecar.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Unsupported"):
        Catalog(output)


def test_unavailable_records_still_require_origins(tmp_path):
    _, output, sidecar, manifest = converted(tmp_path)
    output.write_text(
        output.read_text().replace(
            "b\tCCN\tfixture\tsnapshot\tcatalogued", "b\tCCN\tfixture\tsnapshot\tunavailable"
        )
    )
    manifest["output_sha256"] = digest(output)
    sidecar.write_text(json.dumps(manifest))
    assert Catalog(output).blocks["b"].availability == "unavailable"
    del manifest["record_origins"]["b"]
    sidecar.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="provenance"):
        Catalog(output)


def test_input_drift_during_parsing_fails_before_publication(tmp_path, monkeypatch):
    source = tmp_path / "source.csv"
    source.write_text("id,smiles\na,CCO\n")
    output = tmp_path / "catalog.tsv"
    canonical = assets_module.canonical

    def changing_input(smiles):
        source.write_text("id,smiles\na,CCN\n")
        return canonical(smiles)

    monkeypatch.setattr(assets_module, "canonical", changing_input)
    with pytest.raises(ValueError, match="input changed"):
        convert_catalog([source], output, "fixture", "snapshot")
    assert not output.exists() and not output.with_suffix(".tsv.manifest.json").exists()


def test_normalized_conflicting_ids_fail_before_publication(tmp_path):
    source = tmp_path / "source.csv"
    source.write_text("id,smiles\n a ,CCO\na,CCN\n")
    output = tmp_path / "catalog.tsv"
    with pytest.raises(ValueError, match="Conflicting structures"):
        convert_catalog([source], output, "fixture", "snapshot")
    assert not output.exists() and not output.with_suffix(".tsv.manifest.json").exists()


def test_smiles_logical_indices_exclude_comments_and_blanks(tmp_path):
    source = tmp_path / "source.smi"
    source.write_text("# comment\n\nCCO first\ninvalid bad\nCCN\n")
    output = tmp_path / "catalog.tsv"
    manifest = convert_catalog([source], output, "fixture", "snapshot")
    assert manifest["record_origins"]["first"][0]["record_index"] == 0
    assert manifest["rejected_records"][0]["record"] == 1
    generated = "structure-sha256:" + hashlib.sha256(b"CCN").hexdigest()
    assert manifest["record_origins"][generated][0]["record_index"] == 2


def test_evidence_copies_and_route_audit_export_parity(tmp_path, assets, route):
    from benchmarks.research.qualified_collector import BoundCatalog
    from chembounce_synthesis import catalog as catalog_module, common

    _, output, sidecar, _ = converted(tmp_path, (("acetic", "CC(=O)O"), ("pyridyl", "Nc1ccncc1")))
    catalog = Catalog(output)
    evidence = catalog.evidence("acetic")
    evidence["identifier_provenance"]["source_records"][0]["record_index"] = 999
    assert catalog.evidence("acetic")["identifier_provenance"]["source_records"][0]["record_index"] == 0
    evidence["identifier_provenance"]["source_records"].clear()
    catalog.import_provenance["inputs"][0]["sha256"] = "0" * 64
    verification = verify_route(route, catalog, assets[1])
    audit = audit_route(route, catalog, assets[1])
    bound = BoundCatalog(
        output,
        catalog.sha256,
        {"catalog.py": digest(catalog_module.__file__), "common.py": digest(common.__file__)},
        sidecar_sha256=digest(sidecar),
    )
    assert audit_route(route, bound, assets[1]) == audit
    for verified, audited in zip(verification["building_blocks"], audit["building_blocks"]):
        assert {k: v for k, v in verified.items() if k != "node_id"} == audited["catalog_record"]
        assert audited["catalog_record"]["identifier_provenance"]["status"] == "recorded"
    assert verification["experimental_validation"] is False
    assert catalog.evidence("acetic")["identifier_provenance"]["source_records"][0]["record_index"] == 0
    assert (
        catalog.evidence("acetic")["identifier_provenance"]["source_records"][0]["input_sha256"] != "0" * 64
    )


@pytest.mark.parametrize("change", ["replace", "remove", "add"])
def test_audit_detects_sidecar_freshness_changes(tmp_path, assets, route, change):
    _, output, sidecar, manifest = converted(tmp_path, (("acetic", "CC(=O)O"), ("pyridyl", "Nc1ccncc1")))
    if change == "add":
        sidecar.unlink()
    catalog = Catalog(output)
    audit_route(route, catalog, assets[1])
    if change == "remove":
        sidecar.unlink()
    else:
        sidecar.write_text(json.dumps(manifest) + "\n")
    with pytest.raises(ValueError, match="import manifest changed"):
        audit_route(route, catalog, assets[1])


def test_manifest_hash_binds_the_exact_bytes_parsed(tmp_path, monkeypatch):
    _, output, sidecar, _ = converted(tmp_path)
    original_bytes = sidecar.read_bytes()
    read_bytes = Path.read_bytes

    def replacing_sidecar(path):
        data = read_bytes(path)
        if path == sidecar:
            path.write_bytes(data + b" ")
        return data

    monkeypatch.setattr(Path, "read_bytes", replacing_sidecar)
    catalog = Catalog(output)
    assert (
        catalog.evidence("a")["identifier_provenance"]["manifest_sha256"]
        == hashlib.sha256(original_bytes).hexdigest()
    )
    with pytest.raises(ValueError, match="import manifest changed"):
        catalog._check_import_unchanged()
