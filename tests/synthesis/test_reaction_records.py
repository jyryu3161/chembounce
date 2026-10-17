"""Real adverse ORD source fixture plus explicitly synthetic adversarial mutations.

No model, rule replay, scoring, catalog or campaign data. Mutated records below
are test-only data and never claimed to be public experiments.
"""

import json
import os
import shutil
import subprocess
import sys

import pytest

from chembounce_synthesis import reaction_records as r
from chembounce_synthesis import reaction_record_report as report

PRODUCT = "C(=C/c1ccccc1)\\c1ccccc1"
REACTANTS = ["O=S(=O)(Oc1ccccc1)C(F)(F)F", "C=Cc1ccccc1"]
ID = "ord-268b5629ab3849f890d911d8df2bb683"


@pytest.fixture(scope="module")
def source():
    pytest.importorskip("pyarrow")
    pytest.importorskip("google.protobuf")
    return r.RESOURCE


@pytest.fixture(scope="module")
def bundle(tmp_path_factory, source):
    target = tmp_path_factory.mktemp("ord-real") / "bundle"
    r.import_snapshot(source, target)
    return target


def save_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True) + "\n")


def edit_bundle(tmp_path, bundle, mutation):
    target = tmp_path / "edited"
    shutil.copytree(bundle, target)
    path = target / "bundle.json"
    obj = json.loads(path.read_text())
    mutation(obj)
    save_json(path, obj)
    return target


def synthetic_source(tmp_path, source, mutate, *, duplicate=False, mutate_snapshot=None):
    """Rewrite selected row from public fixture into explicitly synthetic test data."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    from google.protobuf import descriptor_pool, message_factory

    pool = descriptor_pool.DescriptorPool()
    pool.AddSerializedFile((source / "reaction.descriptor.pb").read_bytes())
    cls = message_factory.GetMessageClass(pool.FindMessageTypeByName("ord.Reaction"))
    table = pq.read_table(source / "sample.parquet")
    first = table.slice(0, 1).to_pylist()[0]
    msg = cls.FromString(first["reaction"])
    mutate(msg)
    raw = msg.SerializeToString()
    rows = [{"reaction_id": ID, "reaction": raw}]
    if duplicate:
        rows.append(first)
    target = tmp_path / "synthetic-source"
    shutil.copytree(source, target)
    pq.write_table(
        pa.Table.from_pylist(rows, schema=table.schema), target / "sample.parquet", compression="gzip"
    )
    snapshot = json.loads((target / "snapshot.json").read_text())
    data = (target / "sample.parquet").read_bytes()
    snapshot["dataset"].update(sha256=r.sha(data), bytes=len(data))
    snapshot["source_metadata"]["synthetic_test_only"] = True
    if duplicate:
        snapshot["selectors"].append({"row_group": 0, "row": 1, "reaction_id": ID})
    if mutate_snapshot:
        mutate_snapshot(snapshot)
    save_json(target / "snapshot.json", snapshot)
    return target


def test_real_import_match_render_and_wire_precision(bundle):
    loaded = r.load_bundle(bundle)
    rec = loaded["records"][0]
    assert rec["derived"]["reaction_id"] == ID
    assert rec["raw_sha256"] == "5bdd1c9535bca04c5b649d8350d50befb77ff11ad156a46cb543407841f325ff"
    assert (
        rec["derived"]["projection"]["outcomes"][0]["products"][1]["measurements"][0]["float_value"]["value"]
        == 163541580.0
    )
    assert (
        next(
            x["decoded_value"]
            for x in rec["derived"]["wire_floats"]
            if x["path"] == "/outcomes/0/products/1/measurements/0/float_value/value"
        )
        == 163541584.0
    )
    products = rec["derived"]["projection"]["outcomes"][0]["products"]
    assert sum(p["measurements"][0]["float_value"]["value"] == 0 for p in products) == 4
    presence = rec["derived"]["presence"]
    assert presence["/outcomes/0/products/0/measurements/0/uses_authentic_standard"]["state"] == "present"
    assert presence["/outcomes/0/products/1/measurements/0/uses_authentic_standard"]["state"] == "absent"
    result = r.match_product(bundle, PRODUCT, reactants=REACTANTS, context=r.record_context(rec, 0))
    match = result["result"]["matches"][0]
    assert match["reactants"] == "reported_multiset_matched"
    assert match["context"] == "same_recorded_context"
    assert match["claim_scope"] == "exact_product_recorded_precedent"
    assert any("unknown_input_roles" in x for x in match["unresolved_reasons"])
    assert any("Ligand" in x and "unresolved_identity" in x for x in match["unresolved_reasons"])
    assert match["measurements"][0]["type"] == "AREA"
    assert match["analyses"]["GC-FID"]["is_of_isolated_species"] is False
    text = report.render_captured_report(result)
    assert "AREA" in text and "163541584" in text and "163541580" in text
    assert "does not freshly verify" in text
    assert "source_audited" not in json.dumps(result)
    observation = loaded["runtime"]["parser_observation"]
    assert observation["address_space_cap_bytes"] - observation["baseline_virtual_bytes"] <= 256 * 1024**2


def test_receipt_is_external_and_never_inherited(bundle, tmp_path):
    digest = r.sha((bundle / "snapshot.json").read_bytes())
    assert r.load_bundle(bundle, expected_snapshot_sha256=digest)["receipt"][
        "external_receipt_supplied_and_matched"
    ]
    forged = edit_bundle(
        tmp_path, bundle, lambda b: b["receipt"].update(external_receipt_supplied_and_matched=True)
    )
    assert not r.load_bundle(forged)["receipt"]["external_receipt_supplied_and_matched"]
    with pytest.raises(r.RecordError, match="external_receipt_mismatch"):
        r.load_bundle(bundle, expected_snapshot_sha256="0" * 64)
    with pytest.raises(TypeError):
        r.load_bundle(bundle, verify_hashes=False)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda b: b["records"][0]["derived"]["products"][1]["identity"].update(key="CCC"),
        lambda b: b["records"][0]["derived"]["inputs"][0].update(role="REACTANT"),
        lambda b: b["records"][0]["derived"]["products"][1]["product"]["measurements"][0].update(
            type="YIELD"
        ),
        lambda b: b["records"][0]["derived"]["projection"]["provenance"].update(doi="invented"),
        lambda b: b["records"][0]["selector"].update(row=1),
        lambda b: b["records"][0].update(raw_sha256="0" * 64),
        lambda b: b["records"][0]["derived"]["wire_floats"][0].update(decoded_value=123),
    ],
)
def test_normalized_tampering_rejected(bundle, tmp_path, mutation):
    path = edit_bundle(tmp_path, bundle, mutation)
    with pytest.raises(r.RecordError, match="derived_record_mismatch"):
        r.match_product(path, PRODUCT)


@pytest.mark.parametrize("filename", ["sample.parquet", "reaction.descriptor.pb", "record-0.pb"])
def test_raw_tamper(bundle, tmp_path, filename):
    target = tmp_path / "raw-edit"
    shutil.copytree(bundle, target)
    path = target / filename
    path.write_bytes(path.read_bytes() + b"tampered")
    with pytest.raises(r.RecordError):
        r.load_bundle(target)


@pytest.mark.parametrize(
    "field,value", [("row", 96), ("row_group", 1), ("row", True), ("reaction_id", "ord-" + "0" * 32)]
)
def test_selector_tamper(source, tmp_path, field, value):
    target = tmp_path / "source"
    shutil.copytree(source, target)
    path = target / "snapshot.json"
    snap = json.loads(path.read_text())
    snap["selectors"][0][field] = value
    save_json(path, snap)
    with pytest.raises(r.RecordError):
        r.import_snapshot(target, tmp_path / "bundle")


def test_schema_binding_cannot_be_rewritten(source, tmp_path):
    target = tmp_path / "source"
    shutil.copytree(source, target)
    (target / "reaction.descriptor.pb").write_bytes(b"fake")
    snap = json.loads((target / "snapshot.json").read_text())
    snap["descriptor"]["sha256"] = r.sha(b"fake")
    save_json(target / "snapshot.json", snap)
    with pytest.raises(r.RecordError, match="unsupported_schema"):
        r.import_snapshot(target, tmp_path / "bundle")


def test_duplicate_ids_conflicting_payloads(source, tmp_path):
    target = synthetic_source(
        tmp_path, source, lambda m: setattr(m.provenance, "city", "SYNTHETIC"), duplicate=True
    )
    with pytest.raises(r.RecordError, match="conflicting_record_id"):
        r.import_snapshot(target, tmp_path / "bundle")


@pytest.mark.parametrize(
    "left,right",
    [
        ("F[C@H](Cl)Br", "F[C@@H](Cl)Br"),
        ("F/C=C/F", "F/C=C\\F"),
        ("[13CH4]", "C"),
        ("[NH4+]", "N"),
        ("CC.[Na+]", "CC"),
        ("CC.CC", "CC"),
    ],
)
def test_exact_identity_negative(left, right):
    assert r.smiles_identity(left)["key"] != r.smiles_identity(right)["key"]


@pytest.mark.parametrize("text", ["FC(Cl)Br", "FC=CF", "F[C@H](Cl)Br |o1:1|", "CC |wU:0.0|", "*C"])
def test_unsupported_or_unspecified_identity_abstains(text):
    assert r.smiles_identity(text)["key"] is None


def test_conflicting_smiles_molblock_and_names_only():
    from rdkit import Chem

    identifiers = [
        {"type": "SMILES", "value": "CC"},
        {"type": "MOLBLOCK", "value": Chem.MolToMolBlock(Chem.MolFromSmiles("CO"))},
    ]
    assert "conflicting_structural_identifiers" in r.identity(identifiers)["reasons"]
    assert r.identity([{"type": "NAME", "value": "ethane"}])["key"] is None
    assert r.smiles_identity("[CH3:1][OH:2]")["key"] == r.smiles_identity("CO")["key"]


def test_no_internal_standard_match_and_unknown_query(bundle):
    assert (
        r.match_product(bundle, "COc1cc(OC)cc(OC)c1")["result"]["status"]
        == "search_incomplete_no_assessable_match"
    )
    assert r.match_product(bundle, "FC(Cl)Br")["result"]["status"] == "query_identity_unassessable"
    assert r.match_product(bundle, "CCC")["result"]["status"] == "search_incomplete_no_assessable_match"


def test_multiplicity_and_missing_context(bundle):
    result = r.match_product(bundle, PRODUCT, reactants=REACTANTS + [REACTANTS[0]])["result"]["matches"][0]
    assert result["reactants"] == "different" and result["context"] == "not_compared"
    rec = r.load_bundle(bundle)["records"][0]
    context = r.record_context(rec, 0)
    context["conditions"]["temperature"]["setpoint"]["value"] = 90
    result = r.match_product(bundle, PRODUCT, reactants=REACTANTS, context=context)["result"]["matches"][0]
    assert result["context"] == "different_or_partial"
    assert (
        r.match_product(bundle, PRODUCT, context={})["result"]["matches"][0]["context"]
        == "different_or_partial"
    )


@pytest.mark.parametrize(
    "mutate,code",
    [
        (
            lambda m: setattr(m.outcomes[0].products[1].measurements[0], "analysis_key", "missing"),
            "unresolved_analysis_reference",
        ),
        (
            lambda m: setattr(m.outcomes[0].products[1].measurements[0].float_value, "value", -1),
            "measurement_value_missing_or_negative",
        ),
        (
            lambda m: setattr(m.outcomes[0].products[1].measurements[0], "type", 3),
            "measurement_value_type_mismatch_or_missing",
        ),
    ],
)
def test_typed_measurement_rejections_retained(source, tmp_path, mutate, code):
    target = synthetic_source(tmp_path, source, mutate)
    bundle = tmp_path / "bundle"
    imported = r.import_snapshot(target, bundle)
    assert imported["rejected_records"] == 1 and (bundle / "record-0.pb").exists()
    result = r.match_product(bundle, PRODUCT)["result"]
    assert any(x["code"] == code for x in result["rejected_records"][0]["reasons"])
    assert result["matches"] == []


@pytest.mark.parametrize(
    "mutate,code",
    [
        (
            lambda m: setattr(m.outcomes[0].products[1].measurements[0].float_value, "value", float("nan")),
            "nonfinite_record",
        ),
        (lambda m: setattr(m.inputs["Metal"].components[0], "reaction_role", 999), "unknown_enum"),
    ],
)
def test_nonfinite_unknown_enum(source, tmp_path, mutate, code):
    target = synthetic_source(tmp_path, source, mutate)
    with pytest.raises(r.RecordError, match=code):
        r.import_snapshot(target, tmp_path / "bundle")


def test_unknown_wire_fields(source, tmp_path):
    target = synthetic_source(tmp_path, source, lambda m: m.MergeFromString(b"\xf8\x7f\x01"))
    with pytest.raises(r.RecordError, match="unknown_protobuf_fields"):
        r.import_snapshot(target, tmp_path / "bundle")


def test_unsafe_paths_symlinks_and_file_bound(source, tmp_path):
    for name in ("../snapshot.json", "/etc/passwd", "a/b", "a\\b"):
        with pytest.raises(r.RecordError, match="unsafe_path"):
            r.read_local(source, name)
    (tmp_path / "link").symlink_to(source / "snapshot.json")
    with pytest.raises(r.RecordError):
        r.read_local(tmp_path, "link")
    with pytest.raises(r.RecordError, match="file_bound"):
        r.read_local(source, "sample.parquet", 10)
    with pytest.raises(r.RecordError, match="duplicate_json_key"):
        r.json_read(b'{"a":1,"a":2}')


def test_decompressed_cap_on_small_compressed_source(source, tmp_path):
    target = synthetic_source(tmp_path, source, lambda m: setattr(m.provenance, "city", "X" * (17 * 1024**2)))
    assert (target / "sample.parquet").stat().st_size < r.MAX_DATASET
    with pytest.raises(r.RecordError, match="decompressed_bound"):
        r.import_snapshot(target, tmp_path / "bundle")


def make_saved_run(path):
    path.mkdir()
    nodes = [
        {"id": "a", "building_block_id": "a"},
        {"id": "b", "building_block_id": "b"},
        {"id": "final", "rule_id": "synthetic-test-only", "reactants": ["a", "b"], "product": PRODUCT},
    ]
    route = {
        "route_id": "synthetic-test-route",
        "route": {"root": "final", "nodes": nodes},
        "verification": {
            "product": PRODUCT,
            "building_blocks": [{"id": "a", "smiles": REACTANTS[0]}, {"id": "b", "smiles": REACTANTS[1]}],
            "steps": [
                {
                    "node_id": "final",
                    "rule_id": "synthetic-test-only",
                    "product": PRODUCT,
                    "reactants": REACTANTS,
                }
            ],
        },
    }
    save_json(path / "routes.json", {"routes": [route]})
    save_json(path / "results.json", {"results": [{"route_id": "synthetic-test-route", "smiles": PRODUCT}]})
    return route


def test_saved_route_report_is_bound_and_cannot_change_acceptance(bundle, tmp_path):
    run = tmp_path / "synthetic-run"
    make_saved_run(run)
    before = {p.name: p.read_bytes() for p in run.iterdir()}
    document = report.match_saved_route(bundle, run, "synthetic-test-route")
    assert document["steps"][0]["matches"][0]["claim_scope"] == "exact_product_recorded_precedent"
    assert document["saved_run_binding"]["routes_json_sha256"] == r.sha(before["routes.json"])
    assert document["route_scope"]["reported_product_matches"] == 1
    report.write_report(document, tmp_path / "report")
    assert {p.name: p.read_bytes() for p in run.iterdir()} == before
    assert "synthetic-test-route" in (tmp_path / "report" / "index.html").read_text()
    results = json.loads((run / "results.json").read_text())
    results["results"][0]["smiles"] = "CCC"
    save_json(run / "results.json", results)
    with pytest.raises(r.RecordError, match="saved_final_identity"):
        report.match_saved_route(bundle, run, "synthetic-test-route")


def test_saved_route_multiplicity_and_graph_tampering(bundle, tmp_path):
    run = tmp_path / "run"
    route = make_saved_run(run)
    route["route"]["nodes"][-1]["reactants"] = ["a", "a", "b"]
    save_json(run / "routes.json", {"routes": [route]})
    with pytest.raises(r.RecordError, match="saved_step_binding"):
        report.match_saved_route(bundle, run, "synthetic-test-route")


def test_escape_captured_source_text(bundle):
    document = r.match_product(bundle, PRODUCT)
    document["result"]["matches"][0]["reaction_id"] = '<script>alert("bad")</script>'
    text = report.render_captured_report(document)
    assert "<script>" not in text and "&lt;script&gt;" in text
    assert "Rendering displays captured assertions" in text


def test_cli_actual_flow(source, tmp_path):
    assert report.main(["example", str(tmp_path / "snapshot")]) == 0
    assert report.main(["import", str(tmp_path / "snapshot"), str(tmp_path / "bundle")]) == 0
    assert (
        report.main(["match", str(tmp_path / "bundle"), str(tmp_path / "report"), "--product", PRODUCT]) == 0
    )
    assert (tmp_path / "report" / "index.html").exists()


def test_default_imports_and_capture_render_without_optional_modules(tmp_path):
    script = """
import sys
class Block:
    def find_spec(self, fullname, *args):
        if fullname == "pyarrow" or fullname.startswith("google.protobuf"):
            raise ImportError("test blocks optional dependencies")
sys.meta_path.insert(0, Block())
from chembounce_synthesis import report
from chembounce_synthesis import reaction_record_report as sidecar
from chembounce_synthesis import reaction_records as r
assert "pyarrow" not in sys.modules
assert "Captured public reaction-record report" in sidecar.render_captured_report({"schema":r.VERSION})
try:
    r._dependencies()
except r.RecordError as exc:
    assert exc.code == "optional_dependency_missing"
else:
    raise AssertionError("optional dependency refusal missing")
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
    "kind,value",
    [
        ("CXSMILES", "F[C@@H](Cl)Br"),
        ("CXSMILES", "F[C@H](Cl)Br |o1:1|"),
        ("INCHI", "InChI=1S/CH4/h1H4"),
        ("INCHI_KEY", "VNWKTOKETHGBQD-UHFFFAOYSA-N"),
        ("XYZ", "1\nX\nC 0 0 0"),
        ("AMINO_ACID_SEQUENCE", "G"),
        ("HELM", "PEPTIDE1{G}$$$$"),
    ],
)
def test_additional_structural_identifier_cannot_hide_conflict(kind, value):
    result = r.identity([{"type": "SMILES", "value": "F[C@H](Cl)Br"}, {"type": kind, "value": value}])
    assert result["key"] is None


def test_multiple_outcomes_and_optional_false_remain_distinct(source, tmp_path):
    def mutate(msg):
        added = msg.outcomes.add()
        added.CopyFrom(msg.outcomes[0])
        added.reaction_time.value = 24
        added.products[1].measurements[0].uses_authentic_standard = False

    target = synthetic_source(tmp_path, source, mutate)
    bundle = tmp_path / "bundle"
    r.import_snapshot(target, bundle)
    result = r.match_product(bundle, PRODUCT)["result"]
    assert len(result["matches"]) == 2
    assert "uses_authentic_standard" not in result["matches"][0]["measurements"][0]
    assert result["matches"][1]["measurements"][0]["uses_authentic_standard"] is False
    assert [m["recorded_context"]["outcome"]["reaction_time"]["value"] for m in result["matches"]] == [16, 24]


def test_record_byte_bound(source, tmp_path):
    target = synthetic_source(
        tmp_path, source, lambda m: setattr(m.provenance, "city", "X" * (r.MAX_RECORD + 1))
    )
    with pytest.raises(r.RecordError, match="record_bound"):
        r.import_snapshot(target, tmp_path / "bundle")


def test_snapshot_row_and_selector_bounds(source, tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq

    target = synthetic_source(tmp_path, source, lambda m: None)
    table = pq.read_table(target / "sample.parquet")
    pq.write_table(
        pa.concat_tables([table] * (r.MAX_ROWS + 1)), target / "sample.parquet", compression="gzip"
    )
    snap = json.loads((target / "snapshot.json").read_text())
    raw = (target / "sample.parquet").read_bytes()
    snap["dataset"].update(bytes=len(raw), sha256=r.sha(raw))
    save_json(target / "snapshot.json", snap)
    with pytest.raises(r.RecordError, match="row_bound"):
        r.import_snapshot(target, tmp_path / "bundle")
    snap["selectors"] *= 65
    save_json(target / "snapshot.json", snap)
    with pytest.raises(r.RecordError, match="selector_bound"):
        r.import_snapshot(target, tmp_path / "bundle")


def test_source_unsupported_schema_revision(source, tmp_path):
    target = synthetic_source(
        tmp_path, source, lambda m: None, mutate_snapshot=lambda s: s.update(schema_revision="0" * 40)
    )
    with pytest.raises(r.RecordError, match="unsupported_source_schema"):
        r.import_snapshot(target, tmp_path / "bundle")


def test_rejected_record_is_not_a_no_match_claim(source, tmp_path):
    target = synthetic_source(
        tmp_path, source, lambda m: setattr(m.outcomes[0].products[1].measurements[0], "type", 3)
    )
    r.import_snapshot(target, tmp_path / "bundle")
    result = r.match_product(tmp_path / "bundle", PRODUCT)["result"]
    assert result["status"] == "search_incomplete_no_assessable_match"
    assert result["rejected_records"] and result["search_completeness"] == "incomplete"


def test_multi_step_saved_binding_and_unknown_context(bundle, tmp_path):
    run = tmp_path / "run"
    route = make_saved_run(run)
    route["route"]["nodes"].insert(
        2, {"id": "mid", "rule_id": "synthetic-intermediate", "reactants": ["a"], "product": REACTANTS[0]}
    )
    route["route"]["nodes"][-1]["reactants"][0] = "mid"
    route["verification"]["steps"].insert(
        0,
        {
            "node_id": "mid",
            "rule_id": "synthetic-intermediate",
            "reactants": [REACTANTS[0]],
            "product": REACTANTS[0],
        },
    )
    save_json(run / "routes.json", {"routes": [route]})
    value = report.match_saved_route(bundle, run, "synthetic-test-route")
    assert [s["node_id"] for s in value["steps"]] == ["mid", "final"]
    assert value["route_scope"]["step_count"] == 2
    assert value["route_scope"]["reported_product_matches"] == 1
    assert value["steps"][1]["matches"][0]["context"] == "not_compared"


def test_synthetic_complete_reported_context_stays_report_only(source, tmp_path):
    def mutate(message):
        metal = message.inputs["Metal"].components[0]
        metal.reaction_role = (
            metal.DESCRIPTOR.fields_by_name["reaction_role"].enum_type.values_by_name["CATALYST"].number
        )
        message.inputs["Ligand"].components[0].identifiers[0].value = "P"

    target = synthetic_source(tmp_path, source, mutate)
    path = tmp_path / "bundle"
    r.import_snapshot(target, path)
    context = r.record_context(r.load_bundle(path)["records"][0], 0)
    result = r.match_product(path, PRODUCT, reactants=REACTANTS, context=context)
    match = result["result"]["matches"][0]
    assert match["claim_scope"] == "reported_reaction_with_matching_recorded_context"
    assert "not isolated material" in match["scope_limit"]
    assert "source_audited" not in json.dumps(result)
    context["inputs"]["Base"]["components"][0]["amount"]["moles"]["value"] = 0.2
    changed = r.match_product(path, PRODUCT, reactants=REACTANTS, context=context)["result"]["matches"][0]
    assert changed["claim_scope"] == "exact_product_recorded_precedent"


def test_missing_measurement_is_not_zero(source, tmp_path):
    target = synthetic_source(
        tmp_path, source, lambda m: m.outcomes[0].products[3].measurements[0].float_value.ClearField("value")
    )
    path = tmp_path / "bundle"
    r.import_snapshot(target, path)
    rec = r.load_bundle(path)["records"][0]
    measurement = rec["derived"]["projection"]["outcomes"][0]["products"][3]["measurements"][0]
    assert measurement["float_value"] == {}
    assert (
        rec["derived"]["presence"]["/outcomes/0/products/3/measurements/0/float_value/value"]["state"]
        == "absent"
    )
    assert any(e["code"] == "measurement_value_missing_or_negative" for e in rec["derived"]["rejections"])


@pytest.mark.parametrize(
    "measurement,reason",
    [
        ({"type": "YIELD", "percentage": {"value": 101}}, "percentage_out_of_range"),
        ({"type": "SELECTIVITY"}, "measurement_value_missing"),
        ({"type": "AMOUNT", "amount": {"mass": {"value": 1}}}, "invalid_or_unmeasured_product_amount"),
        (
            {"type": "AREA", "float_value": {"value": 1}, "retention_time": {"value": -1, "units": "SECOND"}},
            "invalid_retention_time",
        ),
    ],
)
def test_invalid_measurement_shapes(measurement, reason):
    assert reason in r._measurement_reasons(measurement)


@pytest.mark.parametrize("value", [[], {"routes": "not a list"}, {"routes": [None]}])
def test_malformed_saved_documents_fail_explicitly(bundle, tmp_path, value):
    run = tmp_path / "run"
    make_saved_run(run)
    save_json(run / "routes.json", value)
    with pytest.raises(r.RecordError, match="saved_documents"):
        report.match_saved_route(bundle, run, "synthetic-test-route")


def unsupported_material(kind):
    from rdkit import Chem

    if kind == "polymer":
        mol = Chem.MolFromSmiles("CC")
        group = Chem.CreateMolSubstanceGroup(mol, "SRU")
        group.AddAtomWithIdx(0)
        group.AddAtomWithIdx(1)
        group.SetProp("LABEL", "n")
        group.SetProp("CONNECT", "HT")
        return "MOLBLOCK", Chem.MolToMolBlock(mol)
    if kind == "query_atom":
        return "MOLBLOCK", Chem.MolToMolBlock(Chem.MolFromSmarts("[#6,#7]-C"))
    if kind == "query_bond":
        return "MOLBLOCK", Chem.MolToMolBlock(Chem.MolFromSmarts("[#6]~[#6]"))
    return "SMILES", "C~C"


@pytest.mark.parametrize("kind", ["polymer", "query_atom", "query_bond", "unspecified_bond"])
def test_unsupported_material_cannot_be_exact(kind):
    identifier_type, value = unsupported_material(kind)
    outcome = r.identity([{"type": identifier_type, "value": value}])
    assert outcome["key"] is None and outcome["status"] == "unassessable"
    assert any("unsupported_" in reason for reason in outcome["reasons"])


@pytest.mark.parametrize("kind", ["polymer", "query_atom", "query_bond", "unspecified_bond"])
def test_raw_backed_unsupported_material_never_matches_ethane(source, tmp_path, kind):
    identifier_type, value = unsupported_material(kind)

    def mutate(message):
        product = message.outcomes[0].products[1]
        product.ClearField("identifiers")
        identifier = product.identifiers.add()
        identifier.type = (
            identifier.DESCRIPTOR.fields_by_name["type"].enum_type.values_by_name[identifier_type].number
        )
        identifier.value = value

    target = synthetic_source(tmp_path, source, mutate)
    bundle = tmp_path / "bundle"
    r.import_snapshot(target, bundle)
    record = r.load_bundle(bundle)["records"][0]
    assert record["derived"]["products"][1]["identity"]["key"] is None
    outcome = r.match_product(bundle, "CC")["result"]
    assert not outcome["matches"]
    assert outcome["status"] == "search_incomplete_no_assessable_match"
    assert any(item["path"] == "/outcomes/0/products/1" for item in outcome["unassessable_products"])


def test_concrete_explicit_hydrogens_still_round_trip():
    from rdkit import Chem

    mol = Chem.AddHs(Chem.MolFromSmiles("F[C@H](Cl)Br"))
    identifier = {"type": "MOLBLOCK", "value": Chem.MolToMolBlock(mol)}
    assert r.identity([identifier])["key"] == r.smiles_identity("F[C@H](Cl)Br")["key"]


@pytest.mark.parametrize(
    "mutation", ["swapped_node_ids", "nested_catalog_smiles", "nested_catalog_id", "duplicate_node_id"]
)
def test_all_provided_saved_leaf_layers_must_agree(bundle, tmp_path, mutation):
    run = tmp_path / "run"
    route = make_saved_run(run)
    blocks = route["verification"]["building_blocks"]
    blocks[0]["node_id"], blocks[1]["node_id"] = "a", "b"
    route["synthesis_evidence"] = {
        "route_audit": {
            "building_blocks": [
                {
                    "node_id": "a",
                    "smiles": REACTANTS[0],
                    "catalog_record": {"id": "a", "smiles": REACTANTS[0]},
                },
                {
                    "node_id": "b",
                    "smiles": REACTANTS[1],
                    "catalog_record": {"id": "b", "smiles": REACTANTS[1]},
                },
            ]
        }
    }
    audit = route["synthesis_evidence"]["route_audit"]["building_blocks"]
    if mutation == "swapped_node_ids":
        blocks[0]["node_id"], blocks[1]["node_id"] = "b", "a"
    elif mutation == "nested_catalog_smiles":
        audit[0]["catalog_record"]["smiles"] = "CCC"
    elif mutation == "nested_catalog_id":
        audit[0]["catalog_record"]["id"] = "b"
    else:
        blocks[1]["node_id"] = "a"
    save_json(run / "routes.json", {"routes": [route]})
    with pytest.raises(r.RecordError, match="saved_blocks"):
        report.match_saved_route(bundle, run, "synthetic-test-route")


def test_legacy_leaf_absence_is_visible_but_explicit_empty_is_not_absence(bundle, tmp_path):
    run = tmp_path / "run"
    route = make_saved_run(run)
    route["synthesis_evidence"] = {
        "route_audit": {
            "building_blocks": [{"node_id": "a", "smiles": REACTANTS[0], "catalog_record": {"id": "a"}}]
        }
    }
    save_json(run / "routes.json", {"routes": [route]})
    value = report.match_saved_route(bundle, run, "synthetic-test-route")
    missing = value["saved_run_binding"]["optional_leaf_fields_absent"]
    assert "verification/building_blocks/0/node_id" in missing
    assert "route_audit/building_blocks/a/catalog_record/smiles" in missing
    route["synthesis_evidence"]["route_audit"]["building_blocks"][0]["catalog_record"]["smiles"] = None
    save_json(run / "routes.json", {"routes": [route]})
    with pytest.raises(r.RecordError, match="saved_blocks"):
        report.match_saved_route(bundle, run, "synthetic-test-route")
