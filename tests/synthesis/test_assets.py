import os
import pickle
import sys
import types

import pytest

from chembounce_synthesis.assets import convert_catalog, convert_rules, read_index_smiles
from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.reactions import Rules


def test_pickle_cannot_execute_external_reduce(tmp_path):
    marker = tmp_path / "executed"

    class Malicious:
        def __reduce__(self):
            return (os.system, (f"touch {marker}",))

    path = tmp_path / "bad.pkl"
    path.write_bytes(pickle.dumps(Malicious()))
    with pytest.raises(ValueError, match="Expected a SynLlama"):
        list(read_index_smiles(path))
    assert not marker.exists()


def test_inert_index_conversion_without_importing_upstream(tmp_path, monkeypatch):
    # Emulate the documented upstream pickle structure using ephemeral classes;
    # remove the modules before reading to ensure the importer never needs them.
    names = ["synllama", "synllama.chem", "synllama.chem.mol", "synllama.chem.fpindex"]
    for name in names:
        monkeypatch.setitem(sys.modules, name, types.ModuleType(name))
    Mol = type("Molecule", (), {"__module__": "synllama.chem.mol", "__getstate__": lambda self: self.smiles})
    Index = type("FingerprintIndex", (), {"__module__": "synllama.chem.fpindex"})
    sys.modules["synllama.chem.mol"].Molecule = Mol
    sys.modules["synllama.chem.fpindex"].FingerprintIndex = Index
    molecule = Mol()
    molecule.smiles = "CCO"
    data = Index()
    data._molecules = (molecule,)
    path = tmp_path / "index.pkl"
    path.write_bytes(pickle.dumps(data))
    for name in names:
        monkeypatch.delitem(sys.modules, name)
    assert list(read_index_smiles(path)) == ["CCO"]
    output = tmp_path / "catalog.tsv"
    result = convert_catalog([path], output, "published snapshot", "sha256:fixture")
    assert result["records"] == 1
    assert next(iter(Catalog(output).blocks)).startswith("structure-sha256:")


def test_conversion_preserves_supplier_ids_and_rejections(tmp_path):
    path = tmp_path / "supplier.csv"
    path.write_text("sku,structure\nvendor-42,CCO\nvendor-43,notSMILES\nvendor-44,CCO.Cl\n")
    output = tmp_path / "catalog.tsv"
    manifest = convert_catalog([path], output, "Test supplier", "2026-10-11", "structure", "sku")
    assert manifest["records"] == 1 and len(manifest["rejected_records"]) == 2
    assert Catalog(output).get("vendor-42").smiles == "CCO"
    assert Catalog(output).import_provenance["rejected_record_count"] == 2
    with pytest.raises(ValueError, match="overwrite"):
        convert_catalog([path], output, "Test supplier", "2026-10-11", "structure", "sku")
    output.write_text(output.read_text().replace("CCO", "CCN"))
    with pytest.raises(ValueError, match="checksum differs"):
        Catalog(output)


def test_rule_import_records_duplicates(assets, tmp_path):
    smarts = assets[1].get("amide").smarts
    path = tmp_path / "templates.txt"
    path.write_text(smarts + "\n" + smarts + "\n")
    out = tmp_path / "converted.json"
    data = convert_rules(path, out, "fixture")
    assert len(Rules(out).rules) == 1
    assert data["duplicate_source_templates"] == [{"source_line": 2, "same_as_rule": "synllama-rule-001"}]


def test_catalog_missing_provenance_or_conflicting_id(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("id\tsmiles\n1\tCCO\n")
    with pytest.raises(ValueError, match="requires"):
        Catalog(path)
    path.write_text("id\tsmiles\tsource\tsnapshot\n1\tCCO\ttest\tdate\n1\tCCN\ttest\tdate\n")
    with pytest.raises(ValueError, match="Duplicate"):
        Catalog(path)


def test_extension_opcode_is_rejected_before_unpickling(tmp_path):
    path = tmp_path / "extension.pkl"
    path.write_bytes(b"\x80\x02\x82\x01.")
    with pytest.raises(ValueError, match="opcode: EXT1"):
        list(read_index_smiles(path))


def test_search_caches_release_discarded_catalogs_and_rule_sets(assets):
    import gc
    import weakref

    catalog = Catalog(assets[0].path)
    rules = Rules(assets[1].path)
    catalog.compatible(rules.get("amide").roles[0])
    rules.apply("amide", ("CC(=O)O", "CCN"))
    refs = [weakref.ref(catalog), weakref.ref(rules)]
    del catalog, rules
    gc.collect()
    assert [ref() is None for ref in refs] == [True, True]
