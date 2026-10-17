"""Optional byte pins; tiny software fixtures are not supplier or synthesis evidence."""

import hashlib
import json
from pathlib import Path

import pytest

from chembounce_synthesis import catalog as catalog_module, workflow
from chembounce_synthesis.assets import convert_catalog
from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.common import digest
from chembounce_synthesis.continuation import continuation_inputs
from chembounce_synthesis.route_audit import audit_route


def snapshot(tmp_path):
    source = tmp_path / "source.csv"
    source.write_text("id,smiles\na,CCO\nb,CCN\n")
    path = tmp_path / "catalog.tsv"
    convert_catalog([source], path, "fixture", "snapshot")
    sidecar = path.with_suffix(".tsv.manifest.json")
    return path, sidecar, {"catalog_sha256": digest(path), "manifest_sha256": digest(sidecar)}


def config_for(path, binding):
    return {
        "schema_version": 1,
        "catalog": str(path),
        "catalog_binding": binding,
        "rules": "rules.json",
        "provider": {"kind": "saved", "path": "proposals.json"},
    }


@pytest.mark.parametrize(
    "binding",
    [
        None,
        False,
        [],
        "bad",
        {},
        {"catalog_sha256": "0" * 64},
        {"manifest_sha256": "0" * 64},
        {"catalog_sha256": "0" * 64, "manifest_sha256": "0" * 64, "extra": 1},
        *[
            {"catalog_sha256": value, "manifest_sha256": "0" * 64}
            for value in [None, True, 0, "", "a" * 63, "A" * 64, "g" * 64]
        ],
        {"catalog_sha256": "0" * 64, "manifest_sha256": False},
    ],
)
def test_malformed_declaration_rejected_before_any_asset_read(tmp_path, binding, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("asset read during declaration validation")

    monkeypatch.setattr(Path, "read_bytes", forbidden)
    with pytest.raises(ValueError, match="catalog_binding"):
        workflow._config(config_for(tmp_path / "missing.tsv", binding))
    # None is the backward-compatible default only for the direct constructor.
    if binding is not None:
        with pytest.raises(ValueError, match="catalog_binding"):
            Catalog(tmp_path / "missing.tsv", binding=binding)


@pytest.mark.parametrize("change", ["missing", "valid_manifest", "tsv"])
def test_mismatch_rejects_before_molecule_construction(tmp_path, monkeypatch, change):
    path, sidecar, binding = snapshot(tmp_path)
    if change == "missing":
        sidecar.unlink()
    elif change == "valid_manifest":
        data = json.loads(sidecar.read_bytes())
        data["record_origins"]["a"][0]["record_index"] = 9
        sidecar.write_text(json.dumps(data))
        assert Catalog(path).evidence("a")["identifier_provenance"]["source_records"][0]["record_index"] == 9
    else:
        path.write_bytes(path.read_bytes().replace(b"CCO", b"CCC"))
    monkeypatch.setattr(catalog_module, "molecule", lambda *_: pytest.fail("molecule before pin match"))
    with pytest.raises(ValueError, match="catalog_binding"):
        Catalog(path, binding=binding)


def test_matching_pair_preserves_evidence_and_freshness(tmp_path):
    path, sidecar, binding = snapshot(tmp_path)
    catalog = Catalog(path, binding=binding)
    assert catalog.evidence("a") == Catalog(path).evidence("a")
    assert catalog.sha256 == binding["catalog_sha256"]
    assert catalog.import_provenance["sha256"] == binding["manifest_sha256"]
    assert catalog.evidence("a")["supplier_sku_verified"] is False
    assert catalog.evidence("a")["current_stock_verified"] is False
    binding["manifest_sha256"] = "0" * 64
    catalog._check_import_unchanged()
    sidecar.write_bytes(sidecar.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="changed since loading"):
        catalog._check_import_unchanged()


def test_explicit_pins_can_bind_legacy_and_unpinned_absence_stays_legacy(tmp_path):
    path, sidecar, binding = snapshot(tmp_path)
    sidecar.write_text(json.dumps({"output_sha256": digest(path)}))
    binding["manifest_sha256"] = digest(sidecar)
    assert Catalog(path, binding=binding).evidence("a")["identifier_provenance"]["status"] == "unknown_legacy"
    sidecar.unlink()
    assert Catalog(path).evidence("a")["identifier_provenance"]["status"] == "unknown_legacy"


def test_captured_tsv_is_hashed_and_parsed_despite_replace_restore_aba(tmp_path, monkeypatch):
    path, sidecar, binding = snapshot(tmp_path)
    original = path.read_bytes()
    altered = original.replace(b"CCO", b"CCC")
    read_bytes = Path.read_bytes
    parse = catalog_module.molecule
    opened = []

    def capture(candidate):
        data = read_bytes(candidate)
        if candidate == path:
            opened.append(candidate)
            path.write_bytes(altered)
        return data

    def restore_at_parse(smiles):
        # Simulate restoration after an alternate disk read would have occurred.
        path.write_bytes(original)
        return parse(smiles)

    monkeypatch.setattr(Path, "read_bytes", capture)
    monkeypatch.setattr(catalog_module, "molecule", restore_at_parse)
    catalog = Catalog(path, binding=binding)
    assert opened == [path]
    assert catalog.get("a").smiles == "CCO"
    assert catalog.sha256 == hashlib.sha256(original).hexdigest()
    assert read_bytes(path) == original


def test_captured_sidecar_is_hashed_and_parsed_from_same_bytes(tmp_path, monkeypatch):
    path, sidecar, binding = snapshot(tmp_path)
    read_bytes = Path.read_bytes
    original = read_bytes(sidecar)
    changed = json.loads(original)
    changed["record_origins"]["a"][0]["record_index"] = 9

    def capture(candidate):
        data = read_bytes(candidate)
        if candidate == sidecar:
            sidecar.write_text(json.dumps(changed))
        return data

    monkeypatch.setattr(Path, "read_bytes", capture)
    catalog = Catalog(path, binding=binding)
    assert catalog.evidence("a")["identifier_provenance"]["source_records"][0]["record_index"] == 0
    assert catalog._manifest_sha256 == binding["manifest_sha256"]
    with pytest.raises(ValueError, match="changed since loading"):
        catalog._check_import_unchanged()


def test_pinned_csv_retains_bom_crlf_and_quoted_multiline_semantics(tmp_path):
    path = tmp_path / "catalog.tsv"
    path.write_bytes(b'\xef\xbb\xbfid\tsmiles\tsource\tsnapshot\r\na\tCCO\t"two\r\nlines"\tsnapshot\r\n')
    sidecar = path.with_suffix(".tsv.manifest.json")
    sidecar.write_text(json.dumps({"output_sha256": digest(path)}))
    binding = {"catalog_sha256": digest(path), "manifest_sha256": digest(sidecar)}
    assert Catalog(path, binding=binding).evidence("a") == Catalog(path).evidence("a")


def test_workflow_pin_failure_precedes_provider_and_forward(tmp_path, assets, monkeypatch):
    catalog, rules = assets
    config = config_for(catalog.path, {"catalog_sha256": catalog.sha256, "manifest_sha256": "0" * 64})
    config["rules"] = str(rules.path)
    for name in ("SavedProposals", "SynLlama", "_evidence_evaluator"):
        monkeypatch.setattr(workflow, name, lambda *a, **k: pytest.fail("provider/forward before pin check"))
    with pytest.raises(ValueError, match="pinned import manifest"):
        workflow.run("CC(=O)Nc1ccccc1", ["CC(=O)Nc1ccncc1"], config, tmp_path / "run", scoring_mode="2d")


def test_workflow_exports_declaration_and_continuation_preserves_it(tmp_path, assets, proposal):
    catalog, rules = assets
    sidecar = catalog.path.with_suffix(".tsv.manifest.json")
    sidecar.write_text(json.dumps({"output_sha256": catalog.sha256}))
    binding = {"catalog_sha256": catalog.sha256, "manifest_sha256": digest(sidecar)}
    config = config_for(catalog.path, binding)
    config["rules"] = str(rules.path)
    proposals = tmp_path / "proposals.json"
    proposals.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": "CC(=O)Nc1ccncc1", "response": proposal}]})
    )
    config["provider"]["path"] = str(proposals)
    config["max_scored_products"] = 1
    out = tmp_path / "run"
    assert workflow.run("CC(=O)Nc1ccccc1", ["CC(=O)Nc1ccncc1"], config, out, scoring_mode="2d")
    saved = json.loads((out / "manifest.json").read_text())
    assert saved["config"]["catalog_binding"] == binding
    assert saved["catalog_import"]["sha256"] == binding["manifest_sha256"]
    _, _, inherited = continuation_inputs(out)
    assert inherited["catalog_binding"] == binding
    copied = tmp_path / "changed.tsv"
    copied.write_bytes(catalog.path.read_bytes().replace(b"CC(=O)O", b"CCC(=O)O"))
    inherited["catalog"] = str(copied)
    with pytest.raises(ValueError, match="catalog_binding.catalog_sha256"):
        workflow.run(
            "CC(=O)Nc1ccccc1", ["CC(=O)Nc1ccncc1"], inherited, tmp_path / "continued", scoring_mode="2d"
        )


def test_route_audit_keeps_postload_tsv_freshness(tmp_path, assets, route):
    catalog, rules = assets
    sidecar = catalog.path.with_suffix(".tsv.manifest.json")
    sidecar.write_text(json.dumps({"output_sha256": catalog.sha256}))
    bound = Catalog(
        catalog.path, binding={"catalog_sha256": catalog.sha256, "manifest_sha256": digest(sidecar)}
    )
    # Freshness is independent of concrete atom-origin metadata support. Older
    # RDKit correctly reports strict_eligible=False but must still reach this guard.
    audited = audit_route(route, bound, rules)
    assert audited["catalog_sha256"] == bound.sha256
    assert audited["product"] == route["nodes"][-1]["product"]
    catalog.path.write_bytes(catalog.path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="Authoritative catalog/rules files changed since loading"):
        audit_route(route, bound, rules)


def test_malformed_binding_recommendation_preflight_does_no_work(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("preflight performed work before rejecting malformed binding")

    for name in ("Catalog", "Rules", "SynLlama", "_evidence_evaluator", "prepare_reference"):
        monkeypatch.setattr(workflow, name, forbidden)
    config = config_for(tmp_path / "missing.tsv", None)
    with pytest.raises(ValueError, match="catalog_binding"):
        workflow.preflight_recommendation("CC(=O)Nc1ccccc1", config)


def test_load_config_retains_pair_without_content_preflight(tmp_path, monkeypatch):
    binding = {"catalog_sha256": "1" * 64, "manifest_sha256": "2" * 64}
    config = config_for(Path("not-present.tsv"), binding)
    filename = tmp_path / "config.json"
    filename.write_text(json.dumps(config))
    monkeypatch.setattr(Path, "read_bytes", lambda *_: pytest.fail("content loaded by load_config"))
    loaded = workflow.load_config(filename)
    assert loaded["catalog"] == str(tmp_path / "not-present.tsv")
    assert loaded["catalog_binding"] == binding
