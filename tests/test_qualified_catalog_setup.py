"""Tiny full catalogs and simulated orchestration; no inference or route search."""

from dataclasses import FrozenInstanceError
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmarks.research import qualified_collector as collector, qualified_runner as runner
from benchmarks.research.protocol import digest, write_json
from tests.test_qualified_runner import protocol


def pinned(tmp_path, name="catalog", sidecar=False):
    path = tmp_path / (name + ".tsv")
    path.write_text("id\tsmiles\tsource\tsnapshot\nethane\tCC\ttest\tv1\n")
    assets = {"catalog": {"artifact": str(path), "sha256": digest(path)}}
    resolved = {"catalog": str(path)}
    method = {"id": name, "arm": "synllama", "assets": assets, "resolved_assets": resolved}
    if sidecar:
        manifest = path.with_suffix(".tsv.manifest.json")
        write_json(manifest, {"output_sha256": digest(path)})
        assets["catalog_manifest"] = {"artifact": str(manifest), "sha256": digest(manifest)}
        resolved["catalog_manifest"] = str(manifest)
    else:
        method["catalog_manifest_absent"] = True
    root = Path(collector.__file__).resolve().parents[2] / "chembounce_synthesis"
    sources = {str(root / name): digest(root / name) for name in ("catalog.py", "common.py")}
    return method, sources


def load(methods, source, termination=None):
    setup = {"catalogs": []}
    result = runner.preload_catalogs(
        {"methods": methods, "source_sha256": source},
        termination or SimpleNamespace(signal_number=None),
        setup,
    )
    return result, setup


def test_once_per_pinned_path_and_no_cross_campaign_cache(tmp_path, monkeypatch):
    first, source = pinned(tmp_path)
    second = dict(first, id="second")
    third, _ = pinned(tmp_path, "other")
    calls = []
    constructor = collector.Catalog

    def counted(path):
        calls.append(path)
        return constructor(path)

    monkeypatch.setattr(collector, "Catalog", counted)
    bound, setup = load([first, second, third], source)
    assert len(calls) == len(setup["catalogs"]) == 2
    assert bound["catalog"] is bound["second"]
    assert bound["catalog"] is not bound["other"]
    assert setup["catalogs"][0]["method_ids"] == ["catalog", "second"]
    assert all(row["wall_seconds"] >= 0 and row["cpu_seconds"] >= 0 for row in setup["catalogs"])
    again, _ = load([first], source)
    assert len(calls) == 3 and again["catalog"] is not bound["catalog"]
    assert bound["catalog"].evidence("ethane") == constructor(first["resolved_assets"]["catalog"]).evidence(
        "ethane"
    )
    with pytest.raises(FrozenInstanceError):
        bound["catalog"].get("ethane").smiles = "CCC"
    with pytest.raises(TypeError):
        bound["catalog"]._records["ethane"] = None
    assert not hasattr(bound["catalog"], "nearest") and not hasattr(bound["catalog"], "fps")


@pytest.mark.parametrize(
    "mutation",
    [
        "path",
        "content",
        "source",
        "runtime",
        "collector_source",
        "sidecar_added",
        "sidecar_changed",
        "sidecar_removed",
    ],
)
def test_snapshot_mutation_is_rejected(tmp_path, monkeypatch, mutation):
    method, source = pinned(tmp_path, sidecar=mutation in {"sidecar_changed", "sidecar_removed"})
    bound = load([method], source)[0]["catalog"]
    path = Path(method["resolved_assets"]["catalog"])
    expected_source = runner.synthesis_sources(source)
    sidecar = path.with_suffix(".tsv.manifest.json")
    if mutation == "path":
        other = tmp_path / "copy.tsv"
        other.write_bytes(path.read_bytes())
        path = other
    elif mutation == "content":
        path.write_text(path.read_text().replace("CC", "CCC"))
    elif mutation == "source":
        expected_source["common.py"] = "0" * 64
    elif mutation == "runtime":
        monkeypatch.setattr(collector.platform, "python_version", lambda: "different")
    elif mutation == "collector_source":
        real = collector.digest
        monkeypatch.setattr(
            collector, "digest", lambda p: "0" * 64 if str(p) == collector.__file__ else real(p)
        )
    elif mutation == "sidecar_removed":
        sidecar.unlink()
    else:
        sidecar.write_text("{}")
    with pytest.raises(ValueError):
        bound.validate(path, expected_source)


@pytest.mark.parametrize(
    "problem",
    [
        "missing_absence",
        "wrong_sidecar_path",
        "conflicting_absence",
        "wrong_hash",
        "unbound_existing_sidecar",
        "invalid_late_row",
    ],
)
def test_preload_requires_complete_valid_pinned_catalog(tmp_path, problem):
    method, source = pinned(tmp_path, sidecar=problem in {"wrong_sidecar_path", "conflicting_absence"})
    path = Path(method["resolved_assets"]["catalog"])
    if problem == "missing_absence":
        del method["catalog_manifest_absent"]
    elif problem == "wrong_sidecar_path":
        method["resolved_assets"]["catalog_manifest"] = str(tmp_path / "other.json")
    elif problem == "conflicting_absence":
        method["catalog_manifest_absent"] = True
    elif problem == "wrong_hash":
        method["assets"]["catalog"]["sha256"] = "0" * 64
    elif problem == "unbound_existing_sidecar":
        path.with_suffix(".tsv.manifest.json").write_text("{}")
    else:
        path.write_text(path.read_text() + "bad\tnot_smiles\ttest\tv1\n")
        method["assets"]["catalog"]["sha256"] = digest(path)
    with pytest.raises(ValueError):
        load([method], source)


def test_empty_and_cancelled_setup_never_constructs(tmp_path, monkeypatch):
    method, source = pinned(tmp_path)
    monkeypatch.setattr(collector, "Catalog", lambda *args: pytest.fail("Unexpected construction"))
    assert load([], source) == ({}, {"catalogs": []})
    with pytest.raises(InterruptedError):
        load([method], source, SimpleNamespace(signal_number=15))


def test_cancellation_during_construction_finishes_validation_then_stops(tmp_path, monkeypatch):
    method, source = pinned(tmp_path)
    termination = SimpleNamespace(signal_number=None)
    constructor = collector.Catalog

    def cancelled(path):
        catalog = constructor(path)
        termination.signal_number = 15
        return catalog

    monkeypatch.setattr(collector, "Catalog", cancelled)
    setup = {"catalogs": []}
    with pytest.raises(InterruptedError):
        runner.preload_catalogs({"methods": [method], "source_sha256": source}, termination, setup)
    assert setup["catalogs"][0]["status"] == "completed"


@pytest.mark.parametrize("failure", ["invalid", "cancelled"])
def test_setup_failure_preserves_full_denominator_and_costs(tmp_path, monkeypatch, failure):
    p = protocol()
    path, output = tmp_path / "protocol.json", tmp_path / "out"
    write_json(path, p)
    monkeypatch.setattr(runner, "resolve_protocol", lambda *args: p)
    monkeypatch.setattr(runner, "attempt", lambda *args, **kwargs: pytest.fail("No attempt allowed"))

    def failing(*args):
        raise InterruptedError("controlled signal") if failure == "cancelled" else ValueError("bad catalog")

    monkeypatch.setattr(runner, "preload_catalogs", failing)
    if failure == "invalid":
        with pytest.raises(ValueError, match="bad catalog"):
            runner.run(path, output)
    else:
        runner.run(path, output)
    data = runner.read(output / "collection.json")
    setup = runner.read(output / "setup.json")
    assert data["attempts"] == []
    assert setup["status"] == ("failed" if failure == "invalid" else "interrupted")
    assert setup["wall_seconds"] >= 0 and setup["cpu_seconds"] >= 0
    assert setup["accounting_share_per_planned_attempt"]["denominator"] == 10
    assert setup["simultaneous_parent_plus_child_peak_rss_mb"] is None
    assert runner.read(output / "analysis.json")["administrative_setup"] == setup
    assert (output / "collection-binding.json").is_file()


def test_old_protocol_cannot_silently_activate_new_policy(tmp_path):
    p = protocol()
    del p["collector_catalog_policy"]
    with pytest.raises(ValueError, match="explicitly freeze"):
        runner.resolve_protocol(p, tmp_path)


def test_runner_forwards_one_snapshot_to_every_attempt(tmp_path, monkeypatch):
    p = protocol()
    method, source = pinned(tmp_path, sidecar=True)
    p["methods"] = [dict(method, id="a"), dict(method, id="b")]
    p["source_sha256"] = source
    path, output = tmp_path / "protocol.json", tmp_path / "out"
    write_json(path, p)
    monkeypatch.setattr(runner, "resolve_protocol", lambda *args: p)
    received = []

    def attempted(method, query, seed, output, plan, termination, *, bound_catalog):
        received.append(bound_catalog)
        assert (output.parent / "setup.json").is_file()
        return {
            "method_id": method["id"],
            "query_id": query["id"],
            "target_id": query["target_id"],
            "query_sha256": runner.canonical_sha256(query),
            "seed": seed,
            "status": "no_candidates",
            "qualified": [],
            "failures": [],
            "costs": {
                "wall_seconds": 1,
                "cpu_seconds": None,
                "gpu_seconds": None,
                "peak_rss_mb": None,
                "peak_gpu_mb": None,
                "native_units": {},
                "common_units": {},
            },
        }

    monkeypatch.setattr(runner, "attempt", attempted)
    collection = runner.run(path, output)
    assert len(received) == 4 and all(value is received[0] for value in received)
    assert len(collection["administrative_setup"]["catalogs"]) == 1
    assert all(row["costs"]["cpu_seconds"] is None for row in collection["attempts"])
    assert collection["administrative_setup"]["accounting_share_per_planned_attempt"]["denominator"] == 4


def test_mutation_during_constructor_rejected(tmp_path, monkeypatch):
    method, source = pinned(tmp_path, sidecar=True)
    constructor = collector.Catalog

    def changed(path):
        catalog = constructor(path)
        path.with_suffix(".tsv.manifest.json").write_text("{}")
        return catalog

    monkeypatch.setattr(collector, "Catalog", changed)
    with pytest.raises(ValueError, match="sidecar"):
        load([method], source)
