"""Real benchmark functions with inert native generation; no DB/model/chemistry work.

Normal NumPy/RDKit and chembounce initialization occurs during imports, including
its empty-fingerprint convention constant. Synthetic profiles are test-only.
"""

import builtins
import copy

import pytest

import chembounce as core
from benchmarks.research import qualified_adapter as adapter, qualified_runner as runner
from chembounce_synthesis import profiles

pytestmark = pytest.mark.usefixtures("synthetic_benchmark_profile_bindings", "deny_science")


def binding(version="v3"):
    name = {
        "v2": profiles.SELECTIVE_PROFILE,
        "v3": profiles.LANDMARK_SELECTIVE_PROFILE,
        "v4": profiles.ALLOCATION_SELECTIVE_PROFILE,
    }[version]
    return {
        "schema": "qualified-benchmark-protocol-" + version,
        "recommendation_profile": name,
        "profile_artifact_sha256": profiles.PROFILE_SHA256[name],
    }


@pytest.fixture(params=["v3", "v4"])
def fixed_library_version(request):
    return request.param


@pytest.fixture
def deny_science(monkeypatch):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from chembounce_synthesis import conformers, selective_scoring, workflow
    from chembounce_synthesis.catalog import Catalog
    from chembounce_synthesis.reactions import Rules
    from chembounce_synthesis.synllama import SynLlama
    from chembounce_synthesis.forward_support import PinnedReactionT5Forward

    def forbidden(*args, **kwargs):
        pytest.fail("Native binding tests cannot perform scientific work")

    hooks = [
        (Chem, ("MolFromSmiles", "SDMolSupplier")),
        (AllChem, ("EmbedMolecule", "EmbedMultipleConfs")),
        (conformers, ("ensemble", "score_candidate", "_overlay")),
        (selective_scoring, ("ensemble", "score_candidate", "_overlay")),
        (workflow, ("run",)),
        (Catalog, ("__init__",)),
        (Rules, ("__init__",)),
        (SynLlama, ("__init__",)),
        (PinnedReactionT5Forward, ("__init__",)),
        (core, ("chembounce",)),
        (runner, ("bounded_method", "preload_catalogs")),
    ]
    if conformers.rdShapeAlign is not None:
        hooks.append((conformers.rdShapeAlign, ("AlignMol",)))
    for target, names in hooks:
        for name in names:
            monkeypatch.setattr(target, name, forbidden)


def native_request(tmp_path, assets=None, options=None, version="v3"):
    return {
        "query": {"id": "metadata", "target_id": "metadata", "smiles": "NOT_PARSED"},
        "seed": 42,
        "settings": {"arm": "full", "generator_options": options or {}},
        "assets": assets or {},
        "output_dir": str(tmp_path),
        **({"benchmark_protocol": binding(version)} if version else {}),
    }


def native_spy(monkeypatch):
    calls = []

    def generate(hit, **options):
        calls.append((hit, options))
        return {"Final structure": ["INERT_DESIGN"]}, {"test_only": True}

    monkeypatch.setattr(core, "chembounce", generate)
    return calls


@pytest.mark.parametrize(
    "libraries,use_fingerprints,accepted",
    [
        ((), None, False),
        ((), False, False),
        (("scaffold_db",), None, True),
        (("scaffold_db",), True, True),
        (("scaffold_db",), False, True),
        (("fingerprint_db",), None, True),
        (("fingerprint_db",), True, True),
        (("fingerprint_db",), False, False),
        (("scaffold_db", "fingerprint_db"), None, True),
        (("scaffold_db", "fingerprint_db"), True, True),
        (("scaffold_db", "fingerprint_db"), False, True),
    ],
)
def test_v3_native_library_modes_forward_only_explicit_inputs(
    tmp_path, monkeypatch, libraries, use_fingerprints, accepted, fixed_library_version
):
    assets = {key: str(tmp_path / key) for key in libraries}
    options = {} if use_fingerprints is None else {"use_fingerprint_db": use_fingerprints}
    req = native_request(tmp_path, assets, options, version=fixed_library_version)
    calls = native_spy(monkeypatch)
    if not accepted:
        with pytest.raises(ValueError, match="bound scaffold DB or enabled fingerprint DB"):
            adapter.native_designs(req)
        assert calls == []
        return
    adapter.native_designs(req)
    assert len(calls) == 1 and calls[0][0] == "NOT_PARSED"
    actual = calls[0][1]
    for key in ("scaffold_db", "fingerprint_db"):
        assert actual.get(key + "_file") == assets.get(key)
    assert actual["overall_max_n"] == 12 and actual["seed"] == 42
    assert actual.get("use_fingerprint_db", True) == options.get("use_fingerprint_db", True)


@pytest.mark.parametrize(
    "option,value",
    [
        ("output_dir", "unbound"),
        ("seed", 97),
        ("overall_max_n", 20),
        ("scaffold_db_file", "unbound"),
        ("fingerprint_db_file", "unbound"),
        ("replace_scaffold_files", ["unbound-replacements.tsv"]),
        ("fragments_DB", ["INLINE_BOUND_BUT_UNSUPPORTED"]),
        ("use_fingerprint_db", 0),
        ("use_fingerprint_db", "false"),
    ],
)
def test_v3_invalid_options_reject_before_generator_import_or_asset_hash(
    tmp_path, monkeypatch, option, value, fixed_library_version
):
    options = {option: value}
    req = native_request(tmp_path, {"scaffold_db": "bound-scaffolds"}, options, version=fixed_library_version)
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name == "chembounce":
            pytest.fail("Invalid native binding reached generator import")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    with pytest.raises(ValueError, match="V3 native"):
        adapter.native_designs(req)

    def no_hash(path):
        pytest.fail("Invalid library mode reached source/asset hashing")

    monkeypatch.setattr(runner, "digest", no_hash)
    p = {
        **binding(fixed_library_version),
        "budget": runner.BUDGET,
        "collector_catalog_policy": runner.CATALOG_POLICY,
        "source_sha256": {},
        "methods": [
            {
                "id": "full",
                "arm": "full",
                "assets": {"scaffold_db": {"artifact": "inert", "sha256": "a" * 64}},
                "settings": {"generator_options": options},
            }
        ],
    }
    with pytest.raises(ValueError, match="V3 native"):
        runner.resolve_protocol(p, tmp_path)


@pytest.mark.parametrize(
    "assets,options", [({}, {}), ({"fingerprint_db": "bound"}, {"use_fingerprint_db": False})]
)
def test_v3_missing_library_mode_rejects_before_any_hash(
    tmp_path, monkeypatch, assets, options, fixed_library_version
):
    monkeypatch.setattr(runner, "digest", lambda *args: pytest.fail("Unexpected hash"))
    p = {
        **binding(fixed_library_version),
        "budget": runner.BUDGET,
        "collector_catalog_policy": runner.CATALOG_POLICY,
        "source_sha256": {},
        "methods": [
            {"id": "full", "arm": "full", "assets": assets, "settings": {"generator_options": options}}
        ],
    }
    with pytest.raises(ValueError, match="bound scaffold DB or enabled fingerprint DB"):
        runner.resolve_protocol(p, tmp_path)


@pytest.mark.parametrize("version", [None, "v2", "v3", "v4"])
@pytest.mark.parametrize("arm", ["full", "baseline", "crem", "mol2mol", "synllama"])
def test_native_binding_is_v3_full_only(tmp_path, version, arm):
    method = {"id": arm, "arm": arm, "resolved_assets": {}, "native_method": {"asset_names": {}}}
    plan = {"benchmark_protocol": binding(version)} if version else {}
    for phase in ("native", "common"):
        req = runner.request_for(
            method, {"id": "q", "target_id": "t", "smiles": "NOT_PARSED"}, 42, tmp_path, 360, phase, plan
        )
        expected = version is not None and (phase == "common" or (version in ("v3", "v4") and arm == "full"))
        assert ("benchmark_protocol" in req) == expected
        if expected:
            assert req["benchmark_protocol"] == binding(version)


@pytest.mark.parametrize(
    "attack",
    [
        None,
        "missing_binding",
        "changed_binding",
        "changed_scaffold",
        "changed_fp",
        "missing_scaffold",
        "mode",
        "options",
        "arm",
    ],
)
def test_v3_native_supply_rechecks_cross_phase_authority(tmp_path, attack, fixed_library_version):
    native = tmp_path / "native"
    native.mkdir()
    req = native_request(
        tmp_path,
        {"scaffold_db": "bound-scaffolds", "fingerprint_db": "bound-fp"},
        version=fixed_library_version,
    )
    prior = copy.deepcopy(req)
    if attack == "missing_binding":
        prior.pop("benchmark_protocol")
    elif attack == "changed_binding":
        prior["benchmark_protocol"]["profile_artifact_sha256"] = "0" * 64
    elif attack in ("changed_scaffold", "changed_fp"):
        prior["assets"]["scaffold_db" if attack == "changed_scaffold" else "fingerprint_db"] = "other"
    elif attack == "missing_scaffold":
        prior["assets"].pop("scaffold_db")
    elif attack == "mode":
        prior["settings"]["generator_options"]["use_fingerprint_db"] = False
    elif attack == "options":
        prior["settings"]["generator_options"]["replace_scaffold_files"] = ["unbound"]
    elif attack == "arm":
        prior["settings"]["arm"] = "synllama"
    adapter.write_json(native / "request.json", prior)
    adapter.write_json(native / "response.json", {"candidates": [{"smiles": "INERT_DESIGN"}]})
    adapter.write_json(native / "design-pool.json", {"designs": ["INERT_DESIGN"]})
    # Rehash the forged request so rejection must come from cross-phase authority.
    req["native"] = {"directory": str(native), "files_sha256": runner.file_manifest(native)}
    if attack:
        with pytest.raises(ValueError, match="common request authority"):
            adapter.native_supply(req, tmp_path)
    else:
        designs, bridge, counts = adapter.native_supply(req, tmp_path)
        assert designs == ["INERT_DESIGN"] and bridge is None and counts["raw_count"] == 1


@pytest.mark.parametrize("version", [None, "v2"])
def test_historical_native_paths_and_unbound_envelope_remain_unchanged(tmp_path, monkeypatch, version):
    # Legacy native requests have no benchmark field even when common uses v2.
    calls = native_spy(monkeypatch)
    req = native_request(tmp_path, options={"replace_scaffold_files": ["historical-setting"]}, version=None)
    adapter.native_designs(req)
    assert calls[0][1]["replace_scaffold_files"] == ["historical-setting"]
    native = tmp_path / "native"
    native.mkdir()
    adapter.write_json(native / "request.json", req)
    adapter.write_json(native / "response.json", {"candidates": [{"smiles": "INERT_DESIGN"}]})
    adapter.write_json(native / "design-pool.json", {"designs": ["INERT_DESIGN"]})
    common = copy.deepcopy(req)
    if version:
        common["benchmark_protocol"] = binding(version)
    common["native"] = {"directory": str(native), "files_sha256": runner.file_manifest(native)}
    assert adapter.native_supply(common, tmp_path)[0] == ["INERT_DESIGN"]
