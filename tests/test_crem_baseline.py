"""No benchmark chemistry: public-API integration plus adversarial provenance tests."""

import copy
import csv
import json
import os
from pathlib import Path
import sys
import time

import pytest
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from benchmarks.external import crem_baseline as cb
from benchmarks.research.training_library import export_library


@pytest.fixture
def panel(tmp_path):
    source = tmp_path / "source.txt"
    source.write_text("Synthetic engineering fixtures; no biological labels or evidence.\n")
    # Simple rings exercise partial-ring replacement at radius 3 without any
    # benchmark input. Public train query is intentionally retained in the DB.
    corpus = [
        "CC1CCCCC1",
        "CC1CCCC1",
        "CC1CCC1",
        "CC1CCCCCC1",
        "CC1CCOCC1",
        "CC1CCNCC1",
        "C1CCCCC1",
        "C1CCCC1",
        "C1CCOCC1",
        "C1CCNCC1",
    ]
    compounds = [
        {
            "id": f"m{i}",
            "smiles": s,
            "target_id": "training",
            "label": "unlabelled",
            "split": "train",
            "source_id": "fixture",
        }
        for i, s in enumerate(corpus)
    ]
    compounds.append(
        {
            "id": "heldout",
            "smiles": "Cc1nccs1",
            "target_id": "other",
            "label": "unlabelled",
            "split": "test",
            "source_id": "fixture",
        }
    )
    data = {
        "schema_version": 1,
        "scope": "synthetic_fixture",
        "dataset_id": "engineering",
        "sources": [
            {
                "id": "fixture",
                "url": "https://example.invalid/synthetic",
                "artifact": str(source),
                "sha256": cb.digest(source),
            }
        ],
        "compounds": compounds,
        "queries": [{"id": "q", "compound_id": "m0"}],
        "split": {"kind": "target"},
    }
    dataset = tmp_path / "panel.json"
    cb.write_json(dataset, data)
    train = tmp_path / "train.tsv"
    export_library(dataset, train)
    return dataset, train, train.with_suffix(".tsv.provenance.json")


def test_training_inventory_preserves_public_train_queries(panel):
    records, audit = cb.train_inventory(*panel)
    assert len(records) == 10
    assert audit["nontrain_partition_family_overlap"] == 0
    assert audit["public_train_query_overlap"][0]["query_id"] == "q"
    assert audit["public_train_query_overlap"][0]["partition"] == "train"
    assert all("label" not in row for row in records)


def input_hashes(panel):
    return dict(zip(("dataset", "train_tsv", "provenance"), map(cb.digest, panel)))


def rewrite_train(panel, transform):
    dataset, train, provenance = panel
    with train.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows, fields = list(reader), reader.fieldnames
    transform(rows, fields)
    with train.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    report = cb.read_json(provenance)
    report.update(library_sha256=cb.digest(train), records=len(rows))
    cb.write_json(provenance, report)


@pytest.mark.parametrize("what", ["library", "dataset", "source", "provenance"])
def test_input_or_source_hash_tamper_rejected(panel, what):
    dataset, train, provenance = panel
    path = {
        "library": train,
        "dataset": dataset,
        "source": dataset.parent / "source.txt",
        "provenance": provenance,
    }[what]
    if what == "provenance":
        data = cb.read_json(path)
        data["labels_exposed"] = True
        cb.write_json(path, data)
    else:
        with path.open("a") as f:
            f.write("\n")
    with pytest.raises(ValueError):
        cb.train_inventory(*panel)


@pytest.mark.parametrize("replacement", ["Cc1nccs1", "COc1ccccc1"])
def test_heldout_or_nonmember_smiles_rejected_even_with_rehashed_tsv(panel, replacement):
    rewrite_train(panel, lambda rows, fields: rows[0].update(smiles=replacement))
    with pytest.raises(ValueError, match="exact member"):
        cb.train_inventory(*panel)


def test_query_partition_family_overlap_fails_closed(panel):
    dataset, train, provenance = panel
    data = cb.read_json(dataset)
    # Query-only roles may share a family in the dataset loader; the DB guard
    # must still reject it. This also tests stereo-equivalent family leakage.
    data["compounds"].append(
        {
            "id": "private-query",
            "smiles": "[13CH3]C1CCCCC1",
            "target_id": "other",
            "split": "query",
            "label": "unlabelled",
            "source_id": "fixture",
        }
    )
    cb.write_json(dataset, data)
    rewrite_train(panel, lambda rows, fields: [r.update(snapshot=cb.digest(dataset)) for r in rows])
    report = cb.read_json(provenance)
    report["dataset_sha256"] = cb.digest(dataset)
    cb.write_json(provenance, report)
    with pytest.raises(ValueError, match="Nontrain partition family leakage"):
        cb.train_inventory(*panel)


def test_label_column_never_accepted(panel):
    def add(rows, fields):
        fields.append("label")
        for row in rows:
            row["label"] = "active"

    rewrite_train(panel, add)
    with pytest.raises(ValueError, match="unlabelled"):
        cb.train_inventory(*panel)


def test_rehashed_member_subset_is_not_the_expected_complete_library(panel):
    rewrite_train(panel, lambda rows, fields: rows.pop())
    with pytest.raises(ValueError, match="exact complete deterministic"):
        cb.train_inventory(*panel)


def test_parent_input_pins_cannot_be_silently_recomputed(panel, tmp_path):
    expected = input_hashes(panel)
    with panel[1].open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="Parent-declared"):
        cb.prepare(*panel, "unused-upstream", tmp_path / "preparation", expected_hashes=expected)
    assert not (tmp_path / "preparation").exists()


@pytest.mark.parametrize(
    "smiles,reason",
    [
        ("CC.O", "connected"),
        ("[*]C1CCCCC1", "dummy"),
        ("not-a-smiles", "Invalid"),
    ],
)
def test_common_connected_valid_structure_gate(smiles, reason):
    with pytest.raises(ValueError, match=reason):
        cb.normalized(smiles)


def test_map_removal_preserves_defined_stereo_charge_and_descriptor():
    a, ai, _ = cb.normalized("[CH3:1][C@H:2](F)c1cc[nH+]cc1")
    b, bi, _ = cb.normalized("C[C@H](F)c1cc[nH+]cc1")
    assert ai.exact_smiles == bi.exact_smiles
    assert "@" in ai.exact_smiles and "+" in ai.exact_smiles
    assert all(atom.GetAtomMapNum() == 0 for atom in a.GetAtoms())
    assert cb.COUNT_FP.GetSparseCountFingerprint(a) == cb.COUNT_FP.GetSparseCountFingerprint(b)


def test_ranking_is_direct_count_no_chirality_and_family_dedup():
    products = [
        "Cc1ccncc1",
        "C[C@H](F)c1ccncc1",
        "C[C@@H](F)c1ccncc1",
        "CCc1ccccc1",
        "Cc1ccccc1",
        "C1CCCC1",
        "bad",
        "CC(F)c1ccncc1",
    ]
    selected, decisions = cb.select_products("Cc1ccccc1", products)
    direct = rdFingerprintGenerator.GetMorganGenerator(radius=2, includeChirality=False)
    query = direct.GetSparseCountFingerprint(Chem.MolFromSmiles("Cc1ccccc1"))
    assert 1 <= len(selected) <= 3
    for row in selected:
        expected = DataStructs.TanimotoSimilarity(
            query, direct.GetSparseCountFingerprint(Chem.MolFromSmiles(row["smiles"]))
        )
        assert row["score_samples"] == [expected]
        assert "route" not in row
        assert row["route_evidence"]["level"] == "no_route_supplied"
    reasons = {r.get("reason") for r in decisions}
    assert {
        "duplicate_family",
        "same_query_family",
        "unchanged_murcko_scaffold",
    } <= reasons
    second, _ = cb.select_products("Cc1ccccc1", list(reversed(products)))
    assert selected == second


def test_unspecified_stereo_is_recorded_and_not_a_two_d_exclusion():
    mol, identity, state = cb.normalized("CC(F)c1ccccc1")
    assert state["unspecified_stereo"] == 1 and "@" not in identity.exact_smiles
    selected, decisions = cb.select_products("CC(F)c1ccccc1", ["CC(F)c1ccncc1"])
    assert len(selected) == 1 and decisions[0]["unspecified_stereo"] == 1


def test_no_arbitrary_size_or_organic_element_filter():
    assert cb.normalized("C" * 61 + "c1ccccc1")[0].GetNumHeavyAtoms() > 60
    assert cb.normalized("[SiH3]c1ccccc1")[0].HasSubstructMatch(Chem.MolFromSmarts("[#14]"))


def test_isotope_query_state_preserved_but_builder_input_explicitly_excluded():
    mol, identity, _ = cb.normalized("[13CH3]C1CCCCC1")
    assert any(a.GetIsotope() == 13 for a in mol.GetAtoms()) and "13C" in identity.exact_smiles
    with pytest.raises(ValueError, match="isotopic_training_structure"):
        cb.normalized("[13CH3]C1CCCCC1", for_training=True)


def request(tmp_path):
    return {
        "schema_version": 1,
        "kind": "external",
        "method_id": "crem",
        "query": {"id": "q", "smiles": "CC1CCCCC1", "target_id": "fixture"},
        "seed": 42,
        "budget": {"wall_seconds": 360, "max_memory_mb": 24000, "max_candidates": 3, "max_evaluations": 256},
        "assets": {
            "crem_db": str(tmp_path / "db"),
            "crem_manifest": str(tmp_path / "manifest"),
            **{key: str(cb.ROOT / rel) for key, rel in cb.SOURCE_ASSETS.items()},
        },
        "settings": {"policy": cb.POLICY, "manifest_sha256": "0" * 64},
        "output_dir": str(tmp_path),
        "deadline_monotonic": time.monotonic() + 360,
    }


@pytest.mark.parametrize(
    "change",
    [
        lambda r: r.update(labels={"active": ["private"]}),
        lambda r: r["query"].update(references=["private"]),
        lambda r: r["assets"].update(dataset="private"),
        lambda r: r["settings"].update(target_smiles_path="private"),
        lambda r: r.update(seed=True),
        lambda r: r.update(seed=45),
        lambda r: r["budget"].update(max_candidates=4),
        lambda r: r["budget"].update(max_evaluations=257),
        lambda r: r["budget"].update(wall_seconds=float("nan")),
        lambda r: r["budget"].update(max_memory_mb=24001),
    ],
)
def test_private_inputs_and_invalid_budgets_rejected(tmp_path, change):
    value = request(tmp_path)
    change(value)
    with pytest.raises(ValueError):
        cb.validate_request(value)


def test_output_protection_before_work(tmp_path, monkeypatch):
    out = tmp_path / "response.json"
    out.write_text("previous")
    monkeypatch.setattr(cb, "verify_database", lambda r: pytest.fail("must not verify/execute"))
    with pytest.raises(ValueError, match="exists"):
        cb.sample(request(tmp_path), out)
    assert out.read_text() == "previous"


def test_wrong_hashseed_fails_before_official_work(monkeypatch):
    monkeypatch.setenv("PYTHONHASHSEED", "random")
    with pytest.raises(ValueError, match="environment"):
        cb.require_environment()


@pytest.fixture(scope="module")
def official():
    root = os.environ.get("CHEMBOUNCE_CREM_SOURCE")
    if not root:
        pytest.skip(
            "Set CHEMBOUNCE_CREM_SOURCE to the exact b84003c official checkout in a dedicated CPU venv"
        )
    cb.source_snapshot(root)
    cb.runtime_snapshot(root)
    return Path(root).resolve()


@pytest.fixture
def real_prepared(panel, official, tmp_path):
    out = tmp_path / "prepared"
    cb.prepare(*panel, official, out, expected_hashes=input_hashes(panel))
    return out / "plan.json"


def test_actual_official_builder_adapter_roundtrip(real_prepared, tmp_path):
    from benchmarks.research.runner import bounded_method, load_protocol
    from benchmarks.research.protocol import load_dataset

    db_dir = tmp_path / "built"
    manifest = cb.build(real_prepared, db_dir, wall_seconds=60, max_memory_mb=1000)
    assert manifest["official_builder_counters"]["Molecules processed"] == 10
    assert manifest["database"]["ring_mapping_rows"] > 0
    serialized = json.dumps(manifest)
    assert "heldout" not in serialized and "Cc1nccs1" not in serialized
    assert manifest["public_train_query_overlap_count"] == 1
    req = request(tmp_path)
    req["assets"].update(
        {
            "crem_db": str(db_dir / "fragments.db"),
            "crem_manifest": str(db_dir / "database-manifest.json"),
        }
    )
    req["settings"]["manifest_sha256"] = cb.digest(req["assets"]["crem_manifest"])
    req["budget"].update(wall_seconds=30, max_memory_mb=1000)
    method = cb.method_fragment(
        req["assets"]["crem_db"],
        req["assets"]["crem_manifest"],
        req["settings"]["manifest_sha256"],
        tmp_path / "method.json",
    )
    assert set(method["assets"]) == {"crem_db", "crem_manifest", *cb.SOURCE_ASSETS}
    dataset, _ = load_dataset(cb.read_json(real_prepared)["inputs"]["dataset"]["path"])
    protocol = {"schema_version": 1, "budget": req["budget"], "seeds": [42], "methods": [method]}
    protocol_path = tmp_path / "protocol.json"
    cb.write_json(protocol_path, protocol)
    assert load_protocol(protocol_path, dataset)["methods"][0]["resolved_assets"] == req["assets"]
    # Independent runner preflight detects compromised adapter bytes even if
    # that copy removes every in-adapter verification. No source file is edited.
    changed = copy.deepcopy(protocol)
    key = next(k for k, p in cb.SOURCE_ASSETS.items() if p.endswith("crem_baseline.py"))
    compromised = tmp_path / "compromised_adapter.py"
    compromised.write_text("raise RuntimeError('all checks removed')\n")
    changed["methods"][0]["assets"][key]["artifact"] = str(compromised)
    cb.write_json(tmp_path / "tampered-protocol.json", changed)
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        load_protocol(tmp_path / "tampered-protocol.json", dataset)
    reports = [bounded_method(method, copy.deepcopy(req), tmp_path / f"run{i}") for i in range(2)]
    assert all(r["status"] == "success" for r in reports), reports
    assert reports[0]["response"]["candidates"] == reports[1]["response"]["candidates"]
    assert 0 < reports[0]["response"]["evaluations"] <= 256
    assert (tmp_path / "run0" / "response.raw.jsonl").read_bytes() == (
        tmp_path / "run1" / "response.raw.jsonl"
    ).read_bytes()
    assert all(
        cb.chemical_identity(r["smiles"]).scaffold_family_smiles != "C1CCCCC1"
        for r in reports[0]["response"]["candidates"]
    )
    # DB bytes are immutable under generation.
    assert cb.digest(db_dir / "fragments.db") == manifest["db_sha256"]
    with (db_dir / "fragments.db").open("ab") as f:
        f.write(b"tampered")
    with pytest.raises(ValueError, match="database changed"):
        cb.verify_database(req)


def test_real_acyclic_db_cannot_masquerade_as_partial_ring_database(official, tmp_path):
    from crem.db import create_db

    path = tmp_path / "acyclic.db"
    create_db(["CC1CCCCC1", "CC1CCCC1"], path, **(cb.BUILD | {"frag_mode": "acyclic"}))
    with pytest.raises(ValueError, match="ring-fragment"):
        cb.database_summary(path)


def test_real_partial_ring_mode_adds_reachability(official, tmp_path):
    from crem.db import create_db
    from crem.crem import mutate_mol

    path = tmp_path / "ring.db"
    # Disubstitution exercises exo side cuts; monosubstituted whole rings are
    # also reachable through ordinary acyclic cuts and cannot distinguish modes.
    create_db(
        ["CC1CCC(CC)CC1", "CC1CCC(CC)CCC1", "CC1CC(CC)CC1", "CC1CCC(CC)OC1", "CC1CCC(CCC)CC1"],
        path,
        **cb.BUILD,
    )
    mol = Chem.MolFromSmiles("CC1CCC(CC)CC1")
    selected = list(mutate_mol(mol, str(path), seed=42, **cb.MUTATION))
    no_cycles = list(mutate_mol(mol, str(path), seed=42, **(cb.MUTATION | {"replace_cycles": "no"})))
    assert set(selected) - set(no_cycles)
    assert selected == list(mutate_mol(mol, str(path), seed=42, **cb.MUTATION))
    assert all(Chem.MolFromSmiles(s) is not None for s in selected)


def test_actual_official_api_preserves_retained_query_isotope(official, tmp_path):
    from crem.db import create_db
    from crem.crem import mutate_mol

    path = tmp_path / "isotope.db"
    create_db(["CC1CCCCC1", "CC1CCCC1", "CC1CCOCC1"], path, **cb.BUILD)
    mol = Chem.MolFromSmiles("[13CH3]C1CCCCC1")
    before = Chem.MolToSmiles(mol)
    products = list(mutate_mol(mol, str(path), seed=42, **cb.MUTATION))
    assert Chem.MolToSmiles(mol) == before
    retained = [s for s in products if "13C" in s]
    assert retained
    selected, rows = cb.select_products(before, retained)
    assert selected and all("13C" in r["smiles"] for r in selected)


def test_wrong_rdkit_is_rejected_before_source_import(monkeypatch):
    monkeypatch.setattr(cb, "require_environment", lambda: None)
    monkeypatch.setattr(cb.rdBase, "rdkitVersion", "2026.09.1")
    with pytest.raises(ValueError, match="RDKit 2025.09.6"):
        cb.runtime_snapshot("unused")


def test_plan_input_tamper_and_overwrite_rejected(real_prepared, panel, official, tmp_path):
    cb.verify_plan(real_prepared)
    with pytest.raises(ValueError, match="exists"):
        cb.prepare(*panel, official, real_prepared.parent, expected_hashes=input_hashes(panel))
    (real_prepared.parent / "train.smi").write_text("C1CCCCC1\n")
    with pytest.raises(ValueError, match="Artifact changed"):
        cb.verify_plan(real_prepared)


def test_manifest_and_upstream_pin_tamper_fail_closed(tmp_path, monkeypatch):
    req = request(tmp_path)
    manifest = {
        "policy": cb.POLICY,
        "complete": True,
        "build_settings": cb.BUILD,
        "mutation_settings": cb.MUTATION,
        "labels_exposed_to_sampler": False,
        "inventory_counts": {"nontrain_partition_family_overlap": 0},
        "db_sha256": "0" * 64,
    }
    cb.write_json(req["assets"]["crem_manifest"], manifest)
    with pytest.raises(ValueError, match="manifest changed"):
        cb.verify_database(req)
    req["settings"]["manifest_sha256"] = cb.digest(req["assets"]["crem_manifest"])
    manifest["complete"] = False
    cb.write_json(req["assets"]["crem_manifest"], manifest)
    req["settings"]["manifest_sha256"] = cb.digest(req["assets"]["crem_manifest"])
    with pytest.raises(ValueError, match="Incomplete"):
        cb.verify_database(req)


@pytest.mark.parametrize("mode", ["timeout", "error"])
def test_partial_outputs_never_accepted_on_failure(official, tmp_path, monkeypatch, mode):
    import crem.crem

    req = request(tmp_path)
    monkeypatch.setattr(cb, "verify_database", lambda r: {"db_sha256": "0" * 64})

    def interrupted(*args, **kwargs):
        yield "CC1CCCC1"
        raise (
            TimeoutError("fixture timeout")
            if mode == "timeout"
            else RuntimeError("fixture generation failure")
        )

    monkeypatch.setattr(crem.crem, "mutate_mol", interrupted)
    out = tmp_path / "response.json"
    with pytest.raises(RuntimeError, match="internal_timeout|generation_error"):
        cb.sample(req, out)
    result = cb.read_json(out)
    assert result["candidates"] == [] and result["evaluations"] == 1
    assert result["metadata"]["completed_within_budget"] is False
    assert cb.read_json(tmp_path / "response.progress.json")["raw_products"] == 1
    assert "CC1CCCC1" in (tmp_path / "response.raw.jsonl").read_text()


def test_unsupported_structure_counts_without_calling_official_generator(official, tmp_path, monkeypatch):
    import crem.crem

    req = request(tmp_path)
    req["query"]["smiles"] = "CC(F)c1ccccc1.[Na+]"
    monkeypatch.setattr(cb, "verify_database", lambda r: {"db_sha256": "0" * 64})
    monkeypatch.setattr(
        crem.crem, "mutate_mol", lambda *a, **kw: pytest.fail("unsupported query must not generate")
    )
    out = tmp_path / "response.json"
    cb.sample(req, out)
    result = cb.read_json(out)
    assert result["evaluations"] == 0 and result["candidates"] == []
    assert result["metadata"]["mutation_calls"] == 0
    assert result["metadata"]["stop_reason"] == "unsupported_input"


def test_official_source_commit_tamper_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(cb, "git", lambda *a: "not-the-pinned-commit")
    with pytest.raises(ValueError, match="pinned upstream"):
        cb.source_snapshot(tmp_path)


def test_api_yield_cap_closes_generator_without_backfill(official, tmp_path, monkeypatch):
    import crem.crem

    req = request(tmp_path)
    req["budget"]["max_evaluations"] = 2
    monkeypatch.setattr(cb, "verify_database", lambda r: {"db_sha256": "0" * 64})
    seen = []

    def many(*args, **kwargs):
        assert kwargs["max_replacements"] == 2 and kwargs["seed"] == 42
        try:
            for value in ["CC1CCCC1", "CC1CCOCC1", "CC1CCNCC1"]:
                seen.append(value)
                yield value
        finally:
            seen.append("closed")

    monkeypatch.setattr(crem.crem, "mutate_mol", many)
    out = tmp_path / "response.json"
    cb.sample(req, out)
    result = cb.read_json(out)
    assert result["evaluations"] == 2 and result["metadata"]["mutation_calls"] == 1
    assert result["metadata"]["stop_reason"] == "api_product_cap"
    assert result["metadata"]["upstream_internal_attempts"] is None
    assert seen == ["CC1CCCC1", "CC1CCOCC1", "closed"]


def test_cold_start_exhaustion_is_failed_without_inference(official, tmp_path, monkeypatch):
    import crem.crem

    req = request(tmp_path)
    req["deadline_monotonic"] = time.monotonic() - 1
    monkeypatch.setattr(cb, "verify_database", lambda r: {"db_sha256": "0" * 64})
    monkeypatch.setattr(crem.crem, "mutate_mol", lambda *a, **kw: pytest.fail("no budget for generation"))
    out = tmp_path / "response.json"
    with pytest.raises(RuntimeError, match="internal_timeout"):
        cb.sample(req, out)
    result = cb.read_json(out)
    assert result["evaluations"] == 0 and result["candidates"] == []
    assert result["metadata"]["mutation_calls"] == 0


def test_hard_timeout_retains_raw_but_no_accepted_candidate(tmp_path):
    from benchmarks.research.runner import bounded_method

    script = tmp_path / "slow.py"
    script.write_text(
        "import json, pathlib, sys, time\n"
        "p=pathlib.Path(sys.argv[2]); p.with_suffix('.raw.jsonl').write_text('fixture only\\n')\n"
        "p.write_text(json.dumps({'candidates':[{'smiles':'C1CCCC1'}], 'evaluations':1, 'evaluation_unit':'fixture'}))\n"
        "time.sleep(5)\n"
    )
    req = request(tmp_path)
    req["budget"]["wall_seconds"] = 0.4
    method = {"kind": "external", "command": [sys.executable, str(script), "{request}", "{response}"]}
    result = bounded_method(method, req, tmp_path / "run")
    assert result["status"] == "timeout" and result["response"] == {}
    assert (tmp_path / "run" / "response.raw.jsonl").read_text() == "fixture only\n"


def test_builder_rejects_changed_plan_before_import(tmp_path):
    req = {"assets": {"plan": str(tmp_path / "plan.json")}, "settings": {"plan_sha256": "0" * 64}}
    cb.write_json(tmp_path / "plan.json", {})
    cb.write_json(tmp_path / "request.json", req)
    with pytest.raises(ValueError, match="plan changed"):
        cb.build_worker(tmp_path / "request.json", tmp_path / "response.json")
