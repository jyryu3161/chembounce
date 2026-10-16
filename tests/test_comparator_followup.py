"""Model-free protocol integration fixtures; no real panel or model is executed."""

import copy
import json
from pathlib import Path
import sys
import subprocess

import pytest

from benchmarks.research import prepare_comparator_followup as prep
from benchmarks.research.runner import load_protocol, public_request


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(prep.encode(value))
    return {"artifact": str(path), "sha256": prep.digest(path)}


def tiny_upstream(root, namespace):
    (root / namespace).mkdir(parents=True)
    (root / namespace / "__init__.py").write_text("# committed synthetic upstream fixture\n")
    (root / "LICENSE.txt").write_text("synthetic fixture\n")
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    subprocess.run(["git", "-C", str(root), "add", "."], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(root),
            "-c",
            "user.name=Fixture",
            "-c",
            "user.email=fixture@example.invalid",
            "commit",
            "-qm",
            "fixture",
        ],
        check=True,
    )
    commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    return commit, {
        str(p.relative_to(root)): prep.digest(p)
        for p in (root / namespace / "__init__.py", root / "LICENSE.txt")
    }


@pytest.fixture
def study(tmp_path, monkeypatch):
    base = tmp_path / "inputs"
    root = tmp_path / "tested-source"
    root.mkdir()
    python = str(Path(sys.executable).absolute())
    monkeypatch.setattr(
        prep, "source_snapshot", lambda r, c, p: {"root": str(r), "commit": c, "files_sha256": {}}
    )
    monkeypatch.setattr(
        prep,
        "runtime_snapshot",
        lambda py, p, rdkit_version="2025.9.6", import_crem_from=None: {
            "python": py,
            "packages": {"rdkit": rdkit_version},
            "rdkit_policy": rdkit_version,
            "metadata_probe_crem_root": str(import_crem_from) if import_crem_from else None,
            "version": "synthetic Python",
            "prefix": str(Path(py).parent.parent),
        },
    )
    fp = base / "fingerprints"
    fp.mkdir(parents=True)
    for name in ("fingerprints.npy", "metadata.json", "smiles.npy"):
        (fp / name).write_text(name)
    fixture_model = base / "medium.prior"
    fixture_model.write_text("synthetic bytes; not a model")
    crem_upstream, mol_upstream = base / "crem-upstream", base / "mol-upstream"
    crem_commit, crem_files = tiny_upstream(crem_upstream, "crem")
    mol_commit, mol_files = tiny_upstream(mol_upstream, "reinvent")
    crem_python, sampler_python = base / "crem-venv/bin/python", base / "sampler-venv/bin/python"
    for interpreter in (crem_python, sampler_python):
        interpreter.parent.mkdir(parents=True)
        interpreter.write_text("synthetic non-executable interpreter marker")
    decision = {
        "source_root": str(tmp_path / "original"),
        "panels": {},
        "crem": {
            "upstream_commit": crem_commit,
            "mutation_proposal": {
                "min_size": 1,
                "max_size": 15,
                "min_rel": 0,
                "max_rel": 1,
                "min_inc": -2,
                "max_inc": 2,
                "max_replacements": 256,
                "replace_cycles": "partial_exo",
                "min_freq": 1,
                "set_names": ["train"],
                "ncores": 1,
                "discard_ring_geometry": True,
            },
        },
        "mol2mol": {"upstream_commit": mol_commit, "official_sha512": prep.digest(fixture_model, "sha512")},
    }
    protocols = {}
    for split in prep.SPLITS:
        panel = {
            "compounds": [
                {
                    "id": f"compound-{i}",
                    "smiles": "c1ccccc1",
                    "target_id": f"target-{i // 2}",
                    "label": "active",
                    "activity": {"value": 7},
                    "split": "test",
                }
                for i in range(prep.COUNTS[split])
            ],
            "queries": [
                {"id": f"query-{i}", "compound_id": f"compound-{i}"} for i in range(prep.COUNTS[split])
            ],
        }
        panel_path = base / f"{split}-panel.json"
        panel_spec = save(panel_path, panel)
        decision["panels"][split] = {"path": panel_path.name, "sha256": panel_spec["sha256"]}
        original_path = base / f"{split}-original.json"
        original = {
            "seeds": prep.SEEDS,
            "methods": [
                {
                    "id": "pipeline",
                    "kind": "pipeline",
                    "budget_overrides": {"max_evaluations": 12},
                    "assets": {},
                    "settings": {"generator_options": {"search_policy": "diverse", "scaffold_top_n": 8}},
                }
            ],
        }
        original_spec = save(original_path, original)
        protocols[split] = (original_path.name, original_spec["sha256"])
    monkeypatch.setattr(prep, "ORIGINAL_PROTOCOLS", protocols)
    decision_path = base / "decision.json"
    decision_sha = save(decision_path, decision)["sha256"]
    monkeypatch.setattr(prep, "DECISION_SHA", decision_sha)
    amendment = {
        "previous_decision": {"sha256": decision_sha},
        "supplementary_grid": {"methods": list(prep.METHODS), "seeds": prep.SEEDS},
        "current_pipeline": {
            "source_worktree": str(root),
            "source_commit": "tested-commit",
            "runtime": python,
            "fingerprint_db": {
                "path": str(fp),
                "files_sha256": {p.name: prep.digest(p) for p in fp.iterdir()},
            },
        },
        "execution_policy": {
            "single_gpu": "GPU-physical-zero",
            "protected_resources": [str(tmp_path / "protected")],
        },
        "required_before_launch": ["independent parent assessment and activation record"],
    }
    amendment_path = base / "amendment.json"
    monkeypatch.setattr(prep, "AMENDMENT_SHA", save(amendment_path, amendment)["sha256"])
    crem_root, mol_root = base / "crem-source", base / "mol-source"
    for source_root, relative_files in (
        (crem_root, prep.CREM_FILES),
        (mol_root, prep.MOL2MOL_FILES.values()),
    ):
        for relative in relative_files:
            path = source_root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# synthetic source only\n")
    fragments = {prep.METHODS[1]: {}}
    for split in prep.SPLITS:
        db = base / f"{split}.db"
        db.write_text("synthetic DB marker; not a database")
        mutation = copy.deepcopy(decision["crem"]["mutation_proposal"])
        mutation["min_rel_size"] = mutation.pop("min_rel")
        mutation["max_rel_size"] = mutation.pop("max_rel")
        mutation.update(
            radius=3, symmetry_fixes=False, return_rxn=False, return_rxn_freq=False, return_mol=False
        )
        manifest = {
            "complete": True,
            "labels_exposed_to_sampler": False,
            "input_hashes": {"dataset": decision["panels"][split]["sha256"]},
            "inventory_counts": {"nontrain_partition_family_overlap": 0},
            "db_sha256": prep.digest(db),
            "mutation_settings": mutation,
            "upstream": {
                "commit": crem_commit,
                "root": str(crem_upstream),
                "tracked_files_sha256": crem_files,
                "generated_files_sha256": {},
            },
            "implementation": {
                "root": str(crem_root),
                "files_sha256": {p: prep.digest(crem_root / p) for p in prep.CREM_FILES},
            },
            "runtime": {
                "executable": str(crem_python),
                "python": "synthetic Python",
                "prefix": str(crem_python.parent.parent),
                "distributions": {"rdkit": "2025.9.6"},
                "origins": {},
            },
            "build_settings": {
                "set_name": "train",
                "radii": [3],
                "ncpu": 1,
                "max_heavy_atoms": 15,
                "keep_stereo": True,
                "frag_mode": "both_optimal",
            },
        }
        manifest_spec = save(base / f"{split}-manifest.json", manifest)
        assets = {"crem_db": {"artifact": str(db), "sha256": prep.digest(db)}, "crem_manifest": manifest_spec}
        assets.update(
            {
                "source__" + p.replace("/", "__").replace(".", "_"): {
                    "artifact": str(crem_root / p),
                    "sha256": prep.digest(crem_root / p),
                }
                for p in prep.CREM_FILES
            }
        )
        method = {
            "id": prep.METHODS[1],
            "kind": "external",
            "command": [
                str(crem_python),
                str(crem_root / prep.CREM_FILES[0]),
                "sample",
                "--request",
                "{request}",
                "--response",
                "{response}",
            ],
            "settings": {
                "policy": "crem-train-only-partial-exo-v1",
                "manifest_sha256": manifest_spec["sha256"],
            },
            "budget_overrides": {"max_evaluations": 256},
            "assets": assets,
        }
        fragments[prep.METHODS[1]][split] = save(base / f"{split}-crem.json", method)
    bundle = {
        "upstream_commit": mol_commit,
        "model_sha256": prep.digest(fixture_model),
        "model_sha512": prep.digest(fixture_model, "sha512"),
        "common_python": python,
        "upstream": str(mol_upstream),
        "source_sha256": mol_files,
        "adapter_root": str(mol_root / "benchmarks/external/mol2mol"),
        "common_root": str(mol_root),
        "adapter_sha256": {
            Path(p).name: prep.digest(mol_root / p)
            for k, p in prep.MOL2MOL_FILES.items()
            if k.startswith("adapter_")
        },
        "common_sha256": {
            p: prep.digest(mol_root / p) for k, p in prep.MOL2MOL_FILES.items() if k.startswith("common_")
        },
        "common_versions": {"rdkit": "2025.9.6"},
        "common_python_sha256": prep.digest(python),
        "sampler_python": str(sampler_python),
        "sampler_versions": {"rdkit": "2023.9.5"},
        "sampler_python_sha256": prep.digest(sampler_python),
    }
    bundle_spec = save(base / "bundle.json", bundle)
    mol_assets = {
        "mol2mol_bundle": bundle_spec,
        "mol2mol_model": {"artifact": str(fixture_model), "sha256": prep.digest(fixture_model)},
    }
    mol_assets.update(
        {
            k: {"artifact": str(mol_root / p), "sha256": prep.digest(mol_root / p)}
            for k, p in prep.MOL2MOL_FILES.items()
        }
    )
    method = {
        "id": prep.METHODS[2],
        "kind": "external",
        "command": [
            python,
            str(mol_root / "benchmarks/external/mol2mol/adapter.py"),
            "--request",
            "{request}",
            "--response",
            "{response}",
        ],
        "settings": prep.MOL2MOL_SETTINGS,
        "budget_overrides": {"max_evaluations": 256},
        "assets": mol_assets,
    }
    fragments[prep.METHODS[2]] = save(base / "mol2mol.json", method)
    fragments_path = base / "fragments.json"
    save(fragments_path, fragments)
    return {
        "decision": decision_path,
        "amendment": amendment_path,
        "fragments": fragments_path,
        "output": tmp_path / "prepared",
        "base": base,
    }


def execute(study):
    return prep.prepare(*(study[k] for k in ("decision", "amendment", "fragments", "output")))


def edit_method(study, change, *, method_id=prep.METHODS[1], split="scaffold"):
    fragments = prep.read_json(study["fragments"])
    old = fragments[method_id][split] if method_id == prep.METHODS[1] else fragments[method_id]
    path = Path(old["artifact"])
    method = prep.read_json(path)
    change(method)
    updated = save(path, method)
    if method_id == prep.METHODS[1]:
        fragments[method_id][split] = updated
    else:
        fragments[method_id] = updated
    save(study["fragments"], fragments)


def test_exact_grid_commands_and_immutability(study):
    manifest = execute(study)
    assert manifest["execution_authorized"] is False
    plan = prep.read_json(study["output"] / "execution-plan.json")
    grid = prep.read_json(study["output"] / "grid.json")
    assert len(grid) == len({(r["split"], r["method_id"], r["query_id"], r["seed"]) for r in grid}) == 270
    assert len(plan["commands"]) == 9
    for method in prep.METHODS:
        assert sum(r["method_id"] == method for r in grid) == 90
    for command in plan["commands"]:
        env = command["environment_overrides"]
        assert env["CUDA_VISIBLE_DEVICES"] == (
            "" if command["method_id"] == prep.METHODS[1] else "GPU-physical-zero"
        )
        assert env["OMP_NUM_THREADS"] == env["OPENBLAS_NUM_THREADS"] == "1"
        assert env["PYTHONPATH"] == command["cwd"]
        assert command["no_retry_or_resume"]
        assert not Path(command["argv"][-1]).exists()
        assert not Path(command["argv"][-1]).is_relative_to(study["output"])
    assert prep.verify_prepared(study["output"], prep.digest(study["output"] / "manifest.json")) == manifest
    with pytest.raises(ValueError, match="immutable"):
        execute(study)


def test_original_pipeline_preserved_except_identity_and_index(study):
    execute(study)
    for split in prep.SPLITS:
        original = prep.read_json(study["base"] / f"{split}-original.json")["methods"][0]
        protocol = prep.read_json(study["output"] / f"protocols/{split}-pipeline_current.json")
        pipeline = protocol["methods"][0]
        assert pipeline["settings"] == original["settings"]
        assert pipeline["budget_overrides"] == {"max_evaluations": 12}
        assert protocol["budget"] == prep.BUDGET
        assert pipeline["assets"]["fingerprint_db"]["directory"] == str(study["base"] / "fingerprints")


def test_private_labels_do_not_reach_requests(study):
    execute(study)
    dataset = prep.read_json(study["base"] / "scaffold-panel.json")
    for method_id in prep.METHODS:
        method = prep.read_json(study["output"] / f"protocols/scaffold-{method_id}.json")["methods"][0]
        method["resolved_assets"] = {
            k: v.get("artifact", v.get("directory")) for k, v in method["assets"].items()
        }
        request = public_request(
            dataset, dataset["queries"][0], method, 42, prep.BUDGET, study["output"] / "unused"
        )
        assert set(request["query"]) == {"id", "smiles", "target_id"}
        assert set(request) == {
            "schema_version",
            "kind",
            "method_id",
            "query",
            "seed",
            "budget",
            "assets",
            "settings",
            "output_dir",
        }
        serialized = json.dumps(request)
        assert "scaffold-panel.json" not in serialized
        assert '"label"' not in serialized and '"activity"' not in serialized
        assert request["budget"]["max_evaluations"] == (12 if method_id == prep.METHODS[0] else 256)
    public = prep.read_json(study["output"] / "public-queries/scaffold.json")
    assert all(set(q) == {"id", "smiles", "target_id"} for q in public)


@pytest.mark.parametrize(
    "change",
    [
        lambda m: m.update(labels={"query": "active"}),
        lambda m: m["assets"].update(reference_panel={"artifact": "private.json", "sha256": "a" * 64}),
        lambda m: m["assets"].update(source__private_labels={"artifact": "private.json", "sha256": "a" * 64}),
        lambda m: m["command"].append("--held-out-labels"),
        lambda m: m.update(budget_overrides={"max_evaluations": 512}),
    ],
)
def test_reject_private_or_undeclared_fragment_fields(study, change):
    edit_method(study, change)
    with pytest.raises(ValueError):
        execute(study)
    assert not study["output"].exists()


def test_reject_changed_sampler_settings(study):
    edit_method(study, lambda m: m["settings"].update(temperature=2), method_id=prep.METHODS[2])
    with pytest.raises(ValueError, match="settings were changed"):
        execute(study)
    assert not study["output"].exists()


def test_missing_final_db_cannot_create_fake_protocol(study):
    (study["base"] / "scaffold.db").unlink()
    with pytest.raises(ValueError, match="finalized file is absent"):
        execute(study)
    assert not study["output"].exists()


def test_cross_split_db_rejected(study):
    fragments = prep.read_json(study["fragments"])
    fragments[prep.METHODS[1]]["time"] = fragments[prep.METHODS[1]]["scaffold"]
    # Fixture panels have identical SMILES but split-specific file bytes must differ.
    time_panel = prep.read_json(study["base"] / "time-panel.json")
    time_panel["split_marker"] = "time"
    save(study["base"] / "time-panel.json", time_panel)
    # Test the isolated fragment check with a different evaluator panel digest.
    with pytest.raises(ValueError, match="does not belong"):
        prep.external_method(
            fragments[prep.METHODS[1]]["time"],
            study["base"],
            prep.METHODS[1],
            "time",
            prep.read_json(study["decision"]),
            "b" * 64,
            prep.Pins(),
        )


def test_input_hash_change_rejected(study):
    (study["base"] / "scaffold-panel.json").write_text("{}")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        execute(study)
    assert not study["output"].exists()


def test_output_and_asset_mutations_detected(study):
    execute(study)
    output = study["output"] / "grid.json"
    original = output.read_text()
    output.write_text("[]")
    with pytest.raises(ValueError, match="Prepared output was modified"):
        prep.verify_prepared(study["output"], prep.digest(study["output"] / "manifest.json"))
    output.write_text(original)
    (study["base"] / "scaffold.db").write_text("changed")
    with pytest.raises(ValueError, match="Prepared input changed"):
        prep.verify_prepared(study["output"], prep.digest(study["output"] / "manifest.json"))


def test_directory_inventory_addition_detected(study):
    execute(study)
    (study["base"] / "fingerprints" / "extra.npy").write_text("new")
    with pytest.raises(ValueError, match="directory asset changed"):
        prep.verify_prepared(study["output"], prep.digest(study["output"] / "manifest.json"))


def test_analysis_plan_retains_failures_and_separate_common_3d(study):
    execute(study)
    analysis = prep.read_json(study["output"] / "analysis-plan.json")
    scoring = prep.read_json(study["output"] / "common-3d-plan.json")
    assert analysis["split_attempts"] == {"scaffold": 108, "time": 108, "target": 54}
    assert "original 630 remain separate" in analysis["design"]
    assert "paired method differences at identical query/seed slots" in analysis["hierarchy"]
    assert "zero partial candidates" in analysis["failures"]
    assert "No route supplied is unassessed" in analysis["route_policy"]
    assert "selection bias" in analysis["route_policy"]
    assert scoring["basic_runner_computes_3d"] is False
    assert scoring["conformers"] == prep.CONFORMERS
    assert scoring["seeds"] == [42, 43, 44]


def test_duplicate_json_keys_rejected(tmp_path):
    path = tmp_path / "duplicate.json"
    path.write_text('{"seeds": [42], "seeds": [43]}')
    with pytest.raises(ValueError, match="Duplicate JSON key"):
        prep.read_json(path)


def test_no_run_directory_reuse(study):
    study["output"].with_name(study["output"].name + "-runs").mkdir()
    with pytest.raises(ValueError, match="no reuse"):
        execute(study)
    assert not study["output"].exists()


def test_external_protocols_load_in_existing_runner(study):
    execute(study)
    dataset = prep.read_json(study["base"] / "scaffold-panel.json")
    for method_id in prep.METHODS[1:]:
        protocol = load_protocol(study["output"] / f"protocols/scaffold-{method_id}.json", dataset)
        method = protocol["methods"][0]
        request = public_request(
            dataset, dataset["queries"][0], method, 42, protocol["budget"], study["output"] / "unused"
        )
        assert request["budget"] == prep.BUDGET
        assert set(request["query"]) == {"id", "smiles", "target_id"}
        assert method["overlap_audit"][0]["status"] == "unknown"


def test_protected_output_and_parent_symlink(study):
    decision, amendment = prep.read_json(study["decision"]), prep.read_json(study["amendment"])
    root = Path(decision["source_root"])
    root.mkdir()
    alias = study["base"] / "root-alias"
    alias.symlink_to(root, target_is_directory=True)
    for output in (
        root / "source/new-plan",
        alias / "source/new-plan",
        Path(amendment["current_pipeline"]["source_worktree"]) / "plan",
    ):
        with pytest.raises(ValueError, match="restricted|protected"):
            prep.validate_output_location(output, decision, amendment)
    allowed = alias / "benchmark-results/research-validation/new-plan"
    assert (
        prep.validate_output_location(allowed, decision, amendment)
        == root / "benchmark-results/research-validation/new-plan"
    )


def test_untracked_upstream_source_cannot_be_ignored(study):
    (study["base"] / "mol-upstream/reinvent/untracked.py").write_text("# shadow source")
    with pytest.raises(ValueError, match="Unpinned upstream Python"):
        execute(study)
    assert not study["output"].exists()


def test_nested_runtime_and_source_inventory_are_pinned(study):
    execute(study)
    plan = prep.read_json(study["output"] / "execution-plan.json")
    sampler_python = str(study["base"] / "sampler-venv/bin/python")
    assert plan["runtimes"][sampler_python]["rdkit_policy"] == "2023.9.5"
    assert sampler_python in plan["input_files_sha256"]
    assert str(study["base"] / "mol-upstream/reinvent/__init__.py") in plan["input_files_sha256"]
    assert len(plan["upstream_commits"]) == 2
    anchor = prep.digest(study["output"] / "manifest.json")
    (study["base"] / "sampler-venv/bin/python").write_text("changed sampler")
    with pytest.raises(ValueError, match="Prepared input changed"):
        prep.verify_prepared(study["output"], anchor)


def test_declared_runtime_distribution_tampering_rejected(study):
    def change(method):
        path = Path(method["assets"]["mol2mol_bundle"]["artifact"])
        bundle = prep.read_json(path)
        bundle["sampler_versions"]["rdkit"] = "wrong"
        method["assets"]["mol2mol_bundle"] = save(path, bundle)

    edit_method(study, change, method_id=prep.METHODS[2])
    with pytest.raises(ValueError, match="distribution map differs"):
        execute(study)
    assert not study["output"].exists()


def test_source_inventory_escape_rejected(tmp_path):
    with pytest.raises(ValueError, match="escapes"):
        prep.pin_inventory(tmp_path, {"../outside.py": "a" * 64}, prep.Pins())


def test_manifest_anchor_catches_coordinated_output_rewrite(study):
    execute(study)
    anchor = prep.digest(study["output"] / "manifest.json")
    manifest = prep.read_json(study["output"] / "manifest.json")
    grid = study["output"] / "grid.json"
    grid.write_text("[]")
    manifest["files_sha256"]["grid.json"] = prep.digest(grid)
    save(study["output"] / "manifest.json", manifest)
    with pytest.raises(ValueError, match="externally retained"):
        prep.verify_prepared(study["output"], anchor)
