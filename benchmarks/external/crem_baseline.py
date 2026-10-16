"""Pinned CReM baseline and separately audited train-only fragment DB preparation.

This module never installs software or downloads data. `prepare` can read private
panel roles; `sample` accepts only a public query and two preverified DB assets.
"""

import argparse
from collections import Counter
import csv
from dataclasses import asdict
import hashlib
import importlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import re
import resource
import sqlite3
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator

from chembounce_synthesis.common import molecule
from chembounce_synthesis.identity import chemical_identity, IDENTITY_POLICY
from fingerprint_db import stable_canonical_molecule


UPSTREAM = "b84003cd8c84bc2afe98b5473570b7328e703ead"
UPSTREAM_URL = "https://github.com/DrrDom/crem"
POLICY = "crem-train-only-partial-exo-v1"
BUILD = {
    "set_name": "train",
    "radii": [3],
    "ncpu": 1,
    "max_heavy_atoms": 15,
    "keep_stereo": True,
    "mode": 0,
    "frag_mode": "both_optimal",
    "chunk_size": 100,
    "flush_every": 1,
    "parallel_shards": 1,
    "prefetch": 1,
    "verbose": True,
    "log_every": 1,
    "timings": True,
    "fragment_error_log": True,
}
MUTATION = {
    "radius": 3,
    "min_size": 1,
    "max_size": 15,
    "min_rel_size": 0,
    "max_rel_size": 1,
    "min_inc": -2,
    "max_inc": 2,
    "max_replacements": 256,
    "replace_cycles": "partial_exo",
    "symmetry_fixes": False,
    "min_freq": 1,
    "set_names": ["train"],
    "return_rxn": False,
    "return_rxn_freq": False,
    "return_mol": False,
    "ncores": 1,
    "discard_ring_geometry": True,
}
ENV = {
    "CUDA_VISIBLE_DEVICES": "",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0",
}
PROJECT_FILES = (
    "benchmarks/external/crem_baseline.py",
    "chembounce_synthesis/common.py",
    "chembounce_synthesis/identity.py",
    "fingerprint_db.py",
    "benchmarks/research/protocol.py",
    "benchmarks/research/runner.py",
    "benchmarks/research/evaluation.py",
    "chemistry.py",
    "chembounce_synthesis/exploration.py",
    "chembounce_synthesis/__init__.py",
    "benchmarks/research/__init__.py",
)
SOURCE_ASSETS = {"source__" + p.replace("/", "__").replace(".", "_"): p for p in PROJECT_FILES}
COUNT_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, includeChirality=False)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_hash(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def read_json(path):
    def bad(value):
        raise ValueError(f"Nonfinite JSON: {value}")

    return json.loads(Path(path).read_text(), parse_constant=bad)


def write_json(path, value):
    """Atomic replacement for this invocation's progress; final outputs use fresh dirs."""
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def artifact(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": digest(path), "bytes": path.stat().st_size}


def verify_artifact(spec):
    path = Path(spec["path"])
    if digest(path) != spec["sha256"] or path.stat().st_size != spec["bytes"]:
        raise ValueError(f"Artifact changed: {path}")
    return path


def require_environment():
    if any(os.environ.get(k) != v for k, v in ENV.items()):
        raise ValueError("Start Python with the declared CPU/thread/PYTHONHASHSEED environment")
    if sys.prefix == sys.base_prefix:
        raise ValueError("A dedicated isolated virtual environment is required")
    config = Path(sys.prefix) / "pyvenv.cfg"
    if not config.is_file() or not re.search(r"include-system-site-packages\s*=\s*false", config.read_text()):
        raise ValueError("Virtual environment must not include system site packages")


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def source_snapshot(upstream):
    upstream = Path(upstream).resolve()
    if git(upstream, "rev-parse", "HEAD") != UPSTREAM:
        raise ValueError("CReM checkout is not the pinned upstream commit")
    if git(upstream, "status", "--porcelain", "--untracked-files=no"):
        raise ValueError("Tracked upstream files were modified")
    files = git(upstream, "ls-files", "-z").split("\0")
    files = [p for p in files if p]
    hashes = {p: digest(upstream / p) for p in files}
    # Editable setuptools_scm creates _version.py. No other untracked importable
    # source is allowed to shadow the pinned implementation.
    extra = sorted(
        str(p.relative_to(upstream))
        for p in (upstream / "crem").rglob("*.py")
        if str(p.relative_to(upstream)) not in hashes
    )
    if set(extra) - {"crem/_version.py"}:
        raise ValueError(f"Unexpected untracked CReM Python source: {extra}")
    return {
        "root": str(upstream),
        "url": UPSTREAM_URL,
        "commit": UPSTREAM,
        "tracked_files_sha256": hashes,
        "generated_files_sha256": {p: digest(upstream / p) for p in extra},
        "code_license": "BSD-3-Clause",
        "license_sha256": hashes["LICENSE.txt"],
    }


def runtime_snapshot(upstream):
    require_environment()
    if rdBase.rdkitVersion != "2025.09.6":
        raise ValueError("This comparison requires common identity/scoring RDKit 2025.09.6")
    upstream = Path(upstream).resolve()
    modules = {}
    for name, expected in (
        ("crem", upstream),
        ("crem.crem", upstream),
        ("crem.db", upstream),
        ("crem.scripts.cremdb_create", upstream),
        ("chembounce_synthesis.identity", ROOT),
        ("fingerprint_db", ROOT),
    ):
        module = importlib.import_module(name)
        path = Path(module.__file__).resolve()
        if not path.is_relative_to(expected):
            raise ValueError(f"Wrong import origin for {name}: {path}")
        modules[name] = artifact(path)
    for name in ("rdkit", "numpy", "pandas", "psutil", "joblib", "tqdm"):
        path = Path(importlib.import_module(name).__file__).resolve()
        if not path.is_relative_to(Path(sys.prefix).resolve()):
            raise ValueError(f"Dependency outside isolated environment: {name}")
        modules[name] = artifact(path)
    frozen_requirements = subprocess.check_output(
        [sys.executable, "-m", "pip", "freeze", "--all"], text=True
    ).splitlines()
    return {
        "python": sys.version,
        "executable": sys.executable,
        "prefix": sys.prefix,
        "platform": platform.platform(),
        "rdkit": rdBase.rdkitVersion,
        "environment": {k: os.environ[k] for k in ENV},
        "pip_freeze_all": frozen_requirements,
        "origins": modules,
        "distributions": dict(
            sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions())
        ),
        "distribution_scope": "Discovered metadata includes setuptools-vendored packages; pip freeze separately records installed requirements",
    }


def implementation():
    return {
        "root": str(ROOT),
        "commit": git(ROOT, "rev-parse", "HEAD"),
        "files_sha256": {p: digest(ROOT / p) for p in PROJECT_FILES},
    }


def verify_implementation(record):
    # A later docs-only commit does not change the algorithm; compare full file set.
    if record["files_sha256"] != implementation()["files_sha256"]:
        raise ValueError("Baseline/project implementation changed")


def normalized(smiles, *, for_training=False):
    mol = molecule(smiles)
    if for_training and any(a.GetIsotope() for a in mol.GetAtoms()):
        raise ValueError("isotopic_training_structure_not_retained_by_official_DB_builder")
    exact, parsed, rounds = stable_canonical_molecule(mol)
    potential = [s for s in Chem.FindPotentialStereo(parsed) if str(s.specified) == "Unspecified"]
    identity = chemical_identity(exact)
    cuts = len(parsed.GetSubstructMatches(Chem.MolFromSmarts("[!#1]!@!=!#[!#1]")))
    return (
        parsed,
        identity,
        {
            "canonical_roundtrips": rounds,
            "unspecified_stereo": len(potential),
            "acyclic_heavy_cut_bonds": cuts,
            "ordinary_fragmentation_maxCutBonds30_exceeded": cuts > 30,
        },
    )


def train_inventory(dataset_path, train_tsv, provenance_path):
    """Parent-side preparation only. Labels/private references never enter sampler assets."""
    from benchmarks.research.protocol import load_dataset

    dataset, audit = load_dataset(dataset_path)
    provenance = read_json(provenance_path)
    if (
        provenance.get("dataset_sha256") != digest(dataset_path)
        or provenance.get("library_sha256") != digest(train_tsv)
        or provenance.get("partition") != "train only"
        or provenance.get("labels_exposed") is not False
    ):
        raise ValueError("Train TSV provenance does not bind this dataset/library")
    allowed, forbidden, public_queries = {}, set(), {}
    for row in dataset["compounds"]:
        ident = dataset.identities[row["id"]]
        if row["split"] == "train":
            allowed.setdefault(ident.exact_smiles, []).append(row["id"])
        else:
            forbidden.add(ident.family_smiles)
    for query in dataset["queries"]:
        row = dataset.compounds_by_id[query["compound_id"]]
        public_queries[query["id"]] = {
            "compound_id": row["id"],
            "partition": row["split"],
            "family": dataset.identities[row["id"]].family_smiles,
        }
    with Path(train_tsv).open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if set(reader.fieldnames or []) != {"id", "smiles", "source", "snapshot"}:
            raise ValueError("Only unlabelled training-library TSV columns are allowed")
        rows = list(reader)
    if not rows or len(rows) != provenance.get("records"):
        raise ValueError("Train TSV row count is inconsistent")
    # Reproduce the *complete* existing training_library.export_library policy,
    # including first representative per family, ordering and IDs. Membership
    # alone would permit a silently reduced library with a rewritten row count.
    original_private = {
        dataset.identities[r["id"]].family_smiles
        for r in dataset["compounds"]
        if r["split"] in {"test", "query"}
    }
    expected_families = {}
    for row in dataset["compounds"]:
        identity = dataset.identities[row["id"]]
        if row["split"] == "train" and identity.family_smiles not in original_private:
            expected_families.setdefault(identity.family_smiles, identity.exact_smiles)
    expected_rows = [
        {
            "id": f"train-{index:06d}",
            "smiles": smiles,
            "source": "Benchmark development compound; not supplier availability",
            "snapshot": audit["dataset_sha256"],
        }
        for index, smiles in enumerate(sorted(expected_families.values()))
    ]
    records, exclusions, seen_ids, seen_exact, seen_family = [], [], set(), set(), set()
    for index, row in enumerate(rows):
        if not row["id"] or row["id"] in seen_ids or row["snapshot"] != audit["dataset_sha256"]:
            raise ValueError("Duplicate/missing training ID or contradictory dataset hash")
        seen_ids.add(row["id"])
        ident = chemical_identity(row["smiles"])
        if ident.exact_smiles not in allowed:
            raise ValueError("TSV molecule is not an exact member of the train partition")
        if ident.family_smiles in forbidden:
            raise ValueError("Nontrain partition family leakage in training library")
        if ident.exact_smiles in seen_exact or ident.family_smiles in seen_family:
            raise ValueError("Training TSV must contain one exact representative per family")
        seen_exact.add(ident.exact_smiles)
        seen_family.add(ident.family_smiles)
        try:
            _, accepted, state = normalized(row["smiles"], for_training=True)
            records.append(
                {
                    "input_index": index,
                    "id": row["id"],
                    "input_smiles": row["smiles"],
                    **asdict(accepted),
                    **state,
                    "dataset_train_ids": allowed[ident.exact_smiles],
                }
            )
        except ValueError as exc:
            exclusions.append(
                {
                    "input_index": index,
                    "id": row["id"],
                    "reason": str(exc),
                    "input_smiles": row["smiles"],
                    **asdict(ident),
                }
            )
    retained = {r["family_smiles"] for r in records}
    overlaps = [
        {"query_id": qid, **row, "family_sha256": json_hash(row["family"])}
        for qid, row in public_queries.items()
        if row["family"] in retained
    ]
    if any(r["partition"] != "train" for r in overlaps):
        raise ValueError("Nontrain public-query family unexpectedly retained")
    if not records:
        raise ValueError("No supported training compounds remain")
    if rows != expected_rows:
        raise ValueError("TSV differs from the exact complete deterministic training-library export")
    summary = {
        "input_records": len(rows),
        "input_exact_inventory_sha256": json_hash(sorted(seen_exact)),
        "input_family_inventory_sha256": json_hash(sorted(seen_family)),
        "retained_records": len(records),
        "unique_exact": len({r["exact_smiles"] for r in records}),
        "unique_families": len(retained),
        "exclusions": exclusions,
        "nontrain_partition_family_overlap": 0,
        "public_train_query_overlap": overlaps,
        "role_policy": "Public train-query families retained; all nontrain partition families excluded",
        "identity_policy": IDENTITY_POLICY,
        "rdkit": rdBase.rdkitVersion,
        "exact_inventory_sha256": json_hash(sorted(r["exact_smiles"] for r in records)),
        "family_inventory_sha256": json_hash(sorted(retained)),
        "scope": "Membership in this retained train TSV/panel only; no general clean-holdout claim",
    }
    return records, summary


def prepare(dataset, train_tsv, provenance, upstream, output, *, expected_hashes):
    output = Path(output).resolve()
    if output.exists():
        raise ValueError("Output already exists; preparation does not overwrite")
    original_inputs = {"dataset": dataset, "train_tsv": train_tsv, "provenance": provenance}
    if set(expected_hashes) != set(original_inputs) or any(
        digest(path) != expected_hashes[name] for name, path in original_inputs.items()
    ):
        raise ValueError("Parent-declared dataset/train TSV/provenance SHA-256 mismatch")
    source = source_snapshot(upstream)
    runtime = runtime_snapshot(upstream)
    records, inventory = train_inventory(dataset, train_tsv, provenance)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "inventory.json", {"records": records, "audit": inventory})
    # Sort exact structures before builder insertion; upstream set/hash iteration
    # is also fixed by process-start PYTHONHASHSEED=0.
    with (output / "train.smi").open("x") as f:
        for row in sorted(records, key=lambda r: r["exact_smiles"]):
            f.write(f"{row['exact_smiles']} {row['id']}\n")
    plan = {
        "schema_version": 1,
        "policy": POLICY,
        "stage": "prepared_not_built",
        "parent_declared_input_sha256": expected_hashes,
        "inputs": {
            "dataset": artifact(dataset),
            "train_tsv": artifact(train_tsv),
            "provenance": artifact(provenance),
            "inventory": artifact(output / "inventory.json"),
            "builder_input": artifact(output / "train.smi"),
        },
        "upstream": source,
        "runtime": runtime,
        "implementation": implementation(),
        "build_settings": BUILD,
        "mutation_settings": MUTATION,
        "inventory_summary": inventory,
        "budget_proposal": {
            "seeds": [42, 43, 44],
            "wall_seconds": 360,
            "max_memory_mb": 24000,
            "max_candidates": 3,
            "max_evaluations": 256,
            "mutation_calls": 1,
        },
        "notes": {
            "ranking": "Unfolded count Morgan radius 2, no chirality; descending Tanimoto then exact SMILES",
            "query_and_output_stereo": "Preserve supplied state; no enumeration or repair; unspecified stereo is recorded, not a 2D exclusion",
            "training_stereo": "Keep source annotations; unspecified train stereo may survive into fragments",
            "isotopes": "Training isotopic records excluded because official builder strips isotope labels; query and product labels supported via official backup/restore",
            "geometry": "Official rigid-ring ETKDGv3 seed 42/maxIterations 4 and torsion scan; exceptions fail open",
            "route": "No synthetic route supplied; CReM edit transformations are not synthesis evidence",
            "work_unit": "API-yielded products, not internal fragment matches; upstream suppresses some failures/duplicates",
        },
    }
    write_json(output / "plan.json", plan)
    return plan


def verify_plan(path):
    plan = read_json(path)
    if (
        plan.get("policy") != POLICY
        or plan.get("stage") != "prepared_not_built"
        or plan.get("build_settings") != BUILD
        or plan.get("mutation_settings") != MUTATION
    ):
        raise ValueError("Unsupported or changed preparation policy")
    for spec in plan["inputs"].values():
        verify_artifact(spec)
    if plan.get("parent_declared_input_sha256") != {
        k: plan["inputs"][k]["sha256"] for k in ("dataset", "train_tsv", "provenance")
    }:
        raise ValueError("Prepared inputs differ from parent-declared original hashes")
    verify_implementation(plan["implementation"])
    if source_snapshot(plan["upstream"]["root"]) != plan["upstream"]:
        raise ValueError("Upstream source changed")
    if runtime_snapshot(plan["upstream"]["root"]) != plan["runtime"]:
        raise ValueError("Runtime changed since preparation")
    return plan


def database_summary(path):
    with sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro&immutable=1", uri=True) as con:
        if con.execute("PRAGMA integrity_check").fetchall() != [("ok",)]:
            raise ValueError("SQLite integrity failure")
        if con.execute("PRAGMA user_version").fetchone()[0] != 2:
            raise ValueError("CReM v2 database is required")
        columns = [r[1] for r in con.execute("PRAGMA table_info(radius3)")]
        if set(columns) != {"env_id", "core_smi_id", "core_num_atoms", "dist2", "train"}:
            raise ValueError("Unexpected fragment sets or radius schema")
        tables = {
            name: con.execute(f'SELECT COUNT(*) FROM "{name}"').fetchone()[0]
            for name in ("envs", "frags", "frags_h", "radius3")
        }
        ring = con.execute(
            "SELECT COUNT(*) FROM radius3 r JOIN frags f USING(core_smi_id) WHERE f.core_smi LIKE '%[1*%'"
        ).fetchone()[0]
        if not tables["radius3"] or not ring:
            raise ValueError("DB has no ring-fragment coverage for partial_exo")
        if con.execute("SELECT COUNT(*) FROM radius3 WHERE train<1 OR core_num_atoms>15").fetchone()[0]:
            raise ValueError("DB violates frequency/core-size policy")
    return {
        "schema_version": 2,
        "table_rows": tables,
        "ring_mapping_rows": ring,
        "bytes": Path(path).stat().st_size,
    }


def build_worker(request_path, response_path):
    request = read_json(request_path)
    if digest(request["assets"]["plan"]) != request["settings"]["plan_sha256"]:
        raise ValueError("Prepared plan changed before worker launch")
    plan = verify_plan(request["assets"]["plan"])
    db = Path(request["assets"]["new_database"])
    if db.exists() or Path(response_path).exists():
        raise ValueError("Builder refuses additive/overwrite operation")
    from crem.db import create_db

    started, cpu = time.monotonic(), time.process_time()
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    create_db(plan["inputs"]["builder_input"]["path"], str(db), **BUILD)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    summary = database_summary(db)
    write_json(
        response_path,
        {
            "candidates": [],
            "evaluations": plan["inventory_summary"]["retained_records"],
            "evaluation_unit": "offline training molecules passed to official DB builder",
            "metadata": {
                "database": summary,
                "worker_wall_seconds": time.monotonic() - started,
                "worker_cpu_seconds": time.process_time() - cpu,
                "completed_children_cpu_seconds": children.ru_utime
                + children.ru_stime
                - before_children.ru_utime
                - before_children.ru_stime,
            },
        },
    )


def build(plan_path, output, wall_seconds=3600, max_memory_mb=24000):
    from benchmarks.research.runner import bounded_method

    plan = verify_plan(plan_path)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    db = output / "fragments.db"
    request = {
        "assets": {"plan": str(Path(plan_path).resolve()), "new_database": str(db)},
        "settings": {"plan_sha256": digest(plan_path)},
        "budget": {
            "wall_seconds": wall_seconds,
            "max_memory_mb": max_memory_mb,
            "max_candidates": 1,
            "max_evaluations": plan["inventory_summary"]["retained_records"],
        },
    }
    method = {
        "kind": "external",
        "command": [
            sys.executable,
            str(Path(__file__).resolve()),
            "build-worker",
            "--request",
            "{request}",
            "--response",
            "{response}",
        ],
    }
    report = bounded_method(method, request, output / "worker")
    report["plan"] = artifact(plan_path)
    report["worker_completed"] = report["status"] == "no_candidates"
    report["successful_database_build"] = False
    write_json(output / "build-report.json", report)
    if not report["worker_completed"]:
        raise RuntimeError("DB build did not complete; partial DB is diagnostic only")
    stdout = (output / "worker" / "stdout.log").read_text()
    counters = {}
    for name in (
        "Molecules processed",
        "Fragments generated",
        "Env/core pairs",
        "MMPA failures",
        "Fragment validation issues",
        "Chunks processed",
        "Chunks skipped",
    ):
        matches = re.findall(r"^" + re.escape(name) + r": (\d+)$", stdout, re.M)
        if len(matches) != 1:
            raise ValueError(f"Official builder counter unavailable: {name}")
        counters[name] = int(matches[0])
    if (
        counters["Molecules processed"] != plan["inventory_summary"]["retained_records"]
        or counters["Chunks skipped"]
    ):
        raise ValueError("Official builder did not process every prepared training molecule")
    database = database_summary(db)
    report["successful_database_build"] = True
    report["official_builder_counters"] = counters
    write_json(output / "build-report.json", report)
    # Public sampler manifest: only hashes/counts of private preparation inputs.
    # No labelled panel path or reference structure inventory is exposed.
    public = {
        "schema_version": 1,
        "policy": POLICY,
        "complete": True,
        "db_sha256": digest(db),
        "database": database,
        "upstream": plan["upstream"],
        "runtime": plan["runtime"],
        "implementation": plan["implementation"],
        "build_settings": BUILD,
        "mutation_settings": MUTATION,
        "official_builder_counters": counters,
        "resources": report["resources"],
        "build_report_sha256": digest(output / "build-report.json"),
        "input_hashes": {k: v["sha256"] for k, v in plan["inputs"].items()},
        "inventory_counts": {
            k: plan["inventory_summary"][k]
            for k in (
                "retained_records",
                "unique_exact",
                "unique_families",
                "nontrain_partition_family_overlap",
            )
        },
        "public_train_query_overlap_count": len(plan["inventory_summary"]["public_train_query_overlap"]),
        "exact_inventory_sha256": plan["inventory_summary"]["exact_inventory_sha256"],
        "family_inventory_sha256": plan["inventory_summary"]["family_inventory_sha256"],
        "build_artifacts_sha256": {
            str(p.relative_to(output)): digest(p)
            for p in (
                output / "worker" / "stdout.log",
                output / "worker" / "stderr.log",
                output / "fragments.errors",
            )
            if p.is_file()
        },
        "fragment_coverage": "Every input line processed; upstream MMPA/fragment exclusions counted, no exhaustive chemical-space claim",
        "upstream_zero_fragment_molecule_ids": None,
        "upstream_zero_fragment_scope": "Public create_db reports aggregate counts, not a per-molecule zero-fragment inventory",
        "labels_exposed_to_sampler": False,
        "scope": "Train membership verified by preparation; no universal leakage or synthesis claim",
    }
    write_json(output / "database-manifest.json", public)
    return public


def validate_request(request):
    allowed = {
        "schema_version",
        "kind",
        "method_id",
        "query",
        "seed",
        "budget",
        "assets",
        "settings",
        "output_dir",
        "deadline_monotonic",
        "finalization_margin_seconds",
    }
    if set(request) - allowed or request.get("schema_version") != 1 or request.get("kind") != "external":
        raise ValueError("Only the external public-request schema is accepted")
    if set(request.get("query", {})) - {"id", "smiles", "target_id", "reference_sdf", "reference_sdf_sha256"}:
        raise ValueError("Query contains private/unknown fields")
    if set(request.get("assets", {})) != {"crem_db", "crem_manifest", *SOURCE_ASSETS}:
        raise ValueError("Sampler assets must be exactly DB, manifest and the project source allowlist")
    settings = request.get("settings", {})
    if set(settings) != {"policy", "manifest_sha256"} or settings["policy"] != POLICY:
        raise ValueError("Unrecognized or mutable sampler settings")
    if not re.fullmatch(r"[0-9a-f]{64}", settings["manifest_sha256"]):
        raise ValueError("Public manifest SHA-256 is required")
    if type(request.get("seed")) is not int or request["seed"] not in (42, 43, 44):
        raise ValueError("Declared seeds are 42, 43, 44")
    budget = request.get("budget", {})
    for key, cap in (("max_candidates", 3), ("max_evaluations", 256)):
        if type(budget.get(key)) is not int or not 1 <= budget[key] <= cap:
            raise ValueError(f"Invalid {key}")
    for key, cap in (("wall_seconds", 360), ("max_memory_mb", 24000)):
        value = budget.get(key)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or not 0 < value <= cap
        ):
            raise ValueError(f"Invalid {key}")
    if not isinstance(request.get("query", {}).get("smiles"), str):
        raise ValueError("Missing query SMILES")


def verify_database(request):
    assets = request["assets"]
    if digest(assets["crem_manifest"]) != request["settings"]["manifest_sha256"]:
        raise ValueError("Public DB manifest changed")
    manifest = read_json(assets["crem_manifest"])
    if (
        manifest.get("policy") != POLICY
        or manifest.get("complete") is not True
        or manifest.get("build_settings") != BUILD
        or manifest.get("mutation_settings") != MUTATION
        or manifest.get("labels_exposed_to_sampler") is not False
        or manifest["inventory_counts"].get("nontrain_partition_family_overlap") != 0
    ):
        raise ValueError("Incomplete or incompatible train-only database")
    if digest(assets["crem_db"]) != manifest["db_sha256"]:
        raise ValueError("CReM database changed")
    for name, relative in SOURCE_ASSETS.items():
        if Path(assets[name]).resolve() != (ROOT / relative).resolve():
            raise ValueError("Pinned project source asset points outside the executing checkout")
        if digest(assets[name]) != manifest["implementation"]["files_sha256"][relative]:
            raise ValueError("Pinned project source asset changed")
    verify_implementation(manifest["implementation"])
    if source_snapshot(manifest["upstream"]["root"]) != manifest["upstream"]:
        raise ValueError("Upstream CReM source changed")
    if runtime_snapshot(manifest["upstream"]["root"]) != manifest["runtime"]:
        raise ValueError("Sampler/build runtime mismatch")
    return manifest


def method_fragment(database, manifest_path, manifest_sha256, output):
    """Emit independent runner asset pins; no query, labels or private inventory is read."""
    output = Path(output).resolve()
    if output.exists():
        raise ValueError("Method fragment output already exists")
    assets = {
        "crem_db": str(Path(database).resolve()),
        "crem_manifest": str(Path(manifest_path).resolve()),
        **{key: str(ROOT / rel) for key, rel in SOURCE_ASSETS.items()},
    }
    settings = {"policy": POLICY, "manifest_sha256": manifest_sha256}
    manifest = verify_database({"assets": assets, "settings": settings})
    method = {
        "id": "crem_train_partial_exo",
        "kind": "external",
        "command": [
            sys.executable,
            str(Path(__file__).resolve()),
            "sample",
            "--request",
            "{request}",
            "--response",
            "{response}",
        ],
        "assets": {key: {"artifact": value, "sha256": digest(value)} for key, value in assets.items()},
        "settings": settings,
        "budget_overrides": {"max_evaluations": 256},
        "source_binding": {
            "scope": "Existing runner independently hashes this exact project-source asset allowlist before launch",
            "relative_project_files": list(PROJECT_FILES),
            "upstream_commit": manifest["upstream"]["commit"],
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, method)
    return method


def select_products(query, products, max_candidates=3):
    """Deterministic postfilter/ranking. No target labels, reference poses or DB reads."""
    query_mol, query_id, _ = normalized(query)
    if not query_id.has_scaffold:
        raise ValueError("query_has_no_ring_scaffold")
    qfp = COUNT_FP.GetSparseCountFingerprint(query_mol)
    rows, accepted = [], []
    for index, value in enumerate(products):
        row = {"index": index, "raw_smiles": value}
        try:
            mol, identity, state = normalized(value)
            row.update(asdict(identity), **state)
            if not identity.has_scaffold:
                raise ValueError("no_ring_scaffold")
            if identity.family_smiles == query_id.family_smiles:
                raise ValueError("same_query_family")
            if identity.scaffold_family_smiles == query_id.scaffold_family_smiles:
                raise ValueError("unchanged_murcko_scaffold")
            row["count_tanimoto"] = DataStructs.TanimotoSimilarity(
                qfp, COUNT_FP.GetSparseCountFingerprint(mol)
            )
            row["status"] = "eligible"
            accepted.append(row)
        except (ValueError, RuntimeError, TypeError) as exc:
            row.update(status="rejected", reason=str(exc))
        rows.append(row)
    accepted.sort(key=lambda r: (-r["count_tanimoto"], r["exact_smiles"], r["index"]))
    families, selected = set(), []
    for row in accepted:
        if row["family_smiles"] in families:
            row.update(status="rejected", reason="duplicate_family")
            continue
        families.add(row["family_smiles"])
        if len(selected) >= max_candidates:
            row["status"] = "eligible_outside_output_cap"
            continue
        row["status"] = "selected"
        selected.append(
            {
                "smiles": row["exact_smiles"],
                "score_samples": [row["count_tanimoto"]],
                "ranking_metric": "Morgan radius2 sparse counts, includeChirality=False",
                "route_evidence": {
                    "level": "no_route_supplied",
                    "reason": "CReM fragment edits are not synthetic routes",
                },
            }
        )
    return selected, rows


def sample(request, response_path):
    """One official mutation call; no test panel or training inventory is opened."""
    response_path = Path(response_path)
    raw_path = response_path.with_name(response_path.stem + ".raw.jsonl")
    progress_path = response_path.with_name(response_path.stem + ".progress.json")
    if any(p.exists() for p in (response_path, raw_path, progress_path)):
        raise ValueError("Sampler output exists; no overwrite or resume")
    validate_request(request)
    started = time.monotonic()
    # The external runner starts this clock before imports. Standalone callers
    # without that clock cannot claim the matched cold-start budget.
    deadline = request.get("deadline_monotonic")
    if isinstance(deadline, bool) or not isinstance(deadline, (int, float)) or not math.isfinite(deadline):
        raise ValueError("A bounded_method monotonic deadline is required")
    if deadline > started + request["budget"]["wall_seconds"] + 0.1:
        raise ValueError("Deadline exceeds the declared wall budget")
    manifest = verify_database(request)
    metadata = {
        "policy": POLICY,
        "request_sha256": json_hash(request),
        "db_sha256": manifest["db_sha256"],
        "upstream_commit": UPSTREAM,
        "mutation_settings": MUTATION,
        "mutation_seed": request["seed"],
        "internal_geometry_seed": 42,
        "internal_geometry_failure_policy": "upstream fail-open",
        "mutation_calls": 0,
        "upstream_internal_attempts": None,
        "upstream_internal_invalid_duplicate_geometry_rejections": None,
        "upstream_count_limitation": "API suppresses these before yielding; unknown is not zero",
        "route_status": "no_route_supplied",
        "raw_products": 0,
        "postfilter_counts": {},
        "stop_reason": None,
    }

    def publish_progress():
        write_json(progress_path, metadata)

    publish_progress()
    raw, success = [], False
    try:
        mol, ident, normalization = normalized(request["query"]["smiles"])
        if not ident.has_scaffold:
            raise ValueError("query_has_no_ring_scaffold")
        metadata["query_normalization"] = {**asdict(ident), **normalization}
    except (ValueError, RuntimeError) as exc:
        metadata.update(stop_reason="unsupported_input", error=str(exc))
        write_json(
            response_path,
            {
                "candidates": [],
                "evaluations": 0,
                "evaluation_unit": "CReM API-yielded products",
                "metadata": metadata,
            },
        )
        publish_progress()
        return
    from crem.crem import mutate_mol

    with raw_path.open("x") as raw_file:
        # Include normalization/ranking inside the same deadline. No cooperative
        # partial acceptance: if work stops for wall time, all products stay raw.
        try:
            if time.monotonic() >= deadline - 2:
                raise TimeoutError("deadline reached during cold-start validation")
            metadata["mutation_calls"] = 1
            settings = MUTATION | {"max_replacements": request["budget"]["max_evaluations"]}
            generator = mutate_mol(mol, request["assets"]["crem_db"], seed=request["seed"], **settings)
            try:
                for product in generator:
                    raw.append(product)
                    raw_file.write(json.dumps({"index": len(raw) - 1, "raw_smiles": product}) + "\n")
                    raw_file.flush()
                    metadata["raw_products"] = len(raw)
                    publish_progress()
                    if time.monotonic() >= deadline - 2:
                        raise TimeoutError("deadline reached during official mutation")
                    if len(raw) >= request["budget"]["max_evaluations"]:
                        metadata["stop_reason"] = "api_product_cap"
                        break
                else:
                    metadata["stop_reason"] = "official_generator_exhausted"
            finally:
                generator.close()
            candidates, decisions = select_products(
                request["query"]["smiles"], raw, request["budget"]["max_candidates"]
            )
            if time.monotonic() >= deadline - 0.2:
                raise TimeoutError("deadline reached during postfilter/ranking")
            metadata["postfilter_counts"] = dict(Counter(r.get("reason", r["status"]) for r in decisions))
            metadata["decisions"] = decisions
            success = True
        except Exception as exc:
            metadata.update(
                stop_reason="internal_timeout" if isinstance(exc, TimeoutError) else "generation_error",
                error=f"{type(exc).__name__}: {exc}",
            )
            candidates = []
        metadata["sampler_wall_seconds_excluding_import"] = time.monotonic() - started
        metadata["completed_within_budget"] = success
        publish_progress()
        write_json(
            response_path,
            {
                "candidates": candidates,
                "evaluations": len(raw),
                "evaluation_unit": "CReM API-yielded products",
                "metadata": metadata,
            },
        )
    if not success:
        # The outer runner retains the files but marks this attempt unsuccessful;
        # even completed partial molecules cannot become accepted candidates.
        raise RuntimeError(metadata["stop_reason"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="action", required=True)
    p = subs.add_parser("prepare")
    for name in ("dataset", "train-tsv", "provenance", "upstream", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    for name in ("dataset", "train-tsv", "provenance"):
        p.add_argument("--" + name + "-sha256", required=True)
    p = subs.add_parser("build")
    p.add_argument("--plan", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--wall-seconds", type=float, default=3600)
    p.add_argument("--max-memory-mb", type=float, default=24000)
    p = subs.add_parser("method")
    for name in ("database", "manifest", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--manifest-sha256", required=True)
    for name in ("sample", "build-worker"):
        p = subs.add_parser(name)
        p.add_argument("--request", type=Path, required=True)
        p.add_argument("--response", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.action == "prepare":
        prepare(
            args.dataset,
            args.train_tsv,
            args.provenance,
            args.upstream,
            args.output,
            expected_hashes={
                "dataset": args.dataset_sha256,
                "train_tsv": args.train_tsv_sha256,
                "provenance": args.provenance_sha256,
            },
        )
    elif args.action == "build":
        if not 0 < args.wall_seconds <= 86400 or not 0 < args.max_memory_mb <= 24000:
            parser.error("Build limits must be positive; maximum 86400 seconds / 24000 MB")
        build(args.plan, args.output, args.wall_seconds, args.max_memory_mb)
    elif args.action == "build-worker":
        build_worker(args.request, args.response)
    elif args.action == "method":
        method_fragment(args.database, args.manifest, args.manifest_sha256, args.output)
    else:
        sample(read_json(args.request), args.response)


if __name__ == "__main__":
    main()
