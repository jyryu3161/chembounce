#!/usr/bin/env python3
"""Frozen-generator, index-only diagnostic. Never updates a campaign or runs models.

Plan creation snapshots references; execution requires a separate explicit command.
This deliberately has no resume/overwrite mode. Interrupted diagnostics remain
partial and require review, rather than silently rerunning selected conditions.
"""

import argparse
import collections
import concurrent.futures
import copy
import csv
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import threading
import time
import traceback


SCHEMA = "frozen-generator-index-sensitivity-v1"
DEPENDENCIES = ("rdkit", "numpy", "pandas", "scipy", "scaffoldgraph", "oddt", "psutil")
THREADS = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "USE_TF": "0"}
OPTIONS = {
    "calculate_shape": False,
    "candidate_thresholds": {},
    "exploration_fraction": 0.25,
    "frag_max_n": 10,
    "overall_max_n": 3,
    "protocol": "research-v1",
    "scaffold_top_n": 8,
    "search_candidate_budget": 512,
    "tqdm_quiet": True,
}
DEFAULTS = {
    "tanimoto_threshold": 0.5,
    "_search_scf_thr_": 0.3,
    "cand_max_n__rplc": None,
    "_merge_structure_top_n_": 100,
    "murcko_frag_itr_rnd": 1,
    "core_smiles": None,
    "fragments": None,
    "replace_scaffold_files": None,
    "low_mem": False,
    "use_fingerprint_db": True,
    "trust_legacy_db": False,
}
TABLE_COLUMNS = [
    "Fragment_no",
    "Original scaffold",
    "Replaced scaffold",
    "Final structure",
    "Standardized final structure",
    "Tanimoto Similarity",
    "Electron shape Similarity",
    "CATS2D dist",
    "QED",
    "SAscore",
    "logP",
    "MW",
    "H_Donors",
    "H_Acceptors",
]
PATTERNS = (
    "run_manifest.json",
    "run_summary.json",
    "overall_result.txt",
    "fragment_info.tsv",
    "fragment_result_*.txt",
    "replace_scaffold_list.fragment_*.tsv",
    "search_decisions.json",
)
LIMITATIONS = [
    "Generator-only diagnostic; no corrected-method activity/synthesis/performance evidence.",
    "Identical selected designs do not establish full-pipeline equality under its 360-second budget.",
    "Canonical spellings can change hash reservoirs and tie order, beyond fingerprint-bit changes.",
    "Generator-complete references never change a failed or incomplete original campaign outcome.",
    "All 126 conditions and 252 arm calls remain in denominators, including failures and unattempted calls.",
    "Wall/RSS observations on a shared host are descriptive; no speedup inference or global optimality claim.",
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def object_sha(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def atomic_json(path, value):
    path = Path(path)
    temp = path.with_name(path.name + ".tmp")
    created = False
    try:
        with temp.open("x") as stream:
            created = True
            json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp, path)
    finally:
        # A signal during finalization must not leave a stale temp blocking the
        # interruption checkpoint; the previously published file remains intact.
        if created:
            temp.unlink(missing_ok=True)


def file_pin(path):
    p = Path(path).resolve()
    return {"path": str(p), "sha256": sha(p), "bytes": p.stat().st_size}


def check_pin(pin):
    require(sha(pin["path"]) == pin["sha256"], "Changed pinned file: " + pin["path"])


def directory_pin(path):
    p = Path(path).resolve()
    require(p.is_dir(), f"Missing directory: {p}")
    files = {str(f.relative_to(p)): sha(f) for f in sorted(p.rglob("*")) if f.is_file()}
    require(files, f"Empty directory: {p}")
    return {"path": str(p), "files_sha256": files}


def check_directory(pin):
    require(directory_pin(pin["path"]) == pin, "Changed pinned directory: " + pin["path"])


def version_tuple(value):
    return tuple(int(v) for v in value.split("."))


def runtime():
    from rdkit import rdBase

    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "process_niceness": os.getpriority(os.PRIO_PROCESS, 0),
        "rdkit_runtime": rdBase.rdkitVersion,
        "dependencies": {n: importlib.metadata.version(n) for n in DEPENDENCIES},
    }


def source_check(root, expected):
    root = Path(root).resolve()
    for relative, digest in expected.items():
        path = (root / relative).resolve()
        require(path.is_relative_to(root), "Source path escapes root")
        require(sha(path) == digest, f"Frozen source mismatch: {relative}")
    # Generator provenance enumerates these files. Extra modules are not silently ignored.
    actual_generator = {str(p.relative_to(root)) for p in root.glob("*.py")} | {
        str(p.relative_to(root)) for p in (root / "chembounce_synthesis").glob("*.py")
    }
    expected_generator = {p for p in expected if "/" not in p or p.startswith("chembounce_synthesis/")}
    require(actual_generator == expected_generator, "Frozen generator source inventory changed")


def install_source(root, expected):
    source_check(root, expected)
    root = str(Path(root).resolve())
    if root in sys.path:
        sys.path.remove(root)
    sys.path.insert(0, root)


def canonical_selection(raw, max_designs=2):
    # Exact frozen workflow.run semantics, without importing/running workflow.
    from chembounce_synthesis.common import canonical

    all_designs = list(dict.fromkeys(canonical(s) for s in raw))
    return {"canonical_deduplicated": all_designs, "downstream_selected": all_designs[:max_designs]}


def validate_index_report(index, scaffold_sha, runtime_record, report):
    metadata = read_json(Path(index["path"]) / "metadata.json")
    require(metadata.get("schema_version") == 2, "Expected directory index schema 2")
    require(metadata.get("source_sha256") == scaffold_sha, "New index source differs")
    require(metadata.get("canonical_sorted_unique") is True, "New index must declare sorted unique strings")
    require(
        metadata.get("canonicalization_policy") == "rdkit-roundtrip-fixed-point-v1"
        and metadata.get("canonicalization_max_roundtrips") == 8,
        "Expected finalized bounded fixed-point builder policy",
    )
    require(
        version_tuple(metadata["rdkit_version"]) == version_tuple(runtime_record["rdkit_runtime"]),
        "New index runtime differs",
    )
    params = {"bitorder": "big", "fp_size": 2048, "include_chirality": False, "radius": 2}
    require(all(metadata.get(k) == v for k, v in params.items()), "Fingerprint parameters changed")
    require(report.get("complete") is True, "Independent exhaustive index audit incomplete")
    require(report.get("input_asset_sha256") == index["files_sha256"], "Audit index hashes differ")
    require(report.get("fingerprint_parameters") == params, "Audit fingerprint parameters differ")
    rows = metadata.get("rows")
    require(isinstance(rows, int) and rows > 0, "Invalid index row count")
    require(
        report.get("rows")
        == report.get("fingerprints_checked")
        == report.get("runtime_canonical_unique_strings")
        == rows,
        "Audit did not cover every unique row",
    )
    require(report.get("raw_strings_strictly_sorted_unique") is True, "Audit strings not sorted unique")
    for field in (
        "invalid_rows",
        "fingerprint_mismatches",
        "canonical_string_changed_rows",
        "canonical_duplicate_rows",
    ):
        require(report.get(field) == [], f"Independent index audit has unassessed or failing {field}")
    for field in ("database_rdkit_version", "runtime_rdkit_version"):
        require(
            version_tuple(report[field]) == version_tuple(runtime_record["rdkit_runtime"]),
            "Audit runtime differs",
        )
    return metadata


def verified_asset(spec, base, directory=False):
    path = (base / spec["directory" if directory else "artifact"]).resolve()
    pin = directory_pin(path) if directory else file_pin(path)
    require(
        pin["files_sha256" if directory else "sha256"] == spec["files_sha256" if directory else "sha256"],
        "Protocol asset hash mismatch",
    )
    return pin


def build_grid(campaigns, enforce_counts=True):
    """Extract public query inputs only; never normalize or pass held-out references."""
    jobs, inputs, asset_cache = {}, [], {}
    explicit_thread_parents, linked_thread_parents = set(), []
    for dataset_path, protocol_path, archive_path in campaigns:
        dataset_path, protocol_path, archive_path = map(
            lambda p: Path(p).resolve(), (dataset_path, protocol_path, archive_path)
        )
        dataset, protocol = read_json(dataset_path), read_json(protocol_path)
        dataset_pin, protocol_pin = file_pin(dataset_path), file_pin(protocol_path)
        require(protocol["seeds"] == [42, 43, 44], "Unexpected frozen seeds")
        budget = protocol["budget"]
        require(
            budget["max_candidates"] == 3
            and budget["wall_seconds"] == 360
            and budget["max_memory_mb"] == 24000,
            "Unexpected frozen outer budgets",
        )
        declared_threads = protocol["campaign"].get("thread_environment")
        if declared_threads is not None:
            require({k: str(v) for k, v in declared_threads.items()} == THREADS, "Unexpected frozen threads")
            explicit_thread_parents.add(protocol["campaign"].get("parent_protocol_sha256"))
        else:
            primary_hash = protocol["campaign"].get("primary_protocol_sha256")
            require(
                isinstance(primary_hash, str) and len(primary_hash) == 64,
                "Missing secondary primary-protocol binding",
            )
            linked_thread_parents.append(primary_hash)
        compounds = {r["id"]: r for r in dataset["compounds"]}
        require(len(compounds) == len(dataset["compounds"]), "Duplicate compound IDs")
        require(len({q["id"] for q in dataset["queries"]}) == len(dataset["queries"]), "Duplicate query IDs")
        if (archive_path / "manifest.json").exists():
            archive_manifest = read_json(archive_path / "manifest.json")
            require(archive_manifest["protocol_sha256"] == protocol_pin["sha256"], "Archive protocol differs")
            require(
                archive_manifest["dataset"]["dataset_sha256"] == dataset_pin["sha256"],
                "Archive dataset differs",
            )
        inputs.append(
            {
                "dataset": dataset_pin,
                "protocol": protocol_pin,
                "archive": str(archive_path),
                "query_count": len(dataset["queries"]),
                "declared_thread_environment": declared_threads,
                "thread_scope": "explicit primary declaration"
                if declared_threads is not None
                else "secondary protocol links primary; no explicit historical thread declaration; diagnostic uses matched primary settings",
            }
        )
        for method in protocol["methods"]:
            if method["kind"] != "pipeline":
                continue
            options = copy.deepcopy(method.get("settings", {}).get("generator_options", {}))
            policy = options.pop("search_policy", None)
            require(options == OPTIONS and policy in {"legacy", "diverse"}, "Unexpected generator options")
            options["search_policy"] = policy
            require(
                not method.get("settings", {}).get("disable_synthesis"),
                "Compound design-only baseline is not this frozen grid",
            )
            assets = {}
            for name in ("fingerprint_db", "scaffold_db", "synthesis_config"):
                key = object_sha([str(protocol_path.parent), method["assets"][name]])
                if key not in asset_cache:
                    asset_cache[key] = verified_asset(
                        method["assets"][name], protocol_path.parent, name == "fingerprint_db"
                    )
                assets[name] = asset_cache[key]
            config = read_json(assets["synthesis_config"]["path"])
            config.update(method.get("settings", {}).get("synthesis_overrides", {}))
            require(config.get("max_designs", 20) == 2, "Frozen downstream design cap changed")
            for query_index, query in enumerate(dataset["queries"]):
                record = compounds[query["compound_id"]]
                public = {"id": query["id"], "smiles": record["smiles"], "target_id": record["target_id"]}
                require(isinstance(public["smiles"], str) and public["smiles"], "Empty query")
                for seed in protocol["seeds"]:
                    key = {
                        "query_smiles": public["smiles"],
                        "seed": seed,
                        "options": options,
                        "scaffold_sha256": assets["scaffold_db"]["sha256"],
                        "old_index_files_sha256": assets["fingerprint_db"]["files_sha256"],
                        "max_designs": 2,
                    }
                    ident = object_sha(key)
                    if ident not in jobs:
                        jobs[ident] = {"id": ident, **key, "members": []}
                    jobs[ident]["members"].append(
                        {
                            "query": public,
                            "query_index": query_index,
                            "seed": seed,
                            "method_id": method["id"],
                            "settings": method.get("settings", {}),
                            "run_directory": str(
                                archive_path / method["id"] / f"query-{query_index:04d}-seed-{seed}"
                            ),
                            "protocol_sha256": protocol_pin["sha256"],
                            "dataset_sha256": dataset_pin["sha256"],
                            "synthesis_config": assets["synthesis_config"],
                            "scaffold": assets["scaffold_db"],
                            "old_index": assets["fingerprint_db"],
                        }
                    )
    jobs = list(jobs.values())
    require(
        all(h in explicit_thread_parents for h in linked_thread_parents),
        "Secondary primary-protocol lineage absent from supplied inputs",
    )
    counts = collections.Counter(j["options"]["search_policy"] for j in jobs)
    if enforce_counts:
        require(
            len(jobs) == 126 and counts == {"diverse": 90, "legacy": 36},
            "Expected exactly 126 frozen generator conditions (90 diverse, 36 legacy)",
        )
        require(sum(len(j["members"]) for j in jobs) == 270, "Expected 270 frozen pipeline grid members")
        require(len({j["query_smiles"] for j in jobs}) == 30, "Expected 30 distinct public query strings")
        require(len({i["dataset"]["sha256"] for i in inputs}) == 3, "Expected three panels")
    require(
        len({object_sha(j["members"][0]["old_index"]) for j in jobs}) == 1,
        "Old index differs between conditions",
    )
    require(
        len({object_sha(j["members"][0]["scaffold"]) for j in jobs}) == 1,
        "Scaffold source differs between conditions",
    )
    return jobs, inputs


def table_view(path, max_designs=2):
    # Cells are lossless strings: full table comparisons never round numeric values.
    rows = list(csv.reader(io.StringIO(Path(path).read_text()), delimiter="\t"))
    require(rows and rows[0] == TABLE_COLUMNS, "Unexpected generator table schema")
    require(all(len(r) == len(TABLE_COLUMNS) for r in rows[1:]), "Malformed generator table row")
    require(len(rows) <= 4, "Generator exceeded three-output cap")
    index = TABLE_COLUMNS.index("Final structure")
    raw = [r[index] for r in rows[1:]]
    require(all(raw), "Missing final design SMILES")
    result = {"columns": rows[0], "rows": rows[1:], "raw_designs": raw, "sha256": sha(path)}
    try:
        result.update(canonical_selection(raw, max_designs))
        result["downstream_selection_status"] = "success"
    except Exception as exc:
        result.update(
            downstream_selection_status="error", downstream_selection_error=f"{type(exc).__name__}: {exc}"
        )
    return result


def check_generator_manifest(manifest, job, index, source_hashes, runtime_record):
    expected = {
        k: v for k, v in source_hashes.items() if "/" not in k or k.startswith("chembounce_synthesis/")
    }
    require(manifest.get("source_code_sha256") == expected, "Archived generator code differs")
    require(manifest.get("source_files_changed_after_import") == [], "Generator source changed after import")
    require(manifest.get("python") == runtime_record["python"], "Generator Python differs")
    require(
        manifest.get("dependencies")
        == {k: runtime_record["dependencies"][k] for k in DEPENDENCIES if k != "psutil"},
        "Generator dependencies differ",
    )
    params = manifest["parameters"]
    for k, value in {
        **DEFAULTS,
        **job["options"],
        "seed": job["seed"],
        "input_smiles": job["query_smiles"],
    }.items():
        require(params.get(k) == value, f"Generator parameter differs: {k}")
    require(
        params.get("effective_limits")
        == {"overall": 3, "per_fragment": 10, "scaffolds": 8, "per_replacement": 10, "beam": 100},
        "Generator effective limits differ",
    )
    require(
        str(Path(params["fingerprint_db_file"]).resolve()) == index["path"], "Generator index path differs"
    )
    require(sha(params["scaffold_db_file"]) == job["scaffold_sha256"], "Generator scaffold source differs")
    require(
        manifest.get("sources") == {index["path"]: index["files_sha256"]}, "Generator index hashes differ"
    )


def inspect_generator(directory, job, index, source_hashes, runtime_record):
    directory = Path(directory)
    manifest, summary = read_json(directory / "run_manifest.json"), read_json(directory / "run_summary.json")
    check_generator_manifest(manifest, job, index, source_hashes, runtime_record)
    require(summary.get("status") in {"success", "no_candidates"}, "Generator did not finish successfully")
    require(
        summary.get("source_files_changed_after_import") == [], "Generator summary reports source mutation"
    )
    table = table_view(directory / "overall_result.txt", job["max_designs"])
    require((summary["status"] == "success") == bool(table["rows"]), "Generator status/table disagreement")
    files = {p.name: sha(p) for pattern in PATTERNS for p in sorted(directory.glob(pattern)) if p.is_file()}
    intermediates = {
        k: v
        for k, v in files.items()
        if k not in {"run_manifest.json", "run_summary.json", "overall_result.txt"}
    }
    normalized_summary = copy.deepcopy(summary)
    normalized_summary.get("database", {}).pop("path", None)
    return {
        "table": table,
        "status": summary["status"],
        "files_sha256": files,
        "intermediates_sha256": intermediates,
        "summary_without_database_path": normalized_summary,
    }


def stable_copy(source, destination):
    source, destination = Path(source), Path(destination)
    before = source.stat()
    data = source.read_bytes()
    after = source.stat()
    require(
        (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns),
        "Reference changed during snapshot",
    )
    digest = hashlib.sha256(data).hexdigest()
    require(sha(source) == digest, "Reference changed during snapshot read")
    with destination.open("xb") as stream:
        stream.write(data)
    destination.chmod(0o444)
    return {"source_path": str(source.resolve()), "sha256": digest, "bytes": len(data)}


def snapshot_references(jobs, output, hashes, runtime_record):
    references = []
    for job in jobs:
        for member in job["members"]:
            directory = Path(member["run_directory"])
            record = {
                "job_id": job["id"],
                "member": member,
                "status": "absent",
                "original_campaign_status": "not_observed",
            }
            evaluation = directory / "evaluation.json"
            if evaluation.is_file():
                data = evaluation.read_bytes()
                record["original_campaign_status"] = json.loads(data)["status"]
                record["original_evaluation_sha256"] = hashlib.sha256(data).hexdigest()
            if not directory.exists():
                references.append(record)
                continue
            record["status"] = "incomplete"
            required = [directory / "request.json"] + [
                directory / "designs" / p
                for p in ("run_manifest.json", "run_summary.json", "overall_result.txt")
            ]
            if not all(p.is_file() for p in required):
                references.append(record)
                continue
            if read_json(directory / "designs/run_summary.json").get("status") not in {
                "success",
                "no_candidates",
            }:
                record["status"] = "generator_failed"
                references.append(record)
                continue
            ident = object_sha([job["id"], str(directory)])
            dest = output / "references" / ident
            (dest / "designs").mkdir(parents=True)
            copied = {"request.json": stable_copy(directory / "request.json", dest / "request.json")}
            for pattern in PATTERNS:
                for path in sorted((directory / "designs").glob(pattern)):
                    if path.is_file():
                        copied["designs/" + path.name] = stable_copy(path, dest / "designs" / path.name)
            request = read_json(dest / "request.json")
            require(
                request["query"]["smiles"] == job["query_smiles"]
                and request["query"]["id"] == member["query"]["id"],
                "Archive query differs",
            )
            require(
                request["seed"] == job["seed"]
                and request["kind"] == "pipeline"
                and request["method_id"] == member["method_id"],
                "Archive condition differs",
            )
            require(request["settings"] == member["settings"], "Archive settings differ")
            require(
                request["budget"]["max_candidates"] == 3 and request["budget"]["wall_seconds"] == 360,
                "Archive budget differs",
            )
            require(
                request["assets"]["fingerprint_db"] == member["old_index"]["path"]
                and request["assets"]["scaffold_db"] == member["scaffold"]["path"],
                "Archive assets differ",
            )
            view = inspect_generator(dest / "designs", job, member["old_index"], hashes, runtime_record)
            synthesis_manifest = directory / "synthesis/manifest.json"
            if synthesis_manifest.is_file():
                copied["synthesis_manifest.json"] = stable_copy(
                    synthesis_manifest, dest / "synthesis_manifest.json"
                )
                actual_designs = read_json(dest / "synthesis_manifest.json")["designs"]
                require(
                    view["table"].get("downstream_selection_status") == "success"
                    and actual_designs == view["table"]["downstream_selected"],
                    "Archived synthesis selection differs from frozen canonical/dedup/cap semantics",
                )
            record.update(
                status="complete_generator_reference",
                snapshot=str(dest.relative_to(output)),
                files=copied,
                view=view,
            )
            references.append(record)
    return references


def execution_environment(root):
    env = os.environ.copy()
    env.update(THREADS)
    env.update(CUDA_VISIBLE_DEVICES="", PYTHONPATH=str(Path(root).resolve()), PYTHONDONTWRITEBYTECODE="1")
    return env


def environment_record(env):
    keys = set(THREADS) | {
        "CUDA_VISIBLE_DEVICES",
        "PYTHONPATH",
        "PYTHONHASHSEED",
        "PYTHONDONTWRITEBYTECODE",
        "LC_ALL",
        "LANG",
        "TZ",
    }
    return {key: env.get(key) for key in sorted(keys)}


def protect_output(output, archives, indices=(), plan_path=None):
    output = Path(output).resolve()
    protected = [Path(p).resolve() for p in archives] + [Path(p).resolve() for p in indices]
    if plan_path is not None:
        protected.append(Path(plan_path).resolve().parent)
    require(
        not any(output.is_relative_to(p) for p in protected),
        "Output would modify an archive, index, or immutable plan snapshot",
    )


def plan(args):
    output = Path(args.output).resolve()
    require(not output.exists(), "Plan output already exists; no overwrite")
    require(
        os.path.abspath(args.python) == os.path.abspath(sys.executable),
        "Create plan using the exact --python interpreter",
    )
    frozen = read_json(args.frozen_manifest)
    hashes = frozen["code_sha256"]
    require(len(hashes) == frozen["file_count"] == 39, "Expected frozen 39-source manifest")
    install_source(args.source_root, hashes)
    current_runtime = runtime()
    jobs, inputs = build_grid(args.campaign)
    index = directory_pin(args.new_index)
    old_index, scaffold = jobs[0]["members"][0]["old_index"], jobs[0]["members"][0]["scaffold"]
    require(index["path"] != old_index["path"], "Counterfactual must use a separate index")
    report = read_json(args.index_validation)
    metadata = validate_index_report(index, scaffold["sha256"], current_runtime, report)
    protect_output(output, [i["archive"] for i in inputs], [index["path"], old_index["path"]])
    output.mkdir(parents=True)
    references = snapshot_references(jobs, output, hashes, current_runtime)
    for name, path in (
        ("frozen-source-manifest.json", args.frozen_manifest),
        ("index-validation.json", args.index_validation),
    ):
        stable_copy(path, output / name)
    record = {
        "schema": SCHEMA,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "helper": file_pin(__file__),
        "source_root": str(Path(args.source_root).resolve()),
        "source_sha256": hashes,
        "frozen_manifest": file_pin(output / "frozen-source-manifest.json"),
        "implementation_commit": frozen["implementation_commit"],
        "python_executable": os.path.abspath(args.python),
        "runtime": current_runtime,
        "environment": environment_record(execution_environment(args.source_root)),
        "indices": {"old": old_index, "new": index},
        "scaffold": scaffold,
        "new_index_metadata": metadata,
        "index_validation": file_pin(output / "index-validation.json"),
        "inputs": inputs,
        "jobs": jobs,
        "references": references,
        "counts": {
            "conditions": len(jobs),
            "arm_calls": 2 * len(jobs),
            "campaign_members": sum(len(j["members"]) for j in jobs),
        },
        "execution": {
            "wall_seconds": 360.0,
            "max_memory_mb": 24000,
            "memory_sample_seconds": 0.05,
            "pair_order": "job-id hash parity in dispatch queue; fresh subprocess per arm; parallel dispatch may overlap paired arms",
            "parallelism": args.workers,
        },
        "limitations": LIMITATIONS,
    }
    atomic_json(output / "plan.json", record)
    atomic_json(output / "plan.sha256.json", {"sha256": sha(output / "plan.json")})
    return {
        "plan": str(output / "plan.json"),
        "sha256": sha(output / "plan.json"),
        "counts": record["counts"],
        "reference_statuses": dict(collections.Counter(r["status"] for r in references)),
    }


def verify_plan(path, expected_sha):
    path = Path(path).resolve()
    require(sha(path) == expected_sha, "Plan hash mismatch")
    data = read_json(path)
    require(data["schema"] == SCHEMA, "Unknown plan schema")
    check_pin(data["helper"])
    require(sha(__file__) == data["helper"]["sha256"], "Executing a different helper")
    install_source(data["source_root"], data["source_sha256"])
    require(os.path.abspath(sys.executable) == data["python_executable"], "Interpreter changed")
    require(runtime() == data["runtime"], "Runtime changed after plan")
    require(
        environment_record(execution_environment(data["source_root"])) == data["environment"],
        "Execution environment changed after plan",
    )
    for pin in (data["frozen_manifest"], data["index_validation"], data["scaffold"]):
        check_pin(pin)
    for index in data["indices"].values():
        check_directory(index)
    validate_index_report(
        data["indices"]["new"],
        data["scaffold"]["sha256"],
        data["runtime"],
        read_json(data["index_validation"]["path"]),
    )
    for item in data["inputs"]:
        check_pin(item["dataset"])
        check_pin(item["protocol"])
    for ref in data["references"]:
        if ref["status"] == "complete_generator_reference":
            for relative, pin in ref["files"].items():
                require(
                    sha(path.parent / ref["snapshot"] / relative) == pin["sha256"],
                    "Archived reference snapshot changed",
                )
    require(
        len(data["jobs"]) == data["counts"]["conditions"] == 126 and data["counts"]["arm_calls"] == 252,
        "Grid count changed",
    )
    return data


def stop_group(process):
    # Kill descendants even when the immediate child already exited.
    cleanup = {
        "term_sent": False,
        "kill_sent": False,
        "reaped": False,
        "term_wait_seconds": 0.5,
        "kill_wait_seconds": 1.0,
    }
    try:
        os.killpg(process.pid, signal.SIGTERM)
        cleanup["term_sent"] = True
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=0.5)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
        cleanup["kill_sent"] = True
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=1.0)
        cleanup.update(reaped=True, status="reaped")
    except subprocess.TimeoutExpired:
        cleanup["status"] = "unreaped_requires_parent_review"
    return cleanup


def run_process(command, cwd, env, directory, wall_seconds, max_memory_mb, sample_seconds=0.05, cancel=None):
    import psutil

    started = time.monotonic()
    samples, peak = 0, 0.0
    status = None
    with (directory / "worker.log").open("xb") as log:
        process = subprocess.Popen(
            command, cwd=cwd, env=env, start_new_session=True, stdout=log, stderr=subprocess.STDOUT
        )
        try:
            while process.poll() is None:
                if cancel is not None and cancel.is_set():
                    status = "infrastructure_interruption"
                    break
                if time.monotonic() - started >= wall_seconds:
                    status = "timeout"
                    break
                try:
                    root = psutil.Process(process.pid)
                    rss = 0
                    for p in [root] + root.children(recursive=True):
                        try:
                            rss += p.memory_info().rss
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            pass
                    peak = max(peak, rss / 1024**2)
                    samples += 1
                    if rss / 1024**2 > max_memory_mb:
                        status = "memory_limit"
                        break
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
                time.sleep(min(sample_seconds, max(0.0, wall_seconds - (time.monotonic() - started))))
            if status is None:
                status = "completed" if process.returncode == 0 else "worker_error"
        finally:
            cleanup = stop_group(process)
    return {
        "status": status if cleanup["reaped"] else "cleanup_unreaped",
        "status_before_cleanup": status,
        "cleanup": cleanup,
        "returncode": process.returncode,
        "resources": {
            "wall_seconds": time.monotonic() - started,
            "sampled_group_peak_rss_mb": peak if samples else None,
            "memory_samples": samples,
            "sampling_interval_seconds": sample_seconds,
            "timing_scope": "subprocess start through cleanup; imports/provenance/frozen generation/finalization included",
            "memory_scope": "sampled sum of live worker/descendant RSS; shared pages may be counted more than once",
        },
    }


def worker(request_path):
    request = read_json(request_path)
    out = Path(request["output"])
    require(not (out / "worker-result.json").exists(), "Worker result already exists")
    try:
        install_source(request["source_root"], request["source_sha256"])
        require(runtime() == request["runtime"], "Worker runtime differs")
        require(environment_record(os.environ) == request["environment"], "Worker environment differs")
        job = request["job"]
        require(job["options"]["calculate_shape"] is False, "3D generation is prohibited")
        from chembounce import chembounce

        require(
            Path(sys.modules["chembounce"].__file__).resolve()
            == Path(request["source_root"]) / "chembounce.py",
            "Wrong generator import",
        )
        options = dict(job["options"])
        options.update(
            output_dir=out / "designs",
            seed=job["seed"],
            overall_max_n=3,
            scaffold_db_file=request["scaffold"]["path"],
            fingerprint_db_file=request["index"]["path"],
        )
        frame, cost = chembounce(job["query_smiles"], **options)
        view = inspect_generator(
            out / "designs", job, request["index"], request["source_sha256"], request["runtime"]
        )
        require(
            list(frame["Final structure"]) == view["table"]["raw_designs"],
            "Returned frame differs from saved designs",
        )
        source_check(request["source_root"], request["source_sha256"])
        prohibited = {
            "torch",
            "transformers",
            "chembounce_synthesis.workflow",
            "chembounce_synthesis.synllama",
            "chembounce_synthesis.evidence",
            "chembounce_synthesis.conformers",
        }
        require(not prohibited.intersection(sys.modules), "Prohibited inference/workflow module imported")
        atomic_json(
            out / "worker-result.json",
            {
                "status": "success",
                "generator_cost": cost,
                "view": view,
                "no_prohibited_modules_imported": True,
            },
        )
        return 0
    except Exception as exc:
        atomic_json(
            out / "worker-result.json",
            {"status": "error", "error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()},
        )
        return 1


def compare_views(left, right):
    a, b = left["table"], right["table"]
    selected_assessed = a["downstream_selection_status"] == b["downstream_selection_status"] == "success"
    return {
        "raw_table_bytes_equal": a["sha256"] == b["sha256"],
        "full_table_cells_equal": a["columns"] == b["columns"] and a["rows"] == b["rows"],
        "raw_ordered_designs_equal": a["raw_designs"] == b["raw_designs"],
        "canonical_ordered_designs_equal": a.get("canonical_deduplicated") == b.get("canonical_deduplicated")
        if selected_assessed
        else None,
        "downstream_first_two_equal": a.get("downstream_selected") == b.get("downstream_selected")
        if selected_assessed
        else None,
        "downstream_comparison_status": "assessed" if selected_assessed else "unassessed",
        "intermediate_files_equal": left["intermediates_sha256"] == right["intermediates_sha256"],
        "summary_excluding_index_path_equal": left["summary_without_database_path"]
        == right["summary_without_database_path"],
    }


def summarize(data, reports):
    lookup = {(r["job_id"], r["arm"]): r for r in reports}
    require(len(lookup) == len(reports), "Duplicate condition/arm report")
    rows = []
    for job in data["jobs"]:
        a, b = lookup.get((job["id"], "old")), lookup.get((job["id"], "new"))
        row = {
            "job_id": job["id"],
            "policy": job["options"]["search_policy"],
            "seed": job["seed"],
            "old_status": a["status"] if a else "unattempted",
            "new_status": b["status"] if b else "unattempted",
            "comparison_status": "unassessed",
        }
        if a and b and a["status"] == b["status"] == "success":
            row.update(comparison_status="assessed", comparison=compare_views(a["view"], b["view"]))
        refs = [
            r
            for r in data["references"]
            if r["job_id"] == job["id"] and r["status"] == "complete_generator_reference"
        ]
        row["archived_reference_count"] = len(refs)
        row["old_replay_status"] = "no_complete_reference" if not refs else "unassessed"
        if refs and a and a["status"] == "success":
            comparisons = [compare_views(a["view"], r["view"]) for r in refs]
            row["old_replay_comparisons"] = comparisons
            row["old_replay_status"] = (
                "matched"
                if all(
                    c["raw_table_bytes_equal"]
                    and c["downstream_first_two_equal"] is True
                    and c["intermediate_files_equal"]
                    and c["summary_excluding_index_path_equal"]
                    for c in comparisons
                )
                else "mismatch_or_unassessed"
            )
        row["asset_only_attribution_status"] = (
            "replay_mismatch_blocks_attribution"
            if row["old_replay_status"] == "mismatch_or_unassessed"
            else "conditional_on_pinned_runtime_and_code"
        )
        rows.append(row)
    return {
        "schema": SCHEMA,
        "planned_conditions": len(data["jobs"]),
        "planned_arm_calls": len(data["jobs"]) * 2,
        "recorded_arm_calls": len(reports),
        "unattempted_arm_calls": len(data["jobs"]) * 2 - len(reports),
        "arm_statuses": dict(collections.Counter(r["status"] for r in reports)),
        "condition_statuses": dict(collections.Counter(r["comparison_status"] for r in rows)),
        "old_replay_statuses": dict(collections.Counter(r["old_replay_status"] for r in rows)),
        "conditions": rows,
        "limitations": LIMITATIONS,
    }


def run(args):
    output = Path(args.output).resolve()
    require(not output.exists(), "Execution output already exists; no resume or overwrite")
    data = verify_plan(args.plan, args.plan_sha256)
    protect_output(
        output,
        [i["archive"] for i in data["inputs"]],
        [i["path"] for i in data["indices"].values()],
        args.plan,
    )
    output.mkdir(parents=True)
    atomic_json(
        output / "run-manifest.json",
        {
            "schema": SCHEMA,
            "plan": file_pin(args.plan),
            "helper": file_pin(__file__),
            "runtime": runtime(),
            "limitations": LIMITATIONS,
        },
    )
    reports = []
    env = execution_environment(data["source_root"])
    cancel = threading.Event()
    dispatch = []
    for job in data["jobs"]:
        arms = ("old", "new") if int(job["id"][-1], 16) % 2 == 0 else ("new", "old")
        dispatch.extend((job, arm) for arm in arms)
    require(data["execution"]["parallelism"] in (1, 2), "Only one or two CPU workers are supported")

    def execute(item):
        job, arm = item
        directory = output / job["id"] / arm
        directory.mkdir(parents=True)
        request = {
            "schema": SCHEMA,
            "job": {k: v for k, v in job.items() if k != "members"},
            "arm": arm,
            "source_root": data["source_root"],
            "source_sha256": data["source_sha256"],
            "runtime": data["runtime"],
            "environment": data["environment"],
            "index": data["indices"][arm],
            "scaffold": data["scaffold"],
            "output": str(directory),
        }
        request_path = directory / "request.json"
        atomic_json(request_path, request)
        command = [
            data["python_executable"],
            str(Path(__file__).resolve()),
            "_worker",
            "--request",
            str(request_path),
        ]
        report = None
        try:
            report = run_process(
                command,
                data["source_root"],
                env,
                directory,
                wall_seconds=data["execution"]["wall_seconds"],
                max_memory_mb=data["execution"]["max_memory_mb"],
                sample_seconds=data["execution"]["memory_sample_seconds"],
                cancel=cancel,
            )
            result_path = directory / "worker-result.json"
            if result_path.exists():
                report["worker_result_sha256"] = sha(result_path)
                result = read_json(result_path)
                if report["status"] == "completed" and result.get("status") == "success":
                    # Independently reread artifacts rather than trusting worker's copied view.
                    view = inspect_generator(
                        directory / "designs",
                        job,
                        data["indices"][arm],
                        data["source_sha256"],
                        data["runtime"],
                    )
                    require(view == result["view"], "Worker result view differs from artifacts")
                    report.update(status="success", view=view, generator_cost=result["generator_cost"])
                else:
                    report["worker_result"] = result
            elif report["status"] == "completed":
                report["status"] = "missing_worker_result"
        except Exception as exc:
            report = report or {"resources": {"wall_seconds": None, "measurement_status": "not_available"}}
            report.update(status="infrastructure_error", error=f"{type(exc).__name__}: {exc}")
        report.update(job_id=job["id"], arm=arm, command=command, request_sha256=sha(request_path))
        atomic_json(directory / "report.json", report)
        return report

    def record(report):
        reports.append(report)
        atomic_json(output / "runs.json", reports)
        atomic_json(output / "summary.json", summarize(data, reports))
        atomic_json(
            output / "progress.json",
            {"completed_arm_calls": len(reports), "planned_arm_calls": 252, "state": "running"},
        )
        print(f"{len(reports)}/252 {report['job_id'][:12]} {report['arm']}: {report['status']}", flush=True)

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=data["execution"]["parallelism"])
    pending, cursor = set(), 0
    previous_sigterm = signal.getsignal(signal.SIGTERM)

    def interrupted(signum, frame):
        raise KeyboardInterrupt(f"Controller received signal {signum}")

    signal.signal(signal.SIGTERM, interrupted)
    try:
        atomic_json(
            output / "progress.json", {"completed_arm_calls": 0, "planned_arm_calls": 252, "state": "running"}
        )
        while cursor < len(dispatch) or pending:
            while cursor < len(dispatch) and len(pending) < data["execution"]["parallelism"]:
                pending.add(executor.submit(execute, dispatch[cursor]))
                cursor += 1
            done, _ = concurrent.futures.wait(
                pending, timeout=0.2, return_when=concurrent.futures.FIRST_COMPLETED
            )
            for future in done:
                pending.remove(future)
                report = future.result()
                record(report)
                if report.get("cleanup", {}).get("reaped") is False:
                    cancel.set()
                    raise RuntimeError(
                        "Worker could not be reaped after bounded TERM/KILL; dispatch stopped for parent review"
                    )
        # Recheck all original inputs/assets/snapshot hashes after the diagnostic.
        verify_plan(args.plan, args.plan_sha256)
        require(
            len(reports) == 252
            and {(r["job_id"], r["arm"]) for r in reports}
            == {(j["id"], a) for j in data["jobs"] for a in ("old", "new")},
            "Incomplete final grid",
        )
        atomic_json(
            output / "progress.json",
            {"state": "complete", "completed_arm_calls": 252, "planned_arm_calls": 252},
        )
    except BaseException as exc:
        cancel.set()
        for future in pending:
            try:
                record(future.result())
            except Exception:
                pass  # Partial output remains reviewable; summary retains unassessed calls.
        atomic_json(
            output / "progress.json",
            {
                "state": "interrupted_requires_review",
                "completed_arm_calls": len(reports),
                "planned_arm_calls": 252,
                "error": f"{type(exc).__name__}: {exc}",
            },
        )
        atomic_json(output / "summary.json", summarize(data, reports))
        raise
    finally:
        executor.shutdown(wait=True)
        signal.signal(signal.SIGTERM, previous_sigterm)
    return {"output": str(output), "recorded_arm_calls": len(reports)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("plan")
    p.add_argument(
        "--campaign", nargs=3, action="append", required=True, metavar=("DATASET", "PROTOCOL", "ARCHIVE")
    )
    for name in ("source-root", "frozen-manifest", "python", "new-index", "index-validation", "output"):
        p.add_argument("--" + name, required=True)
    p.add_argument("--workers", type=int, choices=(1, 2), default=1)
    p = sub.add_parser("run")
    for name in ("plan", "plan-sha256", "output"):
        p.add_argument("--" + name, required=True)
    p = sub.add_parser("_worker", help=argparse.SUPPRESS)
    p.add_argument("--request", required=True)
    args = parser.parse_args()
    if args.action == "_worker":
        return worker(args.request)
    print(json.dumps(plan(args) if args.action == "plan" else run(args), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
