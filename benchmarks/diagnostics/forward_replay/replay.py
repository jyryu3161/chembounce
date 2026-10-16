"""Replay normalized, independently reviewed common-forward statistics; no inference."""

import argparse
import ast
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import sys
from types import ModuleType

HERE = Path(__file__).resolve().parent
HELPER_SHA = "01b4c9468d23bb2870083c0b8f66886fafbf0262b429e08e5e102e683e20d83e"
PLAN_SHA = "385be03fee3d2df6cfc04975cb63b42717088eabe62077bdbd3c6db9ca5db529"
EXTRACTED_SHA = "0bf96138f0a137a6048a700d78bc064f6b278766d8b09da87f3f54ca45739cd1"
EXTRACTION_SHA = "d9c53088f3f5cdf04fdbeecc4851a55a6517dc3088311740db09992178ca2e27"
TRANSPORT_SHA = "9751209528300cdc38250d68c3a6a8748ed873362af1d1dabfa14ecbb2234113"
COMMIT = "d972124d587e647ef48455cf1a038e6be96e5b04"
SPLITS = {"scaffold": 360, "time": 180, "target": 90}
DATA = {
    "attempts.json",
    "normalized-occurrences.json",
    "integrity-validation.json",
    "recorded-analysis.json",
    "recorded-validation.json",
    "audit-protocol.json",
    "audit-progress.json",
}
NODES = {
    "SPLITS",
    "PRIMARY",
    "INFRA",
    "COUNTS",
    "COVERAGE",
    "BOUNDED",
    "CONDITIONAL",
    "METRICS",
    "TERMINAL",
    "ValidationError",
    "require",
    "run_key",
    "mean",
    "bootstrap",
    "summarize",
    "paired_rows",
    "aggregate_split",
}
SCOPE = (
    "CPU statistics from parent-validated normalized original-630 attempts and common-forward occurrences. "
    "No chemistry, route replay, model inference, raw-prediction reinterpretation, or supplementary-270 analysis. "
    "Original validation/provenance is retained, not re-executed. Same forward model as selection gate; "
    "selected-subset agreement is not independent synthesis evidence, activity, yield, or synthesis probability. "
    "No-route-supplied is coverage absence, not chemical synthesis failure. Exact query/seed paired descriptive "
    "target-cluster uncertainty retains few-target, reused-input, and cross-target dependence limitations."
)


def module_bytes(raw, name, path):
    module = ModuleType(name)
    module.__file__ = str(path)
    exec(compile(raw, str(path), "exec"), module.__dict__)
    return module


transport_path = HERE.parent / "generation_replay/replay.py"
transport_raw = transport_path.read_bytes()
if hashlib.sha256(transport_raw).hexdigest() != TRANSPORT_SHA:
    raise ValueError("Reviewed generation-replay transport implementation changed")
io = module_bytes(transport_raw, "forward_replay_transport", transport_path)


def fixed_files():
    return {
        "analyze_forward_historical.py": (HERE / "analyze_forward_historical.py", HELPER_SHA),
        "statistics_extracted.py": (HERE / "statistics_extracted.py", EXTRACTED_SHA),
        "extraction.json": (HERE / "extraction.json", EXTRACTION_SHA),
        "forward-analysis-plan.json": (HERE / "forward-analysis-plan.json", PLAN_SHA),
        "transport_helpers.py": (transport_path, TRANSPORT_SHA),
        io.PROVENANCE: (io.HERE / io.PROVENANCE, io.PROVENANCE_SHA),
        io.ATTRIBUTION: (io.HERE / io.ATTRIBUTION, io.ATTRIBUTION_SHA),
    }


def extract_statistics(raw):
    """Copy only allowlisted original AST node bytes; never execute the full helper."""
    lines = raw.splitlines(keepends=True)
    parts, records, found = [], [], set()
    for node in ast.parse(raw).body:
        label = None
        if isinstance(node, ast.ImportFrom) and node.module == "collections":
            label = "import:collections"
        elif (
            isinstance(node, ast.Import)
            and len(node.names) == 1
            and node.names[0].name in {"statistics", "numpy"}
        ):
            label = "import:" + node.names[0].name
        elif isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in NODES:
            label = node.name
        elif (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id in NODES
        ):
            label = node.targets[0].id
        if label:
            io.require(label not in found, "Duplicate extraction node")
            found.add(label)
            part = b"".join(lines[node.lineno - 1 : node.end_lineno])
            parts.append(part)
            records.append(
                {
                    "node": label,
                    "first_line": node.lineno,
                    "last_line": node.end_lineno,
                    "sha256": io.sha(part),
                }
            )
    io.require(
        found == NODES | {"import:collections", "import:statistics", "import:numpy"},
        "Extraction node set mismatch",
    )
    return b"\n".join(parts), records


def check_fixed(blobs):
    for name, (_, digest) in fixed_files().items():
        io.same(io.sha(blobs[name]), digest, "fixed artifact " + name)
    raw, records = extract_statistics(blobs["analyze_forward_historical.py"])
    io.same(raw, blobs["statistics_extracted.py"], "original/extracted bytes")
    spec = io.parse(blobs["extraction.json"])
    io.same(spec["selected_nodes"], records, "extraction boundaries")
    io.same(spec["historical_sha256"], HELPER_SHA, "extraction original SHA")
    io.same(spec["extracted_sha256"], EXTRACTED_SHA, "extraction SHA")


def content_key(value):
    return io.sha(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode())


def historical_source_hash(mapping, relative):
    """Compare recorded path suffixes only; never dereference historical paths."""
    io.safe_relative(relative)
    matches = [
        digest for path, digest in mapping.items() if path == relative or path.endswith("/" + relative)
    ]
    io.require(len(matches) == 1, "Missing/ambiguous historical hash: " + relative)
    return io.sha_value(matches[0])


def parent_review(blobs):
    review = io.parse(blobs["parent-review.json"])
    for name, value in {
        "schema_version": 1,
        "status": "passed",
        "total_attempts": 630,
        "expected_attempts": SPLITS,
        "analysis_helper_sha256": HELPER_SHA,
        "analysis_plan_sha256": PLAN_SHA,
        "extracted_statistics_sha256": EXTRACTED_SHA,
        "original_source_commit": COMMIT,
        "normalization_source_validation_sha256": io.sha(blobs["recorded-validation.json"]),
        "input_sha256": {n: io.sha(blobs[n]) for n in DATA},
    }.items():
        io.same(review[name], value, "parent review " + name)
    evidence = review["assessment_evidence"]
    io.require(isinstance(evidence, list) and evidence, "Actual final parent assessment evidence required")
    seen = set()
    for item in evidence:
        io.require(set(item) == {"artifact", "sha256"}, "Invalid parent evidence entry")
        name = io.safe_relative(item["artifact"])
        io.require(name not in seen, "Duplicate parent evidence")
        seen.add(name)
        raw = blobs["parent-review-evidence/" + name]
        io.same(io.sha(raw), io.sha_value(item["sha256"]), "parent assessment evidence hash")
        io.require(isinstance(io.parse(raw), dict), "Parent assessment must be a JSON object")
    io.require(
        set(blobs)
        == DATA | set(fixed_files()) | {"parent-review.json"} | {"parent-review-evidence/" + p for p in seen},
        "Bundle artifact membership mismatch",
    )


def validate_normalized(attempts, occurrences, plan):
    io.require(
        isinstance(attempts, list) and isinstance(occurrences, list), "Normalized tables must be arrays"
    )
    attempts_by_split, occurrences_by_split = defaultdict(list), defaultdict(list)
    indexed, caps, targets, expected = {}, {}, {}, set()
    for split, settings in plan["expected_grid"].items():
        io.require(split in SPLITS, "Unexpected split")
        targets[split] = settings["queries"]
        caps[split] = {
            m: cap for spec in settings["groups"].values() for m, cap in spec["output_caps"].items()
        }
        for spec in settings["groups"].values():
            io.same(spec["seeds"], [42, 43, 44], "frozen seeds")
            for method in spec["methods"]:
                io.same(caps[split][method], 3, "frozen output cap")
                for query in targets[split]:
                    for seed in spec["seeds"]:
                        key = split, method, query, seed
                        io.require(key not in expected, "Duplicate planned attempt")
                        expected.add(key)
    terminals = {
        "success",
        "no_candidates",
        "timeout",
        "memory_limit",
        "budget_exceeded",
        "worker_error",
        "missing_response",
        "invalid_response",
        "launch_error",
        "interrupted",
        "infrastructure_interruption",
    }
    for row in attempts:
        io.require(
            set(row) == {"split", "method_id", "query_id", "target_id", "seed", "status"},
            "Attempt fields mismatch",
        )
        io.require(type(row["seed"]) is int, "Generation seed must be integer")
        key = row["split"], row["method_id"], row["query_id"], row["seed"]
        io.require(key in expected and key not in indexed, "Unexpected/duplicate attempt")
        io.same(row["target_id"], targets[row["split"]][row["query_id"]], "attempt target")
        io.require(row["status"] in terminals, "Nonterminal/unknown attempt status")
        indexed[key] = row
        attempts_by_split[row["split"]].append(row)
    io.require(indexed.keys() == expected and len(expected) == 630, "Incomplete original 630-attempt grid")
    io.same({s: len(rows) for s, rows in attempts_by_split.items()}, SPLITS, "attempt split counts")
    ordinal, cached, seen = defaultdict(list), {}, set()
    required = {
        "split",
        "method_id",
        "query_id",
        "target_id",
        "generation_seed",
        "rank",
        "smiles",
        "status",
        "forward_all_steps_matched",
        "matched",
        "evaluated",
        "input_keys",
        "route_key",
    }
    for row in occurrences:
        io.require(isinstance(row, dict) and required <= row.keys(), "Missing normalized occurrence fields")
        io.require(
            type(row["generation_seed"]) is int and type(row["rank"]) is int and row["rank"] > 0,
            "Invalid occurrence seed/ordinal",
        )
        key = row["split"], row["method_id"], row["query_id"], row["generation_seed"]
        io.require(
            key in indexed and indexed[key]["status"] == "success", "Occurrence without successful attempt"
        )
        io.require((*key, row["rank"]) not in seen, "Duplicate normalized occurrence")
        seen.add((*key, row["rank"]))
        ordinal[key].append(row["rank"])
        io.same(row["target_id"], indexed[key]["target_id"], "occurrence target")
        io.require(isinstance(row["smiles"], str) and row["smiles"], "Missing stored candidate identity")
        io.require(
            type(row["matched"]) is bool and type(row["evaluated"]) is bool, "Outcome flags must be booleans"
        )
        io.require(not row["matched"] or row["evaluated"], "Matched outcome must be evaluated")
        io.require(isinstance(row["input_keys"], list), "Input keys must be a list")
        for value in row["input_keys"]:
            io.sha_value(value)
        status = row["status"]
        io.require(status in {"no_route_supplied", "assessed", "audit_failed"}, "Unknown occurrence status")
        if status == "no_route_supplied":
            io.require(
                row["route_key"] is None and "route_cache_key" not in row,
                "No-route promoted to supplied route",
            )
        else:
            io.sha_value(row["route_key"])
            io.same(row["route_cache_key"], row["route_key"], "normalized/original route key")
        if status == "assessed":
            io.same(row["forward_all_steps_matched"], row["matched"], "original matched flag")
            io.same(row["all_steps_evaluated"], row["evaluated"], "original evaluated flag")
            io.require(bool(row["input_keys"]), "Assessed route has no step inputs")
        else:
            io.require(
                not row["matched"]
                and not row["evaluated"]
                and row["input_keys"] == []
                and row["forward_all_steps_matched"] is None
                and row.get("all_steps_evaluated") is None,
                "Unassessed route promoted to evidence",
            )
            if status == "audit_failed":
                io.require(
                    isinstance(row.get("error"), str) and row["error"], "Audit failure explanation missing"
                )
        if row["route_key"]:
            value = {k: row[k] for k in ("status", "smiles", "matched", "evaluated", "input_keys")}
            if row["route_key"] in cached:
                io.same(value, cached[row["route_key"]], "reused route outcomes")
            cached[row["route_key"]] = value
        occurrences_by_split[row["split"]].append(row)
    for key, ranks in ordinal.items():
        io.require(
            sorted(ranks) == list(range(1, len(ranks) + 1)) and len(ranks) <= caps[key[0]][key[1]],
            "Noncontiguous valid-output ordinals or exceeded cap",
        )
    return attempts_by_split, occurrences_by_split, targets, caps


def validate(blobs, calculate=False):
    check_fixed(blobs)
    parent_review(blobs)
    data = {n: io.parse(blobs[n]) for n in DATA}
    plan = io.parse(blobs["forward-analysis-plan.json"])
    io.same(plan["expected_attempts"], SPLITS, "plan count")
    io.same(set(plan["expected_grid"]), set(SPLITS), "plan split set")
    io.same(io.sha(blobs["audit-protocol.json"]), plan["audit_protocol_sha256"], "frozen audit protocol")
    protocol, progress, validated = (
        data[n] for n in ("audit-protocol.json", "audit-progress.json", "recorded-validation.json")
    )
    io.same(protocol["expected_attempts"], SPLITS, "audit protocol counts")
    io.same(progress["complete"], True, "common-forward completion")
    io.same(progress["observed_attempts"], SPLITS, "observed attempt counts")
    for key, expected in {
        "passed": True,
        "total_attempts": 630,
        "analysis_helper_sha256": HELPER_SHA,
        "analysis_plan_sha256": PLAN_SHA,
        "scientific_comparison_emitted": True,
    }.items():
        io.same(validated[key], expected, "original validation " + key)
    io.same(
        validated["analysis_runtime"]["forward_inference_performed"],
        False,
        "original analysis inference flag",
    )
    io.same(validated["audit_model_provenance"], protocol["model"], "original model provenance")
    original_hashes = validated["source_sha256"]
    io.require(isinstance(original_hashes, dict) and original_hashes, "Original source hashes missing")
    for value in original_hashes.values():
        io.sha_value(value)
    for name, value in plan["source_sha256"].items():
        io.same(historical_source_hash(original_hashes, name), value, "plan/original input hash")
    for role in ("protocol", "progress"):
        io.same(
            historical_source_hash(original_hashes, f"common-forward-audit-cpu/{role}.json"),
            io.sha(blobs[f"audit-{role}.json"]),
            "original audit input hash",
        )
    source = io.parse(blobs[io.PROVENANCE])
    roots = validated["resolved_inputs"]
    io.same(roots["frozen_code_file_count"], 39, "frozen source count")
    io.same(roots["frozen_code_inventory_sha256"], content_key(source["code_sha256"]), "frozen source map")
    for name, value in source["code_sha256"].items():
        key = str(PurePosixPath(roots["source_root"]) / name)
        io.same(original_hashes[key], value, "original source hash")
    expected_folders = {
        f"comparison-{group}" if split == "scaffold" else f"comparison-{split}-{group}"
        for split, settings in plan["expected_grid"].items()
        for group in settings["groups"]
    }
    io.same(set(progress["input_run_sha256"]), expected_folders, "audit run hash membership")
    for folder, value in progress["input_run_sha256"].items():
        io.same(historical_source_hash(original_hashes, folder + "/runs.json"), value, "audit run hash")
    attempts, occurrences = data["attempts.json"], data["normalized-occurrences.json"]
    rows, outputs, targets, caps = validate_normalized(attempts, occurrences, plan)
    io.same(progress["candidate_occurrences"], len(occurrences), "occurrence coverage")
    io.same(progress["status_counts"], dict(Counter(r["status"] for r in occurrences)), "occurrence statuses")
    unique_routes = len({r["route_key"] for r in occurrences if r["route_key"]})
    io.same(progress["unique_cached_routes"], unique_routes, "route cache coverage")
    integrity = data["integrity-validation.json"]
    for key, expected in {
        "passed": True,
        "attempts": 630,
        "candidate_occurrences": len(occurrences),
        "unique_routes": unique_routes,
        "outputs_sha256": {n: io.sha(blobs[n]) for n in ("attempts.json", "normalized-occurrences.json")},
    }.items():
        io.same(integrity[key], expected, "normalization integrity " + key)
    integrity_hashes = integrity["files_sha256_verified_before_after"]
    for folder, value in progress["input_run_sha256"].items():
        io.same(
            historical_source_hash(integrity_hashes, folder + "/runs.json"),
            value,
            "normalization raw run hash",
        )
    for role in ("protocol", "progress", "occurrences"):
        relative = f"common-forward-audit-cpu/{role}.json"
        io.same(
            historical_source_hash(integrity_hashes, relative),
            historical_source_hash(original_hashes, relative),
            "normalization/original raw audit input",
        )
    grids = {
        s: {"expected": count, "observed": len(rows[s]), "complete": True, "missing": [], "duplicates": []}
        for s, count in SPLITS.items()
    }
    io.same(validated["grid"], grids, "recorded full grid")
    expected = data["recorded-analysis.json"]
    meta = {
        "analysis_plan": plan,
        "validation": validated,
        **{k: plan[k] for k in ("interpretation", "dependence", "conditional_metric", "infrastructure")},
    }
    io.require(set(expected) == set(meta) | {"splits"}, "Recorded analysis field set mismatch")
    for key, value in meta.items():
        io.same(expected[key], value, "recorded analysis " + key)
    io.same(set(expected["splits"]), set(SPLITS), "analysis split membership")
    if not calculate:
        return None
    # Only now import NumPy and execute the pinned, chemistry-free statistical nodes.
    statistics_module = module_bytes(
        blobs["statistics_extracted.py"], "forward_statistics", HERE / "statistics_extracted.py"
    )
    result = {
        **meta,
        "splits": {
            s: statistics_module.aggregate_split(rows[s], outputs[s], targets[s], caps[s], plan["bootstrap"])
            for s in SPLITS
        },
    }
    io.same(result, expected, "complete recorded analysis", tolerance=1e-12)
    return result


def manifest_for(blobs):
    return {
        "schema_version": 1,
        "kind": "original-630-common-forward-statistical-replay",
        "scope": SCOPE,
        "files": {n: {"sha256": io.sha(v), "bytes": len(v)} for n, v in sorted(blobs.items())},
    }


def pack(inputs, analysis, audit, parent, parent_sha256, output):
    io.require(not io.no_symlinks(output).exists(), "Output already exists")
    mapping = {
        n: io.no_symlinks(inputs) / n
        for n in ("attempts.json", "normalized-occurrences.json", "integrity-validation.json")
    }
    mapping.update(
        {
            "recorded-analysis.json": io.no_symlinks(analysis) / "analysis.json",
            "recorded-validation.json": io.no_symlinks(analysis) / "validation.json",
            "audit-protocol.json": io.no_symlinks(audit) / "protocol.json",
            "audit-progress.json": io.no_symlinks(audit) / "progress.json",
            "parent-review.json": io.no_symlinks(parent),
        }
    )
    mapping.update({n: p for n, (p, _) in fixed_files().items()})
    raw = io.read_bytes(parent)
    io.same(io.sha(raw), io.sha_value(parent_sha256), "external parent review SHA")
    for entry in io.parse(raw)["assessment_evidence"]:
        name = io.safe_relative(entry["artifact"])
        mapping["parent-review-evidence/" + name] = Path(parent).parent / name
    blobs = io.snapshot(mapping)
    io.same(blobs["parent-review.json"], raw, "parent review changed before snapshot")
    validate(blobs, calculate=True)
    io.unchanged(mapping, blobs)
    output = io.new_directory(output)
    for name, value in blobs.items():
        path = output / io.safe_relative(name)
        path.parent.mkdir(parents=True, exist_ok=True)
        io.write_new(path, value)
    io.unchanged(mapping, blobs)
    manifest = io.encoded(manifest_for(blobs))
    io.write_new(output / "manifest.json", manifest)
    return {"status": "packed", "attempts": 630, "manifest_sha256": io.sha(manifest)}


def load_bundle(bundle, expected_manifest_sha256):
    bundle = io.no_symlinks(bundle)
    raw = io.read_bytes(bundle / "manifest.json")
    io.same(io.sha(raw), io.sha_value(expected_manifest_sha256), "external manifest SHA")
    manifest = io.parse(raw)
    io.require(set(manifest) == {"schema_version", "kind", "scope", "files"}, "Manifest fields mismatch")
    io.same(
        {k: v for k, v in manifest.items() if k != "files"},
        {k: v for k, v in manifest_for({}).items() if k != "files"},
        "manifest schema/scope",
    )
    mapping = {io.safe_relative(n): bundle / io.safe_relative(n) for n in manifest["files"]}
    io.require("manifest.json" not in mapping, "Manifest lists itself")
    observed = set()
    for base, dirs, files in os.walk(bundle, followlinks=False):
        for name in dirs + files:
            path = Path(base) / name
            io.no_symlinks(path)
            if name in files:
                observed.add(path.relative_to(bundle).as_posix())
    io.require(observed == set(mapping) | {"manifest.json"}, "Unexpected/missing bundle file")
    blobs = io.snapshot(mapping)
    io.same(manifest["files"], manifest_for(blobs)["files"], "bundle artifact hashes/lengths")
    return {**mapping, "manifest.json": bundle / "manifest.json"}, {**blobs, "manifest.json": raw}, blobs


def check(bundle, manifest_sha256, output, calculate):
    output = io.no_symlinks(output)
    io.require(
        not output.exists() and not output.is_relative_to(io.no_symlinks(bundle)),
        "Output exists or is inside bundle",
    )
    mapping, before, blobs = load_bundle(bundle, manifest_sha256)
    report = validate(blobs, calculate=calculate)
    io.unchanged(mapping, before)
    result = {
        "status": "passed",
        "attempts": 630,
        "scope": SCOPE,
        "comparison_recomputed": calculate,
        "manifest_sha256": manifest_sha256,
        "analysis_helper_sha256": HELPER_SHA,
        "extracted_statistics_sha256": EXTRACTED_SHA,
        "absolute_tolerance": 1e-12,
        "input_before_after_sha256": {n: io.sha(v) for n, v in sorted(before.items())},
    }
    if calculate:
        import numpy as np

        result["runtime"] = {"python": sys.version, "numpy": np.__version__}
        result["recomputed_analysis_sha256"] = io.sha(io.encoded(report))
        output = io.new_directory(output)
        io.write_new(output / "analysis.json", io.encoded(report))
        io.unchanged(mapping, before)
        io.write_new(output / "replay-validation.json", io.encoded(result))
    else:
        io.write_new(output, io.encoded(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    pack_parser = commands.add_parser("pack")
    for name in ("inputs", "analysis", "audit", "parent-review", "output"):
        pack_parser.add_argument("--" + name, type=Path, required=True)
    pack_parser.add_argument("--parent-review-sha256", required=True)
    for command in ("verify", "replay"):
        p = commands.add_parser(command)
        p.add_argument("--bundle", type=Path, required=True)
        p.add_argument("--expected-manifest-sha256", required=True)
        p.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "pack":
        result = pack(
            args.inputs, args.analysis, args.audit, args.parent_review, args.parent_review_sha256, args.output
        )
    else:
        result = check(args.bundle, args.expected_manifest_sha256, args.output, args.command == "replay")
    print(json.dumps({k: result[k] for k in ("status", "attempts", "manifest_sha256")}))


if __name__ == "__main__":
    main()
