"""Byte-preserving, model-free replay of the completed original 630-attempt study."""

import argparse
import ast
from collections import Counter
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
import statistics
import sys

HERE = Path(__file__).resolve().parent
HELPER = "analyze_generation_historical.py"
HELPER_SHA = "cb6992cda28895fa42791ad3969df5cbafd48a71b4ab8362ecc180cfeb4f2cc0"
ATTRIBUTION = "chembl-source-provenance.json"
ATTRIBUTION_SHA = "5d680d5dad8a273989d3548035d1e3352453e4939bd6fa22f041c7ad0da5f94b"
PROVENANCE = "frozen-campaign-code-provenance.json"
PROVENANCE_SHA = "8c8d9a883fe29acc8ffc1d6c3826a0056def13a93fa852508617069f3d48075e"
COMMIT = "d972124d587e647ef48455cf1a038e6be96e5b04"
SPLITS = {"scaffold": (360, 12), "time": (180, 12), "target": (90, 6)}
GROUPS = ("baselines", "gpu0", "gpu1")
SEEDS = [42, 43, 44]
METHODS = {
    "baselines": ["nearest_catalog", "nearest_3d", "reaction_enumeration"],
    "gpu0": ["pipeline", "pipeline_no_3d", "pipeline_no_refinement"],
    "gpu1": [
        "synllama_only",
        "pipeline_legacy_search",
        "pipeline_encounter_allocation",
        "pipeline_rule_only",
    ],
}
STATUSES = {
    "success",
    "no_candidates",
    "timeout",
    "memory_limit",
    "worker_error",
    "missing_response",
    "budget_exceeded",
    "invalid_response",
    "launch_error",
    "infrastructure_interruption",
}
SCOPE = (
    "Stored original 630 generation attempts and common 3D scores only. No generation, chemistry, "
    "3D rescoring, independent forward re-inference, or supplementary 270-study analysis. "
    "Recovery/similarity are proxies, not experimental activity, synthesis feasibility, or yield proof. "
    "Repeated targets across splits are not independent cohorts. Hash-bound review records describe "
    "their own checks; hashes do not establish independent review or reproduce packages/models. "
    "Original paired contrasts pair targets: method-specific success-conditional populations may differ; "
    "no like-for-like conditional-score improvement inference is supported."
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def sha_value(value):
    require(isinstance(value, str) and re.fullmatch("[0-9a-f]{64}", value), "Invalid SHA256")
    return value


def safe_relative(value):
    require(isinstance(value, str) and value and "\\" not in value, "Unsafe artifact path")
    p = PurePosixPath(value)
    require(
        not p.is_absolute() and all(x not in {".", "..", ""} for x in value.split("/")),
        "Unsafe artifact path",
    )
    require("\x00" not in value and ":" not in value, "Unsafe artifact path")
    return value


def no_symlinks(path):
    path = Path(os.path.abspath(path))
    for part in (path, *path.parents):
        require(not part.is_symlink(), f"Symlink forbidden: {part}")
    return path


def read_bytes(path):
    path = no_symlinks(path)
    require(stat.S_ISREG(path.stat().st_mode), "Input must be a regular file")
    with path.open("rb") as stream:
        require(stat.S_ISREG(os.fstat(stream.fileno()).st_mode), "Input must be a regular file")
        return stream.read()


def parse(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def number(value):
        result = float(value)
        require(math.isfinite(result), "Nonfinite JSON number")
        return result

    def constant(value):
        raise ValueError(f"Nonfinite JSON constant: {value}")

    return json.loads(raw, object_pairs_hook=pairs, parse_float=number, parse_constant=constant)


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def same(actual, expected, at="root", tolerance=0.0):
    """JSON type/shape equality; only float values use an absolute tolerance."""
    require(type(actual) is type(expected), f"Type mismatch at {at}")
    if isinstance(expected, dict):
        require(actual.keys() == expected.keys(), f"Field mismatch at {at}")
        for name in expected:
            same(actual[name], expected[name], f"{at}.{name}", tolerance)
    elif isinstance(expected, list):
        require(len(actual) == len(expected), f"Length mismatch at {at}")
        for i, (a, e) in enumerate(zip(actual, expected)):
            same(a, e, f"{at}[{i}]", tolerance)
    elif isinstance(expected, float):
        require(
            math.isfinite(actual) and math.isfinite(expected) and abs(actual - expected) <= tolerance,
            f"Numeric mismatch at {at}",
        )
    else:
        require(actual == expected, f"Value mismatch at {at}")


def paths(split):
    result = {"panel": f"panels/direct-human-{split}/panel.json"}
    for group in GROUPS:
        result[f"protocol_{group}"] = (
            f"main-protocol-{group}.json"
            if split == "scaffold"
            else f"secondary-{split}-protocol-{group}.json"
        )
        folder = f"comparison-{group}" if split == "scaffold" else f"comparison-{split}-{group}"
        result[f"runs_{group}"] = f"{folder}/runs.json"
        result[f"manifest_{group}"] = f"{folder}/manifest.json"
    folder = "common-generated-scores-8c" if split == "scaffold" else f"common-generated-{split}-scores-8c"
    for role in ("protocol", "progress", "generated_scores"):
        result[f"common_{role}"] = f"{folder}/{role}.json"
    result["comparison"] = f"analysis-generation-{split}/comparison.json"
    return result


def input_names():
    return {PROVENANCE, ATTRIBUTION, *(p for split in SPLITS for p in paths(split).values())}


def historical(raw):
    require(sha(raw) == HELPER_SHA, "Historical helper hash mismatch")
    # Execute only the reviewed module, never its filesystem-dependent main().
    spec = importlib.util.spec_from_loader("generation_historical", loader=None)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(HERE / HELPER)
    exec(compile(raw, HELPER, "exec"), module.__dict__)
    return module


def report_text(raw):
    """Read unchanged explanatory literals from the pinned helper's report dictionary."""
    result = {}
    for node in ast.walk(ast.parse(raw)):
        if isinstance(node, ast.Dict):
            for key, value in zip(node.keys, node.values):
                if isinstance(key, ast.Constant) and key.value in {
                    "uncertainty",
                    "scope",
                    "assessment_design_note",
                }:
                    result[key.value] = ast.literal_eval(value)
    require(len(result) == 3, "Historical report literals missing")
    return result


def normalized_protocol(raw, recorded):
    recorded = copy.deepcopy(recorded)
    require(len(recorded["methods"]) == len(raw["methods"]), "Normalized method count mismatch")
    for before, after in zip(raw["methods"], recorded["methods"]):
        # Only loader-added fields absent in the actual raw protocol can be projected out.
        for name in ("resolved_assets", "overlap_audit", "search_space_overlap_audit"):
            if name not in before:
                after.pop(name, None)
    same(recorded, raw, "recorded protocol")


def validate_split(split, blobs, code):
    names = paths(split)
    data = {role: parse(blobs[name]) for role, name in names.items()}
    panel = data["panel"]
    require(panel["split"]["kind"] == split, "Panel split mismatch")
    records = {r["id"]: r for r in panel["compounds"]}
    require(len(records) == len(panel["compounds"]), "Duplicate compound ID")
    queries = {q["id"]: records[q["compound_id"]] for q in panel["queries"]}
    require(len(queries) == len(panel["queries"]) == SPLITS[split][1], "Query count/ID mismatch")
    require(not any(q.get("reference_sdf") for q in panel["queries"]), "Expected ligand-only queries")
    targets = {q: r["target_id"] for q, r in queries.items()}
    require(all(isinstance(t, str) and t for t in targets.values()), "Invalid target")
    source_inventory = {s["id"]: s["sha256"] for s in panel["sources"]}
    require(len(source_inventory) == len(panel["sources"]), "Duplicate panel source ID")
    attribution = parse(blobs[ATTRIBUTION])
    sources = {s["name"]: s for s in attribution["sources"]}
    for s in panel["sources"]:
        same(s["sha256"], sources[s["id"]]["sha256"], "ChEMBL source hash")
        same(s["url"], sources[s["id"]]["url"], "ChEMBL source URL")
    runs, caps, run_hashes = [], {}, {}
    expected = set()
    for group in GROUPS:
        protocol = data[f"protocol_{group}"]
        manifest = data[f"manifest_{group}"]
        same(protocol["seeds"], SEEDS, "generation seeds")
        mids = (
            METHODS[group]
            if split == "scaffold"
            else METHODS[group][:1]
            if group != "baselines"
            else METHODS[group]
        )
        same([m["id"] for m in protocol["methods"]], mids, "method grid")
        budget = protocol["budget"]
        same(
            budget,
            {"max_candidates": 3, "max_evaluations": 1000000, "max_memory_mb": 24000, "wall_seconds": 360.0},
            "common budget",
        )
        same(manifest["protocol_sha256"], sha(blobs[names[f"protocol_{group}"]]), "protocol SHA")
        same(manifest["dataset"]["dataset_sha256"], sha(blobs[names["panel"]]), "dataset SHA")
        same(manifest["dataset"]["source_artifacts"], source_inventory, "dataset source inventory")
        same(manifest["code_sha256"], code, "original source39")
        normalized_protocol(protocol, manifest["protocol"])
        require(
            isinstance(manifest["versions"], dict)
            and all(
                isinstance(manifest["versions"].get(k), str) and manifest["versions"][k]
                for k in ("numpy", "rdkit", "scipy", "psutil")
            ),
            "Missing original package versions",
        )
        require(
            all(isinstance(manifest.get(k), str) and manifest[k] for k in ("python", "platform")),
            "Missing original runtime provenance",
        )
        if "dataset_sha256" in protocol.get("campaign", {}):
            same(protocol["campaign"]["dataset_sha256"], sha(blobs[names["panel"]]), "inline dataset SHA")
        group_grid = {(m, q, s) for m in mids for q in queries for s in SEEDS}
        group_runs = data[f"runs_{group}"]
        actual = [(r["method_id"], r["query_id"], r["seed"]) for r in group_runs]
        require(
            len(actual) == len(set(actual)) and set(actual) == group_grid, "Incomplete/duplicate run grid"
        )
        expected.update(group_grid)
        for method in protocol["methods"]:
            override = method.get("budget_overrides", {})
            require(not set(override) - set(budget), "Unknown budget override")
            require(
                all(type(v) in (int, float) and 0 < v <= budget[k] for k, v in override.items()),
                "Invalid method cap",
            )
            cap = override.get("max_candidates", budget["max_candidates"])
            require(type(cap) is int and cap == 3, "Original output cap must be 3")
            caps[method["id"]] = cap
        for run in group_runs:
            require(type(run["seed"]) is int and run["seed"] in SEEDS, "Invalid generation seed")
            same(run["target_id"], targets[run["query_id"]], "run target")
            require(run["status"] in STATUSES, "Unknown/nonterminal status")
            candidates = run["evaluation"]["candidates"]
            require(isinstance(candidates, list) and len(candidates) <= 3, "Invalid accepted output count")
            if run["status"] != "success":
                require(not candidates, "Failed attempt has accepted candidates")
            metrics = run["evaluation"]["metrics"]
            require(
                type(metrics["valid_candidate_count"]) is int
                and metrics["valid_candidate_count"] == len(candidates),
                "Valid count mismatch",
            )
            wall = run["resources"]["wall_seconds"]
            if wall is None:
                require(
                    run["status"] == "infrastructure_interruption"
                    and run["resources"].get("measurement_status") == "lost_with_runner",
                    "Unexplained missing runtime",
                )
            else:
                require(type(wall) in (int, float) and wall >= 0, "Invalid wall seconds")
            for c in candidates:
                require(
                    all(isinstance(c.get(k), str) and c[k] for k in ("smiles", "family")),
                    "Missing candidate identity",
                )
        runs.extend(group_runs)
        run_hashes[group] = sha(blobs[names[f"runs_{group}"]])
    require(len(runs) == len(expected) == SPLITS[split][0], "Split attempt total mismatch")
    common, progress, scored = (data[f"common_{r}"] for r in ("protocol", "progress", "generated_scores"))
    common_sha = sha(blobs[names["common_protocol"]])
    same(common["dataset_sha256"], sha(blobs[names["panel"]]), "common dataset SHA")
    same(common["seeds"], SEEDS, "common seeds")
    same(common["expected_runs"], SPLITS[split][0], "common planned count")
    same(
        common["code_sha256"],
        {p: code[p] for p in ("benchmarks/research/controls.py", "chembounce_synthesis/conformers.py")},
        "common scorer source",
    )
    same(
        common["conformers"],
        {
            "energy_window": 5.0,
            "max_iterations": 500,
            "num_conformers": 8,
            "pharmacophore_sigma": 1.0,
            "reject_unspecified_stereo": True,
            "strain_penalty": 0.05,
        },
        "original common conformer settings",
    )
    lookup, cached = {}, {}
    for row in scored:
        key = (row["method_id"], row["query_id"], row["generation_seed"], row["rank"], row["smiles"])
        require(
            key not in lookup and type(row["rank"]) is int and type(row["generation_seed"]) is int,
            "Duplicate/invalid scored occurrence",
        )
        same(row["target_id"], targets[row["query_id"]], "scored target")
        require(type(row["all_3d_seeds_succeeded"]) is bool, "Invalid 3D completion flag")
        trials = row["trials"]
        if trials:
            same([t["seed"] for t in trials], SEEDS, "3D trial seeds")
            require(all(t["status"] in {"success", "failed"} for t in trials), "Unknown 3D trial status")
        else:
            require(isinstance(row.get("error"), str) and row["error"], "Unexplained missing trials")
        same(
            row["all_3d_seeds_succeeded"],
            bool(trials) and all(t["status"] == "success" for t in trials),
            "3D completion disagrees with trials",
        )
        for metric in ("hit_shape", "hit_pharmacophore", "score"):
            if row["all_3d_seeds_succeeded"]:
                values = [t["metrics"][metric] for t in trials]
                require(
                    all(type(v) in (int, float) and math.isfinite(v) for v in values), "Invalid trial score"
                )
                same(row["scores"][metric], statistics.median(values), "3D trial median", tolerance=1e-12)
                # Original evaluation.finite_scores coerces stability samples to float.
                values = [float(v) for v in values]
                stability = {
                    "n": len(values),
                    "mean": statistics.mean(values),
                    "median": statistics.median(values),
                    "min": min(values),
                    "max": max(values),
                    "range": max(values) - min(values),
                    "sample_sd": statistics.stdev(values),
                }
                same(row["scores"][metric + "_stability"], stability, "3D trial stability", tolerance=1e-12)
            else:
                require(
                    metric not in row["scores"] and metric + "_stability" not in row["scores"],
                    "Failed 3D seeds must not have 3D summary fields",
                )
        value = {k: row[k] for k in ("scores", "trials", "error", "all_3d_seeds_succeeded")}
        if row["cache_key"] in cached:
            same(value, cached[row["cache_key"]], "reused common score cache")
        cached[row["cache_key"]] = value
        lookup[key] = row
    cache_keys, used = set(), set()
    for run in runs:
        for ordinal, candidate in enumerate(run["evaluation"]["candidates"], 1):
            key = (run["method_id"], run["query_id"], run["seed"], ordinal, candidate["smiles"])
            require(key in lookup, "Missing valid-output ordinal score join")
            used.add(key)
            cache = sha(
                json.dumps(
                    [common_sha, queries[run["query_id"]]["smiles"], candidate["smiles"]],
                    separators=(",", ":"),
                ).encode()
            )
            same(lookup[key]["cache_key"], cache, "common cache identity")
            cache_keys.add(cache)
    require(used == lookup.keys(), "Extraneous scored occurrence")
    same(
        progress,
        {
            "protocol_sha256": common_sha,
            "input_runs_sha256": run_hashes,
            "runs_observed": len(runs),
            "expected_runs": SPLITS[split][0],
            "candidate_occurrences": len(used),
            "unique_query_candidate_pairs": len(cache_keys),
            "assessed_occurrences": len(scored),
            "failed_3d_occurrences": sum(not r["all_3d_seeds_succeeded"] for r in scored),
            "queued_unique_pairs": 0,
            "complete": True,
        },
        "common completion",
    )
    sources = {names[f"{r}_{g}"]: sha(blobs[names[f"{r}_{g}"]]) for g in GROUPS for r in ("protocol", "runs")}
    sources[names["common_generated_scores"]] = sha(blobs[names["common_generated_scores"]])
    audit = {
        "split": split,
        "expected_runs": len(expected),
        "observed_runs": len(runs),
        "missing": [],
        "unexpected": [],
        "duplicates": [],
        "status_counts": dict(Counter(r["status"] for r in runs)),
        "complete": True,
    }
    return (
        runs,
        scored,
        targets,
        caps,
        {
            "split": split,
            "grid_audit": audit,
            "source_sha256": sources,
            "dataset_sha256": sha(blobs[names["panel"]]),
        },
    )


def review_check(review, blobs):
    same(review["schema_version"], 1, "review schema")
    same(review["status"], "passed", "parent status")
    same(review["total_runs"], 630, "parent total")
    same(review["split_counts"], {s: n for s, (n, _) in SPLITS.items()}, "parent split counts")
    same(review["historical_helper_sha256"], HELPER_SHA, "parent helper SHA")
    same(review["original_source_commit"], COMMIT, "parent source commit")
    same(review["input_sha256"], {n: sha(blobs[n]) for n in input_names()}, "parent input snapshot")
    source = parse(blobs[PROVENANCE])
    same(review["source39_sha256"], source["code_sha256"], "parent original source map")
    evidence = review["assessment_evidence"]
    require(isinstance(evidence, list) and evidence, "At least one actual parent audit is required")
    seen = set()
    for item in evidence:
        require(set(item) == {"artifact", "sha256"}, "Invalid parent evidence entry")
        name = safe_relative(item["artifact"])
        require(name not in seen, "Duplicate parent evidence")
        seen.add(name)
        raw = blobs["parent-review-evidence/" + name]
        same(sha(raw), sha_value(item["sha256"]), "parent evidence SHA")
        require(isinstance(parse(raw), dict), "Parent audit must be a JSON object")
    return seen


def validate(blobs, calculate=False):
    same(sha(blobs[HELPER]), HELPER_SHA, "historical helper SHA")
    same(sha(blobs[ATTRIBUTION]), ATTRIBUTION_SHA, "historical ChEMBL attribution SHA")
    same(sha(blobs[PROVENANCE]), PROVENANCE_SHA, "frozen campaign provenance SHA")
    source = parse(blobs[PROVENANCE])
    same(source["implementation_commit"], COMMIT, "frozen commit")
    require(len(source["code_sha256"]) == source["file_count"] == 39, "Original 39-source inventory required")
    evidence = review_check(parse(blobs["parent-review.json"]), blobs)
    require(
        set(blobs)
        == input_names() | {HELPER, "parent-review.json"} | {"parent-review-evidence/" + p for p in evidence},
        "Bundle file membership mismatch",
    )
    prepared = {s: validate_split(s, blobs, source["code_sha256"]) for s in SPLITS}
    literals = report_text(blobs[HELPER])
    analysis_keys = {
        "methods",
        "paired_pipeline_minus_comparator",
        "infrastructure_missing_outcome_sensitivity",
        "complete_case_secondary",
    }
    for split, (_, _, _, _, metadata) in prepared.items():
        report = parse(blobs[paths(split)["comparison"]])
        require(set(report) == set(metadata) | set(literals) | analysis_keys, "Comparison fields mismatch")
        require(all(isinstance(report[k], dict) for k in analysis_keys), "Invalid comparison section")
        for name, value in literals.items():
            same(report[name], value, f"comparison.{split}.{name}")
        for name, value in metadata.items():
            same(report[name], value, f"comparison.{split}.{name}")
    # All three grids and common scoring passes must validate before any analysis begins.
    reports = {}
    if calculate:
        module = historical(blobs[HELPER])
        for split, (runs, scored, targets, caps, metadata) in prepared.items():
            actual = {**metadata, **module.analyze_runs(runs, scored, targets, caps), **literals}
            expected = parse(blobs[paths(split)["comparison"]])
            same(actual, expected, f"comparison.{split}", tolerance=1e-12)
            reports[split] = actual
    return reports


def snapshot(paths_by_name):
    return {name: read_bytes(path) for name, path in paths_by_name.items()}


def unchanged(paths_by_name, before):
    for name, path in paths_by_name.items():
        same(sha(read_bytes(path)), sha(before[name]), f"input changed during operation: {name}")


def write_new(path, raw):
    path = no_symlinks(path)
    with path.open("xb") as stream:
        stream.write(raw)


def new_directory(path):
    path = no_symlinks(path)
    path.mkdir()  # Parent must exist; existing directories and dangling symlinks fail.
    return path


def pack(source, parent_review, parent_sha256, output):
    source = no_symlinks(source)
    parent_review = no_symlinks(parent_review)
    require(not no_symlinks(output).exists(), "Output already exists")
    mapping = {n: source / n for n in input_names()}
    mapping[ATTRIBUTION] = HERE / ATTRIBUTION
    mapping[HELPER] = HERE / HELPER
    mapping["parent-review.json"] = parent_review
    raw = read_bytes(parent_review)
    same(sha(raw), sha_value(parent_sha256), "external parent review SHA")
    for evidence in parse(raw)["assessment_evidence"]:
        name = safe_relative(evidence["artifact"])
        mapping["parent-review-evidence/" + name] = parent_review.parent / name
    blobs = snapshot(mapping)
    same(blobs["parent-review.json"], raw, "parent review changed before snapshot")
    validate(blobs, calculate=True)
    unchanged(mapping, blobs)
    output = new_directory(output)
    for name, content in blobs.items():
        path = output / safe_relative(name)
        path.parent.mkdir(parents=True, exist_ok=True)
        write_new(path, content)
    unchanged(mapping, blobs)
    manifest = {
        "schema_version": 1,
        "kind": "original-630-statistical-replay",
        "scope": SCOPE,
        "historical_helper_sha256": HELPER_SHA,
        "original_source_commit": COMMIT,
        "files": {n: {"sha256": sha(v), "bytes": len(v)} for n, v in sorted(blobs.items())},
    }
    # Write manifest last: a failed/incomplete pack never has a complete manifest.
    write_new(output / "manifest.json", encoded(manifest))
    return {
        "manifest_sha256": sha(encoded(manifest)),
        "files": len(blobs),
        "attempts": 630,
        "input_before_after_sha256": {n: sha(v) for n, v in sorted(blobs.items())},
    }


def load_bundle(bundle, expected_manifest_sha256):
    bundle = no_symlinks(bundle)
    raw = read_bytes(bundle / "manifest.json")
    same(sha(raw), sha_value(expected_manifest_sha256), "externally retained manifest SHA")
    manifest = parse(raw)
    require(
        set(manifest)
        == {"schema_version", "kind", "scope", "historical_helper_sha256", "original_source_commit", "files"},
        "Manifest fields mismatch",
    )
    for key, value in {
        "schema_version": 1,
        "kind": "original-630-statistical-replay",
        "scope": SCOPE,
        "historical_helper_sha256": HELPER_SHA,
        "original_source_commit": COMMIT,
    }.items():
        same(manifest[key], value, "manifest." + key)
    mapping = {safe_relative(n): bundle / safe_relative(n) for n in manifest["files"]}
    require("manifest.json" not in mapping, "Manifest cannot list itself")
    observed = set()
    for base, dirs, files in os.walk(bundle, followlinks=False):
        for name in dirs + files:
            path = Path(base) / name
            no_symlinks(path)
            if name in files:
                observed.add(path.relative_to(bundle).as_posix())
    require(observed == set(mapping) | {"manifest.json"}, "Unexpected/missing bundle file")
    blobs = snapshot(mapping)
    for name, content in blobs.items():
        same(manifest["files"][name], {"sha256": sha(content), "bytes": len(content)}, "artifact " + name)
    mapping["manifest.json"] = bundle / "manifest.json"
    before = {**blobs, "manifest.json": raw}
    return mapping, before, blobs


def verify(bundle, expected_manifest_sha256):
    mapping, before, blobs = load_bundle(bundle, expected_manifest_sha256)
    validate(blobs)
    unchanged(mapping, before)
    return {
        "status": "verified",
        "attempts": 630,
        "comparison_recomputed": False,
        "manifest_sha256": expected_manifest_sha256,
        "scope": SCOPE,
        "input_before_after_sha256": {n: sha(v) for n, v in sorted(before.items())},
    }


def replay(bundle, expected_manifest_sha256, output):
    require(not no_symlinks(output).is_relative_to(no_symlinks(bundle)), "Output cannot be inside the bundle")
    mapping, before, blobs = load_bundle(bundle, expected_manifest_sha256)
    # Reject existing output before analysis; mkdir remains exclusive after validation.
    output = no_symlinks(output)
    require(not output.exists(), "Output already exists")
    reports = validate(blobs, calculate=True)
    unchanged(mapping, before)
    output = new_directory(output)
    for split, report in reports.items():
        write_new(output / f"{split}-comparison.json", encoded(report))
    unchanged(mapping, before)
    import numpy as np

    result = {
        "status": "passed",
        "attempts": 630,
        "absolute_tolerance": 1e-12,
        "manifest_sha256": expected_manifest_sha256,
        "historical_helper_sha256": HELPER_SHA,
        "runtime": {"python": sys.version, "numpy": np.__version__},
        "scope": SCOPE,
        "input_before_after_sha256": {n: sha(v) for n, v in sorted(before.items())},
        "recomputed_sha256": {s: sha(encoded(r)) for s, r in reports.items()},
    }
    write_new(output / "validation.json", encoded(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("pack")
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--parent-review", type=Path, required=True)
    p.add_argument("--parent-review-sha256", required=True)
    p.add_argument("--output", type=Path, required=True)
    for cmd in ("verify", "replay"):
        p = commands.add_parser(cmd)
        p.add_argument("--bundle", type=Path, required=True)
        p.add_argument("--expected-manifest-sha256", required=True)
        p.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "pack":
        result = pack(args.source, args.parent_review, args.parent_review_sha256, args.output)
    elif args.command == "verify":
        require(
            not no_symlinks(args.output).is_relative_to(no_symlinks(args.bundle)),
            "Output cannot be inside the bundle",
        )
        require(not no_symlinks(args.output).exists(), "Output already exists")
        result = verify(args.bundle, args.expected_manifest_sha256)
        write_new(args.output, encoded(result))
    else:
        result = replay(args.bundle, args.expected_manifest_sha256, args.output)
    print(json.dumps({k: result[k] for k in ("status", "attempts", "manifest_sha256") if k in result}))


if __name__ == "__main__":
    main()
