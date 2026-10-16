"""Portable replay of reviewed supplementary-270 stored statistics; no chemistry."""

import argparse
import ast
from collections import Counter
import hashlib
import math
import os
from pathlib import Path
import statistics
import sys
from types import ModuleType, SimpleNamespace

HERE = Path(__file__).resolve().parent
PINS_SHA = "ec13ca5080cc0d814a04736105821f48bd8262f22f1f60f8a567cf75108f032c"
TRANSPORT_SHA = "9751209528300cdc38250d68c3a6a8748ed873362af1d1dabfa14ecbb2234113"
HELPER_SHA = "b9e24075c3940121c0dea6f0a6323a826824f372ca00129fe576194494f4a295"
PREPARED_SHA = "db76ab6355526ce833c06b443ccf9bd9b1222cbb26e253441f2f23ee2d4d9004"
GENERATION_SHA = "7d01e374c6e4e2a983f0c322b3036afc9258e3e0c7d7da334bf0403ab8f4d2f2"
GENERATION_GATE_SHA = "b9d1bc72af859c40cc2b8df67b3a156f71b21ebb4c3a9309eb9dc7b122f4c776"
SOURCE_COMMIT = "2328588d0969369e694d9293a1dccc89392a8203"
SPLITS = {"scaffold": 108, "time": 108, "target": 54}
METHODS = ["pipeline_current", "crem_train_partial_exo", "mol2mol_pretrained"]
SEEDS = [42, 43, 44]
FIELDS = ("split", "method_id", "query_id", "target_id", "seed")
ANALYSIS_FILES = {
    "comparison.json",
    "attempts.json",
    "recomputed-attempts.json",
    "unique-identities.json",
    "raw-accounting-summary.json",
    "scaffold-retention-summary.json",
    "scaffold-retention-occurrences.json",
    "analysis-plan.json",
}
NODES = {
    "assess_historical.py": {
        "import:collections",
        "import:statistics",
        "INFRA",
        "RETENTION_FIELDS",
        "BOUNDED",
        "key",
        "occurrence_key",
        "mean",
        "attempt_metrics",
        "hierarchical",
        "comparisons",
    },
    "evaluation_historical.py": {
        "import:math",
        "import:numpy",
        "finite_scores",
        "score_stability",
        "cluster_interval",
    },
}
SCOPE = (
    "Stored supplementary 270 attempts only, separate from original 630. CPU statistics; no chemistry, "
    "generation, identity normalization, route replay, models, or 3D rescoring. Chemistry-dependent records "
    "and validation are preserved and parent-reviewed, not independently re-established here. Delivery and "
    "structural proxies do not establish activity, synthesis probability, yield, or universal scaffold hopping. "
    "No-route is unassessed coverage. Declared forward evidence shares selection-model bias. "
    "Few-target, repeated-input and cross-split dependence limit descriptive unadjusted intervals."
)


def module_bytes(raw, name):
    module = ModuleType(name)
    module.__file__ = str(HERE / (name + ".py"))
    exec(compile(raw, name, "exec"), module.__dict__)
    return module


transport_path = HERE.parent / "generation_replay/replay.py"
transport_raw = transport_path.read_bytes()
if hashlib.sha256(transport_raw).hexdigest() != TRANSPORT_SHA:
    raise ValueError("Reviewed transport implementation changed")
io = module_bytes(transport_raw, "supplement_transport")
# Transport helpers use HERE only for their unrelated historical replay entrypoints,
# which are never called by this module.


def content_key(value):
    import json

    return io.sha(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode())


def fixed_files():
    raw = io.read_bytes(HERE / "pins.json")
    io.same(io.sha(raw), PINS_SHA, "fixed pins")
    return {
        **{n: (HERE / io.safe_relative(n), h) for n, h in io.parse(raw).items()},
        "pins.json": (HERE / "pins.json", PINS_SHA),
        "transport_helpers.py": (transport_path, TRANSPORT_SHA),
    }


def extraction(blobs):
    parts, records = [], []
    for filename, allowed in NODES.items():
        raw = blobs[filename]
        lines, found = raw.splitlines(keepends=True), set()
        for node in ast.parse(raw).body:
            label = None
            if isinstance(node, ast.ImportFrom):
                label = "import:" + str(node.module)
            elif isinstance(node, ast.Import) and len(node.names) == 1:
                label = "import:" + node.names[0].name
            elif isinstance(node, ast.FunctionDef):
                label = node.name
            elif (
                isinstance(node, ast.Assign)
                and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)
            ):
                label = node.targets[0].id
            if label in allowed:
                io.require(label not in found, "Duplicate selected AST node")
                found.add(label)
                part = b"".join(lines[node.lineno - 1 : node.end_lineno])
                parts.append(part)
                records.append(
                    {
                        "source": filename,
                        "node": label,
                        "first_line": node.lineno,
                        "last_line": node.end_lineno,
                        "sha256": io.sha(part),
                    }
                )
        io.same(found, allowed, "allowlisted extraction nodes")
    return b"\n".join(parts), records


def check_fixed(blobs):
    for name, (_, pin) in fixed_files().items():
        io.same(io.sha(blobs[name]), pin, "fixed artifact " + name)
    extracted, boundaries = extraction(blobs)
    io.same(extracted, blobs["statistics_extracted.py"], "verbatim AST extraction")
    io.same(boundaries, io.parse(blobs["extraction.json"])["selected_nodes"], "AST boundaries")
    io.same(io.sha(blobs["assess_historical.py"]), HELPER_SHA, "historical helper")
    for folder, anchor in (("prepared", PREPARED_SHA), ("generation", GENERATION_SHA)):
        io.same(io.sha(blobs[folder + "/manifest.json"]), anchor, folder + " anchor")
        manifest = io.parse(blobs[folder + "/manifest.json"])
        for name, raw in blobs.items():
            if name.startswith(folder + "/") and not name.endswith("/manifest.json"):
                io.same(io.sha(raw), manifest["files_sha256"][name.split("/", 1)[1]], "frozen metadata hash")
    plan = io.parse(blobs["prepared/analysis-plan.json"])
    grid = io.parse(blobs["prepared/grid.json"])
    common = io.parse(blobs["prepared/common-3d-plan.json"])
    execution = io.parse(blobs["prepared/execution-plan.json"])
    io.same(plan["expected_grid"], grid, "analysis grid")
    io.same(plan["expected_attempts"], 270, "analysis total")
    io.same(plan["split_attempts"], SPLITS, "analysis splits")
    io.same(common["source"], execution["source"], "current source")
    io.same(common["source"]["commit"], SOURCE_COMMIT, "current source commit")
    io.same(len(common["source"]["files_sha256"]), 41, "current source file count")
    io.same(common["seeds"], SEEDS, "fixed scoring seeds")
    io.same(
        common["source"]["files_sha256"]["benchmarks/research/evaluation.py"],
        io.sha(blobs["evaluation_historical.py"]),
        "interval source",
    )
    stats = module_bytes(extracted, "supplement_statistics_only")
    stats.preparation = SimpleNamespace(SPLITS=tuple(SPLITS), METHODS=METHODS)
    return stats, grid, common, execution


def reviewed(blobs):
    review = io.parse(blobs["parent-review.json"])
    dynamic = {
        n: io.sha(v)
        for n, v in blobs.items()
        if n.startswith(("analysis/", "common/"))
        or n in {"normalized-occurrences.json", "integrity-validation.json"}
    }
    expected = {
        "schema_version": 1,
        "status": "passed",
        "complete": True,
        "stage": "final_supplementary_270_completion",
        "total_attempts": 270,
        "split_attempts": SPLITS,
        "analysis_helper_sha256": HELPER_SHA,
        "prepared_manifest_sha256": PREPARED_SHA,
        "generation_manifest_sha256": GENERATION_SHA,
        "generation_parent_gate_sha256": GENERATION_GATE_SHA,
        "source_commit": SOURCE_COMMIT,
        "input_sha256": dynamic,
    }
    for name, value in expected.items():
        io.same(review[name], value, "final parent review " + name)
    evidence, names = {}, set()
    for item in review["assessment_evidence"]:
        io.require(set(item) == {"role", "artifact", "sha256"}, "Parent evidence schema")
        name = io.safe_relative(item["artifact"])
        io.require(name not in names and item["role"] not in evidence, "Duplicate parent evidence")
        names.add(name)
        raw = blobs["parent-review-evidence/" + name]
        io.same(io.sha(raw), io.sha_value(item["sha256"]), "parent evidence SHA")
        io.require(isinstance(io.parse(raw), dict), "Parent evidence must be JSON object")
        evidence[item["role"]] = (item["sha256"], io.parse(raw))
    io.require(
        {"generation_parent_gate", "final_integrity_and_statistics"} <= set(evidence),
        "Required actual reviews missing",
    )
    io.same(evidence["generation_parent_gate"][0], GENERATION_GATE_SHA, "final generation gate")
    io.same(evidence["final_integrity_and_statistics"][1]["passed"], True, "final independent assessment")
    io.same(
        set(blobs),
        set(fixed_files())
        | set(dynamic)
        | {"parent-review.json"}
        | {"parent-review-evidence/" + n for n in names},
        "exact bundle membership",
    )
    return dynamic


def archive(blobs, prefix, expected_stage):
    manifest = io.parse(blobs[prefix + "/manifest.json"])
    io.same(manifest["stage"], expected_stage, "recorded archive stage")
    io.same(manifest["assessment_code_sha256"], HELPER_SHA, "recorded assessment source")
    io.same(manifest["study_manifest_sha256"], PREPARED_SHA, "recorded study")
    io.same(manifest["generation_manifest_sha256"], GENERATION_SHA, "recorded generation snapshot")
    contents = {
        n[len(prefix) + 1 :]: io.sha(v)
        for n, v in blobs.items()
        if n.startswith(prefix + "/") and n != prefix + "/manifest.json"
    }
    io.same(manifest["files_sha256"], contents, "recorded archive contents")
    return manifest


def finite_metric(value, name):
    io.require(
        value is None or (type(value) in {int, float} and math.isfinite(value)),
        "Invalid numeric metric " + name,
    )


def validate_records(blobs, stats, grid, common, execution):
    am = archive(blobs, "analysis", "final_comparison")
    cm = archive(blobs, "common", "common_3d")
    io.same(am["complete_requires_parent_review"], True, "analysis completion gate")
    io.same(am["expected_attempts"], 270, "analysis attempt count")
    io.same(cm["complete"], True, "complete common scoring required")
    io.same(
        am["common_scores_manifest_sha256"], io.sha(blobs["common/manifest.json"]), "common archive anchor"
    )
    io.same(cm["common_plan_sha256"], io.sha(blobs["prepared/common-3d-plan.json"]), "common plan hash")
    io.same(cm["source"], common["source"], "scoring source")
    python = execution["commands"][0]["argv"][0]
    io.same(cm["runtime_python"], python, "historical runtime path")
    io.same(cm["runtime_packages"], execution["runtimes"][python]["packages"], "historical runtime inventory")
    io.same(cm["execution"], "CPU only, sequential, generation budgets excluded", "fixed sequential scoring")
    io.same(set(am["files_sha256"]), ANALYSIS_FILES, "full original analysis membership")
    io.same(
        blobs["analysis/analysis-plan.json"],
        blobs["prepared/analysis-plan.json"],
        "original analysis plan bytes",
    )
    io.same(
        blobs["analysis/recomputed-attempts.json"],
        blobs["common/recomputed-attempts.json"],
        "stored chemical reevaluation bytes",
    )
    attempts = io.parse(blobs["analysis/recomputed-attempts.json"])
    metrics = io.parse(blobs["analysis/attempts.json"])
    occurrences = io.parse(blobs["normalized-occurrences.json"])
    joined = io.parse(blobs["common/occurrences.json"])
    retained = io.parse(blobs["analysis/scaffold-retention-occurrences.json"])
    io.same(len(grid), 270, "fixed grid count")
    io.same(len({stats.key(g) for g in grid}), 270, "unique fixed grid")
    io.same([stats.key(a) for a in attempts], [stats.key(g) for g in grid], "exact attempt key order")
    io.same([stats.key(a) for a in metrics], [stats.key(g) for g in grid], "exact metric key order")
    io.same(dict(Counter(g["split"] for g in grid)), SPLITS, "split counts")
    expected_occurrences = []
    metric_names = set(io.parse(blobs["prepared/analysis-plan.json"])["evaluator_metrics"])
    for a, m, g in zip(attempts, metrics, grid):
        for field in FIELDS:
            io.same(a[field], g[field], "attempt " + field)
            io.same(m[field], g[field], "metric " + field)
        io.require(type(a["seed"]) is int and a["seed"] in SEEDS, "Invalid attempt seed")
        io.require(
            a["method_id"] in METHODS and a["status"] in io.STATUSES | {"interrupted"},
            "Invalid method/status",
        )
        io.same(a["status"], m["status"], "stored status")
        io.same(set(a["evaluation"]["metrics"]), metric_names, "original evaluator metric schema")
        candidates = a["evaluation"]["candidates"]
        io.require(isinstance(candidates, list) and len(candidates) <= 3, "Candidate occurrence cap")
        io.require(a["status"] == "success" or candidates == [], "Failed attempt promoted to delivery")
        ranks = [c["rank"] for c in candidates]
        io.require(
            all(type(r) is int and 1 <= r <= 3 for r in ranks) and ranks == sorted(set(ranks)),
            "Invalid raw ranks",
        )
        for name, value in a["evaluation"]["metrics"].items():
            finite_metric(value, name)
        for name, value in m["metrics"].items():
            finite_metric(value, name)
        expected_occurrences.extend((a, c) for c in candidates)
    io.same(len(occurrences), len(expected_occurrences), "normalized occurrence count")
    io.same(len(joined), len(occurrences), "common occurrence count")
    io.same(len(retained), len(occurrences), "retention occurrence count")
    scores, cache_names = {}, set()
    for o, (a, c), j, r in zip(occurrences, expected_occurrences, joined, retained):
        for field in FIELDS:
            io.same(o[field], a[field], "occurrence " + field)
        io.same(o["candidate"], c, "normalized candidate equals stored evaluation")
        io.same(o["rank"], c["rank"], "raw rank")
        io.same(o["smiles"], c["smiles"], "stored exact candidate")
        io.require(isinstance(o["query_exact"], str) and o["query_exact"], "Stored query identity required")
        io.require(type(o["route_supplied"]) is bool, "Route supplied must be boolean")
        if o["route_supplied"]:
            io.sha_value(o["serialized_route_sha256"])
        else:
            io.same(o["serialized_route_sha256"], None, "Absent route cannot become failed synthesis")
        for field in stats.RETENTION_FIELDS:
            io.require(
                o["scaffold_retention"][field] is None or type(o["scaffold_retention"][field]) is bool,
                "Unassessed scaffold annotation must remain null",
            )
        projection = {k: o[k] for k in (*FIELDS, "rank", "query_exact", "smiles", "scaffold_retention")}
        io.same(r, projection | {"family": c["family"]}, "stored retention occurrence")
        payload = {
            "query_exact": o["query_exact"],
            "candidate_exact": o["smiles"],
            "source": common["source"],
            "conformers": common["conformers"],
            "seeds": SEEDS,
            "identity_policy": "protonation-tautomer-family-v1",
            "rdkit": "2025.09.6",
        }
        key = content_key(payload)
        io.same(
            j,
            {k: o[k] for k in (*FIELDS, "rank", "smiles", "query_exact")} | {"cache_key": key},
            "exact ordered cache join",
        )
        name = "cache/" + key + ".json"
        cache_names.add(name)
        score = io.parse(blobs["common/" + name])
        validate_score(score, payload, c, stats)
        io.require(stats.occurrence_key(o) not in scores, "Duplicate occurrence")
        scores[stats.occurrence_key(o)] = score
    expected_common = {
        "context.json",
        "recomputed-attempts.json",
        "occurrences.json",
        "coverage.json",
    } | cache_names
    io.same(set(cm["files_sha256"]), expected_common, "full common scoring membership")
    coverage = io.parse(blobs["common/coverage.json"])
    for name, value in {
        "attempts": 270,
        "valid_candidate_occurrences": len(occurrences),
        "unique_exact_pairs": len(cache_names),
        "successful_3d_occurrences": sum(s["all_3d_seeds_succeeded"] for s in scores.values()),
    }.items():
        io.same(coverage[name], value, "common coverage " + name)
    io.same(
        stats.attempt_metrics(attempts, occurrences, scores),
        metrics,
        "all stored attempt statistics",
        tolerance=1e-12,
    )
    integrity = io.parse(blobs["integrity-validation.json"])
    expected = {
        "passed": True,
        "total_attempts": 270,
        "candidate_occurrences": len(occurrences),
        "raw_candidate_joins_verified": True,
        "route_flags_verified": True,
        "scaffold_annotations_verified": True,
        "source_sha256": {n: io.sha(v) for n, v in blobs.items() if n.startswith(("analysis/", "common/"))},
        "output_sha256": {"normalized-occurrences.json": io.sha(blobs["normalized-occurrences.json"])},
    }
    for name, value in expected.items():
        io.same(integrity[name], value, "parent normalization integrity " + name)
    return metrics


def validate_score(score, payload, candidate, stats):
    io.same(score["payload"], payload, "exact cache payload")
    io.same(score["cache_key"], content_key(payload), "cache key")
    trials = score["trials"]
    io.same([t["seed"] for t in trials], SEEDS, "all seed outcomes/order")
    for t in trials:
        io.require(t["status"] in {"success", "failed", "unassessed_due_to_pair_error"}, "Unknown 3D outcome")
        if t["status"] == "success":
            for metric in ("hit_shape", "hit_pharmacophore", "score"):
                value = t["metrics"][metric]
                io.require(type(value) in {int, float} and math.isfinite(value), "Invalid successful score")
                if metric != "score":
                    io.require(0 <= value <= 1, "Unbounded shape/pharmacophore")
        else:
            io.require(isinstance(t.get("reason"), str) and t["reason"], "Missing 3D failure reason")
    success = all(t["status"] == "success" for t in trials)
    io.same(score["all_3d_seeds_succeeded"], success, "all-seed success flag")
    if score["pair_error"] is not None:
        io.require(isinstance(score["pair_error"], str) and score["pair_error"], "Missing pair error")
        io.same(score["scores"], {}, "pair error cannot supply scores")
        io.require(
            all(t["status"] == "unassessed_due_to_pair_error" for t in trials), "Pair error seed outcomes"
        )
    else:
        io.require(all(t["status"] != "unassessed_due_to_pair_error" for t in trials), "Missing pair error")
        for k in ("bit", "count"):
            value = candidate["hit_" + k + "_tanimoto"]
            io.require(
                type(value) in {int, float} and math.isfinite(value) and 0 <= value <= 1,
                "Invalid stored fingerprint similarity",
            )
            io.same(
                score["scores"][k + "_tanimoto"],
                candidate["hit_" + k + "_tanimoto"],
                "stored fingerprint score",
                tolerance=1e-12,
            )
    for k in ("hit_shape", "hit_pharmacophore", "score"):
        if success:
            values = [t["metrics"][k] for t in trials]
            io.same(
                score["scores"][k], float(statistics.median(values)), "stored seed median", tolerance=1e-12
            )
            io.same(
                score["scores"][k + "_stability"],
                stats.score_stability(values),
                "stored seed stability",
                tolerance=1e-12,
            )
        else:
            io.require(
                k not in score["scores"] and k + "_stability" not in score["scores"],
                "Partial scoring cannot supply aggregate",
            )


def validate(blobs, calculate=True):
    stats, grid, common, execution = check_fixed(blobs)
    reviewed(blobs)
    metrics = validate_records(blobs, stats, grid, common, execution)
    if not calculate:
        return None
    comparison = stats.comparisons(metrics, grid, stats.cluster_interval)
    io.same(
        comparison, io.parse(blobs["analysis/comparison.json"]), "entire recorded comparison", tolerance=1e-12
    )
    return comparison


def manifest_for(blobs):
    return {
        "schema_version": 1,
        "kind": "supplementary-270-statistical-replay",
        "scope": SCOPE,
        "files": {n: {"sha256": io.sha(v), "bytes": len(v)} for n, v in sorted(blobs.items())},
    }


def pack(inputs, analysis, common, parent, parent_sha256, output):
    output = io.no_symlinks(output)
    io.require(not output.exists(), "Output already exists")
    for root in (inputs, analysis, common, Path(parent).parent):
        io.require(not output.is_relative_to(io.no_symlinks(root)), "Output cannot be inside an input tree")
    raw = io.read_bytes(parent)
    io.same(io.sha(raw), io.sha_value(parent_sha256), "external final parent review SHA")
    review = io.parse(raw)
    io.require(
        review.get("status") == "passed" and review.get("complete") is True,
        "Final parent completion required before reading result archives",
    )
    mapping = {
        "parent-review.json": io.no_symlinks(parent),
        **{
            n: io.no_symlinks(inputs) / n
            for n in ("normalized-occurrences.json", "integrity-validation.json")
        },
        **{n: p for n, (p, _) in fixed_files().items()},
    }
    for prefix, root in (("analysis", analysis), ("common", common)):
        path = io.no_symlinks(root) / "manifest.json"
        manifest_raw = io.read_bytes(path)
        io.same(
            io.sha(manifest_raw),
            review["input_sha256"][prefix + "/manifest.json"],
            "parent-pinned archive before read",
        )
        mapping[prefix + "/manifest.json"] = path
        names = set(io.parse(manifest_raw)["files_sha256"])
        io.require("manifest.json" not in names, "Archive manifest cannot list itself")
        observed = set()
        for base, dirs, files in os.walk(root, followlinks=False):
            for name in dirs + files:
                checked = io.no_symlinks(Path(base) / name)
                if name in files:
                    observed.add(checked.relative_to(io.no_symlinks(root)).as_posix())
        io.same(observed, names | {"manifest.json"}, "input archive membership")
        for name in names:
            mapping[prefix + "/" + io.safe_relative(name)] = io.no_symlinks(root) / name
    for item in review["assessment_evidence"]:
        name = io.safe_relative(item["artifact"])
        mapping["parent-review-evidence/" + name] = Path(parent).parent / name
    blobs = io.snapshot(mapping)
    io.same(blobs["parent-review.json"], raw, "review changed before snapshot")
    validate(blobs)
    io.unchanged(mapping, blobs)
    output = io.new_directory(output)
    for name, value in blobs.items():
        path = output / io.safe_relative(name)
        path.parent.mkdir(parents=True, exist_ok=True)
        io.write_new(path, value)
    io.unchanged(mapping, blobs)
    manifest = io.encoded(manifest_for(blobs))
    io.write_new(output / "manifest.json", manifest)
    return {"status": "packed", "attempts": 270, "manifest_sha256": io.sha(manifest)}


def load_bundle(bundle, anchor):
    root = io.no_symlinks(bundle)
    raw = io.read_bytes(root / "manifest.json")
    io.same(io.sha(raw), io.sha_value(anchor), "external bundle manifest SHA")
    manifest = io.parse(raw)
    io.same(
        {k: v for k, v in manifest.items() if k != "files"},
        {k: v for k, v in manifest_for({}).items() if k != "files"},
        "manifest schema",
    )
    mapping = {io.safe_relative(n): root / io.safe_relative(n) for n in manifest["files"]}
    io.require("manifest.json" not in mapping, "Bundle manifest cannot list itself")
    observed = set()
    for base, dirs, files in os.walk(root, followlinks=False):
        for n in dirs + files:
            p = io.no_symlinks(Path(base) / n)
            if n in files:
                observed.add(p.relative_to(root).as_posix())
    io.same(observed, set(mapping) | {"manifest.json"}, "bundle membership")
    blobs = io.snapshot(mapping)
    io.same(manifest["files"], manifest_for(blobs)["files"], "bundle hashes/lengths")
    return mapping | {"manifest.json": root / "manifest.json"}, blobs | {"manifest.json": raw}, blobs


def check(bundle, anchor, output, calculate=True):
    output = io.no_symlinks(output)
    io.require(
        not output.exists() and not output.is_relative_to(io.no_symlinks(bundle)),
        "Output exists or inside bundle",
    )
    mapping, before, blobs = load_bundle(bundle, anchor)
    comparison = validate(blobs, calculate)
    io.unchanged(mapping, before)
    report = {
        "status": "passed",
        "attempts": 270,
        "scope": SCOPE,
        "comparison_recomputed": calculate,
        "manifest_sha256": anchor,
        "absolute_tolerance": 1e-12,
        "analysis_helper_sha256": HELPER_SHA,
        "input_before_after_sha256": {n: io.sha(v) for n, v in sorted(before.items())},
    }
    if calculate:
        import numpy as np

        report["runtime"] = {"python": sys.version, "numpy": np.__version__}
        output = io.new_directory(output)
        io.write_new(output / "comparison.json", io.encoded(comparison))
        io.unchanged(mapping, before)
        io.write_new(output / "replay-validation.json", io.encoded(report))
    else:
        io.write_new(output, io.encoded(report))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("pack")
    for n in ("inputs", "analysis", "common", "parent-review", "output"):
        p.add_argument("--" + n, type=Path, required=True)
    p.add_argument("--parent-review-sha256", required=True)
    for name in ("verify", "replay"):
        p = commands.add_parser(name)
        p.add_argument("--bundle", type=Path, required=True)
        p.add_argument("--expected-manifest-sha256", required=True)
        p.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "pack":
        result = pack(
            args.inputs,
            args.analysis,
            args.common,
            args.parent_review,
            args.parent_review_sha256,
            args.output,
        )
    else:
        result = check(args.bundle, args.expected_manifest_sha256, args.output, args.command == "replay")
    print(io.encoded({k: result[k] for k in ("status", "attempts", "manifest_sha256")}).decode(), end="")


if __name__ == "__main__":
    main()
