"""Fresh query/seed-blocked qualified benchmark with a complete planned denominator.

Two sequential owned process groups reuse the existing bounded_method controller.
Only fully validated catalog records are shared in measured administrative setup.
No prediction/search caches, retries, partial-native rescue, or budget changes.
"""

import argparse
import copy
from pathlib import Path
import resource
import sys
import time

import psutil

from benchmarks.research.protocol import digest, write_json, synthesis_asset_paths
from benchmarks.research.qualified_adapter import ARMS, read, benchmark_profile
from benchmarks.research.qualified_analysis import canonical_sha256, analyze, validate_plan
from benchmarks.research import qualified_collector
from benchmarks.research.runner import bounded_method, _termination_scope, _write_checkpoint

BUDGET = {"wall_seconds": 960, "native_wall_seconds": 360, "common_wall_seconds": 600, "max_memory_mb": 24000}
CATALOG_POLICY = "runner-preload-catalog-v1"


def parent_memory():
    """Parent point observation and lifetime high-water mark, not a tree peak."""
    try:
        current = psutil.Process().memory_info().rss / (1024 * 1024)
    except psutil.Error:
        current = None
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_mb = peak / (1024 * 1024 if sys.platform == "darwin" else 1024)
    return {"rss_mb": current, "lifetime_peak_rss_mb": peak_mb}


def synthesis_sources(source):
    root = Path(qualified_collector.__file__).resolve().parents[2] / "chembounce_synthesis"
    return {Path(path).name: value for path, value in source.items() if Path(path).parent == root}


def preload_catalogs(protocol, termination, setup):
    """One full catalog per pinned path/checksum, owned only by this run.

    Sidecar provenance is explicitly pinned or explicitly absent. Cancellation
    is cooperative between full validations; no attempt starts after a signal.
    """
    cached, methods = {}, {}
    expected_source = synthesis_sources(protocol["source_sha256"])
    for method in protocol["methods"]:
        if termination.signal_number is not None:
            raise InterruptedError("Interrupted during shared catalog setup")
        assets = method["resolved_assets"]
        path = Path(assets["catalog"]).resolve()
        expected = method["assets"]["catalog"]["sha256"]
        sidecar = path.with_suffix(path.suffix + ".manifest.json")
        if "catalog_manifest" in assets:
            if (
                method.get("catalog_manifest_absent", False) is not False
                or Path(assets["catalog_manifest"]).resolve() != sidecar
            ):
                raise ValueError("Catalog sidecar path/absence binding conflicts")
            sidecar_sha = method["assets"]["catalog_manifest"]["sha256"]
        else:
            if method.get("catalog_manifest_absent") is not True:
                raise ValueError("Catalog requires pinned sidecar or explicit absence binding")
            sidecar_sha = None
        key = str(path), expected, sidecar_sha
        if key not in cached:
            entry = {"catalog_path": str(path), "status": "started", "method_ids": []}
            setup["catalogs"].append(entry)
            started, cpu = time.monotonic(), time.process_time()
            entry["parent_before"] = parent_memory()
            try:
                cached[key] = qualified_collector.BoundCatalog(
                    path, expected, expected_source, sidecar_sha256=sidecar_sha
                )
                entry.update(binding=cached[key].binding(), status="completed")
            except BaseException as exc:
                entry.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
                raise
            finally:
                entry.update(
                    wall_seconds=time.monotonic() - started,
                    cpu_seconds=time.process_time() - cpu,
                    parent_after=parent_memory(),
                )
        else:
            # Recheck even when methods reuse the same preloaded snapshot.
            cached[key].validate(path, expected_source)
        next(e for e in setup["catalogs"] if e.get("binding") == cached[key].binding())["method_ids"].append(
            method["id"]
        )
        methods[method["id"]] = cached[key]
        if termination.signal_number is not None:
            raise InterruptedError("Interrupted during shared catalog setup")
    return methods


def make_plan(protocol, protocol_sha256):
    queries = protocol["queries"]
    if any(set(q) != {"id", "target_id", "smiles"} for q in queries):
        raise ValueError("Only public query identity/SMILES enter the benchmark")
    plan = {
        "schema": "qualified-benchmark-plan-v1",
        "run_id": protocol["run_id"],
        "protocol_sha256": protocol_sha256,
        "methods": [m["id"] for m in protocol["methods"]],
        "queries": [
            {"target_id": q["target_id"], "query_id": q["id"], "query_sha256": canonical_sha256(q)}
            for q in queries
        ],
        "seeds": protocol["seeds"],
        "collector_source_sha256": digest(qualified_collector.__file__),
        "budget": protocol["budget"],
        "source_sha256": protocol["source_sha256"],
        "collector_catalog_policy": protocol.get("collector_catalog_policy"),
    }
    if protocol.get("schema") == "qualified-benchmark-protocol-v2":
        plan["benchmark_protocol"] = {
            key: protocol[key] for key in ("schema", "recommendation_profile", "profile_artifact_sha256")
        }
        benchmark_profile(plan["benchmark_protocol"])
    plan["order_rule"] = "query/seed blocks; cyclic method rotation by block index"
    plan["slot_order"] = []
    for block, (query, seed) in enumerate((q, seed) for q in plan["queries"] for seed in plan["seeds"]):
        offset = block % len(plan["methods"])
        order = plan["methods"][offset:] + plan["methods"][:offset]
        plan["slot_order"].extend(
            {
                "method_id": method,
                "target_id": query["target_id"],
                "query_id": query["query_id"],
                "seed": seed,
            }
            for method in order
        )
    validate_plan(plan)
    return plan


def resolve_protocol(protocol, base):
    protocol = copy.deepcopy(protocol)
    if (
        protocol.get("schema") not in {"qualified-benchmark-protocol-v1", "qualified-benchmark-protocol-v2"}
        or protocol["budget"] != BUDGET
    ):
        raise ValueError("Unknown protocol or unfrozen stage budgets")
    if protocol["schema"] == "qualified-benchmark-protocol-v2":
        benchmark_profile(
            {
                key: protocol.get(key)
                for key in ("schema", "recommendation_profile", "profile_artifact_sha256")
            }
        )
    elif "recommendation_profile" in protocol or "profile_artifact_sha256" in protocol:
        raise ValueError("Explicit profile binding requires benchmark protocol v2")
    if protocol.get("collector_catalog_policy") != CATALOG_POLICY:
        raise ValueError("Protocol must explicitly freeze the shared catalog setup policy")
    source = protocol["source_sha256"]
    root = Path(__file__).resolve().parents[2]
    required = list((root / "chembounce_synthesis").glob("*.py")) + [
        Path(__file__),
        Path(qualified_collector.__file__),
        root / "benchmarks/research/qualified_adapter.py",
        root / "benchmarks/research/qualified_analysis.py",
        root / "benchmarks/research/runner.py",
        root / "chembounce.py",
    ]
    if any(str(path.resolve()) not in source for path in required):
        raise ValueError("Protocol omits required current source pins")
    files_seen, directories_seen = {}, {}

    def file_hash(path):
        path = Path(path).resolve()
        if path not in files_seen:
            files_seen[path] = digest(path)
        return files_seen[path]

    if any(file_hash(path) != value for path, value in source.items()):
        raise ValueError("Protocol source changed")
    for method in protocol["methods"]:
        if method["arm"] not in ARMS:
            raise ValueError("Unknown arm")
        assets = {}
        for name, spec in method["assets"].items():
            if "directory" in spec:
                path = (Path(base) / spec["directory"]).resolve()
                if path not in directories_seen:
                    directories_seen[path] = {
                        str(p.relative_to(path)): file_hash(p) for p in path.rglob("*") if p.is_file()
                    }
                actual = directories_seen[path]
                if actual != spec["files_sha256"] or not actual:
                    raise ValueError("Asset directory changed")
            else:
                path = (Path(base) / spec["artifact"]).resolve()
                if file_hash(path) != spec["sha256"]:
                    raise ValueError("Asset file changed")
            assets[name] = str(path)
        for key in ("synthesis_config", "catalog", "rules"):
            if key not in assets:
                raise ValueError("Common synthesis assets missing")
        for name, path in synthesis_asset_paths(assets["synthesis_config"]).items():
            if path not in assets.values():
                raise ValueError(f"Unpinned common nested asset: {name}")
        from chembounce_synthesis.workflow import load_config

        expected_profile = benchmark_profile(
            {key: protocol[key] for key in ("schema", "recommendation_profile", "profile_artifact_sha256")}
            if protocol["schema"] == "qualified-benchmark-protocol-v2"
            else None
        )
        if load_config(assets["synthesis_config"]).get("recommendation_profile") != expected_profile:
            raise ValueError("Method config profile differs from benchmark protocol")
        method["resolved_assets"] = assets
        if method["arm"] in {"baseline", "crem", "mol2mol"}:
            native = method["native_method"]
            if native["kind"] != "external" or not all(isinstance(s, str) for s in native["command"]):
                raise ValueError("Native comparison requires an explicit pinned external command")
            if not any("{request}" in a for a in native["command"]) or not any(
                "{response}" in a for a in native["command"]
            ):
                raise ValueError("Native command must consume/produce bounded request/response")
            for alias in native["asset_names"].values():
                if alias not in assets:
                    raise ValueError("Native asset alias is missing")
    return protocol


def adapter_command(phase):
    return [
        sys.executable,
        "-m",
        "benchmarks.research.qualified_adapter",
        "--phase",
        phase,
        "--request",
        "{request}",
        "--response",
        "{response}",
    ]


def request_for(method, query, seed, output, seconds, phase, plan):
    native = method.get("native_method", {}) if phase == "native" and method["arm"] != "full" else {}
    assets = method["resolved_assets"]
    if native:
        assets = {name: assets[alias] for name, alias in native["asset_names"].items()}
    cap = 12 if phase == "native" and method["arm"] == "full" else 3
    evaluations = 256 if phase == "native" and method["arm"] in {"crem", "mol2mol"} else 12
    settings = copy.deepcopy(native.get("settings", method.get("settings", {})))
    if not native:
        settings["arm"] = method["arm"]
    # Keep native schemas unchanged. The outer phase manifest binds the plan.
    return {
        "schema_version": 1,
        **(
            {"benchmark_protocol": plan["benchmark_protocol"]}
            if phase == "common" and "benchmark_protocol" in plan
            else {}
        ),
        "kind": native.get("request_kind", "external"),
        "method_id": method["id"],
        "query": query,
        "seed": seed,
        "budget": {
            "wall_seconds": seconds,
            "max_candidates": cap,
            "max_evaluations": evaluations,
            "max_memory_mb": BUDGET["max_memory_mb"],
        },
        "assets": assets,
        "settings": settings,
        "output_dir": str(output.resolve()),
    }


def file_manifest(directory):
    return {
        str(path.relative_to(directory)): digest(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }


def _unit_count(value):
    return type(value) is int and value >= 0


def _accounting_object(value):
    return value if isinstance(value, dict) else {}


def unit_counts(directory):
    """Report observed units; omit unknown totals rather than imputing zero.

    A sum requires an explicit collection and a valid count in every relevant
    row. An empty collection records zero; it does not prove that no unrecorded
    work occurred before interruption. Other known counters remain usable when
    a particular count is absent or malformed. Unreadable JSON still follows
    the attempt's accounting-error path, preserving measured resources and
    the original attempt status. This cannot recover missing measurements.
    """
    out = directory / "synthesis"
    summary = _accounting_object(read(out / "summary.json")) if (out / "summary.json").is_file() else {}
    counts = {}
    for key in ("scored_products", "proposals", "verified_routes", "designs_used", "designs_supplied"):
        if _unit_count(summary.get(key)):
            counts[key] = summary[key]
    selection = _accounting_object(summary.get("design_selection"))
    for key, value in _accounting_object(selection.get("counts")).items():
        if _unit_count(value):
            counts["design_" + key] = value
    for key, value in _accounting_object(summary.get("synthesis_evidence_usage")).items():
        if _unit_count(value):
            counts["forward_" + key] = value
    search_reports = summary.get("search_reports")
    if isinstance(search_reports, list) and all(
        isinstance(r, dict) and r.get("stage") in ("projection", "expansion") for r in search_reports
    ):
        for stage in ("projection", "expansion"):
            reports = [r for r in search_reports if r["stage"] == stage]
            if any(type(r.get("cached", False)) is not bool for r in reports):
                continue
            reports = [r for r in reports if not r.get("cached", False)]
            counts[stage + "_search_reports"] = len(reports)
            for key in (
                "combinations",
                "root_state_visits",
                "root_states_generated",
                "root_states_truncated",
                "root_objective_evaluations",
                "root_objective_cache_hits",
                "reaction_truncations",
                "combination_truncations",
                "beam_truncations",
            ):
                if all(_unit_count(r.get(key)) for r in reports):
                    counts[stage + "_" + key] = sum(r[key] for r in reports)
    refinement = (
        _accounting_object(read(out / "refinement.json")) if (out / "refinement.json").is_file() else {}
    )
    candidates = refinement.get("candidates")
    if isinstance(candidates, list) and all(
        isinstance(r, dict) and _unit_count(r.get("candidate_conformer_requests")) for r in candidates
    ):
        counts["refinement_candidate_conformer_requests"] = sum(
            r["candidate_conformer_requests"] for r in candidates
        )
    limit = _accounting_object(summary.get("refinement_budget")).get("limit")
    if _unit_count(limit):
        counts["refinement_candidate_conformer_limit"] = limit
    proposals = out / "proposals.json"
    if proposals.is_file():
        returned = _accounting_object(read(proposals)).get("proposals")
        if isinstance(returned, list) and all(isinstance(r, dict) for r in returned):
            counts["returned_proposal_records"] = len(returned)
            tokens = [_accounting_object(r.get("sampling")).get("generated_tokens") for r in returned]
            counts["returned_proposals_with_token_counts"] = sum(_unit_count(n) for n in tokens)
            if all(_unit_count(n) for n in tokens):
                counts["returned_proposal_generated_tokens"] = sum(tokens)
    return counts


def attempt(method, query, seed, output, plan, termination, *, bound_catalog=None):
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    outer_deadline = started + BUDGET["wall_seconds"]
    parent_cpu_started, parent_before = time.process_time(), parent_memory()
    identity = {
        "method_id": method["id"],
        "target_id": query["target_id"],
        "query_id": query["id"],
        "query_sha256": canonical_sha256(query),
        "seed": seed,
    }
    result = {
        **identity,
        "status": "error",
        "qualified": [],
        "failures": [],
        "costs": {
            "wall_seconds": None,
            "cpu_seconds": None,
            "gpu_seconds": None,
            "peak_rss_mb": None,
            "peak_gpu_mb": None,
            "native_units": {},
            "common_units": {},
        },
    }
    phases, native_files = {}, None
    write_json(
        output / "attempt-binding.json",
        {
            **identity,
            "plan_sha256": canonical_sha256(plan),
            "outer_started_monotonic": started,
            "outer_deadline_monotonic": outer_deadline,
        },
    )
    try:
        if method["arm"] != "synllama":
            native_dir = output / "native"
            seconds = BUDGET["native_wall_seconds"]
            request = request_for(method, query, seed, native_dir, seconds, "native", plan)
            command = (
                adapter_command("native") if method["arm"] == "full" else method["native_method"]["command"]
            )
            phases["native"] = bounded_method(
                {"kind": "external", "command": command}, request, native_dir, _termination=termination
            )
            write_json(output / "native-report.json", phases["native"])
            if phases["native"]["status"] not in {"success", "no_candidates"}:
                raise RuntimeError("Native phase did not complete; partial supply is ineligible")
            if phases["native"].get("resources", {}).get("wall_seconds", 0) > BUDGET["native_wall_seconds"]:
                raise TimeoutError("Native phase including cleanup exceeded its allocation")
            native_files = file_manifest(native_dir)
        if termination.signal_number is not None:
            result["status"] = "interrupted"
            return result
        common_dir = output / "common"
        common_started = time.monotonic()
        seconds = min(BUDGET["common_wall_seconds"], outer_deadline - common_started)
        if seconds <= 0:
            raise TimeoutError("Outer budget exhausted before common phase")
        request = request_for(method, query, seed, common_dir, seconds, "common", plan)
        if native_files is not None:
            request["native"] = {"directory": str(native_dir.resolve()), "files_sha256": native_files}
        phases["common"] = bounded_method(
            {"kind": "external", "command": adapter_command("common")},
            request,
            common_dir,
            _termination=termination,
        )
        write_json(output / "common-report.json", phases["common"])
        if phases["common"]["status"] not in {"success", "no_candidates"}:
            raise RuntimeError("Common phase did not complete")
        if phases["common"].get("resources", {}).get("wall_seconds", 0) > seconds:
            phases["common"]["allocation_exceeded"] = True
            raise TimeoutError("Common phase including cleanup exceeded its requested allocation")
        termination.phase = "qualified_collector"
        collector_started = time.monotonic()
        if collector_started >= outer_deadline:
            raise TimeoutError("Outer budget exhausted before collector")
        expected_source = synthesis_sources(plan["source_sha256"])
        result["qualified"] = qualified_collector.checked_products(
            common_dir,
            phases["common"]["response"],
            query,
            expected_source,
            deadline=outer_deadline,
            bound_catalog=bound_catalog,
            **({"benchmark_protocol": plan["benchmark_protocol"]} if "benchmark_protocol" in plan else {}),
        )
        phases["collector"] = {
            "wall_seconds": time.monotonic() - collector_started,
            "scope": "offline integrity attestation; no new inference",
        }
        if time.monotonic() >= outer_deadline:
            raise TimeoutError("Outer budget exhausted including collector")
        if termination.signal_number is not None:
            result["qualified"], result["status"] = [], "interrupted"
        else:
            result["status"] = "success" if result["qualified"] else "no_candidates"
            if not result["qualified"]:
                summary = (
                    read(common_dir / "synthesis/summary.json")
                    if (common_dir / "synthesis/summary.json").is_file()
                    else {}
                )
                result["failures"].append(
                    {
                        "stage": "common",
                        "reason": summary.get(
                            "no_candidates_reason", "No jointly qualified product within declared limits"
                        ),
                    }
                )
                if summary.get("query_preflight", {}).get("status") == "rejected":
                    result["status"] = "unsupported"
    except (ValueError, KeyError, TypeError, RuntimeError, OSError, TimeoutError) as exc:
        statuses = [p.get("status") for p in phases.values()]
        result["status"] = (
            "interrupted"
            if termination.signal_number is not None
            else "timeout"
            if isinstance(exc, TimeoutError) or "timeout" in statuses
            else "error"
        )
        result["qualified"] = []
        result["failures"].append(
            {
                "stage": "collector"
                if "common" in phases
                and phases["common"]["status"] in {"success", "no_candidates"}
                and not phases["common"].get("allocation_exceeded", False)
                else "common"
                if "common" in phases
                else "native",
                "reason": str(exc),
            }
        )
    finally:
        result["costs"]["wall_seconds"] = time.monotonic() - started
        resources = [r["resources"] for r in phases.values() if "resources" in r]
        if resources:
            result["costs"]["cpu_seconds"] = sum(r["sampled_cpu_seconds_lower_bound"] for r in resources)
            result["costs"]["peak_rss_mb"] = max(r["sampled_peak_rss_mb"] for r in resources)
        for phase in ("native", "common"):
            if (output / phase).exists():
                try:
                    result["costs"][phase + "_units"] = unit_counts(output / phase)
                except (ValueError, KeyError, TypeError, OSError) as exc:
                    result["failures"].append(
                        {
                            "stage": phase + "_cost_accounting",
                            "reason": "Partial artifacts; native/common units unavailable: " + str(exc),
                        }
                    )
                if phase in phases and phases[phase].get("response"):
                    response = phases[phase]["response"]
                    if _unit_count(response.get("evaluations")):
                        result["costs"][phase + "_units"]["returned_evaluations"] = response["evaluations"]
        write_json(
            output / "phases.json",
            {
                "phases": phases,
                "resources_scope": "CPU sampled lower bound; GPU time/memory unavailable; no zero imputation",
                "native_files_sha256": native_files,
                "parent_resources": {
                    "cpu_seconds": time.process_time() - parent_cpu_started,
                    "before": parent_before,
                    "after": parent_memory(),
                    "scope": "Parent process only, including retained shared catalog; excludes child CPU/RSS. RSS points and lifetime high-water mark are not simultaneous tree peaks.",
                    "simultaneous_parent_plus_child_peak_rss_mb": None,
                },
            },
        )
        result["costs"]["wall_seconds"] = time.monotonic() - started
        if time.monotonic() >= outer_deadline and result["status"] not in {"timeout", "interrupted"}:
            result["status"], result["qualified"] = "timeout", []
            result["failures"].append(
                {"stage": "accounting", "reason": "Outer budget exhausted through final accounting"}
            )
        write_json(output / "collection-attempt.json", result)
    return result


def run(protocol_path, output):
    campaign_started = time.monotonic()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    raw = read(protocol_path)
    plan = make_plan(raw, digest(protocol_path))
    write_json(output / "plan.json", plan)
    collection = {
        "schema": "qualified-benchmark-collection-v1",
        "plan_sha256": canonical_sha256(plan),
        "collector_source_sha256": plan["collector_source_sha256"],
        "attempts": [],
        "administrative_setup": {"status": "not_started"},
        "measurement_policy": {
            "cpu_seconds": "sum of phase sampled process-tree CPU lower bounds; collector CPU excluded",
            "peak_rss_mb": "maximum sampled aggregate process-tree RSS across phases; collector RSS excluded",
            "gpu_seconds": "not measured; null",
            "peak_gpu_mb": "not measured; null",
            "wall_seconds": "per-attempt outer elapsed including native/common cleanup and collector; shared administrative setup reported separately, never free compute",
            "shared_setup": "Full catalog validation/preload plus integrity preflight outside attempt clocks; measured wall/parent CPU/RSS. No prediction/search/rule-result cache shared.",
            "parent_memory": "Setup and per-attempt parent point RSS/lifetime high-water marks are separate from sampled child-tree peaks. Simultaneous parent+child peak unmeasured; separate maxima are not summed.",
            "memory_limit": "Existing 24000MB controller applies to each child process tree, not a measured combined parent+child peak.",
        },
    }

    def checkpoint():
        collection["campaign_wall_seconds_through_checkpoint"] = time.monotonic() - campaign_started
        _write_checkpoint(output / "collection.json", collection)
        report = analyze(plan, collection)
        report["measurement_policy"] = collection["measurement_policy"]
        report["administrative_setup"] = collection["administrative_setup"]
        report["campaign_wall_seconds_through_checkpoint"] = collection[
            "campaign_wall_seconds_through_checkpoint"
        ]
        _write_checkpoint(output / "analysis.json", report)

    checkpoint()
    with _termination_scope() as termination:
        try:
            setup_started, setup_cpu = time.monotonic(), time.process_time()
            setup = {
                "policy": CATALOG_POLICY,
                "status": "started",
                "catalogs": [],
                "parent_before": parent_memory(),
                "scope": "Shared full catalog validation/preload and integrity preflight outside attempt clocks. Common/native children remain independent cold invocations and pay their own catalog construction costs.",
            }
            catalogs = None
            try:
                if termination.signal_number is not None:
                    raise InterruptedError("Interrupted before shared setup")
                protocol = resolve_protocol(raw, Path(protocol_path).resolve().parent)
                catalogs = preload_catalogs(protocol, termination, setup)
                setup["status"] = "completed"
            except InterruptedError as exc:
                setup.update(status="interrupted", error={"type": type(exc).__name__, "message": str(exc)})
            except BaseException as exc:
                setup.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
                raise
            finally:
                setup.update(
                    wall_seconds=time.monotonic() - setup_started,
                    cpu_seconds=time.process_time() - setup_cpu,
                    parent_after=parent_memory(),
                    simultaneous_parent_plus_child_peak_rss_mb=None,
                )
                setup["accounting_share_per_planned_attempt"] = {
                    "denominator": len(plan["slot_order"]),
                    "wall_seconds": setup["wall_seconds"] / len(plan["slot_order"]),
                    "cpu_seconds": setup["cpu_seconds"] / len(plan["slot_order"]),
                    "scope": "Administrative cost allocation only; includes every planned slot and never extends any attempt budget.",
                }
                collection["administrative_setup"] = setup
                write_json(output / "setup.json", setup)
            if catalogs is None or termination.signal_number is not None:
                return collection
            methods = {m["id"]: m for m in protocol["methods"]}
            queries = {(q["target_id"], q["id"]): q for q in protocol["queries"]}
            for index, slot in enumerate(plan["slot_order"]):
                if termination.signal_number is not None:
                    break
                row = attempt(
                    methods[slot["method_id"]],
                    queries[(slot["target_id"], slot["query_id"])],
                    slot["seed"],
                    output / f"attempt-{index:04d}",
                    plan,
                    termination,
                    bound_catalog=catalogs[slot["method_id"]],
                )
                collection["attempts"].append(row)
                checkpoint()
        finally:
            checkpoint()
            write_json(
                output / "collection-binding.json",
                {
                    "collection_sha256": digest(output / "collection.json"),
                    "plan_sha256": canonical_sha256(plan),
                    "termination": termination.record(),
                },
            )
    return collection


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    run(args.protocol, args.output)


if __name__ == "__main__":
    main()
