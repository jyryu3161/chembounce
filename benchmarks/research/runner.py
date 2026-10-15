"""Run matched-budget methods without passing held-out labels to generators.

Usage: python -m benchmarks.research.runner --dataset panel.json --protocol methods.json --output fresh-dir
"""

import argparse
from contextlib import contextmanager
import copy
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import threading
import time

import psutil

from benchmarks.research.evaluation import aggregate_runs, evaluate_candidates
from benchmarks.research.protocol import (
    audit_inventory,
    digest,
    load_dataset,
    verified_artifact,
    write_json,
    synthesis_asset_paths,
)
from chembounce_synthesis.identity import chemical_identity


KINDS = {"nearest_catalog", "nearest_3d", "reaction_enumeration", "synllama_only", "pipeline", "external"}


class BenchmarkInterrupted(BaseException):
    """A signal stopped a campaign; completed and interrupted attempts are saved.

    Library callers may catch this explicitly and inspect ``runs``. Inheriting
    BaseException prevents an ordinary error/retry handler from restarting work.
    """

    def __init__(self, termination, runs):
        self.termination = termination
        self.signal_number = termination["signal_number"]
        self.exit_code = termination["exit_code"]
        self.runs = runs
        super().__init__(f"Benchmark interrupted by {termination['signal_name']}")


class _Termination:
    def __init__(self):
        self.signal_number = None
        self.phase = "setup"
        self.received_phase = None

    def receive(self, signum, frame):
        # Do not raise asynchronously: that could lose a just-created Popen
        # handle or interrupt process-group cleanup / checkpoint replacement.
        if self.signal_number is None:
            self.signal_number = signum
            self.received_phase = self.phase

    def record(self):
        if self.signal_number is None:
            return None
        return {
            "signal_number": self.signal_number,
            "signal_name": signal.Signals(self.signal_number).name,
            "exit_code": 128 + self.signal_number,
            "received_phase": self.received_phase,
        }


@contextmanager
def _termination_scope(existing=None):
    """Temporarily own termination signals on the main thread, including cleanup."""
    if threading.current_thread() is not threading.main_thread():
        raise RuntimeError("Benchmark signal handling requires the main thread")
    if existing is not None:
        yield existing
        return
    state = _Termination()
    previous = {}
    try:
        for signum in (signal.SIGTERM, signal.SIGINT):
            previous[signum] = signal.getsignal(signum)
            signal.signal(signum, state.receive)
        yield state
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def _write_checkpoint(path, value):
    """Readers see either the old complete JSON or its complete replacement."""
    path = Path(path)
    temporary = path.with_name(f".{path.name}.tmp")
    write_json(temporary, value)
    temporary.replace(path)


def _positive(value, name, integer=False):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ValueError(f"{name} must be positive and finite")
    if integer and not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    return value


def load_protocol(path, dataset):
    path = Path(path).resolve()
    data = json.loads(path.read_text())
    if data.get("schema_version") != 1:
        raise ValueError("Benchmark protocol requires schema_version=1")
    budget = data.get("budget", {})
    for key in ("wall_seconds", "max_candidates", "max_evaluations"):
        _positive(budget.get(key), key, integer=key != "wall_seconds")
    if "max_memory_mb" in budget:
        _positive(budget["max_memory_mb"], "max_memory_mb")
    seeds = data.get("seeds")
    if (
        not isinstance(seeds, list)
        or not seeds
        or len(seeds) != len(set(seeds))
        or any(isinstance(s, bool) or not isinstance(s, int) or not 0 <= s < 2**31 for s in seeds)
    ):
        raise ValueError("seeds must be distinct nonnegative integers below 2**31")
    methods, seen, library_audits, verified_directories = [], set(), {}, {}
    for method in data.get("methods", []):
        method = copy.deepcopy(method)
        ident = method.get("id")
        if (
            not isinstance(ident, str)
            or not ident
            or ident in seen
            or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in ident)
        ):
            raise ValueError("Method IDs must be unique safe path components")
        seen.add(ident)
        if method.get("kind") not in KINDS:
            raise ValueError("Unknown method kind")
        overrides = method.get("budget_overrides", {})
        if not isinstance(overrides, dict) or set(overrides) - {"max_evaluations"}:
            raise ValueError("Only max_evaluations may differ between methods; wall/output caps stay matched")
        if "max_evaluations" in overrides:
            _positive(overrides["max_evaluations"], "method max_evaluations", integer=True)
            if overrides["max_evaluations"] > budget["max_evaluations"]:
                raise ValueError("Method evaluation cap exceeds common outer cap")
        if method.get("kind") == "external":
            argv = method.get("command")
            if not isinstance(argv, list) or not argv or not all(isinstance(a, str) for a in argv):
                raise ValueError("External method requires command argument list")
            if not any("{request}" in a for a in argv) or not any("{response}" in a for a in argv):
                raise ValueError("External command must explicitly consume request and produce response")
        assets = {}
        for name, spec in method.get("assets", {}).items():
            if "directory" in spec:
                # Hash all files, not a directory name or a provider version claim.
                directory = (path.parent / spec["directory"]).resolve()
                if directory not in verified_directories:
                    verified_directories[directory] = {
                        str(p.relative_to(directory)): digest(p) for p in directory.rglob("*") if p.is_file()
                    }
                actual = verified_directories[directory]
                if not actual or actual != spec.get("files_sha256"):
                    raise ValueError("Asset directory contents do not match pinned manifest")
                assets[name] = str(directory)
            else:
                assets[name] = str(verified_artifact(spec, path.parent))
        required = {
            "nearest_catalog": {"catalog"},
            "nearest_3d": {"catalog"},
            "reaction_enumeration": {"catalog", "rules"},
            "synllama_only": {"catalog", "rules", "model_directory", "prompt"},
            "pipeline": {"catalog", "rules", "synthesis_config"},
            "external": set(),
        }[method["kind"]]
        if not required <= assets.keys():
            raise ValueError(f"Missing assets for {ident}: {sorted(required - assets.keys())}")
        if method["kind"] == "pipeline":
            if not {"scaffold_db", "fingerprint_db"} & assets.keys():
                raise ValueError("Pipeline requires a pinned scaffold or fingerprint database")
            for name, expected in synthesis_asset_paths(assets["synthesis_config"]).items():
                if assets.get(name) != expected:
                    raise ValueError(f"Pipeline nested asset must be independently pinned: {name}")
            overrides = method.get("settings", {}).get("synthesis_overrides", {})
            if set(overrides) & {
                "provider",
                "catalog",
                "rules",
                "reference_sdf",
                "carry_forward",
            } or overrides.get("synthesis_evidence", {}).get("predictor"):
                raise ValueError("Pipeline overrides cannot replace pinned external assets")
            generator = method.get("settings", {}).get("generator_options", {})
            if set(generator) & {
                "scaffold_db_file",
                "fingerprint_db_file",
                "replace_scaffold_files",
                "fragments_DB",
            }:
                raise ValueError("Generator file inputs must use declared pinned assets")
        method["resolved_assets"] = assets
        inventories = method.get("inventories", [])
        if len({i["kind"] for i in inventories}) != len(inventories):
            raise ValueError("Declare at most one training and one prompt inventory per method")
        if method["kind"] in {"pipeline", "synllama_only", "external"}:
            for kind in ("training", "prompt"):
                if kind not in {i["kind"] for i in inventories}:
                    inventories.append(
                        {"kind": kind, "status": "unknown", "reason": "Inventory not supplied"}
                    )
        method["overlap_audit"] = [audit_inventory(spec, dataset, path.parent) for spec in inventories]
        method["search_space_overlap_audit"] = []
        for name in ("catalog", "scaffold_db"):
            if name in method.get("assets", {}):
                audit_key = method["assets"][name]["sha256"]
                if audit_key not in library_audits:
                    library_audits[audit_key] = audit_inventory(
                        method["assets"][name]
                        | {
                            "kind": "training",
                            "status": "available",
                            "coverage": "complete",
                            "max_records": data.get("search_space_audit_max_records", 5000),
                        },
                        dataset,
                        path.parent,
                    )
                audit = copy.deepcopy(library_audits[audit_key])
                audit.update(
                    kind="search_space",
                    asset=name,
                    interpretation="Presence in a declared retrieval/search library; distinct from pretrained-model training leakage",
                )
                method["search_space_overlap_audit"].append(audit)
        methods.append(method)
    if not methods:
        raise ValueError("No methods configured; unexecuted baselines must not appear as results")
    data["methods"] = methods
    return data


def public_request(dataset, query, method, seed, budget, output):
    record = next(r for r in dataset["compounds"] if r["id"] == query["compound_id"])
    # Whitelist fields: no labelled reference structures, assay labels, split
    # membership, source data paths, or reference target sets enter this request.
    result = {
        "schema_version": 1,
        "kind": method["kind"],
        "method_id": method["id"],
        "query": {"id": query["id"], "smiles": record["smiles"], "target_id": record["target_id"]},
        "seed": seed,
        "budget": copy.deepcopy(budget) | method.get("budget_overrides", {}),
        "assets": method["resolved_assets"],
        "settings": copy.deepcopy(method.get("settings", {})),
        "output_dir": str(Path(output).resolve()),
    }
    if query.get("resolved_reference_sdf"):
        result["query"]["reference_sdf"] = query["resolved_reference_sdf"]
        result["query"]["reference_sdf_sha256"] = query["reference_sdf"]["sha256"]
    return result


def _stop_group(process):
    cleanup = {
        "process_group_id": process.pid,
        "sigterm_sent": False,
        "sigkill_sent": False,
    }
    if process.poll() is None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
            cleanup["sigterm_sent"] = True
            process.wait(timeout=1)
        except (ProcessLookupError, subprocess.TimeoutExpired):
            pass
    # A worker may exit while leaving grandchildren running. Its process group
    # must be killed even if the immediate child has already become terminal.
    try:
        os.killpg(process.pid, signal.SIGKILL)
        cleanup["sigkill_sent"] = True
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=1)
    except subprocess.TimeoutExpired:
        # Even SIGKILL cannot immediately reap an uninterruptible kernel task.
        # Keep finalization bounded and report this fact instead of hanging.
        pass
    cleanup["leader_reaped"] = process.poll() is not None
    cleanup["leader_returncode"] = process.returncode
    return cleanup


def bounded_method(method, request, output, *, _termination=None):
    """Run on the main thread with temporary SIGTERM/SIGINT handling.

    An interrupted worker returns an empty response and its sampled resources.
    ``termination`` records a received signal even if it arrived after response
    validation, when the already completed outcome is preserved. Prior handlers
    are restored on success, interruption and unrelated exceptions alike.
    """
    with _termination_scope(_termination) as termination:
        report = _bounded_method(method, request, output, termination)
    if termination.signal_number is not None:
        report["termination"] = termination.record()
    return report


def _bounded_method(method, request, output, termination):
    termination.phase = "worker_setup"
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    request_path, response_path = output / "request.json", output / "response.json"
    start = time.monotonic()
    request = copy.deepcopy(request)
    request["deadline_monotonic"] = start + request["budget"]["wall_seconds"]
    request["finalization_margin_seconds"] = min(2.0, request["budget"]["wall_seconds"] * 0.05)
    write_json(request_path, request)
    if method["kind"] == "external":
        argv = [
            a.replace("{request}", str(request_path)).replace("{response}", str(response_path))
            for a in method["command"]
        ]
    else:
        argv = [
            sys.executable,
            "-m",
            "benchmarks.research.adapters",
            "--request",
            str(request_path),
            "--response",
            str(response_path),
        ]
    peak, cpu_max, monitor_samples = 0, 0, 0
    status, response, error, cleanup = None, {}, None, None
    process = None
    with (output / "stdout.log").open("w") as stdout, (output / "stderr.log").open("w") as stderr:
        try:
            if termination.signal_number is not None:
                status = "interrupted"
            else:
                process = subprocess.Popen(
                    argv,
                    stdout=stdout,
                    stderr=stderr,
                    start_new_session=True,
                    cwd=Path(__file__).resolve().parents[2],
                )
        except OSError as exc:
            return {
                "status": "launch_error",
                "response": {},
                "error": str(exc),
                "command": argv,
                "resources": {
                    "wall_seconds": time.monotonic() - start,
                    "sampled_peak_rss_mb": 0,
                    "sampled_cpu_seconds_lower_bound": 0,
                    "monitor_samples": 0,
                    "gpu_memory": "not_measured",
                    "includes_import_and_model_loading": True,
                },
            }
        try:
            termination.phase = "worker_monitor"
            while process is not None and process.poll() is None:
                if termination.signal_number is not None:
                    status = "interrupted"
                    break
                elapsed = time.monotonic() - start
                if elapsed >= request["budget"]["wall_seconds"]:
                    status = "timeout"
                    break
                try:
                    parent = psutil.Process(process.pid)
                    processes = [parent] + parent.children(recursive=True)
                    rss, cpu = 0, 0
                    for child in processes:
                        try:
                            rss += child.memory_info().rss
                            times = child.cpu_times()
                            cpu += times.user + times.system
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            pass
                    peak = max(peak, rss / 1024**2)
                    cpu_max = max(cpu_max, cpu)
                    monitor_samples += 1
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
                if peak > request["budget"].get("max_memory_mb", float("inf")):
                    status = "memory_limit"
                    break
                time.sleep(min(0.05, max(0.001, request["budget"]["wall_seconds"] - elapsed)))
            returncode = process.poll() if process is not None else None
            termination.phase = "worker_response_validation"
            if termination.signal_number is not None:
                status = "interrupted"
            if status is None:
                if time.monotonic() - start > request["budget"]["wall_seconds"]:
                    status = "timeout"
                elif returncode != 0:
                    status = "worker_error"
                elif not response_path.is_file():
                    status = "missing_response"
                else:
                    try:
                        if response_path.stat().st_size > 100_000_000:
                            raise ValueError("Response exceeds 100 MB")

                        def reject_constant(value):
                            raise ValueError(f"Nonfinite JSON value: {value}")

                        response = json.loads(response_path.read_text(), parse_constant=reject_constant)
                        if not isinstance(response, dict):
                            raise ValueError("Response must be an object")
                        candidates = response["candidates"]
                        if not isinstance(candidates, list):
                            raise ValueError("candidates must be a list")
                        evaluations = response.get("evaluations")
                        if (
                            isinstance(evaluations, bool)
                            or not isinstance(evaluations, int)
                            or evaluations < 0
                        ):
                            raise ValueError("Nonnegative integer evaluations accounting is required")
                        if (
                            not isinstance(response.get("evaluation_unit"), str)
                            or not response["evaluation_unit"]
                        ):
                            raise ValueError("evaluation_unit must identify counted work")
                        if (
                            len(candidates) > request["budget"]["max_candidates"]
                            or evaluations > request["budget"]["max_evaluations"]
                        ):
                            status = "budget_exceeded"
                        else:
                            status = "success" if candidates else "no_candidates"
                    except (KeyError, TypeError, ValueError) as exc:
                        status, error, response = "invalid_response", str(exc), {}
            # A fully written file is still partial work if the worker was
            # interrupted, or validation had not finished when signalled.
            if termination.signal_number is not None:
                status, response = "interrupted", {}
        finally:
            termination.phase = "worker_cleanup"
            if process is not None:
                cleanup = _stop_group(process)
    return {
        "status": status,
        "response": response,
        "error": error,
        "command": argv,
        "cleanup": cleanup,
        "resources": {
            "wall_seconds": time.monotonic() - start,
            "sampled_peak_rss_mb": peak,
            "sampled_cpu_seconds_lower_bound": cpu_max,
            "monitor_samples": monitor_samples,
            "gpu_memory": "not_measured",
            "includes_import_and_model_loading": True,
        },
    }


def evidence_verifier(assets):
    """Re-run structural route evidence; provider declarations cannot promote it."""
    catalog = rules = None
    if "catalog" in assets:
        from chembounce_synthesis.catalog import Catalog

        catalog = Catalog(assets["catalog"])
    if "rules" in assets:
        from chembounce_synthesis.reactions import Rules

        rules = Rules(assets["rules"])

    def verify(row):
        if row.get("route") and catalog is not None and rules is not None:
            from chembounce_synthesis.routes import verify_route

            try:
                result = verify_route(row["route"], catalog, rules)
                if (
                    chemical_identity(result["product"]).exact_smiles
                    != chemical_identity(row["smiles"]).exact_smiles
                ):
                    raise ValueError("Route product differs from candidate")
                return {
                    "level": "catalog_and_rule_verified",
                    "has_implicit_atom_sources": result["has_implicit_atom_sources"],
                    "experimental_validation": False,
                    "conditions_validated": False,
                    "catalog_sha256": catalog.sha256,
                    "rules_sha256": rules.sha256,
                }
            except (ValueError, KeyError, TypeError, RuntimeError) as exc:
                return {"level": "route_rejected", "reason": str(exc)}
        if catalog is not None and chemical_identity(row["smiles"]).exact_smiles in catalog.by_smiles:
            return {
                "level": "catalog_record",
                "catalog_sha256": catalog.sha256,
                "availability": "snapshot_record_not_live_stock",
            }
        return {"level": "unverified", "reason": "No independently replayable route supplied"}

    return verify


def run_benchmark(dataset_path, protocol_path, output):
    """Run a campaign on the main thread and restore the caller's handlers.

    SIGTERM/SIGINT checkpoint the active attempt, stop scheduling new attempts,
    then raise BenchmarkInterrupted. An already validated worker result remains
    valid if the signal arrives during evaluation or checkpointing. Cancellation
    during setup records a campaign with zero attempted runs.
    """
    output = Path(output).resolve()
    with _termination_scope() as termination:
        output.mkdir(parents=True, exist_ok=False)
        runs = _run_benchmark(dataset_path, protocol_path, output, termination)
        termination.phase = "campaign_checkpoint"
        if not runs:
            _write_checkpoint(output / "runs.json", [])
            _write_checkpoint(output / "summary.json", {})

        def checkpoint_campaign():
            _write_checkpoint(
                output / "campaign.json",
                {
                    "schema_version": 1,
                    "status": "interrupted" if termination.signal_number is not None else "completed",
                    "attempted_runs": len(runs),
                    "termination": termination.record(),
                },
            )

        checkpoint_campaign()
        if termination.signal_number is not None:
            # Also captures a first signal received while writing the final
            # campaign record; repeated signals do not disrupt finalization.
            checkpoint_campaign()
    if termination.signal_number is not None:
        raise BenchmarkInterrupted(termination.record(), runs)
    return runs


def _run_benchmark(dataset_path, protocol_path, output, termination):
    dataset, dataset_audit = load_dataset(dataset_path)
    if termination.signal_number is not None:
        return []
    protocol = load_protocol(protocol_path, dataset)
    if termination.signal_number is not None:
        return []
    root = Path(__file__).resolve().parents[2]
    code_paths = [
        *Path(__file__).parent.glob("*.py"),
        *root.glob("*.py"),
        *(root / "chembounce_synthesis").glob("*.py"),
    ]
    code = {str(p.relative_to(root)): digest(p) for p in code_paths}
    versions = {}
    for package in ("rdkit", "numpy", "scipy", "psutil"):
        versions[package] = importlib.metadata.version(package)
    manifest = {
        "schema_version": 1,
        "dataset": dataset_audit,
        "protocol_sha256": digest(protocol_path),
        "code_sha256": code,
        "python": sys.version,
        "platform": platform.platform(),
        "versions": versions,
        "protocol": protocol,
        "comparison_policy": "Identical outer wall/candidate caps; explicit per-method work caps beneath common maximum; evaluation units reported separately",
        "missing_methods": protocol.get("planned_unexecuted_methods", []),
        "ground_truth_policy": "Held-out labels/structures excluded from method requests; methods are trusted not to read dataset files",
    }
    write_json(output / "manifest.json", manifest)
    runs = []
    lookup = {r["id"]: r for r in dataset["compounds"]}
    for method in protocol["methods"]:
        termination.phase = "method_setup"
        if termination.signal_number is not None:
            return runs
        verifier = evidence_verifier(method["resolved_assets"])
        for query_index, query in enumerate(dataset["queries"]):
            for seed in protocol["seeds"]:
                termination.phase = "attempt_setup"
                if termination.signal_number is not None:
                    return runs
                run_dir = output / method["id"] / f"query-{query_index:04d}-seed-{seed}"
                request = public_request(dataset, query, method, seed, protocol["budget"], run_dir)
                if termination.signal_number is not None:
                    return runs
                report = bounded_method(method, request, run_dir, _termination=termination)
                termination.phase = "evaluation"
                start = time.monotonic()
                report["evaluation"] = evaluate_candidates(
                    dataset, query, report["response"].get("candidates", []), report["status"], verifier
                )
                report["evaluation_seconds"] = time.monotonic() - start
                report.update(
                    method_id=method["id"],
                    query_id=query["id"],
                    seed=seed,
                    target_id=lookup[query["compound_id"]]["target_id"],
                )
                report["response_sha256"] = (
                    digest(run_dir / "response.json") if (run_dir / "response.json").exists() else None
                )
                termination.phase = "attempt_checkpoint"
                _write_checkpoint(run_dir / "evaluation.json", report)
                runs.append(report)
                _write_checkpoint(output / "runs.json", runs)
                _write_checkpoint(output / "summary.json", aggregate_runs(runs))
                print(f"{method['id']} / {query['id']} / {seed}: {report['status']}", flush=True)
                if termination.signal_number is not None:
                    return runs
    return runs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        runs = run_benchmark(args.dataset, args.protocol, args.output)
    except BenchmarkInterrupted as exc:
        print(f"{exc}; checkpoint saved to {args.output}", file=sys.stderr, flush=True)
        return exc.exit_code
    return 0 if all(r["status"] in {"success", "no_candidates"} for r in runs) else 1


if __name__ == "__main__":
    raise SystemExit(main())
