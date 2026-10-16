"""Four-slot, saved-supply development ablation; no fresh benchmark authority.

A plan and separate launch authority must be frozen/reviewed before `run`.
Importing this module performs no chemistry, asset construction or model loading.
"""

import argparse
from contextlib import contextmanager
import signal
from collections import Counter
import copy
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import resource
import sys
import subprocess
import time

# Executed by absolute script path so worker PYTHONPATH can contain ONLY its arm.
if __package__:
    from .polar_capacity_supervisor import write, supervise, workload_check, identity, same
else:
    import polar_capacity_supervisor as _supervisor

    write, supervise, workload_check = _supervisor.write, _supervisor.supervise, _supervisor.workload_check
    identity, same = _supervisor.identity, _supervisor.same

DESIGN_SHA = "352a62807ae320a04d943721641061e12b86ce2804154fe7d07606c5cc7d5cde"
GATE_SHA = "6669792201d0d63bea5f6e7ab654a8152afa229da23e2211ea58a6ad6c727c24"
ARMS = {
    "baseline": "2443720dc74b740a9a1761fc30b53e28b667b3f4",
    "capacity": "2556b14db3e7c239352e37168c857ccb185cc5a1",
}
ORDER = [
    ("attempt-0000", "baseline"),
    ("attempt-0000", "capacity"),
    ("attempt-0009", "capacity"),
    ("attempt-0009", "baseline"),
]
PROFILE = "qualified-selective-engineering-v2"
PROFILE_SHA = "a8b47a409a21518dd482acbefdce38774414ae0b65f28d6ea7cb52ab63e1ef52"
PACKAGES = ("rdkit", "numpy", "scipy", "torch", "transformers", "tokenizers", "safetensors", "psutil")
THREADS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")


def read(path):
    def invalid(value):
        raise ValueError("Nonfinite JSON: " + value)

    return json.loads(Path(path).read_text(), parse_constant=invalid)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def binding(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


class ArchiveCoverageError(ValueError):
    """The fixed historical supply has no response for a requested design."""


class IntegrityError(ValueError):
    """Systemic source/input authority failure; abort the remaining plan."""


def verify_pins(pins):
    if not pins:
        raise ValueError("Unresolved pins")
    for filename, expected in pins.items():
        path = Path(filename)
        if not path.is_absolute() or str(path.resolve()) != filename or digest(path) != expected:
            raise IntegrityError("File binding differs: " + filename)


def tracked_inventory(source):
    """Complete committed checkout authority, including root imports and tests."""
    root = Path(source["path"]).resolve()

    def git(*args):
        return subprocess.check_output(["git", *args], cwd=root, timeout=10).decode().strip()

    if git("rev-parse", "HEAD") != source["commit"] or git("rev-parse", "HEAD^{tree}") != source.get("tree"):
        raise IntegrityError("Source commit/tree differs")
    if git("status", "--porcelain", "--untracked-files=all"):
        raise IntegrityError("Source checkout is not clean")
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=root, timeout=10).decode().split("\0")
    paths = {str((root / name).resolve()) for name in tracked if name}
    if paths != set(source["files"]):
        raise IntegrityError("Complete tracked source inventory is required")
    verify_pins(source["files"])
    return paths


def required_runtime_files():
    """Distribution records from the current context, no scientific imports."""
    paths = {str(Path(sys.executable).resolve())}
    for name in PACKAGES:
        dist = importlib.metadata.distribution(name)
        if not dist.files:
            raise IntegrityError("Missing runtime distribution inventory: " + name)
        for relative in dist.files:
            path = Path(dist.locate_file(relative)).resolve()
            if path.suffix != ".pyc" and path.is_file():
                paths.add(str(path))
    generation = Path(
        importlib.metadata.distribution("transformers").locate_file("transformers/generation/utils.py")
    ).resolve()
    if not generation.is_file():
        raise IntegrityError("Missing physical generation source")
    paths.add(str(generation))
    return paths


@contextmanager
def outer_watchdog(deadline):
    def expired(number, frame):
        raise TimeoutError("Fixed 4680-second outer deadline exhausted")

    previous = signal.signal(signal.SIGALRM, expired)
    previous_timer = signal.setitimer(signal.ITIMER_REAL, max(0.001, deadline - time.monotonic()))
    if previous_timer != (0.0, 0.0):
        signal.setitimer(signal.ITIMER_REAL, *previous_timer)
        signal.signal(signal.SIGALRM, previous)
        raise RuntimeError("Existing timer conflicts with owned outer watchdog")
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def runtime():
    return {
        "python": platform.python_version(),
        "executable": str(Path(sys.executable).resolve()),
        "packages": {name: importlib.metadata.version(name) for name in PACKAGES},
        "cwd": str(Path.cwd().resolve()),
        "PYTHONPATH": os.environ.get("PYTHONPATH"),
    }


def environment(source):
    return {
        **os.environ,
        **dict.fromkeys(THREADS, "1"),
        "CUDA_VISIBLE_DEVICES": "",
        "HF_HUB_OFFLINE": "1",
        "TRANSFORMERS_OFFLINE": "1",
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONPATH": str(Path(source).resolve()),
    }


def cpu_guard():
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "" or any(os.environ.get(k) != "1" for k in THREADS):
        raise ValueError("CPU-only one-thread environment required")
    if os.environ.get("HF_HUB_OFFLINE") != "1" or os.environ.get("TRANSFORMERS_OFFLINE") != "1":
        raise ValueError("Offline environment required")

    def audit(event, args):
        if event in {"socket.connect", "socket.getaddrinfo", "socket.bind"}:
            # Denied even if an import catches its own capability-probe failure.
            raise RuntimeError("Development network operation forbidden: " + event)

    sys.addaudithook(audit)


def validate_config(config, archive):
    """Pure fixed policy checks; `_config` resolution is checked in each context."""
    expected = {
        "schema_version": 1,
        "recommendation_profile": PROFILE,
        "mode": "hop",
        "minimum_scores": {"hit_tanimoto": 0.35, "hit_shape": 0.70, "hit_pharmacophore": 0.50},
        "similarity_constraints": {
            "version": "selective-similarity-constraints-v2",
            "polar_recall": 0.75,
            "landmarks": [],
        },
        "ranking_policy": "original-hit-v1",
        "route_search_policy": "original-hit-v1",
        "route_search_cache_size": 4096,
        "allocation_policy": "balanced",
        "max_design_attempts": 12,
        "max_designs": 3,
        "max_scored_products": 12,
        "max_results": 3,
        "provider": {"kind": "saved", "path": archive},
        "conformers": {
            "scoring_objective": "original-hit-v1",
            "num_conformers": 8,
            "seed": 42,
            "energy_window": 5.0,
            "strain_penalty": 0.05,
            "max_iterations": 500,
            "pharmacophore_sigma": 1.0,
            "reject_unspecified_stereo": True,
        },
        "search": {
            "beam_width": 6,
            "building_block_top_k": 3,
            "max_combinations": 256,
            "max_steps": 12,
            "reaction_product_limit": 64,
        },
        "expansion": {"building_block_top_k": 6, "max_routes": 12, "min_tanimoto": 0.0, "rounds": 1},
        "refinement": {
            "top_n": 6,
            "num_conformers": 8,
            "seeds": [42, 43, 44],
            "threshold_policy": "all_seeds",
            "ranking_statistic": "minimum",
            "max_total_conformers": 144,
            "max_score_range": None,
            "max_metric_range": None,
            "adaptive_max_conformers": None,
            "conformer_growth_factor": 2,
            "borderline_margin": 0.0,
        },
    }
    if any(config.get(k) != v for k, v in expected.items()):
        raise ValueError("Fixed development settings differ")
    if (
        "reference_pose" in config
        or "pocket" in config
        or config.get("reference_sdf") is not None
        or config.get("core_smiles") is not None
        or config.get("carry_forward", [])
        or config.get("preserve_design_scaffold", False)
    ):
        raise ValueError("Unapproved reference/core/continuation settings")
    evidence = config["synthesis_evidence"]
    if evidence.get("profile") != "computational-route-support-v1" or evidence.get("policy") != {
        "require_forward_agreement": True
    }:
        raise ValueError("Strict route evidence required")
    predictor = evidence["predictor"]
    required = {
        "kind": "reactiont5",
        "device": "cpu",
        "top_k": 3,
        "max_new_tokens": 256,
        "max_input_tokens": 512,
        "model_revision": "933114058cb2604dc1bf536dbebdfcefbe83d4fc",
    }
    if any(
        evidence.get(k) != v
        for k, v in {"max_predictor_calls": 24, "max_input_tokens": 12288, "max_output_tokens": 18432}.items()
    ):
        raise ValueError("Aggregate evidence budgets differ")
    if any(predictor.get(k) != v for k, v in required.items()):
        raise ValueError("Forward budget/settings differ")


def validate_archive(query):
    pool = read(query["design_pool"])["designs"]
    records = read(query["proposals"])["proposals"]
    if not isinstance(pool, list) or len(pool) != 12 or not all(isinstance(x, str) for x in pool):
        raise ValueError("All 12 literal ordered designs required")
    indices = [1, 1, 2, 2, 3, 3]
    if (
        len(records) != 6
        or [r["design"] for r in records] != [pool[i] for i in indices]
        or [r.get("sampling", {}).get("seed") for r in records] != [42, 43, 42, 43, 42, 43]
    ):
        raise ValueError("Archived response coverage/order/seeds differ")
    original = read(query["request"])["query"]
    if original != query["query"]:
        raise ValueError("Query differs from bound historical public request")
    return pool


def validate_plan(plan, *, check_files=True):
    if set(plan) != {
        "schema",
        "design",
        "gate",
        "pins",
        "sources",
        "checker",
        "queries",
        "observer",
        "python",
        "output",
        "workloads",
        "administrative_runtime",
        "feature_definition",
    }:
        raise ValueError("Unknown or missing plan fields")
    if plan["schema"] != "polar-capacity-development-plan-v2":
        raise ValueError("Unknown development plan")
    if plan["design"]["sha256"] != DESIGN_SHA or plan["gate"]["sha256"] != GATE_SHA:
        raise ValueError("Unapproved design gate")
    if set(plan["sources"]) != set(ARMS) or set(plan["queries"]) != {"attempt-0000", "attempt-0009"}:
        raise ValueError("Fixed source pair/two-input slice required")
    required = {plan["design"]["path"]: DESIGN_SHA, plan["gate"]["path"]: GATE_SHA}
    required[plan["feature_definition"]["path"]] = plan["feature_definition"]["sha256"]
    if check_files:
        verify_pins(required)
        design = read(plan["design"]["path"])
        for slot in design["development_bindings"]["slots"]:
            query = plan["queries"][slot["historical_slot"]]
            for name, suffix in [
                ("design_pool", "/native/design-pool.json"),
                ("proposals", "/common/synthesis/proposals.json"),
                ("request", "/common/request.json"),
            ]:
                path = next(p for p in slot["bindings"] if p.endswith(suffix))
                if query[name] != path or plan["pins"].get(path) != slot["bindings"][path]:
                    raise IntegrityError("Input differs from reviewed development slice")
            model = slot["model_metadata"]
            if query["config"]["synthesis_evidence"]["predictor"]["model_path"] != model["path"]:
                raise IntegrityError("Reviewed model path differs")
            for name, sha in model["model_files_sha256"].items():
                if plan["pins"].get(str(Path(model["path"]) / name)) != sha:
                    raise IntegrityError("Reviewed model metadata/weights binding differs")
        for asset in ("catalog", "rules"):
            bound = design["assets"][asset]
            if plan["pins"].get(bound["path"]) != bound["recorded_sha256"]:
                raise IntegrityError("Reviewed asset binding differs")
            if any(q["config"][asset] != bound["path"] for q in plan["queries"].values()):
                raise IntegrityError("Reviewed asset path differs")
        sidecar = design["assets"]["catalog_sidecar"]
        if plan["pins"].get(sidecar["path"]) != sidecar["sha256"]:
            raise IntegrityError("Reviewed sidecar differs")
    harness = Path(__file__).resolve().parent
    for filename in (
        "polar_capacity_development.py",
        "polar_capacity_observer.py",
        "polar_capacity_supervisor.py",
        "qualified_collector.py",
        "qualified_adapter.py",
    ):
        path = str(harness / filename)
        if plan["checker"]["files"].get(path) != plan["pins"].get(path) or path not in plan["pins"]:
            raise ValueError("Missing executable tooling pin")
    if str(Path(plan["python"]).resolve()) not in plan["pins"]:
        raise ValueError("Missing interpreter pin")
    for role, source in list(plan["sources"].items()) + [("checker", plan["checker"])]:
        root = Path(source["path"])
        if (
            not root.is_absolute()
            or not source.get("commit")
            or not source.get("runtime")
            or not source.get("files")
            or not source.get("synthesis_source")
        ):
            raise ValueError("Unresolved source/runtime authority")
        if role in ARMS and source["commit"] != ARMS[role]:
            raise ValueError("Source arm changed")
        required.update(source["files"])
        package = root / "chembounce_synthesis"
        if check_files:
            tracked_inventory(source)
            actual = {p.name: digest(p) for p in package.glob("*.py")}
            if source["synthesis_source"] != actual:
                raise IntegrityError("Incomplete package source map")
        for name, value in source["synthesis_source"].items():
            if source["files"].get(str(package / name)) != value:
                raise ValueError("Missing package pin")
        legacy = str(package / "data" / "qualified-development-v1.json")
        if source["files"].get(legacy) != "9d3c60a895912b28c13a29aac0a04f385ddb8a0effc9fcd81d64a50a3d48bdd2":
            raise ValueError("Missing unchanged legacy artifact pin")
        artifact = str(package / "data" / (PROFILE + ".json"))
        if source["files"].get(artifact) != PROFILE_SHA:
            raise ValueError("Missing engineering artifact pin")
    for name in (
        "common.py",
        "identity.py",
        "conformers.py",
        "selective_scoring.py",
        "similarity_constraints.py",
    ):
        values = [
            source["synthesis_source"].get(name) for source in [*plan["sources"].values(), plan["checker"]]
        ]
        if None in values or len(set(values)) != 1:
            raise ValueError("Protected numerical kernel differs")
    for query in plan["queries"].values():
        for field in ("design_pool", "proposals", "request"):
            if query[field] not in plan["pins"]:
                raise ValueError("Missing original input pin")
        validate_config(query["config"], query["proposals"])
        for field in ("catalog", "rules"):
            if query["config"][field] not in plan["pins"]:
                raise ValueError("Missing chemistry asset binding")
        sidecar = query["config"]["catalog"] + ".manifest.json"
        if sidecar not in plan["pins"]:
            raise ValueError("Explicit present catalog sidecar binding required")
        model = Path(query["config"]["synthesis_evidence"]["predictor"]["model_path"])
        if check_files and any(str(p.resolve()) not in plan["pins"] for p in model.iterdir() if p.is_file()):
            raise ValueError("Incomplete local forward bundle binding")
        if check_files:
            validate_archive(query)
    if set(plan["observer"]) != set(ARMS) or any(not v for v in plan["observer"].values()):
        raise ValueError("Complete per-arm observer inventory required")
    for entries in plan["observer"].values():
        for entry in entries:
            required[entry["source_path"]] = entry["source_sha256"]
            for implementation in entry.get("implementation", []):
                if implementation["implementation_module_path"]:
                    required[implementation["implementation_module_path"]] = implementation[
                        "implementation_module_sha256"
                    ]
    if any(plan["pins"].get(p) != value for p, value in required.items()):
        raise ValueError("Authority pins missing/conflicting")
    if not plan["workloads"].get("identities") or not plan["workloads"].get("command_markers"):
        raise ValueError("Root-declared known workload identity/command coverage required")
    if not Path(plan["output"]).is_absolute() or not Path(plan["python"]).is_absolute():
        raise ValueError("Absolute execution paths required")
    if check_files:
        if not required_runtime_files() <= plan["pins"].keys():
            raise IntegrityError("Incomplete scientific/runtime/generation/interpreter pins")
        verify_pins(plan["pins"])
    return plan


def load_plan(path, expected):
    if digest(path) != expected:
        raise ValueError("Plan digest mismatch")
    return validate_plan(read(path))


def context_check(source, tooling=None):
    actual_runtime = runtime()
    if actual_runtime != source["runtime"]:
        raise IntegrityError("Actual cwd/runtime metadata differs from frozen context")
    required = {
        "rdkit": "2025.9.6",
        "numpy": "1.26.4",
        "torch": "2.8.0",
        "transformers": "4.51.3",
        "tokenizers": "0.21.4",
        "safetensors": "0.6.2",
    }
    if actual_runtime.get("python") != "3.10.12" or any(
        actual_runtime.get("packages", {}).get(k) != v for k, v in required.items()
    ):
        raise IntegrityError("Declared development numerical runtime differs")
    for name, module in list(sys.modules.items()):
        if (
            name in {"chembounce", "utils", "chembounce_synthesis"}
            or name.startswith("chembounce_synthesis.")
            or name.startswith("benchmarks.research.")
            or name.startswith("polar_capacity_")
        ):
            path = Path(module.__file__).resolve()
            owner = (
                tooling
                if tooling
                and (
                    name.startswith("polar_capacity_")
                    or name.startswith("benchmarks.research.polar_capacity_")
                )
                else source
            )
            owner_root = Path(owner["path"]).resolve()
            if not path.is_relative_to(owner_root) or owner["files"].get(str(path)) != digest(path):
                raise IntegrityError("Imported module escaped pinned authority: " + name)


def request_for(plan, plan_sha, index):
    historical, arm = ORDER[index]
    return {
        "schema": "polar-capacity-development-request-v1",
        "plan_sha256": plan_sha,
        "index": index,
        "historical_slot": historical,
        "arm": arm,
        "query": plan["queries"][historical]["query"],
        "config": plan["queries"][historical]["config"],
        "source_commit": plan["sources"][arm]["commit"],
        "source_sha256": plan["sources"][arm]["synthesis_source"],
        "archive": {
            key: {
                "path": plan["queries"][historical][key],
                "sha256": plan["pins"][plan["queries"][historical][key]],
            }
            for key in ("design_pool", "proposals", "request")
        },
        "observer_sha256": binding(plan["observer"][arm]),
    }


def authenticate_request(request, plan, plan_sha):
    index = request.get("index")
    if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < 4:
        raise ValueError("Invalid planned slot")
    if request != request_for(plan, plan_sha, index):
        raise ValueError("Development request differs from frozen authority")
    return index


def inventory(directory):
    return {
        str(p.relative_to(directory)): digest(p)
        for p in sorted(Path(directory).rglob("*"))
        if p.is_file() and p.name not in {"stdout.log", "stderr.log", "receipt.json", "terminal.json"}
    }


def receipt_eligible(receipt, allocation):
    """Recompute the allocation/cleanup predicate; labels are insufficient."""

    def finite(value):
        return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)

    return (
        receipt.get("status") == "completed"
        and receipt.get("returncode") == 0
        and receipt.get("allocation_seconds") == allocation
        and finite(receipt.get("wall_seconds"))
        and 0 <= receipt["wall_seconds"] <= allocation
        and finite(receipt.get("parent_observed_phase_seconds"))
        and receipt["wall_seconds"] <= receipt["parent_observed_phase_seconds"] <= allocation
        and receipt.get("cleanup", {}).get("complete") is True
        and receipt["cleanup"].get("remaining") == []
        and receipt.get("identity") is not None
        and isinstance(receipt.get("terminal_sha256"), str)
        and len(receipt["terminal_sha256"]) == 64
        and receipt.get("resources", {}).get("samples", 0) > 0
        and receipt["resources"].get("gaps") == []
        and receipt["resources"].get("sampled_tree_peak_rss_bytes") is not None
    )


def launch_authority(args):
    if not args.authorization or digest(args.authorization) != args.authorization_sha256:
        raise ValueError("Separate frozen launch authorization required")
    authorization = read(args.authorization)
    if (
        authorization.get("schema") != "polar-capacity-development-launch-v1"
        or authorization.get("plan_sha256") != args.plan_sha256
        or authorization.get("source_review_passed") is not True
        or authorization.get("resource_clearance") is not True
    ):
        raise ValueError("Unapproved source/resource execution authority")
    return authorization


def child_authority(args, plan):
    launch_authority(args)
    receipt = read(Path(args.output) / "receipt.json")
    parent = identity(os.getppid())
    if (
        receipt.get("status") != "running"
        or not same(parent, receipt.get("supervisor_identity"))
        or not same(identity(os.getpid()), receipt.get("identity"))
        or args.deadline != receipt["started_monotonic"] + receipt["allocation_seconds"]
        or receipt.get("argv", [])[1:] != [str(Path(__file__).resolve()), *sys.argv[1:]]
        or workload_check(plan["workloads"])
    ):
        raise ValueError("Scientific child requires its live authenticated supervisor/launch context")
    if args.command in {"worker", "assess"} and receipt["allocation_seconds"] != (
        600 if args.command == "worker" else 300
    ):
        raise ValueError("Scientific child allocation differs")


def worker(plan, plan_sha, request, output, deadline):
    index = authenticate_request(request, plan, plan_sha)
    query_id, arm = ORDER[index]
    source = plan["sources"][arm]
    context_check(source, plan["checker"])
    from chembounce_synthesis import workflow
    from chembounce_synthesis.forward_support import binding as scientific_binding

    context_check(source, plan["checker"])
    from rdkit import RDConfig

    feature_path = str((Path(RDConfig.RDDataDir) / "BaseFeatures.fdef").resolve())
    if (
        feature_path != plan["feature_definition"]["path"]
        or digest(feature_path) != plan["feature_definition"]["sha256"]
    ):
        raise IntegrityError("Numerical feature definition differs")
    config = workflow._resolve_paths(workflow._config(copy.deepcopy(request["config"])), Path.cwd())
    if scientific_binding(config) != binding(request["config"]):
        raise ValueError("Plan must contain fully resolved effective config, not partial defaults")
    if workflow.effective_policy(config)["calibration"]["status"] != "engineering_policy_matched":
        raise IntegrityError("Engineering numerical policy/runtime is not matched")
    if __package__:
        from .polar_capacity_observer import Observer, install_scientific
    else:
        from polar_capacity_observer import Observer, install_scientific
    observer = Observer(output / "events.jsonl", expected_inventory=plan["observer"][arm])
    try:
        install_scientific(observer)
        context_check(source, plan["checker"])
        observer.phase = "workflow"
        try:
            workflow.run(
                request["query"]["smiles"],
                validate_archive(plan["queries"][query_id]),
                config,
                output / "synthesis",
                deadline_monotonic=deadline,
            )
        except ValueError as exc:
            if str(exc).startswith("No saved proposals for design "):
                raise ArchiveCoverageError(str(exc)) from exc
            raise
    finally:
        observer.close()
    context_check(source, plan["checker"])
    verify_pins(plan["pins"])
    summary = read(output / "synthesis/summary.json")
    result = read(output / "synthesis/results.json")
    coverage_failures = [f for f in summary.get("failures", []) if "No saved proposals for design" in str(f)]
    return {
        "schema": "polar-capacity-development-worker-v1",
        "request_sha256": binding(request),
        "index": index,
        "status": "returned",
        "coverage_failures": coverage_failures,
        "files_sha256": inventory(output),
        "candidates": [{"smiles": r["smiles"], "route_id": r["route_id"]} for r in result["results"]],
        "runtime": runtime(),
        "source_sha256": source["synthesis_source"],
        "observer_sha256": binding(plan["observer"][arm]),
    }


def forbid_assessor_compute():
    """Install explicit denials before importing tokenizer/scientific checker."""

    def denied(*args, **kwargs):
        raise RuntimeError("Assessor forbids model construction/generation/embedding/alignment")

    import torch
    from transformers.generation.utils import GenerationMixin
    from rdkit.Chem import AllChem, rdMolAlign, rdDistGeom, rdShapeAlign

    torch.nn.Module.__init__ = denied
    torch.load = denied
    GenerationMixin.generate = denied
    import safetensors.torch
    import safetensors

    safetensors.safe_open = denied

    safetensors.torch.load_file = denied
    safetensors.torch.load = denied
    safetensors.torch.safe_open = denied
    AllChem.EmbedMolecule = AllChem.EmbedMultipleConfs = denied
    rdDistGeom.EmbedMolecule = rdDistGeom.EmbedMultipleConfs = denied
    rdShapeAlign.AlignMol = denied
    AllChem.AlignMol = AllChem.GetO3A = AllChem.GetCrippenO3A = denied
    rdMolAlign.AlignMol = rdMolAlign.GetO3A = rdMolAlign.GetCrippenO3A = denied
    from chembounce_synthesis import conformers, selective_scoring, workflow

    conformers.ensemble = selective_scoring.ensemble = denied
    conformers._overlay = selective_scoring._overlay = denied
    conformers.score_candidate = selective_scoring.score_candidate = denied
    workflow.run = denied
    from chembounce_synthesis.forward_support import PinnedReactionT5Forward

    PinnedReactionT5Forward.__init__ = denied
    PinnedReactionT5Forward.predict_raw = denied


def assess(plan, plan_sha, request, worker_dir, worker_receipt, deadline):
    index = authenticate_request(request, plan, plan_sha)
    context_check(plan["checker"])
    if not receipt_eligible(worker_receipt, 600):
        raise ValueError("Only completed in-budget worker evidence can certify outputs")
    if digest(worker_dir / "terminal.json") != worker_receipt.get("terminal_sha256"):
        raise ValueError("Worker terminal differs from parent receipt")
    terminal = read(worker_dir / "terminal.json")
    if (
        terminal.get("status") != "returned"
        or terminal.get("request_sha256") != binding(request)
        or terminal.get("source_sha256") != request["source_sha256"]
        or terminal.get("observer_sha256") != request["observer_sha256"]
        or terminal.get("index") != index
        or terminal["files_sha256"] != inventory(worker_dir)
    ):
        raise ValueError("Worker terminal/output provenance differs")
    events = [json.loads(line) for line in (worker_dir / "events.jsonl").read_text().splitlines()]
    event_summary(events, expected_inventory=plan["observer"][request["arm"]])
    if terminal.get("coverage_failures"):
        raise ValueError("Archive coverage failed; no validated outcome")
    forbid_assessor_compute()
    from benchmarks.research.qualified_collector import checked_workflow_artifacts
    from chembounce_synthesis.workflow import _config

    context_check(plan["checker"])
    expected = _config(copy.deepcopy(request["config"]))
    if expected != request["config"]:
        raise ValueError("Development config is not fully resolved")
    response = {"metadata": {"files_sha256": terminal["files_sha256"]}, "candidates": terminal["candidates"]}
    products = checked_workflow_artifacts(
        worker_dir,
        response,
        request["query"],
        request["source_sha256"],
        expected_config=expected,
        profile_name=PROFILE,
        deadline=deadline,
    )
    saved_rows = read(worker_dir / "synthesis/results.json")["results"]
    products = [
        {**product, "route_id": row["route_id"], "route_support_binding": row["route_support_binding"]}
        for product, row in zip(products, saved_rows)
    ]
    context_check(plan["checker"])
    verify_pins(plan["pins"])
    return {
        "schema": "polar-capacity-development-attestation-v1",
        "index": index,
        "plan_sha256": plan_sha,
        "request_sha256": binding(request),
        "status": "validated",
        "worker_terminal_sha256": digest(worker_dir / "terminal.json"),
        "worker_receipt_sha256": digest(worker_dir / "receipt.json"),
        "checker_source_sha256": plan["checker"]["files"],
        "runtime": runtime(),
        "qualified": products,
        "coverage_failures": terminal["coverage_failures"],
    }


def summarize(rows):
    if len(rows) != 4 or [r.get("index") for r in rows] != list(range(4)):
        raise ValueError("All four ordered planned slots required")
    detail, totals = [], {arm: {"positive": 0, "families": 0, "at_least_three": 0} for arm in ARMS}
    for index, row in enumerate(rows):
        historical, arm = ORDER[index]
        valid = (
            receipt_eligible(row.get("worker", {}), 600)
            and receipt_eligible(row.get("assessor", {}), 300)
            and row.get("attestation", {}).get("status") == "validated"
            and row["attestation"].get("index") == index
            and row["attestation"].get("coverage_failures") == []
        )
        products = row["attestation"]["qualified"] if valid else []
        if len(products) > 3 or len({p["family"] for p in products}) != len(products):
            raise ValueError("Validated top-3 families are not unique/bounded")
        n = len(products)
        totals[arm]["positive"] += bool(n)
        totals[arm]["families"] += n
        totals[arm]["at_least_three"] += n >= 3
        detail.append(
            {
                "index": index,
                "historical_slot": historical,
                "arm": arm,
                "primary": int(n > 0),
                "qualified_families": n,
                "status": row.get("status", "not_started"),
                "identities": [p["smiles"] for p in products],
                "routes": [
                    {k: p[k] for k in ("smiles", "route_id", "route_support_binding")} for p in products
                ],
                "diversity": {
                    "normalized_scaffolds": len({p["scaffold"] for p in products}),
                    "generic_scaffolds": len({p["generic_scaffold"] for p in products}),
                }
                if n >= 2
                else None,
                "diversity_reason": None if n >= 2 else "insufficient_qualified_families",
                "capacity_counters": None if arm == "baseline" else row.get("capacity_counters"),
                "capacity_scope": "not_applicable" if arm == "baseline" else "recorded_only",
            }
        )
    result = {
        arm: {
            "primary": {"numerator": t["positive"], "denominator": 2},
            "family_yield": {"numerator": t["families"], "denominator": 6},
            "at_least_three": {"numerator": t["at_least_three"], "denominator": 2},
        }
        for arm, t in totals.items()
    }
    pairs = []
    for historical in ("attempt-0000", "attempt-0009"):
        pair = {r["arm"]: r for r in detail if r["historical_slot"] == historical}
        before, after = set(pair["baseline"]["identities"]), set(pair["capacity"]["identities"])
        pairs.append(
            {
                "historical_slot": historical,
                "primary_difference": pair["capacity"]["primary"] - pair["baseline"]["primary"],
                "new": sorted(after - before),
                "lost": sorted(before - after),
                "shared": sorted(before & after),
                "baseline_routes": pair["baseline"]["routes"],
                "capacity_routes": pair["capacity"]["routes"],
            }
        )
    return {
        "arms": result,
        "rows": detail,
        "pairs": pairs,
        "primary_difference": {
            "numerator": totals["capacity"]["positive"] - totals["baseline"]["positive"],
            "denominator": 2,
        },
        "scope": "Two-pair conditional downstream development only; no inferential or causal timing claim.",
    }


def call_difference(before, after):
    """Multiplicity-aware physical started/returned work; never helper counts."""

    def counts(events):
        return Counter(
            (e["kind"], e["call_key"])
            for e in events
            if e.get("physical") and e["event"] == "started" and e.get("call_key") is not None
        )

    a, b = counts(before), counts(after)
    return {
        name: [{"kind": k[0], "call_key": k[1], "count": n} for k, n in sorted(value.items())]
        for name, value in [("removed_started", a - b), ("new_started", b - a), ("unchanged_started", a & b)]
    }


def event_summary(events, *, expected_inventory=None):
    installed = [e for e in events if e.get("event") == "installed"]
    closed = [e for e in events if e.get("event") == "closed"]
    if expected_inventory is not None and (
        len(installed) != 1
        or installed[0]["inventory"] != expected_inventory
        or installed[0]["inventory_sha256"] != binding(expected_inventory)
        or len(closed) != 1
        or closed[0]["unexpected_replacements"]
    ):
        raise ValueError("Actual observer inventory/completion differs")
    started, finished = {}, {}
    for event in events:
        if event.get("event") == "started":
            if event["id"] in started:
                raise ValueError("Duplicate observed call start")
            started[event["id"]] = event
        elif event.get("event") in {"returned", "failed"}:
            ident = event["id"]
            if ident not in started or ident in finished or event["call_key"] != started[ident]["call_key"]:
                raise ValueError("Unbound observed call completion")
            finished[ident] = event
    counts = {}
    for ident, event in started.items():
        kind = event["kind"]
        row = counts.setdefault(kind, {"started": 0, "returned": 0, "failed": 0, "nonreturned": 0})
        row["started"] += 1
        row[finished[ident]["event"] if ident in finished else "nonreturned"] += 1
    if expected_inventory is not None and any(i not in finished for i in started):
        raise ValueError("Completed worker has unfinished observer call")
    return {"calls": counts, "scope": "Logical and physical layers separate; nonreturned costs unknown"}


def final_report(ledger, output):
    report = summarize(ledger["rows"])
    events_by_slot = {}
    for row in ledger["rows"]:
        directory = output / f"slot-{row['index']:04d}" / "worker"
        event_file = directory / "events.jsonl"
        if event_file.exists():
            events = []
            partial = False
            for line in event_file.read_text().splitlines():
                try:
                    events.append(json.loads(line))
                except ValueError:
                    partial = True
                    break
            events_by_slot[row["index"]] = events
            report["rows"][row["index"]]["observed_work"] = event_summary(events)
            report["rows"][row["index"]]["observer_partial"] = partial
        else:
            report["rows"][row["index"]]["observed_work"] = None
        summary_path = directory / "synthesis/summary.json"
        if summary_path.exists():
            actual = read(summary_path)
            report["rows"][row["index"]]["capacity_counters"] = (
                actual.get("polar_capacity") if ORDER[row["index"]][1] == "capacity" else None
            )
            report["rows"][row["index"]]["synthesis_evidence_usage"] = actual.get("synthesis_evidence_usage")
        report["rows"][row["index"]]["costs"] = {
            phase: row.get(phase, {}).get("resources") for phase in ("worker", "assessor")
        }
    report["physical_call_differences"] = []
    for before, after in ((0, 1), (3, 2)):
        report["physical_call_differences"].append(
            {
                "baseline_index": before,
                "capacity_index": after,
                "comparison": call_difference(events_by_slot[before], events_by_slot[after])
                if before in events_by_slot and after in events_by_slot
                else None,
            }
        )
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare", "run", "worker", "assess", "admin"])
    parser.add_argument("--plan", required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--request")
    parser.add_argument("--worker-dir")
    parser.add_argument("--deadline", type=float)
    parser.add_argument("--authorization")
    parser.add_argument("--authorization-sha256")
    args = parser.parse_args(argv)
    if args.command == "run":
        return run(args)
    output = Path(args.output)
    if args.command == "prepare" and output.exists():
        raise FileExistsError(output)
    if args.command != "prepare":
        if (output / "terminal.json").exists():
            raise FileExistsError("Terminal artifact already exists")
        receipt = read(output / "receipt.json")
        if receipt.get("status") != "running" or not same(
            identity(os.getppid()), receipt.get("supervisor_identity")
        ):
            raise ValueError("No active owned supervisor; refusing output mutation")
    if args.command != "prepare":
        launch_authority(args)
        if digest(args.plan) != args.plan_sha256:
            raise ValueError("Plan digest mismatch before child preparation")
        declared = read(args.plan)
        if workload_check(declared["workloads"]):
            raise ValueError("Known workload blocks even child preparation")
    cpu_guard()
    terminal = {"status": "failed"}
    try:
        plan = load_plan(args.plan, args.plan_sha256)
        if args.command != "prepare":
            child_authority(args, plan)
        if args.command == "prepare":
            output.mkdir(parents=True, exist_ok=False)
            write(
                output / "ledger.json",
                {
                    "plan_sha256": args.plan_sha256,
                    "rows": [
                        {
                            "index": i,
                            "status": "not_started",
                            "request": request_for(plan, args.plan_sha256, i),
                        }
                        for i in range(4)
                    ],
                },
            )
            terminal = {"status": "prepared_not_authorized", "plan_sha256": args.plan_sha256}
        elif args.command == "admin":
            terminal = {"status": "verified", "plan_sha256": args.plan_sha256}
            context_check(plan["checker"])
            if args.request:
                authority = read(args.request)
                if digest(authority["ledger_path"]) != authority["ledger_sha256"]:
                    raise IntegrityError("Final ledger changed")
                terminal["report"] = final_report(read(authority["ledger_path"]), Path(plan["output"]))
        else:
            if args.deadline is None or time.monotonic() >= args.deadline:
                raise TimeoutError("Missing/expired parent deadline")
            request = read(args.request)
            if args.command == "worker":
                terminal = worker(plan, args.plan_sha256, request, output, args.deadline)
            else:
                wd = Path(args.worker_dir)
                terminal = assess(
                    plan, args.plan_sha256, request, wd, read(wd / "receipt.json"), args.deadline
                )
    except BaseException as exc:
        terminal = {
            "status": "failed",
            "failure_category": "integrity"
            if isinstance(exc, IntegrityError)
            else "archive_coverage_failure"
            if isinstance(exc, ArchiveCoverageError)
            else "execution",
            "error": {"type": type(exc).__name__, "message": str(exc)},
        }
    terminal["process_lifetime_ru_maxrss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if output.exists():
        write(output / "terminal.json", terminal)
    return 1 if terminal["status"] == "failed" else 0


def run(args):
    """Parent: only supervision, hashes, requests and arithmetic; no science."""
    started = time.monotonic()
    with outer_watchdog(started + 4680):
        return _run(args, started, started + 4680)


def _run(args, outer_started, outer_deadline):
    setup_deadline = outer_started + 300
    launch_authority(args)
    if digest(args.plan) != args.plan_sha256:
        raise ValueError("Plan changed")
    plan = validate_plan(read(args.plan), check_files=False)
    if str(Path(args.output).resolve()) != plan["output"]:
        raise ValueError("Unexpected output path")
    if runtime() != plan["administrative_runtime"]:
        raise ValueError("Parent runtime/cwd differs")
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    rows = [
        {"index": i, "status": "not_started", "request": request_for(plan, args.plan_sha256, i)}
        for i in range(4)
    ]
    ledger = {
        "schema": "polar-capacity-development-ledger-v1",
        "plan_sha256": args.plan_sha256,
        "status": "running",
        "outer_started_monotonic": outer_started,
        "outer_deadline_monotonic": outer_deadline,
        "rows": rows,
        "administration": {},
    }
    write(output / "ledger.json", ledger)
    script = str(Path(__file__).resolve())
    base = [plan["python"], script]
    shared = [
        "--plan",
        str(Path(args.plan).resolve()),
        "--plan-sha256",
        args.plan_sha256,
        "--authorization",
        str(Path(args.authorization).resolve()),
        "--authorization-sha256",
        args.authorization_sha256,
    ]
    blocked = plan["workloads"]

    def phase(command, directory, source, allocation, watchdog, extras=()):
        phase_start = time.monotonic()
        finalizing = directory.name == "finalization"
        remaining = outer_deadline - phase_start
        # Do not launch a phase that cannot retain its entire frozen ceiling.
        # Administrative preparation includes parent plan/request overhead.
        if directory.name == "setup":
            allocation = watchdog = min(300, setup_deadline - phase_start)
        needed = watchdog + (30 if command != "admin" else 0) + (0 if finalizing else 300)
        if allocation <= 0 or remaining < needed:
            return {"status": "outer_budget_exhausted", "cleanup": {"complete": True}, "resources": None}
        result = supervise(
            base
            + [command]
            + shared
            + ["--output", str(directory), "--deadline", "{deadline}"]
            + list(extras),
            cwd=source["path"],
            env=environment(source["path"]),
            output=directory,
            allocation=allocation,
            watchdog=watchdog,
            blocked=blocked,
            grace=30 if command != "admin" else 0,
        )
        # Include the supervisor's durable receipt write in the owned allocation.
        result["parent_observed_phase_seconds"] = time.monotonic() - phase_start
        if result["status"] == "completed" and result["parent_observed_phase_seconds"] > allocation:
            result["status"] = "allocation_exceeded"
        write(directory / "receipt.json", result)
        return result

    def checkpoint():
        write(output / "ledger.json", ledger)

    try:
        ledger["administration"]["before"] = phase("admin", output / "setup", plan["checker"], 300, 300)
        checkpoint()
        if ledger["administration"]["before"]["status"] != "completed":
            raise RuntimeError("Administrative setup failed; no worker launched")
        for index, row in enumerate(rows):
            slot = output / f"slot-{index:04d}"
            slot.mkdir()
            write(slot / "request.json", row["request"])
            row["status"] = "worker_running"
            checkpoint()
            source = plan["sources"][row["request"]["arm"]]
            row["worker"] = phase(
                "worker", slot / "worker", source, 600, 630, ["--request", str(slot / "request.json")]
            )
            row["status"] = row["worker"]["status"]
            checkpoint()
            child_terminal = (
                read(slot / "worker/terminal.json") if (slot / "worker/terminal.json").is_file() else {}
            )
            if child_terminal.get("failure_category") == "archive_coverage_failure":
                row["status"] = "archive_coverage_failure"
            if child_terminal.get("failure_category") == "integrity":
                raise IntegrityError("Worker source/input drift")
            if row["worker"]["status"] in {
                "unsafe_cleanup",
                "workload_overlap",
                "resource_verification_failed",
                "supervisor_error",
                "interrupted",
                "outer_budget_exhausted",
            }:
                raise RuntimeError("Unsafe/systemic worker condition")
            if row["worker"]["status"] == "completed":
                row["status"] = "assessor_running"
                checkpoint()
                row["assessor"] = phase(
                    "assess",
                    slot / "assessor",
                    plan["checker"],
                    300,
                    330,
                    ["--request", str(slot / "request.json"), "--worker-dir", str(slot / "worker")],
                )
                row["status"] = row["assessor"]["status"]
                if row["assessor"]["status"] == "completed":
                    row["attestation"] = read(slot / "assessor/terminal.json")
                    if row["attestation"].get("status") != "validated":
                        row["status"] = "unvalidated"
                    else:
                        row["status"] = (
                            "qualified" if row["attestation"]["qualified"] else "checked_no_candidates"
                        )
                check_terminal = (
                    read(slot / "assessor/terminal.json")
                    if (slot / "assessor/terminal.json").is_file()
                    else {}
                )
                if check_terminal.get("failure_category") == "integrity":
                    raise IntegrityError("Assessor source/input drift")
                if row["assessor"]["status"] in {
                    "unsafe_cleanup",
                    "workload_overlap",
                    "resource_verification_failed",
                    "interrupted",
                    "supervisor_error",
                    "outer_budget_exhausted",
                }:
                    raise RuntimeError("Unsafe/systemic assessor condition")
            checkpoint()
        ledger["status"] = "returned"
    except BaseException as exc:
        ledger["status"] = "failed"
        ledger["error"] = {"type": type(exc).__name__, "message": str(exc)}
    finally:
        # Never launch finalization while owned processes are still unsafe.
        unsafe = (
            any(
                r.get(p, {}).get("cleanup", {}).get("complete") is False
                for r in rows
                for p in ("worker", "assessor")
            )
            or ledger["administration"].get("before", {}).get("cleanup", {}).get("complete") is False
        )
        checkpoint()
        final_started = time.monotonic()
        if not unsafe:
            final_request = output / "final-request.json"
            write(
                final_request,
                {"ledger_path": str(output / "ledger.json"), "ledger_sha256": digest(output / "ledger.json")},
            )
            ledger["administration"]["after"] = phase(
                "admin", output / "finalization", plan["checker"], 300, 300, ["--request", str(final_request)]
            )
        ledger["summary"] = summarize(rows)
        ledger["parent_total_wall_seconds"] = time.monotonic() - outer_started
        if ledger["parent_total_wall_seconds"] > 4680:
            ledger["status"] = "outer_budget_exhausted"
        if ledger["administration"].get("after", {}).get("status") != "completed":
            ledger["status"] = "integrity_incomplete"
        else:
            final_terminal = read(output / "finalization/terminal.json")
            if "report" in final_terminal:
                ledger["detailed_report"] = final_terminal["report"]
        checkpoint()
        final_elapsed = time.monotonic() - final_started
        ledger["administration"]["parent_finalization_seconds"] = final_elapsed
        if final_elapsed > 300 or time.monotonic() > outer_deadline:
            ledger["status"] = "integrity_incomplete"
        checkpoint()
    return 0 if ledger["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
