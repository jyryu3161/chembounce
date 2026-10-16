"""Prepare or execute the separately reviewed fixed-example Mol2Mol diagnostic.

Preparation is model-free. Execution requires a content-pinned plan and is NOT
benchmark sampling. This helper is separate from the reviewed baseline adapter.
"""

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import inspect
import os
from pathlib import Path
import random
import signal
import subprocess
import sys
import threading
import time

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from benchmarks.external.mol2mol.observation import bind_outputs, observe_decode, validate_tokens
from benchmarks.external.mol2mol.support import (
    COMMON_RDKIT,
    SETTINGS,
    digest,
    read_json,
    verify_bundle,
    upstream_origins,
    write_json,
)

FIXED_INPUT = "CCOc1ccccc1"
GPU_UUID = "GPU-ea3de395-f78f-9f79-e363-cfac7415b6fd"
PANEL_SHA256 = {
    "e87a3b539b8c41dde4887a7d52ec67b29a7a63f688584b579272398ef6522204",
    "5402b98837c4abba4ab85606f8354736ec2ef20e23e1e0e8c6dadfbd0c8f7269",
    "8d87129ebc096d122fb581eb1721b9893243cec958434b5720a94df83c2d074e",
}


def utc():
    return datetime.now(timezone.utc).isoformat()


@contextmanager
def capture_return(function):
    """Independent token observation of the exact unwrapped official return frame.

    CPU copies may synchronize CUDA. Return values and RNG are untouched. Timing
    under this profiler is diagnostic only. Never interfere with an existing profiler.
    """
    if sys.getprofile() is not None:
        raise RuntimeError("Preflight requires no existing Python profiler")
    code = inspect.unwrap(function).__code__
    captured = []

    def profile(frame, event, arg):
        if event == "return" and frame.f_code is code and arg is not None:
            # Official v4.4.22 multinomial frame: ys is the full generated tensor.
            # Observe without going through vocabulary.decode or changing tensors.
            captured.append(
                {
                    "output_tokens": frame.f_locals["ys"].detach().cpu().tolist(),
                    "input_tokens": frame.f_locals["src"].detach().cpu().tolist(),
                }
            )

    sys.setprofile(profile)
    try:
        yield captured
    finally:
        sys.setprofile(None)


def panel_check(paths):
    from rdkit import rdBase
    from chembounce_synthesis.common import canonical

    if rdBase.rdkitVersion != COMMON_RDKIT:
        raise ValueError("Wrong common RDKit runtime for exact query exclusion")
    if len(paths) != 3 or {digest(p) for p in paths} != PANEL_SHA256:
        raise ValueError("Exact three original public panels required")
    count = 0
    for path in paths:
        data = read_json(path)
        records = {record["id"]: record for record in data["compounds"]}
        for query in data["queries"]:
            count += 1
            if canonical(records[query["compound_id"]]["smiles"]) == canonical(FIXED_INPUT):
                raise ValueError(
                    "Predeclared fixed example occurs in benchmark public queries; do not replace it"
                )
    if count != 30:
        raise ValueError("Expected thirty original public queries")
    return {
        "public_queries_checked": count,
        "exact_input_present": False,
        "scope": "Exact canonical public query exclusion only; not training/heldout membership proof",
    }


def prepare_plan(bundle_path, fragment_path, panels, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError("Use a new diagnostic plan path")
    bundle = read_json(bundle_path)
    verify_bundle(bundle)
    check = panel_check(panels)
    fragment = read_json(fragment_path)
    if fragment["settings"] != SETTINGS:
        raise ValueError("Method fragment differs from fixed baseline")
    plan = {
        "schema": "mol2mol-fixed-preflight-plan-v1",
        "status": "awaiting_parent_execution_review",
        "created_utc": utc(),
        "fixed_input": FIXED_INPUT,
        "seed": 42,
        "batch_size": 32,
        "physical_gpu_index": 0,
        "gpu_uuid": GPU_UUID,
        "logical_device": "cuda:0",
        "bundle": str(bundle_path.resolve()),
        "bundle_sha256": digest(bundle_path),
        "method_fragment": str(fragment_path.resolve()),
        "method_fragment_sha256": digest(fragment_path),
        "helper": str(Path(__file__).resolve()),
        "helper_sha256": digest(__file__),
        "common_python": sys.executable,
        "common_python_sha256": digest(sys.executable),
        "panels": [{"path": str(p.resolve()), "sha256": digest(p)} for p in panels],
        "public_query_exclusion": check,
        "comparison_wall_seconds": 360,
        "comparison_max_memory_mb": 24000,
        "adapter_budget": {
            "wall_seconds": 360,
            "max_memory_mb": 24000,
            "max_candidates": 3,
            "max_evaluations": 256,
        },
        "calls": [
            "official batch32 observerOFF",
            "official batch32 observerON restored RNG",
            "one real bounded adapter attempt up to256 raw sequences",
        ],
        "interpretation": "Fixed nonbenchmark engineering diagnostic; no performance/accuracy estimates",
        "profiler_effect": "Independent return-frame tensor copies can synchronize CUDA; timing not comparable",
        "execution_authorized": False,
    }
    write_json(output, plan)


def checked_plan(path, expected):
    if digest(path) != expected:
        raise ValueError("Preflight plan SHA256 mismatch")
    plan = read_json(path)
    if plan["fixed_input"] != FIXED_INPUT or plan["seed"] != 42 or plan["gpu_uuid"] != GPU_UUID:
        raise ValueError("Fixed preflight identity/settings changed")
    for filename, pin in (
        (plan["helper"], plan["helper_sha256"]),
        (plan["bundle"], plan["bundle_sha256"]),
        (plan["method_fragment"], plan["method_fragment_sha256"]),
        (plan["common_python"], plan["common_python_sha256"]),
    ):
        if digest(filename) != pin:
            raise ValueError(f"Pinned preflight input changed: {filename}")
    if Path(plan["helper"]).resolve() != Path(__file__).resolve():
        raise ValueError("Different diagnostic helper path")
    return plan


def process_identity(pid=None):
    """Linux PID start ticks plus wall create time, independent of model packages."""
    pid = os.getpid() if pid is None else pid
    stat = Path(f"/proc/{pid}/stat").read_text()
    fields = stat[stat.rfind(")") + 2 :].split()
    start_ticks = int(fields[19])
    boot = next(
        int(line.split()[1])
        for line in Path("/proc/stat").read_text().splitlines()
        if line.startswith("btime ")
    )
    return {
        "pid": pid,
        "start_ticks": start_ticks,
        "create_time": boot + start_ticks / os.sysconf("SC_CLK_TCK"),
    }


def require_idle_gpu(snapshot):
    if snapshot["processes"]:
        raise RuntimeError(
            "Approved GPU0 already has a compute process; no inference started and none killed"
        )


def smi_snapshot():
    gpu = subprocess.check_output(
        [
            "nvidia-smi",
            "-i",
            "0",
            "--query-gpu=index,uuid,name,driver_version,pci.bus_id",
            "--format=csv,noheader",
        ],
        text=True,
        timeout=5,
    ).strip()
    fields = [x.strip() for x in gpu.split(",")]
    if fields[0] != "0" or fields[1] != GPU_UUID:
        raise RuntimeError("Physical GPU0 UUID changed")
    apps = subprocess.check_output(
        [
            "nvidia-smi",
            "-i",
            GPU_UUID,
            "--query-compute-apps=pid,gpu_uuid,used_memory",
            "--format=csv,noheader,nounits",
        ],
        text=True,
        timeout=5,
    )
    processes = []
    for line in apps.splitlines():
        values = [v.strip() for v in line.split(",")]
        if len(values) == 3 and values[0].isdigit():
            pid = int(values[0])
            try:
                identity = process_identity(pid)
            except (OSError, StopIteration):
                identity = None
            processes.append(
                {"pid": pid, "gpu_uuid": values[1], "used_memory_mb": values[2], "process_identity": identity}
            )
    return {
        "utc": utc(),
        "gpu": {
            "index": 0,
            "uuid": fields[1],
            "name": fields[2],
            "driver": fields[3],
            "pci_bus_id": fields[4],
        },
        "processes": processes,
    }


def rng_state(torch, np):
    return {
        "python": random.getstate(),
        "numpy": np.random.get_state(),
        "torch_cpu": torch.get_rng_state().clone(),
        "torch_cuda": [state.clone() for state in torch.cuda.get_rng_state_all()],
    }


def restore_rng(state, torch, np):
    random.setstate(state["python"])
    np.random.set_state(state["numpy"])
    torch.set_rng_state(state["torch_cpu"])
    torch.cuda.set_rng_state_all(state["torch_cuda"])


def rng_equal(left, right, torch, np):
    return {
        "python": left["python"] == right["python"],
        "numpy": left["numpy"][0] == right["numpy"][0]
        and np.array_equal(left["numpy"][1], right["numpy"][1])
        and left["numpy"][2:] == right["numpy"][2:],
        "torch_cpu": torch.equal(left["torch_cpu"], right["torch_cpu"]),
        "torch_cuda_all": len(left["torch_cuda"]) == len(right["torch_cuda"])
        and all(torch.equal(a, b) for a, b in zip(left["torch_cuda"], right["torch_cuda"])),
    }


def encoded_rng(state):
    # Plain diagnostic serialization; no pickle reads and no RNG-consuming operation.
    return {
        "python": state["python"],
        "numpy": [state["numpy"][0], state["numpy"][1].tolist(), *state["numpy"][2:]],
        "torch_cpu_hex": bytes(state["torch_cpu"].tolist()).hex(),
        "torch_cuda_hex": [bytes(value.tolist()).hex() for value in state["torch_cuda"]],
    }


def nll_record(values):
    return {"dtype": str(values.dtype), "shape": list(values.shape), "bytes_hex": values.tobytes().hex()}


def comparison_worker(plan, output):
    for name in ("worker-ready.json", "comparison-detail.json", "parent-live-ack.json"):
        path = output / name
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Refusing to replace comparison evidence: {name}")
    from benchmarks.external.mol2mol.worker import official_api

    if os.environ.get("CUDA_VISIBLE_DEVICES") != GPU_UUID:
        raise RuntimeError("Diagnostic must expose only the approved physical GPU0 UUID")
    bundle = read_json(plan["bundle"])
    torch, np, model_type, dataset_type, origins = official_api(bundle)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    random.seed(42)
    np.random.seed(42)
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise RuntimeError("Exactly one approved CUDA device is required")
    device = torch.device("cuda:0")
    torch.cuda.set_device(device)
    torch.manual_seed(42)
    torch.cuda.manual_seed_all(42)
    torch.cuda.reset_peak_memory_stats(device)
    model = model_type.load_from_file(bundle["model"], mode="inference", device=device)
    model.set_temperature(1.0)
    capacity = int(model.network.src_embed[1].pe.shape[1])
    validate_tokens(FIXED_INPUT, model.tokenizer, model.vocabulary, capacity)
    dataset = dataset_type([FIXED_INPUT] * 32, model.vocabulary, model.tokenizer)
    src, mask = dataset_type.collate_fn([dataset[i] for i in range(32)])
    source_tokens = src.tolist()
    src, mask = src.to(device), mask.to(device)
    devices = {
        "parameters": {name: str(value.device) for name, value in model.network.named_parameters()},
        "buffers": {name: str(value.device) for name, value in model.network.named_buffers()},
        "src": str(src.device),
        "src_mask": str(mask.device),
    }
    if not devices["parameters"]:
        raise RuntimeError("Official model has no recorded parameters")
    if (
        not all(
            value == "cuda:0"
            for group in (devices["parameters"], devices["buffers"])
            for value in group.values()
        )
        or devices["src"] != "cuda:0"
        or devices["src_mask"] != "cuda:0"
    ):
        raise RuntimeError("Model or source tensor is not on cuda:0")
    own_smi = smi_snapshot()
    if not any(p["pid"] == os.getpid() and p["gpu_uuid"] == GPU_UUID for p in own_smi["processes"]):
        raise RuntimeError("Worker PID absent from physical GPU0 compute-apps")
    write_json(
        output / "worker-ready.json",
        {
            "pid": os.getpid(),
            "process_group": os.getpgrp(),
            "process_identity": process_identity(),
            "own_smi": own_smi,
            "devices": devices,
        },
    )
    ack_deadline = time.monotonic() + 20
    while not (output / "parent-live-ack.json").exists():
        if time.monotonic() >= ack_deadline:
            raise TimeoutError("Independent parent live GPU observation missing")
        time.sleep(0.05)
    before = rng_state(torch, np)
    decode_original = model.vocabulary.decode
    decode_had_instance = "decode" in vars(model.vocabulary)
    with capture_return(model.sample) as independent:
        off = model.sample(src, mask, "multinomial")
    after_off = rng_state(torch, np)
    if len(independent) != 1:
        raise RuntimeError("Expected exactly one official sample return frame")
    restore_rng(before, torch, np)
    restored = rng_equal(before, rng_state(torch, np), torch, np)
    with observe_decode(model.vocabulary) as observed:
        on = model.sample(src, mask, "multinomial")
    after_on = rng_state(torch, np)
    rows = bind_outputs(
        observed,
        source_tokens,
        *on,
        eos=model.vocabulary["$"],
        bos=model.vocabulary["^"],
        max_length=model.max_sequence_length,
        tokenizer=model.tokenizer,
        vocabulary=model.vocabulary,
    )
    checks = {
        "input_strings": off[0] == on[0] == [FIXED_INPUT] * 32,
        "output_strings": off[1] == on[1],
        "native_nll_bytes_dtype_shape": nll_record(off[2]) == nll_record(on[2]),
        "input_token_arrays": independent[0]["input_tokens"] == observed[:32] == source_tokens,
        "output_token_arrays": independent[0]["output_tokens"] == observed[32:],
        "decode_restored": model.vocabulary.decode == decode_original
        and ("decode" in vars(model.vocabulary)) == decode_had_instance,
        "profiler_restored": sys.getprofile() is None,
        "pre_rng_restored": all(restored.values()),
        "post_rng_equal": all(rng_equal(after_off, after_on, torch, np).values()),
    }
    write_json(
        output / "comparison-detail.json",
        {
            "checks": checks,
            "off": {"inputs": off[0], "outputs": off[1], "nll": nll_record(off[2]), "tokens": independent[0]},
            "on": {
                "inputs": on[0],
                "outputs": on[1],
                "nll": nll_record(on[2]),
                "tokens": observed[32:],
                "bound_rows": rows,
            },
            "rng": {
                "before": encoded_rng(before),
                "after_off": encoded_rng(after_off),
                "after_on": encoded_rng(after_on),
                "post_equal": rng_equal(after_off, after_on, torch, np),
            },
            "input_capacity": capacity,
            "output_max_sequence_length": model.max_sequence_length,
            "devices": devices,
            "upstream_origins_before_load": origins,
            "upstream_origins_after_load_and_sampling": upstream_origins(bundle),
            "pid": os.getpid(),
            "cuda_peak_allocated": torch.cuda.max_memory_allocated(device),
            "cuda_peak_reserved": torch.cuda.max_memory_reserved(device),
            "timing_interpretation": "Profiler synchronization; no performance comparison",
        },
    )
    if not all(checks.values()):
        raise RuntimeError("Official observer OFF/ON preflight mismatch")


def stop_owned(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=1)
    except (ProcessLookupError, subprocess.TimeoutExpired):
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait(timeout=2)


def execute_plan(plan_path, expected, output):
    import psutil
    from benchmarks.research.runner import bounded_method

    plan = checked_plan(plan_path, expected)
    if Path(sys.executable).absolute() != Path(plan["common_python"]).absolute():
        raise ValueError("Execute controller using pinned common interpreter")
    panel_check([Path(record["path"]) for record in plan["panels"]])
    bundle = read_json(plan["bundle"])
    verify_bundle(bundle)
    require_idle_gpu(smi_snapshot())  # Never kill or share an occupied GPU for this diagnostic.
    output.mkdir(parents=True, exist_ok=False)
    comparison = output / "observer-comparison"
    comparison.mkdir()
    command = [
        bundle["sampler_python"],
        str(Path(__file__).resolve()),
        "--comparison-worker",
        str(plan_path),
        "--expected-plan-sha256",
        expected,
        "--output",
        str(comparison),
    ]
    env = dict(
        os.environ,
        CUDA_VISIBLE_DEVICES=GPU_UUID,
        PYTHONNOUSERSITE="1",
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        PYTHONHASHSEED="42",
    )
    snapshots = []
    previous = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}

    def interrupted(signum, frame):
        raise InterruptedError(f"Preflight cancelled by signal {signum}")

    for s in previous:
        signal.signal(s, interrupted)
    process = None
    try:
        with (comparison / "stdout.log").open("w") as stdout, (comparison / "stderr.log").open("w") as stderr:
            start = time.monotonic()
            peak_rss = 0
            process = subprocess.Popen(command, stdout=stdout, stderr=stderr, env=env, start_new_session=True)
            launch_identity = process_identity(process.pid)
            write_json(
                comparison / "controller-launch.json",
                {
                    "utc": utc(),
                    "command": command,
                    "controller": process_identity(),
                    "worker": launch_identity,
                },
            )
            while process.poll() is None:
                if time.monotonic() - start >= plan["comparison_wall_seconds"]:
                    raise TimeoutError("Comparison exceeded outer diagnostic deadline")
                try:
                    parent = psutil.Process(process.pid)
                    family = [parent, *parent.children(recursive=True)]
                except psutil.NoSuchProcess:
                    continue
                rss = 0
                for member in family:
                    try:
                        rss += member.memory_info().rss
                    except psutil.NoSuchProcess:
                        pass
                peak_rss = max(peak_rss, rss)
                if rss > plan["comparison_max_memory_mb"] * 1024**2:
                    raise MemoryError("Comparison exceeded diagnostic host RSS cap")
                snapshot = smi_snapshot()
                snapshots.append(snapshot)
                write_json(comparison / "parent-live-observations.json", snapshots)
                ready = comparison / "worker-ready.json"
                if ready.exists() and not (comparison / "parent-live-ack.json").exists():
                    pid = read_json(ready)["pid"]
                    if (
                        read_json(ready)["process_identity"] == launch_identity
                        and pid == process.pid
                        and any(p["pid"] == pid and p["gpu_uuid"] == GPU_UUID for p in snapshot["processes"])
                    ):
                        write_json(comparison / "parent-live-ack.json", {"pid": pid, "snapshot": snapshot})
                time.sleep(0.25)
            comparison_wall = time.monotonic() - start
            if comparison_wall >= plan["comparison_wall_seconds"]:
                raise TimeoutError("Comparison completed outside diagnostic deadline")
            if process.returncode:
                raise RuntimeError(f"Observer comparison failed: {process.returncode}")
        detail = read_json(comparison / "comparison-detail.json")
        if not all(detail["checks"].values()):
            raise RuntimeError("Comparison checks did not all pass")
    finally:
        try:
            if process is not None:
                stop_owned(process)
        finally:
            for s, handler in previous.items():
                signal.signal(s, handler)
    # Separate actual adapter lifecycle, same nonbenchmark input and unchanged declared budget.
    require_idle_gpu(smi_snapshot())
    fragment = read_json(plan["method_fragment"])
    assets = {name: spec["artifact"] for name, spec in fragment["assets"].items()}
    for spec in fragment["assets"].values():
        if digest(spec["artifact"]) != spec["sha256"]:
            raise ValueError("Pinned adapter asset changed before roundtrip")
    attempt = output / "adapter-roundtrip"
    request = {
        "schema_version": 1,
        "kind": "external",
        "method_id": fragment["id"],
        "query": {"id": "fixed-nonbenchmark-diagnostic", "smiles": FIXED_INPUT, "target_id": "none"},
        "seed": 42,
        "budget": plan["adapter_budget"],
        "settings": fragment["settings"],
        "assets": assets,
        "output_dir": str(attempt),
    }
    stop = threading.Event()
    adapter_observations, monitor_errors = [], []

    def monitor():
        while not stop.is_set():
            try:
                adapter_observations.append(smi_snapshot())
                write_json(output / "adapter-parent-live-observations.json", adapter_observations)
            except Exception as exc:
                monitor_errors.append(str(exc))
            stop.wait(0.25)

    thread = threading.Thread(target=monitor)
    old_visibility = os.environ.get("CUDA_VISIBLE_DEVICES")
    os.environ["CUDA_VISIBLE_DEVICES"] = GPU_UUID
    thread.start()
    try:
        result = bounded_method(fragment, request, attempt)
    finally:
        stop.set()
        thread.join(timeout=10)
        if old_visibility is None:
            os.environ.pop("CUDA_VISIBLE_DEVICES", None)
        else:
            os.environ["CUDA_VISIBLE_DEVICES"] = old_visibility
    write_json(output / "adapter-result.json", result)
    if result["status"] not in {"success", "no_candidates"} or monitor_errors or thread.is_alive():
        raise RuntimeError(f"Real adapter diagnostic failed: {result['status']}; monitor={monitor_errors}")
    raw = read_json(attempt / "mol2mol-raw.json")
    if not any(
        p["pid"] == raw["pid"] and p["gpu_uuid"] == GPU_UUID
        for snapshot in adapter_observations
        for p in snapshot["processes"]
    ):
        raise RuntimeError("Parent did not observe actual adapter sampler PID on approved GPU")
    provenance = read_json(attempt / "mol2mol-response-provenance.json")
    if provenance["response_sha256"] != digest(attempt / "response.json"):
        raise RuntimeError("Actual adapter response hash mismatch")
    checked_plan(plan_path, expected)
    verify_bundle(bundle)
    write_json(
        output / "validation.json",
        {
            "status": "fixed_example_preflight_passed",
            "completed_utc": utc(),
            "plan_sha256": expected,
            "comparison_detail_sha256": digest(comparison / "comparison-detail.json"),
            "adapter_result_sha256": digest(output / "adapter-result.json"),
            "adapter_status": result["status"],
            "comparison_wall_seconds": comparison_wall,
            "comparison_sampled_peak_rss_bytes": peak_rss,
            "raw_sample_count": len(raw["rows"]),
            "candidate_count": len(result["response"]["candidates"]),
            "benchmark_attempts": 0,
            "scientific_accuracy_or_performance_claim": False,
            "scope": "One fixed nonbenchmark query, seed42, model observer equality and actual adapter lifecycle",
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare-plan", action="store_true")
    mode.add_argument("--execute-plan", type=Path)
    mode.add_argument("--comparison-worker", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--method-fragment", type=Path)
    parser.add_argument("--panels", type=Path, nargs=3)
    parser.add_argument("--expected-plan-sha256")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.prepare_plan:
        if not args.bundle or not args.method_fragment or not args.panels:
            parser.error("Preparation requires bundle, method fragment and three pinned panels")
        prepare_plan(args.bundle, args.method_fragment, args.panels, args.output)
    elif args.comparison_worker:
        comparison_worker(checked_plan(args.comparison_worker, args.expected_plan_sha256), args.output)
    else:
        existed_before = args.output.exists()
        try:
            execute_plan(args.execute_plan, args.expected_plan_sha256, args.output)
        except BaseException as exc:
            if not existed_before and args.output.is_dir() and not (args.output / "failure.json").exists():
                write_json(
                    args.output / "failure.json",
                    {
                        "status": "failed",
                        "utc": utc(),
                        "type": type(exc).__name__,
                        "reason": str(exc),
                        "no_benchmark_inference": True,
                    },
                )
            raise


if __name__ == "__main__":
    main()
