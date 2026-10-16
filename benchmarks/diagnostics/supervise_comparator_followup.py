"""Run the separately reviewed nine-campaign plan once, sequentially, without retries."""

import argparse
from collections import Counter
from datetime import datetime, timezone
import fcntl
import hashlib
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import psutil

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from benchmarks.research import prepare_comparator_followup as preparation

GPU_UUID = "GPU-ea3de395-f78f-9f79-e363-cfac7415b6fd"
ACCEPTED = {"success", "no_candidates"}
TERMINAL = ACCEPTED | {
    "timeout",
    "memory_limit",
    "launch_error",
    "worker_error",
    "missing_response",
    "invalid_response",
    "budget_exceeded",
}
read_json, digest, encode = preparation.read_json, preparation.digest, preparation.encode


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_new(path, value):
    with Path(path).open("x") as stream:
        stream.write(encode(value))
        stream.flush()
        os.fsync(stream.fileno())


def atomic_progress(runtime, value):
    temporary = runtime / ".progress.json.tmp"
    with temporary.open("w") as stream:
        stream.write(encode({**value, "updated_utc": utc()}))
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(runtime / "progress.json")


def birth(pid):
    process = psutil.Process(pid)
    # /proc start ticks distinguish reuse without depending only on rounded wall time.
    fields = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
    return {
        "pid": pid,
        "create_time": process.create_time(),
        "start_ticks": int(fields[19]),
        "process_group_id": os.getpgid(pid),
    }


def signal_owned(process, identity, signum):
    if process.poll() is not None:
        return False
    if birth(process.pid) != identity or identity["process_group_id"] != process.pid:
        raise RuntimeError("Refuse to signal a child with changed birth/process-group identity")
    os.killpg(process.pid, signum)
    return True


def stop_owned(process, identity):
    actions = []
    if signal_owned(process, identity, signal.SIGTERM):
        actions.append("SIGTERM")
    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        if signal_owned(process, identity, signal.SIGKILL):
            actions.append("SIGKILL")
        process.wait(timeout=5)
    return actions


def load_plan(plan_path, manifest_sha256):
    plan_path = Path(plan_path).resolve()
    if plan_path.name != "execution-plan.json":
        raise ValueError("Use the prepared execution-plan.json")
    prepared = plan_path.parent
    preparation.verify_prepared(prepared, manifest_sha256)
    plan, grid = read_json(plan_path), read_json(prepared / "grid.json")
    pairs = [(s, m) for s in preparation.SPLITS for m in preparation.METHODS]
    commands = plan["commands"]
    if plan["expected_attempts"] != 270 or len(commands) != 9:
        raise ValueError("Exactly nine campaigns and 270 attempts are required")
    run_root = prepared.with_name(prepared.name + "-runs")
    expected_grid = []
    for index, ((split, method), command) in enumerate(zip(pairs, commands)):
        queries = read_json(prepared / f"public-queries/{split}.json")
        if len(queries) != preparation.COUNTS[split] or len({q["id"] for q in queries}) != len(queries):
            raise ValueError("Public query count/identity differs")
        argv = command["argv"]
        output = run_root / split / method
        protocol = str(prepared / plan["protocols"][f"{split}/{method}"])
        expected_env = {
            **preparation.THREAD_ENV,
            "PYTHONPATH": plan["source"]["root"],
            "CUDA_VISIBLE_DEVICES": "" if method == preparation.METHODS[1] else GPU_UUID,
        }
        if (
            (command["split"], command["method_id"]) != (split, method)
            or command["execution_order"] != index
            or command["no_retry_or_resume"] is not True
            or command["expected_attempts"] != len(queries) * 3
            or command["cwd"] != plan["source"]["root"]
            or command["environment_overrides"] != expected_env
            or len(argv) != 9
            or argv[0] not in plan["runtimes"]
            or argv[1:4] != ["-m", "benchmarks.research.runner", "--dataset"]
            or argv[4] not in plan["input_files_sha256"]
            or argv[5:] != ["--protocol", protocol, "--output", str(output)]
        ):
            raise ValueError("Command, ordering, GPU visibility or fixed settings differ")
        expected_grid.extend(
            {
                "split": split,
                "method_id": method,
                "query_id": query["id"],
                "target_id": query["target_id"],
                "query_index": qi,
                "seed": seed,
                "output_dir": str(output / method / f"query-{qi:04d}-seed-{seed}"),
            }
            for qi, query in enumerate(queries)
            for seed in preparation.SEEDS
        )
    if grid != expected_grid or len(grid) != 270:
        raise ValueError("Require the exact ordered 270-slot grid; no duplicates or backfill")
    return plan, grid, run_root


def validate_completion(command, slots, exit_code, source):
    """Check accounting and immutable joins, without importing chemistry or scoring."""
    if exit_code not in (0, 1):
        raise ValueError(f"Unexpected campaign exit {exit_code}")
    output = Path(command["argv"][-1])
    campaign, rows, manifest = (
        read_json(output / name) for name in ("campaign.json", "runs.json", "manifest.json")
    )
    if (
        campaign.get("status") != "completed"
        or campaign.get("termination") is not None
        or campaign.get("attempted_runs") != len(slots)
        or len(rows) != len(slots)
    ):
        raise ValueError("Incomplete or interrupted campaign; do not relaunch")
    if (
        manifest["protocol_sha256"] != digest(command["argv"][6])
        or manifest["dataset"]["dataset_sha256"] != digest(command["argv"][4])
        or manifest["code_sha256"] != source["files_sha256"]
    ):
        raise ValueError("Campaign source/protocol/dataset provenance differs")
    hashes = {name: digest(output / name) for name in ("campaign.json", "runs.json", "manifest.json")}
    for row, slot in zip(rows, slots):
        if any(row.get(k) != slot[k] for k in ("query_id", "target_id", "method_id", "seed")):
            raise ValueError("Attempt identity/order differs from exact planned slots")
        directory = Path(slot["output_dir"])
        if row.get("status") not in TERMINAL or read_json(directory / "evaluation.json") != row:
            raise ValueError("Missing terminal outcome or differing per-attempt record")
        response = directory / "response.json"
        response_hash = digest(response) if response.is_file() else None
        if row.get("response_sha256") != response_hash:
            raise ValueError("Response hash differs from retained attempt evidence")
        if row["status"] in ACCEPTED:
            if response_hash is None or read_json(response) != row["response"]:
                raise ValueError("Accepted response is not preserved")
        elif row.get("evaluation", {}).get("candidates"):
            raise ValueError("Failed attempt contains accepted candidates")
        relative = str(directory.relative_to(output))
        hashes[relative + "/evaluation.json"] = digest(directory / "evaluation.json")
        if response_hash:
            hashes[relative + "/response.json"] = response_hash
    expected_exit = 0 if all(r["status"] in ACCEPTED for r in rows) else 1
    if exit_code != expected_exit:
        raise ValueError("Campaign exit disagrees with retained outcomes")
    return {
        "validated_attempts": len(rows),
        "status_counts": dict(Counter(r["status"] for r in rows)),
        "completion_files_sha256": hashes,
    }


def fresh_location(runtime, prepared, run_root):
    runtime = Path(runtime).absolute()
    if runtime.exists() or runtime.is_symlink():
        raise FileExistsError("Supervisor runtime must be new; no retry/resume")
    resolved = runtime.resolve()
    for protected in (prepared.resolve(), run_root.resolve()):
        if resolved.is_relative_to(protected) or protected.is_relative_to(resolved):
            raise ValueError("Supervisor runtime must be separate from prepared files and runs")
    if run_root.exists() or run_root.is_symlink():
        raise FileExistsError("Planned runs already exist; no resume/backfill")
    return resolved


def supervise(runtime, lock_fd):
    runtime = Path(runtime).resolve()
    launch = read_json(runtime / "launch.json")
    manifest_path = Path(launch["plan"]).parent / "manifest.json"
    inherited, expected = os.fstat(lock_fd), manifest_path.stat()
    if (inherited.st_dev, inherited.st_ino) != (expected.st_dev, expected.st_ino):
        raise ValueError("Inherited lock does not identify the pinned prepared manifest")
    fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    write_new(runtime / "started.json", {"utc": utc(), "supervisor": birth(os.getpid())})
    state = {
        "status": "starting",
        "supervisor": birth(os.getpid()),
        "launch_sha256": digest(runtime / "launch.json"),
        "completed_campaigns": 0,
        "completed_attempts": 0,
        "campaigns": [],
        "signal": None,
    }
    previous = {}
    process, identity = None, None

    def receive(signum, frame):
        if state["signal"] is None:
            state["signal"] = signum

    try:
        for sig in (signal.SIGTERM, signal.SIGINT):
            previous[sig] = signal.signal(sig, receive)
        atomic_progress(runtime, state)
        if digest(__file__) != launch["supervisor_sha256"] or digest(launch["plan"]) != launch["plan_sha256"]:
            raise ValueError("Supervisor or execution plan changed after launch")
        plan, grid, _ = load_plan(launch["plan"], launch["manifest_sha256"])
        for command in plan["commands"]:
            if state["signal"] is not None:
                break
            # Recheck original pins before every new campaign; never repair a mismatch.
            preparation.verify_prepared(Path(launch["plan"]).parent, launch["manifest_sha256"])
            if state["signal"] is not None:
                break
            output = Path(command["argv"][-1])
            if output.exists() or output.is_symlink():
                raise FileExistsError("Campaign output already exists; do not retry")
            index = command["execution_order"]
            record = {
                "index": index,
                "command": command,
                "command_sha256": hashlib.sha256(encode(command).encode()).hexdigest(),
                "started_utc": utc(),
                "status": "launching",
            }
            state["campaigns"].append(record)
            state["status"] = "running"
            atomic_progress(runtime, state)
            log_names = (f"campaign-{index:02d}.stdout.log", f"campaign-{index:02d}.stderr.log")
            try:
                with (
                    (runtime / log_names[0]).open("xb") as stdout,
                    (runtime / log_names[1]).open("xb") as stderr,
                ):
                    process = subprocess.Popen(
                        command["argv"],
                        cwd=command["cwd"],
                        env={**os.environ, **command["environment_overrides"]},
                        stdout=stdout,
                        stderr=stderr,
                        start_new_session=True,
                    )
                    identity = birth(process.pid)
                    record.update(child=identity, status="running")
                    while process.poll() is None:
                        if state["signal"] is not None:
                            record["stop_signals"] = stop_owned(process, identity)
                            break
                        try:
                            record["observed_attempts_not_yet_validated"] = len(
                                read_json(output / "runs.json")
                            )
                        except (OSError, ValueError, TypeError):
                            pass
                        atomic_progress(runtime, state)
                        time.sleep(0.25)
                    record["exit_code"] = process.wait()
            finally:
                if process is not None and process.poll() is None:
                    record["stop_signals"] = stop_owned(process, identity)
                record["finished_utc"] = utc()
                record["logs_sha256"] = {
                    name: digest(runtime / name) for name in log_names if (runtime / name).is_file()
                }
            if state["signal"] is not None:
                record["status"] = "interrupted"
                break
            slots = [
                g for g in grid if (g["split"], g["method_id"]) == (command["split"], command["method_id"])
            ]
            record.update(validate_completion(command, slots, record["exit_code"], plan["source"]))
            record["status"] = "completed"
            state["completed_campaigns"] += 1
            state["completed_attempts"] += record["validated_attempts"]
            atomic_progress(runtime, state)
        state["status"] = "interrupted" if state["signal"] is not None else "completed"
        if state["status"] == "completed" and (state["completed_campaigns"], state["completed_attempts"]) != (
            9,
            270,
        ):
            raise ValueError("Refuse false full completion")
    except Exception as exc:
        state["status"] = "interrupted" if state["signal"] is not None else "failed"
        state["error"] = f"{type(exc).__name__}: {exc}"
        if state["campaigns"] and state["campaigns"][-1]["status"] != "completed":
            state["campaigns"][-1]["status"] = state["status"]
    finally:
        if state["signal"] is not None:
            state["status"] = "interrupted"
        state["finished_utc"] = utc()
        atomic_progress(runtime, state)
        for sig, handler in previous.items():
            signal.signal(sig, handler)
    return 0 if state["status"] == "completed" else 1


def launch(plan_path, manifest_sha256, runtime, detached=False):
    plan_path = Path(plan_path).resolve()
    # A read-only manifest descriptor provides the same advisory lock for every
    # runtime directory. The detached child inherits it; no immutable file is written.
    with (plan_path.parent / "manifest.json").open("rb") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        _, _, run_root = load_plan(plan_path, manifest_sha256)
        runtime = fresh_location(runtime, plan_path.parent, run_root)
        runtime.mkdir(parents=True, exist_ok=False)
        record = {
            "created_utc": utc(),
            "launcher": birth(os.getpid()),
            "plan": str(plan_path),
            "plan_sha256": digest(plan_path),
            "manifest_sha256": manifest_sha256,
            "supervisor_sha256": digest(__file__),
            "detached": detached,
            "retry_resume_backfill": False,
            "inherited_environment_with_exact_plan_overrides": True,
        }
        write_new(runtime / "launch.json", record)
        if not detached:
            return supervise(runtime, lock.fileno())
        argv = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            str(runtime),
            "--lock-fd",
            str(lock.fileno()),
        ]
        with (
            (runtime / "supervisor.stdout.log").open("xb") as stdout,
            (runtime / "supervisor.stderr.log").open("xb") as stderr,
        ):
            child = subprocess.Popen(
                argv,
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
                pass_fds=(lock.fileno(),),
            )
            write_new(
                runtime / "detached.json", {"argv": argv, "child": birth(child.pid), "created_utc": utc()}
            )
        return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--manifest-sha256")
    parser.add_argument("--runtime", type=Path)
    parser.add_argument("--detach", action="store_true")
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--lock-fd", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker is not None:
        if args.lock_fd is None:
            parser.error("Internal worker requires inherited manifest lock")
        return supervise(args.worker, args.lock_fd)
    if args.plan is None or args.manifest_sha256 is None or args.runtime is None:
        parser.error("--plan, --manifest-sha256 and --runtime are required")
    return launch(args.plan, args.manifest_sha256, args.runtime, args.detach)


if __name__ == "__main__":
    raise SystemExit(main())
