"""Linux owned-process supervision for the fixed development ablation.

No scientific imports. RSS is sampled, not an instantaneous memory guarantee.
The watchdog/cleanup interval never enlarges the successful work allocation.
"""

import ctypes
import errno
import hashlib
import json
import os
from pathlib import Path
import signal
import resource
import subprocess
import time

import psutil


def write(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)
    fd = os.open(path.parent, os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def identity(pid):
    try:
        with open(f"/proc/{pid}/stat") as stream:
            fields = stream.read().rsplit(")", 1)[1].split()
        return {
            "pid": pid,
            "start_ticks": int(fields[19]),
            "ppid": int(fields[1]),
            "pgrp": int(fields[2]),
            "session": int(fields[3]),
            "state": fields[0],
        }
    except OSError as exc:
        # procfs can report either error when a process exits during the read.
        # Permission and unexpected I/O errors remain supervision failures.
        if exc.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise


def same(left, right):
    return left is not None and right is not None and all(left[k] == right[k] for k in ("pid", "start_ticks"))


def reap_owned(process, records, resources):
    """One nonblocking pass; retain each owned child's actual wait4 usage.

    A snapshot is not evidence of reaping. Missing/unwaitable children remain
    visible to the caller's next ownership snapshot and cannot certify cleanup.
    """
    leader_usage = None
    for pid, entry in records.items():
        if entry["ppid"] != os.getpid() or not same(identity(pid), entry):
            continue
        if pid == process.pid and process.returncode is not None:
            continue
        try:
            reaped, status, usage = os.wait4(pid, os.WNOHANG)
        except ChildProcessError:
            continue
        if not reaped:
            continue
        code = os.waitstatus_to_exitcode(status)
        if pid == process.pid:
            process.returncode = code
            leader_usage = usage
        else:
            resources["adopted_child_usage"].append(
                {
                    "identity": entry,
                    "ru_maxrss_kib": usage.ru_maxrss,
                    "cpu_seconds": usage.ru_utime + usage.ru_stime,
                    "returncode": code,
                }
            )
    return leader_usage


def subreaper(value=None):
    libc = ctypes.CDLL(None, use_errno=True)
    if value is None:
        result = ctypes.c_int()
        status = libc.prctl(37, ctypes.byref(result), 0, 0, 0)
        answer = result.value
    else:
        status = libc.prctl(36, int(value), 0, 0, 0)
        answer = int(value)
    if status != 0:
        err = ctypes.get_errno()
        raise OSError(err, os.strerror(err))
    return answer


def pidfd_open(pid):
    # The pinned Python 3.10 build lacks os.pidfd_open. Linux x86_64 syscall
    # numbers are explicit; no racy kill(pid) fallback is permitted.
    if os.uname().machine != "x86_64":
        raise RuntimeError("Pinned supervisor requires Linux x86_64")
    libc = ctypes.CDLL(None, use_errno=True)
    result = libc.syscall(434, ctypes.c_int(pid), ctypes.c_uint(0))
    if result < 0:
        err = ctypes.get_errno()
        raise OSError(err, os.strerror(err))
    return result


def pidfd_send_signal(fd, number):
    libc = ctypes.CDLL(None, use_errno=True)
    result = libc.syscall(424, ctypes.c_int(fd), ctypes.c_int(number), ctypes.c_void_p(0), ctypes.c_uint(0))
    if result < 0:
        err = ctypes.get_errno()
        raise OSError(err, os.strerror(err))


def signal_owned(entry, number):
    """Pin the kernel process object, then recheck identity before signaling."""
    try:
        fd = pidfd_open(entry["pid"])
    except ProcessLookupError:
        return False
    try:
        if not same(identity(entry["pid"]), entry):
            return False
        pidfd_send_signal(fd, number)
        return True
    except ProcessLookupError:
        return False
    finally:
        os.close(fd)


def workload_check(blocked):
    """Only process metadata; never inspect a blocked job's scientific files."""
    if not isinstance(blocked, dict):
        return [entry for entry in blocked if same(identity(entry["pid"]), entry)]
    matches = [entry for entry in blocked["identities"] if same(identity(entry["pid"]), entry)]
    for process in psutil.process_iter(["pid", "cmdline"]):
        try:
            command = process.info["cmdline"] or []
            if any(marker in command for marker in blocked["command_markers"]):
                item = identity(process.pid)
                if item:
                    matches.append(item)
        except psutil.NoSuchProcess:
            continue
        except psutil.AccessDenied as exc:
            raise RuntimeError("Known workload scan is incomplete") from exc
    return matches


def owned_snapshot(leader, known):
    # A subreaper receives daemonized descendants; no other child is launched by
    # this supervisor while a phase is active. Session membership covers an
    # orphan before adoption; tracked identities cover setsid/reparenting.
    records = {}
    supervisor_pid = os.getpid()
    for name in os.listdir("/proc"):
        if name.isdigit():
            item = identity(int(name))
            if item is not None:
                records[item["pid"]] = item
    anchored_session = any(
        same(item, known.get(pid)) and item["session"] == leader["pid"] for pid, item in records.items()
    )
    selected = {
        pid: item
        for pid, item in records.items()
        if same(item, known.get(pid))
        or item["ppid"] == supervisor_pid
        or (anchored_session and item["session"] == leader["pid"])
    }
    changed = True
    while changed:
        changed = False
        for pid, item in records.items():
            if pid not in selected and item["ppid"] in selected:
                selected[pid] = item
                changed = True
    known.update(selected)
    return selected


def rss_sample(records):
    total, gaps = 0, []
    for pid, entry in records.items():
        if entry["state"] == "Z":
            continue
        try:
            if not same(identity(pid), entry):
                continue  # Expected exit/reuse, never read an unrelated PID.
            total += psutil.Process(pid).memory_info().rss
        except psutil.NoSuchProcess:
            continue
        except (psutil.AccessDenied, OSError) as exc:
            gaps.append({"identity": entry, "error": type(exc).__name__})
    return {"rss_bytes": None if gaps else total, "gaps": gaps}


def parent_memory():
    try:
        current = psutil.Process().memory_info().rss
    except psutil.Error:
        current = None
    return {
        "point_rss_bytes": current,
        "lifetime_ru_maxrss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "scope": "parent only; not additive to independently sampled child peaks",
    }


def terminal_status(reason, returncode, elapsed, allocation, cleaned, durable):
    if not cleaned:
        return "unsafe_cleanup"
    if reason is not None:
        return reason
    if elapsed > allocation:
        return "allocation_exceeded"
    if returncode != 0:
        return "child_error"
    if not durable:
        return "missing_terminal"
    return "completed"


def supervise(
    argv,
    *,
    cwd,
    env,
    output,
    allocation,
    watchdog,
    grace=30,
    blocked=(),
    memory_bytes=16 * 1024**3,
    interval=0.2,
    clock=time.monotonic,
    sample=rss_sample,
    check_workload=workload_check,
):
    """Own exactly one fresh child tree; persist starts and terminal receipts.

    argv may contain {deadline}, replaced with the parent absolute allocation.
    Caller must authenticate argv/cwd/environment and returned terminal payload.
    The pure classifier and bounded fake children are the test seams.
    """
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    if not (0 < allocation <= watchdog and 0 <= grace <= 30 and interval > 0):
        raise ValueError("Invalid fixed phase limits")
    previous_subreaper = subreaper()
    subreaper_changed = False
    probe = pidfd_open(os.getpid())
    os.close(probe)
    t0 = clock()
    receipt = {
        "status": "running",
        "started_monotonic": t0,
        "allocation_seconds": allocation,
        "watchdog_seconds": watchdog,
        "kill_grace_seconds": grace,
        "argv": [str(a).replace("{deadline}", str(t0 + allocation)) for a in argv],
        "cwd": str(cwd),
        "parent_before": parent_memory(),
        "supervisor_identity": identity(os.getpid()),
        "identity": None,
        "returncode": None,
        "cleanup": None,
        "resources": {
            "sampled_tree_peak_rss_bytes": None,
            "lifetime_ru_maxrss_kib": None,
            "child_cpu_seconds": None,
            "samples": 0,
            "gaps": [],
            "rss_scope": "owned simultaneous RSS sum; shared pages may be double counted",
            "instantaneous_cap": False,
            "sample_interval_seconds": interval,
            "maximum_observed_sample_gap_seconds": None,
            "last_sample_monotonic": None,
            "adopted_child_usage": [],
            "lifetime_ru_maxrss_scope": "Linux wait4 immediate child; adopted children reported separately, not summed peaks",
        },
    }
    write(output / "receipt.json", receipt)
    process, reason, usage, leader, known = None, None, None, None, {}
    killed_at, terminated_at, last_workload, last_sample = None, None, -float("inf"), -float("inf")
    old_handlers = {}
    interrupted = []
    for number in (signal.SIGTERM, signal.SIGINT):
        old_handlers[number] = signal.signal(number, lambda n, f: interrupted.append(n))

    def stop(records, number):
        for pid, entry in records.items():
            if entry["state"] != "Z" and same(identity(pid), entry):
                try:
                    signal_owned(entry, number)
                except ProcessLookupError:
                    pass

    try:
        if psutil.Process().children(recursive=True):
            reason = "preexisting_children"
        elif check_workload(blocked):
            reason = "workload_overlap"
        else:
            subreaper(1)
            subreaper_changed = previous_subreaper != 1
            with (output / "stdout.log").open("w") as stdout, (output / "stderr.log").open("w") as stderr:
                process = subprocess.Popen(
                    receipt["argv"], cwd=cwd, env=env, stdout=stdout, stderr=stderr, start_new_session=True
                )
            leader = identity(process.pid)
            if leader is None:
                # A zombie still has /proc identity until wait4; disappearance is
                # an ownership verification failure, never permission to signal.
                raise RuntimeError("Lost launch identity")
            known[process.pid] = leader
            receipt["identity"] = leader
            write(output / "receipt.json", receipt)
            while True:
                now = clock()
                records = owned_snapshot(leader, known)
                if now - last_workload >= 1:
                    last_workload = now
                    if check_workload(blocked):
                        reason = "workload_overlap"
                if interrupted:
                    reason = "interrupted"
                if now - last_sample >= interval:
                    previous_sample = receipt["resources"]["last_sample_monotonic"]
                    if previous_sample is not None:
                        receipt["resources"]["maximum_observed_sample_gap_seconds"] = max(
                            receipt["resources"]["maximum_observed_sample_gap_seconds"] or 0,
                            now - previous_sample,
                        )
                    receipt["resources"]["last_sample_monotonic"] = now
                    last_sample = now
                    observation = sample(records)
                    receipt["resources"]["samples"] += 1
                    if observation["gaps"]:
                        receipt["resources"]["gaps"].extend(observation["gaps"])
                        reason = "resource_verification_failed"
                    if observation["rss_bytes"] is not None:
                        peak = receipt["resources"]["sampled_tree_peak_rss_bytes"]
                        receipt["resources"]["sampled_tree_peak_rss_bytes"] = max(
                            peak or 0, observation["rss_bytes"]
                        )
                        if observation["rss_bytes"] > memory_bytes:
                            reason = "memory_limit"
                adopted_before = len(receipt["resources"]["adopted_child_usage"])
                current_usage = reap_owned(process, records, receipt["resources"])
                if current_usage is not None:
                    usage = current_usage
                # One full ownership scan per healthy live poll is sufficient;
                # new descendants are observed on the next poll. Reaping and
                # every exit/error/deadline boundary still get a fresh scan.
                now = clock()
                if (
                    process.returncode is not None
                    or len(receipt["resources"]["adopted_child_usage"]) != adopted_before
                    or reason is not None
                    or interrupted
                    or terminated_at is not None
                    or now - t0 >= watchdog
                ):
                    records = owned_snapshot(leader, known)
                live = {p: i for p, i in records.items() if i["state"] != "Z"}
                # /proc scans and RSS sampling can overrun a polling interval.
                # Never make escalation decisions with the pre-scan timestamp.
                now = clock()
                if interrupted and reason is None:
                    reason = "interrupted"
                if process.returncode is not None and not records:
                    break
                if now - t0 >= watchdog and reason is None:
                    reason = "watchdog_timeout"
                # A returned leader may not leave descendants behind as free work.
                if process.returncode is not None and live and reason is None:
                    reason = "owned_descendants_after_exit"
                if reason and terminated_at is None:
                    terminated_at = now
                    stop(owned_snapshot(leader, known), signal.SIGTERM)
                    now = clock()
                if terminated_at is not None and now >= min(
                    terminated_at + max(0, grace - 0.1), t0 + watchdog + max(0, grace - 0.1)
                ):
                    if killed_at is None:
                        killed_at = now
                        stop(owned_snapshot(leader, known), signal.SIGKILL)
                        now = clock()
                    # Bounded reaping only; uninterruptible work is unsafe, not success.
                    if now >= t0 + watchdog + grace or now - killed_at >= 0.2:
                        break
                time.sleep(min(interval, 0.05))
    except BaseException as exc:
        reason = "supervisor_error"
        receipt["error"] = {"type": type(exc).__name__, "message": str(exc)}
        if leader:
            stop(owned_snapshot(leader, known), signal.SIGKILL)
    finally:
        for number, handler in old_handlers.items():
            signal.signal(number, handler)
        remaining = owned_snapshot(leader, known) if leader else {}
        # A child may become a zombie between wait4(WNOHANG) and snapshot.
        # Observe, reap, then observe again, with at most two nonblocking passes
        # and no sleep/new work allowance. Live or unreapable owned processes
        # still make cleanup unsafe. Any deadline overrun stays in wall time.
        for _ in range(2):
            if not process or not remaining:
                break
            current_usage = reap_owned(process, remaining, receipt["resources"])
            if current_usage is not None:
                usage = current_usage
            remaining = owned_snapshot(leader, known)
        cleaned = not remaining
        if cleaned and subreaper_changed:
            subreaper(previous_subreaper)
        receipt["subreaper"] = {
            "before": previous_subreaper,
            "after": subreaper(),
            "restored_when_safe": cleaned,
            "unsafe_retained_ownership": not cleaned,
        }
        elapsed = clock() - t0
        durable = (output / "terminal.json").is_file()
        receipt["terminal_sha256"] = (
            hashlib.sha256((output / "terminal.json").read_bytes()).hexdigest() if durable else None
        )
        # Account for a signal observed after the loop check, including the
        # boundary where the original handlers are restored before cleanup.
        if interrupted and reason is None:
            reason = "interrupted"
        receipt.update(
            status=terminal_status(
                reason, process.returncode if process else None, elapsed, allocation, cleaned, durable
            ),
            returncode=process.returncode if process else None,
            wall_seconds=elapsed,
            cleanup={
                "complete": cleaned,
                "remaining": list(remaining.values()),
                "term_monotonic": terminated_at,
                "kill_monotonic": killed_at,
            },
        )
        if usage:
            receipt["resources"].update(
                lifetime_ru_maxrss_kib=usage.ru_maxrss, child_cpu_seconds=usage.ru_utime + usage.ru_stime
            )
        receipt["parent_after"] = parent_memory()
        write(output / "receipt.json", receipt)
        # The owning caller also measures through this durable write. If the
        # write itself crossed the success allocation, correct the receipt.
        after_durable = clock() - t0
        if receipt["status"] == "completed" and after_durable > allocation:
            receipt["status"] = "allocation_exceeded"
            receipt["wall_seconds"] = after_durable
            write(output / "receipt.json", receipt)
    return receipt
