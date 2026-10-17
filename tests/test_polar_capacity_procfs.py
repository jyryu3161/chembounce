"""Deterministic disappearing-procfs fixtures; no scientific workloads."""

import errno
import json
import os
from pathlib import Path
import sys

import pytest

from benchmarks.research import polar_capacity_supervisor as supervisor


@pytest.mark.parametrize("number", [errno.ENOENT, errno.ESRCH])
def test_identity_disappearance_is_absent(monkeypatch, number):
    error = OSError(number, "controlled disappearing process")

    def missing(path, *args, **kwargs):
        raise error

    monkeypatch.setattr(Path, "read_text", missing)
    assert supervisor.identity(1073741824) is None


@pytest.mark.parametrize("number", [errno.EACCES, errno.EPERM, errno.EIO, errno.EMFILE, errno.EINVAL])
def test_identity_other_io_errors_remain_visible(monkeypatch, number):
    error = OSError(number, "controlled non-disappearance failure")

    def failed(path, *args, **kwargs):
        raise error

    monkeypatch.setattr(Path, "read_text", failed)
    with pytest.raises(OSError) as caught:
        supervisor.identity(1073741824)
    assert caught.value is error


@pytest.mark.parametrize("number", [errno.ENOENT, errno.ESRCH, errno.EACCES, errno.EIO])
def test_unrelated_procfs_race_during_owned_child(number, tmp_path, monkeypatch):
    missing_pid = 1073741824  # Beyond Linux's maximum PID; never an actual process.
    target = Path(f"/proc/{missing_pid}/stat")
    read_text, listdir = Path.read_text, os.listdir
    snapshot, signal_owned = supervisor.owned_snapshot, supervisor.signal_owned
    before = supervisor.subreaper()
    injected, scans, signals = [], [], []
    release = tmp_path / "release"

    def listing(path):
        values = listdir(path)
        return values + [str(missing_pid)] if str(path) == "/proc" else values

    def disappearing(path, *args, **kwargs):
        if path == target and not injected:
            injected.append(number)
            raise OSError(number, "controlled procfs read race", str(path))
        return read_text(path, *args, **kwargs)

    def scanned(leader, known):
        scans.append(True)
        result = snapshot(leader, known)
        if len(scans) == 1:
            # The child completes only after the first scan returns safely. An
            # unexpected scan error leaves it waiting for owned cleanup.
            release.touch()
        return result

    def signaled(entry, signum):
        signals.append((entry["pid"], signum))
        return signal_owned(entry, signum)

    monkeypatch.setattr(os, "listdir", listing)
    monkeypatch.setattr(Path, "read_text", disappearing)
    monkeypatch.setattr(supervisor, "owned_snapshot", scanned)
    monkeypatch.setattr(supervisor, "signal_owned", signaled)
    code = (
        "import json,time\nfrom pathlib import Path\n"
        "while not Path('release').exists(): time.sleep(0.01)\n"
        "Path('child/terminal.json').write_text(json.dumps({'status':'healthy_fixture_completed'}))\n"
    )
    result = supervisor.supervise(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env={"PATH": "/usr/bin:/bin", "PYTHONDONTWRITEBYTECODE": "1"},
        output=tmp_path / "child",
        allocation=3,
        watchdog=3,
        grace=1,
        interval=0.01,
    )
    assert injected == [number] and result["cleanup"]["complete"]
    assert supervisor.subreaper() == before
    assert not any(pid == missing_pid for pid, _ in signals)
    if number in (errno.ENOENT, errno.ESRCH):
        assert result["status"] == "completed"
        assert result["returncode"] == 0 and signals == []
        assert (
            json.loads((tmp_path / "child/terminal.json").read_text())["status"]
            == "healthy_fixture_completed"
        )
    else:
        assert result["status"] == "supervisor_error"
        assert result["error"]["type"] == type(OSError(number, "fixture")).__name__
        assert not release.exists() and not (tmp_path / "child/terminal.json").exists()
        assert signals and all(pid == result["identity"]["pid"] for pid, _ in signals)


def test_live_polls_scan_once_and_exit_refreshes(tmp_path, monkeypatch):
    """A held real child exposes healthy polls and the actual wait4 exit boundary."""
    scan, reap, sleep = supervisor.owned_snapshot, supervisor.reap_owned, supervisor.time.sleep
    events, polls = [], []
    release = tmp_path / "release"

    def scanned(leader, known):
        events.append("scan")
        return scan(leader, known)

    def reaped(process, records, resources):
        result = reap(process, records, resources)
        events.append("exit" if process.returncode is not None else "live")
        return result

    def slept(seconds):
        polls.append(list(events))
        events.clear()
        if len(polls) == 3:
            release.touch()
        sleep(seconds)

    monkeypatch.setattr(supervisor, "owned_snapshot", scanned)
    monkeypatch.setattr(supervisor, "reap_owned", reaped)
    monkeypatch.setattr(supervisor.time, "sleep", slept)
    code = (
        "import time\nfrom pathlib import Path\n"
        "while not Path('release').exists(): time.sleep(0.01)\n"
        "Path('child/terminal.json').write_text('{}')\n"
    )
    result = supervisor.supervise(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=os.environ.copy(),
        output=tmp_path / "child",
        allocation=3,
        watchdog=3,
        grace=1,
    )
    assert result["status"] == "completed" and result["cleanup"]["complete"]
    assert polls[:3] == [["scan", "live"]] * 3
    # The exit poll refreshes after wait4 and finally scans before success.
    assert events == ["scan", "exit", "scan", "scan"]
    assert result["resources"]["sample_interval_seconds"] == 0.2


@pytest.mark.parametrize("boundary", ["adopted", "interrupted", "watchdog"])
def test_boundary_after_reap_forces_fresh_scan(tmp_path, monkeypatch, boundary):
    """Inject a boundary between the first snapshot and the post-reap decision."""
    import signal
    import time

    scan, reap, sleep = supervisor.owned_snapshot, supervisor.reap_owned, supervisor.time.sleep
    events, fired = [], []
    offset = [0.0]
    release = tmp_path / "release"

    def scanned(leader, known):
        events.append("scan")
        return scan(leader, known)

    def reaped(process, records, resources):
        result = reap(process, records, resources)
        events.append("reap")
        if not fired and process.returncode is None:
            fired.append(True)
            if boundary == "adopted":
                # Actual adopted wait4 accounting is covered independently;
                # this isolates detection of a newly appended usage record.
                resources["adopted_child_usage"].append({"synthetic": True})
            elif boundary == "interrupted":
                signal.getsignal(signal.SIGINT)(signal.SIGINT, None)
            else:
                offset[0] = 3.01
        return result

    def slept(seconds):
        release.touch()
        sleep(seconds)

    monkeypatch.setattr(supervisor, "owned_snapshot", scanned)
    monkeypatch.setattr(supervisor, "reap_owned", reaped)
    monkeypatch.setattr(supervisor.time, "sleep", slept)
    code = (
        "import time\nfrom pathlib import Path\n"
        "while not Path('release').exists(): time.sleep(0.01)\n"
        "Path('child/terminal.json').write_text('{}')\n"
    )
    result = supervisor.supervise(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=os.environ.copy(),
        output=tmp_path / "child",
        allocation=3,
        watchdog=3,
        grace=1,
        clock=lambda: time.monotonic() + offset[0],
    )
    assert fired and events[:3] == ["scan", "reap", "scan"]
    assert result["cleanup"]["complete"]
    assert (
        result["status"]
        == {
            "adopted": "completed",
            "interrupted": "interrupted",
            "watchdog": "watchdog_timeout",
        }[boundary]
    )
