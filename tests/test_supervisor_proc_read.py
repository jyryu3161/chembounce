"""Differential contracts for the two dd07 procfs hot-path edits; no science."""

import copy
import errno
import io
import os
from pathlib import Path

import pytest

from benchmarks.research import polar_capacity_supervisor as s


def previous_identity(pid):
    """Frozen dd07 parsing/error oracle, before direct text open."""
    try:
        fields = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
        return dict(
            pid=pid,
            start_ticks=int(fields[19]),
            ppid=int(fields[1]),
            pgrp=int(fields[2]),
            session=int(fields[3]),
            state=fields[0],
        )
    except OSError as exc:
        if exc.errno in (errno.ENOENT, errno.ESRCH):
            return None
        raise


def outcome(function):
    try:
        return "value", function()
    except Exception as exc:
        return "error", type(exc), str(exc)


@pytest.mark.parametrize("comm", ["normal", "two words", "a)b(c)d", "( ) with ) tail", "nonascii-é"])
def test_identity_text_parser_and_context_parity(monkeypatch, comm):
    fields = ["S", "40", "41", "42"] + ["0"] * 15 + ["123456"] + ["0"] * 20
    raw = f"17 ({comm}) " + " ".join(fields) + "\n"
    opened = []

    def opening(path, *args, **kwargs):
        assert path == "/proc/17/stat" and not args and not kwargs
        stream = io.StringIO(raw)
        opened.append(stream)
        return stream

    monkeypatch.setattr(s, "open", opening, raising=False)
    monkeypatch.setattr(Path, "read_text", lambda path: raw)
    assert (
        s.identity(17)
        == previous_identity(17)
        == dict(pid=17, start_ticks=123456, ppid=40, pgrp=41, session=42, state="S")
    )
    assert len(opened) == 1 and opened[0].closed


@pytest.mark.parametrize(
    "raw",
    [
        "",
        "no delimiter",
        "17 (name) S",
        "17 (name) S x 1 1 " + "0 " * 30,
        "17 (name) S 1 1 1 " + "0 " * 15 + "bad",
    ],
)
def test_parse_failures_match_previous_and_close(monkeypatch, raw):
    streams = []

    def opening(*args):
        stream = io.StringIO(raw)
        streams.append(stream)
        return stream

    monkeypatch.setattr(s, "open", opening, raising=False)
    monkeypatch.setattr(Path, "read_text", lambda path: raw)
    assert outcome(lambda: s.identity(17)) == outcome(lambda: previous_identity(17))
    assert streams[0].closed


@pytest.mark.parametrize(
    "number", [errno.ENOENT, errno.ESRCH, errno.EACCES, errno.EPERM, errno.EIO, errno.EMFILE]
)
def test_read_failure_parity_closes_stream(monkeypatch, number):
    error = OSError(number, "controlled read failure")

    class Failing(io.StringIO):
        def read(self, *args):
            raise error

    stream = Failing()

    def fail(*args):
        raise error

    monkeypatch.setattr(s, "open", lambda path: stream, raising=False)
    monkeypatch.setattr(Path, "read_text", fail)
    assert outcome(lambda: s.identity(17)) == outcome(lambda: previous_identity(17))
    assert stream.closed


def previous_snapshot(leader, known):
    """Frozen dd07 ownership selection oracle; intentionally no PID hoist."""
    records = {}
    for name in os.listdir("/proc"):
        if name.isdigit():
            item = s.identity(int(name))
            if item is not None:
                records[item["pid"]] = item
    anchored_session = any(
        s.same(item, known.get(pid)) and item["session"] == leader["pid"] for pid, item in records.items()
    )
    selected = {
        pid: item
        for pid, item in records.items()
        if s.same(item, known.get(pid))
        or item["ppid"] == os.getpid()
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


@pytest.mark.parametrize("order", ["forward", "reverse"])
@pytest.mark.parametrize("anchor", ["leader", "tracked_child", "none"])
def test_full_discovery_selection_order_pid_reuse_and_closure(monkeypatch, order, anchor):
    def entry(pid, ppid, session, start=10):
        return dict(pid=pid, ppid=ppid, session=session, pgrp=session, start_ticks=start, state="S")

    leader = entry(100, 900, 100)
    records = {
        100: leader,
        120: entry(120, 100, 100),
        130: entry(130, 120, 130),
        140: entry(140, 900, 140),
        150: entry(150, 999, 100),
        160: entry(160, 999, 777),
        170: entry(170, 999, 777, 20),
        180: entry(180, 999, 777),
    }
    known = {100: copy.deepcopy(leader), 160: copy.deepcopy(records[160]), 170: entry(170, 999, 777, 19)}
    if anchor != "leader":
        records.pop(100)
        known.pop(100)
        if anchor == "tracked_child":
            known[120] = copy.deepcopy(records[120])
    names = [str(pid) for pid in records]
    if order == "reverse":
        names.reverse()
    names += ["self", "thread-self", "999999"]
    reads = []
    pid_calls = []
    monkeypatch.setattr(os, "listdir", lambda path: names)
    monkeypatch.setattr(os, "getpid", lambda: pid_calls.append(True) or 900)
    monkeypatch.setattr(s, "identity", lambda pid: reads.append(pid) or records.get(pid))
    prior_known = copy.deepcopy(known)
    expected = previous_snapshot(leader, prior_known)
    reads.clear()
    pid_calls.clear()
    actual = s.owned_snapshot(leader, known)
    assert list(actual.items()) == list(expected.items())
    assert list(known.items()) == list(prior_known.items())
    assert reads == [int(name) for name in names if name.isdigit()]
    assert len(pid_calls) == 1
    assert 140 in actual and 160 in actual and 170 not in actual and 180 not in actual
    assert (150 in actual) == (anchor != "none")
    assert (130 in actual) == (anchor != "none")
