"""Model-free checks of the separate, not-yet-executed CUDA diagnostic helper."""

import copy
from functools import wraps
import random
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from rdkit import rdBase

from benchmarks.diagnostics import mol2mol_preflight as preflight


class PlainTensor:
    def __init__(self, value):
        self.value = value

    def detach(self):
        return self

    def cpu(self):
        return self

    def tolist(self):
        return copy.deepcopy(self.value)


def wrapped(function):
    @wraps(function)
    def call(*args, **kwargs):
        return function(*args, **kwargs)

    return call


@wrapped
def sample_fixture(src):
    ys = PlainTensor([[1, 3, 2, 0], [1, 4, 4, 4]])
    return src, ys


def test_independent_profiler_observes_only_exact_return_frame_and_preserves_rng():
    src = PlainTensor([[1, 3, 2], [1, 3, 2]])
    state = random.getstate()
    with preflight.capture_return(sample_fixture) as records:
        output = sample_fixture(src)
        assert output[0] is src
    assert records == [{"input_tokens": src.value, "output_tokens": output[1].value}]
    assert sys.getprofile() is None and random.getstate() == state
    output[1].value[0][0] = 99
    assert records[0]["output_tokens"][0][0] == 1


def test_profiler_restores_after_original_exception_without_masking_it():
    error = ValueError("original sampler failure")

    def failing(src):
        raise error

    with pytest.raises(ValueError) as raised:
        with preflight.capture_return(failing) as records:
            failing(PlainTensor([]))
    assert raised.value is error and records == [] and sys.getprofile() is None


def test_existing_profiler_is_not_overwritten():
    def existing(frame, event, arg):
        pass

    sys.setprofile(existing)
    try:
        with pytest.raises(RuntimeError, match="existing Python profiler"):
            with preflight.capture_return(sample_fixture):
                pytest.fail("body must not execute")
        assert sys.getprofile() is existing
    finally:
        sys.setprofile(None)


def test_native_nll_comparison_is_bitwise_including_dtype_shape_and_signed_zero():
    original = np.array([0, 1], dtype=np.float32)
    record = preflight.nll_record(original)
    assert record == preflight.nll_record(original.copy())
    assert record != preflight.nll_record(np.array([-0.0, 1], dtype=np.float32))
    assert record != preflight.nll_record(original.astype(np.float64))
    assert record != preflight.nll_record(original.reshape(1, 2))


def test_full_rng_equality_detects_every_visible_cuda_and_other_stream_changes():
    state = {
        "python": random.getstate(),
        "numpy": np.random.RandomState(42).get_state(),
        "torch_cpu": b"cpu",
        "torch_cuda": [b"gpu0", b"gpu1"],
    }
    torch = SimpleNamespace(equal=lambda a, b: a == b)
    assert all(preflight.rng_equal(state, copy.deepcopy(state), torch, np).values())
    changed = copy.deepcopy(state)
    changed["torch_cuda"][1] = b"different second visible device"
    assert not preflight.rng_equal(state, changed, torch, np)["torch_cuda_all"]
    changed = copy.deepcopy(state)
    changed["torch_cpu"] = b"different"
    assert not preflight.rng_equal(state, changed, torch, np)["torch_cpu"]
    changed = copy.deepcopy(state)
    changed["numpy"][1][5] += 1
    assert not preflight.rng_equal(state, changed, torch, np)["numpy"]


def test_plan_tamper_rejected_before_any_runtime_import(tmp_path):
    path = tmp_path / "plan.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="plan SHA256 mismatch"):
        preflight.checked_plan(path, "0" * 64)


def test_pinned_public_panel_check_rejects_replacing_fixed_example(tmp_path, monkeypatch):
    # This synthetic panel fixture checks overlap, independently of the installed RDKit release.
    monkeypatch.setattr(preflight, "COMMON_RDKIT", rdBase.rdkitVersion)
    paths = []
    for i in range(3):
        path = tmp_path / f"panel{i}.json"
        preflight.write_json(
            path,
            {
                "compounds": [{"id": "q", "smiles": preflight.FIXED_INPUT}],
                "queries": [{"id": str(i), "compound_id": "q"}],
            },
        )
        paths.append(path)
    monkeypatch.setattr(preflight, "PANEL_SHA256", {preflight.digest(p) for p in paths})
    with pytest.raises(ValueError, match="do not replace it"):
        preflight.panel_check(paths)


def test_pinned_public_panel_check_rejects_wrong_runtime_before_reading_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(preflight, "COMMON_RDKIT", rdBase.rdkitVersion + ".fixture-wrong")
    with pytest.raises(ValueError, match="Wrong common RDKit runtime"):
        preflight.panel_check([tmp_path / f"absent-{i}.json" for i in range(3)])


def test_existing_diagnostic_evidence_is_not_modified_on_refused_launch(tmp_path, monkeypatch):
    original = tmp_path / "validation.json"
    original.write_text("previous diagnostic evidence")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "helper",
            "--execute-plan",
            str(tmp_path / "missing-plan"),
            "--expected-plan-sha256",
            "0" * 64,
            "--output",
            str(tmp_path),
        ],
    )
    with pytest.raises(FileNotFoundError):
        preflight.main()
    assert list(tmp_path.iterdir()) == [original]
    assert original.read_text() == "previous diagnostic evidence"


def test_comparison_worker_refuses_existing_results_before_official_import(tmp_path):
    (tmp_path / "comparison-detail.json").write_text("immutable")
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        preflight.comparison_worker({}, tmp_path)


def test_smi_queries_only_approved_gpu_and_records_process_start_identity(monkeypatch):
    import os
    import psutil

    calls = []

    def fake_command(argv, **kwargs):
        calls.append(argv)
        assert argv[argv.index("-i") + 1] in {"0", preflight.GPU_UUID}
        if any(arg.startswith("--query-gpu=") for arg in argv):
            return f"0, {preflight.GPU_UUID}, NVIDIA L40, driver-fixture, bus-fixture\n"
        return f"{os.getpid()}, {preflight.GPU_UUID}, 12\n"

    monkeypatch.setattr(preflight.subprocess, "check_output", fake_command)
    record = preflight.smi_snapshot()
    assert len(calls) == 2
    identity = record["processes"][0]["process_identity"]
    assert identity["pid"] == os.getpid() and identity["start_ticks"] > 0
    assert identity["create_time"] == pytest.approx(psutil.Process().create_time(), abs=0.000001)


def test_busy_gpu_refuses_without_killing_any_process(monkeypatch):
    monkeypatch.setattr(preflight.os, "killpg", lambda *args: pytest.fail("Foreign group must not be killed"))
    preflight.require_idle_gpu({"processes": []})
    with pytest.raises(RuntimeError, match="none killed"):
        preflight.require_idle_gpu({"processes": [{"pid": 1234}]})


def test_prepare_plan_preserves_dangling_output_link(tmp_path):
    path = tmp_path / "plan.json"
    path.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        preflight.prepare_plan(None, None, [], path)
    assert path.is_symlink() and not (tmp_path / "missing").exists()
