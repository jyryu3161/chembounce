"""Model-free checks of the separate, not-yet-executed CUDA diagnostic helper."""

import copy
from functools import wraps
import random
import sys
from types import SimpleNamespace

import numpy as np
import pytest

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
