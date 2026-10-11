import threading
import pytest
from cost_estimation import EstDecoMem


def test_monitor_stops_when_wrapped_function_raises():
    @EstDecoMem
    def broken():
        raise ValueError("expected")

    with pytest.raises(ValueError):
        broken()
    assert not any(t.name == "chembounce-resource-monitor" for t in threading.enumerate())


def test_metrics_have_consistent_units():
    result, cost = EstDecoMem(lambda: 42)()
    assert result == 42
    assert cost["elapsed_time"] >= 0
    assert cost["max_memory_mb"] >= cost["initial_memory_mb"]
    assert cost["memory_samples_count"] >= 2
