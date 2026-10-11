"""Process-local timing and sampled RSS, with constant-memory monitoring."""

from functools import wraps
from threading import Event, Thread
from time import perf_counter
import psutil


class EstDecoMem:
    """Return (wrapped result, metrics); exceptions propagate and stop monitoring.

    CPU percent uses process CPU seconds / wall seconds (100% = one CPU).
    RSS is sampled, not an exact peak; it excludes independent child processes.
    """

    def __init__(self, func):
        self.func = func
        wraps(func)(self)

    def __call__(self, *args, **kwargs):
        process = psutil.Process()
        start = perf_counter()
        cpu_before = process.cpu_times()
        initial = process.memory_info().rss / 1024**2
        stats = {"sum": initial, "max": initial, "count": 1}
        stop = Event()

        def sample():
            rss = process.memory_info().rss / 1024**2
            stats["sum"] += rss
            stats["max"] = max(stats["max"], rss)
            stats["count"] += 1
            return rss

        def monitor():
            while not stop.wait(0.01):
                try:
                    sample()
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    return

        thread = Thread(target=monitor, daemon=True, name="chembounce-resource-monitor")
        thread.start()
        try:
            result = self.func(*args, **kwargs)
        finally:
            stop.set()
            thread.join()
        elapsed = perf_counter() - start
        cpu_after = process.cpu_times()
        user, system = cpu_after.user - cpu_before.user, cpu_after.system - cpu_before.system
        final = sample()
        return result, {
            "elapsed_time": elapsed,
            "user_cpu_time": user,
            "system_cpu_time": system,
            "total_cpu_time": user + system,
            "cpu_usage_percent": 100 * (user + system) / elapsed if elapsed else 0,
            "cpu_count": psutil.cpu_count(),
            "initial_memory_mb": initial,
            "final_memory_mb": final,
            "max_memory_mb": stats["max"],
            "avg_memory_mb": stats["sum"] / stats["count"],
            "memory_delta_mb": final - initial,
            "memory_samples_count": stats["count"],
            "memory_sampling_interval_seconds": 0.01,
            "memory_scope": "process RSS; sampled peak",
        }


class EstDeco(EstDecoMem):
    """Backward-compatible resource decorator."""
