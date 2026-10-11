#!/usr/bin/env python3
"""Compare the 1.0 dense search kernel with 2.0 packed chunked search.

Synthetic sparse binary fingerprints isolate search cost. This is not a
molecular-quality or end-to-end scaffold-hopping benchmark. Each mode runs in
its own process; timed searches exclude fixture construction and disk loading.
"""

import argparse
import json
from pathlib import Path
import platform
import psutil
from threading import Event, Thread
import subprocess
import sys
import tempfile
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from fingerprint_db import calculate_bulk_tanimoto


def worker(mode, directory, repeats):
    root = Path(directory)
    process = psutil.Process()
    peak = [process.memory_info().rss]
    stop = Event()

    def monitor():
        while not stop.wait(0.002):
            peak[0] = max(peak[0], process.memory_info().rss)

    thread = Thread(target=monitor, daemon=True)
    thread.start()
    data = np.load(root / f"{mode}.npy", mmap_mode="r" if mode == "packed" else None)
    query = np.load(root / "query.npy")
    if mode == "packed":
        query = np.packbits(query)
    times = []
    for _ in range(repeats):
        start = perf_counter()
        if mode == "legacy":
            intersection = np.sum(data & query, axis=1)
            union = np.sum(data | query, axis=1)
            scores = np.divide(intersection, union, out=np.zeros(len(data)), where=union > 0)
            # Canonical tie policy used for result comparison; v1 reversed argsort
            # did not define equal-score ordering.
            indices = np.lexsort((np.arange(len(scores)), -scores))[:100]
            scores = scores[indices]
        else:
            indices, scores = calculate_bulk_tanimoto(query, data, top_n=100, packed=True)
        times.append(perf_counter() - start)
    peak[0] = max(peak[0], process.memory_info().rss)
    stop.set()
    thread.join()
    return {
        "seconds": times,
        "median_seconds": float(np.median(times)),
        "sampled_peak_rss_mb": peak[0] / 1024**2,
        "fingerprint_bytes": int(data.nbytes),
        "indices": indices.tolist(),
        "scores": scores.tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=int, default=100000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", default="benchmark-results/search.json")
    parser.add_argument("--worker", choices=["legacy", "packed"])
    parser.add_argument("--directory")
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker, args.directory, args.repeats)))
        return
    rng = np.random.default_rng(42)
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        dense = (rng.random((args.rows, 2048), dtype=np.float32) < 0.03).astype(np.uint8)
        np.save(root / "legacy.npy", dense)
        np.save(root / "packed.npy", np.packbits(dense, axis=1))
        np.save(root / "query.npy", dense[0])
        del dense
        results = {}
        for mode in ("legacy", "packed"):
            output = subprocess.check_output(
                [
                    sys.executable,
                    __file__,
                    "--worker",
                    mode,
                    "--directory",
                    directory,
                    "--repeats",
                    str(args.repeats),
                ],
                text=True,
            )
            results[mode] = json.loads(output)
        np.testing.assert_array_equal(results["legacy"]["indices"], results["packed"]["indices"])
        np.testing.assert_allclose(results["legacy"]["scores"], results["packed"]["scores"], atol=0, rtol=0)
        for result in results.values():
            del result["indices"], result["scores"]
        report = {
            "description": __doc__,
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "rows": args.rows,
            "bits": 2048,
            "seed": 42,
            "results_equal": True,
            "results": results,
            "speedup": results["legacy"]["median_seconds"] / results["packed"]["median_seconds"],
        }
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
