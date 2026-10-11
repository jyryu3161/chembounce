#!/usr/bin/env python3
"""Reproducible multi-input workflow evaluation; records zero-result cases too."""

import argparse
import hashlib
import json
import multiprocessing
from pathlib import Path
import sys
import time
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
from rdkit import Chem
from chembounce import chembounce


def evaluate(connection, record, options):
    out = Path(options["output_dir"])
    try:
        result, cost = chembounce(record["smiles"], **options)
        molecules = [Chem.MolFromSmiles(s) for s in result["Final structure"]]
        assert all(m is not None and len(Chem.GetMolFrags(m)) == 1 for m in molecules)
        assert result["Final structure"].is_unique
        assert len(result) <= options["overall_max_n"]
        assert result["Tanimoto Similarity"].ge(0.5).all()
        (out / "resource_cost.json").write_text(json.dumps(cost, indent=2) + "\n")
        summary = json.loads((out / "run_summary.json").read_text())
        connection.send(
            {
                "status": summary["status"],
                "candidates": len(result),
                "seconds": cost["elapsed_time"],
                "sampled_peak_rss_mb": cost["max_memory_mb"],
                "shape_failures": summary["counts"].get("shape_failed", 0),
            }
        )
    except Exception:
        connection.send({"status": "evaluation_error", "traceback": traceback.format_exc()})
    finally:
        connection.close()


def bounded_evaluate(record, options, timeout):
    context = multiprocessing.get_context("spawn")
    parent, child = context.Pipe(duplex=False)
    process = context.Process(target=evaluate, args=(child, record, options))
    started = time.perf_counter()
    process.start()
    child.close()
    try:
        if parent.poll(timeout):
            try:
                report = parent.recv()
            except EOFError:
                report = {"status": "worker_crash"}
        else:
            report = {
                "status": "timeout",
                "timeout_seconds": timeout,
                "message": "Incomplete files are preserved for diagnosis; they are not final candidates.",
            }
    finally:
        if process.is_alive():
            process.terminate()
        process.join()
        parent.close()
    report["worker_wall_seconds_including_import"] = time.perf_counter() - started
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fingerprint-db")
    parser.add_argument("--scaffold-db")
    parser.add_argument("--trust-legacy-db", action="store_true")
    parser.add_argument("--output", default="benchmark-results/workflow")
    parser.add_argument("--shape", action="store_true")
    parser.add_argument("--panel", type=Path, default=Path(__file__).with_name("molecules.json"))
    parser.add_argument("--names", nargs="+", help="Optional subset of panel names")
    parser.add_argument(
        "--timeout", type=float, default=300, help="Per-input wall budget, including worker import"
    )
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("Timeout must be positive")
    root = Path(args.output)
    if root.exists():
        parser.error("Output directory exists; choose a new location to preserve prior evidence")
    panel = json.loads(args.panel.read_text())
    if isinstance(panel, dict):
        panel = panel["molecules"]
    if args.names:
        panel = [record for record in panel if record["name"] in args.names]
        if {record["name"] for record in panel} != set(args.names):
            parser.error("Unknown panel name")
    root.mkdir(parents=True)
    protocol = {
        "panel_sha256": hashlib.sha256(args.panel.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "names": [r["name"] for r in panel],
        "timeout_seconds_including_import": args.timeout,
    }
    (root / "evaluation_protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    reports = []
    for record in panel:
        out = root / record["name"]
        options = dict(
            output_dir=out,
            fingerprint_db_file=args.fingerprint_db,
            scaffold_db_file=args.scaffold_db,
            trust_legacy_db=args.trust_legacy_db,
            candidate_thresholds={},
            frag_max_n=10,
            overall_max_n=20,
            scaffold_top_n=20,
            calculate_shape=args.shape,
            seed=42,
            _merge_structure_top_n_=50,
        )
        report = bounded_evaluate(record, options, args.timeout)
        report.update(name=record["name"], manifest=str(out / "run_manifest.json"))
        out.mkdir(exist_ok=True)
        (out / "evaluation_status.json").write_text(json.dumps(report, indent=2) + "\n")
        reports.append(report)
        pd.DataFrame(reports).to_csv(root / "summary.tsv", sep="\t", index=False)
        print(f"{record['name']}: {report['status']}", flush=True)
    print(pd.DataFrame(reports).to_string(index=False))
    if any(r["status"] not in ("success", "no_candidates") for r in reports):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
