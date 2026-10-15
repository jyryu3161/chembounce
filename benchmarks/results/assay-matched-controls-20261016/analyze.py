"""Independently recompute context-matched ROC/AP and missing-score bounds.

No imports from the framework's ranking or aggregation implementation. The panel
was fixed before these chemical scores; this is a follow-up, not an untouched
confirmatory cohort. Similarities are not calibrated activity probabilities.
"""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import statistics
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

METRICS = ["bit_tanimoto", "count_tanimoto", "hit_shape", "hit_pharmacophore", "score"]
NAMES = {
    "CHEMBL203": "EGFR",
    "CHEMBL205": "CA2",
    "CHEMBL2147": "PIM1",
    "CHEMBL220": "AChE",
    "CHEMBL251": "A2a",
    "CHEMBL2637": "JNK3",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_metrics(rows, metric):
    values = [
        (r["label"] == "active", r["scores"].get(metric) if r["status"] == "success" else None) for r in rows
    ]
    observed = [(a, b) for a, b in values if b is not None]

    def compute(pairs):
        if len({a for a, _ in pairs}) != 2:
            return {"roc_auc": None, "average_precision": None}
        labels, scores = zip(*pairs)
        return {
            "roc_auc": float(roc_auc_score(labels, scores)),
            "average_precision": float(average_precision_score(labels, scores)),
        }

    result = compute(observed)
    bounds = []
    for high_is_active in [False, True]:
        bounds.append(
            compute(
                [
                    (label, score if score is not None else (1e6 if label == high_is_active else -1e6))
                    for label, score in values
                ]
            )
        )
    for metric_name in ["roc_auc", "average_precision"]:
        result[metric_name + "_bounds"] = [r[metric_name] for r in bounds]
    result.update(observed=len(observed), missing=len(values) - len(observed))
    return result


def average(children):
    out = {"children": len(children)}
    for name in ["roc_auc", "average_precision"]:
        vals = [r[name] for r in children if r[name] is not None]
        out[name] = statistics.mean(vals) if vals else None
        out[name + "_estimable_children"] = len(vals)
        out[name + "_bounds"] = [statistics.mean(r[name + "_bounds"][i] for r in children) for i in [0, 1]]
    return out


def clustered(values):
    a = np.asarray([v for v in values if v is not None], dtype=float)
    if not len(a):
        return {"mean": None, "ci95": None, "target_groups": 0}
    boot = np.random.default_rng(20261016).choice(a, (10000, len(a)), replace=True).mean(axis=1)
    return {
        "mean": float(a.mean()),
        "ci95": np.quantile(boot, [0.025, 0.975]).tolist(),
        "target_groups": len(a),
    }


def analysis(panel, rows, metric):
    by_query = defaultdict(list)
    for row in rows:
        by_query[row["query_id"]].append(row)
    queries = [
        {
            **{k: q[k] for k in ["context_id", "assay_id", "target_id"]},
            "query_id": q["id"],
            "metrics": rank_metrics(by_query[q["id"]], metric),
        }
        for q in panel["queries"]
    ]

    def group(children, fields):
        buckets = defaultdict(list)
        for child in children:
            buckets[tuple(child[k] for k in fields)].append(child["metrics"])
        return [{**dict(zip(fields, key)), "metrics": average(vals)} for key, vals in sorted(buckets.items())]

    contexts = group(queries, ["target_id", "assay_id", "context_id"])
    assays = group(contexts, ["target_id", "assay_id"])
    targets = group(assays, ["target_id"])
    return {
        "queries": queries,
        "contexts": contexts,
        "assays": assays,
        "targets": targets,
        "target_macro": average([t["metrics"] for t in targets]),
        "target_bootstrap": {
            k: clustered([t["metrics"][k] for t in targets]) for k in ["roc_auc", "average_precision"]
        },
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    folder = args.input
    args.output.mkdir(parents=True, exist_ok=True)
    read = lambda name: json.loads((folder / name).read_text())
    panel, rows, pairs = read("panel.json"), read("scores.json"), read("control_pairs.json")
    if not read("progress.json")["complete"]:
        raise SystemExit("Chemical scoring incomplete; no comparative report")
    lookup = {r["pair_id"]: r for r in rows}
    assert len(lookup) == len(rows) == len(pairs)
    assert set(lookup) == {p["pair_id"] for p in pairs}
    for p in pairs:
        assert all(lookup[p["pair_id"]][k] == v for k, v in p.items())
    paired = [{**r, "scores": r["scores"] if r["all_3d_seeds_succeeded"] else {}} for r in rows]
    summaries = {}
    mismatches = []
    for view, observations, reference in [
        ("all_available", rows, "framework-summary.json"),
        ("paired_complete_case", paired, "framework-paired-summary.json"),
    ]:
        framework = read(reference)
        summaries[view] = {}
        for metric in METRICS:
            independent = analysis(panel, observations, metric)
            summaries[view][metric] = independent
            for level, ids in [
                ("queries", ["query_id"]),
                ("contexts", ["context_id"]),
                ("assays", ["target_id", "assay_id"]),
                ("targets", ["target_id"]),
            ]:
                expected = {tuple(r[k] for k in ids): r["metrics"] for r in framework[metric][level]}
                assert len(expected) == len(independent[level])
                for r in independent[level]:
                    key = tuple(r[k] for k in ids)
                    for name in [
                        "roc_auc",
                        "average_precision",
                        "roc_auc_bounds",
                        "average_precision_bounds",
                    ]:
                        a, b = r["metrics"][name], expected[key][name]
                        equal = (
                            a is None
                            and b is None
                            or a is not None
                            and b is not None
                            and np.allclose(a, b, atol=1e-12, rtol=0)
                        )
                        if not equal:
                            mismatches.append([view, metric, level, list(key), name, a, b])
            for name in ["roc_auc", "average_precision", "roc_auc_bounds", "average_precision_bounds"]:
                a, b = independent["target_macro"][name], framework[metric]["target_macro"][name]
                equal = (
                    a is None
                    and b is None
                    or a is not None
                    and b is not None
                    and np.allclose(a, b, atol=1e-12, rtol=0)
                )
                if not equal:
                    mismatches.append([view, metric, "target_macro", name, a, b])
    differences = {}
    reference = {
        r["target_id"]: r["metrics"]["roc_auc"]
        for r in summaries["paired_complete_case"]["bit_tanimoto"]["targets"]
    }
    for metric in METRICS[1:]:
        deltas = {
            r["target_id"]: r["metrics"]["roc_auc"] - reference[r["target_id"]]
            for r in summaries["paired_complete_case"][metric]["targets"]
            if r["metrics"]["roc_auc"] is not None and reference[r["target_id"]] is not None
        }
        differences[metric] = {"target_differences": deltas, "bootstrap": clustered(list(deltas.values()))}
    counts = defaultdict(Counter)
    for r in rows:
        counts[r["target_id"]]["pairs"] += 1
        counts[r["target_id"]]["failed_3d_pairs"] += not r["all_3d_seeds_succeeded"]
        counts[r["target_id"]]["unexpected_errors"] += r["status"] != "success"
    reasons = Counter(
        t.get("reason", t["status"]) for r in rows for t in r["trials"] if t["status"] != "success"
    )
    records = {r["id"]: r for r in panel["records"]}
    contexts = []
    for c in panel["contexts"]:
        q = next(q for q in panel["queries"] if q["context_id"] == c["id"])
        metadata = c["assay_metadata"]["record"]
        contexts.append(
            {
                "context_id": c["id"],
                "target_id": c["target_id"],
                "assay_id": c["assay_id"],
                "endpoint": c["endpoint"],
                "units": c["units"],
                "bao_format": metadata.get("bao_format"),
                "description": metadata.get("description"),
                "query_smiles": records[q["record_id"]]["smiles"],
            }
        )
    report = {
        "scope": "Follow-up same-assay controls; fixed after heterogeneous-control results, before these chemical scores; no model fitting or threshold tuning",
        "source_sha256": {
            name: digest(folder / name)
            for name in [
                "panel.json",
                "control_pairs.json",
                "scores.json",
                "protocol.json",
                "framework-summary.json",
                "framework-paired-summary.json",
            ]
        },
        "pairs": len(rows),
        "queries": len(panel["queries"]),
        "targets": len(counts),
        "per_target_counts": dict(counts),
        "distinct_control_exact_structures": len({r["control_smiles"] for r in rows}),
        "distinct_control_families": len({r["control_family"] for r in rows}),
        "failed_3d_pairs": sum(not r["all_3d_seeds_succeeded"] for r in rows),
        "failed_seed_reasons": dict(reasons),
        "contexts": contexts,
        "estimates": summaries,
        "paired_roc_differences_from_bit_tanimoto": differences,
        "independent_metric_check": {
            "method": "scikit-learn per-query ranking, separate hierarchical aggregation and extreme missing-score assignment",
            "tolerance": 1e-12,
            "mismatches": mismatches,
        },
        "interpretation": [
            "Same assay ID narrows context but does not independently verify every condition or wild-type sequence.",
            "Scores are not activity probabilities; no activity or synthesis experiments were performed.",
            "Equal query/context/assay/target weights; missing estimates remain null, every pair enters missing-score bounds.",
            "Target bootstrap uses 10,000 draws, seed 20261016; bounds are not confidence intervals.",
            "The same six source targets recur across panels; shared chemistry and related biology can correlate target groups.",
            "Pair counts and seed counts are not independent target sample sizes; label-capped sampling differs from deployment prevalence.",
        ],
    }
    (args.output / "analysis.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    if mismatches:
        raise SystemExit(f"{len(mismatches)} framework disagreement(s) require review")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, fail) = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1.5, 1]})
    labels = ["2D bit", "2D count", "3D shape", "3D features", "3D joint"]
    for i, m in enumerate(METRICS):
        c = summaries["paired_complete_case"][m]["target_bootstrap"]["roc_auc"]
        if c["mean"] is not None:
            ax.errorbar(
                i,
                c["mean"],
                yerr=[[c["mean"] - c["ci95"][0]], [c["ci95"][1] - c["mean"]]],
                fmt="o",
                capsize=4,
                color="#2563a6",
            )
    ax.axhline(0.5, linestyle="--", color="grey", linewidth=1)
    ax.set(
        xticks=range(5),
        xticklabels=labels,
        ylim=(0, 1),
        ylabel="Target-macro ROC AUC",
        title="Paired complete cases; target bootstrap",
    )
    keys = sorted(counts)
    fail.bar([NAMES.get(k, k) for k in keys], [counts[k]["failed_3d_pairs"] for k in keys], color="#b35d3a")
    for i, k in enumerate(keys):
        fail.text(
            i,
            counts[k]["failed_3d_pairs"] + 0.3,
            f"{counts[k]['failed_3d_pairs']}/{counts[k]['pairs']}",
            ha="center",
            fontsize=9,
        )
    fail.set(ylabel="Pairs with a failed 3D seed", title="All selected controls retained")
    fail.set_ylim(0, max([counts[k]["failed_3d_pairs"] for k in keys] + [1]) + 5)
    fig.suptitle(
        f"Same-assay follow-up controls: {len(panel['queries'])} queries, {len(rows)} pairs, {len(counts)} targets"
    )
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(args.output / f"controls.{ext}", dpi=180)
    print(
        json.dumps(
            {
                "pairs": len(rows),
                "failed_3d_pairs": report["failed_3d_pairs"],
                "paired_roc": {
                    m: summaries["paired_complete_case"][m]["target_macro"]["roc_auc"] for m in METRICS
                },
                "metric_mismatches": len(mismatches),
            }
        )
    )


if __name__ == "__main__":
    main()
