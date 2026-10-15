"""Independent control metric recomputation and target-level descriptive figures."""

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
METRICS = ["bit_tanimoto", "count_tanimoto", "hit_shape", "hit_pharmacophore", "score"]
LABELS = ["2D bit", "2D count", "3D shape", "3D features", "3D joint"]
TARGETS = {
    "CHEMBL203": "EGFR",
    "CHEMBL205": "CA2",
    "CHEMBL2147": "PIM1",
    "CHEMBL220": "AChE",
    "CHEMBL251": "A2a",
    "CHEMBL2637": "JNK3",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def estimates(rows, metric):
    pairs = [
        (r["label"] == "active", r["scores"].get(metric))
        for r in rows
        if r["label"] in ("active", "inactive")
    ]
    observed = [(label, value) for label, value in pairs if value is not None]

    def calc(p):
        labels, values = zip(*p) if p else ([], [])
        if len(set(labels)) < 2:
            return {"roc_auc": None, "average_precision": None}
        return {
            "roc_auc": float(roc_auc_score(labels, values)),
            "average_precision": float(average_precision_score(labels, values)),
        }

    result = calc(observed)
    result["observed"] = len(observed)
    result["missing"] = len(pairs) - len(observed)
    extremes = []
    for optimistic in (False, True):
        filled = [
            (label, value if value is not None else (1e6 if label == optimistic else -1e6))
            for label, value in pairs
        ]
        extremes.append(calc(filled))
    for name in ("roc_auc", "average_precision"):
        result[name + "_bounds"] = [e[name] for e in extremes]
    return result


def cluster(values):
    values = np.array(values, dtype=float)
    if not len(values):
        return {"mean": None, "ci95": None, "targets": 0}
    rng = np.random.default_rng(20261016)
    draws = rng.choice(values, (10000, len(values)), replace=True).mean(axis=1)
    return {
        "mean": float(values.mean()),
        "ci95": np.quantile(draws, [0.025, 0.975]).tolist(),
        "targets": len(values),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=BASE / "control_scores.json")
    parser.add_argument("--summary", type=Path, default=BASE / "control_summary.json")
    parser.add_argument("--output", type=Path, default=BASE / "recomputed-controls")
    args = parser.parse_args()
    SOURCE_SCORES, SOURCE_SUMMARY, OUT = args.scores, args.summary, args.output
    OUT.mkdir(exist_ok=True)
    raw = json.loads((SOURCE_SCORES).read_text())
    summary = json.loads((SOURCE_SUMMARY).read_text())
    by_query = defaultdict(list)
    for r in raw:
        by_query[r["query_id"]].append(r)
    errors = []
    query_results = []
    for q in summary["queries"]:
        rows = by_query[q["query_id"]]
        paired = [r for r in rows if all(r["scores"].get(m) is not None for m in METRICS)]
        independent = {m: estimates(rows, m) for m in METRICS}
        for metric in METRICS:
            for name in ["roc_auc", "average_precision", "roc_auc_bounds", "average_precision_bounds"]:
                actual, expected = independent[metric][name], q["discrimination"][metric][name]
                if expected is None:
                    equal = actual is None
                else:
                    equal = np.allclose(actual, expected, rtol=0, atol=1e-12)
                if not equal:
                    errors.append(
                        {
                            "query": q["query_id"],
                            "metric": metric,
                            "statistic": name,
                            "independent": actual,
                            "framework": expected,
                        }
                    )
        query_results.append(
            {
                "query_id": q["query_id"],
                "target_id": q["target_id"],
                "controls": len(rows),
                "failed_3d": q["failed_3d_controls"],
                "prefilter_scaffold_recall": q["prefilter"]["active_scaffold_recall"],
                "active_scaffolds": q["prefilter"]["known_active_scaffolds"],
                "surviving_scaffolds": q["prefilter"]["surviving_active_scaffolds"],
                "independent": independent,
                "complete_case_count": len(paired),
                "paired_complete_case": {m: estimates(paired, m) for m in METRICS},
            }
        )
    targets = []
    for target in TARGETS:
        queries = [q for q in query_results if q["target_id"] == target]
        target_metrics = {}
        for metric in METRICS:
            item = {}
            for name in ["roc_auc", "average_precision"]:
                vals = [
                    q["independent"][metric][name]
                    for q in queries
                    if q["independent"][metric][name] is not None
                ]
                pair_vals = [
                    q["paired_complete_case"][metric][name]
                    for q in queries
                    if q["paired_complete_case"][metric][name] is not None
                ]
                item[name + "_success_conditional"] = float(np.mean(vals)) if vals else None
                item[name + "_paired_complete_case"] = float(np.mean(pair_vals)) if pair_vals else None
                item[name + "_all_control_bounds"] = np.mean(
                    [q["independent"][metric][name + "_bounds"] for q in queries], axis=0
                ).tolist()
            target_metrics[metric] = item
        targets.append(
            {
                "target_id": target,
                "name": TARGETS[target],
                "queries": len(queries),
                "controls": sum(q["controls"] for q in queries),
                "failed_3d": sum(q["failed_3d"] for q in queries),
                "prefilter_scaffold_recall": float(
                    np.mean([q["prefilter_scaffold_recall"] for q in queries])
                ),
                "metrics": target_metrics,
            }
        )
    macro = {}
    for metric in METRICS:
        item = {}
        for name in ["roc_auc", "average_precision"]:
            for population in ["success_conditional", "paired_complete_case"]:
                key = name + "_" + population
                item[key] = cluster(
                    [r["metrics"][metric][key] for r in targets if r["metrics"][metric][key] is not None]
                )
            key = name + "_all_control_bounds"
            item[key] = np.mean([r["metrics"][metric][key] for r in targets], axis=0).tolist()
        macro[metric] = item
    report = {
        "scope": "Retrospective discrimination and prefilter audit on a capped balanced active/inactive control sample; not deployment prevalence or calibrated activity probability",
        "source_sha256": {p.name: digest(p) for p in [SOURCE_SUMMARY, SOURCE_SCORES]},
        "controls": len(raw),
        "observation_unit": "query-control pair",
        "distinct_compound_ids": len({r["compound_id"] for r in raw}),
        "failed_3d_controls": sum(q["failed_3d"] for q in query_results),
        "independent_metric_check": {
            "method": "sklearn ROC AUC and AP, independent missing-score extreme assignment",
            "errors": errors,
            "tolerance": 1e-12,
        },
        "prefilter_target_macro_recall": cluster([t["prefilter_scaffold_recall"] for t in targets]),
        "target_macro_metrics": macro,
        "targets": targets,
        "queries": query_results,
        "uncertainty": "10,000 target bootstrap resamples, seed 20261016; only six independent target groups. Missing-score bounds are identification bounds, not confidence intervals.",
        "paired_complete_case_scope": "All metrics evaluated on the same controls with all requested 3D seeds successful; empty paired query groups remain explicit. Complete-case performance can still be selection-biased.",
    }
    (OUT / "control-analysis.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.5), layout="constrained")
    x = np.arange(len(targets))
    axes[0].bar(x, [100 * t["prefilter_scaffold_recall"] for t in targets], color="#4c78a8")
    for i, t in enumerate(targets):
        vals = [
            100 * q["prefilter_scaffold_recall"] for q in query_results if q["target_id"] == t["target_id"]
        ]
        axes[0].scatter([i - 0.12, i + 0.12][: len(vals)], vals, color="black", s=15, zorder=3)
    axes[0].set(
        xticks=x,
        xticklabels=[t["name"] for t in targets],
        ylabel="Active scaffold recall (%)",
        title="A  Count-Tanimoto prefilter ≥ 0.5",
    )
    axes[0].tick_params(axis="x", rotation=30)
    vals = [macro[m]["roc_auc_paired_complete_case"] for m in METRICS]
    y = np.array([v["mean"] for v in vals])
    ci = np.array([v["ci95"] for v in vals])
    axes[1].errorbar(
        np.arange(len(METRICS)),
        y,
        yerr=np.stack([y - ci[:, 0], ci[:, 1] - y]),
        fmt="o",
        capsize=4,
        color="#4c78a8",
    )
    axes[1].axhline(0.5, color="gray", ls="--", lw=1)
    axes[1].set(
        xticks=np.arange(len(METRICS)),
        xticklabels=LABELS,
        ylim=(0, 1),
        ylabel="Target macro ROC AUC",
        title="B  Paired complete cases; 95% CI",
    )
    axes[1].tick_params(axis="x", rotation=35)
    axes[2].bar(x, [t["failed_3d"] for t in targets], color="#e07b54")
    axes[2].set(
        xticks=x,
        xticklabels=[t["name"] for t in targets],
        ylim=(0, 40),
        ylabel="Failed pairs / 40 per target",
        title="C  At least one 3D seed failed",
    )
    axes[2].tick_params(axis="x", rotation=30)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"controls.{extension}", dpi=200)
    plt.close(fig)
    print(
        json.dumps(
            {
                "controls": len(raw),
                "observation_unit": "query-control pair",
                "distinct_compound_ids": len({r["compound_id"] for r in raw}),
                "failed_3d": report["failed_3d_controls"],
                "metric_mismatches": len(errors),
                "prefilter": report["prefilter_target_macro_recall"],
                "paired_roc": {m: macro[m]["roc_auc_paired_complete_case"] for m in METRICS},
            },
            indent=2,
        )
    )
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
