"""Present recorded, validated original-630 replay statistics; never recompute estimates."""

import argparse
from collections import Counter
import csv
import io
import math
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import replay
else:
    from . import replay

METRICS = ("attempt_with_hop", "attempt_with_verified_route", "common_3d_output_slot_coverage")
TARGETS = {"scaffold": 6, "time": 6, "target": 3}
TITLES = ("Attempt with scaffold hop", "Attempt with rule-replayed route", "Common 3D output-slot coverage")
PRIMARY_SCOPE = (
    "Recorded target-mean contrast; success-conditional populations may differ, not like-for-like improvement"
)
COMPLETE_SCOPE = "Secondary selection-biased contrast; both outcomes available at matching query/seed slots"
SCOPE = (
    "Original frozen generation study, commit d972124d587e647ef48455cf1a038e6be96e5b04 (39 files). "
    "Post hoc descriptive presentation of recorded statistics, not preregistration. "
    "No new estimates, pooling, imputation, generation, models, chemistry or 3D scoring. "
    "Repeated split targets are not independent cohorts. Rule replay is not experimental synthesis "
    "or the separate common-forward re-inference audit. Coverage is successful common-3D output slots "
    "divided by cap 3, not conditional quality. No activity, yield, probability or superiority claim. "
    "The current-pipeline supplementary 270 study is excluded. " + PRIMARY_SCOPE + "."
)
require, same, parse, sha = replay.require, replay.same, replay.parse, replay.sha


def count(value, maximum=None):
    require(type(value) is int and value >= 0, "Expected nonnegative integer count")
    require(maximum is None or value <= maximum, "Count exceeds planned population")


def number(value, domain=None):
    if value is None:
        return
    require(type(value) in (int, float) and math.isfinite(value), "Expected finite number or null")
    require(domain is None or domain[0] <= value <= domain[1], "Metric outside its domain")


def interval(item, maximum, domain=None):
    require(isinstance(item, dict) and set(item) == {"mean", "ci95", "independent_targets"}, "Malformed CI")
    number(item["mean"], domain)
    n, ci = item["independent_targets"], item["ci95"]
    count(n, maximum)
    require((item["mean"] is None) == (n == 0), "Missing estimate/target count disagreement")
    if n < 2:
        require(ci is None, "CI requires at least two targets")
    else:
        require(isinstance(ci, list) and len(ci) == 2 and None not in ci, "Malformed CI endpoints")
        for x in ci:
            number(x, domain)
        require(ci[0] <= ci[1], "Reversed CI")


def summary(item, queries, targets, names=None):
    require(
        all(isinstance(item.get(k), dict) for k in ("per_query", "per_target", "target_macro")),
        "Missing summary",
    )
    names = set(item["target_macro"]) if names is None else set(names)
    require(names and set(item["target_macro"]) == names, "Metric membership differs")
    require(
        set(item["per_query"]) == queries and set(item["per_target"]) == targets, "Summary population differs"
    )
    for level in ("per_query", "per_target"):
        for values in item[level].values():
            require(isinstance(values, dict) and set(values) == names, "Metric row malformed")
            for metric, value in values.items():
                number(value, (0, 1) if metric in METRICS else None)
    for metric, value in item["target_macro"].items():
        interval(value, len(targets), (0, 1) if metric in METRICS else None)
        same(
            value["independent_targets"],
            sum(v[metric] is not None for v in item["per_target"].values()),
            "Target availability",
        )
    return names


def bounds(item, queries, targets, paired=False):
    domain = (-1, 1) if paired else (0, 1)
    require(
        set(item["per_query"]) == queries and set(item["per_target"]) == targets, "Bound population differs"
    )
    groups = [*item["per_query"].values(), *item["per_target"].values(), item["target_macro"]]
    for values in groups:
        require(set(values) == set(METRICS), "Bound metric membership differs")
        for value in values.values():
            require(set(value) == {"lower", "upper"}, "Malformed identification bounds")
            lo, hi = value["lower"], value["upper"]
            number(lo, domain)
            number(hi, domain)
            require((lo is None) == (hi is None), "Incomplete identification bounds")
            require(lo is None or lo <= hi, "Reversed identification bounds")


def availability(item, queries, targets, attempts):
    require(
        set(item)
        == {
            "attempts",
            "queries_with_observed_attempts",
            "targets_with_observed_attempts",
            "per_query_attempts",
            "per_target_attempts",
        },
        "Malformed availability",
    )
    same(item["attempts"], attempts, "Availability count")
    for level, population, key in (
        ("per_query_attempts", queries, "queries_with_observed_attempts"),
        ("per_target_attempts", targets, "targets_with_observed_attempts"),
    ):
        require(set(item[level]) == population, "Availability population differs")
        for n in item[level].values():
            count(n, attempts)
        same(sum(item[level].values()), attempts, "Availability total")
        same(item[key], sum(n > 0 for n in item[level].values()), "Available units")


def status_counts(item, attempts):
    require(isinstance(item, dict) and set(item) <= replay.STATUSES, "Unknown attempt status")
    for n in item.values():
        count(n, attempts)
    same(sum(item.values()), attempts, "Status denominator")


def validate_report(split, report, pins, literals):
    require(
        set(report)
        == {
            "split",
            "grid_audit",
            "source_sha256",
            "dataset_sha256",
            "methods",
            "paired_pipeline_minus_comparator",
            "infrastructure_missing_outcome_sensitivity",
            "complete_case_secondary",
            *literals,
        },
        "Comparison fields differ",
    )
    same(report["split"], split, "Split")
    for key, value in literals.items():
        same(report[key], value, key)
    total, nq = replay.SPLITS[split]
    grid = report["grid_audit"]
    require(
        set(grid)
        == {
            "split",
            "expected_runs",
            "observed_runs",
            "missing",
            "unexpected",
            "duplicates",
            "status_counts",
            "complete",
        },
        "Malformed grid audit",
    )
    for key, value in {
        "split": split,
        "expected_runs": total,
        "observed_runs": total,
        "missing": [],
        "unexpected": [],
        "duplicates": [],
        "complete": True,
    }.items():
        same(grid[key], value, "Complete grid." + key)
    status_counts(grid["status_counts"], total)
    names = replay.paths(split)
    same(report["dataset_sha256"], pins[names["panel"]], "Dataset SHA")
    sources = {
        names[f"{r}_{g}"]: pins[names[f"{r}_{g}"]] for g in replay.GROUPS for r in ("protocol", "runs")
    }
    sources[names["common_generated_scores"]] = pins[names["common_generated_scores"]]
    same(report["source_sha256"], sources, "Recorded source hashes")
    methods = {
        m
        for g in replay.GROUPS
        for m in (replay.METHODS[g] if split == "scaffold" or g == "baselines" else replay.METHODS[g][:1])
    }
    require(set(report["methods"]) == methods, "Original method membership differs")
    query_ids = set(report["methods"]["pipeline"]["per_query"])
    target_ids = set(report["methods"]["pipeline"]["per_target"])
    require(
        len(query_ids) == nq and len(target_ids) == TARGETS[split], "Original query/target population differs"
    )
    require(all(isinstance(x, str) and x for x in query_ids | target_ids), "Invalid population ID")
    infra, complete = report["infrastructure_missing_outcome_sensitivity"], report["complete_case_secondary"]
    require(
        set(infra)
        == {
            "interpretation",
            "bounded_metric_domains",
            "hierarchy",
            "methods",
            "paired_pipeline_minus_comparator",
            "other_metrics",
        },
        "Malformed infrastructure sensitivity",
    )
    require(
        set(complete) == {"interpretation", "methods", "paired_pipeline_minus_comparator_bounded_metrics"},
        "Malformed complete-case analysis",
    )
    same(infra["bounded_metric_domains"], dict.fromkeys(METRICS, [0, 1]), "Bound domains")
    require(set(infra["methods"]) == set(complete["methods"]) == methods, "Sensitivity methods differ")
    statuses = Counter()
    for method in sorted(methods):
        main, missing, available = (
            report["methods"][method],
            infra["methods"][method],
            complete["methods"][method],
        )
        same(main["attempts"], nq * 3, "Method grid")
        status_counts(main["status_counts"], main["attempts"])
        statuses.update(main["status_counts"])
        metrics = summary(main, query_ids, target_ids)
        require(set(METRICS) <= metrics, "Required figure metrics missing")
        lost = main["infrastructure_interruption_attempts"]
        same(lost, main["status_counts"].get("infrastructure_interruption", 0), "Infrastructure count")
        for key in ("runtime_observed_attempts", "runtime_missing_attempts"):
            count(main[key], main["attempts"])
        same(
            main["runtime_observed_attempts"] + main["runtime_missing_attempts"],
            main["attempts"],
            "Runtime denominator",
        )
        for n in main["runtime_measurement_status_counts"].values():
            count(n)
        same(
            sum(main["runtime_measurement_status_counts"].values()), main["attempts"], "Runtime status total"
        )
        availability(
            main["runtime_observation_availability"], query_ids, target_ids, main["runtime_observed_attempts"]
        )
        for key, value in {
            "attempts": nq * 3,
            "infrastructure_missing_attempts": lost,
            "observed_outcome_attempts": nq * 3 - lost,
        }.items():
            same(missing[key], value, "Infrastructure." + key)
        availability(missing["full_grid_availability"], query_ids, target_ids, nq * 3)
        require(
            all(n == 3 for n in missing["full_grid_availability"]["per_query_attempts"].values()),
            "Seed grid differs",
        )
        availability(missing["observed_outcome_availability"], query_ids, target_ids, nq * 3 - lost)
        bounds(missing, query_ids, target_ids)
        same(available["attempts"], nq * 3 - lost, "Complete-case count")
        same(available["excluded_infrastructure_attempts"], lost, "Complete-case exclusion")
        status_counts(available["status_counts"], available["attempts"])
        require(
            "infrastructure_interruption" not in available["status_counts"],
            "Complete case includes infrastructure loss",
        )
        availability(available["availability"], query_ids, target_ids, available["attempts"])
        count(available["runtime_observations"]["attempts"], available["attempts"])
        availability(
            available["runtime_observations"],
            query_ids,
            target_ids,
            available["runtime_observations"]["attempts"],
        )
        summary(available, query_ids, target_ids, metrics)
    same(dict(statuses), grid["status_counts"], "Overall status counts")
    comparators = methods - {"pipeline"}
    require(
        set(report["paired_pipeline_minus_comparator"])
        == set(infra["paired_pipeline_minus_comparator"])
        == set(complete["paired_pipeline_minus_comparator_bounded_metrics"])
        == comparators,
        "Paired method membership differs",
    )
    for method in sorted(comparators):
        paired = report["paired_pipeline_minus_comparator"][method]
        require(
            set(paired) == set(report["methods"]["pipeline"]["target_macro"]),
            "Paired metric membership differs",
        )
        for metric, value in paired.items():
            interval(value, len(target_ids), (-1, 1) if metric in METRICS else None)
        missing = infra["paired_pipeline_minus_comparator"][method]
        same(missing["paired_attempts"], nq * 3, "Paired grid")
        count(missing["pairs_with_infrastructure_missing_outcome"], nq * 3)
        count(
            missing["pairs_with_both_outcomes_missing"], missing["pairs_with_infrastructure_missing_outcome"]
        )
        bounds(missing, query_ids, target_ids, paired=True)
        available = complete["paired_pipeline_minus_comparator_bounded_metrics"][method]
        same(
            available["excluded_infrastructure_pairs"],
            missing["pairs_with_infrastructure_missing_outcome"],
            "Paired exclusions",
        )
        same(
            available["paired_attempts"],
            nq * 3 - available["excluded_infrastructure_pairs"],
            "Available pairs",
        )
        availability(available["availability"], query_ids, target_ids, available["paired_attempts"])
        # This summary contains signed differences, so validate its scalar values separately.
        require(
            set(available["per_query"]) == query_ids and set(available["per_target"]) == target_ids,
            "Paired population differs",
        )
        for level in ("per_query", "per_target"):
            for values in available[level].values():
                require(set(values) == set(METRICS), "Paired bounded metrics differ")
                for value in values.values():
                    number(value, (-1, 1))
        require(set(available["target_macro"]) == set(METRICS), "Paired macro metrics differ")
        for metric, value in available["target_macro"].items():
            interval(value, len(target_ids), (-1, 1))
            same(
                value["independent_targets"],
                sum(v[metric] is not None for v in available["per_target"].values()),
                "Paired target availability",
            )


def _load_completed(directory, expected_validation_sha256):
    directory = replay.no_symlinks(directory)
    names = {"validation.json", *(f"{s}-comparison.json" for s in replay.SPLITS)}
    require(
        directory.is_dir() and {p.name for p in directory.iterdir()} == names,
        "Require exactly three comparisons and validation",
    )
    paths = {n: directory / n for n in sorted(names)}
    before = replay.snapshot(paths)
    same(
        sha(before["validation.json"]),
        replay.sha_value(expected_validation_sha256),
        "External validation SHA",
    )
    validation = parse(before["validation.json"])
    require(
        set(validation)
        == {
            "status",
            "attempts",
            "absolute_tolerance",
            "manifest_sha256",
            "historical_helper_sha256",
            "runtime",
            "scope",
            "input_before_after_sha256",
            "recomputed_sha256",
        },
        "Replay validation fields differ",
    )
    for key, value in {
        "status": "passed",
        "attempts": 630,
        "absolute_tolerance": 1e-12,
        "historical_helper_sha256": replay.HELPER_SHA,
        "scope": replay.SCOPE,
    }.items():
        same(validation[key], value, "Replay validation." + key)
    replay.sha_value(validation["manifest_sha256"])
    require(
        set(validation["runtime"]) == {"python", "numpy"}
        and all(isinstance(v, str) and v for v in validation["runtime"].values()),
        "Replay runtime missing",
    )
    pins = validation["input_before_after_sha256"]
    required = replay.input_names() | {replay.HELPER, "parent-review.json", "manifest.json"}
    require(
        required < set(pins) and all(n in required or n.startswith("parent-review-evidence/") for n in pins),
        "Replay provenance membership differs",
    )
    for name, value in pins.items():
        replay.safe_relative(name)
        replay.sha_value(value)
    for name, value in {
        replay.HELPER: replay.HELPER_SHA,
        replay.ATTRIBUTION: replay.ATTRIBUTION_SHA,
        replay.PROVENANCE: replay.PROVENANCE_SHA,
        "manifest.json": validation["manifest_sha256"],
    }.items():
        same(pins[name], value, "Replay provenance." + name)
    require(set(validation["recomputed_sha256"]) == set(replay.SPLITS), "Missing split hash")
    helper = replay.read_bytes(replay.HERE / replay.HELPER)
    same(sha(helper), replay.HELPER_SHA, "Local archived helper SHA")
    literals = replay.report_text(helper)
    reports = {}
    for split in replay.SPLITS:
        raw = before[f"{split}-comparison.json"]
        same(sha(raw), replay.sha_value(validation["recomputed_sha256"][split]), "Recomputed report SHA")
        reports[split] = parse(raw)
        validate_report(split, reports[split], pins, literals)
    return paths, before, reports


def load_completed(directory, expected_validation_sha256):
    try:
        return _load_completed(directory, expected_validation_sha256)
    except (KeyError, TypeError, IndexError, AttributeError) as exc:
        raise ValueError(f"Malformed recorded replay schema: {exc}") from exc


def metric_scope(metric):
    if "success_conditional" in metric:
        return "Success-conditional score; populations can differ across methods"
    if metric == "wall_seconds":
        return "Observed runtimes only; missing runtime excluded, not zero"
    return "Recorded conservative realized output"


def table_rows(reports):
    tables = {"metrics": [], "paired_contrasts": [], "infrastructure_bounds": [], "availability": []}

    def counts(split, population, method, item, path="", query="", target=""):
        for key, value in sorted(item.items()):
            label = f"{path}.{key}" if path else key
            if isinstance(value, dict):
                for name, n in sorted(value.items()):
                    counts(
                        split,
                        population,
                        method,
                        {name: n},
                        label,
                        name if key == "per_query_attempts" else query,
                        name if key == "per_target_attempts" else target,
                    )
            else:
                tables["availability"].append(
                    dict(
                        split=split,
                        population=population,
                        method=method,
                        query_id=query,
                        target_id=target,
                        recorded_field=label,
                        count=value,
                    )
                )

    for split, report in reports.items():
        infra, complete = (
            report["infrastructure_missing_outcome_sensitivity"],
            report["complete_case_secondary"],
        )
        for population, group in (
            ("realized_all_attempts", report["methods"]),
            ("complete_case_secondary", complete["methods"]),
        ):
            for method, item in sorted(group.items()):
                available = (
                    infra["methods"][method]["full_grid_availability"]
                    if population == "realized_all_attempts"
                    else item["availability"]
                )
                for level in ("per_query", "per_target", "target_macro"):
                    entries = {"": item[level]} if level == "target_macro" else item[level]
                    for identity, metrics in sorted(entries.items()):
                        attempts = (
                            item["attempts"]
                            if level == "target_macro"
                            else available[level + "_attempts"][identity]
                        )
                        for metric, value in sorted(metrics.items()):
                            stat = (
                                value
                                if level == "target_macro"
                                else {"mean": value, "ci95": None, "independent_targets": None}
                            )
                            ci = stat["ci95"] or [None, None]
                            tables["metrics"].append(
                                dict(
                                    split=split,
                                    population=population,
                                    level=level,
                                    method=method,
                                    query_id=identity if level == "per_query" else "",
                                    target_id=identity if level == "per_target" else "",
                                    metric=metric,
                                    value=stat["mean"],
                                    ci95_low=ci[0],
                                    ci95_high=ci[1],
                                    independent_targets=stat["independent_targets"],
                                    attempts=attempts,
                                    scope=metric_scope(metric)
                                    if population == "realized_all_attempts"
                                    else "Secondary selection-biased available cases; "
                                    + metric_scope(metric),
                                )
                            )
                counts(
                    split,
                    population,
                    method,
                    {k: v for k, v in item.items() if k not in ("per_query", "per_target", "target_macro")},
                )
        for population, group, scope in (
            ("recorded_target_contrast", report["paired_pipeline_minus_comparator"], PRIMARY_SCOPE),
            (
                "complete_case_paired_secondary",
                complete["paired_pipeline_minus_comparator_bounded_metrics"],
                COMPLETE_SCOPE,
            ),
        ):
            for method, item in sorted(group.items()):
                groups = (
                    {"target_macro": {"": item}}
                    if population == "recorded_target_contrast"
                    else {
                        "per_query": item["per_query"],
                        "per_target": item["per_target"],
                        "target_macro": {"": item["target_macro"]},
                    }
                )
                for level, entries in groups.items():
                    for identity, values in sorted(entries.items()):
                        for metric, value in sorted(values.items()):
                            stat = (
                                value
                                if level == "target_macro"
                                else {"mean": value, "ci95": None, "independent_targets": None}
                            )
                            ci = stat["ci95"] or [None, None]
                            tables["paired_contrasts"].append(
                                dict(
                                    split=split,
                                    population=population,
                                    level=level,
                                    comparator=method,
                                    query_id=identity if level == "per_query" else "",
                                    target_id=identity if level == "per_target" else "",
                                    metric=metric,
                                    pipeline_minus_comparator=stat["mean"],
                                    ci95_low=ci[0],
                                    ci95_high=ci[1],
                                    independent_targets=stat["independent_targets"],
                                    paired_attempts=item.get("paired_attempts")
                                    if level == "target_macro"
                                    else item["availability"][level + "_attempts"][identity],
                                    contrast_scope=scope,
                                )
                            )
                if population != "recorded_target_contrast":
                    counts(
                        split,
                        population,
                        method,
                        {
                            k: v
                            for k, v in item.items()
                            if k not in ("per_query", "per_target", "target_macro")
                        },
                    )
        for population, group in (
            ("method_missing_outcome", infra["methods"]),
            ("paired_pipeline_minus_comparator", infra["paired_pipeline_minus_comparator"]),
        ):
            for method, item in sorted(group.items()):
                for level in ("per_query", "per_target", "target_macro"):
                    entries = {"": item[level]} if level == "target_macro" else item[level]
                    for identity, values in sorted(entries.items()):
                        for metric, bound in sorted(values.items()):
                            tables["infrastructure_bounds"].append(
                                dict(
                                    split=split,
                                    population=population,
                                    level=level,
                                    method_or_comparator=method,
                                    query_id=identity if level == "per_query" else "",
                                    target_id=identity if level == "per_target" else "",
                                    metric=metric,
                                    lower=bound["lower"],
                                    upper=bound["upper"],
                                    interval_type="Identification bounds, NOT confidence intervals",
                                )
                            )
                counts(
                    split,
                    population,
                    method,
                    {k: v for k, v in item.items() if k not in ("per_query", "per_target", "target_macro")},
                )
    return tables


def csv_bytes(rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)  # None becomes blank; Python's round-trip numeric representation is preserved.
    return stream.getvalue().encode("utf-8")


def draw(split, report, output, synthetic_fixture=False):
    import matplotlib

    matplotlib.use("Agg", force=True)
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D

    methods = sorted(report["methods"])
    fig, axes = plt.subplots(1, 3, figsize=(17, max(6.8, len(methods) * 0.54 + 3.1)), sharey=True)
    fig.subplots_adjust(left=0.27, right=0.95, bottom=0.28, top=0.80, wspace=0.35)
    fig.suptitle(
        f"Original frozen generation study — {split} split\n"
        f"{report['grid_audit']['observed_runs']} attempts; {TARGETS[split]} independent targets; "
        "frozen d972124 (39 files)",
        fontsize=14,
        y=0.97,
    )
    if synthetic_fixture:
        fig.text(
            0.5, 0.87, "SYNTHETIC FIXTURE — NOT STUDY RESULTS", ha="center", color="#9c2020", fontsize=12
        )
    for ax, metric, title in zip(axes, METRICS, TITLES):
        for y, method in enumerate(methods):
            stat = report["methods"][method]["target_macro"][metric]
            bound = report["infrastructure_missing_outcome_sensitivity"]["methods"][method]["target_macro"][
                metric
            ]
            if stat["ci95"] is not None:
                ax.hlines(y - 0.10, *stat["ci95"], color="#174f83", linewidth=2)
            if stat["mean"] is not None:
                ax.plot(stat["mean"], y - 0.10, "o", color="#174f83", markersize=5)
            else:
                ax.text(0.015, y - 0.10, "missing", color="#666666", va="center", fontsize=9)
            if bound["lower"] is not None:
                ax.hlines(
                    y + 0.15, bound["lower"], bound["upper"], color="#a64b00", linestyle=":", linewidth=1.8
                )
                ax.plot(
                    [bound["lower"], bound["upper"]], [y + 0.15, y + 0.15], "D", color="#a64b00", markersize=3
                )
            ax.text(
                1.02,
                y,
                f"n={stat['independent_targets']}",
                transform=ax.get_yaxis_transform(),
                va="center",
                fontsize=8,
            )
        ax.set_title(title, fontsize=11, pad=14)
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(len(methods) - 0.45, -0.65)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_xlabel("Target-macro mean (0–1)")
        ax.grid(axis="x", color="#dddddd", linewidth=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(
        range(len(methods)),
        [
            f"{m}\nA={report['methods'][m]['attempts']}; infra={report['methods'][m]['infrastructure_interruption_attempts']}"
            for m in methods
        ],
        fontsize=9,
    )
    handles = [
        Line2D([0], [0], color="#174f83", marker="o", label="Target-macro mean and target-bootstrap 95% CI"),
        Line2D(
            [0],
            [0],
            color="#a64b00",
            marker="D",
            linestyle=":",
            label="Infrastructure identification bounds — NOT a CI",
        ),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.115), frameon=False, fontsize=10)
    fig.text(
        0.5,
        0.065,
        "A: all attempts; infra: missing infrastructure outcomes; n: independent targets with this metric.\n"
        "Rule replay ≠ experimental synthesis or common-forward audit. 3D coverage = available successful output slots / cap 3.\n"
        "Descriptive unadjusted CIs; overlap is not a test. Post hoc presentation; repeated split targets are not independent cohorts.",
        ha="center",
        va="center",
        fontsize=9,
    )
    try:
        fig.savefig(output / f"{split}.png", dpi=180, metadata={"Description": SCOPE})
        fig.savefig(
            output / f"{split}.pdf",
            metadata={
                "Title": f"Original frozen generation study: {split}",
                "Subject": SCOPE,
                "CreationDate": None,
                "ModDate": None,
            },
        )
    finally:
        plt.close(fig)


def export(directory, expected_validation_sha256, output, synthetic_fixture=False):
    directory, output = replay.no_symlinks(directory), replay.no_symlinks(output)
    require(not output.exists(), "Output already exists")
    require(
        not output.is_relative_to(directory) and not directory.is_relative_to(output),
        "Input/output trees must be separate",
    )
    paths, before, reports = load_completed(directory, expected_validation_sha256)
    # All gates and table projection complete before importing Matplotlib or creating output.
    rows = table_rows(reports)
    tables = {name + ".csv": csv_bytes(values) for name, values in rows.items()}
    replay.unchanged(paths, before)
    require({p.name for p in directory.iterdir()} == set(paths), "Replay directory membership changed")
    import matplotlib
    import numpy as np

    output = replay.new_directory(output)
    for name, content in tables.items():
        replay.write_new(output / name, content)
    for split, report in reports.items():
        draw(split, report, output, synthetic_fixture=synthetic_fixture)
    replay.unchanged(paths, before)
    require({p.name for p in directory.iterdir()} == set(paths), "Replay directory membership changed")
    files = {p.name: sha(replay.read_bytes(p)) for p in sorted(output.iterdir())}
    manifest = {
        "schema_version": 1,
        "status": "completed",
        "synthetic_fixture": synthetic_fixture,
        "scope": SCOPE,
        "original_source_commit": replay.COMMIT,
        "original_scientific_files": 39,
        "expected_validation_sha256": expected_validation_sha256,
        "input_before_after_sha256": {n: sha(raw) for n, raw in before.items()},
        "exporter_sha256": sha(replay.read_bytes(__file__)),
        "replay_contract_sha256": sha(replay.read_bytes(replay.__file__)),
        "historical_helper_sha256": replay.HELPER_SHA,
        "runtime": {
            "python": sys.version,
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
        },
        "method_order": {s: sorted(r["methods"]) for s, r in reports.items()},
        "independent_targets_by_split": TARGETS,
        "plotted_metrics": list(METRICS),
        "table_rows": {name: len(values) for name, values in rows.items()},
        "null_encoding": "Blank CSV cell; missing plot values are labelled missing, never zero",
        "paired_contrast_scope": PRIMARY_SCOPE,
        "files_sha256": files,
    }
    replay.write_new(output / "manifest.json", replay.encoded(manifest))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--expected-validation-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--synthetic-fixture",
        action="store_true",
        help="Watermark synthetic diagnostics; not an authenticity check",
    )
    args = parser.parse_args()
    result = export(args.replay, args.expected_validation_sha256, args.output, args.synthetic_fixture)
    print(
        f"{result['status']}: {len(result['files_sha256'])} presentation files; no new statistical estimates"
    )


if __name__ == "__main__":
    main()
