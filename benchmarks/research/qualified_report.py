"""Offline HTML review of a trusted captured qualified benchmark collection.

``render_report(plan_bytes, collection_bytes, *, collection_sha256)`` validates
and recomputes with qualified_analysis.analyze, then returns a self-contained HTML
string. The digest must come from a trusted external record, not the input itself.
``write_report`` provides the same contract for paths and refuses to overwrite.
Neither API reads sidecars, runs chemistry, or determines live campaign state.
"""

import argparse
import hashlib
from html import escape
import json
from pathlib import Path

from benchmarks.research import qualified_analysis as analysis


ENDPOINT_LABELS = (
    "At least one family",
    "Unique families at 3 / 3",
    "At least three families",
)


def _text(value):
    return escape(str(value), quote=True)


def _table(headers, rows, caption):
    return (
        '<div class="table-wrap"><table><caption>'
        + _text(caption)
        + "</caption><thead><tr>"
        + "".join(f'<th scope="col">{_text(h)}</th>' for h in headers)
        + "</tr></thead><tbody>"
        + "".join("<tr>" + "".join(f"<td>{_text(v)}</td>" for v in row) + "</tr>" for row in rows)
        + "</tbody></table></div>"
    )


def _number(value):
    return "unknown" if value is None else format(value, ".8g")


def _estimate(value):
    return "unavailable" if value is None else f"{value['fraction']} ({value['value']:.2%})"


def _details(title, body):
    return f"<details><summary>{_text(title)}</summary>{body}</details>"


def render_report(plan_bytes: bytes, collection_bytes: bytes, *, collection_sha256: str) -> str:
    """Return escaped static HTML after digest verification and authoritative analysis.

    Only collector-consistent computational outcomes are summarized. A digest
    authenticates bytes relative to the caller's trusted record; it does not prove
    chemical validity. Invalid inputs raise ValueError; no output is written here.
    """
    actual_digest = hashlib.sha256(collection_bytes).hexdigest()
    if collection_sha256 != actual_digest:
        raise ValueError("Collection byte digest mismatch")
    plan, collection = json.loads(plan_bytes), json.loads(collection_bytes)
    report = analysis.analyze(plan, collection)
    recorded = {tuple(row[key] for key in analysis.SLOT_KEYS): row for row in collection["attempts"]}
    planned, observed = report["planned_attempts"], report["reported_attempts"]
    missing = planned - observed
    parts = [
        '<!doctype html><html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        '<meta http-equiv="Content-Security-Policy" '
        "content=\"default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'\">",
        "<title>Qualified benchmark review</title><style>",
        "body{font:16px/1.5 system-ui,sans-serif;max-width:1440px;margin:2rem auto;padding:0 1rem;"
        "color:#17212b;background:#fafbfc}h1,h2{line-height:1.2}h2{margin-top:2.2rem}"
        ".notice{padding:1rem;background:#eaf1f8;border-left:5px solid #365f85}"
        ".table-wrap{overflow-x:auto;margin:1rem 0}table{border-collapse:collapse;width:100%;background:white}"
        "caption{text-align:left;font-weight:600;padding:.5rem 0}th,td{border:1px solid #cbd3dc;"
        "padding:.45rem .65rem;text-align:left;vertical-align:top;overflow-wrap:anywhere}"
        "th{background:#edf1f5}details{margin:1rem 0}summary{cursor:pointer;font-weight:600}"
        "pre{white-space:pre-wrap;overflow-wrap:anywhere}code{overflow-wrap:anywhere}"
        "@media print{body{max-width:none}.table-wrap{overflow:visible}}",
        "</style></head><body><main><h1>Qualified benchmark review</h1>",
        f"<p>Run: <strong>{_text(report['run_id'])}</strong></p>",
        '<p class="notice"><strong>'
        + ("Partial saved snapshot" if missing else "All planned slots have records")
        + f"</strong>: {observed} recorded / {planned} planned; {missing} missing. "
        "Record completeness does not assert operational campaign completion. "
        "This report describes the supplied snapshot, not live process state.</p>",
        "<p>Trust boundary: externally pinned collector bytes and analyzer-consistent results. "
        "No fresh chemical validation, synthesis experiment, activity evidence, or efficacy claim. "
        "Zero candidates do not establish chemical infeasibility; no significance or winner is inferred.</p>",
        "<h2>Method comparison</h2>",
        f"<p>Planned grid: {len(plan['methods'])} methods × {len(plan['queries'])} target/query pairs "
        f"× {len(plan['seeds'])} seeds = {planned} slots across {report['target_count']} targets. "
        "Exact fractions below are hierarchical estimates, not pooled success counts. "
        "Equal weight to seeds within each query, queries within each target, then targets. "
        "Every planned slot remains in the denominator; missing and failed attempts contribute zero.</p>",
        "<p>Endpoints deduplicate selected outputs by chemical family within each attempt, keeping "
        "the first-ranked representative. Three exact molecular states need not be three families. "
        "No refill occurs after deduplication; the family-count endpoint always divides by three.</p>",
        _table(
            ["Method", "Planned slots", *ENDPOINT_LABELS, "Duplicate-family outputs"],
            [
                [
                    method,
                    data["planned_attempts"],
                    *[_estimate(data["estimates"][key]["estimate"]) for key in analysis.ENDPOINTS],
                    data["duplicate_family_outputs"],
                ]
                for method, data in report["methods"].items()
            ],
            "Descriptive endpoints over the full planned grid",
        ),
        "<h2>Slot status</h2><p>no_candidates is a recorded outcome, separate from operational "
        "statuses. not_started includes both explicit records and analyzer-filled missing slots; "
        "the final two columns distinguish these. Status labels describe the captured records.</p>",
    ]
    statuses = ("success", "no_candidates", "error", "timeout", "unsupported", "interrupted", "not_started")
    status_rows = []
    for method, data in report["methods"].items():
        method_records = [r for r in collection["attempts"] if r["method_id"] == method]
        status_rows.append(
            [
                method,
                *[data["status_counts"][status] for status in statuses],
                sum(r["status"] == "not_started" for r in method_records),
                data["planned_attempts"] - len(method_records),
            ]
        )
    parts.append(
        _table(
            ["Method", *statuses, "Explicit not_started", "Missing records"],
            status_rows,
            "Per-status counts (all planned slots)",
        )
    )
    parts.extend(
        [
            "<h2>Failure records</h2><p>Failure records may overlap within or across stages and may "
            "accompany successful attempts. These are not a sequential funnel. Affected-attempt counts "
            "deduplicate within each stage only; do not sum stages to count failed attempts. "
            "Missing slots contribute the analyzer's collection / missing_planned_attempt record.</p>",
        ]
    )
    for method, data in report["methods"].items():
        stages = _table(
            ["Stage", "Failure records", "Affected attempts"],
            [
                [stage, count, data["failure_stage_attempts"][stage]]
                for stage, count in data["failure_stage_records"].items()
            ],
            "Stage counts",
        )
        reasons = _table(
            ["Stage", "Reason", "Failure records (not attempts)"],
            [[r["stage"], r["reason"], r["count"]] for r in data["failure_reasons"]],
            "Reason counts" if data["failure_reasons"] else "No failure records",
        )
        parts.append(_details(f"{method}: failure detail", stages + reasons))
    parts.extend(
        [
            "<h2>Observed attempt costs</h2><p>Known values include unsuccessful attempts. "
            "Numeric field coverage is record availability, not complete resource measurement: "
            "observed / planned attempts for each field. Unknown GPU time or memory is "
            "not zero. Time sums and memory maxima summarize only observed values; incomplete coverage "
            "does not estimate missing costs. CPU and GPU time are distinct resources. Summed attempt "
            "wall time is not campaign duration; memory maxima are not sums or a combined process peak.</p>",
            "<p>Separate administrative setup, setup allocation shares, and campaign wall-through-checkpoint "
            "are not included here: those auxiliary fields are outside this analyzer's validation. "
            "No setup shares are added to attempt costs. These tables are not total campaign resource accounting.</p>",
            "<p><strong>Current runner reference, not an inferred policy for this snapshot:</strong> "
            "qualified_runner records sampled child-process-tree CPU lower bounds and maximum sampled "
            "child-tree RSS across phases, excluding parent/collector CPU and RSS from those fields. "
            "Historical or arbitrary collections may use different scopes; consult their captured policy below. "
            "The combined parent+child memory peak remains unknown in this view; separate maxima are never summed.</p>",
        ]
    )
    policy = collection.get("measurement_policy")
    policy = policy if isinstance(policy, dict) else {}
    policy_fields = dict.fromkeys((*analysis.COSTS, "parent_memory", *policy))
    parts.append(
        _table(
            ["Field", "Captured policy (recorded provenance; not analyzer-validated)"],
            [
                [
                    field,
                    policy[field]
                    if isinstance(policy.get(field), str) and policy[field].strip()
                    else "unavailable (missing or non-text policy)",
                ]
                for field in policy_fields
            ],
            "Captured measurement policy; no scope is inferred for unavailable fields",
        )
    )
    parts.append(
        _table(
            ["Method", "Field", "Aggregation", "Known value", "Observed / planned", "Numeric field coverage"],
            [
                [
                    method,
                    field,
                    data["costs"][field]["aggregation"],
                    _number(data["costs"][field]["observed_value"]),
                    f"{data['costs'][field]['observed_attempts']} / {data['planned_attempts']}",
                    "all slots numeric" if data["costs"][field]["complete"] else "partial or unavailable",
                ]
                for method, data in report["methods"].items()
                for field in analysis.COSTS
            ],
            "Collector-reported attempt measurements (seconds and MB as named)",
        )
    )
    for method, data in report["methods"].items():
        rows = [
            [kind, name, item["observed_total"], f"{item['observed_attempts']} / {item['planned_attempts']}"]
            for kind in ("native_units", "common_units")
            for name, item in data["costs"][kind].items()
        ]
        parts.append(
            _details(
                f"{method}: typed work units",
                _table(
                    ["Scope", "Unit", "Known total", "Observed / planned"],
                    rows,
                    "Observed work units; absent fields are unknown, not zero. Units are not interchangeable.",
                ),
            )
        )
    parts.extend(
        [
            "<h2>Conditional selected similarity</h2><p>Selected qualified outputs only, pooled across "
            "recorded attempts; not the planned-attempt distribution or the equal-target endpoint aggregation. "
            "Includes duplicate-family selected outputs; no family deduplication for these distributions. "
            "Counts are outputs, not attempts. Empty selections have unavailable statistics. "
            "Margins are collector-reported threshold margins.</p>",
            _table(
                ["Method", "Measure", "Metric", "Output count", "Minimum", "Median", "Maximum"],
                [
                    [
                        method,
                        kind,
                        metric,
                        values["count"],
                        *[
                            _number(values[key]) if values[key] is not None else "unavailable"
                            for key in ("minimum", "median", "maximum")
                        ],
                    ]
                    for method, data in report["methods"].items()
                    for kind in ("metrics", "margins")
                    for metric, values in data["selected_condition"][kind].items()
                ],
                "Conditional output distributions",
            ),
            "<h2>Attempt outcomes</h2><p>Every planned slot is shown in plan order. "
            "Zero selected outputs can be located by method, target, query, and seed, with recorded "
            "failure reasons where available. Missing records describe this snapshot only.</p>",
        ]
    )
    attempt_rows = []
    for row in report["attempts"]:
        original = recorded.get(tuple(row[key] for key in analysis.SLOT_KEYS))
        failure_text = (
            "; ".join(f"{f['stage']}: {f['reason']}" for f in original["failures"])
            if original is not None
            else "collection: missing_planned_attempt"
        )
        attempt_rows.append(
            [
                *[row[key] for key in analysis.SLOT_KEYS],
                "recorded" if original is not None else "missing",
                row["status"],
                len(original["qualified"]) if original is not None else 0,
                row["unique_families"],
                row["duplicate_family_outputs"],
                failure_text or "none recorded",
            ]
        )
    parts.append(
        _table(
            [
                "Method",
                "Target",
                "Query",
                "Seed",
                "Record",
                "Status",
                "Selected outputs",
                "Unique families",
                "Duplicate-family outputs",
                "Failure records",
            ],
            attempt_rows,
            "Planned attempt detail",
        )
    )
    provenance = {
        "plan_file_sha256": hashlib.sha256(plan_bytes).hexdigest(),
        "collection_file_sha256": actual_digest,
        **{
            key: report[key]
            for key in (
                "schema",
                "plan_sha256",
                "collection_sha256",
                "protocol_sha256",
                "collector_source_sha256",
                "scope",
                "inference_scope",
                "aggregation",
                "secondary_identity",
            )
        },
    }
    parts.append(
        _details(
            "Input digests and scope",
            "<p>File digests identify the exact supplied bytes; plan_sha256 and collection_sha256 "
            "are the analyzer's canonical JSON digests. The external collection digest is a caller "
            "trust requirement. This view recomputes from inputs rather than accepting a saved analysis. "
            "No auxiliary files, models, chemistry tools, or network resources were read. "
            "Intervals and paired differences are intentionally not presented in this descriptive view.</p>"
            f"<pre>{_text(json.dumps(provenance, indent=2, ensure_ascii=False))}</pre>",
        )
    )
    parts.append("</main></body></html>\n")
    return "".join(parts)


def write_report(plan: Path, collection: Path, *, collection_sha256: str, output: Path) -> None:
    """Render captured files to a new UTF-8 file; existing output is never overwritten."""
    rendered = render_report(
        Path(plan).read_bytes(),
        Path(collection).read_bytes(),
        collection_sha256=collection_sha256,
    )
    with Path(output).open("x", encoding="utf-8") as handle:
        handle.write(rendered)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--collection", type=Path, required=True)
    parser.add_argument("--collection-sha256", required=True, help="Trusted external collection byte SHA-256")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    write_report(args.plan, args.collection, collection_sha256=args.collection_sha256, output=args.output)


if __name__ == "__main__":
    main()
