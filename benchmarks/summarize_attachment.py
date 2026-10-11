#!/usr/bin/env python3
"""Summarize all attachment outcomes; identity controls are reported separately."""

import argparse
from collections import Counter
import json
from pathlib import Path
import statistics


def summarize(report):
    groups = {}
    for result in report["results"]:
        case = report["cases"][result["case"]]
        kind = "identity_control" if case["choice"] == "identity_control" else "replacement"
        key = (result["version"], result["beam_width"], kind)
        groups.setdefault(key, []).append(result)
    rows = []
    for (version, beam, kind), records in sorted(groups.items()):
        counts = Counter(record["status"] for record in records)
        candidates = [candidate for record in records for candidate in record.get("candidates", [])]
        durations = [r["seconds"] for r in records if "seconds" in r]
        rows.append(
            {
                "version": version,
                "beam_width": beam,
                "kind": kind,
                "cases": len(records),
                "statuses": dict(counts),
                "cases_with_failed_checks": sum(bool(r["checks_failed"]) for r in records),
                "returned_candidates": len(candidates),
                "candidates_meeting_tanimoto_0_5": sum(score >= 0.5 for score, _ in candidates),
                "changed_murcko": sum(r["changed_murcko"] for r in records),
                "median_tanimoto": statistics.median(score for score, _ in candidates)
                if candidates
                else None,
                "median_completed_call_seconds": statistics.median(durations) if durations else None,
            }
        )
    return {
        "scope": "Fixed selected panel; counts are descriptive, not population success estimates. "
        "No search, property filters, 3D, activity or synthesis validation. Timeouts remain in denominators; "
        "completed-call timing excludes imports and cannot substitute for censored runtime analysis.",
        "groups": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    summary = summarize(json.loads(args.report.read_text()))
    text = json.dumps(summary, indent=2) + "\n"
    if args.output:
        if args.output.exists():
            parser.error("Output exists; choose a new file")
        args.output.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
