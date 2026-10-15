"""Post hoc assay-context audit; never modifies frozen benchmark panels or scores."""

import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from benchmarks.research.assay_metadata import load_assay_metadata

BASE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "independent_controls", Path(__file__).resolve().with_name("analyze_controls.py")
)
METRICS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(METRICS)
UNKNOWN = re.compile(r"unknown origin|unspecified origin", re.I)
VARIANT = re.compile(r"\b(?:mutant|mutated|mutation|L858R|T790M|C797S|exon[- ]?19|L861Q|G719[A-Z])\b", re.I)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cluster(vals):
    if len(vals) < 2:
        return {"mean": vals[0] if vals else None, "ci95": None, "targets": len(vals)}
    return METRICS.cluster(vals)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--controls", type=Path, default=BASE / "control_scores.json")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=BASE / "recomputed-assay-heterogeneity.json")
    args = parser.parse_args()
    panel, controls = [json.loads(p.read_text()) for p in [args.panel, args.controls]]
    assays, provenance = load_assay_metadata(args.manifest)
    compounds = {c["id"]: c for c in panel["compounds"]}
    annotations = {}
    for ident, c in compounds.items():
        ids = c["supporting_assay_ids"]
        assert ids and all(i in assays for i in ids)
        descriptions = [assays[i]["record"].get("description") or "" for i in ids]
        annotations[ident] = {
            "target_id": c["target_id"],
            "supporting_assay_ids": ids,
            "all_supporting_bao_formats": sorted(
                {assays[i]["record"].get("bao_label") or "missing" for i in ids}
            ),
            "representative_bao_format": c["assay_assignment"]["record"].get("bao_label") or "missing",
            "representative_unknown_origin_text": bool(
                UNKNOWN.search(c["assay_assignment"]["record"].get("description") or "")
            ),
            "any_support_unknown_origin_text": any(UNKNOWN.search(d) is not None for d in descriptions),
            "any_support_variant_text": any(VARIANT.search(d) is not None for d in descriptions),
            "all_support_single_protein_format": all(
                assays[i]["record"].get("bao_label") == "single protein format" for i in ids
            ),
        }
    counts = {}
    for target in sorted({c["target_id"] for c in compounds.values()}):
        rows = [a for a in annotations.values() if a["target_id"] == target]
        counts[target] = {
            "compound_records": len(rows),
            "representative_bao_formats": dict(Counter(a["representative_bao_format"] for a in rows)),
            "representative_unknown_origin_text": sum(a["representative_unknown_origin_text"] for a in rows),
            "any_support_unknown_origin_text": sum(a["any_support_unknown_origin_text"] for a in rows),
            "any_support_variant_text": sum(a["any_support_variant_text"] for a in rows),
            "multiple_support_formats": sum(len(a["all_supporting_bao_formats"]) > 1 for a in rows),
        }
    sensitivity = {}
    for mode in [
        "shared_assay_id",
        "matched_single_format_no_unknown_origin_text",
        "single_protein_no_unknown_origin_text",
    ]:
        query_reports = []
        for query in panel["queries"]:
            qa = annotations[query["compound_id"]]
            rows = [c for c in controls if c["query_id"] == query["id"]]
            retained = []
            for row in rows:
                ca = annotations[row["compound_id"]]
                if mode == "shared_assay_id":
                    keep = bool(set(qa["supporting_assay_ids"]) & set(ca["supporting_assay_ids"]))
                elif mode == "matched_single_format_no_unknown_origin_text":
                    keep = (
                        len(qa["all_supporting_bao_formats"]) == 1
                        and qa["all_supporting_bao_formats"] == ca["all_supporting_bao_formats"]
                        and not qa["any_support_unknown_origin_text"]
                        and not ca["any_support_unknown_origin_text"]
                    )
                else:
                    keep = (
                        qa["all_support_single_protein_format"]
                        and ca["all_support_single_protein_format"]
                        and not qa["any_support_unknown_origin_text"]
                        and not ca["any_support_unknown_origin_text"]
                    )
                if keep:
                    retained.append(row)
            paired = [r for r in retained if all(r["scores"].get(m) is not None for m in METRICS.METRICS)]
            query_reports.append(
                {
                    "query_id": query["id"],
                    "target_id": qa["target_id"],
                    "controls_before_filter": len(rows),
                    "controls_retained": len(retained),
                    "active_controls_retained": sum(r["label"] == "active" for r in retained),
                    "inactive_controls_retained": sum(r["label"] == "inactive" for r in retained),
                    "paired_complete_cases": len(paired),
                    "retained_compound_ids": [r["compound_id"] for r in retained],
                    "paired_complete_compound_ids": [r["compound_id"] for r in paired],
                    "metrics": {m: METRICS.estimates(paired, m) for m in METRICS.METRICS},
                }
            )
        macro, per_target = {}, {}
        for metric in METRICS.METRICS:
            by_target = defaultdict(list)
            for q in query_reports:
                value = q["metrics"][metric]["roc_auc"]
                if value is not None:
                    by_target[q["target_id"]].append(value)
            per_target[metric] = {t: float(np.mean(v)) for t, v in by_target.items()}
            macro[metric] = cluster(list(per_target[metric].values()))
        sensitivity[mode] = {
            "query_results": query_reports,
            "target_paired_roc": per_target,
            "target_macro_paired_roc": macro,
            "retained_pairs": sum(q["controls_retained"] for q in query_reports),
            "paired_complete_pairs": sum(q["paired_complete_cases"] for q in query_reports),
        }
    report = {
        "design": "Post hoc descriptive sensitivity audit specified after the primary labelled-control results were seen; not pre-registered or used to select thresholds/methods. Frozen generation and source panels remain unchanged.",
        "source_sha256": {
            str(p): digest(p) for p in [args.panel, args.controls, args.manifest, Path(__file__)]
        },
        "compound_records": len(compounds),
        "per_target_context": counts,
        "query_annotations": [
            {"query_id": q["id"], "compound_id": q["compound_id"], **annotations[q["compound_id"]]}
            for q in panel["queries"]
        ],
        "sensitivity": sensitivity,
        "interpretation": [
            "Confidence 9 measures assay-to-target assignment, not assay-condition equivalence.",
            "A null variant_sequence is absence of a variant annotation, not experimental proof of wild-type sequence.",
            "Unknown-origin text is an unresolved conflict/ambiguity with structured human metadata, not proof of a non-human assay.",
            "BAO format may itself be incomplete or inconsistent with assay prose. Format matching does not equal matched experimental conditions.",
            "Shared assay-ID filters can leave too few active/inactive controls. Unestimable query/target groups remain explicit.",
            "All-support counts refer to the retained supporting measurements, not every publication on each molecule.",
        ],
        "documentation": "https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions",
        "assay_source_count": len(provenance["assay_sources"]),
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                "contexts": counts,
                "sensitivity": {
                    k: {
                        x: v[x]
                        for x in ["retained_pairs", "paired_complete_pairs", "target_macro_paired_roc"]
                    }
                    for k, v in sensitivity.items()
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
