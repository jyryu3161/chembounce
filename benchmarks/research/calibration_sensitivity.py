"""Reporting-only perturbations of reviewed, locked development calibration fits.

No fitting, scoring, model loading, or eligibility rescue occurs here. The CLI
requires externally supplied fit/procedure/addendum digests and verifies the
original immutable scoring source and artifact chain before reading features.
"""

import argparse
import json
from pathlib import Path
import platform
import statistics

from benchmarks.research import qualified_calibration as calibration


ADDENDUM_SHA256 = "9ac8bbc99cb6b9a2d2672e265de50a0a9303d30a509ca2a1295015dc4dcbb5a1"
SCHEMA = "qualified-development-reporting-sensitivity-v1"
FACTORS = (0.9, 0.95, 1.05, 1.1)
SCORING_FILES = {
    "panel.json",
    "procedure.json",
    "identity_manifest.json",
    "pairs.json",
    "protocol.json",
    "scores.json",
    "progress.json",
    "parity.json",
}


def validate_rows(rows):
    """Check scalar consistency without rerunning the hash-bound pose validator."""
    if not rows or len({r["pair_id"] for r in rows}) != len(rows):
        raise ValueError("Rows must have unique nonempty pair identities")
    for row in rows:
        if row["label"] not in {"active", "inactive"} or type(row["topology_change"]) is not bool:
            raise ValueError("Invalid fixed label or topology subgroup")
        trials = row["trials"]
        if [t.get("seed") for t in trials] != list(calibration.SEEDS) or any(
            type(t.get("seed")) is not int for t in trials
        ):
            raise ValueError("Every original seed outcome must remain present in order")
        if any(t.get("status") not in {"success", "failed", "unreported"} for t in trials):
            raise ValueError("Invalid seed outcome")
        complete = all(t["status"] == "success" for t in trials)
        if row.get("all_3d_seeds_succeeded") is not complete:
            raise ValueError("Original three-seed completeness mismatch")
        if row.get("status") not in ({"success"} if complete else {"failed", "error"}):
            raise ValueError("Original row status mismatch")
        if complete:
            values = [row["hit_tanimoto"]] + [
                t["info"]["hit"][metric] for t in trials for metric in ("shape", "pharmacophore")
            ]
            if any(not calibration.finite(v) or not 0 <= v <= 1 for v in values):
                raise ValueError("Invalid seed-level feature")
            expected = calibration.features(row)
        else:
            expected = None
        if row.get("features") != expected:
            raise ValueError("Original minimum features or failed-row eligibility changed")


def checked_cutoff(fit):
    cutoff = fit["cutoff"]
    if cutoff is None:
        if fit.get("fit_status") != "infeasible":
            raise ValueError("No-fit status mismatch")
    elif (
        fit.get("fit_status") != "fitted"
        or not isinstance(cutoff, list)
        or len(cutoff) != len(calibration.METRICS)
        or any(not calibration.finite(c) or not 0 < c <= 1 for c in cutoff)
    ):
        raise ValueError("Invalid locked cutoff")
    return cutoff


def predictions(rows, cutoffs, aggregation="minimum", seeds=calibration.SEEDS):
    result = {}
    for row in rows:
        cutoff = cutoffs[row["target_id"]]
        eligible = row["features"] is not None
        values = row["features"]
        if eligible and aggregation != "original_minimum":
            aggregate = statistics.median if aggregation == "median" else min
            values = [row["hit_tanimoto"]] + [
                aggregate(t["info"]["hit"][metric] for t in row["trials"] if t["seed"] in seeds)
                for metric in ("shape", "pharmacophore")
            ]
        result[row["pair_id"]] = bool(
            eligible and cutoff is not None and all(v >= c for v, c in zip(values, cutoff))
        )
    return result


def checked_predictions(expected, actual):
    if expected != actual or any(type(v) is not bool for v in actual.values()):
        raise ValueError("Locked predictions disagree with original minimum rule")


def variant_specs():
    return [
        {"id": "three_seed_median", "aggregation": "median", "seeds": list(calibration.SEEDS)},
        *[
            {
                "id": f"omit_seed{seed}",
                "aggregation": "minimum",
                "seeds": [s for s in calibration.SEEDS if s != seed],
            }
            for seed in calibration.SEEDS
        ],
        *[
            {
                "id": f"{metric}_times_{factor}",
                "aggregation": "original_minimum",
                "seeds": list(calibration.SEEDS),
                "feature": metric,
                "factor": factor,
            }
            for metric in calibration.METRICS
            for factor in FACTORS
        ],
    ]


def perturb_cutoffs(cutoffs, variant):
    effective, provenance = {}, {}
    for target, original in sorted(cutoffs.items()):
        adjusted = list(original) if original is not None else None
        item = {
            "original": original,
            "requested": None,
            "effective": adjusted,
            "clipped": False,
            "unchanged": True,
            "no_fit": original is None,
        }
        if original is not None and "feature" in variant:
            index = calibration.METRICS.index(variant["feature"])
            requested = original[index] * variant["factor"]
            adjusted[index] = min(1.0, requested)
            item.update(
                requested=requested,
                requested_hex=requested.hex(),
                effective_coordinate=adjusted[index],
                effective_coordinate_hex=float(adjusted[index]).hex(),
                clipped=requested > 1.0,
                unchanged=adjusted == original,
            )
        effective[target], provenance[target] = adjusted, item
    return effective, provenance


def summarize(rows, predicted, baseline, cutoffs):
    accepted = sorted(k for k, value in predicted.items() if value)
    lost = sorted(k for k, value in baseline.items() if value and not predicted[k])
    gained = sorted(k for k, value in predicted.items() if value and not baseline[k])
    return {
        "predictions": predicted,
        "accepted_pair_ids": accepted,
        "accepted_count": len(accepted),
        "accepted_to_rejected_pair_ids": lost,
        "accepted_to_rejected_count": len(lost),
        "rejected_to_accepted_pair_ids": gained,
        "rejected_to_accepted_count": len(gained),
        "report": calibration.prediction_report(
            rows, predicted, fitted_targets=[t for t, c in cutoffs.items() if c is not None]
        ),
    }


def cutoff_signature(cutoffs):
    """Numeric-equivalent int/float cutoffs are the same effective rule."""
    return calibration.json_hash(
        {t: [float(v) for v in c] if c is not None else None for t, c in cutoffs.items()}
    )


def application(rows, cutoffs, baseline):
    checked_predictions(predictions(rows, cutoffs), baseline)
    result = {
        "baseline": {"cutoffs_by_target": cutoffs, **summarize(rows, baseline, baseline, cutoffs)},
        "variants": [],
    }
    seen = {cutoff_signature(cutoffs): "baseline"}
    seen_target = {t: {cutoff_signature({t: c}): "baseline"} for t, c in cutoffs.items()}
    for variant in variant_specs():
        effective, provenance = perturb_cutoffs(cutoffs, variant)
        duplicate = None
        if "feature" in variant:
            digest = cutoff_signature(effective)
            duplicate = seen.get(digest)
            seen.setdefault(digest, variant["id"])
            for target, c in effective.items():
                key = cutoff_signature({target: c})
                provenance[target]["duplicate_of"] = seen_target[target].get(key)
                seen_target[target].setdefault(key, variant["id"])
        predicted = predictions(rows, effective, variant["aggregation"], variant["seeds"])
        result["variants"].append(
            {
                **variant,
                "cutoffs_by_target": effective,
                "cutoff_provenance": provenance,
                "duplicate_of": duplicate,
                **summarize(rows, predicted, baseline, effective),
            }
        )
    return result


def build_report(rows, fit):
    """Pure reporting operation; inputs must be verified by the CLI for real use."""
    validate_rows(rows)
    targets = sorted({r["target_id"] for r in rows})
    if fit.get("scope") != calibration.SPEC["scope"]:
        raise ValueError("Unexpected fit scope")
    final_cutoff = checked_cutoff(fit["final"])
    expected_minima = dict(zip(calibration.METRICS, final_cutoff)) if final_cutoff is not None else None
    if fit.get("minimum_scores") != expected_minima:
        raise ValueError("Locked minimum_scores disagree with final cutoff")
    folds = fit["folds"]
    if [f["test_target"] for f in folds] != targets:
        raise ValueError("Missing, duplicate or reordered locked target folds")
    oof_cutoffs, oof_predictions, reports = {}, {}, []
    for fold in folds:
        target = fold["test_target"]
        if fold["training_targets"] != sorted(set(targets) - {target}):
            raise ValueError("Locked fold training/test target mismatch")
        cutoff = checked_cutoff(fold)
        test = [r for r in rows if r["target_id"] == target]
        fold_report = application(test, {target: cutoff}, fold["test_predictions"])
        if fold.get("test_report") != fold_report["baseline"]["report"]:
            raise ValueError("Locked fold rate report mismatch")
        reports.append({"test_target": target, **fold_report})
        oof_cutoffs[target] = cutoff
        oof_predictions.update(fold["test_predictions"])
    checked_predictions(oof_predictions, fit["oof"]["predictions"])
    final = application(rows, dict.fromkeys(targets, final_cutoff), fit["final"]["predictions"])
    oof = application(rows, oof_cutoffs, oof_predictions)
    if any(
        fit[name]["report"] != result["baseline"]["report"]
        for name, result in (("final", final), ("oof", oof))
    ):
        raise ValueError("Locked aggregate rate report mismatch")
    return {
        "schema": SCHEMA,
        "scope": "previously_inspected_development_reporting_only",
        "refitting_performed": False,
        "production_rule_changed": False,
        "eligibility": "Original complete validated three-seed eligibility, unchanged in every variant",
        "counts": calibration.counts(rows),
        "computational_failure_pair_ids": sorted(r["pair_id"] for r in rows if r["features"] is None),
        "final": final,
        "locked_oof": oof,
        "folds": reports,
        "limitations": [
            "No selection or tuning of cutoffs, no confidence interval from variants",
            "Does not validate median or two-seed acceptance; failed original trials remain ineligible",
            "No-fit folds remain abstentions; operational inactive rejection is not geometric discrimination",
            "Previously inspected development data, not independent unseen validation",
        ],
    }


def checked_manifest(directory, required, *, optional=()):
    directory = Path(directory)
    manifest = json.loads((directory / "artifacts.json").read_text())
    if not set(required) <= set(manifest) <= set(required) | set(optional) or any(
        calibration.file_hash(directory / name) != digest for name, digest in manifest.items()
    ):
        raise ValueError("Artifact manifest integrity mismatch")
    return manifest


def load_inputs(
    *,
    scores_dir,
    fit_dir,
    fit_sha256,
    procedure_path,
    procedure_sha256,
    addendum_path,
    addendum_sha256,
    scoring_source_root,
):
    if addendum_sha256 != ADDENDUM_SHA256 or calibration.file_hash(addendum_path) != addendum_sha256:
        raise ValueError("Reporting addendum differs from reviewed frozen digest")
    addendum = json.loads(Path(addendum_path).read_text())
    if addendum["schema"] != SCHEMA or addendum["primary_procedure"]["sha256"] != procedure_sha256:
        raise ValueError("Reporting addendum and primary procedure disagree")
    procedure = calibration.load_procedure(procedure_path, procedure_sha256)
    source, fitted = Path(scores_dir), Path(fit_dir)
    fit_manifest = checked_manifest(fitted, {"fit.json"}, optional={"preset-candidate.json"})
    if fit_manifest["fit.json"] != fit_sha256:
        raise ValueError("Fit differs from externally approved fit digest")
    fit = json.loads((fitted / "fit.json").read_text())
    artifacts = checked_manifest(source, SCORING_FILES)
    provenance = fit["provenance"]
    if provenance["procedure_sha256"] != procedure_sha256 or provenance["scoring_artifacts"] != artifacts:
        raise ValueError("Fit and scoring provenance disagree")
    protocol = json.loads((source / "protocol.json").read_text())
    if (
        artifacts["procedure.json"] != procedure_sha256
        or artifacts["panel.json"] != procedure["panel_sha256"]
        or protocol["calibration_spec"] != calibration.SPEC
        or protocol["procedure_sha256"] != procedure_sha256
        or protocol["pairs_sha256"] != artifacts["pairs.json"]
        or protocol["source_sha256"] != provenance["source_sha256"]
    ):
        raise ValueError("Scoring protocol binding mismatch")
    code = {name: calibration.file_hash(Path(scoring_source_root) / name) for name in calibration.CODE_FILES}
    if code != provenance["source_sha256"] or code[
        "benchmarks/research/qualified_calibration.py"
    ] != calibration.file_hash(calibration.__file__):
        raise ValueError("Immutable scoring source or imported reporting math differs")
    for key in ("python", "packages", "pharmacophore_feature_definition"):
        if protocol["runtime"][key] != provenance["runtime"][key]:
            raise ValueError("Scoring and fit scientific runtime provenance disagree")
    progress = json.loads((source / "progress.json").read_text())
    pairs = json.loads((source / "pairs.json").read_text())
    rows = json.loads((source / "scores.json").read_text())
    if (
        progress.get("complete") is not True
        or progress.get("expected_pairs") != len(pairs)
        or progress.get("reported_pairs") != len(pairs)
        or progress.get("expected_trials") != len(pairs) * len(calibration.SEEDS)
    ):
        raise ValueError("Scoring did not complete the full pair/trial ledger")
    identity_path = calibration.verified_reference(procedure_path, procedure, "identity_manifest")
    if artifacts["identity_manifest.json"] != calibration.file_hash(identity_path):
        raise ValueError("Frozen identity inventory digest mismatch")
    calibration.verify_inventory(pairs, json.loads((source / "identity_manifest.json").read_text()))
    if (
        len(rows) != len(pairs)
        or calibration.counts(rows) != procedure["expected_counts"]
        or provenance["counts"] != procedure["expected_counts"]
        or any(any(row.get(k) != v for k, v in pair.items()) for row, pair in zip(rows, pairs))
    ):
        raise ValueError("Original scoring identities or denominators changed")
    return (
        rows,
        fit,
        {
            "addendum_sha256": addendum_sha256,
            "procedure_sha256": procedure_sha256,
            "fit_sha256": fit_sha256,
            "fit_artifacts": fit_manifest,
            "scoring_artifacts": artifacts,
            "scoring_source_sha256": code,
            "scoring_runtime": protocol["runtime"],
            "reporter_sha256": calibration.file_hash(__file__),
            "reporter_python": platform.python_version(),
            "validation_scope": "Original poses/settings were validated by the hash-bound immutable fit; reporter verifies artifacts, scalar summaries and locked predictions without rescoring",
        },
    )


def report(**kwargs):
    output = Path(kwargs.pop("output"))
    rows, fit, provenance = load_inputs(**kwargs)
    result = build_report(rows, fit)
    result["provenance"] = provenance
    output.mkdir(parents=True, exist_ok=False)
    calibration.save(output / "sensitivity.json", result)
    calibration.save(
        output / "artifacts.json", {"sensitivity.json": calibration.file_hash(output / "sensitivity.json")}
    )
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "scores-dir",
        "fit-dir",
        "fit-sha256",
        "procedure-path",
        "procedure-sha256",
        "addendum-path",
        "addendum-sha256",
        "scoring-source-root",
        "output",
    ):
        parser.add_argument("--" + name, required=True)
    report(**vars(parser.parse_args(argv)))


if __name__ == "__main__":
    main()
