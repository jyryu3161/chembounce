"""Method-independent outcome metrics; failures remain in every query denominator."""

from collections import Counter, defaultdict
import itertools
import math
import statistics

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from chembounce_synthesis.identity import chemical_identity
from chemistry import count_fingerprint


FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)


def finite_scores(values):
    if not isinstance(values, list) or not values:
        return None
    if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in values):
        return None
    return [float(v) for v in values]


def score_stability(samples):
    values = finite_scores(samples)
    if values is None:
        return None
    return {
        "n": len(values),
        "mean": statistics.mean(values),
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
        "range": max(values) - min(values),
        "sample_sd": statistics.stdev(values) if len(values) > 1 else None,
    }


def query_references(dataset, query):
    cache = getattr(dataset, "reference_cache", {})
    if query["id"] in cache:
        return cache[query["id"]]
    identities = getattr(dataset, "identities", {})
    lookup = getattr(dataset, "compounds_by_id", None) or {r["id"]: r for r in dataset["compounds"]}
    hit = lookup[query["compound_id"]]
    identity = identities.get(hit["id"]) or chemical_identity(hit["smiles"])
    records = [r for r in dataset["compounds"] if r["target_id"] == hit["target_id"] and r["split"] == "test"]
    references = [(r, identities.get(r["id"]) or chemical_identity(r["smiles"])) for r in records]
    active_scaffolds = {
        ident.scaffold_family_smiles
        for r, ident in references
        if r["label"] == "active"
        and ident.has_scaffold
        and ident.scaffold_family_smiles != identity.scaffold_family_smiles
        and ident.family_smiles != identity.family_smiles
    }
    result = hit, identity, references, active_scaffolds
    cache[query["id"]] = result
    return result


def evaluate_candidates(dataset, query, candidates, status="success", verify_evidence=None):
    """Scaffold recovery is explicitly a proxy, never a claim of candidate activity.

    Invalid outputs and duplicated chemical states are retained in denominator
    counts. An unsuccessful/over-budget run has zero accepted candidates even if
    it leaves a partial response artifact on disk.
    """
    hit, query_identity, references, active_scaffolds = query_references(dataset, query)
    if status not in {"success", "no_candidates"}:
        candidates = []
    by_label_exact = defaultdict(set)
    by_label_family = defaultdict(set)
    for record, ident in references:
        by_label_exact[record["label"]].add(ident.exact_smiles)
        by_label_family[record["label"]].add(ident.family_smiles)
    valid, failures = [], []
    for index, row in enumerate(candidates):
        try:
            if not isinstance(row, dict):
                raise ValueError("Candidate must be an object")
            ident = chemical_identity(row["smiles"])
            mol = Chem.MolFromSmiles(ident.exact_smiles)
            qmol = Chem.MolFromSmiles(query_identity.exact_smiles)
            same_family = ident.family_smiles == query_identity.family_smiles
            same_scaffold = ident.scaffold_family_smiles == query_identity.scaffold_family_smiles
            evidence = (
                verify_evidence(row)
                if verify_evidence
                else {"level": "unverified", "reason": "No verifier configured"}
            )
            item = {
                "rank": index + 1,
                "smiles": ident.exact_smiles,
                "family": ident.family_smiles,
                "scaffold": ident.scaffold_family_smiles,
                "ring_skeleton": ident.ring_skeleton_smiles,
                "is_hop": ident.has_scaffold and not same_family and not same_scaffold,
                "same_family_as_query": same_family,
                "hit_bit_tanimoto": DataStructs.TanimotoSimilarity(
                    FP.GetFingerprint(qmol), FP.GetFingerprint(mol)
                ),
                "hit_count_tanimoto": DataStructs.TanimotoSimilarity(
                    count_fingerprint(qmol), count_fingerprint(mol)
                ),
                "exact_active_match": ident.exact_smiles in by_label_exact["active"],
                "active_family_match": ident.family_smiles in by_label_family["active"],
                "exact_inactive_match": ident.exact_smiles in by_label_exact["inactive"],
                "inactive_family_match": ident.family_smiles in by_label_family["inactive"],
                "decoy_family_match": ident.family_smiles in by_label_family["decoy"],
                "active_scaffold_match": ident.scaffold_family_smiles in active_scaffolds,
                "score_stability": score_stability(row.get("score_samples")),
                "route_evidence": evidence,
                "declared_route_evidence": row.get("route_evidence"),
                "declared_external_synthesis_evidence": row.get("synthesis_evidence"),
            }
            valid.append(item)
        except (ValueError, KeyError, TypeError, RuntimeError) as exc:
            failures.append({"rank": index + 1, "error": str(exc)})
    unique = {}
    for row in valid:
        unique.setdefault(row["family"], row)
    hops = [r for r in unique.values() if r["is_hop"]]
    recovered = {r["scaffold"] for r in hops if r["active_scaffold_match"]}
    output_count = len(candidates)
    metrics = {
        "run_success": int(status == "success" and bool(valid)),
        "candidate_count": output_count,
        "valid_candidate_count": len(valid),
        "valid_fraction": len(valid) / output_count if output_count else 0.0,
        "unique_families": len(unique),
        "duplicate_chemical_states": len(valid) - len(unique),
        "unchanged_query_families": sum(r["same_family_as_query"] for r in unique.values()),
        "hop_families": len(hops),
        "unique_hop_scaffolds": len({r["scaffold"] for r in hops}),
        "unique_hop_ring_skeletons": len({r["ring_skeleton"] for r in hops}),
        "reference_active_scaffolds": len(active_scaffolds),
        "recovered_active_scaffolds": len(recovered),
        "active_scaffold_recall_proxy": len(recovered) / len(active_scaffolds) if active_scaffolds else None,
        "active_scaffold_precision_proxy": sum(r["active_scaffold_match"] for r in hops) / len(hops)
        if hops
        else 0.0,
        "exact_known_active_matches": sum(r["exact_active_match"] for r in unique.values()),
        "known_inactive_family_matches": sum(r["inactive_family_match"] for r in unique.values()),
        "decoy_family_matches": sum(r["decoy_family_match"] for r in unique.values()),
        "verified_rule_routes": sum(
            r["route_evidence"]["level"] == "catalog_and_rule_verified" for r in unique.values()
        ),
    }
    return {"metrics": metrics, "candidates": valid, "invalid_candidates": failures}


def binary_ranking_metrics(labelled_scores):
    """Discrimination on active/inactive controls, excluding unlabelled and decoys.

    Arbitrary similarity scores are not probabilities: no calibration or Brier
    claim is made without an independently fitted probability model.
    """
    rows = [
        (int(label == "active"), float(score))
        for label, score in labelled_scores
        if label in {"active", "inactive"} and math.isfinite(score)
    ]
    positives = [score for label, score in rows if label]
    negatives = [score for label, score in rows if not label]
    report = {
        "active_count": len(positives),
        "inactive_count": len(negatives),
        "roc_auc": None,
        "average_precision": None,
        "probability_calibration": "not_applicable_to_similarity_scores",
    }
    if not positives or not negatives:
        report["reason"] = "Both independently labelled active and inactive controls are required"
        return report
    report["roc_auc"] = sum((a > b) + 0.5 * (a == b) for a in positives for b in negatives) / (
        len(positives) * len(negatives)
    )
    # Threshold-wise AP handles ties without depending on input row order.
    count, tp, previous_recall, ap = 0, 0, 0.0, 0.0
    for score in sorted({s for _, s in rows}, reverse=True):
        at = [label for label, value in rows if value == score]
        count += len(at)
        tp += sum(at)
        recall = tp / len(positives)
        ap += (recall - previous_recall) * tp / count
        previous_recall = recall
    report["average_precision"] = ap
    return report


def binary_ranking_with_failures(labelled_scores):
    """Sharp ranking bounds for missing active/inactive scores, plus observed AUC.

    Each missing active can rank below every measured value while each missing
    inactive ranks above it (lower bound), or vice versa (upper bound). No failed
    control is dropped from the bound denominators.
    """
    selected = [(label, score) for label, score in labelled_scores if label in {"active", "inactive"}]
    measured = [
        (label, float(score))
        for label, score in selected
        if score is not None and isinstance(score, (int, float)) and math.isfinite(score)
    ]
    report = binary_ranking_metrics(measured)
    missing = {
        label: sum(
            s is None or not isinstance(s, (int, float)) or not math.isfinite(s)
            for lab, s in selected
            if lab == label
        )
        for label in ("active", "inactive")
    }
    report["missing_scores"] = missing
    report["total_active_controls"] = sum(label == "active" for label, _ in selected)
    report["total_inactive_controls"] = sum(label == "inactive" for label, _ in selected)
    report["success_conditional_metrics"] = True
    report["roc_auc_bounds"] = report["average_precision_bounds"] = None
    if not report["total_active_controls"] or not report["total_inactive_controls"]:
        return report
    finite = [score for _, score in measured]
    low, high = (min(finite) - 1, max(finite) + 1) if finite else (-1, 1)
    extremes = []
    for optimistic in (False, True):
        complete = measured + [
            (label, high if (label == "active") == optimistic else low)
            for label, n in missing.items()
            for _ in range(n)
        ]
        extremes.append(binary_ranking_metrics(complete))
    for metric in ("roc_auc", "average_precision"):
        report[metric + "_bounds"] = [value[metric] for value in extremes]
    return report


def cluster_interval(values, seed=20261015, resamples=2000):
    """Resample target means; candidates and repeat seeds are not independent n."""
    if not values:
        return {"mean": None, "ci95": None, "independent_targets": 0}
    array = np.asarray(values, dtype=float)
    result = {"mean": float(array.mean()), "ci95": None, "independent_targets": len(array)}
    if len(array) < 2:
        result["reason"] = "Fewer than two targets; between-target uncertainty is unidentifiable"
        return result
    rng = np.random.default_rng(seed)
    draws = rng.choice(array, size=(resamples, len(array)), replace=True).mean(axis=1)
    result["ci95"] = [float(x) for x in np.quantile(draws, [0.025, 0.975])]
    return result


def aggregate_runs(runs):
    """Equal weight per target after averaging queries and repeated seeds."""
    result = {}
    for method in sorted({r["method_id"] for r in runs}):
        selected = [r for r in runs if r["method_id"] == method]
        target_rows = defaultdict(list)
        for row in selected:
            target_rows[row["target_id"]].append(row)
        metric_names = sorted({k for r in selected for k in r["evaluation"]["metrics"]})
        per_target = {}
        for target, rows in sorted(target_rows.items()):
            per_target[target] = {
                name: statistics.mean(vals)
                if (
                    vals := [
                        r["evaluation"]["metrics"][name]
                        for r in rows
                        if r["evaluation"]["metrics"][name] is not None
                    ]
                )
                else None
                for name in metric_names
            }
        repeated = defaultdict(list)
        for row in selected:
            repeated[row["query_id"]].append(row)
        seed_stability = {}
        for query_id, rows in repeated.items():
            sets = [{c["family"] for c in r["evaluation"]["candidates"]} for r in rows]
            similarities = [
                len(a & b) / len(a | b) if a | b else 0.0 for a, b in itertools.combinations(sets, 2)
            ]
            recalls = [
                r["evaluation"]["metrics"]["active_scaffold_recall_proxy"]
                for r in rows
                if r["evaluation"]["metrics"]["active_scaffold_recall_proxy"] is not None
            ]
            seed_stability[query_id] = {
                "seeds_attempted": len(rows),
                "mean_family_jaccard": statistics.mean(similarities) if similarities else None,
                "active_scaffold_recall": score_stability(recalls),
                "empty_pairs_jaccard_policy": "zero; two failures are not stable successful recovery",
            }
        result[method] = {
            "attempts": len(selected),
            "status_counts": dict(Counter(r["status"] for r in selected)),
            "per_target": per_target,
            "target_macro_metrics": {
                name: cluster_interval([p[name] for p in per_target.values() if p[name] is not None])
                for name in metric_names
            },
            "seed_stability": seed_stability,
            "wall_seconds_total": sum(r["resources"]["wall_seconds"] for r in selected),
            "sampled_peak_rss_mb": max(
                (r["resources"].get("sampled_peak_rss_mb", 0) for r in selected), default=0
            ),
            "uncertainty_unit": "target (queries and seeds averaged within each target)",
        }
    return result
