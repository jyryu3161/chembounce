"""Common similarity assessment on independent labelled controls and generated rows.

The same fixed scoring protocol can assess every method's outputs. Control
selection is seed-fixed before scoring, failures remain listed, and an inactive
measurement is distinct from a property-matched decoy or an unlabelled molecule.
"""

import argparse
from collections import defaultdict
from dataclasses import replace
import json
from pathlib import Path
import random
import statistics
import time

from rdkit import Chem, DataStructs

from benchmarks.research.evaluation import (
    FP,
    cluster_interval,
    query_references,
    score_stability,
    binary_ranking_with_failures,
)
from benchmarks.research.protocol import digest, load_dataset, write_json
from chemistry import count_fingerprint


def similarity_scores(hit, product, seeds, conformers=None, reference_sdf=None):
    qmol, mol = Chem.MolFromSmiles(hit), Chem.MolFromSmiles(product)
    scores = {
        "bit_tanimoto": DataStructs.TanimotoSimilarity(FP.GetFingerprint(qmol), FP.GetFingerprint(mol)),
        "count_tanimoto": DataStructs.TanimotoSimilarity(count_fingerprint(qmol), count_fingerprint(mol)),
    }
    trials = []
    if conformers is not None:
        from chembounce_synthesis.conformers import ConformerSettings, score_candidate, reference_from_sdf

        settings = ConformerSettings(**conformers)
        for seed in seeds:
            started = time.monotonic()
            try:
                seeded = replace(settings, seed=seed)
                hit_pose = reference_from_sdf(reference_sdf, hit, seeded) if reference_sdf else None
                info, _, _ = score_candidate(product, hit, hit, seeded, hit_pose)
                metrics = {
                    "score": info["score"],
                    "hit_shape": info["hit"]["shape"],
                    "hit_pharmacophore": info["hit"]["pharmacophore"],
                }
                trials.append(
                    {
                        "seed": seed,
                        "status": "success",
                        "metrics": metrics,
                        "seconds": time.monotonic() - started,
                    }
                )
            except (ValueError, RuntimeError) as exc:
                trials.append(
                    {
                        "seed": seed,
                        "status": "failed",
                        "reason": str(exc),
                        "seconds": time.monotonic() - started,
                    }
                )
        if all(t["status"] == "success" for t in trials):
            for metric in ("hit_shape", "hit_pharmacophore", "score"):
                values = [t["metrics"][metric] for t in trials]
                scores[metric] = statistics.median(values)
                scores[metric + "_stability"] = score_stability(values)
    return scores, trials


def assess_controls(
    dataset_path,
    output,
    max_per_label=20,
    seeds=(42, 43, 44),
    conformers=None,
    count_tanimoto_threshold=0.5,
    run_results=None,
):
    if max_per_label <= 0 or not seeds or len(set(seeds)) != len(seeds):
        raise ValueError("Positive control cap and unique seeds required")
    if not 0 <= count_tanimoto_threshold <= 1:
        raise ValueError("Tanimoto threshold must be in [0,1]")
    dataset, audit = load_dataset(dataset_path)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    records, per_query = [], []
    for query in dataset["queries"]:
        hit, identity, references, _ = query_references(dataset, query)
        pools = defaultdict(list)
        for row, ident in references:
            if (
                ident.family_smiles != identity.family_smiles
                and ident.scaffold_family_smiles != identity.scaffold_family_smiles
            ):
                pools[row["label"]].append((row, ident))
        selected, population = [], []
        # Count-based prefilter recall uses every held-out active scaffold, not a
        # convenient scored subset. Three-dimensional controls are separately capped.
        for row, ident in pools["active"]:
            count = DataStructs.TanimotoSimilarity(
                count_fingerprint(Chem.MolFromSmiles(hit["smiles"])),
                count_fingerprint(Chem.MolFromSmiles(row["smiles"])),
            )
            population.append(
                {
                    "compound_id": row["id"],
                    "scaffold": ident.scaffold_family_smiles,
                    "count_tanimoto": count,
                    "passes_prefilter": count >= count_tanimoto_threshold,
                }
            )
        for label in sorted(pools):
            pool = sorted(pools[label], key=lambda item: item[0]["id"])
            random.Random(seeds[0]).shuffle(pool)
            selected.extend(pool[:max_per_label])
        local = []
        for row, ident in selected:
            scores, trials = similarity_scores(
                hit["smiles"], row["smiles"], seeds, conformers, query.get("resolved_reference_sdf")
            )
            result = {
                "query_id": query["id"],
                "target_id": hit["target_id"],
                "compound_id": row["id"],
                "label": row["label"],
                "source_id": row["source_id"],
                "scores": scores,
                "trials": trials,
                "all_3d_seeds_succeeded": all(t["status"] == "success" for t in trials) if trials else None,
            }
            records.append(result)
            local.append(result)
            write_json(output / "control_scores.json", records)
        metric_names = ["bit_tanimoto", "count_tanimoto"] + (
            ["hit_shape", "hit_pharmacophore", "score"] if conformers is not None else []
        )
        metrics = {
            metric: binary_ranking_with_failures([(r["label"], r["scores"].get(metric)) for r in local])
            for metric in metric_names
        }
        all_scaffolds = {r["scaffold"] for r in population}
        surviving_scaffolds = {r["scaffold"] for r in population if r["passes_prefilter"]}
        per_query.append(
            {
                "query_id": query["id"],
                "target_id": hit["target_id"],
                "selected_controls": len(local),
                "available_by_label": {label: len(pool) for label, pool in pools.items()},
                "failed_3d_controls": sum(r["all_3d_seeds_succeeded"] is False for r in local),
                "discrimination": metrics,
                "prefilter": {
                    "threshold": count_tanimoto_threshold,
                    "known_active_molecules": len(population),
                    "known_active_scaffolds": len(all_scaffolds),
                    "surviving_active_scaffolds": len(surviving_scaffolds),
                    "active_scaffold_recall": len(surviving_scaffolds) / len(all_scaffolds)
                    if all_scaffolds
                    else None,
                    "population": population,
                },
            }
        )
    target_means = defaultdict(list)
    for query in per_query:
        recall = query["prefilter"]["active_scaffold_recall"]
        if recall is not None:
            target_means[query["target_id"]].append(recall)
    summary = {
        "dataset_audit": audit,
        "protocol": {
            "max_per_label": max_per_label,
            "seeds": list(seeds),
            "conformers": conformers,
            "count_tanimoto_threshold": count_tanimoto_threshold,
        },
        "queries": per_query,
        "prefilter_target_macro_recall": cluster_interval(
            [statistics.mean(v) for v in target_means.values()]
        ),
        "interpretation": "Similarity discrimination is retrospective; it is not activity probability calibration. ROC/AP on successful scores may be biased when controls fail; failures and denominator counts are retained.",
    }
    if run_results is not None:
        runs = json.loads(Path(run_results).read_text())
        query_lookup = {q["id"]: q for q in dataset["queries"]}
        generated = []
        for run in runs:
            query = query_lookup[run["query_id"]]
            hit, _, _, _ = query_references(dataset, query)
            for row in run["evaluation"]["candidates"]:
                scores, trials = similarity_scores(
                    hit["smiles"], row["smiles"], seeds, conformers, query.get("resolved_reference_sdf")
                )
                generated.append(
                    {
                        "method_id": run["method_id"],
                        "query_id": run["query_id"],
                        "generation_seed": run["seed"],
                        "smiles": row["smiles"],
                        "scores": scores,
                        "trials": trials,
                    }
                )
                write_json(output / "generated_scores.json", generated)
        summary["run_results_sha256"] = digest(run_results)
        summary["generated_candidates_assessed"] = len(generated)
    write_json(output / "summary.json", summary)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-per-label", type=int, default=20)
    parser.add_argument("--seeds", nargs="+", type=int, default=[42, 43, 44])
    parser.add_argument(
        "--conformers", type=Path, help="Fixed ConformerSettings JSON; omit for 2D-only audit"
    )
    parser.add_argument("--count-tanimoto-threshold", type=float, default=0.5)
    parser.add_argument(
        "--run-results",
        type=Path,
        help="Optional runner runs.json; all methods get the same scoring protocol",
    )
    args = parser.parse_args(argv)
    assess_controls(
        args.dataset,
        args.output,
        args.max_per_label,
        args.seeds,
        json.loads(args.conformers.read_text()) if args.conformers else None,
        args.count_tanimoto_threshold,
        args.run_results,
    )


if __name__ == "__main__":
    main()
