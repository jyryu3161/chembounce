"""Executable baseline adapters. Read only public query inputs, never test labels.

A common wall budget is the primary computational comparison. Counted objective
calls have their declared units recorded: a fingerprint comparison, a 3D
optimization, and a model forward pass are not equivalent computational costs.
"""

import argparse
import itertools
import json
from pathlib import Path
import random
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rdkit import Chem, DataStructs

from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.common import FP, molecule
from chembounce_synthesis.identity import chemical_identity
from chembounce_synthesis.reactions import Rules
from benchmarks.research.protocol import write_json


def remaining_seconds(request):
    return request.get("deadline_monotonic", float("inf")) - time.monotonic()


def time_available(request, next_work_seconds=0.0):
    margin = request.get("finalization_margin_seconds", 2.0)
    return remaining_seconds(request) > margin + max(0.0, next_work_seconds)


def budget_stop_metadata(request, attempted, cap):
    return {
        "anytime_stop": request.get("_anytime_stop_reason")
        or (
            "wall_margin"
            if not time_available(request)
            else "evaluation_cap"
            if attempted >= cap
            else "search_exhausted"
        ),
        "remaining_seconds": max(0, remaining_seconds(request)) if "deadline_monotonic" in request else None,
        "finalization_margin_seconds": request.get("finalization_margin_seconds"),
    }


def nearest_catalog(request):
    """Scan the declared library under the outer wall budget, then select outputs."""
    catalog = Catalog(request["assets"]["catalog"])
    query = request["query"]["smiles"]
    qfp = FP.GetFingerprint(molecule(query))
    family = chemical_identity(query).family_smiles
    pool = sorted(catalog.mols)
    random.Random(request["seed"]).shuffle(pool)
    budget = request["budget"]
    rows, evaluated = [], 0
    for smiles in pool[: budget["max_evaluations"]]:
        if not time_available(request, max(0.5, 0.05 * budget["max_candidates"])):
            request["_anytime_stop_reason"] = "wall_margin"
            break
        score = DataStructs.TanimotoSimilarity(qfp, catalog.fps[smiles])
        evaluated += 1
        rows.append((score, smiles))
    rows.sort(key=lambda item: (-item[0], item[1]))
    # Expensive family normalization is needed only for potentially returned
    # neighbours; do not spend the whole budget normalizing an entire library.
    selected, families, rejected = [], {family}, []
    for score, smiles in rows:
        if len(selected) >= budget["max_candidates"] or remaining_seconds(request) <= min(
            0.2, request.get("finalization_margin_seconds", 2.0)
        ):
            break
        try:
            parent = chemical_identity(smiles).family_smiles
        except ValueError as exc:
            rejected.append({"smiles": smiles, "reason": str(exc)})
            continue
        if parent in families:
            continue
        families.add(parent)
        selected.append(
            {
                "smiles": smiles,
                "score_samples": [score],
                "route_evidence": {"level": "catalog_record", "catalog_sha256": catalog.sha256},
            }
        )
    return {
        "candidates": selected,
        "evaluations": evaluated,
        "evaluation_unit": "Morgan fingerprint comparisons",
        "metadata": {
            "library_size": len(catalog.mols),
            "search_subset_size": evaluated,
            "library_sha256": catalog.sha256,
            "subset_policy": "seeded shuffle without replacement; complete scan if budget permits",
            "family_rejections": rejected,
            **budget_stop_metadata(request, evaluated, budget["max_evaluations"]),
            "scope": "Declared-library retrieval; no synthetic route or biological activity inference",
        },
    }


def nearest_3d(request):
    """3D retrieval from a bounded seed-fixed pool of declared training ligands."""
    from benchmarks.research.controls import similarity_scores

    catalog = Catalog(request["assets"]["catalog"])
    pool = sorted(catalog.mols)
    random.Random(request["seed"]).shuffle(pool)
    family = chemical_identity(request["query"]["smiles"]).family_smiles
    rows, failures, attempted, longest_screen = [], [], 0, 1.0
    for smiles in pool:
        if attempted >= request["budget"]["max_evaluations"] or not time_available(
            request, 1.5 * longest_screen
        ):
            break
        attempted += 1
        try:
            candidate_family = chemical_identity(smiles).family_smiles
        except ValueError as exc:
            failures.append({"smiles": smiles, "stage": "chemical_identity", "reason": str(exc)})
            continue
        if candidate_family == family:
            continue
        started = time.monotonic()
        scores, trials = similarity_scores(
            request["query"]["smiles"],
            smiles,
            [request["seed"]],
            request["settings"].get("conformers", {"num_conformers": 6}),
            request["query"].get("reference_sdf"),
        )
        longest_screen = max(longest_screen, time.monotonic() - started)
        if "score" not in scores:
            failures.append({"smiles": smiles, "trials": trials})
            continue
        rows.append(
            {
                "smiles": smiles,
                "score_samples": [scores["score"]],
                "route_evidence": {"level": "catalog_record", "catalog_sha256": catalog.sha256},
            }
        )
    rows.sort(key=lambda r: (-r["score_samples"][0], r["smiles"]))
    return {
        "candidates": rows[: request["budget"]["max_candidates"]],
        "evaluations": attempted,
        "evaluation_unit": "3D candidate screens",
        "metadata": {
            "failures": failures,
            "library_size": len(catalog.mols),
            **budget_stop_metadata(request, attempted, request["budget"]["max_evaluations"]),
            "scope": "3D retrieval from declared development library; no new synthesis claim",
        },
    }


def reaction_enumeration(request):
    """One-step forward enumeration with a hard global reactant-tuple budget.

    Rules and reactant pools are seed-shuffled; enumeration is deliberately simple
    and fully reported. It is not an exhaustive synthesis planner or RuSH proxy.
    """
    started = time.monotonic()
    catalog = Catalog(request["assets"]["catalog"])
    rules = Rules(request["assets"]["rules"])
    budget = request["budget"]
    rng = random.Random(request["seed"])
    qfp = FP.GetFingerprint(molecule(request["query"]["smiles"]))
    family = chemical_identity(request["query"]["smiles"]).family_smiles
    rule_ids = sorted(rules.rules)
    rng.shuffle(rule_ids)
    startup_seconds = time.monotonic() - started
    proposals = [(rule_id, None) for rule_id in rule_ids]
    products, combinations, scored, failures = {}, 0, 0, []
    visited, initialized, pool_cache, preparation = [], [], {}, []
    combinations_by_rule = {}

    def compatible_pool(role):
        if role in pool_cache:
            return pool_cache[role]
        query = Chem.MolFromSmarts(role)
        if query is None:
            raise ValueError("Invalid reaction-role SMARTS")
        prep_start, values, inspected = time.monotonic(), [], 0
        for smiles, mol in catalog.mols.items():
            # Unlike a single full-library compatible() call, this can finish
            # with explicit incomplete preparation before the outer hard kill.
            if inspected % 256 == 0 and not time_available(request):
                break
            inspected += 1
            if mol.HasSubstructMatch(query, useChirality=True):
                values.append(smiles)
        complete = inspected == len(catalog.mols)
        preparation.append(
            {
                "role_smarts": role,
                "catalog_records_inspected": inspected,
                "compatible_records": len(values),
                "complete": complete,
                "seconds": time.monotonic() - prep_start,
            }
        )
        if not complete:
            request["_anytime_stop_reason"] = "wall_margin_during_rule_preparation"
            return None
        rng.shuffle(values)
        pool_cache[role] = values
        return values

    # Initialize a rule only on its first turn, immediately evaluate one tuple,
    # then visit the next rule. Completed iterators and role pools are reused.
    while proposals and combinations < budget["max_evaluations"] and time_available(request):
        next_round = []
        for rule_id, tuples in proposals:
            if combinations >= budget["max_evaluations"] or not time_available(request):
                break
            if tuples is None:
                visited.append(rule_id)
                pools = []
                for role in rules.get(rule_id).roles:
                    if not time_available(request):
                        break
                    values = compatible_pool(role)
                    if values is None:
                        break
                    pools.append(values)
                if len(pools) != len(rules.get(rule_id).roles):
                    break
                initialized.append(rule_id)
                if not all(pools):
                    continue
                tuples = iter(itertools.product(*pools))
            try:
                reactants = next(tuples)
            except StopIteration:
                continue
            next_round.append((rule_id, tuples))
            combinations += 1
            combinations_by_rule[rule_id] = combinations_by_rule.get(rule_id, 0) + 1
            generated, capped = rules.apply(
                rule_id, tuple(reactants), request["settings"].get("reaction_product_limit", 64)
            )
            for smiles in generated:
                if smiles in products:
                    continue
                # Count attempted validation/scoring, including standardization
                # failures, so one difficult product cannot erase other results.
                if scored >= budget["max_evaluations"]:
                    break
                scored += 1
                try:
                    candidate_family = chemical_identity(smiles).family_smiles
                except ValueError as exc:
                    failures.append({"smiles": smiles, "stage": "chemical_identity", "reason": str(exc)})
                    continue
                if candidate_family == family:
                    continue
                score = DataStructs.TanimotoSimilarity(qfp, FP.GetFingerprint(molecule(smiles)))
                nodes = [
                    {"id": f"b{i}", "building_block_id": sorted(catalog.by_smiles[s])[0]}
                    for i, s in enumerate(reactants)
                ]
                nodes.append(
                    {"id": "r0", "rule_id": rule_id, "reactants": [n["id"] for n in nodes], "product": smiles}
                )
                products[smiles] = {
                    "smiles": smiles,
                    "score_samples": [score],
                    "route": {"schema_version": 1, "root": "r0", "nodes": nodes},
                    "route_evidence": {"level": "rule_replay_claim", "enumeration_truncated": capped},
                }
        proposals = next_round
    rows = sorted(products.values(), key=lambda r: (-r["score_samples"][0], r["smiles"]))
    return {
        "candidates": rows[: budget["max_candidates"]],
        "evaluations": max(combinations, scored),
        "evaluation_unit": "maximum of reaction tuples and unique product fingerprint comparisons",
        "metadata": {
            "reaction_combinations": combinations,
            "reaction_combinations_by_rule": combinations_by_rule,
            "scored_products": scored,
            "startup_seconds": startup_seconds,
            "rule_preparation_seconds": sum(p["seconds"] for p in preparation),
            "visited_rules": visited,
            "initialized_rules": initialized,
            "unvisited_rules": [rule_id for rule_id in rule_ids if rule_id not in visited],
            "rule_pool_preparation": preparation,
            "rule_allocation": "seeded round-robin; lazy pool initialization followed immediately by one tuple",
            "failures": failures,
            **budget_stop_metadata(request, combinations, budget["max_evaluations"]),
            "rules_sha256": rules.sha256,
            "catalog_sha256": catalog.sha256,
            "scope": "One-step forward-rule baseline; no condition or activity validation",
        },
    }


def synllama_only(request):
    """Direct SynLlama projection of the query, without ChemBounce designs or 3D ranking.

    Exact rule/catalog replay remains required; raw model product claims are
    retained in proposal artifacts but never reported as verified candidates.
    """
    from chembounce_synthesis.routes import SearchLimits, resolve_proposal, verify_route
    from chembounce_synthesis.synllama import SynLlama

    settings = dict(request["settings"])
    budget = request["budget"]
    settings["model_path"] = request["assets"]["model_directory"]
    settings["prompt_path"] = request["assets"]["prompt"]
    search = settings.pop("search", {})
    settings["seed"] = request["seed"]
    provider = SynLlama(**settings)
    catalog, rules = Catalog(request["assets"]["catalog"]), Rules(request["assets"]["rules"])
    proposals = provider.propose(
        request["query"]["smiles"],
        deadline_monotonic=(
            request["deadline_monotonic"] - request.get("finalization_margin_seconds", 2.0)
            if "deadline_monotonic" in request
            else None
        ),
    )
    write_json(
        Path(request["output_dir"]) / "proposals.json",
        {"provider": provider.provenance, "proposals": proposals},
    )
    candidates, attempts, failures = {}, 0, []
    qfp = FP.GetFingerprint(molecule(request["query"]["smiles"]))
    query_family = chemical_identity(request["query"]["smiles"]).family_smiles
    for proposal in proposals:
        if attempts >= budget["max_evaluations"] or not time_available(request):
            break
        limits = SearchLimits(**(search | {"max_combinations": budget["max_evaluations"] - attempts}))
        try:
            routes, report = resolve_proposal(proposal["response"], catalog, rules, limits)
            attempts += report.combinations
            for resolved in routes:
                route = resolved["route"]
                evidence = verify_route(route, catalog, rules, limits.reaction_product_limit)
                smiles = evidence["product"]
                try:
                    candidate_family = chemical_identity(smiles).family_smiles
                except ValueError as exc:
                    failures.append(f"Chemical identity for {smiles}: {exc}")
                    continue
                if candidate_family == query_family:
                    continue
                score = DataStructs.TanimotoSimilarity(qfp, FP.GetFingerprint(molecule(smiles)))
                candidates.setdefault(
                    smiles,
                    {
                        "smiles": smiles,
                        "route": route,
                        "score_samples": [score],
                        "route_evidence": {"level": evidence["status"]},
                    },
                )
        except (ValueError, RuntimeError, KeyError) as exc:
            failures.append(str(exc))
    ranked = sorted(candidates.values(), key=lambda r: (-r["score_samples"][0], r["smiles"]))
    return {
        "candidates": ranked[: budget["max_candidates"]],
        "evaluations": attempts,
        "evaluation_unit": "reaction tuples in model-route reconstruction",
        "metadata": {
            "provider": provider.provenance,
            "failures": failures,
            **budget_stop_metadata(request, attempts, budget["max_evaluations"]),
            "ranking": "Morgan Tanimoto to public query; unchanged chemical family excluded",
            "scope": "Direct SynLlama projection; no design generator, 3D ranking, or expansion",
        },
    }


def pipeline(request):
    """Full ChemBounce generation from the public hit, followed by synthesis.

    No held-out design/reference is accepted. ``generator_options`` controls
    documented feature ablations; each method/configuration has its own hash.
    """
    from chembounce import chembounce
    from chembounce_synthesis.workflow import load_config, run

    settings = request["settings"]
    budget = request["budget"]
    options = dict(settings.get("generator_options", {}))
    if settings.get("disable_3d", False) and options.get("calculate_shape", True):
        raise ValueError(
            "disable_3d requires generator_options.calculate_shape=false in the paired parent and variant; otherwise generation still performs 3D work"
        )
    options.update(
        output_dir=Path(request["output_dir"]) / "designs",
        seed=request["seed"],
        overall_max_n=budget["max_candidates"],
    )
    if "fingerprint_db" in request["assets"]:
        options["fingerprint_db_file"] = request["assets"]["fingerprint_db"]
    if "scaffold_db" in request["assets"]:
        options["scaffold_db_file"] = request["assets"]["scaffold_db"]
    designs, cost = chembounce(request["query"]["smiles"], **options)
    if settings.get("disable_synthesis", False):
        # A compound baseline: both synthesis and all 3D evaluation are removed.
        rows = [
            {"smiles": s, "score_samples": [float(t)]}
            for s, t in zip(designs["Final structure"], designs["Tanimoto Similarity"])
        ]
        return {
            "candidates": rows[: budget["max_candidates"]],
            "evaluations": len(rows),
            "evaluation_unit": "returned design scores (internal generator work also in resource artifact)",
            "metadata": {
                "generator_cost": cost,
                "baseline": "design_only_2d",
                "single_factor_ablation": False,
                "removed_components": [
                    "synthesis_reconstruction",
                    "route_evidence",
                    "3d_scoring",
                    "refinement",
                ],
                "budget_accounting": "partial",
            },
        }
    config = load_config(request["assets"]["synthesis_config"])
    config.update(settings.get("synthesis_overrides", {}))
    config["catalog"], config["rules"] = request["assets"]["catalog"], request["assets"]["rules"]
    if request["query"].get("reference_sdf"):
        config["reference_sdf"] = request["query"]["reference_sdf"]
    config["max_scored_products"] = budget["max_evaluations"]
    config["max_results"] = budget["max_candidates"]
    config.setdefault("conformers", {})["seed"] = request["seed"]
    if config["provider"].get("kind") == "synllama":
        config["provider"]["seed"] = request["seed"]
    out = Path(request["output_dir"]) / "synthesis"
    rows = run(
        request["query"]["smiles"],
        list(designs["Final structure"]),
        config,
        out,
        scoring_mode="2d" if settings.get("disable_3d", False) else "3d",
        deadline_monotonic=(
            request["deadline_monotonic"] - request.get("finalization_margin_seconds", 2.0)
            if "deadline_monotonic" in request
            else None
        ),
    )
    summary = json.loads((out / "summary.json").read_text())
    route_data = json.loads((out / "routes.json").read_text())
    route_records = route_data.get("routes", route_data) if isinstance(route_data, dict) else route_data
    lookup = {r["route_id"]: r for r in route_records}
    refinement_path = out / "refinement.json"
    refinement = json.loads(refinement_path.read_text()) if refinement_path.is_file() else {}
    refinement_rows = refinement.get("candidates", []) if isinstance(refinement, dict) else []
    refinements = {
        (r.get("product"), r.get("design")): r for r in refinement_rows if r.get("status") == "success"
    }
    candidates = []
    for row in rows:
        route_record = lookup.get(row["route_id"], {})
        trials = refinements.get((row["smiles"], row["design"]), {}).get("trials", [])
        samples = [trial["score"] for trial in trials if trial.get("status") == "success"] or [row["score"]]
        candidates.append(
            {
                "smiles": row["smiles"],
                "score_samples": samples,
                "score_summary": {k: v for k, v in row.items() if k.startswith("seed_score_")},
                "route": route_record.get("route"),
                "route_evidence": route_record.get("verification"),
                "synthesis_evidence": route_record.get("synthesis_evidence"),
                "design": row["design"],
            }
        )
    return {
        "candidates": candidates,
        "evaluations": summary["scored_products"],
        "evaluation_unit": "synthesis product/design 2D screens"
        if settings.get("disable_3d", False)
        else "synthesis product/design 3D screens (refinement separately reported)",
        "metadata": {
            "generator_cost": cost,
            "ablation": "3d_scoring_disabled" if settings.get("disable_3d", False) else None,
            "synthesis_summary": summary,
            "refinement": refinement,
            "budget_accounting": "partial; generator and refinement costs require wall-time comparison",
        },
    }


ADAPTERS = {
    "nearest_catalog": nearest_catalog,
    "nearest_3d": nearest_3d,
    "reaction_enumeration": reaction_enumeration,
    "synllama_only": synllama_only,
    "pipeline": pipeline,
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--response", type=Path, required=True)
    args = parser.parse_args(argv)
    request = json.loads(args.request.read_text())
    response = ADAPTERS[request["kind"]](request)
    write_json(args.response, response)


if __name__ == "__main__":
    main()
