"""Posthoc saved-product allocation comparator, not a production screening policy.

Selection reads only exact product identities and saved original-hit 2D scores.
Existing three-seed labels are joined after both fixed-quota selections exist.
No geometry, model, catalog, route generation, or forward evaluation is run.
"""

import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path

SEEDS = (42, 43, 44)
QUOTA = 12
POLICY = "original-hit-tanimoto-desc-first-eligible-route-v1"


def select(rows, original):
    """2D-only ranking; ties retain saved first eligible route order."""
    products = [r["product"] for r in rows]
    if len(set(products)) != len(products) or len(original) != QUOTA or len(set(original)) != QUOTA:
        raise ValueError("Allocation requires distinct exact products and12original selections")
    if not set(original) <= set(products):
        raise ValueError("Original allocation outside eligible population")
    for row in rows:
        score = row["hit_tanimoto"]
        if type(score) not in (int, float) or not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError("Invalid saved original-hit similarity")
    ranked = sorted(rows, key=lambda r: (-r["hit_tanimoto"], r["first_eligible_route_index"], r["product"]))
    return {"original_saved": list(original), POLICY: [r["product"] for r in ranked[:QUOTA]]}


def read_observation(binding, candidate):
    """Consume pinned saved results only; missing/inconsistent evidence is unknown."""
    try:
        raw = Path(binding["path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != binding["sha256"]:
            raise ValueError("Saved condition digest mismatch")
        condition = json.loads(raw)
        for field in ("product", "hit", "design"):
            if condition[field] != candidate[field]:
                raise ValueError("Saved condition identity/representative mismatch")
        if condition["seed"] != binding["seed"]:
            raise ValueError("Saved condition seed mismatch")
        observed = condition["baseline"] if binding["origin"] == "paired_baseline" else condition
        status = observed["status"]
        if status == "feasible":
            if not isinstance(observed.get("info"), dict):
                raise ValueError("Missing saved feasible result")
            return {"status": "feasible"}
        if status in {"abstained", "bounded_infeasible"}:
            audit = observed.get("audit", {})
            if not isinstance(audit, dict):
                raise ValueError("Malformed saved pair audit")
            pairs = audit.get("pairs", [])
            if not isinstance(pairs, list) or any(not isinstance(p, dict) for p in pairs):
                raise ValueError("Malformed saved audit pairs")
            if (
                observed.get("reason") == "no_jointly_feasible_hit_pose"
                and pairs
                and all(p.get("status") == "rejected" for p in pairs)
            ):
                return {"status": "bounded_infeasible"}
        return {"status": "unknown", "reason": "Saved result is not a complete bounded assessment"}
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return {"status": "unknown", "reason": f"{type(exc).__name__}: {exc}"}


def join_outcomes(candidate):
    seeds = []
    for seed in SEEDS:
        bindings = [b for b in candidate["observations"] if b["seed"] == seed]
        observations = [read_observation(b, candidate) for b in bindings]
        statuses = {o["status"] for o in observations}
        status = next(iter(statuses)) if len(statuses) == 1 and "unknown" not in statuses else "unknown"
        seeds.append({"seed": seed, "status": status, "observations": observations})
    complete = all(s["status"] != "unknown" for s in seeds)
    return {
        "seeds": seeds,
        "assessment_complete": complete,
        "all_three_feasible": all(s["status"] == "feasible" for s in seeds) if complete else None,
    }


def describe(smiles):
    # Existing2D identity normalization only. No conformer/scorer imports or calls.
    from chembounce_synthesis.identity import chemical_identity

    return asdict(chemical_identity(smiles))


def diversity(products, descriptors, hit):
    known = [descriptors[p] for p in products if "error" not in descriptors[p]]
    result = dict(
        products=len(products),
        identity_observed=len(known),
        identity_unknown=len(products) - len(known),
        observed_families=len({v["family_smiles"] for v in known}),
        observed_atom_aware_scaffolds=len({v["scaffold_family_smiles"] for v in known if v["has_scaffold"]}),
        observed_generic_ring_graphs=len({v["ring_skeleton_smiles"] for v in known if v["has_scaffold"]}),
        observed_acyclic_products=sum(not v["has_scaffold"] for v in known),
    )
    if "error" not in hit:
        result.update(
            observed_atom_aware_changes=sum(
                v["scaffold_family_smiles"] != hit["scaffold_family_smiles"] for v in known
            ),
            observed_generic_graph_changes=sum(
                v["ring_skeleton_smiles"] != hit["ring_skeleton_smiles"] for v in known
            ),
        )
    return result


def compare(rows, original):
    allocations = select(rows, original)  # Deliberately before any label or2Didentity read.
    outcomes = {r["product"]: join_outcomes(r) for r in rows}
    descriptors = {}
    for smiles in dict.fromkeys([rows[0]["hit"], *[r["product"] for r in rows]]):
        try:
            descriptors[smiles] = describe(smiles)
        except Exception as exc:
            descriptors[smiles] = {"error": f"{type(exc).__name__}: {exc}"}
    arms = {}
    for name, products in allocations.items():
        feasible = [p for p in products if outcomes[p]["all_three_feasible"] is True]
        unknown = [p for p in products if not outcomes[p]["assessment_complete"]]
        arms[name] = dict(
            products=products,
            quota=QUOTA,
            eligible_products=len(rows),
            exact_product_coverage=len(products) / len(rows),
            all_three_feasible_observed=len(feasible),
            feasible_products=feasible,
            unknown_products=unknown,
            assessment_complete=not unknown,
            selected_diversity=diversity(products, descriptors, descriptors[rows[0]["hit"]]),
            feasible_diversity=diversity(feasible, descriptors, descriptors[rows[0]["hit"]]),
        )
    return dict(
        arms=arms,
        outcomes=outcomes,
        identities=descriptors,
        full_population_assessment_complete=all(o["assessment_complete"] for o in outcomes.values()),
        shared_products=sorted(set(allocations["original_saved"]) & set(allocations[POLICY])),
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    raw = args.plan.read_bytes()
    if hashlib.sha256(raw).hexdigest() != args.plan_sha256:
        raise ValueError("Plan byte mismatch")
    plan = json.loads(raw)
    if plan["policy"] != POLICY or plan["quota"] != QUOTA or plan["seeds"] != list(SEEDS):
        raise ValueError("Unsupported fixed allocation comparator")
    root = Path(__file__).resolve().parents[2]
    for path, expected in plan["source_sha256"].items():
        if hashlib.sha256((root / path).read_bytes()).hexdigest() != expected:
            raise ValueError("Source byte mismatch: " + path)
    from rdkit import rdBase

    if rdBase.rdkitVersion != plan["rdkit_version"]:
        raise ValueError("2Didentity runtime version mismatch")
    args.output.mkdir(parents=True, exist_ok=False)
    result = {
        "plan_sha256": args.plan_sha256,
        "policy": POLICY,
        "scope": "Posthoc saved fixed-population comparator; not fresh method performance or production policy",
        "inputs": {},
    }
    for item in plan["inputs"]:
        result["inputs"][item["attempt"]] = compare(item["population"], item["original_products"])
    (args.output / "comparison.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
