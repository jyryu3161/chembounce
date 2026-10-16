"""Bounded original-hit objective for complete regenerated route products only."""

from collections import OrderedDict

from rdkit import DataStructs

from .common import FP, canonical, molecule, positive


def search_policy(config):
    default = (
        "original-hit-v1" if config.get("recommendation", {}).get("mode") == "qualified" else "design-v1"
    )
    version = config.get("route_search_policy", default)
    if version not in ("original-hit-v1", "design-v1"):
        raise ValueError("route_search_policy must be original-hit-v1 or design-v1")
    size = config.get("route_search_cache_size", 4096)
    positive(size, "route_search_cache_size")
    guided = version == "original-hit-v1"
    return {
        "version": version,
        "cache_size_per_design": size if guided else 0,
        "audit_sample_size_per_design": min(size, 256) if guided else 0,
        "retain_seed_only": guided and config.get("expansion", {}).get("rounds", 0) > 0,
        "scope": "Complete regenerated root products only; intermediate and building-block targets unchanged",
        "ranking": (
            "eligible first, original-hit binary Morgan, design binary Morgan, deterministic route key; distinct products within each eligibility tier"
            if guided
            else "Legacy node/model target binary Morgan for projection and rebuild; design target for expansion frontier; no original-hit objective or seed-only distinction"
        ),
        "fingerprint": "binary Morgan radius=2 bits=2048 chirality=true; not sparse-count Morgan",
        "final_acceptance": "Unchanged final recommendation contract; design-v1 is an explicit search ablation",
    }


class RootProductObjective:
    def __init__(self, hit, design, minimum_scores, check_constraints, policy):
        self.hit, self.design = canonical(hit), canonical(design)
        self.hit_fp = FP.GetFingerprint(molecule(self.hit))
        self.design_fp = FP.GetFingerprint(molecule(self.design))
        self.minimum = {
            key: value for key, value in minimum_scores.items() if key in {"hit_tanimoto", "design_tanimoto"}
        }
        self.check_constraints = check_constraints
        self.policy = dict(policy)
        self.cache, self.sample = OrderedDict(), {}
        self.counts = dict(
            requests=0,
            evaluations=0,
            cache_hits=0,
            evictions=0,
            eligible_evaluations=0,
            rejected_evaluations=0,
        )

    def assess(self, product, report=None):
        self.counts["requests"] += 1
        if product in self.cache:
            self.counts["cache_hits"] += 1
            if report is not None:
                report.root_objective_cache_hits += 1
            self.cache.move_to_end(product)
            return self.cache[product]
        self.counts["evaluations"] += 1
        if report is not None:
            report.root_objective_evaluations += 1
        fp = FP.GetFingerprint(molecule(product))
        scores = {
            "hit_tanimoto": DataStructs.TanimotoSimilarity(self.hit_fp, fp),
            "design_tanimoto": DataStructs.TanimotoSimilarity(self.design_fp, fp),
        }
        reasons = [
            f"{key} below configured minimum" for key, value in self.minimum.items() if scores[key] < value
        ]
        try:
            self.check_constraints(product)
        except ValueError as exc:
            reasons.append(str(exc))
        decision = {"product": product, "eligible": not reasons, **scores, "reasons": reasons}
        self.counts["rejected_evaluations" if reasons else "eligible_evaluations"] += 1
        self.cache[product] = decision
        if len(self.cache) > self.policy["cache_size_per_design"]:
            self.cache.popitem(last=False)
            self.counts["evictions"] += 1
        if len(self.sample) < self.policy["audit_sample_size_per_design"]:
            self.sample.setdefault(product, decision)
        return decision

    def snapshot(self):
        return {
            "hit": self.hit,
            "design": self.design,
            "policy": self.policy,
            "minimum_scores": self.minimum,
            "counts": dict(self.counts),
            "current_cached_unique_products": len(self.cache),
            "costs": {
                "product_fingerprint_calculations": self.counts["evaluations"],
                "constraint_checks": self.counts["evaluations"],
                "reference_fingerprint_calculations": 2,
            },
            "decision_sample": list(self.sample.values()),
            "audit_scope": "First bounded sample of exact products; not a complete product inventory. Evaluations count cache misses and may repeat products after eviction; requests include repeated pruning visits.",
            "sample_at_capacity": len(self.sample) == self.policy["audit_sample_size_per_design"],
        }


def prioritize_products(items, key, product):
    """Distinct exact products first, then alternate routes, within a ranked tier."""
    unique, alternatives = {}, []
    for item in sorted(items, key=key):
        smiles = product(item)
        if smiles in unique:
            alternatives.append(item)
        else:
            unique[smiles] = item
    return list(unique.values()) + alternatives
