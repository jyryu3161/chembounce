"""Cheap necessary condition for the existing typed polar-recall constraint.

This does not predict geometric feasibility or replace final pose/route evidence.
Only complete final products may be filtered; intermediates can gain features.
"""

from collections import Counter, OrderedDict
import copy
import time

from .common import canonical, molecule, positive
from .similarity_constraints import POLAR, feature_records


class CapacityRejection(ValueError):
    def __init__(self, decision):
        self.decision = copy.deepcopy(decision)
        super().__init__(decision["reason"])


class PolarCapacity:
    def __init__(self, similarity, *, cache_size=4096, deadline_reached=lambda: False):
        positive(cache_size, "route_search_cache_size")
        self.similarity = similarity
        self.deadline_reached = deadline_reached
        self.cache_size = cache_size
        self.sample_size = min(cache_size, 256)
        self.cache, self.sample = OrderedDict(), {}
        self.counts = dict(
            requests=0,
            evaluations=0,
            cache_hits=0,
            evictions=0,
            preflight_failures=0,
            capacity_rejections=0,
            deadline_events=0,
        )
        self.elapsed_seconds = 0.0
        reference_started = time.monotonic()
        self.reference_ready = not self.deadline_reached()
        self.hit_counts = Counter()
        perceptions = 0
        if self.reference_ready:
            self.hit_counts.update(
                f["family"] for f in feature_records(molecule(similarity.hit)) if f["family"] in POLAR
            )
            perceptions = 1
            self.reference_ready = not self.deadline_reached()
        self.reference_elapsed_seconds = time.monotonic() - reference_started
        self.policy = {
            "version": "typed-polar-capacity-v1",
            "hit": similarity.hit,
            "constraint_policy_sha256": similarity.sha256,
            "minimum_polar_reference_recall": similarity.policy["minimum_polar_reference_recall"],
            "reference_counts": dict(sorted(self.hit_counts.items())),
            "reference_feature_perceptions": perceptions,
            "reference_status": "complete" if self.reference_ready else "not_evaluated_deadline",
            "cache_size_per_run": cache_size,
            "sample_size": self.sample_size,
            "scope": "Complete root products only; necessary condition, never geometric acceptance. Optional landmarks remain exclusively subject to final joint validation.",
        }

    def _deadline(self, decision):
        self.counts["deadline_events"] += 1
        return {
            **decision,
            "eligible": False,
            "status": "not_evaluated_deadline",
            "reason": "Deadline reached before polar-capacity decision could be used; completed work is retained",
            "preflight_completed": decision.get("status")
            in {"preflight_rejected", "possible", "not_applicable", "capacity_rejected"},
            "capacity_computed": "recall_upper_bound" in decision,
        }

    def assess(self, product):
        started = time.monotonic()
        try:
            self.counts["requests"] += 1
            product = canonical(product)
            decision = {"product": product, "constraint_policy_sha256": self.similarity.sha256}
            if not self.reference_ready or self.deadline_reached():
                return self._deadline(decision)
            if product in self.cache:
                self.counts["cache_hits"] += 1
                self.cache.move_to_end(product)
                return copy.deepcopy(self.cache[product])
            self.counts["evaluations"] += 1
            try:
                features = self.similarity.preflight_candidate(product)
            except ValueError as exc:
                self.counts["preflight_failures"] += 1
                decision.update(eligible=False, status="preflight_rejected", reason=str(exc))
            else:
                counts = Counter(f["family"] for f in features if f["family"] in POLAR)
                denominator = sum(self.hit_counts.values())
                numerator = sum(min(n, counts[family]) for family, n in self.hit_counts.items())
                # Same float division/strict comparison as the final recall.
                # Exact Fraction comparison would incorrectly reject e.g. 1/10
                # against the binary floating-point representation of 0.1.
                upper = numerator / denominator if denominator else None
                minimum = self.policy["minimum_polar_reference_recall"]
                eligible = upper is None or upper >= minimum
                self.counts["capacity_rejections"] += int(not eligible)
                decision.update(
                    eligible=eligible,
                    status="not_applicable"
                    if upper is None
                    else "possible"
                    if eligible
                    else "capacity_rejected",
                    reason=None if eligible else "polar_reference_capacity_below_minimum",
                    reference_counts=copy.deepcopy(self.policy["reference_counts"]),
                    candidate_counts=dict(sorted(counts.items())),
                    numerator=numerator,
                    denominator=denominator,
                    recall_upper_bound=upper,
                    minimum_polar_reference_recall=minimum,
                )
            if self.deadline_reached():
                return self._deadline(decision)
            self.cache[product] = copy.deepcopy(decision)
            if len(self.cache) > self.cache_size:
                self.cache.popitem(last=False)
                self.counts["evictions"] += 1
            if len(self.sample) < self.sample_size:
                self.sample.setdefault(product, copy.deepcopy(decision))
            return decision
        finally:
            self.elapsed_seconds += time.monotonic() - started

    def require(self, product):
        decision = self.assess(product)
        if not decision["eligible"]:
            raise CapacityRejection(decision)
        return decision

    def snapshot(self):
        return {
            "policy": copy.deepcopy(self.policy),
            "counts": dict(self.counts),
            "elapsed_seconds": self.elapsed_seconds,
            "reference_elapsed_seconds": self.reference_elapsed_seconds,
            "current_cached_products": len(self.cache),
            "decision_sample": copy.deepcopy(list(self.sample.values())),
            "sample_at_capacity": len(self.sample) == self.sample_size,
            "audit_scope": "Evaluations count cache misses, including repeat products after eviction; rejection counts are evaluation events, not unique products. First bounded sample only. No model, embedding, alignment or route-result cache.",
        }
