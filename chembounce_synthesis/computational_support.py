"""Fresh concrete route replay plus pinned per-step forward selection support."""

import copy
import math
import time

from .evidence import EvidencePolicy, SynthesisEvidenceEvaluator, _json_key
from .forward_support import PROFILE, PinnedReactionT5Forward, binding, strict_smiles

STRICT_POLICY = {"require_forward_agreement": True}


def validate_config(config):
    from .forward_support import validate_settings

    validate_settings(config.get("predictor"))
    if config.get("product_matching", "exact") != "exact":
        raise ValueError(f"{PROFILE} requires exact full-product matching")
    policy = EvidencePolicy(**config.get("policy", STRICT_POLICY))
    if not policy.require_forward_agreement:
        raise ValueError(f"{PROFILE} requires every-step forward agreement")
    if policy.require_complete_context or policy.require_observed_success:
        raise ValueError("Verified physical context/empirical evidence is unavailable in this profile")
    if config.get("observations"):
        raise ValueError("Supplied observation claims cannot establish computational route support")


class ComputationalRouteSupport(SynthesisEvidenceEvaluator):
    def __init__(self, predictor=None, **kwargs):
        kwargs.setdefault("policy", STRICT_POLICY)
        super().__init__(predictor=predictor, **kwargs)
        self.provenance.update(
            profile=PROFILE,
            interpretation="Conditional computational selection support; not independent post-selection or experimental validation; training independence unknown",
            physical_context="Temperature/time are not model inputs; missing physical context remains conditional",
        )

    def evaluate(self, verification, step_context=None, deadline_monotonic=None):
        raise ValueError("Strict support requires evaluate_route with the actual route, catalog and rules")

    def _raw_prediction(self, reactants, context, deadline, *, cached_only=False):
        if self.predictor is None:
            return {"status": "unknown", "reason": "Pinned forward model unavailable"}
        if deadline is not None and time.monotonic() >= deadline:
            return {"status": "unknown", "reason": "Forward prediction deadline exhausted"}
        reagents = [strict_smiles(s) for s in context.get("reagents", [])]
        try:
            prepared = self.predictor.prepare(reactants, reagents)
        except (ValueError, RuntimeError) as exc:
            return {"status": "unknown", "reason": str(exc)}
        key = binding({"input": prepared, "predictor": self.predictor.provenance})
        if key in self.cache:
            if not cached_only:
                self.usage["cache_hits"] += 1
            return copy.deepcopy(self.cache[key])
        if cached_only:
            raise ValueError("Strict evidence has no actual model-call cache binding")
        if deadline is not None and time.monotonic() >= deadline:
            return {"status": "unknown", "reason": "Forward prediction deadline exhausted"}
        reserve = 3 * 256
        if (
            self.usage["predictor_calls"] >= self.limits["predictor_calls"]
            or self.usage["input_tokens"] + prepared["input_tokens"] > self.limits["input_tokens"]
            or self.usage["output_tokens"] + reserve > self.limits["output_tokens"]
        ):
            return {"status": "unknown", "reason": "Forward prediction budget exhausted"}
        self.usage["predictor_calls"] += 1
        self.usage["input_tokens"] += prepared["input_tokens"]
        try:
            raw = self.predictor.predict_raw(prepared)
            self.predictor.validate_raw(raw, prepared)
            tokens = raw["output_tokens"]
            if isinstance(tokens, bool) or not isinstance(tokens, int) or not 0 <= tokens <= reserve:
                raise ValueError("Invalid strict output-token count")
            self.usage["output_tokens"] += tokens
            result = {"status": "completed_call", "raw": raw}
        except (ValueError, RuntimeError, KeyError, TypeError) as exc:
            self.usage["output_tokens"] += reserve
            result = {"status": "failed", "reason": str(exc)}
        self.cache[key] = copy.deepcopy(result)
        return result

    def _report(self, audit, deadline=None, *, cached_only=False):
        reports, reasons = [], []
        for step in audit["steps"]:
            reactants = [slot["smiles"] for slot in step["reactant_slots"]]
            product = step["product"]
            context = copy.deepcopy(self.reaction_contexts.get(_json_key([sorted(reactants), product]), {}))
            prediction = self._raw_prediction(reactants, context, deadline, cached_only=cached_only)
            matches = []
            if prediction["status"] == "completed_call":
                raw = prediction["raw"]
                prepared = self.predictor.prepare(
                    reactants, [strict_smiles(s) for s in context.get("reagents", [])]
                )
                predictions = self.predictor.validate_raw(raw, prepared)
                matches = [
                    item["rank"] for item in predictions if item["complete"] and item["product"] == product
                ]
                status = (
                    "matched"
                    if matches
                    else "mismatched"
                    if all(item["complete"] and item["product"] is not None for item in predictions)
                    else "partially_evaluated"
                )
                prediction.update(
                    status=status,
                    matched_rank=min(matches) if matches else None,
                    matching_policy="exact",
                    predictions=predictions,
                )
            rejected = []
            if not step["strict_eligible"]:
                rejected.append("Concrete atom/stereo audit did not establish route integrity")
            if not matches:
                rejected.append("No completed exact full-product hypothesis in bounded top three")
            reasons.extend(f"{step['node_id']}: {reason}" for reason in rejected)
            reports.append(
                {
                    "node_id": step["node_id"],
                    "reactants": reactants,
                    "product": product,
                    "step_binding_sha256": step["step_binding_sha256"],
                    "rule_replay": "matched",
                    "forward_prediction": prediction,
                    "context": {
                        "status": "supplied_unvalidated" if context else "incomplete",
                        "supplied": context,
                        "validated": False,
                        "used_by_predictor": ["reagents"] if "reagents" in context else [],
                    },
                    "atom_sources": {
                        "status": "concrete_replay_verified" if step["strict_eligible"] else "unsupported"
                    },
                    "observations": {"status": "unavailable", "source_audited": False},
                    "eligibility": {"passed": not rejected, "unmet_requirements": rejected},
                }
            )
        if not audit["strict_eligible"]:
            reasons.append("Route audit did not establish complete concrete atom/stereo integrity")
        return {
            "schema_version": 2,
            "profile": PROFILE,
            "evidence_level": "conditional_computational_selection_support",
            "product": audit["product"],
            "route_binding_sha256": audit["route_binding_sha256"],
            "route_audit": audit,
            "steps": reports,
            "eligibility": {"passed": not reasons, "reasons": reasons},
            "policy": self.provenance["policy"],
            "predictor_provenance": self.provenance["predictor"],
            "interpretation": self.provenance["interpretation"],
        }

    def evaluate_route(self, route, catalog, rules, product_limit=64, deadline_monotonic=None):
        from .route_audit import audit_route

        if deadline_monotonic is not None and (
            isinstance(deadline_monotonic, bool)
            or not isinstance(deadline_monotonic, (int, float))
            or not math.isfinite(deadline_monotonic)
        ):
            raise ValueError("Evidence deadline must be a finite monotonic timestamp")
        report = self._report(audit_route(route, catalog, rules, product_limit), deadline_monotonic)
        return {**report, "usage": self.usage.copy()}

    def validate_route(self, evidence, route, catalog, rules, product_limit=64):
        """Fresh replay and raw-token comparison, without another model generation.

        The saved report must reproduce the evaluator's actual input-bound calls.
        Supplied match flags, route hashes or plausible tokens alone are insufficient.
        """
        from .route_audit import audit_route

        expected = self._report(audit_route(route, catalog, rules, product_limit), cached_only=True)
        if not expected["eligibility"]["passed"]:
            raise ValueError("Strict route support did not pass fresh final validation")
        if {key: value for key, value in evidence.items() if key != "usage"} != expected:
            raise ValueError("Strict evidence differs from fresh route/input/raw-prediction validation")
        return True


def from_config(config, *, load_predictor=True):
    validate_config(config)
    values = copy.deepcopy(config)
    values.pop("profile", None)
    provider = values.pop("predictor")
    predictor = (
        PinnedReactionT5Forward(**{k: v for k, v in provider.items() if k != "kind"})
        if load_predictor
        else None
    )
    return ComputationalRouteSupport(predictor=predictor, **values)
