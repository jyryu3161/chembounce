"""Independent, explicitly conditional evidence for synthesis route proposals.

Reaction-rule replay, forward-model agreement, supplied chemical context, and
reported observations are orthogonal evidence. None is a synthesis probability.
The forward predictor never receives the claimed product or reaction template.
"""

import copy
import importlib.metadata
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from rdkit import Chem

from .common import absolute_smiles, parse_absolute_smiles, canonical, digest, positive


@dataclass(frozen=True)
class EvidencePolicy:
    require_forward_agreement: bool = False
    require_complete_context: bool = False
    require_explicit_atom_sources: bool = False
    require_observed_success: bool = False

    def __post_init__(self):
        for key, value in asdict(self).items():
            if not isinstance(value, bool):
                raise ValueError(f"{key} must be true or false")


def evidence_requirements(prediction_status, context_complete, atom_source_status, observed_status):
    """Shared policy checks for evaluation and final recommendation qualification."""
    return {
        "require_forward_agreement": prediction_status == "matched",
        "require_complete_context": context_complete,
        "require_explicit_atom_sources": atom_source_status
        in {"explicit_in_reaction_rule", "supplied_mapped_attribution"},
        "require_observed_success": observed_status == "success",
    }


def matching_prediction_ranks(predictions, product, product_matching):
    """Exact canonical prediction matches used by evaluation and qualification."""
    return [
        p["rank"]
        for p in predictions
        if p.get("product") == product
        or (product_matching == "component" and product in p.get("components", []))
    ]


def _json_key(value):
    try:
        return json.dumps(value, sort_keys=True, allow_nan=False, separators=(",", ":"))
    except (ValueError, TypeError) as exc:
        raise ValueError("Evidence context must contain finite JSON values") from exc


def _component_smiles(smiles):
    """Canonicalize a reagent, allowing explicitly supplied salt components."""
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Reagents require nonempty SMILES")
    mol = parse_absolute_smiles(smiles)
    if mol is None or not mol.GetNumAtoms() or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError(f"Invalid reagent SMILES: {smiles!r}")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return absolute_smiles(mol, isomericSmiles=True)


def _smiles_list(values, name, *, reagents=False):
    if not isinstance(values, list):
        raise ValueError(f"{name} must be a list of SMILES")
    convert = _component_smiles if reagents else canonical
    return [convert(value) for value in values]


def _context(value):
    if value is None:
        value = {}
    if not isinstance(value, dict):
        raise ValueError("Each step context must be an object")
    allowed = {"reagents", "conditions", "source", "atom_sources"}
    if set(value) - allowed:
        raise ValueError(f"Unknown step-context fields: {sorted(set(value) - allowed)}")
    result = copy.deepcopy(value)
    if "reagents" in result:
        _smiles_list(result["reagents"], "reagents", reagents=True)
        # Atom-source indices refer to the supplied representation. Canonicalize
        # only comparison/model inputs, never reorder these indexed source atoms.
        result["reagents"] = [smiles.strip() for smiles in result["reagents"]]
    if "conditions" in result:
        if not isinstance(result["conditions"], dict):
            raise ValueError("conditions must be an object")
        _json_key(result["conditions"])
    if "source" in result and (not isinstance(result["source"], str) or not result["source"].strip()):
        raise ValueError("Context source must be a nonempty locator or attribution")
    if "atom_sources" in result and not isinstance(result["atom_sources"], list):
        raise ValueError("atom_sources must be a list")
    return result


def _reaction_key(reactants, product, context):
    """Multiset identity: preserve duplicate reactants and exact chemical states."""
    return _json_key(
        {
            "reactants": sorted(reactants),
            "product": product,
            "reagents": (
                sorted(_component_smiles(s) for s in context["reagents"]) if "reagents" in context else None
            ),
            "conditions": context.get("conditions"),
        }
    )


def _reaction_contexts(records):
    if not isinstance(records, list):
        raise ValueError("reaction_contexts must be a list")
    indexed = {}
    for record in records:
        if not isinstance(record, dict) or not {"reactants", "product"} <= set(record):
            raise ValueError("Reaction context requires exact reactants and product selectors")
        reactants = _smiles_list(record["reactants"], "context reactants")
        if not reactants:
            raise ValueError("Reaction context requires reactants")
        key = _json_key([sorted(reactants), canonical(record["product"])])
        context = _context(
            {key: value for key, value in record.items() if key not in {"reactants", "product"}}
        )
        if key in indexed:
            raise ValueError("Multiple contexts match the same exact reaction; choose one explicitly")
        indexed[key] = context
    return indexed


def _observations(records):
    if not isinstance(records, list):
        raise ValueError("observations must be a list")
    indexed = {}
    for record in records:
        required = {
            "kind",
            "source",
            "record_id",
            "reactants",
            "product",
            "reagents",
            "conditions",
            "outcome",
        }
        if not isinstance(record, dict) or not required <= set(record):
            raise ValueError(
                "Observation records require kind/source/record_id and exact reaction/context/outcome"
            )
        if record["kind"] not in {"literature", "experiment"} or record["outcome"] not in {
            "success",
            "failure",
        }:
            raise ValueError("Observation kind or outcome is invalid")
        for name in ("source", "record_id"):
            if not isinstance(record[name], str) or not record[name].strip():
                raise ValueError(f"Observation {name} must be nonempty")
        context = _context({"reagents": record["reagents"], "conditions": record["conditions"]})
        if not context["conditions"]:
            raise ValueError("Observation records require specified conditions")
        reactants = _smiles_list(record["reactants"], "observation reactants")
        if not reactants:
            raise ValueError("Observation records require reactants")
        product = canonical(record["product"])
        normalized = {**record, "reactants": reactants, "product": product, **context}
        _json_key(normalized)
        indexed.setdefault(_reaction_key(reactants, product, context), []).append(normalized)
    return indexed


def _atom_source_evidence(step, context):
    implicit = step.get("unmapped_product_template_atoms", [])
    if not implicit:
        return {"status": "explicit_in_reaction_rule", "unattributed_template_atoms": [], "attributions": []}
    expected = {atom["atom_index"]: atom["element"] for atom in implicit}
    mapped, used, accepted = set(), set(), []
    for attribution in context.get("atom_sources", []):
        required = {"template_atom_index", "reagent_index", "reagent_atom_index", "source", "record_id"}
        if not isinstance(attribution, dict) or set(attribution) != required:
            raise ValueError("Atom-source attribution requires atom indices, source and record_id")
        for key in ("template_atom_index", "reagent_index", "reagent_atom_index"):
            if (
                isinstance(attribution[key], bool)
                or not isinstance(attribution[key], int)
                or attribution[key] < 0
            ):
                raise ValueError("Atom-source indices must be nonnegative integers")
        for key in ("source", "record_id"):
            if not isinstance(attribution[key], str) or not attribution[key].strip():
                raise ValueError("Atom-source attributions require source and record_id")
        target = attribution["template_atom_index"]
        origin = (attribution["reagent_index"], attribution["reagent_atom_index"])
        if target not in expected or target in mapped or origin in used:
            raise ValueError("Atom-source attribution duplicates or references an unknown template atom")
        reagents = context.get("reagents", [])
        if origin[0] >= len(reagents):
            raise ValueError("Atom-source attribution references an absent reagent")
        reagent = Chem.MolFromSmiles(reagents[origin[0]])
        if (
            origin[1] >= reagent.GetNumAtoms()
            or reagent.GetAtomWithIdx(origin[1]).GetSymbol() != expected[target]
        ):
            raise ValueError("Atom-source attribution has an invalid reagent atom or element")
        mapped.add(target)
        used.add(origin)
        accepted.append(attribution)
    missing = sorted(set(expected) - mapped)
    return {
        "status": "unknown" if missing else "supplied_mapped_attribution",
        "unattributed_template_atoms": missing,
        "attributions": accepted,
        "source_audited": False,
    }


class SynthesisEvidenceEvaluator:
    """Reuse independent predictions across routes under explicit global budgets.

    A predictor exposes ``provenance``, ``token_budget(reactants, reagents,
    conditions)`` and ``predict(...)``. Prediction inputs never include the
    candidate product. Call/input/output limits apply to this evaluator's whole
    lifetime. A cache hit does not consume inference budget.
    """

    def __init__(
        self,
        predictor=None,
        policy=None,
        max_predictor_calls=128,
        max_input_tokens=32768,
        max_output_tokens=32768,
        observations=None,
        product_matching="exact",
        reaction_contexts=None,
    ):
        self.predictor = predictor
        if policy is not None and not isinstance(policy, (EvidencePolicy, dict)):
            raise ValueError("Evidence policy must be an object")
        self.policy = policy if isinstance(policy, EvidencePolicy) else EvidencePolicy(**(policy or {}))
        if product_matching not in {"exact", "component"}:
            raise ValueError("product_matching must be exact or component")
        self.product_matching = product_matching
        self.limits = {
            "predictor_calls": max_predictor_calls,
            "input_tokens": max_input_tokens,
            "output_tokens": max_output_tokens,
        }
        for key, value in self.limits.items():
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"Evidence {key} budget must be a nonnegative integer")
        self.usage = {"predictor_calls": 0, "input_tokens": 0, "output_tokens": 0, "cache_hits": 0}
        self.cache = {}
        self.observations = _observations([] if observations is None else observations)
        self.reaction_contexts = _reaction_contexts([] if reaction_contexts is None else reaction_contexts)
        self.provenance = {
            "profile": "legacy-conditional-v1",
            "predictor": copy.deepcopy(predictor.provenance) if predictor is not None else None,
            "policy": asdict(self.policy),
            "product_matching": product_matching,
            "limits": self.limits.copy(),
            "observation_record_count": sum(len(records) for records in self.observations.values()),
            "reaction_context_count": len(self.reaction_contexts),
            "interpretation": "Conditional product agreement and supplied evidence; not synthesis probability",
        }

    def _predict(self, reactants, context, deadline_monotonic=None):
        if self.predictor is None:
            return {
                "status": "unknown",
                "reason": "No independent forward predictor configured",
                "cached": False,
            }
        reagents, conditions = context.get("reagents", []), context.get("conditions", {})
        key = _json_key([sorted(reactants), sorted(_component_smiles(s) for s in reagents), conditions])
        if key in self.cache:
            self.usage["cache_hits"] += 1
            return {**copy.deepcopy(self.cache[key]), "cached": True}
        if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
            return {
                "status": "unknown",
                "reason": "Independent predictor deadline exhausted",
                "cached": False,
            }
        if self.usage["predictor_calls"] >= self.limits["predictor_calls"]:
            return {
                "status": "unknown",
                "reason": "Independent predictor call budget exhausted",
                "cached": False,
            }
        try:
            budget = self.predictor.token_budget(reactants, reagents, conditions)
            for name in ("input_tokens", "max_output_tokens"):
                positive(budget[name], name)
        except (ValueError, RuntimeError, KeyError, TypeError) as exc:
            return {"status": "failed", "reason": str(exc), "cached": False}
        if self.usage["input_tokens"] + budget["input_tokens"] > self.limits["input_tokens"]:
            return {
                "status": "unknown",
                "reason": "Independent predictor input-token budget exhausted",
                "cached": False,
            }
        if self.usage["output_tokens"] + budget["max_output_tokens"] > self.limits["output_tokens"]:
            return {
                "status": "unknown",
                "reason": "Independent predictor output-token budget exhausted",
                "cached": False,
            }
        if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
            return {
                "status": "unknown",
                "reason": "Independent predictor deadline exhausted",
                "cached": False,
            }
        self.usage["predictor_calls"] += 1
        self.usage["input_tokens"] += budget["input_tokens"]
        try:
            predicted = self.predictor.predict(reactants, reagents, conditions)
            tokens = predicted["output_tokens"]
            if (
                isinstance(tokens, bool)
                or not isinstance(tokens, int)
                or not 0 <= tokens <= budget["max_output_tokens"]
            ):
                raise ValueError("Forward predictor returned an invalid output-token count")
            products, complete = predicted["products"], predicted["complete"]
            if (
                not isinstance(products, list)
                or not isinstance(complete, list)
                or len(products) != len(complete)
            ):
                raise ValueError("Forward predictor requires products and per-product completion flags")
            predictions = []
            for rank, (raw, finished) in enumerate(zip(products, complete), start=1):
                if not isinstance(raw, str) or not isinstance(finished, bool):
                    raise ValueError("Forward predictor returned malformed products or completion flags")
                item = {"rank": rank, "raw_product": raw, "complete": finished, "product": None}
                try:
                    if not finished:
                        raise ValueError("Product decoding exhausted its output-token limit")
                    item["product"] = _component_smiles(raw)
                    item["components"] = sorted(
                        absolute_smiles(fragment, isomericSmiles=True)
                        for fragment in Chem.GetMolFrags(Chem.MolFromSmiles(item["product"]), asMols=True)
                    )
                    if len(item["components"]) > 1 and self.product_matching == "exact":
                        item["comparison_note"] = "unsupported_mixture_for_exact_matching"
                except ValueError as exc:
                    item["invalid_reason"] = str(exc)
                predictions.append(item)
            self.usage["output_tokens"] += tokens
            result = {
                "status": "predicted" if any(item["product"] for item in predictions) else "unknown",
                "predictions": predictions,
                "input_tokens": budget["input_tokens"],
                "output_tokens": tokens,
                "cached": False,
            }
        except (ValueError, RuntimeError, KeyError, TypeError) as exc:
            # Failed calls may have generated tokens before raising. Charge their
            # reserved upper bound so repeated failures cannot evade the budget.
            self.usage["output_tokens"] += budget["max_output_tokens"]
            result = {"status": "failed", "reason": str(exc), "cached": False}
        self.cache[key] = copy.deepcopy(result)
        return result

    def evaluate(self, verification, step_context=None, deadline_monotonic=None):
        if deadline_monotonic is not None and (
            isinstance(deadline_monotonic, bool)
            or not isinstance(deadline_monotonic, (int, float))
            or not math.isfinite(deadline_monotonic)
        ):
            raise ValueError("Evidence deadline must be a finite monotonic timestamp")
        if not isinstance(verification, dict) or verification.get("status") != "catalog_and_rule_verified":
            raise ValueError("Independent evidence requires completed catalog/rule verification")
        steps = verification.get("steps")
        if not isinstance(steps, list) or not steps:
            raise ValueError("Verified route requires reaction steps")
        contexts = {} if step_context is None else step_context
        if not isinstance(contexts, dict) or set(contexts) - {step["node_id"] for step in steps}:
            raise ValueError("Step context must reference reaction node IDs in this route")
        reports, reasons = [], []
        for step in steps:
            node_id = step["node_id"]
            reactants = _smiles_list(step["reactants"], "verified reactants")
            product = canonical(step["product"])
            context = _context(
                contexts[node_id]
                if node_id in contexts
                else self.reaction_contexts.get(_json_key([sorted(reactants), product]))
            )
            context_complete = "reagents" in context and bool(context.get("conditions"))
            prediction = self._predict(reactants, context, deadline_monotonic)
            if prediction["status"] == "predicted":
                matches = matching_prediction_ranks(prediction["predictions"], product, self.product_matching)
                prediction["status"] = "matched" if matches else "mismatched"
                prediction["matched_rank"] = min(matches) if matches else None
                prediction["matching_policy"] = self.product_matching
                prediction["interpretation"] = (
                    "Exact full-product agreement"
                    if self.product_matching == "exact"
                    else "Exact expected component present; other predicted components remain explicit"
                )
            prediction["input_coverage"] = {
                "reagents_supplied": "reagents" in context,
                "conditions_supplied": bool(context.get("conditions")),
                "conditions_used_by_predictor": bool(
                    self.predictor is not None and self.predictor.provenance.get("uses_conditions", False)
                ),
                "incomplete_context": not context_complete,
            }
            atom_sources = _atom_source_evidence(step, context)
            observations = self.observations.get(_reaction_key(reactants, product, context), [])
            outcomes = {record["outcome"] for record in observations}
            observed_status = (
                "conflicting" if len(outcomes) > 1 else next(iter(outcomes)) if outcomes else "unknown"
            )
            requirements = evidence_requirements(
                prediction["status"], context_complete, atom_sources["status"], observed_status
            )
            rejected = [
                name for name, passed in requirements.items() if getattr(self.policy, name) and not passed
            ]
            reasons.extend(f"{node_id}: {name} not satisfied" for name in rejected)
            reports.append(
                {
                    "node_id": node_id,
                    "reactants": reactants,
                    "product": product,
                    "rule_replay": "matched",
                    "forward_prediction": prediction,
                    "context": {
                        "status": "supplied_unvalidated" if context_complete else "incomplete",
                        "validated": False,
                        "supplied": context,
                    },
                    "atom_sources": atom_sources,
                    "observations": {
                        "status": observed_status,
                        "records": copy.deepcopy(observations),
                        "source_audited": False,
                        "interpretation": "Exact reaction/context matches to supplied observation records",
                    },
                    "eligibility": {"passed": not rejected, "unmet_requirements": rejected},
                }
            )
        return {
            "schema_version": 1,
            "profile": "legacy-conditional-v1",
            "evidence_level": "catalog_rule_replay_with_separately_reported_external_evidence",
            "product": canonical(verification["product"]),
            "steps": reports,
            "eligibility": {"passed": not reasons, "reasons": reasons},
            "policy": asdict(self.policy),
            "predictor_provenance": copy.deepcopy(self.provenance["predictor"]),
            "usage": self.usage.copy(),
            "interpretation": self.provenance["interpretation"],
        }


def evaluate_synthesis_evidence(
    verification, predictor=None, step_context=None, policy=None, observations=None
):
    """Convenience API; reuse SynthesisEvidenceEvaluator to share cache/budgets."""
    return SynthesisEvidenceEvaluator(predictor, policy, observations=observations).evaluate(
        verification, step_context
    )


def validate_evidence_config(config):
    """Validate an inline workflow configuration without loading model weights."""
    if not isinstance(config, dict):
        raise ValueError("synthesis_evidence must be an object")
    allowed = {
        "profile",
        "predictor",
        "policy",
        "max_predictor_calls",
        "max_input_tokens",
        "max_output_tokens",
        "product_matching",
        "observations",
        "reaction_contexts",
    }
    if set(config) - allowed:
        raise ValueError(f"Unknown synthesis_evidence keys: {sorted(set(config) - allowed)}")
    values = copy.deepcopy(config)
    profile = values.pop("profile", "legacy-conditional-v1")
    if profile not in {"legacy-conditional-v1", "computational-route-support-v1"}:
        raise ValueError("Unknown synthesis evidence profile")
    if profile == "computational-route-support-v1":
        from .computational_support import validate_config

        validate_config(config)
    provider = values.pop("predictor", None)
    if provider is not None:
        permitted = {
            "kind",
            "model_path",
            "model_revision",
            "device",
            "top_k",
            "max_new_tokens",
            "max_input_tokens",
        }
        if (
            not isinstance(provider, dict)
            or set(provider) - permitted
            or provider.get("kind") != "reactiont5"
        ):
            raise ValueError("Evidence predictor requires kind=reactiont5 and supported local-model settings")
        if not isinstance(provider.get("model_path"), str) or not provider["model_path"].strip():
            raise ValueError("Evidence predictor requires a local model_path")
        for key in ("top_k", "max_new_tokens", "max_input_tokens"):
            if key in provider:
                positive(provider[key], key)
        if "device" in provider and (
            not isinstance(provider["device"], str) or not provider["device"].strip()
        ):
            raise ValueError("Evidence predictor device must be nonempty")
        if provider.get("model_revision") is not None and (
            not isinstance(provider["model_revision"], str) or not provider["model_revision"].strip()
        ):
            raise ValueError("Evidence model_revision must be a nonempty declared revision")
    try:
        SynthesisEvidenceEvaluator(**values)
    except TypeError as exc:
        raise ValueError(f"Invalid synthesis evidence settings: {exc}") from exc
    return copy.deepcopy(config)


def evaluator_from_config(config, *, load_predictor=True):
    """Build the optional predictor once and share its budgets/cache across routes."""
    values = validate_evidence_config(config)
    if values.get("profile") == "computational-route-support-v1":
        from .computational_support import from_config

        return from_config(values, load_predictor=load_predictor)
    values.pop("profile", None)
    provider = values.pop("predictor", None)
    predictor = (
        ReactionT5Forward(**{key: value for key, value in provider.items() if key != "kind"})
        if provider is not None and load_predictor
        else None
    )
    return SynthesisEvidenceEvaluator(predictor=predictor, **values)


def _require_reactiont5_version():
    try:
        version = importlib.metadata.version("transformers")
    except importlib.metadata.PackageNotFoundError as exc:
        raise RuntimeError(
            "Install chembounce[reactiont5] for the compatible forward-prediction runtime"
        ) from exc
    try:
        major, minor = (int(part) for part in version.split(".")[:2])
    except ValueError as exc:
        raise RuntimeError(f"Unrecognized Transformers version for ReactionT5: {version}") from exc
    if not (4, 51) <= (major, minor) < (5, 0):
        raise RuntimeError(
            f"ReactionT5 requires transformers>=4.51,<5; found {version}. "
            "Transformers 5 changed this checkpoint's tied-weight loading and produced incompatible predictions. "
            "Install chembounce[reactiont5] in a compatible environment; SynLlama's runtime is separate."
        )
    return version


class ReactionT5Forward:
    """Local, deterministic ReactionT5v2 forward inference using documented inputs.

    The model consumes reactants and reagents, not temperature/time conditions.
    Those conditions are retained in evidence but never silently encoded using
    an invented model prompt. Revision attribution is declared separately from
    independently recorded local file hashes.
    """

    def __init__(
        self, model_path, model_revision=None, device="cpu", top_k=5, max_new_tokens=128, max_input_tokens=512
    ):
        for name, value in (
            ("top_k", top_k),
            ("max_new_tokens", max_new_tokens),
            ("max_input_tokens", max_input_tokens),
        ):
            positive(value, name)
        path = Path(model_path).resolve()
        if not path.is_dir() or not (path / "config.json").is_file() or not list(path.glob("*.safetensors")):
            raise ValueError("ReactionT5 requires local config, tokenizer and safetensors weights")
        if model_revision is not None and (not isinstance(model_revision, str) or not model_revision.strip()):
            raise ValueError("model_revision must be a nonempty declared revision")
        transformers_version = _require_reactiont5_version()
        try:
            import torch
            from transformers import AutoModelForSeq2SeqLM, AutoTokenizer
        except ImportError as exc:
            raise RuntimeError(
                "Install the reactiont5 optional dependency group for compatible ReactionT5 inference"
            ) from exc
        self.torch, self.device = torch, device
        self.top_k, self.max_new_tokens, self.max_input_tokens = top_k, max_new_tokens, max_input_tokens
        self.tokenizer = AutoTokenizer.from_pretrained(
            str(path), local_files_only=True, trust_remote_code=False
        )
        self.model = (
            AutoModelForSeq2SeqLM.from_pretrained(
                str(path), local_files_only=True, trust_remote_code=False, use_safetensors=True
            )
            .to(device)
            .eval()
        )
        self.provenance = {
            "provider": "reactiont5v2_forward_local",
            "model_path": str(path),
            "declared_revision": model_revision,
            "model_files_sha256": {p.name: digest(p) for p in sorted(path.iterdir()) if p.is_file()},
            "device": device,
            "top_k": top_k,
            "max_new_tokens": max_new_tokens,
            "max_input_tokens": max_input_tokens,
            "do_sample": False,
            "uses_conditions": False,
            "source": "https://huggingface.co/sagawa/ReactionT5v2-forward",
            "publication": "https://doi.org/10.1186/s13321-025-01075-4",
            "declared_license": "MIT (upstream model card)",
            "input_format": "REACTANT:<dot-separated SMILES>REAGENT:<dot-separated SMILES>",
            "torch_version": importlib.metadata.version("torch"),
            "transformers_version": transformers_version,
            "runtime_compatibility": "transformers>=4.51,<5; v5 checkpoint loading is unsupported",
        }

    def _inputs(self, reactants, reagents, conditions):
        # Sort multisets for cache/order invariance; never deduplicate stoichiometry.
        reactants = sorted(_smiles_list(list(reactants), "reactants"))
        if not reactants:
            raise ValueError("ReactionT5 requires at least one reactant")
        reagents = sorted(_smiles_list(list(reagents), "reagents", reagents=True))
        _json_key(conditions)
        prompt = "REACTANT:" + ".".join(reactants) + "REAGENT:" + ".".join(reagents)
        tokens = self.tokenizer(prompt, return_tensors="pt", truncation=False)
        size = int(tokens["input_ids"].shape[1])
        if size > self.max_input_tokens:
            raise ValueError(f"ReactionT5 input has {size} tokens, exceeding limit {self.max_input_tokens}")
        return tokens, size

    def token_budget(self, reactants, reagents, conditions):
        _, size = self._inputs(reactants, reagents, conditions)
        return {"input_tokens": size, "max_output_tokens": self.top_k * self.max_new_tokens}

    def predict(self, reactants, reagents, conditions):
        inputs, input_tokens = self._inputs(reactants, reagents, conditions)
        with self.torch.inference_mode():
            output = self.model.generate(
                **inputs.to(self.device),
                do_sample=False,
                num_beams=self.top_k,
                num_return_sequences=self.top_k,
                max_new_tokens=self.max_new_tokens,
                return_dict_in_generate=True,
                output_scores=False,
            )
        products, completed, output_tokens = [], [], 0
        eos, pad = self.tokenizer.eos_token_id, self.tokenizer.pad_token_id
        for sequence in output.sequences.tolist():
            generated = sequence[1:]  # T5 decoder start token is not generated.
            if eos is not None and eos in generated:
                generated = generated[: generated.index(eos) + 1]
                complete = True
            else:
                while generated and generated[-1] == pad:
                    generated.pop()
                complete = False
            output_tokens += len(generated)
            products.append(
                self.tokenizer.decode(sequence, skip_special_tokens=True).replace(" ", "").rstrip(".")
            )
            completed.append(complete)
        return {
            "products": products,
            "complete": completed,
            "input_tokens": input_tokens,
            "output_tokens": output_tokens,
        }
