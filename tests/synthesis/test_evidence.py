import copy

import pytest

from chembounce_synthesis.common import canonical
from chembounce_synthesis.evidence import (
    EvidencePolicy,
    ReactionT5Forward,
    SynthesisEvidenceEvaluator,
    evaluator_from_config,
    evaluate_synthesis_evidence,
    validate_evidence_config,
)
from chembounce_synthesis.routes import verify_route


class Predictor:
    provenance = {"provider": "independent_test_predictor", "uses_conditions": False}

    def __init__(self, products, *, complete=None, fail=False):
        self.products = products
        self.complete = [True] * len(products) if complete is None else complete
        self.fail, self.calls, self.inputs = fail, 0, []

    def token_budget(self, reactants, reagents, conditions):
        return {"input_tokens": 5, "max_output_tokens": 16}

    def predict(self, reactants, reagents, conditions):
        self.calls += 1
        self.inputs.append((list(reactants), list(reagents), copy.deepcopy(conditions)))
        if self.fail:
            raise RuntimeError("simulated independent prediction failure")
        return {"products": self.products, "complete": self.complete, "input_tokens": 5, "output_tokens": 8}


def verified(assets, route):
    return verify_route(route, *assets)


def test_external_prediction_does_not_receive_product_or_template(assets, route):
    replay = verified(assets, route)
    predictor = Predictor([replay["product"]])
    engine = SynthesisEvidenceEvaluator(predictor, {"require_forward_agreement": True})
    result = engine.evaluate(replay)
    assert result["eligibility"]["passed"]
    step = result["steps"][0]
    assert step["forward_prediction"]["status"] == "matched"
    assert step["forward_prediction"]["matched_rank"] == 1
    assert step["forward_prediction"]["input_coverage"]["incomplete_context"]
    assert step["context"]["status"] == "incomplete"
    assert predictor.inputs == [(replay["steps"][0]["reactants"], [], {})]
    assert result["predictor_provenance"] == predictor.provenance
    assert "not synthesis probability" in result["interpretation"]


def test_product_mismatch_reuses_input_only_cache_without_false_match(assets, route):
    replay = verified(assets, route)
    predictor = Predictor([replay["product"]])
    engine = SynthesisEvidenceEvaluator(predictor, {"require_forward_agreement": True})
    assert engine.evaluate(replay)["eligibility"]["passed"]
    mismatched = copy.deepcopy(replay)
    mismatched["product"] = mismatched["steps"][0]["product"] = "CCC(=O)Nc1ccncc1"
    result = engine.evaluate(mismatched)
    assert not result["eligibility"]["passed"]
    assert result["steps"][0]["forward_prediction"]["status"] == "mismatched"
    assert result["steps"][0]["forward_prediction"]["cached"]
    assert predictor.calls == 1
    assert engine.usage == {"predictor_calls": 1, "input_tokens": 5, "output_tokens": 8, "cache_hits": 1}


def test_unknown_external_evidence_cannot_satisfy_a_requested_requirement(assets, route):
    result = evaluate_synthesis_evidence(
        verified(assets, route),
        policy={
            "require_forward_agreement": True,
            "require_complete_context": True,
            "require_observed_success": True,
        },
    )
    assert not result["eligibility"]["passed"]
    assert len(result["eligibility"]["reasons"]) == 3
    assert result["steps"][0]["forward_prediction"]["status"] == "unknown"
    assert result["steps"][0]["observations"]["status"] == "unknown"


def test_supplied_context_is_preserved_and_never_claimed_validated(assets, route):
    replay = verified(assets, route)
    context = {
        "r0": {"reagents": ["CCN(CC)CC"], "conditions": {"temperature_C": 25}, "source": "test specification"}
    }
    predictor = Predictor([replay["product"]])
    result = evaluate_synthesis_evidence(
        replay, predictor, context, {"require_forward_agreement": True, "require_complete_context": True}
    )
    assert result["eligibility"]["passed"]
    step = result["steps"][0]
    assert step["context"]["status"] == "supplied_unvalidated"
    assert step["context"]["validated"] is False
    assert step["context"]["supplied"] == context["r0"]
    assert not step["forward_prediction"]["input_coverage"]["conditions_used_by_predictor"]
    assert predictor.inputs[0][1:] == (["CCN(CC)CC"], {"temperature_C": 25})


@pytest.mark.parametrize(
    "budget,reason",
    [
        ({"max_predictor_calls": 0}, "call budget"),
        ({"max_input_tokens": 4}, "input-token"),
        ({"max_output_tokens": 15}, "output-token"),
    ],
)
def test_inference_budgets_fail_closed_before_call(assets, route, budget, reason):
    replay = verified(assets, route)
    predictor = Predictor([replay["product"]])
    engine = SynthesisEvidenceEvaluator(predictor, {"require_forward_agreement": True}, **budget)
    result = engine.evaluate(replay)
    assert not result["eligibility"]["passed"]
    assert reason in result["steps"][0]["forward_prediction"]["reason"]
    assert predictor.calls == 0


def test_failures_charge_reserved_tokens_and_are_cached(assets, route):
    predictor = Predictor([], fail=True)
    engine = SynthesisEvidenceEvaluator(predictor, {"require_forward_agreement": True})
    replay = verified(assets, route)
    assert not engine.evaluate(replay)["eligibility"]["passed"]
    result = engine.evaluate(replay)
    assert result["steps"][0]["forward_prediction"]["status"] == "failed"
    assert predictor.calls == 1
    assert engine.usage["output_tokens"] == 16


def test_context_and_reactant_multiplicity_are_part_of_cache_key(assets, route):
    replay = verified(assets, route)
    predictor = Predictor([replay["product"]])
    engine = SynthesisEvidenceEvaluator(predictor)
    engine.evaluate(replay, {"r0": {"reagents": [], "conditions": {"temperature_C": 20}}})
    engine.evaluate(replay, {"r0": {"reagents": [], "conditions": {"temperature_C": 30}}})
    changed = copy.deepcopy(replay)
    changed["steps"][0]["reactants"].append(changed["steps"][0]["reactants"][0])
    engine.evaluate(changed, {"r0": {"reagents": [], "conditions": {"temperature_C": 30}}})
    assert predictor.calls == 3


def test_valid_mixtures_are_preserved_with_explicit_matching_policy(assets, route):
    replay = verified(assets, route)
    mixture = replay["product"] + ".[Na+]"
    strict = SynthesisEvidenceEvaluator(Predictor([mixture]), {"require_forward_agreement": True})
    exact = strict.evaluate(replay)
    assert not exact["eligibility"]["passed"]
    prediction = exact["steps"][0]["forward_prediction"]["predictions"][0]
    assert "invalid_reason" not in prediction
    assert prediction["comparison_note"] == "unsupported_mixture_for_exact_matching"
    assert prediction["product"].count(".") == 1
    component = SynthesisEvidenceEvaluator(
        Predictor([mixture]), {"require_forward_agreement": True}, product_matching="component"
    ).evaluate(replay)
    assert component["eligibility"]["passed"]
    prediction = component["steps"][0]["forward_prediction"]
    assert prediction["matching_policy"] == "component"
    assert "[Na+]" in prediction["predictions"][0]["components"]


def test_component_matching_ignores_atom_map_labels_without_removing_coproducts(assets, route):
    replay = verified(assets, route)
    mixture = replay["product"].replace("CC(=O)", "[CH3:1]C(=O)") + ".[Na+]"
    result = SynthesisEvidenceEvaluator(
        Predictor([mixture]), {"require_forward_agreement": True}, product_matching="component"
    ).evaluate(replay)
    assert result["eligibility"]["passed"]
    assert "[Na+]" in result["steps"][0]["forward_prediction"]["predictions"][0]["components"]


def test_invalid_or_truncated_outputs_cannot_match(assets, route):
    replay = verified(assets, route)
    for predictor in (Predictor([replay["product"]], complete=[False]), Predictor(["not_smiles.[Na+]"])):
        result = evaluate_synthesis_evidence(replay, predictor, policy={"require_forward_agreement": True})
        assert not result["eligibility"]["passed"]
        assert result["steps"][0]["forward_prediction"]["status"] == "unknown"
        assert result["steps"][0]["forward_prediction"]["predictions"][0]["invalid_reason"]


def test_implicit_atom_sources_need_mapped_source_not_element_presence(assets, route):
    replay = verified(assets, route)
    replay["steps"][0]["unmapped_product_template_atoms"] = [
        {"atom_index": 3, "element": "C", "map_number": 0}
    ]
    engine = SynthesisEvidenceEvaluator(policy={"require_explicit_atom_sources": True})
    # Methanol includes carbon but element presence alone is not provenance.
    context = {"r0": {"reagents": ["OC"], "conditions": {"temperature_C": 20}}}
    assert not engine.evaluate(replay, context)["eligibility"]["passed"]
    context["r0"]["atom_sources"] = [
        {
            "template_atom_index": 3,
            "reagent_index": 0,
            "reagent_atom_index": 1,
            "source": "curated fixture",
            "record_id": "atom-map-1",
        }
    ]
    result = engine.evaluate(replay, context)
    assert result["eligibility"]["passed"]
    assert result["steps"][0]["atom_sources"]["status"] == "supplied_mapped_attribution"
    assert not result["steps"][0]["atom_sources"]["source_audited"]
    context["r0"]["atom_sources"][0]["reagent_atom_index"] = 0
    with pytest.raises(ValueError, match="element"):
        engine.evaluate(replay, context)


def observation(replay):
    return {
        "kind": "literature",
        "source": "https://example.org/test-fixture-not-an-experiment",
        "record_id": "fixture-1",
        "reactants": list(reversed(replay["steps"][0]["reactants"])),
        "product": replay["product"],
        "reagents": ["CCN(CC)CC"],
        "conditions": {"temperature_C": 20},
        "outcome": "success",
    }


def test_observation_requires_exact_multiset_product_and_context(assets, route):
    replay = verified(assets, route)
    record = observation(replay)
    engine = SynthesisEvidenceEvaluator(policy={"require_observed_success": True}, observations=[record])
    assert not engine.evaluate(replay)["eligibility"]["passed"]
    context = {"r0": {"reagents": record["reagents"], "conditions": record["conditions"]}}
    result = engine.evaluate(replay, context)
    assert result["eligibility"]["passed"]
    assert result["steps"][0]["observations"]["records"][0]["record_id"] == "fixture-1"
    assert not result["steps"][0]["observations"]["source_audited"]
    for alteration in ("reactants", "product", "conditions"):
        changed, changed_context = copy.deepcopy(replay), copy.deepcopy(context)
        if alteration == "reactants":
            changed["steps"][0]["reactants"].append(changed["steps"][0]["reactants"][0])
        elif alteration == "product":
            changed["steps"][0]["product"] = canonical("CCC(=O)Nc1ccncc1")
        else:
            changed_context["r0"]["conditions"]["temperature_C"] = 30
        assert not engine.evaluate(changed, changed_context)["eligibility"]["passed"]


def test_conflicting_observations_do_not_satisfy_success_gate(assets, route):
    replay = verified(assets, route)
    success = observation(replay)
    failure = {**success, "record_id": "fixture-2", "outcome": "failure"}
    engine = SynthesisEvidenceEvaluator(
        policy={"require_observed_success": True}, observations=[success, failure]
    )
    result = engine.evaluate(
        replay, {"r0": {"reagents": success["reagents"], "conditions": success["conditions"]}}
    )
    assert not result["eligibility"]["passed"]
    assert result["steps"][0]["observations"]["status"] == "conflicting"


def test_context_validation_is_strict(assets, route):
    replay = verified(assets, route)
    engine = SynthesisEvidenceEvaluator()
    for context in (
        [],
        {"absent-node": {}},
        {"r0": {"reagents": "CCN"}},
        {"r0": {"conditions": {"temperature": float("nan")}}},
    ):
        with pytest.raises(ValueError):
            engine.evaluate(replay, context)
    with pytest.raises(ValueError, match="true or false"):
        EvidencePolicy(require_forward_agreement=1)
    with pytest.raises(ValueError, match="observations must be a list"):
        SynthesisEvidenceEvaluator(observations={})
    with pytest.raises(ValueError, match="local"):
        ReactionT5Forward("/nonexistent/reaction-t5")


def test_configured_context_matches_exact_reaction_and_not_reused_node_id(assets, route):
    replay = verified(assets, route)
    context = {
        "reactants": list(reversed(replay["steps"][0]["reactants"])),
        "product": replay["product"],
        "reagents": [],
        "conditions": {"temperature_C": 20},
    }
    engine = evaluator_from_config(
        {"policy": {"require_complete_context": True}, "reaction_contexts": [context]}
    )
    assert engine.evaluate(replay)["eligibility"]["passed"]
    expanded = copy.deepcopy(replay)
    expanded["steps"][0]["reactants"][0] = "CCC(=O)O"
    assert not engine.evaluate(expanded)["eligibility"]["passed"]
    with pytest.raises(ValueError, match="Multiple contexts"):
        evaluator_from_config({"reaction_contexts": [context, context]})


def test_configuration_validation_does_not_load_model(monkeypatch):
    def forbidden(**kwargs):
        raise AssertionError("Configuration validation must not load model weights")

    monkeypatch.setattr("chembounce_synthesis.evidence.ReactionT5Forward", forbidden)
    config = {
        "predictor": {"kind": "reactiont5", "model_path": "relative/model", "top_k": 3},
        "policy": {"require_forward_agreement": True},
        "max_predictor_calls": 10,
    }
    assert validate_evidence_config(config) == config
    for invalid in (
        {"invented": True},
        {"predictor": {"kind": "reactiont5", "model_path": "model", "top_k": True}},
        {"policy": {"invented": True}},
        {"max_output_tokens": -1},
        {"product_matching": "strip_salts_silently"},
    ):
        with pytest.raises(ValueError):
            validate_evidence_config(invalid)


def test_reactiont5_inputs_follow_published_format_without_product_or_conditions():
    class Tokens:
        shape = (1, 9)

    class Tokenizer:
        def __init__(self):
            self.prompts = []

        def __call__(self, prompt, **kwargs):
            self.prompts.append((prompt, kwargs))
            return {"input_ids": Tokens()}

    adapter = ReactionT5Forward.__new__(ReactionT5Forward)
    adapter.tokenizer, adapter.max_input_tokens = Tokenizer(), 10
    adapter.top_k, adapter.max_new_tokens = 3, 20
    assert adapter.token_budget(["CN", "CC(=O)Cl", "CN"], ["[Cl-].[Na+]"], {"temperature_C": 20}) == {
        "input_tokens": 9,
        "max_output_tokens": 60,
    }
    prompt, kwargs = adapter.tokenizer.prompts[0]
    assert prompt == "REACTANT:CC(=O)Cl.CN.CNREAGENT:[Cl-].[Na+]"
    assert "20" not in prompt
    assert kwargs["truncation"] is False
    adapter.max_input_tokens = 8
    with pytest.raises(ValueError, match="exceeding limit"):
        adapter.token_budget(["CN"], [], {})


@pytest.mark.parametrize("version", ["4.50.3", "5.0.0", "6.1.0"])
def test_reactiont5_rejects_incompatible_transformers_runtime(monkeypatch, version):
    from chembounce_synthesis.evidence import _require_reactiont5_version

    monkeypatch.setattr("chembounce_synthesis.evidence.importlib.metadata.version", lambda name: version)
    with pytest.raises(RuntimeError, match="requires transformers>=4.51,<5"):
        _require_reactiont5_version()


def test_reactiont5_accepts_tested_transformers_runtime(monkeypatch):
    from chembounce_synthesis.evidence import _require_reactiont5_version

    monkeypatch.setattr("chembounce_synthesis.evidence.importlib.metadata.version", lambda name: "4.51.3")
    assert _require_reactiont5_version() == "4.51.3"


def test_expired_deadline_skips_inference_and_preserves_unknown_status(assets, route, monkeypatch):
    replay = verified(assets, route)
    predictor = Predictor([replay["product"]])
    engine = SynthesisEvidenceEvaluator(predictor, {"require_forward_agreement": True})
    monkeypatch.setattr("chembounce_synthesis.evidence.time.monotonic", lambda: 20.0)
    result = engine.evaluate(replay, deadline_monotonic=19.0)
    assert not result["eligibility"]["passed"]
    assert "deadline exhausted" in result["steps"][0]["forward_prediction"]["reason"]
    assert predictor.calls == 0
    # Expiration is not cached as chemical evidence and can be retried later.
    assert engine.evaluate(replay, deadline_monotonic=21.0)["eligibility"]["passed"]
