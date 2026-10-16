"""Synthetic runtime tensors exercise strict decoding; these tests run no model.

Known amide token vectors reproduce the observed EOS-filled tensor layout.
Ancestry below is deliberately synthetic, not reconstructed historical evidence.
"""

import copy
from contextlib import nullcontext
from types import SimpleNamespace

import numpy as np
import pytest

from chembounce_synthesis.computational_support import ComputationalRouteSupport
from chembounce_synthesis.forward_support import PinnedReactionT5Forward


def tensors():
    return [[0, 70, 70, 1, 1], [0, 70, 82, 70, 1], [0, 70, 70, 70, 1]], [
        [0, 1, 2, -1],
        [0, 1, 2, 2],
        [0, 0, 1, 2],
    ]


def test_realistic_known_amide_vectors_use_only_witnessed_generated_prefix(strict_predictor):
    tokens = {
        0: "<pad>",
        1: "</s>",
        2: "<unk>",
        3: "C",
        4: "c",
        5: ")",
        6: "(",
        7: "1",
        9: "O",
        12: "N",
        15: "=",
        19: "ccc",
        21: "F",
        26: "CC",
        45: "O",
        61: "ccnc",
    }
    strict_predictor.tokenizer.decode = lambda ids, **kwargs: "".join(tokens[i] for i in ids)
    generated = [
        [26, 6, 15, 9, 5, 12, 4, 7, 61, 4, 7, 1],
        [45, 15, 3, 6, 3, 12, 4, 7, 61, 4, 7, 5, 12, 4, 7, 61, 4, 7, 1],
        [45, 15, 3, 6, 3, 12, 4, 7, 61, 4, 7, 5, 4, 7, 19, 6, 21, 5, 4, 4, 7, 1],
    ]
    width = max(map(len, generated))
    sequences = [[0] + ids + [1] * (width - len(ids)) for ids in generated]
    ancestry = [[0] * len(ids) + [-1] * (width - len(ids)) for ids in generated]
    predictions, charged = strict_predictor.inspect_sequences(sequences, ancestry)
    assert charged == 53
    assert predictions[0]["complete"] and predictions[0]["product"] == "CC(=O)Nc1ccncc1"
    assert predictions[0]["original_decoded_text"].endswith("</s>" * 11)
    assert predictions[0]["padding_token_ids"] == [1] * 10
    assert predictions[0]["generated_token_ids"] == generated[0]
    assert predictions[0]["beam_indices"] == ancestry[0]
    assert all(p["product"] is not None for p in predictions)


@pytest.mark.parametrize(
    "fault",
    [
        "missing",
        "rows",
        "boolean",
        "range",
        "negative",
        "hole",
        "empty",
        "width",
        "ragged",
        "wrong_fill",
        "hidden_content",
    ],
)
def test_malformed_ancestry_or_unwitnessed_padding_fails_closed(strict_predictor, fault):
    sequences, ancestry = tensors()
    if fault == "missing":
        ancestry = None
    elif fault == "rows":
        ancestry.pop()
    elif fault == "boolean":
        ancestry[0][0] = True
    elif fault == "range":
        ancestry[0][0] = 3
    elif fault == "negative":
        ancestry[0][0] = -2
    elif fault == "hole":
        ancestry[0] = [0, -1, 2, -1]
    elif fault == "empty":
        ancestry[0] = [-1] * 4
    elif fault == "width":
        ancestry[0].pop()
    elif fault == "ragged":
        sequences[0].pop()
    elif fault == "wrong_fill":
        sequences[0][-1] = 0
    else:
        sequences[0][-1] = 70
    with pytest.raises(ValueError):
        strict_predictor.inspect_sequences(sequences, ancestry)


def test_eos_fill_never_completes_unfinished_prefix(strict_predictor):
    sequences, ancestry = tensors()
    ancestry[0] = [0, 1, -1, -1]
    predictions, charged = strict_predictor.inspect_sequences(sequences, ancestry)
    assert not predictions[0]["complete"] and predictions[0]["product"] is None
    assert predictions[0]["generated_token_ids"] == [70, 70]
    assert predictions[0]["padding_token_ids"] == [1, 1]
    assert charged == 10


@pytest.mark.parametrize("kind", ["extra_generated_eos", "internal_eos", "unknown", "generated_pad"])
def test_generated_special_content_is_never_erased_into_a_match(strict_predictor, kind):
    sequences, ancestry = tensors()
    if kind == "extra_generated_eos":
        ancestry[0] = [0, 1, 2, 2]
    elif kind == "internal_eos":
        sequences[0] = [0, 70, 1, 70, 1]
        ancestry[0] = [0, 1, 2, 2]
    elif kind == "unknown":
        sequences[0][1] = 2
    else:
        sequences[0][1] = 0
    predictions, charged = strict_predictor.inspect_sequences(sequences, ancestry)
    assert predictions[0]["product"] is None
    assert "invalid_reason" in predictions[0]
    assert charged == (12 if kind in {"extra_generated_eos", "internal_eos"} else 11)


def mock_generation(predictor, monkeypatch, sequences, ancestry):
    class Inputs(dict):
        def to(self, device):
            assert device == "cpu"
            return self

    prepared = predictor.prepare(["CC"], [])
    inputs = Inputs(input_ids=np.array([prepared["input_token_ids"]]))
    monkeypatch.setattr(predictor, "_inputs", lambda *args: (inputs, prepared["input_tokens"]))
    predictor.torch = SimpleNamespace(inference_mode=nullcontext)
    calls = []

    def generate(**kwargs):
        calls.append(kwargs)
        assert kwargs["output_scores"] is False and kwargs["return_dict_in_generate"] is True
        return SimpleNamespace(
            sequences=np.array(sequences), beam_indices=None if ancestry is None else np.array(ancestry)
        )

    predictor.model.generate = generate
    return prepared, calls


def test_provider_retains_actual_returned_ancestry_and_validates_again(strict_predictor, monkeypatch):
    sequences, ancestry = tensors()
    prepared, calls = mock_generation(strict_predictor, monkeypatch, sequences, ancestry)
    raw = PinnedReactionT5Forward.predict_raw(strict_predictor, prepared)
    assert len(calls) == 1
    assert raw["sequences"] == sequences and raw["beam_indices"] == ancestry
    assert raw["output_tokens"] == 11
    assert strict_predictor.validate_raw(raw, prepared)[0]["product"] == "CC"
    forged = copy.deepcopy(raw)
    forged["beam_indices"][0] = [0, 1, -1, -1]
    with pytest.raises(ValueError, match="evidence mismatch"):
        strict_predictor.validate_raw(forged, prepared)


def test_missing_ancestry_preserved_and_failed_call_charges_full_reservation(strict_predictor, monkeypatch):
    sequences, _ = tensors()
    _, calls = mock_generation(strict_predictor, monkeypatch, sequences, None)
    monkeypatch.setattr(
        strict_predictor,
        "predict_raw",
        lambda prepared: PinnedReactionT5Forward.predict_raw(strict_predictor, prepared),
    )
    evaluator = ComputationalRouteSupport(strict_predictor, max_output_tokens=768)
    prediction = evaluator._raw_prediction(["CC"], {}, None)
    assert prediction["status"] == "failed" and "beam_indices" in prediction["reason"]
    assert prediction["raw"]["sequences"] == sequences
    assert prediction["raw"]["beam_indices"] is None
    assert evaluator.usage["output_tokens"] == 768 and len(calls) == 1
    assert evaluator._raw_prediction(["CO"], {}, None)["status"] == "unknown"
    assert len(calls) == 1
