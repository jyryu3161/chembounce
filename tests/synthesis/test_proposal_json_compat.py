"""Synthetic controls for the optional, audited proposal JSON dialect."""

import copy
import json
import sys
import time
from contextlib import nullcontext
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

from chembounce_synthesis.proposal_json import (
    POLICY,
    bind_row,
    decode_response,
    finalize_timing,
    provider_identity_sha256,
    verify_saved_row,
)
from chembounce_synthesis import proposal_json
from chembounce_synthesis import synllama
from chembounce_synthesis.common import digest
from chembounce_synthesis.routes import parse_response
from chembounce_synthesis.synllama import SavedProposals, SynLlama, load_prompt


def _proposal(smiles="F/C=C\\Cl"):
    return {
        "reactions": [{"reaction_template": "fixture", "reactants": ["C/C=C\\N"], "product": smiles}],
        "building_blocks": ["C/C=C\\N"],
    }


def _raw_invalid(field="product", follow="C"):
    raw = json.dumps(_proposal())
    if field == "product":
        return raw.replace("C\\\\Cl", "C\\" + follow + "l", 1)
    if field == "reactants":
        return raw.replace("C\\\\N", "C\\N", 1)
    if field == "building_blocks":
        head, tail = raw.rsplit("C\\\\N", 1)
        return head + "C\\N" + tail
    raise AssertionError(field)


def _provider(seed=1):
    return {
        "provider": "synllama_local",
        "response_json_dialect": POLICY,
        "prompt_sha256": "a" * 64,
        "response_decoder_source_sha256": digest(Path(proposal_json.__file__)),
        "model_files_sha256": {"fixture": "b" * 64},
        "samples": 1,
        "seed": seed,
    }


def test_strict_json_retains_existing_value_and_raw_string():
    raw = json.dumps(_proposal())
    row = decode_response(raw)
    assert row["response"] == raw
    assert row["compatibility"]["status"] == "strict_valid"
    assert parse_response(row["response"]) == _proposal()
    assert "repairs" not in row["compatibility"]


@pytest.mark.parametrize("field", ["product", "reactants", "building_blocks"])
def test_allowed_literal_bond_backslash_is_exact_and_audited(field):
    raw = _raw_invalid(field)
    row = decode_response(raw)
    assert row["raw_response"] == raw
    assert row["compatibility"]["status"] == "compatibility_decoded"
    assert row["compatibility"]["repairs"][0]["path"] == (
        ["reactions", 0, field]
        if field == "product"
        else ["reactions", 0, "reactants", 0]
        if field == "reactants"
        else ["building_blocks", 0]
    )
    assert row["response"] == _proposal()
    assert json.loads(row["compatibility"]["derived_json_text"]) == _proposal()


@pytest.mark.parametrize("ambiguous_follow", ["n", "b"])
def test_rejects_unknown_escape_and_ambiguous_valid_escape_in_same_response(ambiguous_follow):
    assert decode_response(_raw_invalid(follow="q"))["compatibility"]["status"] == "rejected"
    ambiguous = _raw_invalid().replace("C/C=C\\\\N", "C/C=C\\" + ambiguous_follow, 1)
    assert decode_response(ambiguous)["compatibility"]["status"] == "rejected"
    assert (
        decode_response(_raw_invalid().replace("C\\Cl", "C\\\\\\Cl"))["compatibility"]["status"] == "rejected"
    )


def test_correct_unicode_escape_keeps_json_meaning_on_fallback():
    raw = _raw_invalid().replace("C/C=C\\\\N", "C/C=\\u0043\\\\N", 1)
    row = decode_response(raw)
    assert row["compatibility"]["status"] == "compatibility_decoded"
    assert row["response"]["reactions"][0]["reactants"][0] == "C/C=C\\N"


@pytest.mark.parametrize(
    "raw",
    [
        '{"metadata":{"product":"F/C=C\\Cl"},"reactions":[],"building_blocks":[]}',
        '{"reactions":[{"reaction_template":"F/C=C\\Cl","reactants":[],"product":"CC"}],"building_blocks":[]}',
        '{"reactions":[{"product":"F/C=C\\Cl","reactants":[]}],"building_blocks":[],"extra":NaN}',
        '{"reactions":[{"product":"F/C=C\\Cl","reactants":[]}],"building_blocks":[],"extra":"\\ud800"}',
        '{"reactions":[{"product":"F/C=C\\Cl","reactants":[]}],"building_blocks":[],"extra":"unterminated}',
        '{"reactions":[{"product":"F/C=C\\Cl","reactants":[]}],"building_blocks":[]} trailing',
        '{"reactions":[{"product":"F/C=C\\Cl","reactants":[]}],"building_blocks":[],"reactions":[]}',
        '{"reactions":[{"product":"F/C=C\\Cl","\\u0070roduct":"CC","reactants":[]}],"building_blocks":[]}',
    ],
)
def test_rejects_other_syntax_and_duplicate_decoded_keys(raw):
    assert decode_response(raw)["compatibility"]["status"] == "rejected"


def test_fenced_json_supported_but_trailing_text_rejected():
    raw = "```json\n" + _raw_invalid() + "\n```"
    assert decode_response(raw)["compatibility"]["status"] == "compatibility_decoded"
    assert decode_response(raw + " trailing")["compatibility"]["status"] == "rejected"


def test_repair_count_and_deadline_bound():
    raw = _raw_invalid().replace('"product": "F/C=C\\Cl"', '"product": "' + "C\\c" * 257 + '"')
    assert decode_response(raw)["compatibility"]["status"] == "rejected"
    row = bind_row(
        _raw_invalid(),
        {"seed": 1, "sample_index": 0},
        "a" * 64,
        _provider()["response_decoder_source_sha256"],
        design="CCO",
        provider_identity=provider_identity_sha256(_provider()),
        deadline_monotonic=time.monotonic() - 1,
    )
    assert row["response"] is None and row["raw_response"] == _raw_invalid()
    assert row["compatibility"]["status"] == "not_evaluated_deadline"
    finalize_timing(row, 1.0, True)
    assert verify_saved_row(row, _provider(), design="CCO") == row
    with pytest.raises(ValueError, match="differs from raw rederivation"):
        verify_saved_row(row, _provider(), design="CCC")


def test_opt_in_depth_and_length_bounds_precede_strict_json():
    for depth in (130, 1200):
        raw = '{"reactions":[],"building_blocks":[],"extra":' + "[" * depth + "0" + "]" * depth + "}"
        row = decode_response(raw)
        assert row["raw_response"] == raw
        assert row["compatibility"]["status"] == "rejected"
        assert "nesting limit" in row["compatibility"]["reason"]

    class Oversized(str):
        def __iter__(self):
            raise AssertionError("Oversized input was traversed")

        def encode(self, *args, **kwargs):
            raise AssertionError("Oversized input was hashed")

    with pytest.raises(ValueError, match="oversized"):
        decode_response(Oversized(" " * 250_001))


def test_rejected_and_deadline_responses_cannot_enter_route_parser(tmp_path):
    deep = '{"reactions":[],"building_blocks":[],"extra":' + "[" * 130 + "0" + "]" * 130 + "}"
    for row in (
        decode_response(deep),
        decode_response('{"reactions":[],"building_blocks":[]}', deadline_monotonic=0),
    ):
        assert row["raw_response"]
        assert row["response"] is None
        with pytest.raises(ValueError):
            parse_response(row["response"])
    rejected = bind_row(
        deep,
        {"seed": 1, "sample_index": 0},
        "a" * 64,
        _provider()["response_decoder_source_sha256"],
        design="CCO",
        provider_identity=provider_identity_sha256(_provider()),
    )
    finalize_timing(rejected, 1.0, False)
    archive = tmp_path / "deep-rejected.json"
    archive.write_text(
        json.dumps(
            {"schema_version": 1, "provider": _provider(), "proposals": [{"design": "CCO", **rejected}]}
        )
    )
    replayed = SavedProposals(archive).propose("CCO")[0]
    assert replayed["compatibility"]["status"] == "rejected"
    with pytest.raises(ValueError):
        parse_response(replayed["response"])


def test_oversized_generated_row_is_terminal_and_bounded_on_replay():
    raw = " " * 250_001
    row = bind_row(
        raw,
        {"seed": 1, "sample_index": 0},
        "a" * 64,
        _provider()["response_decoder_source_sha256"],
        design="CCO",
        provider_identity=provider_identity_sha256(_provider()),
    )
    finalize_timing(row, 1.0, False)
    assert row["response"] is None
    assert row["raw_response"] == raw
    assert row["compatibility"]["strict_parser_outcome"] == "not_parsed_input_limit"
    assert verify_saved_row(row, _provider(), design="CCO") == row
    long_row = bind_row(
        raw + "more",
        {"seed": 1, "sample_index": 0},
        "a" * 64,
        _provider()["response_decoder_source_sha256"],
        design="CCO",
        provider_identity=provider_identity_sha256(_provider()),
    )
    assert len(long_row["raw_response"]) == 250_001
    assert long_row["compatibility"]["raw_omitted_chars"] == 4
    finalize_timing(long_row, 1.0, False)
    assert verify_saved_row(long_row, _provider(), design="CCO") == long_row


def test_replay_rejects_oversized_deadline_and_missing_timing_seal():
    row = bind_row(
        _raw_invalid(),
        {"seed": 1, "sample_index": 0},
        "a" * 64,
        _provider()["response_decoder_source_sha256"],
        design="CCO",
        provider_identity=provider_identity_sha256(_provider()),
        deadline_monotonic=0,
    )
    finalize_timing(row, 1.0, True)
    oversized = copy.deepcopy(row)
    oversized["raw_response"] = " " * 250_001
    with pytest.raises(ValueError, match="oversized"):
        verify_saved_row(oversized, _provider(), design="CCO")
    stripped = copy.deepcopy(row)
    for key in ("elapsed_seconds", "deadline_reached"):
        stripped["sampling"].pop(key)
    for key in ("sample_elapsed_seconds", "sample_deadline_reached"):
        stripped["compatibility"].pop(key)
    with pytest.raises(ValueError, match="timing seal"):
        verify_saved_row(stripped, _provider(), design="CCO")


def test_saved_replay_rederives_raw_and_rejects_tampering(tmp_path):
    row = bind_row(
        _raw_invalid(),
        {"seed": 1, "sample_index": 0},
        "a" * 64,
        _provider()["response_decoder_source_sha256"],
        design="CCO",
        provider_identity=provider_identity_sha256(_provider()),
    )
    finalize_timing(row, 1.0, False)
    assert verify_saved_row(row, _provider(), design="CCO") == row
    archive = tmp_path / "proposals.json"
    archive.write_text(
        json.dumps({"schema_version": 1, "provider": _provider(), "proposals": [{"design": "CCO", **row}]})
    )
    saved = SavedProposals(archive)
    assert saved.propose("CCO")[0] == row
    archive.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "provider": saved.provenance,
                "proposals": [{"design": "CCO", **saved.propose("CCO")[0]}],
            }
        )
    )
    assert SavedProposals(archive).propose("CCO")[0] == row
    downgraded = copy.deepcopy(saved.provenance)
    for key in ("response_json_dialect", "prompt_sha256", "response_decoder_source_sha256"):
        downgraded.pop(key)
    stripped = {key: value for key, value in row.items() if key not in ("raw_response", "compatibility")}
    archive.write_text(
        json.dumps(
            {"schema_version": 1, "provider": downgraded, "proposals": [{"design": "CCO", **stripped}]}
        )
    )
    with pytest.raises(ValueError, match="stripped or changed"):
        SavedProposals(archive)
    archive.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "provider": _provider(),
                "proposals": [{"design": "CCO", **row}, {"design": "CCO", **row}],
            }
        )
    )
    with pytest.raises(ValueError, match="Duplicate saved compatibility sample index"):
        SavedProposals(archive)
    changed_provider = _provider()
    changed_provider["model_files_sha256"] = {"other": "c" * 64}
    archive.write_text(
        json.dumps(
            {"schema_version": 1, "provider": changed_provider, "proposals": [{"design": "CCO", **row}]}
        )
    )
    with pytest.raises(ValueError, match="raw rederivation"):
        SavedProposals(archive)
    deeply_wrapped = _provider()
    for _ in range(33):
        deeply_wrapped = {
            "provider": "saved_proposals",
            "response_json_dialect": POLICY,
            "prompt_sha256": "a" * 64,
            "response_decoder_source_sha256": _provider()["response_decoder_source_sha256"],
            "recorded_provider": deeply_wrapped,
        }
    archive.write_text(
        json.dumps({"schema_version": 1, "provider": deeply_wrapped, "proposals": [{"design": "CCO", **row}]})
    )
    with pytest.raises(ValueError, match="chain exceeds limit"):
        SavedProposals(archive)
    archive.write_text(
        json.dumps(
            {"schema_version": 1, "provider": saved.provenance, "proposals": [{"design": "CCO", **row}]}
        )
    )
    wrong_source = json.loads(archive.read_text())
    wrong_source["provider"]["response_decoder_source_sha256"] = "0" * 64
    archive.write_text(json.dumps(wrong_source))
    with pytest.raises(ValueError, match="source"):
        SavedProposals(archive)
    archive.write_text(
        json.dumps(
            {"schema_version": 1, "provider": saved.provenance, "proposals": [{"design": "CCO", **row}]}
        )
    )
    for changed in (
        {"raw_response": _raw_invalid().replace("Cl", "Nl")},
        {"response": {"reactions": [], "building_blocks": []}},
        {"compatibility": {**row["compatibility"], "repairs": []}},
        {"compatibility": None},
    ):
        tampered = copy.deepcopy(row)
        tampered.update(changed)
        with pytest.raises(ValueError):
            verify_saved_row(tampered, _provider(), design="CCO")
    stripped = json.loads(archive.read_text())
    del stripped["proposals"][0]["compatibility"]
    archive.write_text(json.dumps(stripped))
    with pytest.raises(ValueError, match="missing raw response or audit"):
        SavedProposals(archive)


def test_legacy_saved_proposals_keep_existing_dict_transport(tmp_path):
    archive = tmp_path / "legacy.json"
    archive.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": "CCO", "response": _proposal()}]})
    )
    saved = SavedProposals(archive)
    assert saved.propose("CCO")[0] == {"response": _proposal(), "sampling": {}}
    archive.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "provider": saved.provenance,
                "proposals": [{"design": "CCO", **saved.propose("CCO")[0]}],
            }
        )
    )
    assert SavedProposals(archive).propose("CCO")[0] == saved.propose("CCO")[0]


def test_prompt_opt_in_requires_exact_named_policy(tmp_path):
    prompt = tmp_path / "prompt.json"
    base = {"instruction": "fixture", "input": "SMILES_STRING"}
    prompt.write_text(json.dumps(base))
    assert "response_json_dialect" not in load_prompt(prompt)
    for bad in ("unknown", 1, None, True):
        prompt.write_text(json.dumps({**base, "response_json_dialect": bad}))
        with pytest.raises(ValueError, match="response_json_dialect"):
            load_prompt(prompt)
        with pytest.raises(ValueError, match="response_json_dialect"):
            SynLlama(tmp_path / "missing-model", prompt)
    prompt.write_text(json.dumps({**base, "response_json_dialect": POLICY}))
    assert load_prompt(prompt)["response_json_dialect"] == POLICY


def test_provider_acquisition_uses_template_policy_without_model_runtime(monkeypatch):
    class Batch(dict):
        def to(self, device):
            return self

        @property
        def input_ids(self):
            return SimpleNamespace(shape=(1, 1))

    class Tokenizer:
        eos_token_id, unk_token_id, pad_token_id = 1, 2, 0

        def __call__(self, *args, **kwargs):
            return Batch(input_ids=[99])

        def convert_tokens_to_ids(self, token):
            return 2

        def decode(self, tokens, **kwargs):
            return _raw_invalid()

    class Generated:
        def __getitem__(self, key):
            return [7]

    class Model:
        config = SimpleNamespace(max_position_embeddings=100)

        def generate(self, **kwargs):
            return Generated()

    fake_torch = SimpleNamespace(
        random=SimpleNamespace(
            fork_rng=lambda **kwargs: nullcontext(),
            default_generator=SimpleNamespace(manual_seed=lambda seed: None),
        ),
        inference_mode=lambda: nullcontext(),
    )
    provider = SynLlama.__new__(SynLlama)
    provider.template = {"instruction": "fixture", "input": "SMILES_STRING"}
    provider.tokenizer = Tokenizer()
    provider.model = Model()
    provider.torch = fake_torch
    provider.device = "cpu"
    provider.max_new_tokens = 4
    provider.samples = 1
    provider.seed = 17
    provider.response_json_dialect = "strict-json-v1"
    provider.provenance = _provider(seed=17)
    strict = provider.propose("CCO")[0]
    assert strict["response"] == _raw_invalid() and "compatibility" not in strict
    provider.response_json_dialect = POLICY
    compatible = provider.propose("CCO")[0]
    assert compatible["raw_response"] == _raw_invalid()
    assert compatible["response"] == _proposal()
    assert compatible["sampling"]["sample_index"] == 0
    assert verify_saved_row(compatible, _provider(seed=17), design="CCO") == compatible

    provider.samples = 2
    provider.provenance["samples"] = 2
    responses = iter((_raw_invalid(), " " * 250_001))
    provider.tokenizer.decode = lambda *args, **kwargs: next(responses)
    completed, oversized = provider.propose("CCO")
    assert completed["raw_response"] == _raw_invalid()
    assert oversized["compatibility"]["status"] == "rejected"
    assert oversized["raw_response"] == " " * 250_001
    assert oversized["response"] is None
    provider.samples = 1
    provider.provenance["samples"] = 1
    provider.tokenizer.decode = lambda *args, **kwargs: _raw_invalid()

    class Clock:
        now = 0.0

        def monotonic(self):
            return self.now

    clock = Clock()
    monkeypatch.setattr(synllama, "time", clock)
    monkeypatch.setattr(proposal_json, "time", clock)
    transforms = ModuleType("transformers")
    transforms.StoppingCriteria = type("StoppingCriteria", (), {})
    transforms.StoppingCriteriaList = list
    monkeypatch.setitem(sys.modules, "transformers", transforms)
    original_hash = proposal_json._sha

    def costly_hash(value):
        result = original_hash(value)
        clock.now += 2
        return result

    monkeypatch.setattr(proposal_json, "_sha", costly_hash)

    def generated_after_one_second(**kwargs):
        clock.now += 1
        return Generated()

    provider.model.generate = generated_after_one_second
    late = provider.propose("CCO", deadline_monotonic=5)[0]
    assert late["compatibility"]["status"] == "not_evaluated_deadline"
    assert late["sampling"]["deadline_reached"] is True
    assert late["sampling"]["elapsed_seconds"] >= clock.now
    assert late["compatibility"]["sample_elapsed_seconds"] == late["sampling"]["elapsed_seconds"]
    assert verify_saved_row(late, _provider(seed=17), design="CCO") == late
