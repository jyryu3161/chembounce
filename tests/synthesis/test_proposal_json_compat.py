"""Synthetic controls for the optional, audited proposal JSON dialect."""

import copy
import json
import time
from contextlib import nullcontext
from pathlib import Path
from types import SimpleNamespace

import pytest

from chembounce_synthesis.proposal_json import POLICY, bind_row, decode_response, verify_saved_row
from chembounce_synthesis import proposal_json
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


def _provider():
    return {
        "response_json_dialect": POLICY,
        "prompt_sha256": "a" * 64,
        "response_decoder_source_sha256": digest(Path(proposal_json.__file__)),
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
        deadline_monotonic=time.monotonic() - 1,
    )
    assert row["response"] == row["raw_response"]
    assert row["compatibility"]["status"] == "not_evaluated_deadline"
    assert verify_saved_row(row, _provider(), design="CCO") == row
    with pytest.raises(ValueError, match="differs from raw rederivation"):
        verify_saved_row(row, _provider(), design="CCC")


def test_saved_replay_rederives_raw_and_rejects_tampering(tmp_path):
    row = bind_row(
        _raw_invalid(),
        {"seed": 1, "sample_index": 0},
        "a" * 64,
        _provider()["response_decoder_source_sha256"],
        design="CCO",
    )
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
    wrong_source = json.loads(archive.read_text())
    wrong_source["provider"]["response_decoder_source_sha256"] = "0" * 64
    archive.write_text(json.dumps(wrong_source))
    with pytest.raises(ValueError, match="decoder source"):
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
    assert SavedProposals(archive).propose("CCO")[0] == {"response": _proposal(), "sampling": {}}


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


def test_provider_acquisition_uses_template_policy_without_model_runtime():
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
    provider.provenance = {"prompt_sha256": "a" * 64, **_provider()}
    strict = provider.propose("CCO")[0]
    assert strict["response"] == _raw_invalid() and "compatibility" not in strict
    provider.response_json_dialect = POLICY
    compatible = provider.propose("CCO")[0]
    assert compatible["raw_response"] == _raw_invalid()
    assert compatible["response"] == _proposal()
    assert compatible["sampling"]["sample_index"] == 0
    assert verify_saved_row(compatible, _provider(), design="CCO") == compatible
