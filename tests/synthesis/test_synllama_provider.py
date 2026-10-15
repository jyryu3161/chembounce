"""Optional local-runtime contract tests; these tiny weights are not SynLlama."""

import json
import time

import pytest

TORCH = pytest.importorskip("torch", reason="Optional local model runtime")
pytest.importorskip("transformers", reason="Optional local model runtime")


@pytest.mark.parametrize("weight_format", ["safetensors", "pytorch_bin"])
def test_local_inference_formats_capture_raw_outputs_and_restore_rng(tmp_path, weight_format):
    from tokenizers import Tokenizer
    from tokenizers.models import WordLevel
    from tokenizers.pre_tokenizers import Whitespace
    from transformers import LlamaConfig, LlamaForCausalLM, PreTrainedTokenizerFast

    from chembounce_synthesis.synllama import SynLlama

    model_dir = tmp_path / "tiny-fixture"
    model_dir.mkdir()
    tokenizer = Tokenizer(
        WordLevel(
            {"[PAD]": 0, "[EOS]": 1, "[UNK]": 2, "fixture": 3, "CCO": 4, "C": 5, "N": 6, "O": 7},
            unk_token="[UNK]",
        )
    )
    tokenizer.pre_tokenizer = Whitespace()
    fast = PreTrainedTokenizerFast(
        tokenizer_object=tokenizer, pad_token="[PAD]", eos_token="[EOS]", unk_token="[UNK]"
    )
    fast.save_pretrained(model_dir)
    config = LlamaConfig(
        vocab_size=8,
        hidden_size=16,
        intermediate_size=32,
        num_hidden_layers=1,
        num_attention_heads=2,
        num_key_value_heads=1,
        max_position_embeddings=256,
        bos_token_id=None,
        eos_token_id=1,
        pad_token_id=0,
    )
    model = LlamaForCausalLM(config)
    if weight_format == "safetensors":
        model.save_pretrained(model_dir)
    else:
        config.save_pretrained(model_dir)
        TORCH.save(model.state_dict(), model_dir / "pytorch_model.bin")
    prompt = tmp_path / "prompt.json"
    prompt.write_text(json.dumps({"instruction": "fixture", "input": "SMILES_STRING"}))
    provider = SynLlama(model_dir, prompt, samples=1, max_new_tokens=4, seed=17)
    before = TORCH.random.get_rng_state().clone()
    rows = provider.propose("CCO")
    after = TORCH.random.get_rng_state()
    assert TORCH.equal(before, after)
    assert len(rows) == 1 and isinstance(rows[0]["response"], str)
    assert rows[0]["sampling"]["generated_tokens"] <= 4
    assert rows[0]["sampling"]["seed"] == 17
    assert provider.provenance["provider"] == "synllama_local"
    assert provider.provenance["use_cache"] is True
    assert provider.provenance["model_files_sha256"]
    assert provider.propose("CCO", deadline_monotonic=time.monotonic() - 1) == []
    timed = provider.propose("CCO", deadline_monotonic=time.monotonic() + 10)
    assert len(timed) == 1 and timed[0]["sampling"]["elapsed_seconds"] >= 0
    assert TORCH.equal(before, TORCH.random.get_rng_state())
