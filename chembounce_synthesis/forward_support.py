"""Pinned local ReactionT5 adapter for conditional computational selection support.

Only molecular reactant/reagent multisets enter the model. The raw token stream
is retained; tokenization spaces alone may be removed before strict SMILES parsing.
"""

import copy
import hashlib
import importlib.metadata
import json
import os
import platform
from pathlib import Path

from rdkit import Chem

from .common import digest
from .evidence import ReactionT5Forward

GENERATION_SOURCE_SHA256 = "dc9044d10c00abe7d341f2e777fb623925ff59dda2c1d814ecc7d6a685ad2ee1"
PADDING_POLICY = (
    "transformers-4.51.3-beam-ancestry-v1; generated-prefix indices 0..2; -1 suffix filled with EOS1"
)
PROFILE = "computational-route-support-v1"
REVISION = "933114058cb2604dc1bf536dbebdfcefbe83d4fc"
SETTINGS = {
    "kind": "reactiont5",
    "model_revision": REVISION,
    "device": "cpu",
    "top_k": 3,
    "max_new_tokens": 256,
    "max_input_tokens": 512,
}
FILES = {
    "model.safetensors": "69c23ca364e63e0e2be22d114652aad9f3df4e161d43936eee0dc3ce74a87fb7",
    "config.json": "a80eac6938e51ccabfc6a8b281497c6553626145e57b2268101f718fbaa034a2",
    "generation_config.json": "aed028489d41ae6a673920da9d103f2f58a650c3fbc50bd880b35e30a481dba6",
    "tokenizer.json": "511dfd8bb25c19cccae8d7d996c1ff0f1b5263207be4cdd275f23b75da8113b6",
    "tokenizer_config.json": "ac38d16067df56f00dfa2c2e0d1546a6ff08aa37814662a6fc38dea66cd13d9d",
    "special_tokens_map.json": "0f385cf62b07c1b63aaa84d4c86a313c751df7039d958a7d7bd7e723d933bcd9",
}
RUNTIME = {
    "torch": "2.8.0",
    "transformers": "4.51.3",
    "tokenizers": "0.21.4",
    "safetensors": "0.6.2",
    "rdkit": "2025.9.6",
}


def binding(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def strict_smiles(value, *, tokenization_spaces=False):
    if not isinstance(value, str) or not value:
        raise ValueError("Molecular input must be nonempty SMILES")
    text = value.replace(" ", "") if tokenization_spaces else value
    if (
        not text
        or any(c.isspace() for c in text)
        or text.startswith(".")
        or text.endswith(".")
        or ".." in text
    ):
        raise ValueError("Whitespace, empty components or trailing punctuation are unsupported")
    params = Chem.SmilesParserParams()
    params.parseName = False
    params.allowCXSMILES = False
    mol = Chem.MolFromSmiles(text, params)
    if mol is None or not mol.GetNumAtoms() or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError("Invalid or unsupported molecular SMILES")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def molecular_input(reactants, reagents):
    if not isinstance(reactants, list) or not reactants or not isinstance(reagents, list):
        raise ValueError("Forward input requires molecular reactant/reagent lists")
    data = {
        "reactants": sorted(strict_smiles(s) for s in reactants),
        "reagents": sorted(strict_smiles(s) for s in reagents),
    }
    data["prompt"] = "REACTANT:" + ".".join(data["reactants"]) + "REAGENT:" + ".".join(data["reagents"])
    return data


def validate_settings(provider):
    if (
        not isinstance(provider, dict)
        or not isinstance(provider.get("model_path"), str)
        or not provider["model_path"]
    ):
        raise ValueError(f"{PROFILE} requires an explicit local forward model_path")
    for key, expected in SETTINGS.items():
        if provider.get(key, expected) != expected:
            raise ValueError(f"{PROFILE} fixes predictor {key}={expected!r}")
    if set(provider) - set(SETTINGS) - {"model_path"}:
        raise ValueError("Unsupported strict forward-model settings")


def validate_bundle(model_path):
    """Validate numerical assets without loading a model or tokenizer.

    Hugging Face's local-download metadata is not a model search location. Only
    its known directory layout and metadata for this bundle are allowed; nested
    model/config/code files and symlinks remain unsupported.
    """
    path = Path(model_path).resolve()
    if not path.is_dir():
        raise ValueError("Pinned ReactionT5 model directory is missing")
    auxiliary = {"README.md", "trainer_state.json"}
    if set(p.name for p in path.iterdir()) - set(FILES) - auxiliary - {".cache"}:
        raise ValueError("Unexpected files in pinned ReactionT5 bundle")
    ignored = {}
    cache = path / ".cache"
    if cache.exists() or cache.is_symlink():
        if cache.is_symlink() or not cache.is_dir():
            raise ValueError("Unsupported ReactionT5 metadata cache")
        directories = {"huggingface", "huggingface/download"}
        metadata = {"huggingface/.gitignore"} | {
            "huggingface/download/" + name + ".metadata" for name in set(FILES) | auxiliary
        }
        for entry in sorted(cache.rglob("*")):
            relative = entry.relative_to(cache).as_posix()
            if entry.is_symlink() or (
                relative not in directories
                if entry.is_dir()
                else relative not in metadata or not entry.is_file()
            ):
                raise ValueError(f"Unexpected file in ReactionT5 metadata cache: {relative}")
            if entry.is_file():
                ignored[".cache/" + relative] = digest(entry)
    for name, expected in FILES.items():
        if not (path / name).is_file() or digest(path / name) != expected:
            raise ValueError(f"Pinned ReactionT5 asset mismatch: {name}")
    return ignored


class PinnedReactionT5Forward(ReactionT5Forward):
    def __init__(self, model_path, **kwargs):
        validate_settings({"model_path": model_path, **kwargs})
        path = Path(model_path).resolve()
        ignored_metadata = validate_bundle(path)
        runtime = {name: importlib.metadata.version(name) for name in RUNTIME}
        if platform.python_version() != "3.10.12" or runtime != RUNTIME:
            raise ValueError(f"Pinned ReactionT5 runtime mismatch: {runtime}")
        generation_source = Path(
            importlib.metadata.distribution("transformers").locate_file("transformers/generation/utils.py")
        )
        if digest(generation_source) != GENERATION_SOURCE_SHA256:
            raise ValueError("Pinned Transformers generation implementation mismatch")
        super().__init__(model_path, **{k: v for k, v in SETTINGS.items() if k != "kind"})
        self.torch.set_num_threads(1)
        self.provenance.update(
            profile=PROFILE,
            ignored_non_model_metadata_sha256=ignored_metadata,
            python=platform.python_version(),
            platform=platform.platform(),
            torch_threads=self.torch.get_num_threads(),
            thread_environment={
                key: os.environ.get(key)
                for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")
            },
            runtime=runtime,
            generation_source={"path": str(generation_source), "sha256": GENERATION_SOURCE_SHA256},
            padding_policy=PADDING_POLICY,
            adapter_sha256=digest(Path(__file__)),
            base_adapter_sha256=digest(Path(__file__).with_name("evidence.py")),
            normalization="canonical-isomeric-unmapped-v1; ASCII tokenization spaces only; no component or punctuation repair",
            independence="Selection support; training independence unknown; not post-selection or experimental validation",
        )

    def prepare(self, reactants, reagents):
        data = molecular_input(reactants, reagents)
        tokens, size = self._inputs(data["reactants"], data["reagents"], {})
        ids = tokens["input_ids"].tolist()[0]
        # These two registered special tokens are the pinned molecular field markers.
        allowed = {self.tokenizer.eos_token_id, 266, 267}
        if (
            not ids
            or ids[0] != 266
            or ids.count(266) != 1
            or ids.count(267) != 1
            or ids[-1] != self.tokenizer.eos_token_id
        ):
            raise ValueError("Forward input does not match the pinned molecular field serialization")
        if any(i in set(self.tokenizer.all_special_ids) - allowed for i in ids):
            raise ValueError("Forward input contains unknown or unsupported special tokens")
        return {**data, "input_token_ids": ids, "input_tokens": size}

    def inspect_sequences(self, sequences, beam_indices):
        if not isinstance(sequences, list) or len(sequences) != 3:
            raise ValueError("Strict forward prediction requires exactly three raw beams")
        if not isinstance(beam_indices, list) or len(beam_indices) != 3:
            raise ValueError("Strict forward prediction requires actual beam_indices metadata")
        results, output_tokens = [], 0
        eos, pad = self.tokenizer.eos_token_id, self.tokenizer.pad_token_id
        start = self.model.config.decoder_start_token_id
        if (start, pad, eos) != (0, 0, 1):
            raise ValueError("Decoder token identities differ from the pinned runtime")
        width = len(sequences[0]) if isinstance(sequences[0], list) else 0
        for rank, (sequence, ancestry) in enumerate(zip(sequences, beam_indices), 1):
            if (
                not isinstance(sequence, list)
                or not sequence
                or sequence[0] != start
                or any(isinstance(i, bool) or not isinstance(i, int) or i < 0 for i in sequence)
                or len(sequence) != width
                or not 2 <= width <= self.max_new_tokens + 1
            ):
                raise ValueError("Malformed, nonrectangular or oversized raw decoder sequence")
            if (
                not isinstance(ancestry, list)
                or len(ancestry) != width - 1
                or any(
                    isinstance(i, bool) or not isinstance(i, int) or i not in {-1, 0, 1, 2} for i in ancestry
                )
            ):
                raise ValueError("Malformed or misaligned beam_indices metadata")
            generated_length = next((i for i, value in enumerate(ancestry) if value == -1), len(ancestry))
            if not generated_length or any(value != -1 for value in ancestry[generated_length:]):
                raise ValueError("Beam ancestry requires a nonempty contiguous generated prefix")
            generated = sequence[1 : 1 + generated_length]
            padding = sequence[1 + generated_length :]
            # The pinned 4.51.3 implementation uses `pad_token_id or eos_id`:
            # zero PAD falls through to EOS. Only actual -1 ancestry witnesses
            # may identify these positions as fill; EOS itself is no witness.
            if any(token != eos for token in padding):
                raise ValueError("Beam suffix differs from pinned runtime EOS fill")
            complete = generated[-1] == eos and generated.count(eos) == 1
            content = generated[:-1] if complete else list(generated)
            output_tokens += generated_length
            raw = self.tokenizer.decode(
                sequence, skip_special_tokens=False, clean_up_tokenization_spaces=False
            )
            decoded = self.tokenizer.decode(
                content, skip_special_tokens=False, clean_up_tokenization_spaces=False
            )
            item = {
                "rank": rank,
                "token_ids": list(sequence),
                "beam_indices": list(ancestry),
                "generated_token_ids": generated,
                "padding_token_ids": padding,
                "generated_length": generated_length,
                "content_token_ids": content,
                "original_decoded_text": raw,
                "raw_product": decoded,
                "complete": complete,
                "product": None,
            }
            try:
                if not complete:
                    raise ValueError("Incomplete or inconsistent generated EOS decoding")
                if any(i in self.tokenizer.all_special_ids for i in content):
                    raise ValueError("Unknown or unexpected special-token content")
                item["product"] = strict_smiles(decoded, tokenization_spaces=True)
            except ValueError as exc:
                item["invalid_reason"] = str(exc)
            results.append(item)
        return results, output_tokens

    def predict_raw(self, prepared):
        inputs, _ = self._inputs(prepared["reactants"], prepared["reagents"], {})
        if inputs["input_ids"].tolist()[0] != prepared["input_token_ids"]:
            raise ValueError("Forward input token identity changed")
        with self.torch.inference_mode():
            output = self.model.generate(
                **inputs.to("cpu"),
                do_sample=False,
                num_beams=3,
                num_return_sequences=3,
                max_new_tokens=256,
                return_dict_in_generate=True,
                output_scores=False,
            )
        sequences = output.sequences.tolist()
        ancestry = getattr(output, "beam_indices", None)
        beam_indices = ancestry.tolist() if ancestry is not None else None
        raw = {
            "input": copy.deepcopy(prepared),
            "sequences": sequences,
            "beam_indices": beam_indices,
            "predictor_binding_sha256": binding(self.provenance),
        }
        try:
            predictions, tokens = self.inspect_sequences(sequences, beam_indices)
            raw.update(predictions=predictions, output_tokens=tokens)
        except ValueError as exc:
            # Retain actual malformed tensors; validation will reject and the
            # evaluator charges the reserved allowance, without another call.
            raw.update(predictions=[], output_tokens=3 * self.max_new_tokens, decoding_error=str(exc))
        return raw

    def validate_raw(self, raw, prepared):
        if raw.get("input") != prepared or raw.get("predictor_binding_sha256") != binding(self.provenance):
            raise ValueError("Forward prediction input/model binding mismatch")
        predictions, tokens = self.inspect_sequences(raw.get("sequences"), raw.get("beam_indices"))
        if raw.get("predictions") != predictions or raw.get("output_tokens") != tokens:
            raise ValueError("Forward raw token/decoded evidence mismatch")
        return predictions
