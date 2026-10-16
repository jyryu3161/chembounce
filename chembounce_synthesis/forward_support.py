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


class PinnedReactionT5Forward(ReactionT5Forward):
    def __init__(self, model_path, **kwargs):
        validate_settings({"model_path": model_path, **kwargs})
        path = Path(model_path).resolve()
        if not path.is_dir():
            raise ValueError("Pinned ReactionT5 model directory is missing")
        if set(p.name for p in path.iterdir()) - set(FILES) - {"README.md", "trainer_state.json"}:
            raise ValueError("Unexpected files in pinned ReactionT5 bundle")
        for name, expected in FILES.items():
            if not (path / name).is_file() or digest(path / name) != expected:
                raise ValueError(f"Pinned ReactionT5 asset mismatch: {name}")
        runtime = {name: importlib.metadata.version(name) for name in RUNTIME}
        if platform.python_version() != "3.10.12" or runtime != RUNTIME:
            raise ValueError(f"Pinned ReactionT5 runtime mismatch: {runtime}")
        super().__init__(model_path, **{k: v for k, v in SETTINGS.items() if k != "kind"})
        self.torch.set_num_threads(1)
        self.provenance.update(
            profile=PROFILE,
            python=platform.python_version(),
            platform=platform.platform(),
            torch_threads=self.torch.get_num_threads(),
            thread_environment={
                key: os.environ.get(key)
                for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")
            },
            runtime=runtime,
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

    def inspect_sequences(self, sequences):
        if not isinstance(sequences, list) or len(sequences) != 3:
            raise ValueError("Strict forward prediction requires exactly three raw beams")
        results, output_tokens = [], 0
        eos, pad = self.tokenizer.eos_token_id, self.tokenizer.pad_token_id
        start = self.model.config.decoder_start_token_id
        for rank, sequence in enumerate(sequences, 1):
            if (
                not isinstance(sequence, list)
                or not sequence
                or sequence[0] != start
                or any(isinstance(i, bool) or not isinstance(i, int) or i < 0 for i in sequence)
                or len(sequence) > self.max_new_tokens + 1
            ):
                raise ValueError("Malformed or oversized raw decoder sequence")
            generated = sequence[1:]
            complete = eos in generated
            if complete:
                end = generated.index(eos)
                content = generated[:end]
                tail_valid = all(i == pad for i in generated[end + 1 :])
                charged = end + 1
            else:
                # Without EOS every returned position was generated. PAD tokens
                # can be model output, so neither hide them nor refund their cost.
                content = list(generated)
                tail_valid, charged = True, len(generated)
            output_tokens += charged
            raw = self.tokenizer.decode(
                sequence, skip_special_tokens=False, clean_up_tokenization_spaces=False
            )
            decoded = self.tokenizer.decode(
                content, skip_special_tokens=False, clean_up_tokenization_spaces=False
            )
            item = {
                "rank": rank,
                "token_ids": list(sequence),
                "content_token_ids": content,
                "original_decoded_text": raw,
                "raw_product": decoded,
                "complete": complete,
                "product": None,
            }
            try:
                if not complete:
                    raise ValueError("Incomplete decoding")
                if not tail_valid or any(i in self.tokenizer.all_special_ids for i in content):
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
        predictions, tokens = self.inspect_sequences(sequences)
        return {
            "input": copy.deepcopy(prepared),
            "sequences": sequences,
            "predictions": predictions,
            "output_tokens": tokens,
            "predictor_binding_sha256": binding(self.provenance),
        }

    def validate_raw(self, raw, prepared):
        if raw.get("input") != prepared or raw.get("predictor_binding_sha256") != binding(self.provenance):
            raise ValueError("Forward prediction input/model binding mismatch")
        predictions, tokens = self.inspect_sequences(raw.get("sequences"))
        if raw.get("predictions") != predictions or raw.get("output_tokens") != tokens:
            raise ValueError("Forward raw token/decoded evidence mismatch")
        return predictions
