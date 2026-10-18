"""Optional local SynLlama inference; never execute upstream Python code.

Official .bin weights use Torch's restricted weights-only loader (Torch >=2.6).
General-purpose pickle objects and remote model code are never loaded.
"""

import ast
import importlib.metadata
import json
import logging
import math
import time
from pathlib import Path

from .common import canonical, digest, positive
from .proposal_json import (
    POLICY as SMILES_JSON_POLICY,
    bind_row,
    compatibility_origin,
    finalize_timing,
    provider_identity_sha256,
    verify_saved_row,
)

LOGGER = logging.getLogger(__name__)


def load_prompt(filename):
    """Read an upstream vars.py TEMPLATE literal or a user-supplied JSON template."""
    path = Path(filename)
    if path.suffix == ".json":
        template = json.loads(path.read_text())
    else:
        tree = ast.parse(path.read_text())
        assignments = [
            n
            for n in tree.body
            if isinstance(n, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id == "TEMPLATE" for t in n.targets)
        ]
        if len(assignments) != 1:
            raise ValueError("Prompt source must contain one literal TEMPLATE assignment")
        template = ast.literal_eval(assignments[0].value)
    if (
        not isinstance(template, dict)
        or not all(isinstance(template.get(k), str) for k in ("instruction", "input"))
        or "SMILES_STRING" not in template["input"]
    ):
        raise ValueError("Prompt requires instruction and input containing SMILES_STRING")
    if "response_json_dialect" in template and (
        path.suffix != ".json" or template["response_json_dialect"] != SMILES_JSON_POLICY
    ):
        raise ValueError(f"JSON prompt response_json_dialect must be {SMILES_JSON_POLICY}")
    return template


class SynLlama:
    def __init__(self, model_path, prompt_path, device="cpu", samples=4, max_new_tokens=1600, seed=42):
        positive(samples, "samples")
        positive(max_new_tokens, "max_new_tokens")
        if not isinstance(seed, int) or not 0 <= seed <= 2**31 - 1:
            raise ValueError("Invalid inference seed")
        self.model_path = Path(model_path).resolve()
        self.prompt_path = Path(prompt_path).resolve()
        self.template = load_prompt(self.prompt_path)
        self.response_json_dialect = self.template.get("response_json_dialect", "strict-json-v1")
        safe_weights = sorted(self.model_path.glob("*.safetensors"))
        weights = safe_weights or sorted(self.model_path.glob("pytorch_model*.bin"))
        if not weights or not (self.model_path / "config.json").is_file():
            raise ValueError("SynLlama requires local Hugging Face weights and config.json")
        self.samples, self.max_new_tokens, self.seed, self.device = samples, max_new_tokens, seed, device
        self.provenance = {
            "provider": "synllama_local",
            "model_path": str(self.model_path),
            "model_files_sha256": {
                p.name: digest(p) for p in sorted(self.model_path.iterdir()) if p.is_file()
            },
            "prompt_path": str(self.prompt_path),
            "prompt_sha256": digest(self.prompt_path),
            "device": device,
            "samples": samples,
            "max_new_tokens": max_new_tokens,
            "seed": seed,
            "use_cache": True,
            "response_json_dialect": self.response_json_dialect,
        }
        if self.response_json_dialect == SMILES_JSON_POLICY:
            self.provenance["response_decoder_source_sha256"] = digest(
                Path(__file__).with_name("proposal_json.py")
            )
        try:
            import torch
            from transformers import AutoModelForCausalLM, AutoTokenizer
        except ImportError as exc:
            raise RuntimeError(
                "Install the synllama optional dependency group to run model inference"
            ) from exc
        self.torch = torch
        if not safe_weights and tuple(int(p) for p in torch.__version__.split(".")[:2]) < (2, 6):
            raise RuntimeError(
                "Official .bin weights require Torch >=2.6 for restricted weights-only loading"
            )
        self.tokenizer = AutoTokenizer.from_pretrained(
            str(self.model_path), local_files_only=True, trust_remote_code=False
        )
        dtype = torch.float16 if device.startswith("cuda") else torch.float32
        self.model = (
            AutoModelForCausalLM.from_pretrained(
                str(self.model_path),
                torch_dtype=dtype,
                local_files_only=True,
                trust_remote_code=False,
                use_safetensors=bool(safe_weights),
                weights_only=True,
            )
            .to(device)
            .eval()
        )
        if self.tokenizer.pad_token_id is None:
            self.tokenizer.pad_token_id = self.tokenizer.eos_token_id
        for package in ("torch", "transformers", "safetensors"):
            self.provenance[package + "_version"] = importlib.metadata.version(package)
        self.provenance["dtype"] = str(dtype)
        self.provenance["hardware"] = (
            torch.cuda.get_device_name(torch.device(device)) if device.startswith("cuda") else "cpu"
        )

    def propose(self, smiles, deadline_monotonic=None):
        if deadline_monotonic is not None and (
            isinstance(deadline_monotonic, bool)
            or not isinstance(deadline_monotonic, (int, float))
            or not math.isfinite(deadline_monotonic)
        ):
            raise ValueError("deadline_monotonic must be finite")
        smiles = canonical(smiles)
        prompt = f"### Instruction:\n{self.template['instruction']}\n\n### Input:\n{self.template['input'].replace('SMILES_STRING', smiles)}\n\n### Response: \n"
        inputs = self.tokenizer(prompt, return_tensors="pt", return_token_type_ids=False).to(self.device)
        if inputs.input_ids.shape[1] + self.max_new_tokens > self.model.config.max_position_embeddings:
            raise ValueError("Prompt plus output budget exceeds model context length")
        modes = [(0.1, 0.1), (0.6, 0.5), (1.0, 0.7), (1.5, 0.9)]
        eos = [self.tokenizer.eos_token_id]
        eot = self.tokenizer.convert_tokens_to_ids("<|eot_id|>")
        if eot is not None and eot != self.tokenizer.unk_token_id and eot not in eos:
            eos.append(eot)
        responses = []
        devices = [self.torch.device(self.device).index or 0] if self.device.startswith("cuda") else []
        stopping = {}
        if deadline_monotonic is not None:
            from transformers import StoppingCriteria, StoppingCriteriaList

            class DeadlineStop(StoppingCriteria):
                def __call__(self, input_ids, scores, **kwargs):
                    return time.monotonic() >= deadline_monotonic

            stopping["stopping_criteria"] = StoppingCriteriaList([DeadlineStop()])
        longest_sample = 0.0
        provider_identity = None
        for index in range(self.samples):
            if deadline_monotonic is not None and deadline_monotonic - time.monotonic() <= max(
                0.1, longest_sample * 1.2
            ):
                break
            sample_started = time.monotonic()
            if self.response_json_dialect == SMILES_JSON_POLICY and provider_identity is None:
                provider_identity = provider_identity_sha256(self.provenance)
            temperature, top_p = modes[index % len(modes)]
            LOGGER.info("SynLlama sample %s/%s for %s", index + 1, self.samples, smiles)
            with self.torch.random.fork_rng(devices=devices), self.torch.inference_mode():
                self.torch.random.default_generator.manual_seed(self.seed + index)
                if devices:
                    with self.torch.cuda.device(self.device):
                        self.torch.cuda.manual_seed(self.seed + index)
                generated = self.model.generate(
                    **inputs,
                    use_cache=True,
                    do_sample=True,
                    temperature=temperature,
                    top_p=top_p,
                    max_new_tokens=self.max_new_tokens,
                    eos_token_id=eos,
                    pad_token_id=self.tokenizer.pad_token_id,
                    **stopping,
                )
            tokens = generated[0, inputs.input_ids.shape[1] :]
            response = self.tokenizer.decode(tokens, skip_special_tokens=True)
            sampling = {
                "temperature": temperature,
                "top_p": top_p,
                "seed": self.seed + index,
                "generated_tokens": len(tokens),
                "ended_with_eos": int(tokens[-1]) in eos if len(tokens) else False,
            }
            if self.response_json_dialect == SMILES_JSON_POLICY:
                sampling["sample_index"] = index
                row = bind_row(
                    response,
                    sampling,
                    self.provenance["prompt_sha256"],
                    self.provenance["response_decoder_source_sha256"],
                    design=smiles,
                    provider_identity=provider_identity,
                    deadline_monotonic=deadline_monotonic,
                )
                if (
                    deadline_monotonic is not None
                    and time.monotonic() >= deadline_monotonic
                    and row["compatibility"]["status"] != "not_evaluated_deadline"
                ):
                    row = bind_row(
                        response,
                        sampling,
                        self.provenance["prompt_sha256"],
                        self.provenance["response_decoder_source_sha256"],
                        design=smiles,
                        provider_identity=provider_identity,
                        deadline_monotonic=time.monotonic(),
                    )
                finished = time.monotonic()
                finalize_timing(
                    row,
                    finished - sample_started,
                    deadline_monotonic is not None and finished >= deadline_monotonic,
                )
                responses.append(row)
            else:
                sampling["deadline_reached"] = (
                    deadline_monotonic is not None and time.monotonic() >= deadline_monotonic
                )
                sampling["elapsed_seconds"] = time.monotonic() - sample_started
                responses.append({"response": response, "sampling": sampling})
            longest_sample = max(longest_sample, time.monotonic() - sample_started)
        return responses


class SavedProposals:
    """Replay auditable JSON model outputs without installing an LLM runtime."""

    def __init__(self, filename):
        path = Path(filename).resolve()
        data = json.loads(path.read_text())
        if data.get("schema_version") != 1 or not isinstance(data.get("proposals"), list):
            raise ValueError("Saved proposals require schema_version=1 and a proposals list")
        self.proposals = {}
        recorded_provider = data.get("provider", {})
        origin = compatibility_origin(recorded_provider)
        dialect = SMILES_JSON_POLICY if origin is not None else None
        if dialect == SMILES_JSON_POLICY and not all(
            isinstance(recorded_provider.get(key), str)
            and len(recorded_provider[key]) == 64
            and all(ch in "0123456789abcdef" for ch in recorded_provider[key])
            for key in ("prompt_sha256", "response_decoder_source_sha256")
        ):
            raise ValueError("Saved compatibility provider is missing policy provenance")
        if dialect == SMILES_JSON_POLICY and recorded_provider["response_decoder_source_sha256"] != digest(
            Path(__file__).with_name("proposal_json.py")
        ):
            raise ValueError("Saved compatibility decoder source does not match this runtime")
        if dialect == SMILES_JSON_POLICY and (
            type(origin.get("samples")) is not int
            or origin["samples"] < 1
            or type(origin.get("seed")) is not int
        ):
            raise ValueError("Saved compatibility origin is missing sample identity bounds")
        seen_samples = {}
        for row in data["proposals"]:
            design = canonical(row["design"])
            if dialect == SMILES_JSON_POLICY:
                proposal = verify_saved_row(row, recorded_provider, design=design)
                index = proposal["sampling"]["sample_index"]
                if (
                    index >= origin["samples"]
                    or type(proposal["sampling"].get("seed")) is not int
                    or proposal["sampling"]["seed"] != origin["seed"] + index
                ):
                    raise ValueError("Saved compatibility sample identity differs from origin provider")
                slots = seen_samples.setdefault(design, set())
                if index in slots:
                    raise ValueError("Duplicate saved compatibility sample index for design")
                slots.add(index)
            else:
                if any(key in row for key in ("raw_response", "compatibility")):
                    raise ValueError("Saved proposal audit requires declared compatibility policy")
                proposal = {"response": row["response"], "sampling": row.get("sampling", {})}
            self.proposals.setdefault(design, []).append(proposal)
        self.provenance = {
            "provider": "saved_proposals",
            "path": str(path),
            "sha256": digest(path),
            "recorded_provider": recorded_provider,
        }
        if dialect == SMILES_JSON_POLICY:
            self.provenance.update(
                {
                    key: recorded_provider[key]
                    for key in ("response_json_dialect", "prompt_sha256", "response_decoder_source_sha256")
                }
            )

    def propose(self, smiles, deadline_monotonic=None):
        if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
            return []
        key = canonical(smiles)
        if key not in self.proposals:
            raise ValueError(f"No saved proposals for design {key}")
        return self.proposals[key]
