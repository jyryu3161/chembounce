"""Scoped observation, never a replacement cache or scientific decision maker."""

from collections import Counter
from dataclasses import fields, is_dataclass
import functools
import hashlib
import inspect
import json
import marshal
import math
import os
from pathlib import Path
import sys
import time


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compact(value):
    """Identity/coordinate serialization only; no perception, RNG or mutation."""
    if isinstance(value, float) and not math.isfinite(value):
        return {"nonfinite_float": repr(value)}
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if is_dataclass(value):
        return {field.name: compact(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, dict):
        return {str(k): compact(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [compact(v) for v in value]
    if hasattr(value, "GetConformers"):
        raw = [[c.GetId(), c.GetPositions().tolist()] for c in value.GetConformers()]
        return {
            "type": "saved_molecule_coordinates",
            "coordinates_sha256": binding(compact(raw)),
            "atoms": value.GetNumAtoms(),
        }
    if hasattr(value, "tolist"):
        return value.tolist()
    if hasattr(value, "policy") and hasattr(value, "sha256"):
        return {"policy_sha256": value.sha256}
    return {"type": type(value).__module__ + "." + type(value).__qualname__}


def binding(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def implementation_identity(value):
    """Bind wrapper chains as actually installed, not functools.wraps metadata."""
    chain = []
    seen = set()
    while value is not None:
        if id(value) in seen:
            raise ValueError("Callable wrapper cycle")
        seen.add(id(value))
        code = getattr(value, "__code__", None)
        module = sys.modules.get(getattr(value, "__module__", ""))
        path = getattr(module, "__file__", None)
        chain.append(
            {
                "type": type(value).__module__ + "." + type(value).__qualname__,
                "code_sha256": hashlib.sha256(marshal.dumps(code)).hexdigest() if code else None,
                "defaults": compact(getattr(value, "__defaults__", None)),
                "implementation_module_path": str(Path(path).resolve()) if path else None,
                "implementation_module_sha256": digest(path) if path else None,
                "native_limit": None
                if code
                else "No Python code object; bound module/binary and callable type only",
            }
        )
        value = getattr(value, "__wrapped__", None)
    return chain


class Observer:
    def __init__(self, path, *, expected_inventory=None, module_prefix="chembounce_synthesis"):
        self.path = Path(path)
        self.stream = self.path.open("x")
        self.expected = expected_inventory
        self.prefix = module_prefix
        self.inventory, self.restores, self.events = [], [], Counter()
        self.phase = "startup"
        self.sequence = 0
        self.active = []
        self.contexts = []

    def emit(self, row):
        self.stream.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")
        self.stream.flush()
        os.fsync(self.stream.fileno())

    def install(self, owner, name, *, kind, physical=False):
        original = getattr(owner, name)
        if hasattr(original, "__capacity_observer__"):
            raise ValueError("Unexpected already installed observer")
        aliases = [(owner, name)]
        for modname, module in sorted(sys.modules.copy().items()):
            if module is None or not modname.startswith(self.prefix):
                continue
            for attr, value in list(vars(module).items()):
                if value is original and (module, attr) not in aliases:
                    aliases.append((module, attr))
        owner_module = sys.modules[owner.__module__] if inspect.isclass(owner) else owner
        source = Path(owner_module.__file__).resolve()
        try:
            signature = str(inspect.signature(original))
        except (ValueError, TypeError):
            signature = "builtin_signature_unavailable"
        entry = {
            "kind": kind,
            "physical": physical,
            "owner": owner.__name__,
            "attribute": name,
            "source_path": str(source),
            "source_sha256": digest(source),
            "signature": signature,
            "callable": getattr(original, "__qualname__", name),
            "implementation": implementation_identity(original),
            "aliases": [o.__name__ + "." + a for o, a in aliases],
            "wrapper_sha256": digest(__file__),
            "install_index": len(self.inventory),
        }
        if self.expected is not None and (
            len(self.expected) <= len(self.inventory) or self.expected[len(self.inventory)] != entry
        ):
            raise ValueError("Observer inventory differs from frozen authority")
        self.inventory.append(entry)
        observer = self

        @functools.wraps(original)
        def observed(*args, **kwargs):
            observer.sequence += 1
            ident = observer.sequence
            key_args = compact(args[1:] if inspect.isclass(owner) else args)
            # Absolute deadline is recorded but excluded from comparable scientific keys.
            key_kwargs = {k: v for k, v in kwargs.items() if k != "deadline_monotonic"}
            scientific_context = [
                item
                for item in observer.contexts
                if item["kind"] in {"ensemble", "score_candidate", "reference_request", "refinement"}
            ]
            comparable = kind not in {"embedding", "overlay"} or bool(scientific_context)
            key = (
                binding(
                    {
                        "kind": kind,
                        "args": key_args,
                        "kwargs": compact(key_kwargs),
                        "scientific_context": scientific_context,
                    }
                )
                if comparable
                else None
            )
            start, cpu = time.monotonic(), time.process_time()
            cache = original.cache_info()._asdict() if hasattr(original, "cache_info") else None
            row = {
                "event": "started",
                "id": ident,
                "kind": kind,
                "physical": physical,
                "phase": "screening"
                if kind == "score_candidate" and observer.phase == "workflow"
                else observer.phase,
                "parent": observer.active[-1] if observer.active else None,
                "call_key": key,
                "comparable": comparable,
                "scientific_context": scientific_context,
                "arguments": key_args,
                "kwargs": compact(kwargs),
                "monotonic": start,
                "cache_before": cache,
            }
            observer.emit(row)
            observer.events[kind + ":started"] += 1
            observer.active.append(ident)
            observer.contexts.append({"kind": kind, "args": key_args, "kwargs": compact(key_kwargs)})
            previous_phase = observer.phase
            if kind == "refinement":
                observer.phase = "final_refinement"
            try:
                result = original(*args, **kwargs)
            except BaseException as exc:
                observer.emit(
                    {
                        "event": "failed",
                        "id": ident,
                        "kind": kind,
                        "call_key": key,
                        "error": {"type": type(exc).__name__, "message": str(exc)},
                        "wall_seconds": time.monotonic() - start,
                        "cpu_seconds": time.process_time() - cpu,
                    }
                )
                observer.events[kind + ":failed"] += 1
                raise
            else:
                metadata = {}
                if kind == "embedding":
                    metadata = {"returned_conformers": len(result)}
                elif kind == "generation":
                    sequences = getattr(result, "sequences", result)
                    metadata = {"sequences": compact(sequences)}
                if hasattr(original, "cache_info"):
                    metadata["cache_after"] = original.cache_info()._asdict()
                observer.emit(
                    {
                        "event": "returned",
                        "id": ident,
                        "kind": kind,
                        "call_key": key,
                        "wall_seconds": time.monotonic() - start,
                        "cpu_seconds": time.process_time() - cpu,
                        **metadata,
                    }
                )
                observer.events[kind + ":returned"] += 1
                return result
            finally:
                observer.active.pop()
                observer.contexts.pop()
                observer.phase = previous_phase

        observed.__capacity_observer__ = self
        # The original LRU object still owns state; this is a transparent call alias.
        for attr in ("cache_info", "cache_parameters"):
            if hasattr(original, attr):
                setattr(observed, attr, getattr(original, attr))
        for target, attr in aliases:
            self.restores.append((target, attr, original, observed))
            setattr(target, attr, observed)
        return entry

    def seal(self):
        if self.expected is None or self.inventory != self.expected:
            raise ValueError("Complete frozen observer inventory is required")
        self.emit(
            {
                "event": "installed",
                "inventory": self.inventory,
                "inventory_sha256": binding(self.inventory),
                "runtime_identities": [
                    {"original": id(o), "wrapper": id(w), "target": t.__name__ + "." + a}
                    for t, a, o, w in self.restores
                ],
            }
        )

    def close(self):
        changed = []
        for target, attr, original, wrapper in reversed(self.restores):
            if getattr(target, attr) is not wrapper:
                changed.append(target.__name__ + "." + attr)
            setattr(target, attr, original)
        self.emit({"event": "closed", "counts": dict(self.events), "unexpected_replacements": changed})
        self.stream.close()
        if changed:
            raise ValueError("Observed callables changed during execution")


def install_scientific(observer):
    """Called only by a separately authorized worker, after module-path checks.

    No model/asset constructors are invoked to install these wrappers. Complete
    dependency imports and inventory are pinned prospectively in the literal.
    """
    from rdkit.Chem import AllChem
    from transformers.generation.utils import GenerationMixin
    from chembounce_synthesis import catalog, reactions, conformers, selective_scoring, refinement, workflow
    from chembounce_synthesis import forward_support

    for owner, name, kind, physical in (
        (catalog.Catalog, "__init__", "catalog_startup", False),
        (reactions.Rules, "__init__", "rules_startup", False),
        (forward_support.PinnedReactionT5Forward, "__init__", "model_startup", False),
        (forward_support.PinnedReactionT5Forward, "predict_raw", "predict_raw", False),
        (GenerationMixin, "generate", "generation", True),
        (conformers, "prepare_reference", "reference_request", False),
        (refinement, "refine_candidate", "refinement", False),
        (conformers, "ensemble", "ensemble", False),
        (AllChem, "EmbedMultipleConfs", "embedding", True),
        (conformers, "_overlay", "overlay", True),
        (selective_scoring, "score_candidate", "score_candidate", False),
    ):
        observer.install(owner, name, kind=kind, physical=physical)
    # Imported aliases in workflow/refinement/selective_scoring are included by
    # identity discovery; import them before installation, never after it.
    _ = refinement, workflow
    observer.seal()
