import json

import pytest

from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.reactions import Rules

AMIDE = "[C:1](=[O:2])[O;H1:3].[N;H1,H2:4]>>[C:1](=[O:2])[N:4]"
HIT = "CC(=O)Nc1ccccc1"
DESIGN = "CC(=O)Nc1ccncc1"


@pytest.fixture
def assets(tmp_path):
    catalog_path = tmp_path / "catalog.tsv"
    rows = [
        ("acetic", "CC(=O)O"),
        ("propionic", "CCC(=O)O"),
        ("butyric", "CCCC(=O)O"),
        ("benzoic", "O=C(O)c1ccccc1"),
        ("pyridyl", "Nc1ccncc1"),
        ("pyridyl-2", "Nc1ccccn1"),
        ("pyridyl-3", "Nc1cccnc1"),
        ("phenyl", "Nc1ccccc1"),
        ("ethylamine", "CCN"),
        ("glycine", "NCC(=O)O"),
        ("alanine_s", "N[C@@H](C)C(=O)O"),
        ("alanine_r", "N[C@H](C)C(=O)O"),
    ]
    catalog_path.write_text(
        "id\tsmiles\tsource\tsnapshot\n"
        + "".join(f"{ident}\t{s}\ttest fixture; not a vendor catalog\t2026-10-11\n" for ident, s in rows)
    )
    rules_path = tmp_path / "rules.json"
    rules_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "rules": [
                    {
                        "id": "amide",
                        "smarts": AMIDE,
                        "source": "unit-test reaction fixture; conditions unspecified",
                    }
                ],
            }
        )
    )
    return Catalog(catalog_path), Rules(rules_path)


@pytest.fixture
def proposal():
    return {
        "reactions": [
            {
                "reaction_number": 1,
                "reaction_template": f"<rxn>{AMIDE}</rxn>",
                "reactants": ["Nc1ccncc1", "CC(=O)O"],
                "product": DESIGN,
            }
        ],
        "building_blocks": ["<bb>CC(=O)O</bb>", "<bb>Nc1ccncc1</bb>"],
    }


@pytest.fixture
def route():
    return {
        "schema_version": 1,
        "root": "r0",
        "nodes": [
            {"id": "b0", "building_block_id": "acetic"},
            {"id": "b1", "building_block_id": "pyridyl"},
            {"id": "r0", "rule_id": "amide", "reactants": ["b0", "b1"], "product": DESIGN},
        ],
    }


@pytest.fixture
def strict_predictor(monkeypatch):
    """Controlled raw beams; real strict decoder/comparison/audit code remains active."""
    import copy
    from types import SimpleNamespace

    from chembounce_synthesis import computational_support
    from chembounce_synthesis.forward_support import PinnedReactionT5Forward, binding, molecular_input

    class Tokenizer:
        eos_token_id, pad_token_id = 1, 0
        all_special_ids = [0, 1, 2]

        def decode(self, ids, **kwargs):
            special = {0: "<pad>", 1: "</s>", 2: "<unk>"}
            return "".join(special[i] if i in special else chr(i - 3) for i in ids)

    class Controlled(PinnedReactionT5Forward):
        def __init__(self):
            self.tokenizer = Tokenizer()
            self.model = SimpleNamespace(config=SimpleNamespace(decoder_start_token_id=0))
            self.max_new_tokens = 256
            self.provenance = {
                "profile": "computational-route-support-v1",
                "fixture": "controlled raw beams; not model evidence",
            }
            self.calls = []
            self.products = None

        def prepare(self, reactants, reagents):
            data = molecular_input(reactants, reagents)
            ids = [ord(c) + 3 for c in data["prompt"]] + [1]
            return {**data, "input_token_ids": ids, "input_tokens": len(ids)}

        def predict_raw(self, prepared):
            self.calls.append(copy.deepcopy(prepared))
            amine = next((s for s in prepared["reactants"] if s.startswith("N")), "N")
            products = self.products if self.products is not None else ["CC(=O)" + amine, "CC", "CO"]
            if callable(products):
                products = products(prepared)
            sequences = [
                ([0] + [ord(c) + 3 for c in value] + [1]) if isinstance(value, str) else value
                for value in products
            ]
            predictions, tokens = self.inspect_sequences(sequences)
            return {
                "input": copy.deepcopy(prepared),
                "sequences": sequences,
                "predictions": predictions,
                "output_tokens": tokens,
                "predictor_binding_sha256": binding(self.provenance),
            }

    predictor = Controlled()
    monkeypatch.setattr(computational_support, "PinnedReactionT5Forward", lambda **kwargs: predictor)
    return predictor
