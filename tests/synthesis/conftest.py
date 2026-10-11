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
