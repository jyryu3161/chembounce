"""Allowlisted reaction rules and bounded, independently sanitized replay."""

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

from .common import canonical, digest, molecule, positive


@dataclass
class ReactionRule:
    id: str
    smarts: str
    source: str
    reaction: object
    validation_warnings: int = 0

    @property
    def roles(self):
        return tuple(Chem.MolToSmarts(m) for m in self.reaction.GetReactants())

    @property
    def unmapped_product_template_atoms(self):
        """Heavy atoms lacking an explicit reactant atom-map source in the rule."""
        reactant_maps = {
            a.GetAtomMapNum() for m in self.reaction.GetReactants() for a in m.GetAtoms() if a.GetAtomMapNum()
        }
        return [
            {"atom_index": a.GetIdx(), "element": a.GetSymbol(), "map_number": a.GetAtomMapNum()}
            for m in self.reaction.GetProducts()
            for a in m.GetAtoms()
            if a.GetAtomicNum() > 1 and a.GetAtomMapNum() not in reactant_maps
        ]


class Rules:
    def __init__(self, filename):
        self.path = Path(filename).resolve()
        self.sha256 = digest(self.path)
        data = json.loads(self.path.read_text())
        if data.get("schema_version") != 1 or not isinstance(data.get("rules"), list):
            raise ValueError("Rules require schema_version=1 and a rules list")
        self.rules, self.by_smarts = {}, {}
        for row in data["rules"]:
            ident, smarts, source = row["id"], row["smarts"], row["source"]
            if not all(isinstance(x, str) and x.strip() for x in (ident, smarts, source)):
                raise ValueError("Rules require nonempty id, smarts, source")
            if ident in self.rules or smarts in self.by_smarts:
                raise ValueError("Duplicate reaction rule ID or SMARTS")
            rxn = rdChemReactions.ReactionFromSmarts(smarts)
            if (
                rxn is None
                or not 1 <= rxn.GetNumReactantTemplates() <= 3
                or rxn.GetNumProductTemplates() != 1
            ):
                raise ValueError(f"Rule {ident}: require 1–3 reactants and one principal product template")
            warnings, errors = rxn.Validate(silent=True)
            if errors:
                raise ValueError(f"Invalid reaction rule: {ident}")
            # Template warnings are retained in route evidence; avoid repeating
            # the same known warnings for every registry load on stderr.
            with rdBase.BlockLogs():
                rxn.Initialize()
            # A whole-molecule role check is insufficient when only one of several
            # sites has the required stereochemistry. Enforce it during replay too.
            rxn.GetSubstructParams().useChirality = True
            self.rules[ident] = ReactionRule(ident, smarts, source, rxn, warnings)
            self.by_smarts[smarts] = ident
        if not self.rules:
            raise ValueError("Reaction allowlist is empty")
        self.apply = lru_cache(maxsize=8192)(self.apply)

    def get(self, ident):
        if ident not in self.rules:
            raise ValueError(f"Unknown reaction rule: {ident}")
        return self.rules[ident]

    def resolve(self, text):
        text = text.removeprefix("<rxn>").removesuffix("</rxn>").strip()
        if text in self.rules:
            return text
        if text not in self.by_smarts:
            raise ValueError("Model proposed a reaction outside the allowlist")
        return self.by_smarts[text]

    def apply(self, ident, reactants, product_limit=64):
        """Reactants are in rule-slot order; return all unique valid products up to cap.

        Count truncation explicitly. Multi-product rules are rejected on import.
        Stereochemistry comes from RDKit reaction execution, never from model text.
        """
        positive(product_limit, "reaction_product_limit")
        rule = self.get(ident)
        if len(reactants) != len(rule.roles):
            raise ValueError(f"Wrong reactant count for {ident}")
        mols = tuple(molecule(s) for s in reactants)
        if not all(
            m.HasSubstructMatch(q, useChirality=True) for m, q in zip(mols, rule.reaction.GetReactants())
        ):
            return (), False
        raw = rule.reaction.RunReactants(mols, maxProducts=product_limit + 1)
        products = set()
        for entry in raw[:product_limit]:
            try:
                mol = entry[0]
                Chem.SanitizeMol(mol)
                Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                products.add(canonical(Chem.MolToSmiles(mol, isomericSmiles=True)))
            except (ValueError, RuntimeError):
                continue
        return tuple(sorted(products)), len(raw) > product_limit
