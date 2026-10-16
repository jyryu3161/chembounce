"""Concrete, bounded reaction witnesses; computational support, not experiments.

Every retained product heavy atom and explicit isotopic hydrogen must originate
in one concrete reactant slot/atom. Ordinary hydrogen/protons, lost atoms and
unmodeled byproducts are not balanced. No free-text reagent attribution is used.
"""

import ast
import hashlib
import json

from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

from .common import canonical, digest, molecule, positive
from .routes import _topology, verify_route


SCHEMA = "concrete-route-audit-v1"
H = ("implicit_hydrogen",)
TETRA = {Chem.ChiralType.CHI_TETRAHEDRAL_CW: 1, Chem.ChiralType.CHI_TETRAHEDRAL_CCW: -1}
DOUBLE = {
    Chem.BondStereo.STEREOCIS: 1,
    Chem.BondStereo.STEREOZ: 1,
    Chem.BondStereo.STEREOTRANS: -1,
    Chem.BondStereo.STEREOE: -1,
}


def binding(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def stereo_state(mol):
    """Potential tetrahedral/double stereo, including unspecified route states."""
    entries = []
    for item in Chem.FindPotentialStereo(mol):
        entries.append(
            {
                "type": str(item.type),
                "index": item.centeredOn,
                "specified": str(item.specified),
                "descriptor": str(item.descriptor),
                "controlling_atoms": list(item.controllingAtoms),
            }
        )
    return {
        "centers": entries,
        "unspecified": [x for x in entries if x["specified"] != "Specified"],
        "unsupported_types": [x for x in entries if x["type"] not in {"Atom_Tetrahedral", "Bond_Double"}],
        "fully_specified": all(
            x["specified"] == "Specified" and x["type"] in {"Atom_Tetrahedral", "Bond_Double"}
            for x in entries
        ),
    }


def template_maps(reaction):
    maps, collisions = {}, []
    for slot, mol in enumerate(reaction.GetReactants()):
        for atom in mol.GetAtoms():
            number = atom.GetAtomMapNum()
            if number:
                if number in maps:
                    collisions.append(number)
                maps[number] = (slot, atom.GetIdx())
    product_maps = [
        a.GetAtomMapNum() for m in reaction.GetProducts() for a in m.GetAtoms() if a.GetAtomMapNum()
    ]
    duplicate_products = sorted({n for n in product_maps if product_maps.count(n) > 1})
    gaps = [
        {
            "template_atom_index": a.GetIdx(),
            "element": a.GetSymbol(),
            "isotope": a.GetIsotope(),
            "map_number": a.GetAtomMapNum(),
        }
        for m in reaction.GetProducts()
        for a in m.GetAtoms()
        if (a.GetAtomicNum() > 1 or a.GetIsotope()) and a.GetAtomMapNum() not in maps
    ]
    return maps, {
        "reactant_map_collisions": sorted(set(collisions)),
        "product_map_collisions": duplicate_products,
        "unbound_product_atoms": gaps,
    }


def atom_origins(product, reactants, maps):
    """Read RDKit's concrete execution metadata before any SMILES round trip."""
    origins, records, problems, used = {}, [], [], set()
    for atom in product.GetAtoms():
        index = atom.GetIdx()
        # Explicit isotope-H is included; ordinary implicit/explicit H is outside scope.
        if atom.GetAtomicNum() == 1 and not atom.GetIsotope():
            continue
        item = {
            "product_atom_index": index,
            "element": atom.GetSymbol(),
            "isotope": atom.GetIsotope(),
            "old_map_number": atom.GetIntProp("old_mapno") if atom.HasProp("old_mapno") else None,
        }
        reasons = []
        if not atom.HasProp("react_atom_idx") or not atom.HasProp("react_idx"):
            reasons.append("missing_concrete_reactant_origin")
        else:
            slot, source_index = atom.GetIntProp("react_idx"), atom.GetIntProp("react_atom_idx")
            item.update(reactant_slot=slot, reactant_atom_index=source_index)
            if not 0 <= slot < len(reactants) or not 0 <= source_index < reactants[slot].GetNumAtoms():
                reasons.append("origin_index_out_of_range")
            else:
                source = reactants[slot].GetAtomWithIdx(source_index)
                item.update(source_element=source.GetSymbol(), source_isotope=source.GetIsotope())
                if (source.GetAtomicNum(), source.GetIsotope()) != (atom.GetAtomicNum(), atom.GetIsotope()):
                    reasons.append("element_or_isotope_transmutation")
                number = item["old_map_number"]
                if number is not None and (number not in maps or maps[number][0] != slot):
                    reasons.append("template_map_origin_slot_mismatch")
                if (slot, source_index) in used:
                    reasons.append("duplicate_donor_atom_within_slot")
                used.add((slot, source_index))
                origins[index] = (slot, source_index)
        item["problems"] = reasons
        records.append(item)
        problems.extend({"product_atom_index": index, "reason": reason} for reason in reasons)
    lost = [
        {
            "reactant_slot": slot,
            "reactant_atom_index": a.GetIdx(),
            "element": a.GetSymbol(),
            "isotope": a.GetIsotope(),
        }
        for slot, mol in enumerate(reactants)
        for a in mol.GetAtoms()
        if (a.GetAtomicNum() > 1 or a.GetIsotope()) and (slot, a.GetIdx()) not in used
    ]
    return origins, {
        "status": "complete" if not problems else "unsupported",
        "atoms": records,
        "problems": problems,
        "lost_input_atoms": lost,
        "lost_input_atoms_scope": "Absent from principal product; byproducts and total mass balance unassessed",
    }


def neighbor_keys(mol, index, slot, origins=None, exclude=None):
    atom = mol.GetAtomWithIdx(index)
    keys = []
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() == exclude:
            continue
        keys.append((slot, neighbor.GetIdx()) if origins is None else origins.get(neighbor.GetIdx()))
    if atom.GetTotalNumHs() == 1:
        keys.append(H)
    return keys


def correspondence(source, product):
    """Preserve common neighbors; one leaving/new substituent may correspond."""
    if (
        len(source) != len(product)
        or None in product
        or len(set(source)) != len(source)
        or len(set(product)) != len(product)
    ):
        return None, None
    missing, added = set(source) - set(product), set(product) - set(source)
    if len(missing) > 1 or len(added) > 1:
        return None, None
    replacement = None
    converted = list(product)
    if missing:
        old, new = next(iter(missing)), next(iter(added))
        converted = [old if x == new else x for x in product]
        replacement = {"leaving_neighbor": list(old), "incoming_neighbor": list(new)}
    if set(converted) != set(source):
        return None, None
    return converted, replacement


def permutation_sign(source, reordered):
    indices = [source.index(item) for item in reordered]
    inversions = sum(indices[i] > indices[j] for i in range(len(indices)) for j in range(i + 1, len(indices)))
    return -1 if inversions % 2 else 1


def tetra_audit(product, atom_index, reactants, origins):
    atom = product.GetAtomWithIdx(atom_index)
    origin = origins.get(atom_index)
    result = {
        "type": "tetrahedral",
        "product_atom_index": atom_index,
        "classification": "unresolved",
        "requires_stereo_prediction": True,
    }
    if origin is None or atom.GetChiralTag() not in TETRA:
        return result
    slot, source_index = origin
    source = reactants[slot].GetAtomWithIdx(source_index)
    result["origin"] = {"reactant_slot": slot, "reactant_atom_index": source_index}
    result["rdkit_inversion_flag"] = (
        atom.GetIntProp("molInversionFlag") if atom.HasProp("molInversionFlag") else None
    )
    if source.GetChiralTag() not in TETRA:
        result["classification"] = "new_or_template_asserted"
        return result
    old = neighbor_keys(reactants[slot], source_index, slot)
    new = neighbor_keys(product, atom_index, slot, origins)
    aligned, replacement = correspondence(old, new)
    if aligned is None or len(old) != 4:
        result["reason"] = "insufficient_mapped_neighbor_correspondence"
        return result
    sign = permutation_sign(old, aligned)
    retained = TETRA[atom.GetChiralTag()] * sign == TETRA[source.GetChiralTag()]
    result.update(
        classification="retained" if retained else "inverted",
        requires_stereo_prediction=not retained or replacement is not None,
        source_neighbor_order=[list(k) for k in old],
        product_neighbor_order=[list(k) for k in new],
        correspondence_parity=sign,
        substituent_replacement=replacement,
        comparison="Chiral neighbor-order parity, not CIP R/S labels",
    )
    return result


def double_audit(product, bond_index, reactants, origins):
    bond = product.GetBondWithIdx(bond_index)
    ends = [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
    origin = [origins.get(i) for i in ends]
    result = {
        "type": "double_bond",
        "product_bond_index": bond_index,
        "classification": "unresolved",
        "requires_stereo_prediction": True,
    }
    if any(x is None for x in origin) or bond.GetStereo() not in DOUBLE:
        return result
    if origin[0][0] != origin[1][0]:
        result["classification"] = "new_or_template_asserted"
        return result
    slot = origin[0][0]
    old_bond = reactants[slot].GetBondBetweenAtoms(origin[0][1], origin[1][1])
    if (
        old_bond is None
        or old_bond.GetBondType() != Chem.BondType.DOUBLE
        or old_bond.GetStereo() not in DOUBLE
    ):
        result["classification"] = "new_or_template_asserted"
        return result
    old_refs = dict(zip([old_bond.GetBeginAtomIdx(), old_bond.GetEndAtomIdx()], old_bond.GetStereoAtoms()))
    new_refs = dict(zip(ends, bond.GetStereoAtoms()))
    if len(old_refs) != 2 or len(new_refs) != 2:
        return result
    flips, sides = 1, []
    for i, end in enumerate(ends):
        source_index, other_source = origin[i][1], origin[1 - i][1]
        old = neighbor_keys(reactants[slot], source_index, slot, exclude=other_source)
        new = neighbor_keys(product, end, slot, origins, exclude=ends[1 - i])
        aligned, replacement = correspondence(old, new)
        selected = origins.get(new_refs[end])
        if aligned is None or len(old) != 2 or selected not in new:
            result["reason"] = "insufficient_mapped_double_bond_neighbor_correspondence"
            return result
        selected_aligned = aligned[new.index(selected)]
        flip = 1 if selected_aligned == (slot, old_refs[source_index]) else -1
        flips *= flip
        sides.append(
            {
                "source_endpoint": source_index,
                "product_endpoint": end,
                "selected_neighbor_flip": flip,
                "substituent_replacement": replacement,
            }
        )
    retained = DOUBLE[bond.GetStereo()] * flips == DOUBLE[old_bond.GetStereo()]
    result.update(
        classification="retained" if retained else "inverted",
        requires_stereo_prediction=not retained or any(s["substituent_replacement"] for s in sides),
        sides=sides,
        comparison="Relative mapped substituent geometry, not E/Z priority labels",
    )
    return result


def product_stereo(product, reactants, origins):
    state = stereo_state(product)
    changes = []
    for item in state["centers"]:
        if item["specified"] != "Specified":
            changes.append({"type": item["type"], "index": item["index"], "classification": "unresolved"})
        elif item["type"] == "Atom_Tetrahedral":
            changes.append(tetra_audit(product, item["index"], reactants, origins))
        elif item["type"] == "Bond_Double":
            changes.append(double_audit(product, item["index"], reactants, origins))
        else:
            changes.append({"type": item["type"], "index": item["index"], "classification": "unresolved"})
    return {
        "state": state,
        "changes": changes,
        "status": "specified_correspondence"
        if state["fully_specified"] and all(c["classification"] != "unresolved" for c in changes)
        else "unsupported",
        "requires_stereo_prediction": any(c.get("requires_stereo_prediction", True) for c in changes),
        "change_scope": "Product stereocenters; complete reactant stereo states are recorded separately, including centers lost from the principal product",
        "selectivity": "Not experimentally established; new/inverted/template-asserted stereo requires exact external support",
    }


def normalize_product(product):
    """Keep concrete metadata while aligning origins to canonical product indices."""
    mol = Chem.Mol(product)
    for atom in mol.GetAtoms():
        atom.SetIntProp("_audit_raw_index", atom.GetIdx())
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    mol = Chem.RemoveHs(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    exact = canonical(smiles)
    if smiles != exact:
        raise ValueError("Canonical product indexing changed during strict molecule normalization")
    order = ast.literal_eval(mol.GetProp("_smilesAtomOutputOrder"))
    parsed = molecule(exact)
    if len(order) != parsed.GetNumAtoms() or any(
        (mol.GetAtomWithIdx(old).GetAtomicNum(), mol.GetAtomWithIdx(old).GetIsotope())
        != (parsed.GetAtomWithIdx(new).GetAtomicNum(), parsed.GetAtomWithIdx(new).GetIsotope())
        for new, old in enumerate(order)
    ):
        raise ValueError("Canonical product atom correspondence is unavailable")
    return mol, exact, {old: new for new, old in enumerate(order)}


def replay_witnesses(rule, reactants, target, product_limit):
    reaction = rdChemReactions.ReactionFromSmarts(rule.smarts)
    with rdBase.BlockLogs():
        reaction.Initialize()
    reaction.GetSubstructParams().useChirality = True
    maps, template = template_maps(reaction)
    raw = reaction.RunReactants(tuple(reactants), maxProducts=product_limit + 1)
    witnesses, invalid = [], []
    for index, entry in enumerate(raw[:product_limit]):
        try:
            product, exact, canonical_indices = normalize_product(entry[0])
        except (ValueError, RuntimeError) as exc:
            invalid.append({"enumeration_index": index, "reason": str(exc)})
            continue
        if exact != target:
            continue
        origins, atoms = atom_origins(product, reactants, maps)
        for atom in atoms["atoms"]:
            atom["canonical_product_atom_index"] = canonical_indices[atom["product_atom_index"]]
            atom["raw_product_atom_index"] = product.GetAtomWithIdx(atom["product_atom_index"]).GetIntProp(
                "_audit_raw_index"
            )
        stereo = product_stereo(product, reactants, origins)
        witness = {
            "enumeration_index": index,
            "product": exact,
            "atom_origins": atoms,
            "stereo": stereo,
            "canonical_atom_order": [canonical_indices[i] for i in range(product.GetNumAtoms())],
        }
        witness["witness_sha256"] = binding(witness)
        witnesses.append(witness)
    return (
        witnesses,
        template,
        {
            "requested_product_limit": product_limit,
            "raw_entries_observed": len(raw),
            "examined_entries": min(len(raw), product_limit),
            "truncated": len(raw) > product_limit,
            "exact_matching_witnesses": len(witnesses),
            "invalid_entries": invalid,
        },
    )


def audit_route(route, catalog, rules, product_limit=64):
    """Recompute exact route and concrete origin/stereo witnesses from real assets."""
    positive(product_limit, "reaction_product_limit")
    if product_limit > 64:
        raise ValueError("Concrete route audit permits at most 64 reaction products")
    if (
        not isinstance(route, dict)
        or set(route) != {"schema_version", "root", "nodes"}
        or type(route.get("schema_version")) is not int
        or route["schema_version"] != 1
    ):
        raise ValueError("Authoritative route requires only schema_version=1, root and nodes")
    if digest(catalog.path) != catalog.sha256 or digest(rules.path) != rules.sha256:
        raise ValueError("Authoritative catalog/rules files changed since loading")
    catalog._check_import_unchanged()
    verification = verify_route(route, catalog, rules, product_limit)
    order = _topology(route["nodes"], route["root"])
    values, parents, blocks, steps = {}, {}, [], []
    for node in order:
        ident = node["id"]
        if "building_block_id" in node:
            block = catalog.get(node["building_block_id"])
            values[ident] = canonical(block.smiles)
            record = {
                "node_id": ident,
                "smiles": values[ident],
                "catalog_record": catalog.evidence(block.id),
                "stereo": stereo_state(molecule(values[ident])),
            }
            record["building_block_binding_sha256"] = binding(record)
            parents[ident] = record["building_block_binding_sha256"]
            blocks.append(record)
            continue
        slots = [
            {"slot": slot, "node_id": ref, "smiles": values[ref], "parent_binding_sha256": parents[ref]}
            for slot, ref in enumerate(node["reactants"])
        ]
        reactants = [molecule(s["smiles"]) for s in slots]
        product = canonical(node["product"])
        rule = rules.get(node["rule_id"])
        witnesses, template, enumeration = replay_witnesses(rule, reactants, product, product_limit)
        if not witnesses:
            raise ValueError(f"Reaction {ident} lacks a fresh exact concrete replay witness")
        atom_ok = all(w["atom_origins"]["status"] == "complete" for w in witnesses) and not any(
            template.values()
        )
        reactant_stereo = [stereo_state(m) for m in reactants]
        stereo_ok = all(s["fully_specified"] for s in reactant_stereo) and all(
            w["stereo"]["status"] == "specified_correspondence" for w in witnesses
        )
        signatures = {
            binding(
                [
                    {
                        "product_atom_index": a["canonical_product_atom_index"],
                        "reactant_slot": a.get("reactant_slot"),
                        "reactant_atom_index": a.get("reactant_atom_index"),
                    }
                    for a in sorted(
                        w["atom_origins"]["atoms"], key=lambda a: a["canonical_product_atom_index"]
                    )
                ]
            )
            for w in witnesses
        }
        step = {
            "node_id": ident,
            "rule_id": rule.id,
            "smarts": rule.smarts,
            "rule_source": rule.source,
            "rule_smarts_sha256": hashlib.sha256(rule.smarts.encode()).hexdigest(),
            "reactant_slots": slots,
            "product": product,
            "template_audit": template,
            "reactant_stereo": reactant_stereo,
            "enumeration": enumeration,
            "witnesses": witnesses,
            "origin_mapping_ambiguous": len(signatures) > 1,
            "origin_mapping_distinct_witnesses": len(signatures),
            "atom_origin_status": "complete" if atom_ok else "unsupported",
            "stereo_status": "specified_correspondence" if stereo_ok else "unsupported",
            "requires_stereo_prediction": any(w["stereo"]["requires_stereo_prediction"] for w in witnesses),
            "strict_eligible": atom_ok and stereo_ok,
            "catalog_sha256": catalog.sha256,
            "rules_sha256": rules.sha256,
            "product_limit": product_limit,
        }
        step["step_binding_sha256"] = binding(step)
        values[ident], parents[ident] = product, step["step_binding_sha256"]
        steps.append(step)
    result = {
        "schema": SCHEMA,
        "root": route["root"],
        "product": values[route["root"]],
        "verification": verification,
        "catalog_sha256": catalog.sha256,
        "rules_sha256": rules.sha256,
        "rdkit_version": rdBase.rdkitVersion,
        "product_limit": product_limit,
        "enumeration_truncated": any(s["enumeration"]["truncated"] for s in steps),
        "building_blocks": blocks,
        "steps": steps,
        "atom_origin_status": "complete"
        if all(s["atom_origin_status"] == "complete" for s in steps)
        else "unsupported",
        "stereo_status": "specified_correspondence"
        if all(s["stereo_status"] == "specified_correspondence" for s in steps)
        else "unsupported",
        "strict_eligible": all(s["strict_eligible"] for s in steps),
        "requires_stereo_prediction": any(s["requires_stereo_prediction"] for s in steps),
        "limitations": [
            "Bounded concrete computational witnesses, not experimental atom-source or selectivity evidence",
            "Heavy atoms and explicit isotope-H audited; ordinary H/protons and total mass balance unassessed",
            "Lost input atoms and unmodeled byproducts are reported without inventing reagents",
            "Symmetry-equivalent origins and enumeration truncation remain explicit; witness enumeration is not exhaustive beyond the cap",
            "Stepwise canonical atom bindings permit lineage tracing; shared DAG inputs in different slots are separate equivalents",
            "RDKit must expose concrete react_idx/react_atom_idx metadata; runtimes lacking slot origins fail closed",
            "Conditions, kinetics, yield, practical applicability and empirical validation are unassessed",
        ],
    }
    result["route_binding_sha256"] = binding(result)
    return result
