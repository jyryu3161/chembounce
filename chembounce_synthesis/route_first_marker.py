"""Bounded set-valued departing-site marker projection; never edits a reactant."""

import re
import time

from rdkit import Chem, DataStructs

from chembounce_synthesis.common import FP
from chembounce_synthesis.route_audit import template_maps

HALOGENS = frozenset((9, 17, 35, 53))
MAX_MATCHES = 65  # The 65th match is a truncation sentinel.


class UnknownMarker(ValueError):
    pass


class Meter:
    def __init__(
        self, deadline, *, fp_limit=20, match_limit=20, parse_limit=20,
        score_limit=20, units_limit=8192,
    ):
        self.deadline = deadline
        self.fp_limit = fp_limit
        self.match_limit = match_limit
        self.parse_limit = parse_limit
        self.score_limit = score_limit
        self.units_limit = units_limit
        self.fingerprints = 0
        self.match_requests = 0
        self.parses = 0
        self.explicit_parses = 0
        self.internal_canonical_parse_reservations = 0
        self.tanimoto_attempts = 0
        self.projected_serialization_attempts = 0
        self.structural_units = 0

    def check(self):
        if time.monotonic() >= self.deadline:
            raise UnknownMarker("deadline")

    def charge_structure(self, *mols):
        self.check()
        units = 3 * sum(mol.GetNumAtoms() + mol.GetNumBonds() for mol in mols)
        if any(mol.GetNumAtoms() > 256 or mol.GetNumBonds() > 512 for mol in mols):
            raise UnknownMarker("molecule_size_cap")
        if self.structural_units + units > self.units_limit:
            raise UnknownMarker("structural_unit_cap")
        self.structural_units += units

    def parse(self, smiles):
        self.check()
        if self.parses >= self.parse_limit:
            raise UnknownMarker("parse_cap")
        self.parses += 1
        self.explicit_parses += 1
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise UnknownMarker("parse_failed")
        self.charge_structure(mol)
        return mol

    def reserve_internal_canonical_parse(self):
        """Conservatively reserve a parse for typed_reverse_hint's canonical()."""
        self.check()
        if self.parses >= self.parse_limit:
            raise UnknownMarker("parse_cap")
        self.parses += 1
        self.internal_canonical_parse_reservations += 1

    def fingerprint(self, mol):
        self.check()
        if self.fingerprints >= self.fp_limit:
            raise UnknownMarker("fingerprint_cap")
        self.fingerprints += 1  # Charge the attempted generator call before it runs.
        try:
            fp = FP.GetFingerprint(mol)
        except (ValueError, RuntimeError) as exc:
            raise UnknownMarker("projected_fingerprint_failed") from exc
        self.check()
        return fp


def _halogen_query(atom):
    if not atom.HasQuery():
        return atom.GetAtomicNum() in HALOGENS
    lines = atom.DescribeQuery().splitlines()
    roots, stack = [], []
    for line in lines:
        match = re.fullmatch(r"( *)(.+)", line)
        if match is None or len(match[1]) % 2:
            return False
        depth = len(match[1]) // 2
        if depth > len(stack):
            return False
        node = [match[2], []]
        stack = stack[:depth]
        (stack[-1][1] if stack else roots).append(node)
        stack.append(node)

    def proved(node):
        label, children = node
        if label == "AtomOr" and children:
            return all(proved(child) for child in children)
        if not children:
            item = re.fullmatch(r"AtomType (\d+) = val", label)
            return item is not None and int(item[1]) in HALOGENS
        return False

    return len(roots) == 1 and proved(roots[0])


def _fixed_element(atom):
    """Prove an endpoint's element from a small supported Boolean query tree."""
    if not atom.HasQuery():
        return atom.GetAtomicNum() or None
    roots, stack = [], []
    for line in atom.DescribeQuery().splitlines():
        match = re.fullmatch(r"( *)(.+)", line)
        if match is None or len(match[1]) % 2:
            return None
        depth = len(match[1]) // 2
        if depth > len(stack):
            return None
        node = [match[2], []]
        stack = stack[:depth]
        (stack[-1][1] if stack else roots).append(node)
        stack.append(node)

    def fixed(node):
        label, children = node
        if label in ("AtomAnd", "AtomOr") and children:
            values = [fixed(child) for child in children]
            proven = {value for value in values if value is not None}
            if label == "AtomOr" and len(proven) != 1:
                return None
            return next(iter(proven)) if len(proven) == 1 and (
                label == "AtomAnd" or all(value is not None for value in values)
            ) else None
        if children:
            return None
        match = re.fullmatch(r"Atom(?:AtomicNum|Type) (\d+) = val", label)
        if match is None:
            return None
        value = int(match[1])
        return value - 1000 if value in (1006, 1007) else value if value < 1000 else None

    return fixed(roots[0]) if len(roots) == 1 else None


def departure_spec(reaction, slot):
    """Return (template atom, mapped attachment), or remain unknown."""
    if not 0 <= slot < reaction.GetNumReactantTemplates() or reaction.GetNumReactantTemplates() != 2:
        raise UnknownMarker("unsupported_arity")
    if reaction.GetNumProductTemplates() != 1:
        raise UnknownMarker("unsupported_products")
    maps, issues = template_maps(reaction)
    if any(issues.values()):
        raise UnknownMarker("ambiguous_maps")
    product = reaction.GetProducts()[0]
    retained = {atom.GetAtomMapNum() for atom in product.GetAtoms() if atom.GetAtomMapNum()}
    role = reaction.GetReactants()[slot]
    departing = [atom for atom in role.GetAtoms() if not atom.GetAtomMapNum()]
    if len(departing) != 1 or not _halogen_query(departing[0]):
        raise UnknownMarker("unproved_single_halogen_template")
    atom = departing[0]
    if atom.GetDegree() != 1:
        raise UnknownMarker("nonterminal_template_atom")
    bond = atom.GetBonds()[0]
    center = bond.GetOtherAtom(atom)
    if (
        center.GetAtomMapNum() not in retained
        or maps.get(center.GetAtomMapNum()) != (slot, center.GetIdx())
        or bond.GetBondType() != Chem.BondType.SINGLE
        or bond.GetIsAromatic()
        or bond.GetStereo() != Chem.BondStereo.STEREONONE
        or bond.GetBondDir() != Chem.BondDir.NONE
    ):
        raise UnknownMarker("unproved_departing_attachment")
    cross = [
        bond for bond in product.GetBonds()
        if (a := bond.GetBeginAtom().GetAtomMapNum())
        and (b := bond.GetEndAtom().GetAtomMapNum())
        and maps[a][0] != maps[b][0]
    ]
    if len(cross) != 1 or center.GetAtomMapNum() not in {
        cross[0].GetBeginAtom().GetAtomMapNum(), cross[0].GetEndAtom().GetAtomMapNum()
    }:
        raise UnknownMarker("unproved_single_cross_bond")
    for end in (cross[0].GetBeginAtom(), cross[0].GetEndAtom()):
        source_slot, source_index = maps[end.GetAtomMapNum()]
        source_atom = reaction.GetReactants()[source_slot].GetAtomWithIdx(source_index)
        source_element, product_element = _fixed_element(source_atom), _fixed_element(end)
        if source_element is None or source_element != product_element:
            raise UnknownMarker("endpoint_element_not_proven_same")
    return atom.GetIdx(), center.GetIdx()


def _state(mol):
    atoms = tuple(
        (
            atom.GetAtomicNum(), atom.GetIsotope(), atom.GetFormalCharge(),
            atom.GetNumExplicitHs(), atom.GetNumImplicitHs(), atom.GetNoImplicit(),
            atom.GetIsAromatic(), atom.GetNumRadicalElectrons(), int(atom.GetChiralTag()),
            atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None,
        )
        for atom in mol.GetAtoms()
    )
    bonds = tuple(
        (
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), str(bond.GetBondType()),
            bond.GetIsAromatic(), int(bond.GetStereo()), int(bond.GetBondDir()),
            tuple(bond.GetStereoAtoms()),
        )
        for bond in mol.GetBonds()
    )
    return atoms, bonds


def _site(mol, leaving, attachment):
    if not 0 <= leaving < mol.GetNumAtoms() or not 0 <= attachment < mol.GetNumAtoms():
        raise UnknownMarker("site_out_of_range")
    atom, center = mol.GetAtomWithIdx(leaving), mol.GetAtomWithIdx(attachment)
    bond = mol.GetBondBetweenAtoms(leaving, attachment)
    if (
        atom.GetAtomicNum() not in HALOGENS
        or atom.GetDegree() != 1
        or atom.GetIsotope()
        or atom.GetFormalCharge()
        or atom.GetNumExplicitHs()
        or atom.GetNumImplicitHs()
        or atom.GetNumRadicalElectrons()
        or atom.GetIsAromatic()
        or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        or bond is None
        or bond.GetBondType() != Chem.BondType.SINGLE
        or bond.GetIsAromatic()
        or center.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    ):
        raise UnknownMarker("unsafe_concrete_departure")
    if any(
        edge.GetStereo() != Chem.BondStereo.STEREONONE
        or edge.GetBondDir() != Chem.BondDir.NONE
        or leaving in edge.GetStereoAtoms()
        for edge in center.GetBonds()
    ):
        raise UnknownMarker("adjacent_specified_stereo")


def marked_view(mol, leaving, attachment, meter):
    meter.charge_structure(mol)
    _site(mol, leaving, attachment)
    before_atoms, before_bonds = _state(mol)
    view = Chem.Mol(mol)
    view.GetAtomWithIdx(leaving).SetAtomicNum(0)
    if Chem.SanitizeMol(view, catchErrors=True) != Chem.SanitizeFlags.SANITIZE_NONE:
        raise UnknownMarker("marker_sanitization_failed")
    after_atoms, after_bonds = _state(view)
    if before_bonds != after_bonds or any(
        before != after
        for index, (before, after) in enumerate(zip(before_atoms, after_atoms))
        if index != leaving
    ):
        raise UnknownMarker("marker_changed_other_atoms_or_bonds")
    if (0,) + before_atoms[leaving][1:] != after_atoms[leaving]:
        raise UnknownMarker("marker_changed_atom_state")
    meter.check()
    return view


def source_view(mol, reaction, slot, meter):
    leaving_query, center_query = departure_spec(reaction, slot)
    meter.charge_structure(mol)
    if meter.match_requests >= meter.match_limit:
        raise UnknownMarker("match_request_cap")
    meter.match_requests += 1
    matches = mol.GetSubstructMatches(
        reaction.GetReactants()[slot], useChirality=True, uniquify=False, maxMatches=MAX_MATCHES
    )
    if not matches:
        raise UnknownMarker("source_role_ineligible")
    if len(matches) >= MAX_MATCHES:
        raise UnknownMarker("65_match_sentinel")
    sites = {(row[leaving_query], row[center_query]) for row in matches}
    if len(sites) != 1:
        raise UnknownMarker("distinct_departing_sites")
    meter.check()
    return marked_view(mol, *next(iter(sites)), meter)


def source_views(mol, reaction, slot, meter):
    """Return every safe site and view; a single bad or missing view rejects all."""
    leaving_query, center_query = departure_spec(reaction, slot)
    meter.charge_structure(mol)
    if meter.match_requests >= meter.match_limit:
        raise UnknownMarker("match_request_cap")
    meter.match_requests += 1
    matches = mol.GetSubstructMatches(
        reaction.GetReactants()[slot], useChirality=True, uniquify=False, maxMatches=MAX_MATCHES
    )
    if not matches:
        raise UnknownMarker("source_role_ineligible")
    if len(matches) >= MAX_MATCHES:
        raise UnknownMarker("65_match_sentinel")
    sites = sorted({(row[leaving_query], row[center_query]) for row in matches})
    if len(sites) > 64:
        raise UnknownMarker("site_cap")
    views = {}
    associations = []
    for leaving, attachment in sites:
        meter.check()
        marked = marked_view(mol, leaving, attachment, meter)
        meter.projected_serialization_attempts += 1
        try:
            smiles = Chem.MolToSmiles(marked, isomericSmiles=True)
        except (ValueError, RuntimeError) as exc:
            raise UnknownMarker("projected_serialization_failed") from exc
        if smiles not in views:
            views[smiles] = meter.fingerprint(marked)
        associations.append({
            "departing_atom_index": leaving,
            "attachment_atom_index": attachment,
            "projected_smiles": smiles,
        })
    meter.check()
    return associations, views


def max_marked_similarity(guide_fp, views, meter):
    """One fixed optimistic score; it never chooses a reaction site."""
    if guide_fp is None or not views:
        raise UnknownMarker("missing_complete_projected_views")
    scores = []
    for fp in views.values():
        meter.check()
        if meter.tanimoto_attempts >= meter.score_limit:
            raise UnknownMarker("tanimoto_cap")
        meter.tanimoto_attempts += 1
        scores.append(DataStructs.TanimotoSimilarity(guide_fp, fp))
        meter.check()
    return max(scores)


def guide_view(raw, clean, reaction, slot, tags, meter):
    leaving_query, center_query = departure_spec(reaction, slot)
    meter.charge_structure(raw, clean)
    if raw.GetNumAtoms() != clean.GetNumAtoms() or raw.GetNumBonds() != clean.GetNumBonds():
        raise UnknownMarker("raw_clean_graph_changed")
    tag = tags.get((slot, leaving_query))
    if tag is None:
        raise UnknownMarker("template_departure_tag_missing")
    tagged = [
        atom for atom in raw.GetAtoms()
        if atom.HasProp("old_mapno") and atom.GetIntProp("old_mapno") == tag
    ]
    if len(tagged) != 1 or tagged[0].GetAtomicNum() != 0:
        raise UnknownMarker("unproved_originless_representative")
    atom = tagged[0]
    if atom.HasProp("react_idx") or atom.HasProp("react_atom_idx"):
        raise UnknownMarker("representative_has_parent_origin")
    if atom.GetDegree() != 1:
        raise UnknownMarker("nonterminal_raw_representative")
    neighbor = atom.GetNeighbors()[0]
    if not neighbor.HasProp("old_mapno") or neighbor.GetIntProp("old_mapno") != (
        reaction.GetReactants()[slot].GetAtomWithIdx(center_query).GetAtomMapNum()
    ):
        raise UnknownMarker("representative_wrong_attachment")
    for index in range(raw.GetNumAtoms()):
        if index == atom.GetIdx():
            continue
        a, b = raw.GetAtomWithIdx(index), clean.GetAtomWithIdx(index)
        if (a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge(), int(a.GetChiralTag())) != (
            b.GetAtomicNum(), b.GetIsotope(), b.GetFormalCharge(), int(b.GetChiralTag())
        ):
            raise UnknownMarker("raw_clean_other_atom_changed")
    for bond in raw.GetBonds():
        same = clean.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if same is None or (
            str(bond.GetBondType()), bond.GetIsAromatic(), int(bond.GetStereo()), int(bond.GetBondDir()),
            tuple(bond.GetStereoAtoms()),
        ) != (
            str(same.GetBondType()), same.GetIsAromatic(), int(same.GetStereo()), int(same.GetBondDir()),
            tuple(same.GetStereoAtoms()),
        ):
            raise UnknownMarker("raw_clean_bond_changed")
    meter.check()
    if meter.match_requests >= meter.match_limit:
        raise UnknownMarker("match_request_cap")
    meter.match_requests += 1
    matches = clean.GetSubstructMatches(
        reaction.GetReactants()[slot], useChirality=True, uniquify=False, maxMatches=MAX_MATCHES
    )
    if len(matches) >= MAX_MATCHES:
        raise UnknownMarker("65_match_sentinel")
    if not any(
        row[leaving_query] == atom.GetIdx() and row[center_query] == neighbor.GetIdx()
        for row in matches
    ):
        raise UnknownMarker("tagged_site_not_role_eligible")
    return marked_view(clean, atom.GetIdx(), neighbor.GetIdx(), meter)
