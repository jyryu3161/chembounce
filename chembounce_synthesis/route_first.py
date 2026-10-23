"""Opt-in, bounded route-first retrieval. Reverse templates are hints, never evidence."""

import hashlib
import heapq
import itertools
import math
import re
import time
from dataclasses import asdict, dataclass

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdChemReactions

from .common import FP, absolute_smiles, canonical, molecule, positive
from .route_audit import audit_route, stereo_state, template_maps
from .routes import verify_route

VERSION = "catalog-rule-route-first-v2"
REVERSE_ROLE_REPRESENTATIVES = (
    (8, 0, 1, "OH"),
    (9, 0, 0, "F"),
    (17, 0, 0, "Cl"),
    (35, 0, 0, "Br"),
    (53, 0, 0, "I"),
    (8, -1, 0, "O-"),
)
MAX_REVERSE_ROLE_TRIALS_PER_OUTCOME = 16
MAX_REVERSE_ROLE_WORK_PER_SEARCH = 1024
MAX_REVERSE_ROLE_MATCHES = 65
MAX_OUTER_CENTER_METADATA_CHECKS = 32768
MAX_REPRESENTATIVE_ORIGIN_VISITS = 32768
MAX_REPRESENTATIVE_ORIGIN_GUIDES = 128
MAX_REPRESENTATIVE_ORIGIN_ATOMS = 256
MAX_REPRESENTATIVE_ORIGIN_BONDS = 512
MAX_TOPOLOGY_MATCH_REQUESTS = 128
MAX_TOPOLOGY_MATCHES = 65  # 64 plus a truncation sentinel.
MAX_TOPOLOGY_SPEC_GRAPH_VISITS = 262144
MAX_CATALOG_STEREO_PERCEPTIONS = 2048
MAX_CATALOG_STEREO_ATOM_BOND_CHARGES = 131072
MAX_CATALOG_STEREO_CACHE_ENTRIES = 2048
MAX_CATALOG_STEREO_EXCLUSION_SAMPLES = 16
MAX_POOL_QUERY_ROLE_PAIRS = 128
MAX_POOL_WINDOW_MATCH_REQUESTS = 128
MAX_POOL_SITE_MATCH_REQUESTS = 384
MAX_POOL_MATCHES = 65
MAX_POOL_QUERY_ATOMS = 256
MAX_POOL_QUERY_BONDS = 512
MAX_POOL_GRAPH_VISITS = 262144
MAX_ROLE_INDEX_MATCH_CHECKS = 500000
MAX_FALLBACK_HIT_SECOND_TANIMOTOS = 65536
MAX_FALLBACK_HIT_SECOND_SAMPLES = 16
MAX_CONTINUATION_ROLE_MATCHES = 65  # 64 plus a truncation sentinel.
MAX_CONTINUATION_ROLE_ATOMS = 256
MAX_CONTINUATION_ROLE_BONDS = 512
MAX_CONTINUATION_PRECURSORS = 4
MAX_CONTINUATION_ORIGIN_VISITS = 8192
MAX_CONTINUATION_RAW_SNAPSHOTS = 1024
MAX_CONTINUATION_SNAPSHOT_UNITS = 262144
MAX_CONTINUATION_ROLE_CONTEXTS = 128
MAX_CONTINUATION_REGISTRIES = 512
MAX_CONTINUATION_RANKED_ROWS = 65536
MAX_CONTINUATION_ROLE_CHECKS_PER_QUERY = 8192
MAX_CONTINUATION_ROLE_CHECKS_PER_SEARCH = 65536
MAX_CONTINUATION_OUTER_ROLE_QUERIES = 16
MAX_CONTINUATION_SELECTED_CONTEXTS = 128


def _pool_raw_counts(mol, atoms=None):
    if atoms is None:
        atoms = set(range(mol.GetNumAtoms()))
    return (
        sum(mol.GetAtomWithIdx(i).GetAtomicNum() != 1 for i in atoms),
        sum(mol.GetAtomWithIdx(i).GetAtomicNum() not in (1, 6) for i in atoms),
        sum(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms),
        sum(set(ring) <= atoms for ring in mol.GetRingInfo().AtomRings()),
    )


def _pool_bins(counts):
    return counts[0] // 4, counts[1] // 2, counts[2] // 4, min(counts[3], 8)


def _normalized_role_smarts(role):
    query = Chem.MolFromSmarts(role)
    if query is None:
        raise ValueError("Invalid catalog retrieval role SMARTS")
    for atom in query.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmarts(query)


def _pool_ring_systems(mol):
    systems = []
    for ring in (set(item) for item in mol.GetRingInfo().AtomRings()):
        merging = [index for index, system in enumerate(systems) if system & ring]
        for index in reversed(merging):
            ring.update(systems.pop(index))
        systems.append(ring)
    changed = True
    while changed:
        changed = False
        for first in range(len(systems)):
            for second in range(first + 1, len(systems)):
                if systems[first] & systems[second]:
                    systems[first].update(systems.pop(second))
                    changed = True
                    break
            if changed:
                break
    return systems


def topology_rule_specification(reaction):
    """Describe a supported product cut without asserting chemical feasibility."""
    if len(reaction.GetReactants()) != 2 or len(reaction.GetProducts()) != 1:
        return None, "not_two_reactants_one_product"
    maps, issues = template_maps(reaction)
    if any(issues.values()):
        return None, "ambiguous_or_unbound_maps"
    product = reaction.GetProducts()[0]
    retained = {atom.GetAtomMapNum() for atom in product.GetAtoms() if atom.GetAtomMapNum()}
    product_edges = {
        frozenset((first, second))
        for bond in product.GetBonds()
        if (first := bond.GetBeginAtom().GetAtomMapNum()) and (second := bond.GetEndAtom().GetAtomMapNum())
    }
    for reactant in reaction.GetReactants():
        for bond in reactant.GetBonds():
            first, second = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if first in retained and second in retained and frozenset((first, second)) not in product_edges:
                return None, "retained_reactant_edge_deleted"
    cross = []
    for bond in product.GetBonds():
        first, second = bond.GetBeginAtom(), bond.GetEndAtom()
        a, b = first.GetAtomMapNum(), second.GetAtomMapNum()
        if a and b and maps[a][0] != maps[b][0]:
            cross.append((first.GetIdx(), second.GetIdx()))
    if not cross:
        return None, "no_cross_slot_product_bond"
    return (product, maps, tuple(cross)), None


def reverse_partition_issue(outcome, target_atoms):
    """Require each target atom to belong to exactly one reverse precursor."""
    origins = set()
    for precursor in outcome:
        for atom in precursor.GetAtoms():
            has_index = atom.HasProp("react_atom_idx")
            has_slot = atom.HasProp("react_idx")
            if not has_index and not has_slot:
                # A reverse template may introduce a leaving group as a hint.
                continue
            if not has_index or not has_slot:
                return "incomplete_origin_metadata"
            if atom.GetIntProp("react_idx") != 0:
                return "wrong_parent_slot"
            index = atom.GetIntProp("react_atom_idx")
            if not 0 <= index < target_atoms:
                return "origin_out_of_range"
            if index in origins:
                return "reused_parent_atom"
            origins.add(index)
    if len(origins) != target_atoms:
        return "missing_parent_atoms"
    return None


def continuation_role_origin(parent, outcome, role, *, deadline_monotonic=None, max_origin_visits=None):
    """Describe an outer-role embedding in one raw inner reverse outcome.

    This is diagnostic only. Local atom/bond identity is not proof that an
    arbitrary SMARTS predicate, or a future concrete catalog leaf, preserves
    the role. Call before reverse-hint normalization/SMILES serialization.
    """
    result = {
        "classification": "unknown",
        "reason": None,
        "source_slot": None,
        "match_count": None,
        "match_limit": MAX_CONTINUATION_ROLE_MATCHES,
        "origin_visits": 0,
        "matched_parent_atoms": None,
        "owner_slots": None,
        "atom_changes": [],
        "bond_changes": [],
        "full_smarts_preservation_proven": False,
        "concrete_leaf_requirement_proven": False,
    }

    def unknown(reason):
        result["reason"] = reason
        return result

    if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
        return unknown("deadline")
    if (
        parent.GetNumAtoms() > MAX_CONTINUATION_ROLE_ATOMS
        or parent.GetNumBonds() > MAX_CONTINUATION_ROLE_BONDS
        or role.GetNumAtoms() > MAX_CONTINUATION_ROLE_ATOMS
        or role.GetNumBonds() > MAX_CONTINUATION_ROLE_BONDS
    ):
        return unknown("parent_size_limit")
    if not outcome or len(outcome) > MAX_CONTINUATION_PRECURSORS:
        return unknown("precursor_size_limit")
    if any(
        part.GetNumAtoms() > MAX_CONTINUATION_ROLE_ATOMS or part.GetNumBonds() > MAX_CONTINUATION_ROLE_BONDS
        for part in outcome
    ):
        return unknown("precursor_size_limit")

    def charge(amount=1):
        if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
            return False, "deadline"
        ceiling = (
            MAX_CONTINUATION_ORIGIN_VISITS
            if max_origin_visits is None
            else min(MAX_CONTINUATION_ORIGIN_VISITS, max_origin_visits)
        )
        if result["origin_visits"] + amount > ceiling:
            return False, "origin_visit_limit"
        result["origin_visits"] += amount
        return True, None

    try:
        ok, reason = charge(sum(part.GetNumAtoms() for part in outcome))
        if not ok:
            return unknown(reason)
        issue = reverse_partition_issue(outcome, parent.GetNumAtoms())
        if issue is not None:
            return unknown(f"partition_{issue}")
        ok, reason = charge(parent.GetNumAtoms() + parent.GetNumBonds() + role.GetNumAtoms())
        if not ok:
            return unknown(reason)
        matches = parent.GetSubstructMatches(
            role,
            useChirality=True,
            uniquify=False,
            maxMatches=MAX_CONTINUATION_ROLE_MATCHES,
        )
        result["match_count"] = len(matches)
        if len(matches) >= MAX_CONTINUATION_ROLE_MATCHES:
            return unknown("role_match_truncated")
        if not matches:
            return unknown("role_not_in_parent")
        if len(matches) != 1:
            return unknown("ambiguous_role_embeddings")
        match = matches[0]
        result["matched_parent_atoms"] = list(match)
        stereo_sensitive = any(
            parent.GetAtomWithIdx(index).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            or any(
                bond.GetStereo() != Chem.BondStereo.STEREONONE or bond.GetBondDir() != Chem.BondDir.NONE
                for bond in parent.GetAtomWithIdx(index).GetBonds()
            )
            for index in match
        )
        owners = {}
        for slot, precursor in enumerate(outcome):
            ok, reason = charge(precursor.GetNumAtoms() + precursor.GetNumBonds())
            if not ok:
                return unknown(reason)
            for atom in precursor.GetAtoms():
                if atom.HasProp("react_atom_idx"):
                    owners[atom.GetIntProp("react_atom_idx")] = (slot, atom.GetIdx())
        slots = sorted({owners[index][0] for index in match})
        result["owner_slots"] = slots
        if len(slots) == 1:
            result["source_slot"] = slots[0]

        def atom_state(atom):
            return (
                atom.GetAtomicNum(),
                atom.GetIsotope(),
                atom.GetFormalCharge(),
                atom.GetIsAromatic(),
                atom.GetNumExplicitHs(),
                atom.GetNoImplicit(),
                int(atom.GetChiralTag()),
                atom.GetNumRadicalElectrons(),
            )

        def incident(mol, atom, origin):
            signatures = []
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                other = origin(neighbor)
                signatures.append(
                    (
                        other,
                        str(bond.GetBondType()),
                        bond.GetIsAromatic(),
                        int(bond.GetStereo()),
                        int(bond.GetBondDir()),
                    )
                )
            return sorted(signatures, key=repr)

        for index in match:
            slot, raw_index = owners[index]
            precursor = outcome[slot]
            ok, reason = charge(
                parent.GetAtomWithIdx(index).GetDegree() + precursor.GetAtomWithIdx(raw_index).GetDegree() + 1
            )
            if not ok:
                return unknown(reason)
            original = parent.GetAtomWithIdx(index)
            raw = precursor.GetAtomWithIdx(raw_index)
            if atom_state(original) != atom_state(raw):
                result["atom_changes"].append(index)
            before = incident(parent, original, lambda atom: atom.GetIdx())
            after = incident(
                precursor,
                raw,
                lambda atom: (
                    atom.GetIntProp("react_atom_idx")
                    if atom.HasProp("react_atom_idx") and atom.HasProp("react_idx")
                    else None
                ),
            )
            if before != after:
                result["bond_changes"].append(index)
        if stereo_sensitive:
            return unknown("specified_stereo_mapping_unproven")
        if len(slots) != 1:
            result["classification"] = "split_across_sources"
        elif result["atom_changes"] or result["bond_changes"]:
            result["classification"] = "single_source_edited"
        else:
            result["classification"] = "single_source_local_intact"
        return result
    except (ValueError, RuntimeError, KeyError, AttributeError):
        return unknown("origin_or_match_error")


def _normalized_origin_site_unchanged(raw, clean, target_indices):
    """Require the same target-derived local site through hint normalization."""

    def origin(atom):
        if atom.HasProp("react_atom_idx") and atom.HasProp("react_idx"):
            return atom.GetIntProp("react_idx"), atom.GetIntProp("react_atom_idx")
        return None

    def state(atom):
        return (
            atom.GetAtomicNum(),
            atom.GetIsotope(),
            atom.GetFormalCharge(),
            atom.GetIsAromatic(),
            atom.GetNumExplicitHs(),
            atom.GetNoImplicit(),
            int(atom.GetChiralTag()),
            atom.GetNumRadicalElectrons(),
        )

    def around(atom):
        return sorted(
            [
                (
                    origin(neighbor := bond.GetOtherAtom(atom)),
                    state(neighbor),
                    str(bond.GetBondType()),
                    bond.GetIsAromatic(),
                    int(bond.GetStereo()),
                    int(bond.GetBondDir()),
                )
                for bond in atom.GetBonds()
            ],
            key=repr,
        )

    for target_index in target_indices:
        key = (0, target_index)
        before = [atom for atom in raw.GetAtoms() if origin(atom) == key]
        after = [atom for atom in clean.GetAtoms() if origin(atom) == key]
        if len(before) != 1 or len(after) != 1:
            return False
        if state(before[0]) != state(after[0]) or around(before[0]) != around(after[0]):
            return False
    return True


class ReverseHintError(ValueError):
    def __init__(self, reason):
        self.reason = reason
        super().__init__(reason)


def _role_forces_aliphatic_element(atom, atomic_number):
    """Prove the positive aliphatic atom type through RDKit's Boolean query tree."""
    if not atom.HasQuery():
        return atom.GetAtomicNum() == atomic_number and not atom.GetIsAromatic()
    roots, stack = [], []
    for line in atom.DescribeQuery().splitlines():
        match = re.fullmatch(r"( *)(.+)", line)
        if match is None or len(match[1]) % 2:
            return False
        depth = len(match[1]) // 2
        if depth > len(stack):
            return False
        node = [match[2], []]
        stack = stack[:depth]
        if stack:
            stack[-1][1].append(node)
        else:
            roots.append(node)
        stack.append(node)
    if len(roots) != 1:
        return False

    def proves(node):
        label, children = node
        if label == f"AtomType {atomic_number} = val" and not children:
            return True
        if label == "AtomAnd" and children:
            return any(proves(child) for child in children)
        if label == "AtomOr" and children:
            return all(proves(child) for child in children)
        return False

    return proves(roots[0])


def _restore_reactant_aromatic_query_flags(clean, parent, rule):
    """Drop stale product aromatic atom flags only where a mapped role proves them stale."""
    if not any(atom.GetIsAromatic() for atom in clean.GetAtoms()):
        return 0
    _, issues = template_maps(rule.reaction)
    if any(issues.values()):
        return 0
    reactant_atoms, product_atoms = {}, {}
    for template in rule.reaction.GetReactants():
        for atom in template.GetAtoms():
            if atom.GetAtomMapNum():
                reactant_atoms.setdefault(atom.GetAtomMapNum(), []).append(atom)
    for template in rule.reaction.GetProducts():
        for atom in template.GetAtoms():
            if atom.GetAtomMapNum():
                product_atoms.setdefault(atom.GetAtomMapNum(), []).append(atom)
    # Ring perception is read from this private copy, before sanitization can
    # mutate aromatic flags or bond types.
    Chem.GetSymmSSSR(clean)
    changes = 0
    for atom in clean.GetAtoms():
        if (
            not atom.GetIsAromatic()
            or atom.IsInRing()
            or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            or not atom.HasProp("old_mapno")
            or not atom.HasProp("react_idx")
            or not atom.HasProp("react_atom_idx")
            or atom.GetIntProp("react_idx") != 0
        ):
            continue
        number = atom.GetIntProp("old_mapno")
        source_index = atom.GetIntProp("react_atom_idx")
        role_atoms = reactant_atoms.get(number, ())
        product = product_atoms.get(number, ())
        if (
            len(role_atoms) != 1
            or len(product) != 1
            or source_index < 0
            or source_index >= parent.GetNumAtoms()
            or not _role_forces_aliphatic_element(role_atoms[0], atom.GetAtomicNum())
            or not product[0].GetIsAromatic()
            or product[0].GetAtomicNum() != atom.GetAtomicNum()
            or any(
                bond.GetIsAromatic()
                or bond.GetBondType() == Chem.BondType.AROMATIC
                or bond.GetStereo() != Chem.BondStereo.STEREONONE
                or bond.GetBondDir() != Chem.BondDir.NONE
                for bond in atom.GetBonds()
            )
        ):
            continue
        source = parent.GetAtomWithIdx(source_index)
        if (
            not source.GetIsAromatic()
            or source.GetAtomicNum() != atom.GetAtomicNum()
            or source.GetIsotope() != atom.GetIsotope()
            or source.GetChiralTag() != atom.GetChiralTag()
        ):
            continue
        atom.SetIsAromatic(False)
        changes += 1
    return changes


def _reverse_hint_invariants(clean):
    atoms = tuple(
        (
            atom.GetAtomicNum(),
            atom.GetIsotope(),
            atom.GetFormalCharge(),
            atom.GetChiralTag(),
            atom.GetIsAromatic(),
            atom.GetNumExplicitHs(),
            atom.GetNoImplicit(),
        )
        for atom in clean.GetAtoms()
    )
    bonds = tuple(
        (
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond.GetBondType(),
            bond.GetIsAromatic(),
            bond.GetStereo(),
            bond.GetBondDir(),
            tuple(bond.GetStereoAtoms()),
        )
        for bond in clean.GetBonds()
    )
    return atoms, bonds


def _sanitize_reverse_hint(clean, parent, rule, *, aromatic_changes=None):
    """Sanitize one copied query without changing its target-derived atom identity."""
    if clean.GetNumAtoms() == 0 or any(atom.GetAtomicNum() == 0 for atom in clean.GetAtoms()):
        raise ReverseHintError("remaining_dummy")
    if len(Chem.GetMolFrags(clean)) != 1:
        raise ReverseHintError("disconnected_hint")
    changes = _restore_reactant_aromatic_query_flags(clean, parent, rule)
    if aromatic_changes is not None:
        aromatic_changes["reset_product_aromatic_query_flags"] = changes
    preserved = _reverse_hint_invariants(clean) if changes else None
    with rdBase.BlockLogs():
        issue = Chem.SanitizeMol(clean, catchErrors=True)
    resets = 0
    if issue != Chem.SanitizeFlags.SANITIZE_NONE:
        if issue != Chem.SanitizeFlags.SANITIZE_PROPERTIES:
            raise ReverseHintError("hint_sanitization_failed")
        if preserved is not None and _reverse_hint_invariants(clean) != preserved:
            raise ReverseHintError("aromatic_query_cleanup_changed_identity_or_bonds")
        role_atoms = {
            atom.GetAtomMapNum(): atom
            for role in rule.reaction.GetReactants()
            for atom in role.GetAtoms()
            if atom.GetAtomMapNum()
        }
        for atom in clean.GetAtoms():
            if (
                atom.GetAtomicNum() != 7
                or atom.GetIsAromatic()
                or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                or atom.GetFormalCharge() != 0
                or atom.GetNumExplicitHs() == 0
                or atom.GetDegree() + atom.GetNumExplicitHs() <= 3
                or not atom.HasProp("react_atom_idx")
                or not atom.HasProp("old_mapno")
            ):
                continue
            source = parent.GetAtomWithIdx(atom.GetIntProp("react_atom_idx"))
            role_atom = role_atoms.get(atom.GetIntProp("old_mapno"))
            if (
                role_atom is None
                or role_atom.GetAtomicNum() != 7
                or role_atom.GetNumExplicitHs() != 0
                or source.GetAtomicNum() != 7
                or source.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                or source.GetIsotope() != atom.GetIsotope()
                or source.GetFormalCharge() != atom.GetFormalCharge()
                or source.GetNumExplicitHs() != 0
                or source.GetDegree() <= atom.GetDegree()
            ):
                continue
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(False)
            resets += 1
        if preserved is not None and resets:
            # The existing N-H repair is separately guarded by parent and
            # reactant-template provenance above. Everything else remains fixed.
            preserved = _reverse_hint_invariants(clean)
        with rdBase.BlockLogs():
            retry_issue = Chem.SanitizeMol(clean, catchErrors=True)
        if not resets or retry_issue != Chem.SanitizeFlags.SANITIZE_NONE:
            raise ReverseHintError("hint_sanitization_failed")
    if preserved is not None and _reverse_hint_invariants(clean) != preserved:
        raise ReverseHintError("aromatic_query_cleanup_changed_identity_or_bonds")
    return resets


def normalize_reverse_hint(raw, parent, rule, *, capture_clean=None):
    """Make a copied reverse outcome usable for retrieval, never route evidence."""
    editable = Chem.RWMol(Chem.Mol(raw))
    edits = {
        "removed_originless_terminal_dummy": 0,
        "reset_template_n_hydrogens": 0,
        "reset_product_aromatic_query_flags": 0,
    }
    for index in reversed(range(editable.GetNumAtoms())):
        atom = editable.GetAtomWithIdx(index)
        if atom.GetAtomicNum() != 0:
            continue
        if atom.HasProp("react_atom_idx") or atom.HasProp("react_idx") or atom.GetDegree() != 1:
            raise ReverseHintError("unsafe_dummy")
        bond = atom.GetBonds()[0]
        neighbor = bond.GetOtherAtom(atom)
        if (
            neighbor.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            or bond.GetStereo() != Chem.BondStereo.STEREONONE
            or any(
                adjacent.GetStereo() != Chem.BondStereo.STEREONONE
                or adjacent.GetBondDir() != Chem.BondDir.NONE
                or index in adjacent.GetStereoAtoms()
                for adjacent in neighbor.GetBonds()
            )
        ):
            raise ReverseHintError("unsafe_dummy_stereo")
        editable.RemoveAtom(index)
        edits["removed_originless_terminal_dummy"] += 1
    clean = editable.GetMol()
    edits["reset_template_n_hydrogens"] = _sanitize_reverse_hint(clean, parent, rule, aromatic_changes=edits)
    try:
        smiles = canonical(absolute_smiles(clean, isomericSmiles=True))
    except (ValueError, RuntimeError) as exc:
        raise ReverseHintError("hint_serialization_failed") from exc
    if capture_clean is not None:
        try:
            capture_clean(clean, smiles)
        except Exception:
            # Diagnostics cannot change which retrieval hint is accepted.
            pass
    return smiles, edits


def tagged_reverse_reaction(rule):
    """Tag only a local reverse reaction; the allowlisted rule stays untouched."""
    lhs, rhs = rule.smarts.split(">>")
    reverse = rdChemReactions.ReactionFromSmarts(rhs + ">>" + lhs)
    next_map = max(
        (
            atom.GetAtomMapNum()
            for template in itertools.chain(rule.reaction.GetReactants(), rule.reaction.GetProducts())
            for atom in template.GetAtoms()
        ),
        default=0,
    )
    tags = {}
    for slot, template in enumerate(reverse.GetProducts()):
        for atom in template.GetAtoms():
            if atom.GetAtomMapNum() == 0:
                next_map += 1
                atom.SetAtomMapNum(next_map)
                tags[(slot, atom.GetIdx())] = next_map
    with rdBase.BlockLogs():
        reverse.Initialize()
    return reverse, tags


def reverse_role_issue(clean, role, tags, slot):
    """Check the original role at every locally tagged or mapped atom position."""
    positions_by_map = {}
    for atom in clean.GetAtoms():
        if atom.HasProp("old_mapno"):
            positions_by_map.setdefault(atom.GetIntProp("old_mapno"), []).append(atom.GetIdx())
    expected = {}
    for atom in role.GetAtoms():
        number = atom.GetAtomMapNum() or tags.get((slot, atom.GetIdx()))
        if number is None or len(positions_by_map.get(number, ())) != 1:
            return "missing_or_ambiguous_role_origin"
        expected[atom.GetIdx()] = positions_by_map[number][0]
    matches = clean.GetSubstructMatches(
        role, useChirality=True, uniquify=False, maxMatches=MAX_REVERSE_ROLE_MATCHES
    )
    if len(matches) >= MAX_REVERSE_ROLE_MATCHES:
        return "role_match_cap"
    if not any(all(match[index] == position for index, position in expected.items()) for match in matches):
        return "role_mismatch_at_mapped_positions"
    return None


def typed_reverse_hint(
    raw,
    parent,
    rule,
    tags,
    slot,
    representative=None,
    count_match=None,
    capture_clean=None,
    aromatic_changes=None,
):
    """One role-conforming retrieval query from a raw copied precursor."""
    editable = Chem.RWMol(Chem.Mol(raw))
    dummies = [atom for atom in editable.GetAtoms() if atom.GetAtomicNum() == 0]
    if representative is None:
        if dummies:
            raise ReverseHintError("unsupported_placeholders")
    else:
        if (
            len(dummies) != 1
            or dummies[0].GetDegree() != 1
            or dummies[0].HasProp("react_atom_idx")
            or dummies[0].HasProp("react_idx")
            or not dummies[0].HasProp("old_mapno")
        ):
            raise ReverseHintError("unsupported_placeholders")
        atom = dummies[0]
        bond = atom.GetBonds()[0]
        neighbor = bond.GetOtherAtom(atom)
        if (
            neighbor.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            or bond.GetStereo() != Chem.BondStereo.STEREONONE
            or any(
                adjacent.GetStereo() != Chem.BondStereo.STEREONONE
                or adjacent.GetBondDir() != Chem.BondDir.NONE
                or atom.GetIdx() in adjacent.GetStereoAtoms()
                for adjacent in neighbor.GetBonds()
            )
        ):
            raise ReverseHintError("unsupported_placeholder_stereo")
        atomic_number, charge, hydrogens, _label = representative
        atom.SetAtomicNum(atomic_number)
        atom.SetFormalCharge(charge)
        atom.SetNumExplicitHs(hydrogens)
        atom.SetNoImplicit(True)
    clean = editable.GetMol()
    resets = _sanitize_reverse_hint(clean, parent, rule, aromatic_changes=aromatic_changes)
    for atom in raw.GetAtoms():
        if atom.HasProp("react_atom_idx"):
            same = clean.GetAtomWithIdx(atom.GetIdx())
            if (atom.GetAtomicNum(), atom.GetIsotope(), atom.GetFormalCharge(), atom.GetChiralTag()) != (
                same.GetAtomicNum(),
                same.GetIsotope(),
                same.GetFormalCharge(),
                same.GetChiralTag(),
            ):
                raise ReverseHintError("target_atom_changed")
    for bond in raw.GetBonds():
        same = clean.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if (
            same is None
            or bond.GetBondType() != same.GetBondType()
            or bond.GetStereo() != same.GetStereo()
            or bond.GetBondDir() != same.GetBondDir()
            or tuple(bond.GetStereoAtoms()) != tuple(same.GetStereoAtoms())
        ):
            raise ReverseHintError("bond_or_stereo_changed")
    if count_match is not None:
        count_match()
    issue = reverse_role_issue(clean, rule.reaction.GetReactants()[slot], tags, slot)
    if issue is not None:
        raise ReverseHintError(issue)
    try:
        smiles = canonical(absolute_smiles(clean, isomericSmiles=True))
    except (ValueError, RuntimeError) as exc:
        raise ReverseHintError("hint_serialization_failed") from exc
    if capture_clean is not None:
        try:
            capture_clean(clean, smiles)
        except Exception:
            pass
    return smiles, resets


class HitStereoAnchors:
    """Bounded local hit anchors; compile only representable, unique source mappings."""

    def __init__(self, hit):
        self.source = molecule(hit)
        if self.source.GetNumAtoms() > 128:
            raise ValueError("Hit stereo anchor molecule exceeds bounded comparison size")
        tetra = [
            (atom.GetIdx(),)
            for atom in self.source.GetAtoms()
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        ]
        bonds = [
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            for bond in self.source.GetBonds()
            if bond.GetStereo() != Chem.BondStereo.STEREONONE
        ]
        if len(tetra) + len(bonds) > 16:
            raise ValueError("Too many specified hit stereo anchors")
        if any(
            self.source.GetBondBetweenAtoms(*pair).GetStereo()
            not in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ)
            for pair in bonds
        ):
            raise ValueError("Unsupported specified hit bond stereo")
        self.queries = [self._compile(roots) for roots in tetra + bonds]

    def _compile(self, roots):
        included, frontier = set(roots), set(roots)
        rings = self.source.GetRingInfo().AtomRings()
        for radius in range(1, 9):
            previous_included = set(included)
            frontier = {
                neighbor.GetIdx()
                for index in frontier
                for neighbor in self.source.GetAtomWithIdx(index).GetNeighbors()
            } - included
            included.update(frontier)
            for _ in range(len(rings) + 1):
                previous = len(included)
                for ring in rings:
                    if included.intersection(ring):
                        included.update(ring)
                if len(included) == previous:
                    break
            frontier = included - previous_included
            if len(included) > 64:
                break
            marked = Chem.Mol(self.source)
            for position, index in enumerate(roots):
                marked.GetAtomWithIdx(index).SetAtomMapNum(1000000 + position)
            fragment = Chem.MolFragmentToSmiles(
                marked, atomsToUse=sorted(included), canonical=True, isomericSmiles=True
            )
            query = Chem.MolFromSmiles(fragment)
            if query is None:
                continue
            positions = [
                next(
                    (atom.GetIdx() for atom in query.GetAtoms() if atom.GetAtomMapNum() == 1000000 + i), None
                )
                for i in range(len(roots))
            ]
            if None in positions:
                continue
            if len(roots) == 1:
                defined = query.GetAtomWithIdx(positions[0]).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            else:
                bond = query.GetBondBetweenAtoms(*positions)
                defined = bond is not None and bond.GetStereo() in (
                    Chem.BondStereo.STEREOE,
                    Chem.BondStereo.STEREOZ,
                )
            if not defined:
                continue
            for position in positions:
                query.GetAtomWithIdx(position).SetAtomMapNum(0)
            matches = self.source.GetSubstructMatches(
                query, useChirality=False, uniquify=False, maxMatches=65
            )
            locations = {tuple(match[position] for position in positions) for match in matches}
            chiral_matches = self.source.GetSubstructMatches(
                query, useChirality=True, uniquify=False, maxMatches=65
            )
            chiral_locations = {tuple(match[position] for position in positions) for match in chiral_matches}
            if (
                len(matches) < 65
                and len(chiral_matches) < 65
                and locations == {tuple(roots)}
                and chiral_locations == {tuple(roots)}
            ):
                return query, tuple(positions), radius
            if not frontier:
                break
        raise ValueError("Hit stereo anchor is not uniquely representable within bounded local context")

    def check(self, product):
        candidate = molecule(product)
        if candidate.GetNumAtoms() > 128:
            raise ValueError("Hit stereo anchor molecule exceeds bounded comparison size")
        used_tetra, used_bonds = set(), set()
        for query, positions, _radius in self.queries:
            achiral_matches = candidate.GetSubstructMatches(
                query, useChirality=False, uniquify=False, maxMatches=65
            )
            if len(achiral_matches) >= 65:
                raise ValueError("Hit stereo anchor match bound reached")
            locations = {tuple(match[position] for position in positions) for match in achiral_matches}
            if len(locations) != 1:
                raise ValueError("Hit stereo anchor is missing or has ambiguous correspondence")
            chiral_matches = candidate.GetSubstructMatches(
                query, useChirality=True, uniquify=False, maxMatches=65
            )
            if len(chiral_matches) >= 65:
                raise ValueError("Hit stereo anchor match bound reached")
            if locations != {tuple(match[position] for position in positions) for match in chiral_matches}:
                raise ValueError("Hit stereo anchor is inverted")
            mapped = next(iter(locations))
            used = used_tetra if len(positions) == 1 else used_bonds
            mapped = mapped if len(positions) == 1 else tuple(sorted(mapped))
            if mapped in used:
                raise ValueError("Two hit stereo anchors map to one product anchor")
            used.add(mapped)


def check_hit_stereo_anchors(hit, product):
    """Conservative standalone control; workflows compile anchors only once."""
    HitStereoAnchors(hit).check(product)


@dataclass(frozen=True)
class RouteFirstLimits:
    version: str = VERSION
    max_rule_inspections: int = 512
    max_rule_scans_per_target: int = 128
    max_rules: int = 16
    max_reverse_outcomes: int = 4
    max_outer_guides: int = 8
    max_catalog_index_rows: int = 250000
    max_catalog_rows_per_query: int = 8192
    max_catalog_row_visits: int = 262144
    max_fingerprint_calls: int = 65536
    building_block_top_k: int = 2
    max_combinations: int = 512
    max_products: int = 1024
    max_routes: int = 64
    beam_width: int = 12
    max_steps: int = 2
    max_inner_rules: int = 4
    fallback_rules: int = 2
    reaction_product_limit: int = 16

    def __post_init__(self):
        if self.version != VERSION:
            raise ValueError(f"route_first.version must be {VERSION}")
        for key, value in asdict(self).items():
            if key != "version":
                positive(value, f"route_first.{key}")
        if self.max_steps not in (1, 2):
            raise ValueError("route_first.max_steps must be 1 or 2")
        if self.max_rules > self.max_rule_inspections or self.fallback_rules > self.max_rules:
            raise ValueError("route_first rule budgets are inconsistent")
        ceilings = {
            "max_rule_inspections": 1024,
            "max_reverse_outcomes": 16,
            "max_catalog_index_rows": 500000,
            "max_catalog_rows_per_query": 32768,
            "max_catalog_row_visits": 524288,
            "max_fingerprint_calls": 262144,
            "max_combinations": 4096,
            "max_products": 8192,
            "max_routes": 256,
            "beam_width": 64,
        }
        for key, ceiling in ceilings.items():
            if getattr(self, key) > ceiling:
                raise ValueError(f"route_first.{key} exceeds the implementation ceiling")


class RouteFirstSearch:
    def __init__(
        self,
        hit,
        catalog,
        rules,
        limits,
        deadline_monotonic=None,
        accept=None,
        hit_anchors=None,
        *,
        diagnose_continuation=False,
        prefer_continuation_second=False,
        preserve_typed_prefix=False,
        correct_typed_prefix_size=False,
    ):
        self.hit = canonical(hit)
        self.hit_mol = molecule(self.hit)
        self.hit_fp = FP.GetFingerprint(self.hit_mol)
        self.catalog, self.rules, self.limits = catalog, rules, limits
        self.deadline = deadline_monotonic
        self.accept = accept or (lambda product: None)
        self._admission_fp_reserve = min(2 * limits.max_products, limits.max_fingerprint_calls // 2)
        self._generation_fp_ceiling = limits.max_fingerprint_calls - self._admission_fp_reserve
        self._depth_one_fp_ceiling = min(limits.max_fingerprint_calls // 2, self._generation_fp_ceiling)
        self.report = {
            "version": VERSION,
            "limits": asdict(limits),
            "catalog_sha256": catalog.sha256,
            "rules_sha256": rules.sha256,
            "rule_inspections": 0,
            "rule_template_checks": 0,
            "rule_template_exclusions": [],
            "topology_limits": {
                "match_requests_per_search": MAX_TOPOLOGY_MATCH_REQUESTS,
                "matches_per_request": MAX_TOPOLOGY_MATCHES - 1,
                "spec_and_graph_visits_per_search": MAX_TOPOLOGY_SPEC_GRAPH_VISITS,
            },
            "topology_match_requests": 0,
            "topology_spec_graph_visits": 0,
            "topology_static_checks": 0,
            "topology_static_cache_hits": 0,
            "topology_assessment_cache_hits": 0,
            "topology_pruned_rules": 0,
            "topology_unknown_reasons": {},
            "reverse_nonring_query_compiles": 0,
            "reverse_nonring_compile_failures": 0,
            "reverse_nonring_query_calls": 0,
            "reverse_nonring_empty": 0,
            "rules_selected": [],
            "reverse_attempts": 0,
            "reverse_failures": 0,
            "reverse_partition_checks": 0,
            "reverse_partition_rejections": {},
            "reverse_hint_normalizations": {},
            "reverse_hint_rejections": {},
            "reverse_role_query_limits": {
                "trials_per_raw_outcome": MAX_REVERSE_ROLE_TRIALS_PER_OUTCOME,
                "trial_and_match_work_per_search": MAX_REVERSE_ROLE_WORK_PER_SEARCH,
                "matches_per_query": MAX_REVERSE_ROLE_MATCHES - 1,
                "representative_order": [row[3] for row in REVERSE_ROLE_REPRESENTATIVES],
            },
            "reverse_role_work": 0,
            "reverse_role_trials": 0,
            "reverse_role_match_checks": 0,
            "reverse_role_rejections": {},
            "reverse_hint_quality_counts": {
                "typed": 0,
                "approximate": 0,
                "explicit_fallback": 0,
                "implicit_fallback": 0,
            },
            "reverse_local_template_compiles": 0,
            "outer_center_metadata_check_limit": MAX_OUTER_CENTER_METADATA_CHECKS,
            "outer_center_metadata_checks": 0,
            "outer_center_metadata_incomplete": 0,
            "representative_origin_limits": {
                "structural_units": MAX_REPRESENTATIVE_ORIGIN_VISITS,
                "guide_entries": MAX_REPRESENTATIVE_ORIGIN_GUIDES,
                "atoms_per_part": MAX_REPRESENTATIVE_ORIGIN_ATOMS,
                "bonds_per_part": MAX_REPRESENTATIVE_ORIGIN_BONDS,
            },
            "representative_origin_structural_units": 0,
            "representative_origin_checks": 0,
            "representative_origin_guide_entries": 0,
            "representative_origin_known_guides": 0,
            "representative_origin_legacy_reasons": {},
            "representative_origin_score_corrections": 0,
            "outer_plan_candidates": 0,
            "outer_distinct_cut_centers": 0,
            "phase2_plan_trace": [],
            "phase2_plan_trace_truncated": 0,
            "reverse_outcomes": 0,
            "reverse_fallbacks": 0,
            "retrieval_guides_normalized": 0,
            "retrieval_guides_unusable": 0,
            "catalog_index_rows_scanned": 0,
            "catalog_index_hashes": 0,
            "anchor_index_queries": min(2, len(hit_anchors.queries)) if hit_anchors else 0,
            "anchor_index_queries_unindexed": max(0, len(hit_anchors.queries) - 2) if hit_anchors else 0,
            "anchor_index_match_checks": 0,
            "anchor_index_matches": [],
            "anchor_index_postings_retained": [],
            "anchor_index_postings_truncated": [],
            "anchor_index_incomplete": False,
            "guide_anchor_checks": 0,
            "anchor_posting_metadata_visits": 0,
            "anchor_role_checks": 0,
            "anchor_role_compatible": 0,
            "catalog_row_visits": 0,
            "catalog_stereo_limits": {
                "perceptions": MAX_CATALOG_STEREO_PERCEPTIONS,
                "atom_bond_charges": MAX_CATALOG_STEREO_ATOM_BOND_CHARGES,
                "cache_entries": MAX_CATALOG_STEREO_CACHE_ENTRIES,
                "exclusion_samples": MAX_CATALOG_STEREO_EXCLUSION_SAMPLES,
            },
            "catalog_stereo_perceptions": 0,
            "catalog_stereo_atom_bond_charges": 0,
            "catalog_stereo_cache_hits": 0,
            "catalog_stereo_exclusions": 0,
            "catalog_stereo_unknowns": {},
            "catalog_stereo_exclusion_samples": [],
            "catalog_stereo_exclusion_samples_truncated": 0,
            "fallback_hit_second_limit": MAX_FALLBACK_HIT_SECOND_TANIMOTOS,
            "fallback_hit_second_sample_limit": MAX_FALLBACK_HIT_SECOND_SAMPLES,
            "fallback_hit_second_tanimotos": 0,
            "fallback_hit_second_queries": 0,
            "fallback_hit_second_replacements": 0,
            "fallback_hit_second_legacy_reasons": {},
            "fallback_hit_second_samples": [],
            "fallback_hit_second_samples_omitted": 0,
            "catalog_bucket_rankings": 0,
            "reaction_pool_limits": {
                "query_role_pairs": MAX_POOL_QUERY_ROLE_PAIRS,
                "window_match_requests": MAX_POOL_WINDOW_MATCH_REQUESTS,
                "site_match_requests": MAX_POOL_SITE_MATCH_REQUESTS,
                "matches_per_request": MAX_POOL_MATCHES - 1,
                "query_atoms": MAX_POOL_QUERY_ATOMS,
                "query_bonds": MAX_POOL_QUERY_BONDS,
                "combined_graph_visits": MAX_POOL_GRAPH_VISITS,
                "extra_catalog_role_match_checks": MAX_ROLE_INDEX_MATCH_CHECKS,
                "total_pool_rows": self.limits.max_catalog_rows_per_query,
            },
            "reaction_pool_selected_roles": [],
            "reaction_pool_role_index_checks": 0,
            "reaction_pool_role_index_complete": False,
            "reaction_pool_role_index_failure": None,
            "reaction_pool_role_index_sizes": [],
            "reaction_pool_window_match_requests": 0,
            "reaction_pool_site_match_requests": 0,
            "reaction_pool_graph_visits": 0,
            "reaction_pool_query_pairs": 0,
            "reaction_pool_metadata_cache_hits": 0,
            "reaction_pool_metadata_unknowns": {},
            "reaction_pool_index_choices": {"exact": 0, "same_site": 0, "legacy": 0},
            "reaction_pool_index_applied": 0,
            "reaction_pool_index_disabled_incomplete": 0,
            "reaction_pool_caller_only_sites": 0,
            "reaction_pool_index_only_sites": 0,
            "reaction_pool_indexed_bucket_rankings": 0,
            "reaction_pool_ordinary_bucket_rankings": 0,
            "reaction_pool_indexed_rows_selected": 0,
            "reaction_pool_ordinary_rows_selected": 0,
            "fingerprint_calls": 1,
            "hit_fingerprint_calls": 1,
            "generation_fingerprint_calls": 0,
            "admission_rank_fingerprint_calls": 0,
            "admission_fingerprint_calls": 0,
            "admission_fingerprint_reserve": self._admission_fp_reserve,
            "generation_fingerprint_ceiling": self._generation_fp_ceiling,
            "depth_one_fingerprint_ceiling": self._depth_one_fp_ceiling,
            "admission_pending_unique_products": 0,
            "admission_reserved_cost_units": 0,
            "admission_ranked_products": 0,
            "admission_attempted_products": 0,
            "admission_fp_unassessed_products": 0,
            "admission_callback_unattempted_products": 0,
            "admission_rank_selection_policy": "cached fingerprint first, then canonical SMILES; reserve one possible callback fingerprint plus zero or one new rank fingerprint per selected product",
            "role_queries": 0,
            "combinations": 0,
            "products_generated": 0,
            "inner_role_checks": 0,
            "inner_role_compatible": 0,
            "inner_role_rejected": 0,
            "route_audit_attempts": 0,
            "routes_verified": 0,
            "product_rejections": {},
            "limits_reached": [],
            "selection_policy": "Product-template hit match first; bounded reverse precursor guides; explicit hit-role fallback",
            "retrieval_scope": "One bounded catalog feature index; each guide ranks charged bucket metadata and scans at most the nearest 8192 catalog rows; exact target ID is additionally considered",
            "evidence_scope": "Reverse results only guide retrieval. Forward rule replay and later route audit/evidence are authoritative.",
        }
        self._nearest_cache = {}
        self._rule_selection_cache = {}
        self._topology_rule_cache = {}
        self._topology_assessment_cache = {}
        self._eligible_rules = None
        self._reverse_cache = {}
        self._diagnose_continuation = diagnose_continuation or prefer_continuation_second
        self._continuation_capture_all_phases = diagnose_continuation
        self._prefer_continuation_second = prefer_continuation_second
        self._preserve_typed_prefix = preserve_typed_prefix
        self._correct_typed_prefix_size = correct_typed_prefix_size
        self._guide_representative_heavy = {} if correct_typed_prefix_size else None
        if preserve_typed_prefix:
            self.report["mixed_typed_prefix_hints"] = 0
            self.report["mixed_typed_prefix_parts"] = 0
        self._continuation_raw = {} if self._diagnose_continuation else None
        self._continuation_context_cache = {} if self._diagnose_continuation else None
        self._continuation_snapshots = 0
        self._continuation_snapshot_units = 0
        self._continuation_role_contexts = 0
        self._continuation_ranked = {} if prefer_continuation_second else None
        self._continuation_selected = {} if prefer_continuation_second else None
        self._continuation_outer_roles = {} if prefer_continuation_second else None
        if prefer_continuation_second:
            self.report["continuation_preference"] = {
                "ranked_rows_limit": MAX_CONTINUATION_RANKED_ROWS,
                "role_checks_per_query_limit": MAX_CONTINUATION_ROLE_CHECKS_PER_QUERY,
                "role_checks_per_search_limit": MAX_CONTINUATION_ROLE_CHECKS_PER_SEARCH,
                "outer_role_queries_limit": MAX_CONTINUATION_OUTER_ROLE_QUERIES,
                "selected_contexts_limit": MAX_CONTINUATION_SELECTED_CONTEXTS,
                "ranked_rows_retained": 0,
                "role_checks": 0,
                "role_matches": 0,
                "stereo_rejections": 0,
                "stereo_unknowns": 0,
                "outer_role_queries": 0,
                "replacements": 0,
                "legacy_reasons": {},
                "samples": [],
                "samples_omitted": 0,
                "bridge_positive_work_items": 0,
                "bridge_unknown_work_items": 0,
                "bridge_unknown_reasons": {},
                "bridge_snapshot_attempts": 0,
                "bridge_structural_units": 0,
                "bridge_role_contexts": 0,
            }
        self._tagged_reverse_cache = {}
        self._nonring_unavailable_rules = set()
        self._guide_quality = {}
        self._guide_missing_role_slots = {}
        self._guide_centers = {}
        self._active_phase2_slot = None
        self._retrieval_target_cache = {}
        self._fp_cache = {self.hit: self.hit_fp}
        self._catalog_index = None
        self._catalog_pools = {}
        self._pool_last_complete = True
        self._pool_metadata_cache = {}
        self._selected_role_indexes = ()
        self._role_index_queries = {}
        self._role_index_buckets = {}
        self._role_ordinary_buckets = {}
        self._catalog_stereo_cache = {}
        self._anchor_queries = tuple(item[0] for item in hit_anchors.queries[:2]) if hit_anchors else ()
        self._anchor_postings = ()
        self._anchor_scores = {}
        self._guide_anchor_cache = {}
        self._phase = 1
        self.routes = []
        self._route_keys = set()
        self._pending = []

    def _stop(self, name=None):
        if self.deadline is not None and time.monotonic() >= self.deadline:
            name = "deadline_monotonic"
        if name is not None and name not in self.report["limits_reached"]:
            self.report["limits_reached"].append(name)
        return name is not None

    def _rank(self, smiles, target=None):
        if self._stop():
            return 0
        fp = self._fingerprint(smiles)
        reference = self.hit_fp if target is None else self._fingerprint(target)
        return (
            -DataStructs.TanimotoSimilarity(reference, fp) if fp is not None and reference is not None else 0
        )

    def _fingerprint(self, smiles, *, purpose="generation"):
        if smiles in self._fp_cache:
            return self._fp_cache[smiles]
        if purpose == "admission":
            phase_limit = self.limits.max_fingerprint_calls
            limit_reason = "admission_fingerprint_budget"
        elif purpose == "generation":
            depth_one = self._phase == 1 and self.limits.max_steps == 2
            phase_limit = self._depth_one_fp_ceiling if depth_one else self._generation_fp_ceiling
            limit_reason = "depth_one_fingerprint_budget" if depth_one else "generation_fingerprint_budget"
        else:
            raise ValueError("Unknown route-first fingerprint purpose")
        if self.report["fingerprint_calls"] + getattr(self.accept, "fingerprint_calls", 0) >= phase_limit:
            self._stop(limit_reason)
            return None
        self.report["fingerprint_calls"] += 1
        self.report[
            "admission_rank_fingerprint_calls" if purpose == "admission" else "generation_fingerprint_calls"
        ] += 1
        fp = FP.GetFingerprint(molecule(smiles))
        self._fp_cache[smiles] = fp
        return fp

    @staticmethod
    def _features(mol):
        heavy = mol.GetNumHeavyAtoms()
        hetero = sum(atom.GetAtomicNum() not in (1, 6) for atom in mol.GetAtoms())
        aromatic = sum(atom.GetIsAromatic() for atom in mol.GetAtoms())
        rings = mol.GetRingInfo().NumRings()
        return (heavy // 4, hetero // 2, aromatic // 4, min(rings, 8))

    def _pool_unknown(self, reason):
        counts = self.report["reaction_pool_metadata_unknowns"]
        counts[reason] = counts.get(reason, 0) + 1

    def _pool_charge_graph(self, visits):
        if self._stop() or self.report["reaction_pool_graph_visits"] + visits > MAX_POOL_GRAPH_VISITS:
            return False
        self.report["reaction_pool_graph_visits"] += visits
        return True

    def _pool_matches(self, mol, query, kind):
        field = f"reaction_pool_{kind}_match_requests"
        ceiling = MAX_POOL_WINDOW_MATCH_REQUESTS if kind == "window" else MAX_POOL_SITE_MATCH_REQUESTS
        if self.report[field] >= ceiling:
            return None, f"{kind}_match_request_cap"
        if not self._pool_charge_graph(mol.GetNumAtoms() + mol.GetNumBonds()):
            return None, "graph_visit_cap_or_deadline"
        self.report[field] += 1
        try:
            matches = mol.GetSubstructMatches(
                query, useChirality=True, uniquify=True, maxMatches=MAX_POOL_MATCHES
            )
        except (ValueError, RuntimeError, TypeError):
            return None, f"{kind}_match_error"
        if len(matches) >= MAX_POOL_MATCHES:
            return None, f"{kind}_match_truncated"
        if not self._pool_charge_graph(sum(map(len, matches))):
            return None, "graph_visit_cap_or_deadline"
        return matches, None

    def _pool_window_counts(self, mol, role_query, full_counts):
        matches, reason = self._pool_matches(mol, role_query, "window")
        if reason or not matches:
            return full_counts, reason or "no_role_match", reason is None
        matched = {index for match in matches for index in match}
        if not self._pool_charge_graph(mol.GetNumAtoms() + 3 * mol.GetNumBonds()):
            return full_counts, "graph_visit_cap_or_deadline", False
        stereo = {
            atom.GetIdx() for atom in mol.GetAtoms() if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        }
        for bond in mol.GetBonds():
            if bond.GetStereo() != Chem.BondStereo.STEREONONE or bond.GetBondDir() != Chem.BondDir.NONE:
                stereo.update((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
                stereo.update(index for index in bond.GetStereoAtoms() if index >= 0)
        stereo.update(
            neighbor.GetIdx()
            for index in tuple(stereo)
            for neighbor in mol.GetAtomWithIdx(index).GetNeighbors()
        )
        core = set(matched)
        for index in matched:
            core.update(atom.GetIdx() for atom in mol.GetAtomWithIdx(index).GetNeighbors())
        core.update(stereo)
        rings = mol.GetRingInfo().AtomRings()
        if not self._pool_charge_graph(sum(map(len, rings)) + len(core) + len(rings) ** 2):
            return full_counts, "graph_visit_cap_or_deadline", False
        systems = _pool_ring_systems(mol)
        core.update(index for system in systems if system & core for index in system)
        retained = set(core)
        for index in core:
            retained.update(atom.GetIdx() for atom in mol.GetAtomWithIdx(index).GetNeighbors())
        if not self._pool_charge_graph(len(retained) + sum(map(len, rings))):
            return full_counts, "graph_visit_cap_or_deadline", False
        return _pool_raw_counts(mol, retained), None, True

    def _pool_site_choice(self, mol, caller, normalized):
        if normalized in self._selected_role_indexes:
            return normalized, "exact", True
        caller_matches, reason = self._pool_matches(mol, caller, "site")
        if reason:
            return None, reason, False
        caller_sites = {frozenset(match) for match in caller_matches}
        selected = None
        index_sites = set()
        for indexed_role in self._selected_role_indexes:
            matches, reason = self._pool_matches(mol, self._role_index_queries[indexed_role], "site")
            if reason:
                return None, reason, False
            sites = {frozenset(match) for match in matches}
            index_sites.update(sites)
            if selected is None and caller_sites & sites:
                selected = indexed_role
        self.report["reaction_pool_caller_only_sites"] += len(caller_sites - index_sites)
        self.report["reaction_pool_index_only_sites"] += len(index_sites - caller_sites)
        return selected, "same_site" if selected is not None else "legacy", True

    def _pool_metadata(self, target, guide, role):
        key = (target, guide, role)
        if key in self._pool_metadata_cache:
            self.report["reaction_pool_metadata_cache_hits"] += 1
            return self._pool_metadata_cache[key]
        full_mol = molecule(target)
        full_features = self._features(full_mol)
        if full_mol.GetNumAtoms() > MAX_POOL_QUERY_ATOMS or full_mol.GetNumBonds() > MAX_POOL_QUERY_BONDS:
            self._pool_unknown("query_size_cap")
            return full_features, None, False
        if self.report["reaction_pool_query_pairs"] >= MAX_POOL_QUERY_ROLE_PAIRS:
            self._pool_unknown("query_pair_cache_cap")
            return full_features, None, False
        self.report["reaction_pool_query_pairs"] += 1
        try:
            guide_mol = molecule(guide)
            caller = Chem.MolFromSmarts(role)
            normalized = _normalized_role_smarts(role)
        except (ValueError, RuntimeError, TypeError):
            self._pool_unknown("query_or_role_parse")
            return full_features, None, False
        if (
            caller is None
            or guide_mol.GetNumAtoms() > MAX_POOL_QUERY_ATOMS
            or guide_mol.GetNumBonds() > MAX_POOL_QUERY_BONDS
        ):
            self._pool_unknown("query_size_or_role_parse")
            return full_features, None, False
        if not self._pool_charge_graph(
            full_mol.GetNumAtoms()
            + full_mol.GetNumBonds()
            + sum(map(len, full_mol.GetRingInfo().AtomRings()))
        ):
            self._pool_unknown("graph_visit_cap_or_deadline")
            return full_features, None, False
        full_counts = _pool_raw_counts(full_mol)
        local_counts, reason, window_complete = self._pool_window_counts(guide_mol, caller, full_counts)
        if reason is not None:
            self._pool_unknown(reason)
        if not window_complete:
            return full_features, None, False
        midpoint = tuple((full + local + 1) // 2 for full, local in zip(full_counts, local_counts))
        features = _pool_bins(midpoint)
        selected, choice, site_complete = self._pool_site_choice(guide_mol, caller, normalized)
        if not site_complete:
            self._pool_unknown(choice)
            return full_features, None, False
        self.report["reaction_pool_index_choices"][choice] += 1
        complete = window_complete and site_complete
        outcome = features, selected if site_complete else None, complete
        if complete:
            self._pool_metadata_cache[key] = outcome
        return outcome

    def _select_role_indexes(self, plans):
        counts = {}
        for rule, _guide, _fallback, _keys, _kind, _quality in plans:
            for role in rule.roles:
                normalized = _normalized_role_smarts(role)
                counts[normalized] = counts.get(normalized, 0) + 1
        selected = tuple(
            role for role, _count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:2]
        )
        self._selected_role_indexes = selected
        self._role_index_queries = {role: Chem.MolFromSmarts(role) for role in selected}
        self.report["reaction_pool_selected_roles"] = [
            {"normalized_smarts": role, "occurrences": counts[role]} for role in selected
        ]

    def _index(self):
        if self._catalog_index is not None:
            return self._catalog_index
        buckets = {}
        indexed = {role: {} for role in self._selected_role_indexes}
        ordinary = {role: {} for role in self._selected_role_indexes}
        role_index_complete = bool(self._selected_role_indexes)
        reserve = min(256, max(1, self.limits.max_catalog_rows_per_query // 4))
        postings = [[] for _ in self._anchor_queries]
        match_totals = [0 for _ in self._anchor_queries]
        for smiles, mol in itertools.islice(self.catalog.mols.items(), self.limits.max_catalog_index_rows):
            if self._stop():
                break
            self.report["catalog_index_rows_scanned"] += 1
            self.report["catalog_index_hashes"] += 1
            key = self._features(mol)
            score = int.from_bytes(hashlib.sha256(smiles.encode()).digest()[:8], "big")
            buckets.setdefault(key, []).append((score, smiles))
            if role_index_complete:
                for role in self._selected_role_indexes:
                    if self._stop():
                        role_index_complete = False
                        self.report["reaction_pool_role_index_failure"] = "deadline"
                        break
                    if self.report["reaction_pool_role_index_checks"] >= MAX_ROLE_INDEX_MATCH_CHECKS:
                        role_index_complete = False
                        self.report["reaction_pool_role_index_failure"] = "match_check_cap"
                        break
                    self.report["reaction_pool_role_index_checks"] += 1
                    try:
                        matches_role = mol.HasSubstructMatch(
                            self._role_index_queries[role], useChirality=True
                        )
                    except (ValueError, RuntimeError, TypeError):
                        role_index_complete = False
                        self.report["reaction_pool_role_index_failure"] = "role_match_error"
                        break
                    group = indexed if matches_role else ordinary
                    group[role].setdefault(key, []).append((score, smiles))
            for anchor_index, anchor in enumerate(self._anchor_queries):
                if self.report["anchor_index_match_checks"] >= 500000:
                    self.report["anchor_index_incomplete"] = True
                    break
                self.report["anchor_index_match_checks"] += 1
                if mol.HasSubstructMatch(anchor, useChirality=True):
                    match_totals[anchor_index] += 1
                    heap = postings[anchor_index]
                    item = (-score, smiles)
                    if len(heap) < reserve:
                        heapq.heappush(heap, item)
                    elif item > heap[0]:
                        heapq.heapreplace(heap, item)
        if len(self.catalog.mols) > self.limits.max_catalog_index_rows:
            self._stop("max_catalog_index_rows")
        for rows in buckets.values():
            rows.sort()
        if role_index_complete and self.report["catalog_index_rows_scanned"] == len(self.catalog.mols):
            for group in (indexed, ordinary):
                for role_buckets in group.values():
                    for rows in role_buckets.values():
                        rows.sort()
            self._role_index_buckets = indexed
            self._role_ordinary_buckets = ordinary
            self.report["reaction_pool_role_index_complete"] = True
            self.report["reaction_pool_role_index_sizes"] = [
                sum(map(len, indexed[role].values())) for role in self._selected_role_indexes
            ]
        else:
            self._role_index_buckets = {}
            self._role_ordinary_buckets = {}
            if self._selected_role_indexes and self.report["reaction_pool_role_index_failure"] is None:
                self.report["reaction_pool_role_index_failure"] = "incomplete_catalog_index"
        self._anchor_postings = tuple(
            tuple(smiles for _score, smiles in sorted(heap, reverse=True)) for heap in postings
        )
        self._anchor_scores = {
            smiles: -negative_score for heap in postings for negative_score, smiles in heap
        }
        self.report["anchor_index_matches"] = match_totals
        self.report["anchor_index_postings_retained"] = [len(posting) for posting in self._anchor_postings]
        self.report["anchor_index_postings_truncated"] = [
            total > len(posting) for total, posting in zip(match_totals, self._anchor_postings)
        ]
        if self.report["catalog_index_rows_scanned"] < min(
            len(self.catalog.mols), self.limits.max_catalog_index_rows
        ):
            self.report["anchor_index_incomplete"] = True
        self._catalog_index = buckets
        self.report["catalog_index_size"] = sum(map(len, buckets.values()))
        self.report["catalog_index_buckets"] = len(buckets)
        return buckets

    def _pool(self, target, *, guide=None, role=None):
        buckets = self._index()
        if self._stop():
            self._pool_last_complete = False
            return ()
        if guide is not None and role is not None:
            features, indexed_role, metadata_complete = self._pool_metadata(target, guide, role)
        else:
            features, indexed_role, metadata_complete = self._features(molecule(target)), None, True
        if not self.report["reaction_pool_role_index_complete"]:
            if indexed_role is not None:
                self.report["reaction_pool_index_disabled_incomplete"] += 1
            indexed_role = None
        elif indexed_role is not None:
            self.report["reaction_pool_index_applied"] += 1
        cache_key = (target, features, indexed_role)
        self._pool_last_complete = metadata_complete
        if cache_key in self._catalog_pools:
            return self._catalog_pools[cache_key]
        phase_limit = (
            self.limits.max_catalog_row_visits // 2
            if self._phase == 1 and self.limits.max_steps == 2
            else self.limits.max_catalog_row_visits
        )

        def distance(key):
            return (
                4 * abs(key[0] - features[0])
                + 3 * abs(key[1] - features[1])
                + 2 * abs(key[2] - features[2])
                + 2 * abs(key[3] - features[3]),
                key,
            )

        selected = []

        def take(group, label):
            if self._stop() or self.report["catalog_row_visits"] + len(group) > phase_limit:
                self._stop("max_catalog_row_visits")
                return False
            self.report["catalog_row_visits"] += len(group)
            self.report["catalog_bucket_rankings"] += len(group)
            self.report[f"reaction_pool_{label}_bucket_rankings"] += len(group)
            before = len(selected)
            for key in sorted(group, key=distance):
                remaining = self.limits.max_catalog_rows_per_query - len(selected)
                if remaining <= 0:
                    break
                selected.extend(smiles for _, smiles in group[key][:remaining])
            self.report[f"reaction_pool_{label}_rows_selected"] += len(selected) - before
            return True

        if indexed_role is not None:
            if not take(self._role_index_buckets[indexed_role], "indexed"):
                self._pool_last_complete = False
                return ()
            if len(selected) < self.limits.max_catalog_rows_per_query and not take(
                self._role_ordinary_buckets[indexed_role], "ordinary"
            ):
                self._pool_last_complete = False
                return ()
        elif not take(buckets, "ordinary"):
            self._pool_last_complete = False
            return ()
        pool = tuple(selected)
        if metadata_complete:
            self._catalog_pools[cache_key] = pool
        return pool

    def _eligible_rule_ids(self):
        if self._eligible_rules is not None:
            return self._eligible_rules
        eligible = []
        bounded_ids = sorted(self.rules.rules)[: self.limits.max_rule_scans_per_target]
        for ident in bounded_ids:
            if self._stop():
                break
            if self.report["rule_inspections"] >= self.limits.max_rule_inspections:
                self._stop("max_rule_inspections")
                break
            self.report["rule_inspections"] += 1
            self.report["rule_template_checks"] += 1
            rule = self.rules.get(ident)
            _, issues = template_maps(rule.reaction)
            if any(issues.values()):
                self.report["rule_template_exclusions"].append({"rule_id": ident, "template_audit": issues})
            else:
                eligible.append(ident)
        # A partial precheck is usable, but never admits uninspected rules.
        self._eligible_rules = tuple(eligible)
        return self._eligible_rules

    def _topology_unknown(self, reason):
        counts = self.report["topology_unknown_reasons"]
        counts[reason] = counts.get(reason, 0) + 1
        return "unknown"

    def _topology_visit(self, count=1):
        if self._stop():
            return False
        if self.report["topology_spec_graph_visits"] + count > MAX_TOPOLOGY_SPEC_GRAPH_VISITS:
            self._stop("max_topology_spec_graph_visits")
            return False
        self.report["topology_spec_graph_visits"] += count
        return True

    def _topology_spec(self, rule):
        if rule.id in self._topology_rule_cache:
            self.report["topology_static_cache_hits"] += 1
            return self._topology_rule_cache[rule.id]
        cost = sum(
            part.GetNumAtoms() + part.GetNumBonds()
            for part in itertools.chain(rule.reaction.GetReactants(), rule.reaction.GetProducts())
        )
        if not self._topology_visit(cost):
            return None, "spec_budget_or_deadline"
        self.report["topology_static_checks"] += 1
        result = topology_rule_specification(rule.reaction)
        # At most the fixed eligible rule prefix can enter this cache.
        self._topology_rule_cache[rule.id] = result
        return result

    def _topology_assess(self, target, parent, rule):
        key = (target, rule.id)
        if key in self._topology_assessment_cache:
            self.report["topology_assessment_cache_hits"] += 1
            return self._topology_assessment_cache[key]
        specification, reason = self._topology_spec(rule)
        if reason:
            status = self._topology_unknown(reason)
        elif self.report["topology_match_requests"] >= MAX_TOPOLOGY_MATCH_REQUESTS:
            self._stop("max_topology_match_requests")
            status = self._topology_unknown("match_request_budget")
        elif self._stop():
            status = self._topology_unknown("deadline")
        else:
            self.report["topology_match_requests"] += 1
            query, maps, cross = specification
            matches = parent.GetSubstructMatches(
                query, useChirality=True, uniquify=False, maxMatches=MAX_TOPOLOGY_MATCHES
            )
            if len(matches) >= MAX_TOPOLOGY_MATCHES:
                status = self._topology_unknown("match_truncation")
            elif not matches:
                # The existing product match was positive. A disagreement is
                # uncertainty, never evidence that the rule is impossible.
                status = self._topology_unknown("no_enumerated_match")
            else:
                status = "impossible"
                for match in matches:
                    removed = {frozenset((match[a], match[b])) for a, b in cross}
                    if not self._topology_visit(parent.GetNumAtoms()):
                        status = self._topology_unknown("graph_budget_or_deadline")
                        break
                    neighbors = [[] for _ in range(parent.GetNumAtoms())]
                    complete = True
                    for bond in parent.GetBonds():
                        if not self._topology_visit():
                            complete = False
                            break
                        first, second = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                        if frozenset((first, second)) not in removed:
                            neighbors[first].append(second)
                            neighbors[second].append(first)
                    if not complete:
                        status = self._topology_unknown("graph_budget_or_deadline")
                        break
                    components = {}
                    component = 0
                    for origin in range(parent.GetNumAtoms()):
                        if origin in components:
                            continue
                        stack = [origin]
                        components[origin] = component
                        while stack:
                            if not self._topology_visit():
                                complete = False
                                break
                            for neighbor in neighbors[stack.pop()]:
                                if not self._topology_visit():
                                    complete = False
                                    break
                                if neighbor not in components:
                                    components[neighbor] = component
                                    stack.append(neighbor)
                            if not complete:
                                break
                        if not complete:
                            break
                        component += 1
                    if not complete:
                        status = self._topology_unknown("graph_budget_or_deadline")
                        break
                    slots = (set(), set())
                    for atom in query.GetAtoms():
                        if not self._topology_visit():
                            complete = False
                            break
                        number = atom.GetAtomMapNum()
                        if number:
                            slots[maps[number][0]].add(components[match[atom.GetIdx()]])
                    if not complete:
                        status = self._topology_unknown("graph_budget_or_deadline")
                        break
                    if not slots[0].isdisjoint(slots[1]):
                        continue
                    status = "feasible"
                    break
        # Only bounded completed assessments are cached. A budget miss remains
        # a legacy-eligible rule without consuming a result cache slot.
        if (
            self.report["topology_match_requests"] <= MAX_TOPOLOGY_MATCH_REQUESTS
            and len(self._topology_assessment_cache) < MAX_TOPOLOGY_MATCH_REQUESTS
            and status in ("feasible", "impossible")
        ):
            self._topology_assessment_cache[key] = status
        return status

    def _select_rules(self, target, maximum, fallback=0):
        key = (target, maximum, fallback, self._phase == 2)
        if key in self._rule_selection_cache:
            return self._rule_selection_cache[key]
        mol = molecule(target)
        primary, secondary = [], []
        eligible_ids = self._eligible_rule_ids()
        # For inner searches, this is the same ordering used by the original
        # full scan's primary sort, so the first maximum matches are final.
        if not fallback:
            eligible_ids = tuple(
                sorted(
                    eligible_ids,
                    key=lambda ident: (
                        -self.rules.get(ident).reaction.GetProducts()[0].GetNumHeavyAtoms(),
                        ident,
                    ),
                )
            )
        complete = True
        for ident in eligible_ids:
            if self._stop():
                complete = False
                break
            if self.report["rule_inspections"] >= self.limits.max_rule_inspections:
                self._stop("max_rule_inspections")
                complete = False
                break
            self.report["rule_inspections"] += 1
            rule = self.rules.get(ident)
            product_template = rule.reaction.GetProducts()[0]
            if mol.HasSubstructMatch(product_template, useChirality=True):
                if (
                    self._phase == 2
                    and not fallback
                    and self._topology_assess(target, mol, rule) == "impossible"
                ):
                    self.report["topology_pruned_rules"] += 1
                    continue
                primary.append(rule)
                if not fallback and len(primary) >= maximum:
                    break
            elif fallback and any(
                mol.HasSubstructMatch(q, useChirality=True) for q in rule.reaction.GetReactants()
            ):
                secondary.append(rule)
        primary.sort(key=lambda rule: (-rule.reaction.GetProducts()[0].GetNumHeavyAtoms(), rule.id))
        secondary.sort(key=lambda rule: rule.id)
        selected = [(rule, False) for rule in primary[:maximum]] + [
            (rule, True) for rule in secondary[: max(0, min(fallback, maximum - len(primary)))]
        ]
        if complete and not self._stop():
            self._rule_selection_cache[key] = selected
        return selected

    def _role_work(self, *, trial=False):
        if self._stop():
            return False
        if self.report["reverse_role_work"] >= MAX_REVERSE_ROLE_WORK_PER_SEARCH:
            self._stop("reverse_role_work_budget")
            return False
        self.report["reverse_role_work"] += 1
        if trial:
            self.report["reverse_role_trials"] += 1
        return True

    def _role_rejection(self, reason):
        counts = self.report["reverse_role_rejections"]
        counts[reason] = counts.get(reason, 0) + 1

    def _count_role_match(self):
        self.report["reverse_role_match_checks"] += 1

    def _outer_center_work(self):
        if self._stop():
            return False
        if self.report["outer_center_metadata_checks"] >= MAX_OUTER_CENTER_METADATA_CHECKS:
            self._stop("outer_center_metadata_budget")
            return False
        self.report["outer_center_metadata_checks"] += 1
        return True

    def _reverse_cut_bonds(self, outcome, parent):
        """Cross-precursor bonds of original parent atoms, after partition validation."""
        owners = {}
        for slot, precursor in enumerate(outcome):
            for atom in precursor.GetAtoms():
                if not self._outer_center_work():
                    self.report["outer_center_metadata_incomplete"] += 1
                    return None
                if atom.HasProp("react_atom_idx"):
                    owners[atom.GetIntProp("react_atom_idx")] = slot
        cuts = []
        for bond in parent.GetBonds():
            if not self._outer_center_work():
                self.report["outer_center_metadata_incomplete"] += 1
                return None
            first, second = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if owners[first] != owners[second]:
                cuts.append((min(first, second), max(first, second)))
        return tuple(sorted(cuts))

    def _representative_origin_count(self, outcome, clean_parts, typed_prefix_count, parent_atom_count):
        """Count only surviving typed representatives in the emitted guide.

        Parent-derived atoms and bonds must still be identifiable in every
        accepted clean part.  This metadata never changes guide acceptance.
        """
        if self._stop():
            return None, "deadline"
        if (
            clean_parts is None
            or len(clean_parts) != len(outcome)
            or any(part is None for part in clean_parts)
        ):
            return None, "clean_part_unavailable"
        if typed_prefix_count > len(outcome):
            return None, "typed_slot_mismatch"
        parts = tuple(outcome) + tuple(clean_parts)
        if any(
            part.GetNumAtoms() > MAX_REPRESENTATIVE_ORIGIN_ATOMS
            or part.GetNumBonds() > MAX_REPRESENTATIVE_ORIGIN_BONDS
            for part in parts
        ):
            return None, "part_size_limit"
        # Precharge three structural passes, including failed proofs.  These
        # are bounded work units rather than literal Python loop iterations.
        units = 3 * sum(part.GetNumAtoms() + part.GetNumBonds() for part in parts)
        if self.report["representative_origin_structural_units"] + units > MAX_REPRESENTATIVE_ORIGIN_VISITS:
            return None, "structural_unit_limit"
        # Precharge even a failed proof.  This meter is separate from reverse,
        # role, catalog and fingerprint budgets.
        self.report["representative_origin_structural_units"] += units
        self.report["representative_origin_checks"] += 1

        def parent_atoms(part):
            result = {}
            for atom in part.GetAtoms():
                has_index = atom.HasProp("react_atom_idx")
                has_slot = atom.HasProp("react_idx")
                if has_index != has_slot:
                    raise ValueError("partial_parent_origin")
                if not has_index:
                    continue
                if atom.GetIntProp("react_idx") != 0:
                    raise ValueError("wrong_parent_slot")
                index = atom.GetIntProp("react_atom_idx")
                if index in result:
                    raise ValueError("duplicate_parent_origin")
                result[index] = (
                    atom.GetAtomicNum(),
                    atom.GetIsotope(),
                    atom.GetFormalCharge(),
                    int(atom.GetChiralTag()),
                )
            return result

        try:
            count = 0
            covered_parent_atoms = set()
            for slot, (raw, clean) in enumerate(zip(outcome, clean_parts)):
                raw_parent_atoms = parent_atoms(raw)
                if raw_parent_atoms != parent_atoms(clean):
                    return None, "parent_atom_changed"
                if covered_parent_atoms.intersection(raw_parent_atoms):
                    return None, "duplicate_parent_origin"
                covered_parent_atoms.update(raw_parent_atoms)

                def parent_bonds(part):
                    result = []
                    for bond in part.GetBonds():
                        first, second = bond.GetBeginAtom(), bond.GetEndAtom()
                        if not (
                            first.HasProp("react_idx")
                            and first.HasProp("react_atom_idx")
                            and second.HasProp("react_idx")
                            and second.HasProp("react_atom_idx")
                        ):
                            continue
                        endpoints = sorted(
                            (first.GetIntProp("react_atom_idx"), second.GetIntProp("react_atom_idx"))
                        )
                        result.append((tuple(endpoints), str(bond.GetBondType()), int(bond.GetStereo())))
                    return sorted(result)

                if parent_bonds(raw) != parent_bonds(clean):
                    return None, "parent_bond_changed"
                if slot >= typed_prefix_count:
                    continue
                if clean.GetNumAtoms() != raw.GetNumAtoms() or clean.GetNumBonds() != raw.GetNumBonds():
                    return None, "typed_graph_changed"
                for atom in raw.GetAtoms():
                    if atom.GetAtomicNum() != 0:
                        continue
                    if (
                        atom.HasProp("react_atom_idx")
                        or atom.HasProp("react_idx")
                        or not atom.HasProp("old_mapno")
                    ):
                        return None, "dummy_origin_unclear"
                    replacement = clean.GetAtomWithIdx(atom.GetIdx())
                    if replacement.GetAtomicNum() <= 1 or replacement.HasProp("react_atom_idx"):
                        return None, "representative_not_retained"
                    count += 1
                for bond in raw.GetBonds():
                    same = clean.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                    if (
                        same is None
                        or bond.GetBondType() != same.GetBondType()
                        or bond.GetIsAromatic() != same.GetIsAromatic()
                        or bond.GetStereo() != same.GetStereo()
                    ):
                        return None, "typed_bond_changed"
            if covered_parent_atoms != set(range(parent_atom_count)):
                return None, "parent_coverage_incomplete"
            if self._stop():
                return None, "deadline"
            return count, None
        except (ValueError, RuntimeError, KeyError, AttributeError):
            return None, "origin_proof_error"

    @staticmethod
    def _merge_representative_origin(counts, unknown, guide, count, reason):
        if reason is not None:
            unknown[guide] = reason
            counts.pop(guide, None)
        elif guide not in unknown:
            previous = counts.setdefault(guide, count)
            if previous != count:
                counts.pop(guide, None)
                unknown[guide] = "raw_alternatives_disagree"

    def _typed_reverse_outcome(
        self, outcome, parent, rule, tags, *, capture_parts=None, capture_partial=None
    ):
        parts, resets, aromatic_resets, instantiated, trials = [], 0, 0, 0, 0
        for slot, raw in enumerate(outcome):
            if self._stop():
                return None
            dummies = [atom for atom in raw.GetAtoms() if atom.GetAtomicNum() == 0]
            if dummies:
                selected = None
                for representative in REVERSE_ROLE_REPRESENTATIVES:
                    captured = []
                    aromatic_changes = {}
                    if trials >= MAX_REVERSE_ROLE_TRIALS_PER_OUTCOME:
                        self._role_rejection("outcome_trial_budget")
                        break
                    if not self._role_work(trial=True):
                        self._role_rejection("search_trial_and_match_budget")
                        break
                    trials += 1
                    try:
                        selected = typed_reverse_hint(
                            raw,
                            parent,
                            rule,
                            tags,
                            slot,
                            representative,
                            count_match=self._count_role_match,
                            aromatic_changes=aromatic_changes,
                            **(
                                {"capture_clean": lambda mol, smiles: captured.append((mol, smiles))}
                                if capture_parts is not None
                                else {}
                            ),
                        )
                    except ReverseHintError as exc:
                        self._role_rejection(exc.reason)
                        continue
                    except (ValueError, RuntimeError):
                        self._role_rejection("typed_query_failed")
                        continue
                    break
                if selected is None:
                    return None
                part, reset_count = selected
                instantiated += 1
            else:
                captured = []
                aromatic_changes = {}
                if trials >= MAX_REVERSE_ROLE_TRIALS_PER_OUTCOME:
                    self._role_rejection("outcome_trial_budget")
                    return None
                if not self._role_work(trial=True):
                    self._role_rejection("search_trial_and_match_budget")
                    return None
                trials += 1
                try:
                    part, reset_count = typed_reverse_hint(
                        raw,
                        parent,
                        rule,
                        tags,
                        slot,
                        count_match=self._count_role_match,
                        aromatic_changes=aromatic_changes,
                        **(
                            {"capture_clean": lambda mol, smiles: captured.append((mol, smiles))}
                            if capture_parts is not None
                            else {}
                        ),
                    )
                except ReverseHintError as exc:
                    self._role_rejection(exc.reason)
                    return None
                except (ValueError, RuntimeError):
                    self._role_rejection("typed_query_failed")
                    return None
            parts.append(part)
            if capture_parts is not None:
                capture_parts.append(
                    captured[0][0] if len(captured) == 1 and captured[0][1] == part else None
                )
            resets += reset_count
            aromatic_resets += aromatic_changes.get("reset_product_aromatic_query_flags", 0)
            if capture_partial is not None:
                capture_partial.append(
                    (
                        part,
                        {
                            "typed_originless_placeholder_instantiated": int(bool(dummies)),
                            "typed_reset_template_n_hydrogens": reset_count,
                            "typed_reset_product_aromatic_query_flags": aromatic_changes.get(
                                "reset_product_aromatic_query_flags", 0
                            ),
                        },
                    )
                )
        return tuple(parts), {
            "typed_originless_placeholder_instantiated": instantiated,
            "typed_reset_template_n_hydrogens": resets,
            "typed_reset_product_aromatic_query_flags": aromatic_resets,
        }

    def _capture_continuation_raw(self, registry, guide, outcome, clean_parts):
        """Attach one emitted guide to its exact raw alternative; never affect it."""
        if not self._diagnose_continuation:
            return
        if len(clean_parts) != len(outcome) or any(part is None for part in clean_parts):
            registry["incomplete"] = "normalized_molecule_unavailable"
            return
        if len(outcome) > MAX_CONTINUATION_PRECURSORS:
            registry["incomplete"] = "snapshot_precursor_count_limit"
            return
        parts = tuple(outcome) + tuple(clean_parts)
        if any(
            part.GetNumAtoms() > MAX_CONTINUATION_ROLE_ATOMS
            or part.GetNumBonds() > MAX_CONTINUATION_ROLE_BONDS
            for part in parts
        ):
            registry["incomplete"] = "snapshot_part_size_limit"
            return
        units = sum(part.GetNumAtoms() + part.GetNumBonds() for part in parts)
        if self._continuation_snapshots >= MAX_CONTINUATION_RAW_SNAPSHOTS or (
            self._continuation_snapshot_units + units > MAX_CONTINUATION_SNAPSHOT_UNITS
        ):
            registry["incomplete"] = "snapshot_budget"
            return
        # Reserve attempted copy work before any allocation. Failed partial
        # copies still consume the diagnostic budget and cannot be retried free.
        self._continuation_snapshots += 1
        self._continuation_snapshot_units += units
        try:
            captured = (
                tuple(Chem.Mol(part) for part in outcome),
                tuple(Chem.Mol(part) for part in clean_parts),
            )
        except Exception:
            registry["incomplete"] = "snapshot_copy_failed"
            return
        registry["by_guide"].setdefault(guide, []).append(captured)

    def continuation_guide_origin(self, rule, target, fallback, guide, outer_rule, outer_slot):
        """Opt-in exact-guide observation; unknown never affects source selection."""
        result = {
            "classification": "unknown",
            "reason": None,
            "source_slot": None,
            "raw_alternatives": 0,
            "full_role_on_normalized_guide_same_site": False,
            "concrete_leaf_requirement_proven": False,
        }

        def unknown(reason):
            result["reason"] = reason
            return result

        if not self._diagnose_continuation:
            return unknown("diagnostic_disabled")
        if self.deadline is not None and time.monotonic() >= self.deadline:
            return unknown("deadline")
        cache_key = (rule.id, target, fallback)
        context_key = (cache_key, guide, outer_rule.id, outer_slot)
        if context_key in self._continuation_context_cache:
            return dict(self._continuation_context_cache[context_key])
        if self._continuation_role_contexts >= MAX_CONTINUATION_ROLE_CONTEXTS:
            return unknown("role_context_budget")
        self._continuation_role_contexts += 1
        registry = self._continuation_raw.get(cache_key)
        if registry is None:
            return unknown("raw_provenance_missing_or_cached_before_capture")
        if registry["incomplete"] is not None:
            return unknown(registry["incomplete"])
        if registry["raw_cap_saturated"]:
            return unknown("raw_outcome_cap_saturated")
        alternatives = registry["by_guide"].get(guide, ())
        result["raw_alternatives"] = len(alternatives)
        if not alternatives:
            return unknown("guide_has_no_raw_provenance")
        try:
            parent = molecule(target)
            role = outer_rule.reaction.GetReactants()[outer_slot]
            owners = set()
            for raw, clean in alternatives:
                remaining = MAX_CONTINUATION_SNAPSHOT_UNITS - self._continuation_snapshot_units
                if remaining <= 0:
                    return unknown("structural_work_budget")
                origin = continuation_role_origin(
                    parent,
                    raw,
                    role,
                    deadline_monotonic=self.deadline,
                    max_origin_visits=remaining,
                )
                self._continuation_snapshot_units += origin["origin_visits"]
                if origin["classification"] != "single_source_local_intact":
                    return unknown(f"raw_{origin['classification']}_{origin['reason'] or 'observed'}")
                slot = origin["source_slot"]
                if slot is None or slot >= len(clean):
                    return unknown("source_slot_missing")
                normalized = clean[slot]
                if normalized.GetNumAtoms() > MAX_CONTINUATION_ROLE_ATOMS:
                    return unknown("normalized_part_size_limit")
                normalized_units = (
                    (normalized.GetNumAtoms() + raw[slot].GetNumAtoms())
                    * (len(origin["matched_parent_atoms"]) + 1)
                    + normalized.GetNumBonds()
                    + raw[slot].GetNumBonds()
                    + role.GetNumAtoms()
                    + role.GetNumBonds()
                )
                if self._continuation_snapshot_units + normalized_units > MAX_CONTINUATION_SNAPSHOT_UNITS:
                    return unknown("structural_work_budget")
                self._continuation_snapshot_units += normalized_units
                positions = []
                for target_index in origin["matched_parent_atoms"]:
                    found = [
                        atom.GetIdx()
                        for atom in normalized.GetAtoms()
                        if atom.HasProp("react_atom_idx")
                        and atom.HasProp("react_idx")
                        and atom.GetIntProp("react_atom_idx") == target_index
                        and atom.GetIntProp("react_idx") == 0
                    ]
                    if len(found) != 1:
                        return unknown("normalized_origin_missing_or_ambiguous")
                    positions.append(found[0])
                if self.deadline is not None and time.monotonic() >= self.deadline:
                    return unknown("deadline")
                matches = normalized.GetSubstructMatches(
                    role, useChirality=True, uniquify=False, maxMatches=MAX_CONTINUATION_ROLE_MATCHES
                )
                if len(matches) >= MAX_CONTINUATION_ROLE_MATCHES:
                    return unknown("normalized_role_match_truncated")
                if len(matches) != 1:
                    return unknown("normalized_role_embedding_ambiguous_or_absent")
                if tuple(positions) != matches[0]:
                    return unknown("normalized_same_site_role_mismatch")
                if not _normalized_origin_site_unchanged(
                    raw[slot], normalized, origin["matched_parent_atoms"]
                ):
                    return unknown("normalization_edited_role_site")
                owners.add(slot)
            if len(owners) != 1:
                return unknown("raw_alternative_owner_conflict")
            if self.deadline is not None and time.monotonic() >= self.deadline:
                return unknown("deadline")
            result["classification"] = "single_source_local_intact_normalized_role"
            result["source_slot"] = owners.pop()
            result["full_role_on_normalized_guide_same_site"] = True
            self._continuation_context_cache[context_key] = dict(result)
            return result
        except (ValueError, RuntimeError, AttributeError, IndexError, KeyError):
            return unknown("role_or_origin_check_failed")

    def _reverse(self, rule, target, fallback, *, record_centers=False):
        key = (rule.id, target, fallback)
        if key in self._reverse_cache:
            if record_centers:
                # Metadata is optional. A guide already cached by an inner
                # call cannot justify rerunning reverse chemistry; retain it
                # with explicitly unknown cut provenance instead.
                for guide in self._reverse_cache[key]:
                    self._guide_centers.setdefault((rule.id, target, guide), (None,))
            return self._reverse_cache[key]
        guides = []
        qualities = {}
        centers = {}
        representative_counts = {}
        representative_unknown = {}
        raw_seen = False
        reverse_failed = False
        nonring_query = False
        registry = (
            {"by_guide": {}, "incomplete": None, "raw_cap_saturated": False}
            if self._diagnose_continuation
            and (self._phase == 2 or self._continuation_capture_all_phases)
            and len(self._continuation_raw) < MAX_CONTINUATION_REGISTRIES
            else None
        )
        if not fallback and not self._stop():
            self.report["reverse_attempts"] += 1
            try:
                restriction = None
                if rule.id not in self._nonring_unavailable_rules:
                    specification, reason = self._topology_spec(rule)
                    if reason is None and len(specification[2]) == 1:
                        restriction = specification[2][0]
                if self._stop():
                    return []
                template_key = (rule.id, restriction is not None)
                if template_key not in self._tagged_reverse_cache:
                    new_template = True
                    reverse, tags = tagged_reverse_reaction(rule)
                    if restriction is not None:
                        try:
                            bond = reverse.GetReactants()[0].GetBondBetweenAtoms(*restriction)
                            if bond is None:
                                raise ValueError("Missing local reverse cross bond")
                            bond.ExpandQuery(Chem.BondFromSmarts("!@"), Chem.CompositeQueryType.COMPOSITE_AND)
                            with rdBase.BlockLogs():
                                reverse.Initialize()
                        except (ValueError, RuntimeError, AttributeError):
                            # Failure to build an optional local query must
                            # not delete the original reverse search path.
                            self.report["reverse_nonring_compile_failures"] += 1
                            self._nonring_unavailable_rules.add(rule.id)
                            restriction = None
                            template_key = (rule.id, False)
                            if template_key in self._tagged_reverse_cache:
                                reverse, tags = self._tagged_reverse_cache[template_key]
                                new_template = False
                            else:
                                reverse, tags = tagged_reverse_reaction(rule)
                        else:
                            self.report["reverse_nonring_query_compiles"] += 1
                    if new_template:
                        self._tagged_reverse_cache[template_key] = reverse, tags
                        self.report["reverse_local_template_compiles"] += 1
                reverse, tags = self._tagged_reverse_cache[template_key]
                nonring_query = restriction is not None
                parent = molecule(target)
                if self._stop():
                    # Static support checks, local query compilation and
                    # parent parsing can consume the remaining deadline.
                    # No chemistry call or fallback may start after expiry.
                    return []
                if nonring_query:
                    self.report["reverse_nonring_query_calls"] += 1
                with rdBase.BlockLogs():
                    outcomes = reverse.RunReactants((parent,), maxProducts=self.limits.max_reverse_outcomes)
                if nonring_query and not outcomes:
                    self.report["reverse_nonring_empty"] += 1
                raw_seen = bool(outcomes)
                if registry is not None:
                    registry["raw_cap_saturated"] = len(outcomes) >= self.limits.max_reverse_outcomes
                for outcome in outcomes:
                    if self._stop():
                        if registry is not None:
                            registry["incomplete"] = "reverse_interrupted"
                        break
                    if len(outcome) != len(rule.roles):
                        self._role_rejection("wrong_precursor_count")
                        continue
                    self.report["reverse_partition_checks"] += 1
                    issue = reverse_partition_issue(outcome, parent.GetNumAtoms())
                    if issue is not None:
                        counts = self.report["reverse_partition_rejections"]
                        counts[issue] = counts.get(issue, 0) + 1
                        continue
                    cut_bonds = self._reverse_cut_bonds(outcome, parent) if record_centers else None
                    clean_parts = (
                        []
                        if registry is not None
                        or self._preserve_typed_prefix
                        or self._correct_typed_prefix_size
                        else None
                    )
                    partial_typed = [] if self._preserve_typed_prefix else None
                    typed = self._typed_reverse_outcome(
                        outcome,
                        parent,
                        rule,
                        tags,
                        capture_parts=clean_parts,
                        capture_partial=partial_typed,
                    )
                    if self._stop():
                        if registry is not None:
                            registry["incomplete"] = "reverse_interrupted"
                        break
                    if typed is not None:
                        guide, changes = typed
                        quality = "typed"
                        typed_slots = len(guide)
                    else:
                        use_prefix = (
                            bool(partial_typed)
                            and clean_parts is not None
                            and len(clean_parts) == len(partial_typed)
                            and all(part is not None for part in clean_parts)
                        )
                        if clean_parts is not None and not use_prefix:
                            clean_parts.clear()
                        typed_slots = len(partial_typed) if use_prefix else 0
                        # A reverse guide is a search hint. Keep the old safe
                        # approximation when exact role reconstruction fails.
                        # This never bypasses the parent-origin partition check.
                        try:
                            parts, changes = [], {}
                            if use_prefix:
                                for part, edits in partial_typed:
                                    parts.append(part)
                                    for label, count in edits.items():
                                        changes[label] = changes.get(label, 0) + count
                            for mol in outcome[len(parts) :]:
                                captured = []
                                part, edits = normalize_reverse_hint(
                                    mol,
                                    parent,
                                    rule,
                                    capture_clean=(lambda clean, smiles: captured.append((clean, smiles)))
                                    if clean_parts is not None
                                    else None,
                                )
                                parts.append(part)
                                if clean_parts is not None:
                                    clean_parts.append(
                                        captured[0][0]
                                        if len(captured) == 1 and captured[0][1] == part
                                        else None
                                    )
                                for label, count in edits.items():
                                    changes[label] = changes.get(label, 0) + count
                            guide = tuple(parts)
                            quality = "approximate"
                            if use_prefix:
                                self.report["mixed_typed_prefix_hints"] += 1
                                self.report["mixed_typed_prefix_parts"] += len(partial_typed)
                        except ReverseHintError as exc:
                            counts = self.report["reverse_hint_rejections"]
                            counts[exc.reason] = counts.get(exc.reason, 0) + 1
                            self.report["reverse_failures"] += 1
                            continue
                        except (ValueError, RuntimeError):
                            self.report["reverse_failures"] += 1
                            continue
                    guides.append(guide)
                    if self._correct_typed_prefix_size and record_centers:
                        if guide not in representative_counts and guide not in representative_unknown:
                            if (
                                self.report["representative_origin_guide_entries"]
                                >= MAX_REPRESENTATIVE_ORIGIN_GUIDES
                            ):
                                representative_unknown[guide] = "guide_entry_limit"
                            else:
                                self.report["representative_origin_guide_entries"] += 1
                        if guide not in representative_unknown:
                            count, reason = self._representative_origin_count(
                                outcome, clean_parts, typed_slots, parent.GetNumAtoms()
                            )
                            self._merge_representative_origin(
                                representative_counts, representative_unknown, guide, count, reason
                            )
                    if registry is not None:
                        try:
                            self._capture_continuation_raw(registry, guide, outcome, clean_parts)
                        except Exception:
                            registry["incomplete"] = "snapshot_unexpected_error"
                    if record_centers:
                        centers.setdefault(guide, set()).add(cut_bonds)
                    if quality == "typed" or guide not in qualities:
                        qualities[guide] = quality
                    for label, count in changes.items():
                        if count:
                            edits = self.report["reverse_hint_normalizations"]
                            edits[label] = edits.get(label, 0) + count
            except (ValueError, RuntimeError, AttributeError):
                reverse_failed = True
                self.report["reverse_failures"] += 1
                if registry is not None:
                    registry["incomplete"] = "reverse_failed"
        if (
            not guides
            and not fallback
            and not self._stop()
            and not raw_seen
            and not reverse_failed
            and not nonring_query
            and any(
                molecule(target).HasSubstructMatch(q, useChirality=True) for q in rule.reaction.GetReactants()
            )
        ):
            fallback = True
            self.report["reverse_fallbacks"] += 1
            implicit = True
        else:
            implicit = False
        if fallback:
            # Coarse, explicitly capped fallback. Unmatched roles use a small
            # retrieval query so a pair of nearly whole-hit blocks is not forced.
            mol = molecule(target)
            role_matches = tuple(
                mol.HasSubstructMatch(q, useChirality=True) for q in rule.reaction.GetReactants()
            )
            guides = [tuple(target if matched else "C" for matched in role_matches)]
            qualities = {guides[0]: "implicit_fallback" if implicit else "explicit_fallback"}
            if not implicit:
                self._guide_missing_role_slots[(rule.id, target, guides[0], True)] = tuple(
                    slot for slot, matched in enumerate(role_matches) if not matched
                )
            if record_centers:
                centers = {guides[0]: {None}}
            if registry is not None:
                registry["incomplete"] = "fallback_has_no_raw_origin"
        guides = sorted(set(guides))[: self.limits.max_reverse_outcomes]
        self.report["reverse_outcomes"] += len(guides)
        for guide in guides:
            quality = qualities[guide]
            self._guide_quality[(rule.id, target, guide)] = quality
            if record_centers:
                self._guide_centers[(rule.id, target, guide)] = tuple(
                    sorted(centers[guide], key=lambda value: (value is None, value or ()))
                )
            self.report["reverse_hint_quality_counts"][quality] += 1
            if self._correct_typed_prefix_size and record_centers:
                guide_key = (rule.id, target, guide)
                if guide in representative_counts and guide not in representative_unknown:
                    self._guide_representative_heavy[guide_key] = representative_counts[guide]
                    self.report["representative_origin_known_guides"] += 1
                else:
                    reason = representative_unknown.get(guide, "no_raw_origin")
                    reasons = self.report["representative_origin_legacy_reasons"]
                    reasons[reason] = reasons.get(reason, 0) + 1
        self._reverse_cache[key] = guides
        if registry is not None:
            self._continuation_raw[key] = registry
        return guides

    def _retrieval_target(self, target):
        """Erase reverse-SMARTS dummy/radical artifacts for retrieval only.

        The original rule, catalog block and forward replay are never changed.
        An unconstructive guide falls back to its original literal query.
        """
        if target in self._retrieval_target_cache:
            return self._retrieval_target_cache[target]
        try:
            source = molecule(target)
            edited = Chem.RWMol(source)
            changed = False
            for atom in edited.GetAtoms():
                if atom.GetNumRadicalElectrons():
                    atom.SetNumRadicalElectrons(0)
                    atom.SetNumExplicitHs(0)
                    atom.SetNoImplicit(False)
                    changed = True
            for index in reversed(range(edited.GetNumAtoms())):
                if edited.GetAtomWithIdx(index).GetAtomicNum() == 0:
                    edited.RemoveAtom(index)
                    changed = True
            if edited.GetNumAtoms() == 0:
                raise ValueError("Guide has no concrete atoms")
            clean = edited.GetMol()
            Chem.SanitizeMol(clean)
            if len(Chem.GetMolFrags(clean)) != 1:
                raise ValueError("Guide is disconnected after placeholder removal")
            result = canonical(Chem.MolToSmiles(clean, isomericSmiles=True))
            if changed:
                self.report["retrieval_guides_normalized"] += 1
        except (ValueError, RuntimeError):
            self.report["retrieval_guides_unusable"] += 1
            result = target
        self._retrieval_target_cache[target] = result
        return result

    def _outer_size_delta(self, rule, guide, normalized, artifacts):
        normalized_heavy = sum(molecule(part).GetNumHeavyAtoms() for part in normalized)
        representative_count = (
            self._guide_representative_heavy.get((rule.id, self.hit, guide))
            if self._correct_typed_prefix_size and not artifacts
            else None
        )
        if representative_count:
            self.report["representative_origin_score_corrections"] += 1
        elif self._correct_typed_prefix_size and artifacts:
            reasons = self.report["representative_origin_legacy_reasons"]
            reasons["retrieval_target_changed"] = reasons.get("retrieval_target_changed", 0) + 1
        return abs(normalized_heavy - (representative_count or 0) - self.hit_mol.GetNumHeavyAtoms())

    def _guide_anchors(self, target):
        if target in self._guide_anchor_cache:
            return self._guide_anchor_cache[target]
        matched = []
        if self._stop():
            return ()
        mol = molecule(target)
        for index, query in enumerate(self._anchor_queries):
            if self._stop():
                break
            if self.report["guide_anchor_checks"] >= 2048:
                self._stop("guide_anchor_check_budget")
                break
            self.report["guide_anchor_checks"] += 1
            if mol.HasSubstructMatch(query, useChirality=True):
                matched.append(index)
        result = tuple(matched)
        self._guide_anchor_cache[target] = result
        return result

    def _catalog_block_stereo_eligible(self, smiles):
        """A strict route requires every concrete source reactant to be specified."""
        if smiles in self._catalog_stereo_cache:
            self.report["catalog_stereo_cache_hits"] += 1
            return self._catalog_stereo_cache[smiles]
        if self._stop():
            reason = "deadline"
        elif self.report["catalog_stereo_perceptions"] >= MAX_CATALOG_STEREO_PERCEPTIONS:
            self._stop("max_catalog_stereo_perceptions")
            reason = "perception_budget"
        elif len(self._catalog_stereo_cache) >= MAX_CATALOG_STEREO_CACHE_ENTRIES:
            self._stop("max_catalog_stereo_cache_entries")
            reason = "cache_budget"
        else:
            mol = self.catalog.mols[smiles]
            cost = mol.GetNumAtoms() + mol.GetNumBonds()
            if self.report["catalog_stereo_atom_bond_charges"] + cost > MAX_CATALOG_STEREO_ATOM_BOND_CHARGES:
                self._stop("max_catalog_stereo_atom_bond_charges")
                reason = "atom_bond_budget"
            else:
                self.report["catalog_stereo_perceptions"] += 1
                self.report["catalog_stereo_atom_bond_charges"] += cost
                try:
                    eligible = stereo_state(mol)["fully_specified"]
                except (ValueError, RuntimeError, AttributeError):
                    reason = "perception_failed"
                else:
                    if self._stop():
                        reason = "deadline"
                    else:
                        self._catalog_stereo_cache[smiles] = eligible
                        return eligible
        counts = self.report["catalog_stereo_unknowns"]
        counts[reason] = counts.get(reason, 0) + 1
        return None  # Unknown remains eligible at its legacy rank.

    def _continuation_second(self, base_key, legacy, context):
        """Prefer a role-matching ordinary second from the already ranked pool."""
        if context is None or not self._prefer_continuation_second or self._phase != 2:
            return legacy
        outer_rule, outer_slot, source_slot, inner_rule_id, inner_target, inner_guide = context
        key = (base_key, outer_rule.id, outer_slot, inner_rule_id, inner_target, inner_guide, source_slot)
        report = self.report["continuation_preference"]

        def old(reason):
            counts = report["legacy_reasons"]
            counts[reason] = counts.get(reason, 0) + 1
            return legacy

        if self.deadline is not None and time.monotonic() >= self.deadline:
            return old("deadline")
        if key in self._continuation_selected:
            return self._continuation_selected[key]
        if len(self._continuation_selected) >= MAX_CONTINUATION_SELECTED_CONTEXTS:
            return old("selected_context_cache_budget")
        stored = self._continuation_ranked.get(base_key)
        if stored is None:
            return old("complete_ranked_pool_unavailable")
        ordinary, ordinary_second = stored
        if self.limits.building_block_top_k != 2 or len(legacy) != 2 or not ordinary_second:
            return old("no_ordinary_second")
        role_key = (outer_rule.id, outer_slot)
        query = self._continuation_outer_roles.get(role_key)
        if query is None:
            if len(self._continuation_outer_roles) >= MAX_CONTINUATION_OUTER_ROLE_QUERIES:
                return old("outer_role_query_budget")
            try:
                query = outer_rule.reaction.GetReactants()[outer_slot]
                if (
                    query.GetNumAtoms() > MAX_CONTINUATION_ROLE_ATOMS
                    or query.GetNumBonds() > MAX_CONTINUATION_ROLE_BONDS
                ):
                    return old("outer_role_query_size")
            except (IndexError, AttributeError):
                return old("outer_role_query_unavailable")
            if self.deadline is not None and time.monotonic() >= self.deadline:
                return old("deadline")
            self._continuation_outer_roles[role_key] = query
            report["outer_role_queries"] += 1
        checked = 0
        for row in ordinary:
            smiles = row[2]
            if smiles == legacy[0].smiles:
                continue
            if self.deadline is not None and time.monotonic() >= self.deadline:
                return old("deadline")
            if checked >= MAX_CONTINUATION_ROLE_CHECKS_PER_QUERY:
                return old("role_check_query_budget")
            if report["role_checks"] >= MAX_CONTINUATION_ROLE_CHECKS_PER_SEARCH:
                return old("role_check_search_budget")
            checked += 1
            report["role_checks"] += 1
            try:
                if not self.catalog.mols[smiles].HasSubstructMatch(query, useChirality=True):
                    continue
            except (ValueError, RuntimeError, KeyError):
                return old("outer_role_match_failed")
            report["role_matches"] += 1
            if smiles not in self._catalog_stereo_cache:
                mol = self.catalog.mols[smiles]
                cost = mol.GetNumAtoms() + mol.GetNumBonds()
                if (
                    self.report["catalog_stereo_perceptions"] >= MAX_CATALOG_STEREO_PERCEPTIONS
                    or len(self._catalog_stereo_cache) >= MAX_CATALOG_STEREO_CACHE_ENTRIES
                    or self.report["catalog_stereo_atom_bond_charges"] + cost
                    > MAX_CATALOG_STEREO_ATOM_BOND_CHARGES
                ):
                    return old("stereo_budget")
            eligible = self._catalog_block_stereo_eligible(smiles)
            if eligible is None:
                report["stereo_unknowns"] += 1
                return old("stereo_unknown")
            if not eligible:
                report["stereo_rejections"] += 1
                continue
            if self.deadline is not None and time.monotonic() >= self.deadline:
                return old("deadline")
            replacement = self.catalog.blocks[sorted(self.catalog.by_smiles[smiles])[0]]
            selected = (legacy[0], replacement)
            if selected != legacy:
                report["replacements"] += 1
                if len(report["samples"]) < MAX_FALLBACK_HIT_SECOND_SAMPLES:
                    report["samples"].append(
                        {
                            "outer_rule_id": outer_rule.id,
                            "outer_slot": outer_slot,
                            "inner_rule_id": inner_rule_id,
                            "inner_target": inner_target,
                            "source_slot": source_slot,
                            "first": legacy[0].smiles,
                            "old_second": legacy[1].smiles,
                            "new_second": smiles,
                            "role_checks": checked,
                        }
                    )
                else:
                    report["samples_omitted"] += 1
            self._continuation_selected[key] = selected
            return selected
        if self.deadline is not None and time.monotonic() >= self.deadline:
            return old("deadline")
        self._continuation_selected[key] = legacy
        return old("no_role_and_stereo_eligible_ordinary")

    def _nearest(
        self,
        target,
        role,
        *,
        rule_id=None,
        role_slot=None,
        fallback_missing_role=False,
        continuation_context=None,
    ):
        # This context changes only the second block. Keep its completed
        # selection separate from the ordinary nearest result for this query.
        key = (target, role, fallback_missing_role)
        trace = self._active_phase2_slot
        if trace is not None:
            trace["pool_queries"] += 1
        if key in self._nearest_cache:
            if trace is not None:
                trace["pool_cache_hits"] += 1
                trace["pool_selected_blocks"] += len(self._nearest_cache[key])
            return self._continuation_second(key, self._nearest_cache[key], continuation_context)
        self.report["role_queries"] += 1
        query = Chem.MolFromSmarts(role)
        retrieval_target = self._retrieval_target(target)
        # Radical hydrogens in reversed templates can retain useful local
        # positional features in the fingerprint. Dummy atoms cannot, so only
        # those use the normalized fingerprint; both remain retrieval hints.
        fingerprint_target = (
            retrieval_target
            if any(atom.GetAtomicNum() == 0 for atom in molecule(target).GetAtoms())
            else target
        )
        fp = self._fingerprint(fingerprint_target)
        if fp is None:
            return ()
        ranked = []
        pool = self._pool(retrieval_target, guide=target, role=role)
        if trace is not None:
            trace["pool_bucket_candidates"] += len(pool)
        active_anchors = self._guide_anchors(target)
        anchor_strength = {}
        for index in active_anchors:
            for smiles in self._anchor_postings[index]:
                anchor_strength[smiles] = anchor_strength.get(smiles, 0) + 1
        complete = self._pool_last_complete
        phase_limit = (
            self.limits.max_catalog_row_visits // 2
            if self._phase == 1 and self.limits.max_steps == 2
            else self.limits.max_catalog_row_visits
        )
        if self.report["catalog_row_visits"] + len(anchor_strength) > phase_limit:
            self._stop("max_catalog_row_visits")
            return ()
        self.report["catalog_row_visits"] += len(anchor_strength)
        self.report["anchor_posting_metadata_visits"] += len(anchor_strength)
        ordered_anchors = sorted(
            anchor_strength,
            key=lambda smiles: (-anchor_strength[smiles], self._anchor_scores[smiles], smiles),
        )
        chosen_anchors = []
        posting_attempts = 0
        posting_limit = min(256, self.limits.max_catalog_rows_per_query)
        for smiles in ordered_anchors[:posting_limit]:
            if self._stop():
                complete = False
                break
            if self.report["catalog_row_visits"] >= phase_limit:
                self._stop("max_catalog_row_visits")
                complete = False
                break
            self.report["catalog_row_visits"] += 1
            self.report["anchor_role_checks"] += 1
            if trace is not None:
                trace["pool_role_checks"] += 1
            posting_attempts += 1
            mol = self.catalog.mols[smiles]
            if mol.HasSubstructMatch(query, useChirality=True):
                self.report["anchor_role_compatible"] += 1
                chosen_anchors.append(smiles)
        remaining = self.limits.max_catalog_rows_per_query - posting_attempts
        base_pool = tuple(smiles for smiles in pool if smiles not in chosen_anchors)[:remaining]
        if target in self.catalog.mols and target not in chosen_anchors and target not in base_pool:
            base_pool = (target,) + base_pool[: max(0, remaining - 1)] if remaining else ()
        for smiles in itertools.chain(chosen_anchors, base_pool):
            if self._stop():
                complete = False
                break
            if smiles not in chosen_anchors:
                if self.report["catalog_row_visits"] >= phase_limit:
                    self._stop("max_catalog_row_visits")
                    complete = False
                    break
                self.report["catalog_row_visits"] += 1
                mol = self.catalog.mols[smiles]
                if trace is not None:
                    trace["pool_role_checks"] += 1
                if not mol.HasSubstructMatch(query, useChirality=True):
                    continue
            if trace is not None:
                trace["pool_role_compatible"] += 1
            candidate_fp = self._fingerprint(smiles)
            if candidate_fp is None:
                complete = False
                break
            score = DataStructs.TanimotoSimilarity(fp, candidate_fp)
            ranked.append((-score, -anchor_strength.get(smiles, 0), smiles))
        anchored = sorted((row for row in ranked if row[1] < 0), key=lambda row: (row[1], row[0], row[2]))
        ordinary = sorted((row for row in ranked if row[1] == 0), key=lambda row: (row[0], row[2]))
        selected = []
        stereo_complete = True

        def eligible_row(row):
            nonlocal stereo_complete
            smiles = row[2]
            eligible = self._catalog_block_stereo_eligible(smiles)
            if eligible is False:
                self.report["catalog_stereo_exclusions"] += 1
                samples = self.report["catalog_stereo_exclusion_samples"]
                if len(samples) < MAX_CATALOG_STEREO_EXCLUSION_SAMPLES:
                    samples.append(
                        {
                            "building_block_id": sorted(self.catalog.by_smiles[smiles])[0],
                            "smiles": smiles,
                            "target": target,
                            "rule_id": rule_id,
                            "role_slot": role_slot,
                            "selection_phase": self._phase,
                            "reason": "strict_route_reactant_stereo_not_fully_specified",
                        }
                    )
                else:
                    self.report["catalog_stereo_exclusion_samples_truncated"] += 1
                return False
            if eligible is None:
                stereo_complete = False
            return True

        def first_eligible(rows):
            for index, row in enumerate(rows):
                if eligible_row(row):
                    return row, rows[index + 1 :]
            return None, []

        # Preserve the original one-anchor/one-ordinary reservation after
        # screening. A rejected ordinary leader must not be replaced by a
        # second anchor before the next ordinary candidate is considered.
        anchored_tail, ordinary_tail = anchored, ordinary
        if anchored:
            first, anchored_tail = first_eligible(anchored)
            if first is not None:
                selected.append(first)
        if ordinary and len(selected) < self.limits.building_block_top_k:
            first, ordinary_tail = first_eligible(ordinary)
            if first is not None:
                selected.append(first)
        for row in anchored_tail + ordinary_tail:
            if len(selected) >= self.limits.building_block_top_k:
                break
            if eligible_row(row):
                selected.append(row)
        hit_second_complete = True
        if fallback_missing_role and self.limits.building_block_top_k == 2:
            self.report["fallback_hit_second_queries"] += 1

            def use_legacy(reason):
                nonlocal hit_second_complete
                hit_second_complete = False
                counts = self.report["fallback_hit_second_legacy_reasons"]
                counts[reason] = counts.get(reason, 0) + 1

            if not complete or not stereo_complete or self._stop():
                use_legacy("incomplete_original_selection")
            elif len(selected) < 2 or selected[1][1] != 0:
                # Preserve the existing anchor/ordinary reservation and first.
                use_legacy("no_ordinary_second")
            else:
                candidates = [(rank, row) for rank, row in enumerate(ordinary) if row[2] != selected[0][2]]
                remaining = MAX_FALLBACK_HIT_SECOND_TANIMOTOS - self.report["fallback_hit_second_tanimotos"]
                if len(candidates) > self.limits.max_catalog_rows_per_query:
                    use_legacy("ranked_query_cap")
                elif len(candidates) > remaining:
                    use_legacy("tanimoto_budget")
                elif any(row[2] not in self._fp_cache for _rank, row in candidates):
                    use_legacy("missing_cached_fingerprint")
                else:
                    scored = []
                    for rank, row in candidates:
                        if self._stop():
                            use_legacy("deadline_during_scoring")
                            break
                        self.report["fallback_hit_second_tanimotos"] += 1
                        score = DataStructs.TanimotoSimilarity(self.hit_fp, self._fp_cache[row[2]])
                        scored.append((-score, rank, row[2], row))
                    if hit_second_complete:
                        hit_scores = {smiles: -score for score, _rank, smiles, _row in scored}
                        for _score, _rank, smiles, row in sorted(scored):
                            if self._stop():
                                use_legacy("deadline_during_stereo_screen")
                                break
                            if smiles not in self._catalog_stereo_cache:
                                mol = self.catalog.mols[smiles]
                                cost = mol.GetNumAtoms() + mol.GetNumBonds()
                                if (
                                    self.report["catalog_stereo_perceptions"]
                                    >= MAX_CATALOG_STEREO_PERCEPTIONS
                                    or len(self._catalog_stereo_cache) >= MAX_CATALOG_STEREO_CACHE_ENTRIES
                                    or self.report["catalog_stereo_atom_bond_charges"] + cost
                                    > MAX_CATALOG_STEREO_ATOM_BOND_CHARGES
                                ):
                                    use_legacy("stereo_budget")
                                    break
                            eligible = eligible_row(row)
                            if not stereo_complete:
                                use_legacy("unknown_source_stereo")
                                break
                            if eligible:
                                if self._stop():
                                    use_legacy("deadline_after_stereo_screen")
                                    break
                                if row != selected[1]:
                                    self.report["fallback_hit_second_replacements"] += 1
                                    samples = self.report["fallback_hit_second_samples"]
                                    if len(samples) < MAX_FALLBACK_HIT_SECOND_SAMPLES:
                                        samples.append(
                                            {
                                                "rule_id": rule_id,
                                                "role_slot": role_slot,
                                                "placeholder_target": target,
                                                "first": selected[0][2],
                                                "old_second": selected[1][2],
                                                "new_second": smiles,
                                                "old_second_hit_tanimoto": hit_scores[selected[1][2]],
                                                "new_second_hit_tanimoto": hit_scores[smiles],
                                                "comparison_count": len(scored),
                                            }
                                        )
                                    else:
                                        self.report["fallback_hit_second_samples_omitted"] += 1
                                selected[1] = row
                                break
        blocks = tuple(
            self.catalog.blocks[sorted(self.catalog.by_smiles[smiles])[0]]
            for _score, _anchors, smiles in selected
        )
        if trace is not None:
            trace["pool_selected_blocks"] += len(blocks)
        if complete and stereo_complete and hit_second_complete:
            self._nearest_cache[key] = blocks
            if continuation_context is not None and self._prefer_continuation_second and self._phase == 2:
                report = self.report["continuation_preference"]
                if (
                    len(ordinary) <= self.limits.max_catalog_rows_per_query
                    and report["ranked_rows_retained"] + len(ordinary) <= MAX_CONTINUATION_RANKED_ROWS
                ):
                    self._continuation_ranked[key] = (
                        tuple(ordinary),
                        len(selected) >= 2 and selected[1][1] == 0,
                    )
                    report["ranked_rows_retained"] += len(ordinary)
        return self._continuation_second(key, blocks, continuation_context)

    def _stage(
        self,
        rule,
        guide,
        prefix,
        replace=None,
        budget=None,
        product_budget=None,
        fallback_missing_slots=(),
        continuation_context=None,
    ):
        options = []
        for slot, (target, role) in enumerate(zip(guide, rule.roles)):
            if replace is not None and replace[0] == slot:
                options.append(replace[1])
            else:
                options.append(
                    tuple(
                        {
                            "root": f"{prefix}b{slot}",
                            "product": block.smiles,
                            "nodes": [
                                {
                                    "id": f"{prefix}b{slot}",
                                    "building_block_id": block.id,
                                }
                            ],
                        }
                        for block in self._nearest(
                            target,
                            role,
                            rule_id=rule.id,
                            role_slot=slot,
                            fallback_missing_role=slot in fallback_missing_slots,
                            continuation_context=(
                                continuation_context if continuation_context[2] == slot else None
                            )
                            if continuation_context is not None
                            else None,
                        )
                    )
                )
        states = []
        for combination in itertools.product(*options):
            if self._stop():
                break
            if self.report["combinations"] >= self.limits.max_combinations or (
                budget is not None and self.report["combinations"] >= budget
            ):
                self._stop("max_combinations" if budget is None else "depth_one_combination_budget")
                break
            self.report["combinations"] += 1
            reactants = tuple(item["product"] for item in combination)
            products, capped = self.rules.apply(rule.id, reactants, self.limits.reaction_product_limit)
            if capped:
                self._stop("reaction_product_limit")
            for product in products:
                if self.report["products_generated"] >= (
                    product_budget if product_budget is not None else self.limits.max_products
                ):
                    self._stop("depth_one_product_budget" if product_budget is not None else "max_products")
                    break
                self.report["products_generated"] += 1
                nodes = [node for item in combination for node in item["nodes"]]
                nodes.append(
                    {
                        "id": f"{prefix}r",
                        "rule_id": rule.id,
                        "reactants": [item["root"] for item in combination],
                        "product": product,
                    }
                )
                states.append({"root": f"{prefix}r", "product": product, "nodes": nodes})
            if self.report["products_generated"] >= (
                product_budget if product_budget is not None else self.limits.max_products
            ):
                break
        return states

    def _compatible_inner(self, states, outer_rule, slot):
        role = outer_rule.reaction.GetReactants()[slot]
        compatible = []
        for state in states:
            if self._stop():
                break
            self.report["inner_role_checks"] += 1
            if molecule(state["product"]).HasSubstructMatch(role, useChirality=True):
                self.report["inner_role_compatible"] += 1
                compatible.append(state)
            else:
                self.report["inner_role_rejected"] += 1
        return compatible

    def _retain(self, states):
        # Keep the bounded forward outcomes until both search phases finish.
        # Product admission and route audit then spend their budget on the
        # closest eligible exact products, rather than the first rule's output.
        for state in states:
            if self._stop():
                return
            self._pending.append(state)

    def _admit(self):
        unique_products = sorted({state["product"] for state in self._pending})
        self.report["admission_pending_unique_products"] = len(unique_products)
        self.report["admission_fp_unassessed_products"] = len(unique_products)
        if self._stop() or not unique_products:
            return
        available = max(
            0,
            self.limits.max_fingerprint_calls
            - self.report["fingerprint_calls"]
            - getattr(self.accept, "fingerprint_calls", 0),
        )
        product_ranks = {}
        # Select a bounded subset before similarity ordering. Each selected
        # product gets its possible new rank fingerprint and one callback
        # fingerprint reserved; an unranked product has no artificial score.
        for product in sorted(unique_products, key=lambda item: (item not in self._fp_cache, item)):
            if self._stop():
                break
            required = (product not in self._fp_cache) + 1
            if available < required:
                continue
            available -= required
            self.report["admission_reserved_cost_units"] += required
            fp = self._fingerprint(product, purpose="admission")
            if fp is None:
                continue
            product_ranks[product] = -DataStructs.TanimotoSimilarity(self.hit_fp, fp)
        self.report["admission_ranked_products"] = len(product_ranks)
        self.report["admission_fp_unassessed_products"] = len(unique_products) - len(product_ranks)
        if self.report["admission_fp_unassessed_products"] and not self._stop():
            self._stop("admission_fingerprint_reserve_exhausted")
        eligible = set()
        for product in sorted(product_ranks, key=lambda item: (product_ranks[item], item)):
            if self._stop():
                break
            self.report["admission_attempted_products"] += 1
            before_accept = getattr(self.accept, "fingerprint_calls", 0)
            try:
                self.accept(product)
            except (ValueError, KeyError, TypeError) as exc:
                reason = f"{type(exc).__name__}: {exc}"
                self.report["product_rejections"][reason] = (
                    self.report["product_rejections"].get(reason, 0) + 1
                )
                continue
            finally:
                after_accept = getattr(self.accept, "fingerprint_calls", 0)
                self.report["admission_fingerprint_calls"] = after_accept
                if (
                    not 0 <= after_accept - before_accept <= 1
                    or self.report["fingerprint_calls"] + after_accept > self.limits.max_fingerprint_calls
                ):
                    raise RuntimeError("Route-first admission callback exceeded its one-fingerprint contract")
            eligible.add(product)
        self.report["admission_callback_unattempted_products"] = (
            len(product_ranks) - self.report["admission_attempted_products"]
        )
        if self._stop():
            return
        states = sorted(
            (state for state in self._pending if state["product"] in eligible),
            key=lambda state: (
                product_ranks[state["product"]],
                state["product"],
                sum("rule_id" in node for node in state["nodes"]),
                str(state["nodes"]),
            ),
        )
        first, alternatives, seen_products = [], [], set()
        for state in states:
            if state["product"] in seen_products:
                alternatives.append(state)
            else:
                first.append(state)
                seen_products.add(state["product"])
        for state in itertools.chain(first, alternatives):
            if self._stop():
                return
            if self.report["route_audit_attempts"] >= self.limits.max_routes:
                self._stop("max_routes")
                return
            route = {"schema_version": 1, "root": state["root"], "nodes": state["nodes"]}
            key = str(route)
            if key in self._route_keys:
                continue
            self.report["route_audit_attempts"] += 1
            try:
                verification = verify_route(
                    route, self.catalog, self.rules, self.limits.reaction_product_limit
                )
                audit = audit_route(route, self.catalog, self.rules, self.limits.reaction_product_limit)
                if not audit["strict_eligible"]:
                    raise ValueError("Full route atom/stereo audit is not strict eligible")
            except (ValueError, KeyError, TypeError) as exc:
                reason = f"{type(exc).__name__}: {exc}"
                self.report["product_rejections"][reason] = (
                    self.report["product_rejections"].get(reason, 0) + 1
                )
                continue
            self._route_keys.add(key)
            self.routes.append(
                {
                    "route": route,
                    "verification": verification,
                    "generation": {
                        "policy": VERSION,
                        "reference": "self_product_design",
                        "route_audit_binding_sha256": audit["route_binding_sha256"],
                    },
                }
            )
            self.report["routes_verified"] += 1

    def run(self):
        outer = self._select_rules(self.hit, self.limits.max_rules, fallback=self.limits.fallback_rules)
        self.report["rules_selected"] = [
            {"rule_id": rule.id, "fallback": fallback} for rule, fallback in outer
        ]
        depth_one_budget = (
            self.limits.max_combinations
            if self.limits.max_steps == 1
            else max(1, self.limits.max_combinations // 2)
        )
        plan_candidates = []
        for rule, fallback in outer:
            if self._stop():
                break
            for guide in self._reverse(rule, self.hit, fallback, record_centers=True):
                if self._stop():
                    break
                normalized = tuple(self._retrieval_target(part) for part in guide)
                artifacts = sum(a != b for a, b in zip(guide, normalized))
                size_delta = self._outer_size_delta(rule, guide, normalized, artifacts)
                centers = self._guide_centers[(rule.id, self.hit, guide)]
                proven = tuple(center for center in centers if center)
                if proven:
                    center_keys = tuple(("cut", center) for center in proven)
                    center_kind = "cut"
                else:
                    center_kind = "unknown" if None in centers else "no_cut"
                quality = self._guide_quality[(rule.id, self.hit, guide)]
                if quality in ("explicit_fallback", "implicit_fallback"):
                    center_kind = quality
                if not proven:
                    center_keys = ((center_kind, rule.id, guide),)
                plan_candidates.append(
                    (rule, guide, fallback, artifacts, size_delta, center_keys, center_kind, quality)
                )
        self.report["outer_plan_candidates"] = len(plan_candidates)
        self.report["outer_distinct_cut_centers"] = len(
            {
                key
                for _rule, _guide, _fallback, _artifacts, _size_delta, keys, _kind, _quality in plan_candidates
                for key in keys
                if key[0] == "cut"
            }
        )
        plan_candidates.sort(
            key=lambda item: (
                item[2],
                {"typed": 0, "approximate": 1, "implicit_fallback": 2, "explicit_fallback": 2}[item[7]],
                {
                    "cut": 0,
                    "no_cut": 1,
                    "unknown": 2,
                    "explicit_fallback": 3,
                    "implicit_fallback": 3,
                }[item[6]],
                -item[0].reaction.GetProducts()[0].GetNumHeavyAtoms(),
                item[3],
                item[4],
                item[0].id,
                item[1],
            )
        )
        plans = []
        covered_centers = set()
        for only_new_centers in (True, False):
            for rule, guide, fallback, _artifacts, _size_delta, keys, kind, quality in plan_candidates:
                if len(plans) >= self.limits.max_outer_guides:
                    break
                if only_new_centers and all(key in covered_centers for key in keys):
                    continue
                if any(item[0] == rule and item[1] == guide for item in plans):
                    continue
                plans.append((rule, guide, fallback, keys, kind, quality))
                covered_centers.update(keys)
        if len(plan_candidates) > len(plans):
            self._stop("max_outer_guides")
        self.report["outer_plans_selected"] = [
            {
                "rule_id": rule.id,
                "guide": list(guide),
                "fallback": fallback,
                "hint_quality": quality,
                "center_kind": kind,
                "cut_bond_centers": [list(map(list, key[1])) for key in keys if key[0] == "cut"],
            }
            for rule, guide, fallback, keys, kind, quality in plans
        ]
        self._select_role_indexes(plans)
        depth_one_product_budget = (
            self.limits.max_products if self.limits.max_steps == 1 else max(1, self.limits.max_products // 2)
        )
        for rule, guide, fallback, _keys, _kind, quality in plans:
            if self._stop() or self.report["combinations"] >= depth_one_budget:
                break
            states = self._stage(
                rule,
                guide,
                "o",
                budget=depth_one_budget,
                product_budget=depth_one_product_budget,
                fallback_missing_slots=(
                    self._guide_missing_role_slots.get((rule.id, self.hit, guide, True), ())
                    if fallback and quality == "explicit_fallback"
                    else ()
                ),
            )
            self._retain(states)
        if self.limits.max_steps == 2 and not self._stop():
            self._phase = 2
            work_by_plan = []
            plan_traces = []
            for plan_index, (outer_rule, guide, _fallback, _keys, _kind, _quality) in enumerate(plans):
                plan_trace = None
                if len(self.report["phase2_plan_trace"]) < 8:
                    plan_trace = {
                        "plan_index": plan_index,
                        "rule_id": outer_rule.id,
                        "slots": [],
                        "prepared_work_items": 0,
                        "typed_prepared": 0,
                        "approximate_prepared": 0,
                        "fallback_prepared": 0,
                        "work_items_attempted": 0,
                        "typed_attempted": 0,
                        "approximate_attempted": 0,
                        "fallback_attempted": 0,
                        "work_items_deferred": 0,
                        "defer_reason": None,
                        "preparation_limit_observed": None,
                        "rule_inspections_preparation": 0,
                        "reverse_role_work_preparation": 0,
                        "catalog_rows_execution": 0,
                        "catalog_exhausted_attempts": 0,
                    }
                    self.report["phase2_plan_trace"].append(plan_trace)
                else:
                    self.report["phase2_plan_trace_truncated"] += 1
                plan_traces.append(plan_trace)
                work = []
                work_by_plan.append(work)
                before_reverse_role_work = self.report["reverse_role_work"]
                for slot, target in enumerate(guide):
                    if self._stop() or self.report["combinations"] >= self.limits.max_combinations:
                        break
                    slot_trace = None
                    if plan_trace is not None:
                        slot_trace = {
                            "slot": slot,
                            "target": target,
                            "inner_rules_selected": [],
                            "inner_reverse_guides": 0,
                            "inner_stage_calls": 0,
                            "inner_stage_states": 0,
                            "inner_role_compatible_states": 0,
                            "inner_beam_states": 0,
                            "outer_stage_calls": 0,
                            "outer_stage_states": 0,
                            "pool_queries": 0,
                            "pool_cache_hits": 0,
                            "pool_bucket_candidates": 0,
                            "pool_role_checks": 0,
                            "pool_role_compatible": 0,
                            "pool_selected_blocks": 0,
                            "rule_inspections_preparation": 0,
                            "prepared_work_items": 0,
                            "work_items_attempted": 0,
                            "typed_attempted": 0,
                            "approximate_attempted": 0,
                            "fallback_attempted": 0,
                            "catalog_rows_execution": 0,
                        }
                        plan_trace["slots"].append(slot_trace)
                    before_rules = self.report["rule_inspections"]
                    inner_rules = self._select_rules(target, self.limits.max_inner_rules)
                    inspections = self.report["rule_inspections"] - before_rules
                    if plan_trace is not None:
                        plan_trace["rule_inspections_preparation"] += inspections
                        slot_trace["rule_inspections_preparation"] += inspections
                        slot_trace["inner_rules_selected"] = [rule.id for rule, _ in inner_rules]
                    for inner_rule, _ in inner_rules:
                        if self._stop():
                            break
                        inner_guides = self._reverse(inner_rule, target, False)
                        if slot_trace is not None:
                            slot_trace["inner_reverse_guides"] += len(inner_guides)
                        for inner_guide in inner_guides:
                            quality = self._guide_quality[(inner_rule.id, target, inner_guide)]
                            work.append((slot, target, inner_rule, inner_guide, quality, slot_trace))
                            if plan_trace is not None:
                                plan_trace["prepared_work_items"] += 1
                                slot_trace["prepared_work_items"] += 1
                                key = (
                                    f"{quality}_prepared"
                                    if quality in ("typed", "approximate")
                                    else "fallback_prepared"
                                )
                                plan_trace[key] += 1
                if plan_trace is not None:
                    plan_trace["reverse_role_work_preparation"] = (
                        self.report["reverse_role_work"] - before_reverse_role_work
                    )
                    if self._stop():
                        plan_trace["preparation_limit_observed"] = "deadline_monotonic"
                    elif self.report["combinations"] >= self.limits.max_combinations:
                        plan_trace["preparation_limit_observed"] = "max_combinations"
                    elif "max_rule_inspections" in self.report["limits_reached"]:
                        plan_trace["preparation_limit_observed"] = "max_rule_inspections"
                    elif "reverse_role_work_budget" in self.report["limits_reached"]:
                        plan_trace["preparation_limit_observed"] = "reverse_role_work_budget"
                # Python's stable sort preserves original slot/rule/guide
                # order within a quality class. This does not add outcomes.
                work.sort(key=lambda item: item[4] != "typed")

            offsets = [0] * len(plans)
            while any(offset < len(work) for offset, work in zip(offsets, work_by_plan)):
                if self._stop():
                    defer_reason = "deadline_monotonic"
                    break
                if self.report["combinations"] >= self.limits.max_combinations:
                    defer_reason = "max_combinations"
                    break
                for plan_index, (outer_rule, guide, _fallback, _keys, _kind, _quality) in enumerate(plans):
                    if offsets[plan_index] >= len(work_by_plan[plan_index]):
                        continue
                    if self._stop() or self.report["combinations"] >= self.limits.max_combinations:
                        break
                    slot, target, inner_rule, inner_guide, quality, slot_trace = work_by_plan[plan_index][
                        offsets[plan_index]
                    ]
                    offsets[plan_index] += 1
                    plan_trace = plan_traces[plan_index]
                    before_rows = self.report["catalog_row_visits"]
                    if plan_trace is not None:
                        plan_trace["work_items_attempted"] += 1
                        slot_trace["work_items_attempted"] += 1
                        key = (
                            f"{quality}_attempted"
                            if quality in ("typed", "approximate")
                            else "fallback_attempted"
                        )
                        plan_trace[key] += 1
                        slot_trace[key] += 1
                        if before_rows >= self.limits.max_catalog_row_visits:
                            plan_trace["catalog_exhausted_attempts"] += 1
                    self._active_phase2_slot = slot_trace
                    try:
                        source_slot = None
                        if self._prefer_continuation_second:
                            origin = self.continuation_guide_origin(
                                inner_rule,
                                target,
                                False,
                                inner_guide,
                                outer_rule,
                                slot,
                            )
                            if origin["classification"] == "single_source_local_intact_normalized_role":
                                source_slot = origin["source_slot"]
                                self.report["continuation_preference"]["bridge_positive_work_items"] += 1
                            else:
                                report = self.report["continuation_preference"]
                                report["bridge_unknown_work_items"] += 1
                                reason = origin["reason"] or origin["classification"]
                                counts = report["bridge_unknown_reasons"]
                                counts[reason] = counts.get(reason, 0) + 1
                        continuation = (
                            (outer_rule, slot, source_slot, inner_rule.id, target, inner_guide)
                            if source_slot is not None
                            else None
                        )
                        inner = self._stage(inner_rule, inner_guide, "i", continuation_context=continuation)
                        if slot_trace is not None:
                            slot_trace["inner_stage_calls"] += 1
                            slot_trace["inner_stage_states"] += len(inner)
                        inner = self._compatible_inner(inner, outer_rule, slot)
                        if slot_trace is not None:
                            slot_trace["inner_role_compatible_states"] += len(inner)
                        if self._stop():
                            continue
                        inner = sorted(
                            inner,
                            key=lambda state: (
                                self._rank(state["product"], target),
                                state["product"],
                            ),
                        )[: self.limits.beam_width]
                        if slot_trace is not None:
                            slot_trace["inner_beam_states"] += len(inner)
                        if not inner:
                            continue
                        states = self._stage(
                            outer_rule,
                            guide,
                            "o",
                            replace=(slot, inner),
                            fallback_missing_slots=(
                                self._guide_missing_role_slots.get((outer_rule.id, self.hit, guide, True), ())
                                if _fallback and _quality == "explicit_fallback"
                                else ()
                            ),
                        )
                        if slot_trace is not None:
                            slot_trace["outer_stage_calls"] += 1
                            slot_trace["outer_stage_states"] += len(states)
                        self._retain(states)
                    finally:
                        self._active_phase2_slot = None
                        if plan_trace is not None:
                            rows = self.report["catalog_row_visits"] - before_rows
                            plan_trace["catalog_rows_execution"] += rows
                            slot_trace["catalog_rows_execution"] += rows
            else:
                defer_reason = None
            if defer_reason is not None:
                for plan_index, work in enumerate(work_by_plan):
                    plan_trace = plan_traces[plan_index]
                    if plan_trace is not None and offsets[plan_index] < len(work):
                        plan_trace["work_items_deferred"] = len(work) - offsets[plan_index]
                        plan_trace["defer_reason"] = defer_reason

        self._admit()
        if self._prefer_continuation_second:
            report = self.report["continuation_preference"]
            report["bridge_snapshot_attempts"] = self._continuation_snapshots
            report["bridge_structural_units"] = self._continuation_snapshot_units
            report["bridge_role_contexts"] = self._continuation_role_contexts
        self.report["total_fingerprint_calls"] = (
            self.report["fingerprint_calls"] + self.report["admission_fingerprint_calls"]
        )
        if self.report["total_fingerprint_calls"] != (
            self.report["hit_fingerprint_calls"]
            + self.report["generation_fingerprint_calls"]
            + self.report["admission_rank_fingerprint_calls"]
            + self.report["admission_fingerprint_calls"]
        ):
            raise RuntimeError("Route-first fingerprint purpose accounting mismatch")
        self.report["distinct_exact_products"] = len(
            {item["verification"]["product"] for item in self.routes}
        )
        self.report["depth_counts"] = {
            str(depth): sum(len(item["verification"]["steps"]) == depth for item in self.routes)
            for depth in (1, 2)
        }
        return self.routes, self.report


def generate(hit, catalog, rules, config, *, deadline_monotonic=None, accept=None, hit_anchors=None):
    if not isinstance(config, dict):
        raise ValueError("route_first must be an object")
    if set(config) - set(RouteFirstLimits.__dataclass_fields__):
        raise ValueError("Unknown route_first option")
    limits = RouteFirstLimits(**config)
    if deadline_monotonic is not None and (
        isinstance(deadline_monotonic, bool)
        or not isinstance(deadline_monotonic, (int, float))
        or not math.isfinite(deadline_monotonic)
    ):
        raise ValueError("Invalid route-first deadline")
    return RouteFirstSearch(hit, catalog, rules, limits, deadline_monotonic, accept, hit_anchors).run()
