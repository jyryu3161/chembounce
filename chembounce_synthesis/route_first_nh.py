"""Conservative mapped-N/H retrieval hint; strict route audit remains authoritative."""

from rdkit import Chem


SENTINEL = 65


def _safe_n(atom):
    if (
        atom.GetAtomicNum() != 7
        or atom.GetIsAromatic()
        or atom.GetFormalCharge()
        or atom.GetIsotope()
        or atom.GetNumRadicalElectrons()
        or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        or any(bond.GetBondType() != Chem.BondType.SINGLE for bond in atom.GetBonds())
    ):
        return False
    if any(
        neighbor.GetAtomicNum() == 1 for neighbor in (bond.GetOtherAtom(atom) for bond in atom.GetBonds())
    ):
        # An explicit H, especially an isotope-labelled H, needs its own audit.
        return False
    return atom.GetTotalValence() == 3


def prove_parent_n(reaction, parent):
    """Return a uniquely mapped product-N site, or an explicit unknown reason."""
    reactants, products = reaction.GetReactants(), reaction.GetProducts()
    if len(products) != 1:
        return None, "product_template_count"
    mapped = {}
    for slot, template in enumerate(reactants):
        for atom in template.GetAtoms():
            number = atom.GetAtomMapNum()
            if number:
                mapped.setdefault(number, []).append((slot, atom))
    product = products[0]
    product_mapped = {}
    for atom in product.GetAtoms():
        number = atom.GetAtomMapNum()
        if number:
            product_mapped.setdefault(number, []).append(atom)
    candidates = []
    for number, atoms in product_mapped.items():
        if len(atoms) != 1 or len(mapped.get(number, ())) != 1:
            continue
        atom = atoms[0]
        source_slot, source_atom = mapped[number][0]
        if (
            atom.GetAtomicNum() != 7
            or source_atom.GetAtomicNum() != 7
            or atom.GetFormalCharge()
            or source_atom.GetFormalCharge()
            or atom.GetNumExplicitHs()
            or source_atom.GetNumExplicitHs()
            or atom.GetIsotope()
            or source_atom.GetIsotope()
            or atom.GetIsAromatic()
            or source_atom.GetIsAromatic()
            or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            or source_atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            or atom.GetNumRadicalElectrons()
            or source_atom.GetNumRadicalElectrons()
            or atom.DescribeQuery().strip() != "AtomType 7 = val"
        ):
            continue
        before = {
            (bond.GetOtherAtom(source_atom).GetAtomMapNum(), bond.GetBondType())
            for bond in source_atom.GetBonds()
        }
        after = {(bond.GetOtherAtom(atom).GetAtomMapNum(), bond.GetBondType()) for bond in atom.GetBonds()}
        if any(other == 0 for other, _ in before | after):
            continue
        added = after - before
        if before - after or len(added) != 1:
            continue
        other, order = next(iter(added))
        if (
            order != Chem.BondType.SINGLE
            or len(product_mapped.get(other, ())) != 1
            or len(mapped.get(other, ())) != 1
            or mapped[other][0][0] == source_slot
        ):
            continue
        new_bond = atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), product_mapped[other][0].GetIdx())
        if new_bond.GetStereo() != Chem.BondStereo.STEREONONE or new_bond.GetBondDir() != Chem.BondDir.NONE:
            continue
        candidates.append((number, source_slot, atom.GetIdx(), product_mapped[other][0].GetIdx()))
    if len(candidates) != 1:
        return None, "reactive_n_template_ambiguous_or_unsupported"
    number, source_slot, template_index, other_template_index = candidates[0]
    matches = parent.GetSubstructMatches(product, useChirality=True, uniquify=False, maxMatches=SENTINEL)
    if not matches or len(matches) >= SENTINEL:
        return None, "product_template_match_absent_or_truncated"
    positions = {(match[template_index], match[other_template_index]) for match in matches}
    if len(positions) != 1:
        return None, "parent_reactive_n_or_new_bond_ambiguous"
    n_index, other_index = next(iter(positions))
    atom = parent.GetAtomWithIdx(n_index)
    if not _safe_n(atom) or atom.GetDegree() != 2 or atom.GetTotalNumHs() != 1:
        return None, "parent_n_state_unsupported"
    if parent.GetBondBetweenAtoms(n_index, other_index).GetBondType() != Chem.BondType.SINGLE:
        return None, "parent_new_bond_not_single"
    return {
        "map_number": number,
        "source_slot": source_slot,
        "parent_n": n_index,
        "parent_added_neighbor": other_index,
        "parent_h": 1,
        "required_source_h": 2,
    }, None


def prove_raw_n(outcomes, parent, site, max_outcomes=4):
    """Require all actual raw alternatives to map the same site and one lost bond."""
    if not outcomes or len(outcomes) >= max_outcomes:
        return False, "raw_absent_or_cap_saturated"
    parent_n = parent.GetAtomWithIdx(site["parent_n"])
    expected = sorted(
        (bond.GetOtherAtomIdx(site["parent_n"]), str(bond.GetBondType()))
        for bond in parent_n.GetBonds()
        if bond.GetOtherAtomIdx(site["parent_n"]) != site["parent_added_neighbor"]
    )
    for outcome in outcomes:
        if site["source_slot"] >= len(outcome):
            return False, "raw_source_slot_missing"
        found = []
        for slot, precursor in enumerate(outcome):
            for atom in precursor.GetAtoms():
                if (
                    atom.HasProp("react_idx")
                    and atom.HasProp("react_atom_idx")
                    and atom.HasProp("old_mapno")
                    and atom.GetIntProp("react_idx") == 0
                    and atom.GetIntProp("react_atom_idx") == site["parent_n"]
                    and atom.GetIntProp("old_mapno") == site["map_number"]
                ):
                    found.append((slot, atom))
        if len(found) != 1 or found[0][0] != site["source_slot"]:
            return False, "raw_n_origin_missing_ambiguous_or_other_slot"
        atom = found[0][1]
        if (
            atom.GetAtomicNum() != 7
            or atom.GetFormalCharge() != parent_n.GetFormalCharge()
            or atom.GetIsotope() != parent_n.GetIsotope()
            or atom.GetIsAromatic() != parent_n.GetIsAromatic()
            or atom.GetNumRadicalElectrons() != parent_n.GetNumRadicalElectrons()
            or atom.GetChiralTag() != parent_n.GetChiralTag()
        ):
            return False, "raw_n_identity_changed"
        actual = []
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if (
                not neighbor.HasProp("react_idx")
                or not neighbor.HasProp("react_atom_idx")
                or neighbor.GetIntProp("react_idx") != 0
            ):
                return False, "raw_n_originless_or_unmapped_neighbor"
            actual.append((neighbor.GetIntProp("react_atom_idx"), str(bond.GetBondType())))
        if sorted(actual) != expected:
            return False, "raw_n_local_bonds_changed_beyond_one_added_bond"
    return True, None


def candidate_role_n_site(mol, role, map_number, required_h=2):
    """Classify one original-role source; unknown is never skipped in ranking."""
    if mol.GetNumAtoms() > 256 or mol.GetNumBonds() > 512:
        return None, "candidate_size_limit", None
    matches = mol.GetSubstructMatches(role, useChirality=True, uniquify=False, maxMatches=SENTINEL)
    if len(matches) >= SENTINEL:
        return None, "role_match_truncated", None
    if not matches:
        return False, "role_absent", None
    mapped = [atom.GetIdx() for atom in role.GetAtoms() if atom.GetAtomMapNum() == map_number]
    if len(mapped) != 1:
        return None, "mapped_role_n_missing_or_ambiguous", None
    # The mapped N site must be unique over every bounded role embedding.
    positions = {match[mapped[0]] for match in matches}
    if len(positions) != 1:
        return None, "reactive_n_role_ambiguous", None
    index = next(iter(positions))
    atom = mol.GetAtomWithIdx(index)
    if not _safe_n(atom):
        return None, "candidate_n_state_unsupported", None
    return atom.GetTotalNumHs() == required_h, None, index


def candidate_role_n(mol, role, map_number, required_h=2):
    eligible, reason, _ = candidate_role_n_site(mol, role, map_number, required_h)
    return eligible, reason


def same_n_origin_and_bonds(raw, clean, parent_n):
    """Prove the accepted hint retains one parent-derived N and its bond origins.

    Hydrogen state in a reverse hint is deliberately not used to prove source
    eligibility; the concrete catalog role and H count are checked separately.
    """

    def origin(atom):
        if not atom.HasProp("react_idx") or not atom.HasProp("react_atom_idx"):
            return None
        return atom.GetIntProp("react_idx"), atom.GetIntProp("react_atom_idx")

    def site(mol):
        found = [atom for atom in mol.GetAtoms() if origin(atom) == (0, parent_n)]
        return found[0] if len(found) == 1 else None

    before, after = site(raw), site(clean)
    if before is None or after is None:
        return False
    identity = lambda atom: (
        atom.GetAtomicNum(),
        atom.GetIsotope(),
        atom.GetFormalCharge(),
        atom.GetIsAromatic(),
        atom.GetNumRadicalElectrons(),
        int(atom.GetChiralTag()),
    )
    if identity(before) != identity(after) or before.GetAtomicNum() != 7:
        return False

    def bonds(atom):
        neighbors = []
        for bond in atom.GetBonds():
            key = origin(bond.GetOtherAtom(atom))
            if key is None:
                return None
            neighbors.append((key, str(bond.GetBondType()), bond.GetIsAromatic()))
        return sorted(neighbors)

    old_bonds, new_bonds = bonds(before), bonds(after)
    return old_bonds is not None and old_bonds == new_bonds
