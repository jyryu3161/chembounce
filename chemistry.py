"""Attachment-aware scaffold replacement with a deterministic bounded beam.

Only the first deterministic substructure match is replaced. All bonds crossing
that match are retained, and RDKit molzip preserves stereo at retained atoms.
The beam is a heuristic, not an exhaustive enumeration of substitution isomers.
"""

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator

GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2)


def count_fingerprint(mol):
    return GENERATOR.GetSparseCountFingerprint(mol)


def cut_scaffold(target, pattern):
    """Return the outside components, labelled at each severed bond."""
    matched = target.GetSubstructMatch(pattern, useChirality=True)
    if not matched:
        raise ValueError("Fragment is not a substructure of the input molecule")
    match = set(matched)
    work = Chem.Mol(target)
    for atom in work.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetBoolProp("_cb_outside", atom.GetIdx() not in match)
    bonds = [
        b.GetIdx() for b in work.GetBonds() if (b.GetBeginAtomIdx() in match) != (b.GetEndAtomIdx() in match)
    ]
    if not bonds:
        return []
    if any(work.GetBondWithIdx(index).GetIsAromatic() for index in bonds):
        raise ValueError("Fragment boundary crosses aromatic bonds; select the entire aromatic ring system")
    cut = Chem.FragmentOnBonds(work, bonds, dummyLabels=[(i, i) for i in range(1, len(bonds) + 1)])
    outside = []
    for part in Chem.GetMolFrags(cut, asMols=True, sanitizeFrags=False):
        if any(a.HasProp("_cb_outside") and a.GetBoolProp("_cb_outside") for a in part.GetAtoms()):
            Chem.SanitizeMol(part)
            outside.append(part)
    return sorted(outside, key=Chem.MolToSmiles)


def _attach(state, site_idx, label):
    edit = Chem.RWMol(state)
    endpoint = next(a for a in edit.GetAtoms() if a.GetAtomicNum() == 0 and a.GetIsotope() == label)
    bond_type = endpoint.GetBonds()[0].GetBondType()
    endpoint.SetAtomMapNum(label)
    endpoint.SetIsotope(0)
    site = edit.GetAtomWithIdx(site_idx)
    # Explicit hydrogen on aromatic [nH] or bracket atoms is replaced by the bond.
    order = int(bond_type)
    if site.GetNumExplicitHs():
        site.SetNumExplicitHs(max(0, site.GetNumExplicitHs() - order))
    dummy = Chem.Atom(0)
    dummy.SetAtomMapNum(label)
    idx = edit.AddAtom(dummy)
    edit.AddBond(site_idx, idx, bond_type)
    product = Chem.molzip(edit.GetMol())
    Chem.SanitizeMol(product)
    Chem.AssignStereochemistry(product, cleanIt=True, force=True)
    return product


def prepare_replacement(target, pattern):
    """Cache invariant cutting and query fingerprint work once per fragment."""
    return tuple(cut_scaffold(target, pattern)), count_fingerprint(target)


def _insufficient_capacity(replacement, outside):
    """Conservatively reject impossible neutral C/N/O replacements before search.

    Count bond-order demand against an upper bound after replacing all removable
    hydrogens. N receives a deliberately generous valence bound of five because
    sanitization can convert hypervalent nitro-like representations to charged
    forms. Other elements, charged scaffolds and unusual bonds use normal assembly;
    no capacity assumption is made for them. Restrict retained atoms to C/N/O/H
    too, avoiding metal/dative-bond cleanup that could change the bond-order demand.
    """
    maximum = {6: 4, 7: 5, 8: 2}
    capacity = 0
    for atom in replacement.GetAtoms():
        if atom.GetAtomicNum() not in maximum or atom.GetFormalCharge():
            return False
        capacity += maximum[atom.GetAtomicNum()] - atom.GetTotalValence() + atom.GetTotalNumHs()
    ordinary = {Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE}
    if any(b.GetBondType() not in ordinary | {Chem.BondType.AROMATIC} for b in replacement.GetBonds()):
        return False
    demand = 0
    for part in outside:
        for atom in part.GetAtoms():
            if atom.GetAtomicNum() not in (0, 1, 6, 7, 8):
                return False
            if atom.GetAtomicNum() == 0:
                bond = atom.GetBonds()[0]
                if bond.GetBondType() not in ordinary:
                    return False
                demand += bond.GetBondTypeAsDouble()
    return demand > capacity


def replace(target, pattern, replacement, max_n=10, beam_width=100, *, prepared=None):
    if max_n < 1 or beam_width < 1:
        raise ValueError("Candidate and beam limits must be positive")
    if replacement is None or not replacement.GetNumAtoms():
        raise ValueError("Replacement scaffold must be a nonempty molecule")
    if len(Chem.GetMolFrags(replacement)) != 1 or any(a.GetAtomicNum() == 0 for a in replacement.GetAtoms()):
        raise ValueError("Replacement scaffold must be connected and contain no dummy atoms")
    outside, query_fp = prepared if prepared is not None else prepare_replacement(target, pattern)
    if _insufficient_capacity(replacement, outside):
        return []
    seed = Chem.Mol(replacement)
    for atom in seed.GetAtoms():
        atom.SetAtomMapNum(0)
        atom.SetBoolProp("_cb_site", True)
    for part in outside:
        seed = Chem.CombineMols(seed, part)
    labels = sorted(a.GetIsotope() for a in seed.GetAtoms() if a.GetAtomicNum() == 0)
    states = [seed]
    for label in labels:
        candidates = {}
        for state in states:
            sites = [a.GetIdx() for a in state.GetAtoms() if a.HasProp("_cb_site")]
            for site in sites:
                try:
                    with rdBase.BlockLogs():
                        product = _attach(state, site, label)
                    smiles = Chem.MolToSmiles(product, isomericSmiles=True)
                    if smiles not in candidates:
                        score = DataStructs.TanimotoSimilarity(query_fp, count_fingerprint(product))
                        candidates[smiles] = (score, product)
                except (ValueError, RuntimeError):
                    # Invalid valence at a proposed attachment is an expected rejection.
                    continue
        ranked = sorted(candidates.items(), key=lambda x: (-x[1][0], x[0]))[:beam_width]
        states = [value[1] for _, value in ranked]
        if not states:
            break
    final = {}
    for mol in states:
        if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
            continue
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        final[smiles] = DataStructs.TanimotoSimilarity(query_fp, count_fingerprint(mol))
    return [(score, smiles) for smiles, score in sorted(final.items(), key=lambda x: (-x[1], x[0]))[:max_n]]
