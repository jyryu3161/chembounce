"""Separate exact scoring structures from chemical-family novelty identities.

Family identities deliberately ignore stereochemistry, isotopes, reversible
protonation and RDKit-enumerated tautomers. They are *not* representative states
at any pH and must never replace a structure supplied to scoring or reaction
replay. The supported input remains a connected molecule, not a salt mixture.
"""

from dataclasses import dataclass
from functools import lru_cache

from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

from .common import canonical, molecule

IDENTITY_POLICY = "protonation-tautomer-family-v1"
TAUTOMER_LIMIT = 1000
TRANSFORM_LIMIT = 1000
IDENTITY_CACHE_SIZE = 32768


@dataclass(frozen=True)
class ChemicalIdentity:
    """Exact state and increasingly coarse, explicitly labelled identities."""

    exact_smiles: str
    family_smiles: str
    scaffold_family_smiles: str
    ring_skeleton_smiles: str
    has_scaffold: bool
    policy_version: str = IDENTITY_POLICY
    rdkit_version: str = rdBase.rdkitVersion


def _protonation_parent(mol):
    """Remove reversible proton charges, retaining non-proton charge changes.

    Older supported RDKit Uncharger versions can add/remove hydride as well as
    protons. Accept a local charge change only when its hydrogen-count change
    is equal: proton addition changes both by +1, removal changes both by -1.
    This also preserves permanent ammonium and carbocation/borohydride charges
    on the minimum supported RDKit without version-specific constructor flags.
    ``force`` permits a protonation family to have a permanent net charge; it
    does not authorize removal of permanent charges.
    """
    result = rdMolStandardize.Uncharger(canonicalOrder=True, force=True).uncharge(mol)
    for original, atom in zip(mol.GetAtoms(), result.GetAtoms()):
        charge_delta = atom.GetFormalCharge() - original.GetFormalCharge()
        hydrogen_delta = atom.GetTotalNumHs() - original.GetTotalNumHs()
        if charge_delta != hydrogen_delta:
            atom.SetFormalCharge(original.GetFormalCharge())
            atom.SetNumExplicitHs(original.GetNumExplicitHs())
            atom.SetNoImplicit(original.GetNoImplicit())
            atom.SetNumRadicalElectrons(original.GetNumRadicalElectrons())
    Chem.SanitizeMol(result)
    return result


def _canonical_tautomer(mol):
    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(TAUTOMER_LIMIT)
    enumerator.SetMaxTransforms(TRANSFORM_LIMIT)
    result = enumerator.Enumerate(mol)
    if result.status != rdMolStandardize.TautomerEnumeratorStatus.Completed:
        # A truncated enumeration can depend on the input tautomer and must not
        # silently assign a new-family identity to an equivalent chemical state.
        raise ValueError(f"Chemical-family tautomer enumeration incomplete: {result.status}")
    return enumerator.PickCanonical(result)


def _family_parent(mol):
    parent = Chem.Mol(mol)
    Chem.RemoveStereochemistry(parent)
    for atom in parent.GetAtoms():
        atom.SetIsotope(0)
    parent = Chem.RemoveHs(parent)
    return _canonical_tautomer(_protonation_parent(parent))


@dataclass(frozen=True)
class _IdentityFailure:
    reason: str


@lru_cache(maxsize=IDENTITY_CACHE_SIZE)
def _canonical_identity_input(smiles):
    return canonical(smiles)


@lru_cache(maxsize=IDENTITY_CACHE_SIZE)
def _scaffold_identity(scaffold_smiles, normalization_signature):
    # Canonicalize independently of decorations. Many dataset analogues share
    # exactly this expensive normalization, even when full molecules differ.
    scaffold = _family_parent(molecule(scaffold_smiles))
    # Exocyclic double-bonded atoms belong to the atom-aware Murcko scaffold but
    # not the generic ring/linker graph after making every bond single.
    skeleton = MurckoScaffold.GetScaffoldForMol(MurckoScaffold.MakeScaffoldGeneric(scaffold))
    return Chem.MolToSmiles(scaffold, isomericSmiles=False), Chem.MolToSmiles(skeleton, isomericSmiles=False)


@lru_cache(maxsize=IDENTITY_CACHE_SIZE)
def _identity(exact_smiles, normalization_signature):
    try:
        parent = _family_parent(molecule(exact_smiles))
        scaffold = MurckoScaffold.GetScaffoldForMol(parent)
        has_scaffold = bool(scaffold.GetNumAtoms())
        scaffold_smiles, skeleton_smiles = (
            _scaffold_identity(Chem.MolToSmiles(scaffold, isomericSmiles=False), normalization_signature)
            if has_scaffold
            else ("", "")
        )
        return ChemicalIdentity(
            exact_smiles=exact_smiles,
            family_smiles=Chem.MolToSmiles(parent, isomericSmiles=False),
            scaffold_family_smiles=scaffold_smiles,
            ring_skeleton_smiles=skeleton_smiles,
            has_scaffold=has_scaffold,
        )
    except ValueError as exc:
        # Failed normalization is deterministic too. Store the reason, never a
        # partial identity, so repeated dataset passes still fail closed without
        # repeating an exhausted tautomer search. Limits/version are cache keys.
        return _IdentityFailure(str(exc))


def chemical_identity(smiles):
    """Return a deterministic family key while retaining the exact input state.

    The policy and RDKit version accompany each identity because the underlying
    tautomer rule set is versioned. An incomplete normalization fails closed.
    """
    if not isinstance(smiles, str):
        canonical(smiles)  # Preserve the public validation error for invalid types.
    signature = (IDENTITY_POLICY, rdBase.rdkitVersion, TAUTOMER_LIMIT, TRANSFORM_LIMIT)
    result = _identity(_canonical_identity_input(smiles), signature)
    if isinstance(result, _IdentityFailure):
        raise ValueError(result.reason)
    return result


def identity_cache_info():
    """Expose bounded cache statistics for reproducible performance diagnostics."""
    return {
        "exact_inputs": _canonical_identity_input.cache_info()._asdict(),
        "identities_including_failures": _identity.cache_info()._asdict(),
        "scaffolds": _scaffold_identity.cache_info()._asdict(),
    }


def clear_identity_caches():
    """Release all cached identities and normalization errors without policy changes."""
    _canonical_identity_input.cache_clear()
    _identity.cache_clear()
    _scaffold_identity.cache_clear()


def family_key(smiles):
    """Identity for counting distinct chemical families, not scoring states."""
    return chemical_identity(smiles).family_smiles


def scaffold_family_key(smiles):
    """Atom/bond-aware, state-normalized Murcko identity; empty for acyclic inputs."""
    return chemical_identity(smiles).scaffold_family_smiles


def classify_change(reference, candidate):
    """Return ``unchanged``, ``state_only``, ``decoration`` or ``scaffold``.

    ``state_only`` includes stereo/isotope changes for scaffold-novelty purposes,
    not an assertion that the two structures have identical biological effects.
    Heteroatom substitutions are scaffold changes even with the same topology.
    Two acyclic molecules have no Murcko scaffold and are decoration changes.
    """
    original, proposed = chemical_identity(reference), chemical_identity(candidate)
    if original.exact_smiles == proposed.exact_smiles:
        return "unchanged"
    if original.family_smiles == proposed.family_smiles:
        return "state_only"
    if original.scaffold_family_smiles == proposed.scaffold_family_smiles:
        return "decoration"
    return "scaffold"
