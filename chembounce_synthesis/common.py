"""Strict, stereochemistry-aware identities and portable provenance."""

import hashlib
import json
import re
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)


def require_absolute_stereo(mol):
    """Reject stereo information ordinary SMILES cannot faithfully retain.

    OR/AND groups and atropisomeric bonds are unsupported. Unconsumed file
    wedging is rejected as unsupported/ambiguous, not classified as an atropisomer.
    """
    tetra = {Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW}
    for group in mol.GetStereoGroups():
        if group.GetGroupType() != Chem.StereoGroupType.STEREO_ABSOLUTE:
            raise ValueError(
                "Unsupported nonabsolute enhanced stereochemistry: "
                + str(group.GetGroupType())
                + "; supply an explicitly defined absolute state, not an OR/AND group"
            )
        if (
            not group.GetAtoms()
            or any(a.GetChiralTag() not in tetra for a in group.GetAtoms())
            or (hasattr(group, "GetBonds") and group.GetBonds())
        ):
            raise ValueError("Unsupported ABS group: only specified tetrahedral atom groups are supported")
    for bond in mol.GetBonds():
        if str(bond.GetStereo()) in {"STEREOATROPCW", "STEREOATROPCCW"}:
            raise ValueError("Unsupported atropisomeric bond stereochemistry; no axial state is selected")
        # Older RDKit SDF readers can retain the wedge property while losing the
        # axial interpretation. A wedge not accounted for by a specified
        # tetrahedral origin is not safe to flatten, regardless of runtime age.
        wedged = (
            str(bond.GetBondDir()) in {"BEGINWEDGE", "BEGINDASH"}
            or (bond.HasProp("_MolFileBondStereo") and bond.GetIntProp("_MolFileBondStereo") in {1, 6})
            or (bond.HasProp("_MolFileBondCfg") and bond.GetIntProp("_MolFileBondCfg") in {1, 3})
        )
        if wedged and bond.GetBeginAtom().GetChiralTag() not in tetra:
            raise ValueError(
                "Unsupported or unresolved SDF/molecule bond wedging outside specified tetrahedral stereo"
            )


def parse_absolute_smiles(smiles):
    """Exact-state text: ordinary SMILES or one supported tetrahedral ABS block.

    Reject other CX annotations *before* parsing: older RDKit can silently discard
    wU/wD axial stereo. Names/extra blocks are not part of this exact-state format.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Expected a nonempty SMILES string")
    text = smiles.strip()
    absolute_indices = None
    if "|" in text:
        if re.search(r"(?:[|,])(?:o|&)\d+:", text):
            raise ValueError("Unsupported nonabsolute enhanced stereochemistry: OR/AND CX group")
        match = re.fullmatch(r"([^\s|]+)\s+\|a:(\d+(?:,\d+)*)\|", text)
        if match is None:
            raise ValueError(
                "Unsupported CXSMILES annotation: only one tetrahedral ABS-only block is supported; names/extra fields are not accepted"
            )
        absolute_indices = [int(value) for value in match.group(2).split(",")]
    elif any(c.isspace() for c in text):
        raise ValueError("Exact-state SMILES does not accept trailing names or embedded whitespace")
    params = Chem.SmilesParserParams()
    params.parseName = False
    params.allowCXSMILES = absolute_indices is not None
    # CX atom indices refer to the input graph, before ordinary explicit H is
    # removed. Validate that graph first; default parsing may renumber centers.
    if absolute_indices is not None:
        params.removeHs = False
    mol = Chem.MolFromSmiles(text, params)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    require_absolute_stereo(mol)
    if absolute_indices is not None:
        retained = {a.GetIdx() for group in mol.GetStereoGroups() for a in group.GetAtoms()}
        if len(set(absolute_indices)) != len(absolute_indices) or set(absolute_indices) != retained:
            raise ValueError("ABS atom references are invalid, duplicated, or not retained by the parser")
        # Reparse validated text with the established default H handling, rather
        # than approximating it with a separately configured RemoveHs operation.
        params.removeHs = True
        mol = Chem.MolFromSmiles(text, params)
        if mol is None or mol.GetNumAtoms() == 0:
            raise ValueError(f"Invalid SMILES after ordinary hydrogen normalization: {smiles!r}")
        require_absolute_stereo(mol)
        if sum(len(group.GetAtoms()) for group in mol.GetStereoGroups()) != len(absolute_indices):
            raise ValueError("ABS centers were not retained during ordinary hydrogen normalization")
    return mol


def absolute_smiles(mol, **kwargs):
    """Serialize only after checking stereo information ordinary SMILES cannot keep."""
    require_absolute_stereo(mol)
    return Chem.MolToSmiles(mol, **kwargs)


def molecule(smiles):
    mol = parse_absolute_smiles(smiles)
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError(f"A connected molecule without dummy atoms is required: {smiles!r}")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol


def canonical(smiles):
    return Chem.MolToSmiles(molecule(smiles), isomericSmiles=True)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temp.replace(path)


def positive(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value
