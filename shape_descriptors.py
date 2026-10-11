"""ElectroShape descriptor using RDKit coordinates and Gasteiger charges.

Method: Armstrong et al. (2010), doi:10.1007/s10822-010-9374-0.
Five reference points in (x,y,z,25*charge), each with mean distance,
standard deviation, and cube root of the third central moment.
The implementation is numerically compared with ODDT in the regression suite.
"""

import numpy as np
from rdkit.Chem import AllChem, rdMolDescriptors


def electroshape(mol):
    if not mol.GetNumConformers():
        raise ValueError("ElectroShape needs a 3D conformer")
    AllChem.ComputeGasteigerCharges(mol)
    charges = np.array([float(atom.GetProp("_GasteigerCharge")) for atom in mol.GetAtoms()])
    if not np.isfinite(charges).all():
        raise ValueError("Non-finite Gasteiger partial charges")
    coordinates = np.asarray(mol.GetConformer().GetPositions())
    points = np.column_stack((coordinates, 25 * charges))
    center = points.mean(axis=0)
    farthest = points[np.linalg.norm(points - center, axis=1).argmax()]
    opposite = points[np.linalg.norm(points - farthest, axis=1).argmax()]
    a, b = farthest - center, opposite - center
    normal = np.cross(a[:3], b[:3])
    magnitude = np.linalg.norm(normal)
    if magnitude < 1e-12:
        raise ValueError("ElectroShape is undefined for degenerate reference points")
    offset = center[:3] + normal * np.linalg.norm(a) / (2 * magnitude)
    reference = [
        center,
        farthest,
        opposite,
        np.append(offset, 25 * charges.max()),
        np.append(offset, 25 * charges.min()),
    ]
    features = []
    for point in reference:
        distances = np.linalg.norm(points - point, axis=1)
        mean = distances.mean()
        features.extend((mean, distances.std(), np.cbrt(np.mean((distances - mean) ** 3))))
    return np.asarray(features)


def usrcat(mol):
    return np.asarray(rdMolDescriptors.GetUSRCAT(mol))


def similarity(left, right):
    if left.shape != right.shape or left.size not in (15, 60):
        raise ValueError("Mismatched shape descriptor dimensions")
    return float(1 / (1 + np.abs(left - right).mean()))
