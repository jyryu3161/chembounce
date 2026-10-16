"""Explicit geometric engineering constraints on one aligned original-hit pose.

Typed feature coverage is not biological essentiality or receptor interaction
evidence. Feature selectors refer to atom indices in canonical original-hit SMILES.
"""

import copy
import hashlib
import json
import math

import numpy as np
from rdkit import Chem
from scipy.optimize import linear_sum_assignment

from .common import canonical, molecule
from .conformers import FEATURES, FAMILIES, _heavy_view

VERSION = "selective-similarity-constraints-v2"
MAX_FEATURES = 256
MAX_LANDMARKS = 32
POLAR = frozenset({"Donor", "Acceptor", "PosIonizable", "NegIonizable"})
MINIMUM_SCORES = {"hit_tanimoto": 0.35, "hit_shape": 0.70, "hit_pharmacophore": 0.50}


def binding(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _positive(value, name, upper=None):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
        or (upper is not None and value > upper)
    ):
        raise ValueError(f"{name} must be finite and positive" + (f" and <= {upper}" if upper else ""))


def indexed_topology(mol):
    return binding(
        {
            "atoms": [(a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge()) for a in mol.GetAtoms()],
            "bonds": sorted(
                (
                    min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                    max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                    str(b.GetBondType()),
                )
                for b in mol.GetBonds()
            ),
        }
    )


def feature_records(mol, cid=None):
    """Use exactly the existing global feature family/lumped-hydrophobe semantics."""
    indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    heavy = _heavy_view(mol)
    raw = [
        f
        for f in FEATURES.GetFeaturesForMol(heavy, confId=-1 if cid is None else cid)
        if f.GetFamily() in FAMILIES
    ]
    lumped = {atom for f in raw if f.GetFamily() == "LumpedHydrophobe" for atom in f.GetAtomIds()}
    records = []
    for feature in raw:
        family, atoms = feature.GetFamily(), feature.GetAtomIds()
        if family == "Hydrophobe" and set(atoms) <= lumped:
            continue
        record = {
            "family": "Hydrophobe" if family == "LumpedHydrophobe" else family,
            "atom_ids": sorted(indices[a] for a in atoms),
        }
        if cid is not None:
            xyz = [float(x) for x in feature.GetPos()]
            if len(xyz) != 3 or not all(math.isfinite(x) for x in xyz):
                raise ValueError("Nonfinite feature coordinates")
            record["position"] = xyz
        records.append(record)
    if len(records) > MAX_FEATURES:
        raise ValueError(f"Similarity feature bound exceeded ({MAX_FEATURES})")
    records.sort(key=lambda r: (r["family"], r["atom_ids"]))
    keys = [(r["family"], tuple(r["atom_ids"])) for r in records]
    if len(keys) != len(set(keys)):
        raise ValueError("Ambiguous duplicate feature identity")
    return records


def pose_binding(mol, cid):
    if not mol.GetConformer(cid).Is3D():
        raise ValueError("Similarity constraints require a 3D pose")
    positions = np.asarray(mol.GetConformer(cid).GetPositions(), dtype=float)
    if positions.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(positions).all():
        raise ValueError("Nonfinite or malformed pose coordinates")
    value = {
        "smiles": canonical(Chem.MolToSmiles(mol, isomericSmiles=True)),
        "atom_order_sha256": indexed_topology(mol),
        "conformer_id": cid,
        "coordinates": positions.tolist(),
    }
    return {**value, "sha256": binding(value)}


class BoundSimilarity:
    """Bind selectors and query identity before any generation or embedding."""

    def __init__(self, hit, minimum_scores=None, *, polar_recall=0.75, core_smiles=None, landmarks=None):
        self.hit = canonical(hit)
        self._hit_mol = molecule(self.hit)
        self._check_state(self._hit_mol)
        self._reference_features = feature_records(self._hit_mol)
        if not self._reference_features:
            raise ValueError("Selective similarity requires nonempty original-hit features")
        minimum = copy.deepcopy(MINIMUM_SCORES if minimum_scores is None else minimum_scores)
        if not isinstance(minimum, dict) or set(minimum) != set(MINIMUM_SCORES):
            raise ValueError("Selective similarity requires the three original-hit minima")
        for name, value in minimum.items():
            _positive(value, name, 1)
        _positive(polar_recall, "polar_recall", 1)
        core = canonical(core_smiles) if core_smiles is not None else None
        if core is not None and not self._hit_mol.HasSubstructMatch(molecule(core), useChirality=True):
            raise ValueError("Required core is not a chiral substructure of the original hit")
        requests = copy.deepcopy([] if landmarks is None else landmarks)
        if not isinstance(requests, list) or len(requests) > MAX_LANDMARKS:
            raise ValueError(f"At most {MAX_LANDMARKS} required landmarks are supported")
        ids, selectors = set(), set()
        for request in requests:
            if not isinstance(request, dict) or set(request) - {
                "id",
                "family",
                "hit_atom_ids",
                "max_distance",
            }:
                raise ValueError("Invalid required landmark fields")
            ident, family, atoms = request.get("id"), request.get("family"), request.get("hit_atom_ids")
            if not isinstance(ident, str) or not ident or ident in ids:
                raise ValueError("Landmark IDs must be unique nonempty strings")
            if not isinstance(family, str) or family not in FAMILIES - {"LumpedHydrophobe"}:
                raise ValueError("Landmark family must name an existing normalized feature type")
            if (
                not isinstance(atoms, list)
                or not atoms
                or any(type(a) is not int or a < 0 for a in atoms)
                or atoms != sorted(set(atoms))
            ):
                raise ValueError("Landmark atom IDs must be sorted unique canonical-hit indices")
            matches = [
                f for f in self._reference_features if f["family"] == family and f["atom_ids"] == atoms
            ]
            key = (family, tuple(atoms))
            if len(matches) != 1 or key in selectors:
                raise ValueError("Landmark must select one distinct existing original-hit feature")
            request.setdefault("max_distance", 1.0)
            _positive(request["max_distance"], "landmark.max_distance")
            ids.add(ident)
            selectors.add(key)
        self._policy = {
            "version": VERSION,
            "hit": self.hit,
            "minimum_scores": minimum,
            "minimum_polar_reference_recall": polar_recall,
            "sigma": 1.0,
            "assignment_policy": "joint-required-landmarks-max-polar-weight-v1",
            "feature_atom_indexing": "Canonical original-hit full atom order; isotope H excluded from feature perception only",
            "core_smiles": core,
            "landmarks": requests,
            "max_features": MAX_FEATURES,
            "max_landmarks": MAX_LANDMARKS,
            "scope": "Engineering geometric constraints; typed coverage does not identify biological essentiality",
        }
        self.sha256 = binding(self._policy)

    @staticmethod
    def _check_state(mol):
        if any(str(s.specified) != "Specified" for s in Chem.FindPotentialStereo(mol)):
            raise ValueError("Selective similarity requires fully specified stereochemistry")

    @property
    def policy(self):
        return copy.deepcopy(self._policy)

    def preflight_candidate(self, smiles):
        mol = molecule(canonical(smiles))
        self._check_state(mol)
        features = feature_records(mol)
        if not features:
            raise ValueError("Selective similarity requires nonempty candidate features")
        core = self._policy["core_smiles"]
        if core is not None and not mol.HasSubstructMatch(molecule(core), useChirality=True):
            raise ValueError("Candidate lost required core/anchor")
        return features

    def assess(self, reference, aligned, ref_cid, candidate_cid):
        if canonical(Chem.MolToSmiles(reference)) != self.hit or indexed_topology(
            reference
        ) != indexed_topology(self._hit_mol):
            raise ValueError("Original-hit pose identity or canonical atom order changed")
        self.preflight_candidate(Chem.MolToSmiles(aligned))
        ref, cand = feature_records(reference, ref_cid), feature_records(aligned, candidate_cid)
        if not ref or not cand:
            raise ValueError("Empty global pharmacophore features")
        ref_pose, cand_pose = pose_binding(reference, ref_cid), pose_binding(aligned, candidate_cid)
        polar_indices = [i for i, f in enumerate(ref) if f["family"] in POLAR]
        requests = self._policy["landmarks"]
        required = {}
        for request in requests:
            index = next(
                i
                for i, f in enumerate(ref)
                if f["family"] == request["family"] and f["atom_ids"] == request["hit_atom_ids"]
            )
            required[index] = request
        # One joint witness: mandatory landmarks cannot steal a feature already
        # counted toward polar recall. Nonpolar requirements add no objective
        # weight; optional polar rows may take a zero-weight dummy column.
        reference_indices = sorted(set(polar_indices) | set(required))
        costs = np.full((len(reference_indices), len(cand) + len(reference_indices)), np.inf)
        edges = {}
        for row, index in enumerate(reference_indices):
            rf, request = ref[index], required.get(index)
            if request is None:
                costs[row, len(cand) :] = 0.0
            for col, cf in enumerate(cand):
                if cf["family"] != rf["family"]:
                    continue
                delta = np.asarray(rf["position"]) - cf["position"]
                distance = float(np.linalg.norm(delta))
                if not math.isfinite(distance):
                    raise ValueError("Nonfinite feature distance")
                weight = math.exp(-float(delta @ delta) / 2)
                if request is not None and distance > request["max_distance"]:
                    continue
                costs[row, col] = -weight if index in polar_indices else 0.0
                edges[row, col] = (distance, weight)
        try:
            rows, cols = linear_sum_assignment(costs)
            feasible = True
        except ValueError:
            rows, cols, feasible = [], [], False
        assignments = []
        for row, col in zip(rows, cols):
            index = reference_indices[row]
            if col >= len(cand):
                continue
            distance, weight = edges[row, col]
            assignments.append(
                {
                    "reference_feature": index,
                    "candidate_feature": int(col),
                    "distance": distance,
                    "weight": weight,
                    "counts_toward_polar_recall": index in polar_indices,
                    "required_landmark": required[index]["id"] if index in required else None,
                }
            )
        recall = (
            (sum(a["weight"] for a in assignments if a["counts_toward_polar_recall"]) / len(polar_indices))
            if polar_indices and feasible
            else None
        )
        polar_passed = feasible and (
            recall is None or recall >= self._policy["minimum_polar_reference_recall"]
        )
        matched = {a["required_landmark"]: a for a in assignments if a["required_landmark"] is not None}
        landmarks = [
            {"id": r["id"], "passed": r["id"] in matched, "assignment": matched.get(r["id"])}
            for r in requests
        ]
        reasons = []
        if not feasible:
            reasons.append("required_landmark_joint_assignment_infeasible")
        elif not polar_passed:
            reasons.append("polar_reference_recall_below_minimum")
        return {
            "policy_sha256": self.sha256,
            "passed": not reasons,
            "reasons": reasons,
            "reference_pose": ref_pose,
            "candidate_pose_sha256": cand_pose["sha256"],
            "reference_features": ref,
            "candidate_features": cand,
            "joint_assignment": {
                "scope": "Polar recall and required landmarks; global pharmacophore overlap has its own full-feature assignment on this same pose",
                "feasible": feasible,
                "assignments": assignments,
            },
            "polar": {
                "status": "unassessed"
                if not feasible
                else "not_applicable"
                if not polar_indices
                else "passed"
                if polar_passed
                else "failed",
                "reference_count": len(polar_indices),
                "candidate_count": sum(f["family"] in POLAR for f in cand),
                "recall": recall,
            },
            "landmarks": landmarks,
        }

    def validate_witness(self, witness, aligned, candidate_cid):
        """Recompute the joint geometry from saved reference and actual output pose.

        This is a consistency check, not independent provenance for an externally
        supplied pose or evidence that an ensemble search was actually executed.
        """
        try:
            saved = witness["reference_pose"]
            cid = saved["conformer_id"]
            if type(cid) is not int or cid < 0:
                raise ValueError("Invalid reference conformer ID")
            reference = molecule(self.hit)
            xyz = np.asarray(saved["coordinates"], dtype=float)
            if xyz.shape != (reference.GetNumAtoms(), 3) or not np.isfinite(xyz).all():
                raise ValueError("Invalid saved reference coordinates")
            conformer = Chem.Conformer(reference.GetNumAtoms())
            conformer.SetId(cid)
            conformer.Set3D(True)
            for atom, position in enumerate(xyz):
                conformer.SetAtomPosition(atom, position.tolist())
            reference.AddConformer(conformer, assignId=False)
            actual = self.assess(reference, aligned, cid, candidate_cid)
            if not actual["passed"] or binding(actual) != binding(witness):
                raise ValueError("Selective constraint witness mismatch or failed joint geometry")
            return reference
        except (KeyError, TypeError, AttributeError, RuntimeError) as exc:
            raise ValueError(f"Incomplete selective constraint witness: {exc}") from exc
