"""Opt-in bounded rigid proposals for explicitly requested landmarks only.

Correspondences propose poses; the unchanged joint evaluator decides feasibility.
No embedding, reflection, conformer search, or new alignment optimization occurs.
"""

import copy
import hashlib
import math
from pathlib import Path
import time

import numpy as np
from rdkit import Chem

from . import conformers
from .similarity_constraints import binding, indexed_topology, pose_binding

VERSION = "landmark-anchored-rigid-v1"
LIMITS = {
    "alternatives_per_pair": 2,
    "alternatives_per_call": 8,
    "scoremol_per_call": 16,
    "mapping_edge_visits_per_pair": 128,
    "complete_maps_per_pair": 4,
    "helper_edge_visits_per_pair": 4096,
    "cooperative_seconds_per_call": 1.0,
}
PARITY_TOLERANCE = 1e-5
SOURCE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def require_capability():
    api = conformers.rdShapeAlign
    if (
        api is None
        or not callable(getattr(api, "ScoreMol", None))
        or not callable(getattr(api, "ShapeInputOptions", None))
    ):
        raise ValueError(
            "Landmark pose proposals require fixed-pose rdShapeAlign.ScoreMol and ShapeInputOptions"
        )


def declaration():
    require_capability()
    return {
        "version": VERSION,
        "limits": dict(LIMITS),
        "parity_tolerance": PARITY_TOLERANCE,
        "source_sha256": SOURCE_SHA256,
        "scope": "Opt-in rigid proposals only; correspondences are not a final assignment; no calibration or efficacy claim",
    }


def proper_transform(source, target):
    """Pivot-preserving proper fit; deterministic shortest rotation for rank one."""
    p, q = np.asarray(source, dtype=float), np.asarray(target, dtype=float)
    if (
        p.shape != q.shape
        or p.ndim != 2
        or p.shape[1] != 3
        or not len(p)
        or not np.isfinite(p).all()
        or not np.isfinite(q).all()
    ):
        raise ValueError("Invalid proposal anchor positions")
    a, b = p[1:] - p[0], q[1:] - q[0]
    valid = [i for i in range(len(a)) if np.linalg.norm(a[i]) > 1e-10 and np.linalg.norm(b[i]) > 1e-10]
    rotation = np.eye(3)
    if (
        len(valid) >= 2
        and np.linalg.matrix_rank(a[valid], tol=1e-10) >= 2
        and np.linalg.matrix_rank(b[valid], tol=1e-10) >= 2
    ):
        u, _, vt = np.linalg.svd(a[valid].T @ b[valid])
        fix = np.eye(3)
        fix[-1, -1] = np.linalg.det(vt.T @ u.T)
        rotation = vt.T @ fix @ u.T
    elif valid:
        x, y = a[valid[0]], b[valid[0]]
        x, y = x / np.linalg.norm(x), y / np.linalg.norm(y)
        cross, cosine = np.cross(x, y), float(np.clip(x @ y, -1, 1))
        if np.linalg.norm(cross) <= 1e-10:
            if cosine < 0:
                basis = np.eye(3)[int(np.argmin(np.abs(x)))]
                axis = np.cross(x, basis)
                axis /= np.linalg.norm(axis)
                rotation = 2 * np.outer(axis, axis) - np.eye(3)
        else:
            sine = float(np.linalg.norm(cross))
            axis = cross / sine
            k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
            rotation = np.eye(3) + sine * k + (1 - cosine) * (k @ k)
    if (
        not np.isfinite(rotation).all()
        or not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-8, rtol=0)
        or not math.isclose(float(np.linalg.det(rotation)), 1, abs_tol=1e-8)
    ):
        raise ValueError("Proposal is not a proper rigid rotation")
    return rotation, q[0] - rotation @ p[0]


class ProposalSearch:
    def __init__(self, policy, *, deadline=None, setup_started=None, clock=time.monotonic):
        if policy != declaration():
            raise ValueError("Landmark proposal policy/source mismatch")
        self.clock, self.deadline = clock, deadline
        self.spent = 0.0 if setup_started is None else max(0.0, clock() - setup_started)
        self.active = None
        self.policy = copy.deepcopy(policy)
        self.audit = {
            "policy": self.policy,
            "alternatives_attempted": 0,
            "alternatives_returned": 0,
            "scoremol_attempted": 0,
            "scoremol_returned": 0,
            "pairs": [],
            "truncated": False,
            "reasons": [],
            "elapsed_seconds": self.spent,
            "scope": "Added work only; cooperative allowance, native calls may overrun; late proposals discarded",
        }

    def charge(self, seconds):
        self.spent += max(0.0, seconds)
        self.audit["elapsed_seconds"] = self.spent
        self.audit["cooperative_overrun"] = self.spent > LIMITS["cooperative_seconds_per_call"]

    def elapsed(self):
        return self.spent + (0 if self.active is None else max(0, self.clock() - self.active))

    def available(self):
        reason = None
        if self.deadline is not None and self.clock() >= self.deadline:
            reason = "outer_deadline_exhausted"
        elif self.elapsed() >= LIMITS["cooperative_seconds_per_call"]:
            reason = "proposal_allowance_exhausted"
        if reason:
            self.audit["truncated"] = True
            if reason not in self.audit["reasons"]:
                self.audit["reasons"].append(reason)
        return reason is None

    def _score(self, reference, candidate, rid, cid):
        if not self.available():
            return None
        if self.audit["scoremol_attempted"] >= LIMITS["scoremol_per_call"]:
            self.audit["truncated"] = True
            self.audit["reasons"].append("scoremol_cap")
            return None
        ref, cand = Chem.Mol(conformers._heavy_view(reference)), Chem.Mol(conformers._heavy_view(candidate))
        before = pose_binding(ref, rid), pose_binding(cand, cid)
        opts = [conformers.rdShapeAlign.ShapeInputOptions() for _ in range(2)]
        for option in opts:
            option.normalize, option.useColors = False, True
        if not self.available():
            return None
        self.audit["scoremol_attempted"] += 1
        scores = conformers.rdShapeAlign.ScoreMol(ref, cand, *opts, rid, cid)
        self.audit["scoremol_returned"] += 1
        if before != (pose_binding(ref, rid), pose_binding(cand, cid)):
            raise ValueError("Fixed-pose scoring changed coordinates or identity")
        if len(scores) != 2 or not all(math.isfinite(x) and 0 <= x <= 1 for x in scores):
            raise ValueError("Invalid fixed-pose shape/color score")
        return tuple(float(x) for x in scores) if self.available() else None

    def _maps(self, ref, cand, requests, row):
        required = [
            next(
                i
                for i, f in enumerate(ref)
                if f["family"] == r["family"] and f["atom_ids"] == r["hit_atom_ids"]
            )
            for r in requests
        ]
        required.sort()
        choices = {
            i: sorted(
                (j for j, f in enumerate(cand) if f["family"] == ref[i]["family"]),
                key=lambda j: (sum((a - b) ** 2 for a, b in zip(ref[i]["position"], cand[j]["position"])), j),
            )
            for i in required
        }
        maps = []

        def walk(pairs, used):
            if not self.available():
                return
            if len(pairs) == len(required):
                maps.append(list(pairs))
                return
            i = required[len(pairs)]
            for j in choices[i]:
                if (
                    len(maps) >= LIMITS["complete_maps_per_pair"]
                    or row["mapping_edge_visits"] >= LIMITS["mapping_edge_visits_per_pair"]
                ):
                    row["mapping_truncated"] = self.audit["truncated"] = True
                    return
                row["mapping_edge_visits"] += 1
                if j not in used:
                    walk(pairs + [(i, j)], used | {j})

        walk([], set())
        row["complete_maps"] = len(maps)
        return maps

    def _anchors(self, mapping, ref, cand, row):
        yield mapping
        if len(mapping) != 1:
            return
        i, j = mapping[0]
        # Prefer a remote typed anchor; ties use feature identity order, not scores.
        helpers = []
        for a, rf in enumerate(ref):
            for b, cf in enumerate(cand):
                if row["helper_edge_visits"] >= LIMITS["helper_edge_visits_per_pair"]:
                    row["helper_truncated"] = self.audit["truncated"] = True
                    break
                row["helper_edge_visits"] += 1
                if not self.available():
                    return
                if a == i or b == j or rf["family"] != cf["family"]:
                    continue
                dr = float(np.linalg.norm(np.asarray(rf["position"]) - ref[i]["position"]))
                dc = float(np.linalg.norm(np.asarray(cf["position"]) - cand[j]["position"]))
                if min(dr, dc) > 1e-10:
                    helpers.append((-dr, abs(dr - dc), a, b))
        if helpers:
            _, _, a, b = min(helpers)
            yield mapping + [(a, b)]

    def propose(self, reference, aligned, rid, cid, baseline_info, witness, similarity, settings):
        self.active = self.clock()
        row = {
            "reference_conformer_id": rid,
            "candidate_conformer_id": cid,
            "status": "not_searched",
            "mapping_edge_visits": 0,
            "mapping_truncated": False,
            "helper_edge_visits": 0,
            "helper_truncated": False,
            "proposals": [],
        }
        self.audit["pairs"].append(row)
        output = []
        try:
            if not self.available():
                return output
            if self.audit["alternatives_attempted"] >= LIMITS["alternatives_per_call"]:
                self.audit["truncated"] = True
                row["status"] = "proposal_cap"
                return output
            baseline_pose, reference_pose = pose_binding(aligned, cid), pose_binding(reference, rid)
            row["baseline_pose_sha256"] = baseline_pose["sha256"]
            row["reference_pose_sha256"] = reference_pose["sha256"]
            parity = self._score(reference, aligned, rid, cid)
            if parity is None:
                row["status"] = "budget_exhausted"
                return output
            row["parity"] = {
                "shape": abs(parity[0] - baseline_info["shape"]),
                "color": abs(parity[1] - baseline_info["alignment_color"]),
            }
            if max(row["parity"].values()) > PARITY_TOLERANCE:
                row["status"] = "fixed_score_parity_unavailable"
                return output
            ref, cand = witness["reference_features"], witness["candidate_features"]
            seen = set()
            for mapping in self._maps(ref, cand, similarity.policy["landmarks"], row):
                for anchors in self._anchors(mapping, ref, cand, row):
                    if (
                        len(row["proposals"]) >= LIMITS["alternatives_per_pair"]
                        or self.audit["alternatives_attempted"] >= LIMITS["alternatives_per_call"]
                    ):
                        row["proposal_truncated"] = self.audit["truncated"] = True
                        return output
                    if not self.available():
                        return output
                    rotation, translation = proper_transform(
                        [cand[j]["position"] for _, j in anchors], [ref[i]["position"] for i, _ in anchors]
                    )
                    transform = {"rotation": rotation.tolist(), "translation": translation.tolist()}
                    key = binding(transform)
                    if key in seen:
                        continue
                    seen.add(key)
                    entry = {
                        "index": len(row["proposals"]),
                        "status": "started",
                        "transform": transform,
                        "mandatory_correspondences": [
                            {"reference_feature": i, "candidate_feature": j} for i, j in mapping
                        ],
                        "proposal_anchors": [
                            {"reference_feature": i, "candidate_feature": j} for i, j in anchors
                        ],
                        "scope": "Proposal correspondences only; final joint assignment is independently recomputed",
                    }
                    row["proposals"].append(entry)
                    self.audit["alternatives_attempted"] += 1
                    candidate = Chem.Mol(aligned)
                    coords = np.asarray(aligned.GetConformer(cid).GetPositions()) @ rotation.T + translation
                    for a, point in enumerate(coords):
                        candidate.GetConformer(cid).SetAtomPosition(a, point.tolist())
                    saved = pose_binding(candidate, cid)
                    if saved["smiles"] != baseline_pose["smiles"] or indexed_topology(
                        candidate
                    ) != indexed_topology(aligned):
                        raise ValueError("Rigid proposal changed exact molecular identity")
                    scores = self._score(reference, candidate, rid, cid)
                    if scores is None:
                        entry["status"] = "budget_exhausted"
                        return output
                    pharm, counts = conformers.pharmacophore_overlap(
                        reference, candidate, rid, cid, settings.pharmacophore_sigma
                    )
                    final = similarity.assess(reference, candidate, rid, cid)
                    joint = math.sqrt(scores[0] * pharm)
                    info = {
                        "shape": scores[0],
                        "alignment_color": scores[1],
                        "pharmacophore": pharm,
                        "joint": joint,
                        "reference_adjusted_joint": joint
                        * math.exp(-settings.strain_penalty * baseline_info["reference_relative_energy"]),
                        "reference_conformer_id": rid,
                        "reference_relative_energy": baseline_info["reference_relative_energy"],
                        **counts,
                    }
                    self.audit["alternatives_returned"] += 1
                    entry.update(
                        status="evaluated",
                        candidate_pose_sha256=saved["sha256"],
                        reference_pose_sha256=reference_pose["sha256"],
                        shape=scores[0],
                        pharmacophore=pharm,
                        witness_reasons=list(final["reasons"]),
                        polar_recall=final["polar"]["recall"],
                    )
                    if not self.available():
                        entry["status"] = "discarded_after_allowance"
                        return output
                    if (
                        pose_binding(aligned, cid) != baseline_pose
                        or pose_binding(reference, rid) != reference_pose
                    ):
                        raise ValueError("Proposal mutated baseline/reference")
                    info["landmark_proposal"] = {
                        "version": VERSION,
                        "pair_index": len(self.audit["pairs"]) - 1,
                        "proposal_index": entry["index"],
                        "transform": transform,
                        "baseline_candidate_pose": baseline_pose,
                        "candidate_pose_sha256": saved["sha256"],
                        "reference_pose_sha256": reference_pose["sha256"],
                        "proposal_policy_sha256": binding(self.policy),
                    }
                    output.append((info, candidate, final))
            row["status"] = (
                "searched"
                if row["proposals"]
                else ("mapping_truncated" if row["mapping_truncated"] else "no_complete_correspondence")
            )
            return output
        except (ValueError, RuntimeError, OverflowError) as exc:
            row.update(status="unavailable", error=f"{type(exc).__name__}: {exc}")
            return []
        finally:
            if not self.available() and row["status"] == "not_searched":
                row["status"] = "budget_exhausted"
            self.spent = self.elapsed()
            self.active = None
            self.audit["elapsed_seconds"] = self.spent
            self.audit["cooperative_overrun"] = self.spent > LIMITS["cooperative_seconds_per_call"]


def validate_selected(info, aligned, similarity):
    """Check saved proposal association/rigid transform; no new scoring or search."""
    from .similarity_constraints import restore_pose

    policy = similarity.policy["landmark_pose_proposals"]
    audit = info["pose_selection"]["landmark_pose_proposals"]
    if policy != declaration() or audit["policy"] != policy:
        raise ValueError("Landmark proposal policy/source mismatch")
    for field, cap in (
        ("alternatives_attempted", LIMITS["alternatives_per_call"]),
        ("scoremol_attempted", LIMITS["scoremol_per_call"]),
    ):
        if type(audit[field]) is not int or not 0 <= audit[field] <= cap:
            raise ValueError("Landmark proposal call bound mismatch")
    cid, rid = info["candidate_conformer_id"], info["hit"]["reference_conformer_id"]
    pairs = [
        p
        for p in info["pose_selection"]["pairs"]
        if p["candidate_conformer_id"] == cid and p["reference_conformer_id"] == rid
    ]
    if len(pairs) != 1 or pairs[0]["status"] != "feasible":
        raise ValueError("Missing selected conformer pair")
    expected_decision = {
        "status": "feasible",
        "reasons": [],
        "score": info["score"],
        "shape": info["hit"]["shape"],
        "pharmacophore": info["hit"]["pharmacophore"],
        "polar_recall": info["similarity_constraints"]["witness"]["polar"]["recall"],
    }
    if any(pairs[0].get(key) != value for key, value in expected_decision.items()):
        raise ValueError("Selected pair decision/metric mismatch")
    actual = pose_binding(aligned, cid)
    if pairs[0]["selected_pose_sha256"] != actual["sha256"]:
        raise ValueError("Selected pair pose mismatch")
    saved = info["hit"].get("landmark_proposal")
    if saved is None:
        if pairs[0]["selected_origin"] != "baseline":
            raise ValueError("Missing selected proposal origin")
        return
    if (
        pairs[0]["selected_origin"] != "landmark_proposal"
        or saved["version"] != VERSION
        or saved["proposal_policy_sha256"] != binding(policy)
    ):
        raise ValueError("Selected proposal origin/policy mismatch")
    if (
        type(saved["pair_index"]) is not int
        or type(saved["proposal_index"]) is not int
        or min(saved["pair_index"], saved["proposal_index"]) < 0
    ):
        raise ValueError("Invalid selected proposal index")
    if saved["pair_index"] >= len(audit["pairs"]):
        raise ValueError("Selected proposal pair index out of range")
    row = audit["pairs"][saved["pair_index"]]
    if saved["proposal_index"] >= len(row["proposals"]):
        raise ValueError("Selected proposal index out of range")
    entry = row["proposals"][saved["proposal_index"]]
    reference = info["similarity_constraints"]["witness"]["reference_pose"]
    if (
        row["reference_conformer_id"] != rid
        or row["candidate_conformer_id"] != cid
        or entry["status"] != "evaluated"
        or entry["candidate_pose_sha256"] != actual["sha256"]
        or saved["candidate_pose_sha256"] != actual["sha256"]
        or saved["reference_pose_sha256"] != reference["sha256"]
        or entry["reference_pose_sha256"] != reference["sha256"]
        or entry["transform"] != saved["transform"]
    ):
        raise ValueError("Selected proposal does not bind actual saved pose")
    baseline = restore_pose(saved["baseline_candidate_pose"])
    before = pose_binding(baseline, cid)
    if (
        before["sha256"] != row["baseline_pose_sha256"]
        or before["smiles"] != actual["smiles"]
        or before["atom_order_sha256"] != actual["atom_order_sha256"]
    ):
        raise ValueError("Proposal baseline identity mismatch")
    rotation = np.asarray(saved["transform"]["rotation"], dtype=float)
    translation = np.asarray(saved["transform"]["translation"], dtype=float)
    if (
        rotation.shape != (3, 3)
        or translation.shape != (3,)
        or not np.isfinite(rotation).all()
        or not np.isfinite(translation).all()
        or not np.allclose(rotation.T @ rotation, np.eye(3), rtol=0, atol=1e-8)
        or not math.isclose(float(np.linalg.det(rotation)), 1, abs_tol=1e-8)
    ):
        raise ValueError("Saved proposal is not a proper rigid transform")
    expected = np.asarray(before["coordinates"]) @ rotation.T + translation
    if not np.allclose(expected, actual["coordinates"], rtol=0, atol=1e-8):
        raise ValueError("Saved proposal transform/coordinates mismatch")
    features = info["similarity_constraints"]["witness"]
    ref, cand = features["reference_features"], features["candidate_features"]
    mandatory = entry["mandatory_correspondences"]
    anchors = entry["proposal_anchors"]
    required = {(r["family"], tuple(r["hit_atom_ids"])) for r in similarity.policy["landmarks"]}

    def identities(items):
        left, right = [], []
        for pair in items:
            i, j = pair["reference_feature"], pair["candidate_feature"]
            if (
                type(i) is not int
                or type(j) is not int
                or not 0 <= i < len(ref)
                or not 0 <= j < len(cand)
                or ref[i]["family"] != cand[j]["family"]
            ):
                raise ValueError("Invalid typed proposal correspondence")
            left.append(i)
            right.append(j)
        if len(left) != len(set(left)) or len(right) != len(set(right)):
            raise ValueError("Noninjective proposal correspondence")
        return {(ref[i]["family"], tuple(ref[i]["atom_ids"])) for i in left}

    if (
        identities(mandatory) != required
        or not required <= identities(anchors)
        or anchors[: len(mandatory)] != mandatory
    ):
        raise ValueError("Proposal does not cover every mandatory landmark")
    if entry["shape"] != info["hit"]["shape"] or entry["pharmacophore"] != info["hit"]["pharmacophore"]:
        raise ValueError("Selected proposal metric association mismatch")
