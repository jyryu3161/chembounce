"""Bounded joint-feasible original-hit pose selection, separate from v1 kernels.

These are explicit engineering constraints, not activity or synthesis estimates.
This scorer does not itself certify a scaffold hop or a synthesis route.
"""

import math
import time
from dataclasses import asdict, replace

from rdkit import Chem, DataStructs

from .common import FP, canonical, molecule
from .conformers import ConformerSettings, _overlay, ensemble, pharmacophore_overlap, require_alignment
from .similarity_constraints import BoundSimilarity, feature_records

VERSION = "joint-feasible-hit-pose-v2"
MAX_CONFORMERS = 8
MAX_HIT_PAIRS = 64
DEFAULT_SETTINGS = ConformerSettings(num_conformers=8, scoring_objective="original-hit-v1")


class SimilarityAbstention(ValueError):
    """No supported result within the specified, bounded scoring protocol."""

    def __init__(self, reason, audit):
        super().__init__(reason)
        self.reason = reason
        self.audit = audit


def _validate_ensemble(value, smiles):
    if canonical(Chem.MolToSmiles(value.mol)) != smiles:
        raise ValueError("Ensemble chemical identity changed")
    keys = set(value.relative_energies)
    if (
        not keys
        or len(keys) > MAX_CONFORMERS
        or keys != set(value.energies)
        or keys != {c.GetId() for c in value.mol.GetConformers()}
        or value.bound_pose
    ):
        raise ValueError("Invalid or oversized generated ensemble")
    for cid in keys:
        energy, relative = value.energies[cid], value.relative_energies[cid]
        if (
            not all(
                isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)
                for x in (energy, relative)
            )
            or relative < 0
        ):
            raise ValueError("Invalid ensemble energy")


def score_candidate(
    smiles,
    hit,
    design,
    settings=DEFAULT_SETTINGS,
    *,
    similarity=None,
    deadline_monotonic=None,
):
    """Enumerate at most 8x8 hit pairs, filter jointly, then maximize v1 scalar.

    The scalar includes both reference and candidate relative-energy penalties.
    Strict ``>`` preserves the first pair on ties. Real-design diagnostics use the
    same selected candidate conformer and never influence its hit-pose selection.
    Reference SDF/pocket overrides are intentionally outside this version's scope.
    Deadline checks bound further work, but cannot interrupt an individual RDKit call.
    """
    started = time.monotonic()
    if (
        not isinstance(settings, ConformerSettings)
        or settings.num_conformers != MAX_CONFORMERS
        or settings.scoring_objective != "original-hit-v1"
        or settings.pharmacophore_sigma != 1.0
        or not settings.reject_unspecified_stereo
    ):
        raise ValueError(
            "Selective scoring requires fixed 8 conformers, original-hit-v1, sigma 1 and strict stereo"
        )
    if deadline_monotonic is not None and (
        isinstance(deadline_monotonic, bool)
        or not isinstance(deadline_monotonic, (int, float))
        or not math.isfinite(deadline_monotonic)
    ):
        raise ValueError("Invalid scoring deadline")
    smiles, hit, design = canonical(smiles), canonical(hit), canonical(design)
    similarity = BoundSimilarity(hit) if similarity is None else similarity
    if not isinstance(similarity, BoundSimilarity) or similarity.hit != hit:
        raise ValueError("Selective constraint policy must bind the exact original hit")
    similarity.preflight_candidate(smiles)
    # Design availability remains a secondary-diagnostic prerequisite. This
    # does not require the design itself to satisfy candidate core/feature gates.
    BoundSimilarity._check_state(molecule(design))
    feature_records(molecule(design))
    policy = similarity.policy
    audit = {
        "policy": VERSION,
        "constraint_policy_sha256": similarity.sha256,
        "settings": asdict(settings),
        "requested_conformers_per_ensemble": MAX_CONFORMERS,
        "max_hit_pairs": MAX_HIT_PAIRS,
        "ensemble_requests": [],
        "hit_pairs_attempted": 0,
        "hit_pairs_completed": 0,
        "hit_pairs_feasible": 0,
        "design_overlay_attempts": 0,
        "pairs": [],
        "status": "running",
        "scope": "Bounded generated ensembles; not an exhaustive conformational search. Cached ensemble requests are recorded, not asserted to be fresh embeddings.",
    }

    def finish(status, reason=None):
        audit.update(status=status, reason=reason, elapsed_seconds=time.monotonic() - started)

    def abstain(reason):
        finish("abstained", reason)
        raise SimilarityAbstention(reason, audit)

    def check_deadline():
        if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
            abstain("scoring_deadline_exhausted")

    check_deadline()
    tanimoto = float(
        DataStructs.TanimotoSimilarity(FP.GetFingerprint(molecule(hit)), FP.GetFingerprint(molecule(smiles)))
    )
    audit["hit_tanimoto"] = tanimoto
    if tanimoto < policy["minimum_scores"]["hit_tanimoto"]:
        abstain("hit_tanimoto_below_minimum")
    require_alignment()

    def get_ensemble(role, identity):
        check_deadline()
        requested_at = time.monotonic()
        request = {"role": role, "smiles": identity, "status": "requested"}
        audit["ensemble_requests"].append(request)
        try:
            value = ensemble(identity, settings)
            _validate_ensemble(value, identity)
        except (ValueError, RuntimeError) as exc:
            request.update(status="failed", elapsed_seconds=time.monotonic() - requested_at)
            abstain(f"{role}_ensemble_failed: {exc}")
        request.update(
            status="returned",
            elapsed_seconds=time.monotonic() - requested_at,
            retained_conformers=len(value.relative_energies),
        )
        check_deadline()
        return value

    candidate = get_ensemble("candidate", smiles)
    hit_ref = get_ensemble("hit", hit)
    best = None
    for cid, candidate_energy in candidate.relative_energies.items():
        for rid, reference_energy in hit_ref.relative_energies.items():
            check_deadline()
            pair = {"candidate_conformer_id": cid, "reference_conformer_id": rid}
            audit["pairs"].append(pair)
            audit["hit_pairs_attempted"] += 1
            # A single-reference view prevents the old scalar-optimal reference
            # from hiding another jointly feasible reference in the SAME ensemble.
            reference = replace(hit_ref, relative_energies={rid: reference_energy})
            try:
                hit_info, aligned = _overlay(reference, candidate, cid, settings)
                if (
                    hit_info["reference_conformer_id"] != rid
                    or not all(
                        math.isfinite(hit_info[k]) and 0 <= hit_info[k] <= 1
                        for k in ("shape", "pharmacophore")
                    )
                    or canonical(Chem.MolToSmiles(aligned)) != smiles
                ):
                    raise ValueError("Invalid pair overlay or chemical identity")
                witness = similarity.assess(hit_ref.mol, aligned, rid, cid)
            except (ValueError, RuntimeError) as exc:
                pair.update(status="failed", reasons=[str(exc)])
                check_deadline()
                continue
            audit["hit_pairs_completed"] += 1
            reasons = list(witness["reasons"])
            for name, field in (("hit_shape", "shape"), ("hit_pharmacophore", "pharmacophore")):
                if hit_info[field] < policy["minimum_scores"][name]:
                    reasons.append(name + "_below_minimum")
            score = (
                math.sqrt(hit_info["shape"] * hit_info["pharmacophore"])
                * math.exp(-settings.strain_penalty * reference_energy)
                * math.exp(-settings.strain_penalty * candidate_energy)
            )
            pair.update(
                status="rejected" if reasons else "feasible",
                reasons=reasons,
                score=score,
                shape=hit_info["shape"],
                pharmacophore=hit_info["pharmacophore"],
                polar_recall=witness["polar"]["recall"],
            )
            if not reasons:
                audit["hit_pairs_feasible"] += 1
            check_deadline()
            if reasons:
                continue
            if best is None or score > best[0]:
                best = (score, cid, hit_info, aligned, witness)
    if best is None:
        abstain("no_jointly_feasible_hit_pose")
    score, cid, hit_info, aligned, witness = best
    design_ref = get_ensemble("design", design)
    check_deadline()
    audit["design_overlay_attempts"] = 1
    try:
        design_info, design_aligned = _overlay(design_ref, candidate, cid, settings)
    except (ValueError, RuntimeError) as exc:
        abstain(f"selected_conformer_design_diagnostic_failed: {exc}")
    check_deadline()
    finish("selected")
    info = {
        "score": score,
        "scoring_objective": settings.scoring_objective,
        "pose_selection_policy": VERSION,
        "candidate_smiles": smiles,
        "hit_smiles": hit,
        "design_smiles": design,
        "candidate_conformer_id": cid,
        "candidate_energy": candidate.energies[cid],
        "candidate_relative_energy": candidate.relative_energies[cid],
        "force_field": candidate.force_field,
        "hit": hit_info,
        "design": design_info,
        "ensemble_sizes": {
            "candidate": len(candidate.relative_energies),
            "hit": len(hit_ref.relative_energies),
            "design": len(design_ref.relative_energies),
        },
        "ensemble_diagnostics": {
            "candidate": candidate.failures,
            "hit": hit_ref.failures,
            "design": design_ref.failures,
        },
        "settings": asdict(settings),
        "similarity_constraints": {"policy": policy, "witness": witness},
        "pose_selection": audit,
    }
    return info, aligned, design_aligned


def validate_result(result, product, settings, *, hit, design, similarity):
    """Bind joint feasibility to actual exported coordinates and expected policy.

    Reuses v1 validation and recomputes the added geometric witness/global
    pharmacophore metric. This does not attest an external claimed search log.
    """
    from .qualification import validate_result as validate_base

    validate_base(result, product, settings, hit=hit, design=design)
    try:
        info, aligned, _ = result
        if (
            not isinstance(similarity, BoundSimilarity)
            or similarity.hit != hit
            or info["pose_selection_policy"] != VERSION
            or info["similarity_constraints"]["policy"] != similarity.policy
            or settings.num_conformers != MAX_CONFORMERS
            or settings.scoring_objective != "original-hit-v1"
            or settings.pharmacophore_sigma != 1.0
            or not settings.reject_unspecified_stereo
        ):
            raise ValueError("Selective result policy mismatch")
        witness = info["similarity_constraints"]["witness"]
        cid, rid = info["candidate_conformer_id"], info["hit"]["reference_conformer_id"]
        reference = similarity.validate_witness(witness, aligned, cid)
        if rid != witness["reference_pose"]["conformer_id"]:
            raise ValueError("Hit metric and joint witness reference pose differ")
        pharm, _ = pharmacophore_overlap(reference, aligned, rid, cid, settings.pharmacophore_sigma)
        if not math.isclose(pharm, info["hit"]["pharmacophore"], abs_tol=1e-12, rel_tol=1e-12):
            raise ValueError("Hit pharmacophore metric does not match selected coordinates")
        actual = {
            "hit_tanimoto": float(
                DataStructs.TanimotoSimilarity(
                    FP.GetFingerprint(molecule(hit)), FP.GetFingerprint(molecule(product))
                )
            ),
            "hit_shape": info["hit"]["shape"],
            "hit_pharmacophore": pharm,
        }
        if any(actual[k] < floor for k, floor in similarity.policy["minimum_scores"].items()):
            raise ValueError("Selective original-hit minimum failed")
        audit = info["pose_selection"]
        if (
            audit["policy"] != VERSION
            or audit["constraint_policy_sha256"] != similarity.sha256
            or audit["status"] != "selected"
            or audit["settings"] != asdict(settings)
            or audit["hit_tanimoto"] != actual["hit_tanimoto"]
        ):
            raise ValueError("Selective search/result provenance mismatch")
    except (KeyError, TypeError, AttributeError, RuntimeError) as exc:
        raise ValueError(f"Incomplete selective scoring result: {exc}") from exc
