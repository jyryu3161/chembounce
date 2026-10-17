"""Bounded joint-feasible original-hit pose selection, separate from v1 kernels.

These are explicit engineering constraints, not activity or synthesis estimates.
This scorer does not itself certify a scaffold hop or a synthesis route.
"""

import math
import time
from dataclasses import asdict, replace

from rdkit import DataStructs

from .common import absolute_smiles, FP, canonical, molecule
from .conformers import ConformerSettings, _overlay, ensemble, pharmacophore_overlap, require_alignment
from .similarity_constraints import VERSION as CONSTRAINT_VERSION
from .similarity_constraints import BoundSimilarity, binding, feature_records, pose_binding, restore_pose

VERSION = "joint-feasible-hit-pose-v2"
MAX_CONFORMERS = 8
MAX_HIT_PAIRS = 64
DEFAULT_SETTINGS = ConformerSettings(num_conformers=8, scoring_objective="original-hit-v1")


def configured_constraints(config, hit=None):
    """Validate opt-in orchestration before assets; bind query selectors when known."""
    from .refinement import RefinementSettings

    if "similarity_constraints" not in config:
        return None
    request = config["similarity_constraints"]
    if (
        not isinstance(request, dict)
        or set(request) - {"version", "polar_recall", "landmarks", "landmark_pose_proposals"}
        or request.get("version") != CONSTRAINT_VERSION
    ):
        raise ValueError("similarity_constraints requires the selective-similarity-constraints-v2 schema")
    enabled = request.get("landmark_pose_proposals", False)
    if type(enabled) is not bool:
        raise ValueError("landmark_pose_proposals must be true or false")
    if enabled:
        if not request.get("landmarks"):
            raise ValueError("Landmark pose proposals require explicit nonempty landmarks")
        from .landmark_proposals import require_capability

        require_capability()
    settings = ConformerSettings(**config.get("conformers", {}))
    refinement = RefinementSettings(**config.get("refinement", {}))
    if (
        config.get("recommendation", {}).get("mode") != "qualified"
        or config.get("mode", "hop") != "hop"
        or settings.num_conformers != 8
        or settings.seed != 42
        or settings.scoring_objective != "original-hit-v1"
        or settings.pharmacophore_sigma != 1.0
        or not settings.reject_unspecified_stereo
        or refinement.top_n <= 0
        or refinement.num_conformers != 8
        or tuple(refinement.seeds) != (42, 43, 44)
        or refinement.threshold_policy != "all_seeds"
        or refinement.adaptive_max_conformers is not None
        or config.get("reference_sdf")
        or "pocket" in config
    ):
        raise ValueError(
            "Selective constraints require qualified hop mode, fixed8 seed42 screening and complete42/43/44 refinement, original-hit-v1/sigma1/strict stereo, with no adaptive, reference-SDF or pocket override"
        )
    minimum = {
        key: config.get("minimum_scores", {}).get(key)
        for key in ("hit_tanimoto", "hit_shape", "hit_pharmacophore")
    }
    if hit is not None:
        return BoundSimilarity(
            hit,
            minimum,
            polar_recall=request.get("polar_recall", 0.75),
            core_smiles=config.get("core_smiles"),
            landmarks=request.get("landmarks", []),
            landmark_pose_proposals=enabled,
        )
    return None


def failure_record(exc, *, phase, result=None):
    audit = getattr(exc, "audit", None)
    if audit is None and result is not None:
        audit = result[0].get("pose_selection")
    return {
        "phase": phase,
        "reason": str(exc),
        "reason_code": getattr(exc, "reason", None),
        "ensemble_request_status": "none_requested"
        if phase == "preflight" or (isinstance(audit, dict) and audit.get("ensemble_requests") == [])
        else "requested"
        if audit is not None
        else "unknown",
        "scorer_audit": audit,
    }


def result_fields(info):
    witness = info["similarity_constraints"]["witness"]
    return {
        "pose_selection_policy": VERSION,
        "similarity_constraint_policy": CONSTRAINT_VERSION,
        "similarity_constraint_sha256": witness["policy_sha256"],
        "similarity_witness_sha256": binding(witness),
        "polar_reference_recall": witness["polar"]["recall"],
    }


def preflight_design(design):
    """Check secondary-reference eligibility without imposing candidate cores."""
    mol = molecule(canonical(design))
    BoundSimilarity._check_state(mol)
    feature_records(mol)


def restore_result(info):
    """Materialize saved seed coordinates for consistency checking, never scoring."""
    try:
        poses = info["selective_pose_evidence"]
        return info, restore_pose(poses["hit_candidate"]), restore_pose(poses["design_candidate"])
    except (KeyError, TypeError) as exc:
        raise ValueError(f"Missing retained selective seed geometry: {exc}") from exc


class SimilarityAbstention(ValueError):
    """No supported result within the specified, bounded scoring protocol."""

    def __init__(self, reason, audit):
        super().__init__(reason)
        self.reason = reason
        self.audit = audit


def _validate_ensemble(value, smiles):
    if canonical(absolute_smiles(value.mol)) != smiles:
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
    rejected_pose_diagnostics=None,
    rejected_pose_diagnostic_phase="direct",
):
    """Enumerate at most 8x8 hit pairs, filter jointly, then maximize v1 scalar.

    The scalar includes both reference and candidate relative-energy penalties.
    Strict ``>`` preserves the first pair on ties. Real-design diagnostics use the
    same selected candidate conformer and never influence its hit-pose selection.
    Reference SDF/pocket overrides are intentionally outside this version's scope.
    Deadline checks bound further work, but cannot interrupt an individual RDKit call.
    """
    started = time.monotonic()
    if rejected_pose_diagnostics is not None:
        diagnostic_dispatch_started = time.monotonic()
        from .rejected_pose_diagnostics import RejectedPoseDiagnostics

        if rejected_pose_diagnostic_phase not in {"direct", "screening", "refinement"}:
            raise ValueError("Unknown rejected-pose diagnostic phase")
        if not isinstance(rejected_pose_diagnostics, RejectedPoseDiagnostics):
            raise ValueError("rejected_pose_diagnostics requires a shared diagnostic collector")
        rejected_pose_diagnostics.charge_overhead(time.monotonic() - diagnostic_dispatch_started)
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
    preflight_design(design)
    policy = similarity.policy
    proposals = None
    if "landmark_pose_proposals" in policy:
        proposal_setup = time.monotonic()
        from .landmark_proposals import ProposalSearch

        proposals = ProposalSearch(
            policy["landmark_pose_proposals"], deadline=deadline_monotonic, setup_started=proposal_setup
        )
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

    if proposals is not None:
        audit["landmark_pose_proposals"] = proposals.audit
        audit["landmark_pose_proposals"]["bound_policy_preflight_seconds"] = (
            similarity.landmark_proposal_preflight_seconds
        )
        audit["landmark_pose_proposals"]["bound_policy_preflight_scope"] = (
            "Previously incurred once per bound object; reported here, not incurred again per seed"
        )
    diagnostic_dispatch_started = time.monotonic() if rejected_pose_diagnostics is not None else None
    diagnostic_before = rejected_pose_diagnostics._elapsed if rejected_pose_diagnostics is not None else 0
    diagnostic = (
        rejected_pose_diagnostics.start_call(
            {"product": smiles, "hit": hit, "design": design},
            settings,
            policy,
            deadline_monotonic,
            rejected_pose_diagnostic_phase,
        )
        if rejected_pose_diagnostics is not None
        else None
    )

    if rejected_pose_diagnostics is not None:
        rejected_pose_diagnostics.charge_overhead(
            max(
                0.0,
                time.monotonic()
                - diagnostic_dispatch_started
                - (rejected_pose_diagnostics._elapsed - diagnostic_before),
            )
        )

    def finish(status, reason=None):
        audit.update(status=status, reason=reason, elapsed_seconds=time.monotonic() - started)

    def abstain(reason):
        # Scientific abstention is decided before optional serialization.
        if diagnostic is not None:
            audit["rejected_pose_diagnostic"] = diagnostic.finish(reason)
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
                    or canonical(absolute_smiles(aligned)) != smiles
                ):
                    raise ValueError("Invalid pair overlay or chemical identity")
                witness = similarity.assess(hit_ref.mol, aligned, rid, cid)
            except (ValueError, RuntimeError) as exc:
                pair.update(status="failed", reasons=[str(exc)])
                check_deadline()
                continue
            audit["hit_pairs_completed"] += 1
            proposal_work_started = time.monotonic() if proposals is not None else None
            proposal_spent_before = proposals.spent if proposals is not None else 0
            options = [(hit_info, aligned, witness)]
            if (
                proposals is not None
                and "required_landmark_joint_assignment_infeasible" in witness["reasons"]
            ):
                options += proposals.propose(
                    hit_ref.mol, aligned, rid, cid, hit_info, witness, similarity, settings
                )
                check_deadline()
            pair_best = None
            selected_outcome = None
            outcomes = []
            for option_index, (option_info, option_pose, option_witness) in enumerate(options):
                reasons = list(option_witness["reasons"])
                for name, field in (("hit_shape", "shape"), ("hit_pharmacophore", "pharmacophore")):
                    if option_info[field] < policy["minimum_scores"][name]:
                        reasons.append(name + "_below_minimum")
                score = (
                    math.sqrt(option_info["shape"] * option_info["pharmacophore"])
                    * math.exp(-settings.strain_penalty * reference_energy)
                    * math.exp(-settings.strain_penalty * candidate_energy)
                )
                outcome = dict(
                    status="rejected" if reasons else "feasible",
                    reasons=reasons,
                    score=score,
                    shape=option_info["shape"],
                    pharmacophore=option_info["pharmacophore"],
                    polar_recall=option_witness["polar"]["recall"],
                )
                outcomes.append(outcome)
                if reasons and diagnostic is not None:
                    diagnostic_pair = {**pair, **outcome}
                    if option_index:
                        diagnostic_pair["landmark_proposal"] = option_info["landmark_proposal"]
                    diagnostic.consider(
                        len(audit["pairs"]) - 1,
                        diagnostic_pair,
                        option_pose,
                        option_witness,
                        option_info,
                        candidate_energy,
                        reference_energy,
                        hit_tanimoto=tanimoto,
                    )
                if not reasons and (pair_best is None or score > pair_best[0]):
                    pair_best = (score, cid, option_info, option_pose, option_witness)
                    selected_outcome = outcome
            pair.update(outcomes[0])
            if proposals is not None:
                pair["baseline"] = dict(outcomes[0])
                pair["proposal_outcomes"] = outcomes[1:]
            if proposals is not None:
                proposals.charge(
                    time.monotonic() - proposal_work_started - (proposals.spent - proposal_spent_before)
                )
                if (
                    pair_best is not None
                    and "landmark_proposal" in pair_best[2]
                    and not proposals.available()
                ):
                    pair["proposal_selection_status"] = "discarded_after_allowance"
                    pair_best = None
            if pair_best is not None:
                audit["hit_pairs_feasible"] += 1
                if proposals is not None:
                    pair.update(selected_outcome)
                    pair["selected_pose_sha256"] = pair_best[4]["candidate_pose_sha256"]
                    pair["selected_origin"] = (
                        "landmark_proposal" if "landmark_proposal" in pair_best[2] else "baseline"
                    )
                if best is None or pair_best[0] > best[0]:
                    best = pair_best
            check_deadline()

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
    if diagnostic is not None:
        diagnostic.finish("selected")
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
        "selective_pose_evidence": {
            "hit_candidate": pose_binding(aligned, cid),
            "design_candidate": pose_binding(design_aligned, cid),
            "design_reference": pose_binding(design_ref.mol, design_info["reference_conformer_id"]),
        },
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
        if "landmark_pose_proposals" in similarity.policy:
            from .landmark_proposals import validate_selected

            validate_selected(info, aligned, similarity)
        elif "landmark_proposal" in info["hit"] or "landmark_pose_proposals" in info["pose_selection"]:
            raise ValueError("Unrequested landmark proposal evidence")
        for role, pose in (("hit_candidate", aligned), ("design_candidate", result[2])):
            if binding(info["selective_pose_evidence"][role]) != binding(pose_binding(pose, cid)):
                raise ValueError("Retained selective seed geometry does not describe the actual pose")
        reference = similarity.validate_witness(witness, aligned, cid)
        if rid != witness["reference_pose"]["conformer_id"]:
            raise ValueError("Hit metric and joint witness reference pose differ")
        pharm, _ = pharmacophore_overlap(reference, aligned, rid, cid, settings.pharmacophore_sigma)
        if not math.isclose(pharm, info["hit"]["pharmacophore"], abs_tol=1e-12, rel_tol=1e-12):
            raise ValueError("Hit pharmacophore metric does not match selected coordinates")
        design_reference = restore_pose(info["selective_pose_evidence"]["design_reference"])
        if canonical(absolute_smiles(design_reference)) != design:
            raise ValueError("Retained design reference identity changed")
        design_pharm, _ = pharmacophore_overlap(
            design_reference,
            result[2],
            info["design"]["reference_conformer_id"],
            cid,
            settings.pharmacophore_sigma,
        )
        if not math.isclose(design_pharm, info["design"]["pharmacophore"], abs_tol=1e-12, rel_tol=1e-12):
            raise ValueError("Design pharmacophore metric does not match selected coordinates")
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
