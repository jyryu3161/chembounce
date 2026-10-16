"""Energy-filtered ensembles and joint shape/pharmacophore pose selection."""

import math
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, ChemicalFeatures, rdMolAlign, rdMolTransforms

try:
    from rdkit.Chem import rdShapeAlign
except ImportError:
    rdShapeAlign = None
from scipy.optimize import linear_sum_assignment

from .common import canonical, molecule, positive

FEATURES = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
FAMILIES = {"Donor", "Acceptor", "Aromatic", "PosIonizable", "NegIonizable", "Hydrophobe", "LumpedHydrophobe"}


def require_alignment():
    if rdShapeAlign is None:
        raise ValueError(
            "Joint 3D alignment requires an RDKit build with rdShapeAlign; install chembounce[synthesis] (RDKit >=2025.9)"
        )


@dataclass(frozen=True)
class ConformerSettings:
    num_conformers: int = 12
    seed: int = 42
    energy_window: float = 5.0
    strain_penalty: float = 0.05
    max_iterations: int = 500
    pharmacophore_sigma: float = 1.0
    reject_unspecified_stereo: bool = True
    scoring_objective: str = "joint-v1"

    def __post_init__(self):
        if self.scoring_objective not in ("joint-v1", "original-hit-v1"):
            raise ValueError("scoring_objective must be joint-v1 or original-hit-v1")
        if not isinstance(self.reject_unspecified_stereo, bool):
            raise ValueError("reject_unspecified_stereo must be true or false")
        positive(self.num_conformers, "num_conformers")
        positive(self.max_iterations, "max_iterations")
        if isinstance(self.seed, bool) or not isinstance(self.seed, int) or not 0 <= self.seed <= 2**31 - 1:
            raise ValueError("seed must be an integer in [0, 2**31-1]")
        for key in ("energy_window", "strain_penalty", "pharmacophore_sigma"):
            value = getattr(self, key)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
                or (key == "pharmacophore_sigma" and value == 0)
            ):
                raise ValueError(f"Invalid {key}")


@dataclass
class Ensemble:
    mol: object
    energies: dict
    relative_energies: dict
    force_field: str
    failures: dict
    bound_pose: bool = False


def _check_stereo(mol, settings):
    unspecified = [s for s in Chem.FindPotentialStereo(mol) if str(s.specified) == "Unspecified"]
    if unspecified and settings.reject_unspecified_stereo:
        raise ValueError("3D scoring requires specified stereochemistry; enumerate stereoisomers explicitly")
    return len(unspecified)


@lru_cache(maxsize=128)
def ensemble(smiles, settings=ConformerSettings()):
    """Cache by full isomeric/charged SMILES and settings; states are never neutralized."""
    mol = molecule(smiles)
    unspecified = _check_stereo(mol, settings)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = settings.seed
    params.numThreads = 1
    params.pruneRmsThresh = 0.5
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=settings.num_conformers, params=params))
    if not ids:
        raise ValueError("ETKDG produced no conformers")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        force_field = "MMFF94s"
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=force_field)
        make_ff = lambda cid: AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
    elif AllChem.UFFHasAllMoleculeParams(mol):
        force_field = "UFF"
        make_ff = lambda cid: AllChem.UFFGetMoleculeForceField(mol, confId=cid)
    else:
        raise ValueError("No complete MMFF or UFF parameters")
    energies = {}
    failed = 0
    for cid in ids:
        ff = make_ff(cid)
        if ff is None or ff.Minimize(maxIts=settings.max_iterations) != 0:
            failed += 1
            continue
        energy = float(ff.CalcEnergy())
        if math.isfinite(energy):
            energies[cid] = energy
        else:
            failed += 1
    if not energies:
        raise ValueError("No converged finite-energy conformers")
    minimum = min(energies.values())
    kept = {cid: e for cid, e in energies.items() if e - minimum <= settings.energy_window}
    for cid in ids:
        if cid not in kept:
            mol.RemoveConformer(cid)
    mol = Chem.RemoveHs(mol)
    return Ensemble(
        mol,
        kept,
        {cid: e - minimum for cid, e in kept.items()},
        force_field,
        {
            "requested": settings.num_conformers,
            "embedded": len(ids),
            "optimization_failed": failed,
            "outside_energy_window": len(energies) - len(kept),
            "unspecified_stereo": unspecified,
        },
    )


def reference_from_sdf(filename, expected_smiles, settings=ConformerSettings()):
    entries = list(Chem.SDMolSupplier(str(filename), removeHs=False))
    if len(entries) != 1 or entries[0] is None:
        raise ValueError("Reference SDF must contain exactly one valid molecule")
    mol = Chem.RemoveHs(entries[0])
    if canonical(Chem.MolToSmiles(mol)) != canonical(expected_smiles):
        raise ValueError("Reference pose must match the hit's exact chemical state and stereochemistry")
    _check_stereo(mol, settings)
    if not mol.GetNumConformers() or not mol.GetConformer().Is3D():
        raise ValueError("Reference SDF must contain 3D coordinates")
    if not np.isfinite(mol.GetConformer().GetPositions()).all():
        raise ValueError("Reference pose contains nonfinite coordinates")
    # Preserve the experimental pose, including its strain. Do not minimize it.
    cid = mol.GetConformer().GetId()
    return Ensemble(mol, {cid: None}, {cid: 0.0}, "supplied_pose_not_minimized", {}, True)


def prepare_reference(hit, settings=ConformerSettings(), reference_sdf=None):
    """Check fixed 3D prerequisites without embedding or modifying chemical state.

    Return the supplied bound pose, if any, so the workflow can reuse its exact
    coordinates. SDF references retain their existing exact-state/stereo checks;
    generated references use the same stereo policy as ``ensemble``. This does
    not guarantee that later embedding, optimization or alignment will succeed.
    """
    require_alignment()
    if reference_sdf:
        return reference_from_sdf(reference_sdf, hit, settings)
    _check_stereo(molecule(hit), settings)
    return None


def _heavy_view(mol):
    # Ordinary Hs are already removed from optimized ensembles, but RDKit keeps
    # isotopic Hs by default. Explicit D/T must not create extra pharmacophore
    # features or alter shape/color perception solely through graph degree.
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        return mol
    heavy = Chem.RemoveAllHs(mol)
    if not heavy.GetNumAtoms():
        raise ValueError("3D comparison requires at least one heavy atom")
    return heavy


def _features(mol, cid):
    mol = _heavy_view(mol)
    positions = mol.GetConformer(cid).GetPositions()
    result = []
    raw = [f for f in FEATURES.GetFeaturesForMol(mol, confId=cid) if f.GetFamily() in FAMILIES]
    lumped_atoms = {a for f in raw if f.GetFamily() == "LumpedHydrophobe" for a in f.GetAtomIds()}
    for feature in raw:
        family, atoms = feature.GetFamily(), feature.GetAtomIds()
        if family == "Hydrophobe" and set(atoms) <= lumped_atoms:
            continue
        family = "Hydrophobe" if family == "LumpedHydrophobe" else family
        xyz = np.array(feature.GetPos(), dtype=float)
        normal = None
        if family == "Aromatic" and len(atoms) >= 3:
            plane = positions[list(atoms)]
            _, _, vt = np.linalg.svd(plane - plane.mean(axis=0))
            normal = vt[-1]
        result.append((family, xyz, normal))
    return result


def pharmacophore_overlap(reference, probe, ref_cid, probe_cid, sigma=1.0):
    """Typed one-to-one Gaussian feature matching on the SAME aligned pose.

    Aromatic normals are compared without sign. H-bond donor/acceptor positions
    are used; directional H-bond geometry and receptor interactions are not scored.
    """
    ref, cand = _features(reference, ref_cid), _features(probe, probe_cid)
    if not ref or not cand:
        return (1.0 if not ref and not cand else 0.0), {
            "reference_features": len(ref),
            "candidate_features": len(cand),
            "matched_weight": 0.0,
        }
    weights = np.zeros((len(ref), len(cand)))
    for i, (rf, rp, rn) in enumerate(ref):
        for j, (pf, pp, pn) in enumerate(cand):
            if rf == pf:
                value = math.exp(-float(np.sum((rp - pp) ** 2)) / (2 * sigma**2))
                if rn is not None and pn is not None:
                    value *= min(1.0, float(np.dot(rn, pn)) ** 2)
                weights[i, j] = value
    rows, cols = linear_sum_assignment(-weights)
    matched = float(weights[rows, cols].sum())
    return matched / (len(ref) + len(cand) - matched), {
        "reference_features": len(ref),
        "candidate_features": len(cand),
        "matched_weight": matched,
    }


def _overlay(reference, candidate, cid, settings):
    best = None
    reference_view = _heavy_view(reference.mol)
    for rid in reference.relative_energies:
        probe = Chem.Mol(candidate.mol)
        probe_view = _heavy_view(probe)
        shape, color = rdShapeAlign.AlignMol(
            reference_view,
            probe_view,
            refConfId=rid,
            probeConfId=cid,
            useColors=True,
            opt_param=0.5,
            max_preiters=20,
            max_postiters=50,
        )
        pharm, counts = pharmacophore_overlap(
            reference_view, probe_view, rid, cid, settings.pharmacophore_sigma
        )
        if probe_view is not probe:
            # Transfer the rigid alignment to every original atom, including D/T,
            # so exported coordinates retain the actual candidate's identity.
            indices = [a.GetIdx() for a in probe.GetAtoms() if a.GetAtomicNum() != 1]
            rmsd, transform = rdMolAlign.GetAlignmentTransform(
                probe,
                probe_view,
                prbCid=cid,
                refCid=cid,
                atomMap=[(old, new) for new, old in enumerate(indices)],
                reflect=False,
            )
            if rmsd > 1e-3:
                raise ValueError("Could not transfer the heavy-atom alignment to the full molecule")
            rdMolTransforms.TransformConformer(probe.GetConformer(cid), transform)
        if not all(math.isfinite(x) for x in (shape, color, pharm)):
            continue
        joint = math.sqrt(max(0.0, shape * pharm))
        adjusted = joint * math.exp(-settings.strain_penalty * reference.relative_energies[rid])
        info = {
            "shape": float(shape),
            "alignment_color": float(color),
            "pharmacophore": pharm,
            "joint": joint,
            "reference_adjusted_joint": adjusted,
            "reference_conformer_id": rid,
            "reference_relative_energy": reference.relative_energies[rid],
            **counts,
        }
        if best is None or adjusted > best[0]["reference_adjusted_joint"]:
            best = (info, probe)
    if best is None:
        raise ValueError("No finite joint shape/pharmacophore alignment")
    return best


def score_candidate(smiles, hit, design, settings=ConformerSettings(), hit_pose=None, pose_constraints=None):
    """Select one candidate conformer under the versioned strain-adjusted objective.

    Each reference overlay computes shape and pharmacophore on one common pose.
    joint-v1 uses both references. original-hit-v1 selects solely against the
    original hit, then overlays the actual design on that selected conformation.
    Both retain reference and candidate energy penalties and first-on-tie choice.
    Rigid-body orientation can differ. Scores are heuristics, not probabilities.
    Identical reference ensemble objects share an alignment calculation; their
    returned comparison dictionaries and pose molecules remain independent.
    """
    require_alignment()
    if pose_constraints is not None:
        pose_constraints.validate_reference(hit_pose)
    candidate = ensemble(canonical(smiles), settings)
    hit_ref = hit_pose if hit_pose is not None else ensemble(canonical(hit), settings)
    design_ref = ensemble(canonical(design), settings)
    shared_reference = hit_ref is design_ref
    hit_only = settings.scoring_objective == "original-hit-v1"
    best, rejected_poses = None, []
    for cid in candidate.relative_energies:
        hit_info, hit_aligned = _overlay(hit_ref, candidate, cid, settings)
        if shared_reference and not hit_only:
            # Equal SMILES do not imply equal bound/generated reference poses.
            # Copy before pocket assessment so either exported pose can be
            # annotated or modified independently, as with two overlays.
            design_info, design_aligned = dict(hit_info), Chem.Mol(hit_aligned)
        pocket = pose_constraints.evaluate(hit_aligned, cid) if pose_constraints is not None else None
        if pocket is not None and not pocket["passed"]:
            rejected_poses.append({"candidate_conformer_id": cid, **pocket})
            continue
        if not shared_reference and not hit_only:
            design_info, design_aligned = _overlay(design_ref, candidate, cid, settings)
        relative_energy = candidate.relative_energies[cid]
        objective = (
            hit_info["reference_adjusted_joint"]
            if hit_only
            else math.sqrt(hit_info["reference_adjusted_joint"] * design_info["reference_adjusted_joint"])
        )
        score = objective * math.exp(-settings.strain_penalty * relative_energy)
        info = {
            "score": score,
            "scoring_objective": settings.scoring_objective,
            "candidate_smiles": canonical(smiles),
            "hit_smiles": canonical(hit),
            "design_smiles": canonical(design),
            "candidate_conformer_id": cid,
            "candidate_energy": candidate.energies[cid],
            "candidate_relative_energy": relative_energy,
            "force_field": candidate.force_field,
            "hit": hit_info,
            **({"design": design_info} if not hit_only else {}),
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
        }
        if pocket is not None:
            info["pocket"] = pocket
        if best is None or score > best[0]["score"]:
            best = (info, hit_aligned, design_aligned if not hit_only else None)
    if best is None and pose_constraints is not None:
        raise ValueError(
            f"No candidate conformation satisfies bound-pocket constraints ({len(rejected_poses)} rejected)"
        )
    if hit_only and best is not None:
        info, hit_aligned, _ = best
        if shared_reference:
            design_info, design_aligned = dict(info["hit"]), Chem.Mol(hit_aligned)
        else:
            design_info, design_aligned = _overlay(
                design_ref, candidate, info["candidate_conformer_id"], settings
            )
        info["design"] = design_info
        best = (info, hit_aligned, design_aligned)
    if pose_constraints is not None:
        best[0]["pocket_selection"] = {
            "candidate_conformations_evaluated": len(candidate.relative_energies),
            "rejected_conformations": rejected_poses,
            "provenance": pose_constraints.provenance,
        }
    return best
