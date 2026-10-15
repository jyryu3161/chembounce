import math

import numpy as np
import pytest
from rdkit import Chem

pytest.importorskip(
    "rdkit.Chem.rdShapeAlign", reason="Optional synthesis alignment requires a newer RDKit wheel"
)

from chembounce_synthesis.conformers import (
    ConformerSettings,
    ensemble,
    pharmacophore_overlap,
    reference_from_sdf,
    score_candidate,
)

SETTINGS = ConformerSettings(num_conformers=4)


def test_identity_joint_pose_and_determinism():
    smiles = "CC(=O)Nc1ccncc1"
    one, mol1, _ = score_candidate(smiles, smiles, smiles, SETTINGS)
    two, mol2, _ = score_candidate(smiles, smiles, smiles, SETTINGS)
    assert one == two
    assert one["score"] > 0.99
    assert one["hit"]["shape"] > 0.99 and one["hit"]["pharmacophore"] > 0.99
    cid = one["candidate_conformer_id"]
    assert np.allclose(mol1.GetConformer(cid).GetPositions(), mol2.GetConformer(cid).GetPositions())
    ref = ensemble(smiles, SETTINGS)
    score, _ = pharmacophore_overlap(ref.mol, mol1, one["hit"]["reference_conformer_id"], cid)
    assert score == pytest.approx(one["hit"]["pharmacophore"])


def test_pharmacophore_depends_on_shared_pose_not_graph_identity():
    ens = ensemble("CC(=O)Nc1ccccc1", SETTINGS)
    mol = Chem.Mol(ens.mol)
    cid = next(iter(ens.energies))
    before, _ = pharmacophore_overlap(ens.mol, mol, cid, cid)
    conf = mol.GetConformer(cid)
    for index in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(index, conf.GetAtomPosition(index) + Chem.rdGeometry.Point3D(20, 0, 0))
    after, _ = pharmacophore_overlap(ens.mol, mol, cid, cid)
    assert before > 0.999 and after < 1e-10


def test_energy_window_and_mutation_isolation():
    ens = ensemble("CCCCCCO", ConformerSettings(num_conformers=12, energy_window=1))
    assert all(0 <= e <= 1 for e in ens.relative_energies.values())
    original = ens.mol.GetConformer().GetPositions().copy()
    score_candidate("CCCCCCO", "CCCCCO", "CCCCCCO", ConformerSettings(num_conformers=12, energy_window=1))
    assert np.array_equal(original, ens.mol.GetConformer().GetPositions())


def test_unspecified_stereo_is_an_explicit_failure():
    with pytest.raises(ValueError, match="stereochemistry"):
        ensemble("CC(O)CC", SETTINGS)


def test_bound_pose_is_preserved_and_identity_checked(tmp_path):
    smiles = "CC(=O)Nc1ccccc1"
    ref = ensemble(smiles, SETTINGS)
    path = tmp_path / "pose.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(ref.mol)
    loaded = reference_from_sdf(path, smiles)
    assert loaded.bound_pose
    assert list(loaded.energies.values()) == [None]
    assert np.allclose(
        loaded.mol.GetConformer().GetPositions(), ref.mol.GetConformer().GetPositions(), atol=1e-4
    )
    with pytest.raises(ValueError, match="exact chemical state"):
        reference_from_sdf(path, "CC(=O)Nc1ccncc1")


@pytest.mark.parametrize(
    "kwargs",
    [
        {"num_conformers": 0},
        {"energy_window": math.nan},
        {"pharmacophore_sigma": 0},
        {"seed": -1},
        {"seed": True},
        {"energy_window": True},
        {"pharmacophore_sigma": "1"},
    ],
)
def test_invalid_ensemble_settings(kwargs):
    with pytest.raises(ValueError):
        ConformerSettings(**kwargs)


def test_both_references_must_use_one_candidate_conformation(monkeypatch):
    import chembounce_synthesis.conformers as module

    hit = module.Ensemble(Chem.MolFromSmiles("CCO"), {0: 0.0}, {0: 0.0}, "fixture", {})
    design = module.Ensemble(Chem.MolFromSmiles("CCN"), {0: 0.0}, {0: 0.0}, "fixture", {})
    candidate = module.Ensemble(Chem.MolFromSmiles("CCC"), {0: 0.0, 1: 0.0}, {0: 0.0, 1: 0.0}, "fixture", {})
    monkeypatch.setattr(
        module, "ensemble", lambda s, settings: {"CCO": hit, "CCN": design, "CCC": candidate}[s]
    )

    def overlay(reference, probe, cid, settings):
        # Opposite conformers individually fit each reference perfectly. Taking
        # independent best scores would incorrectly report a joint score of one.
        score = 1.0 if ((reference is hit and cid == 0) or (reference is design and cid == 1)) else 0.25
        return {"reference_adjusted_joint": score}, Chem.Mol(probe.mol)

    monkeypatch.setattr(module, "_overlay", overlay)
    result, _, _ = module.score_candidate("CCC", "CCO", "CCN", SETTINGS)
    assert result["score"] == pytest.approx(0.5)


@pytest.mark.parametrize("isotope", [0, 2])
def test_explicit_hydrogens_do_not_create_features_and_isotopes_survive_export(monkeypatch, isotope):
    import chembounce_synthesis.conformers as module
    from chembounce_synthesis.common import canonical

    smiles = "Cc1csc(NC(=O)c2cccnc2N2CCOCC2)n1"
    reference = ensemble(smiles, ConformerSettings(num_conformers=2))
    full = Chem.AddHs(reference.mol, addCoords=True)
    for neighbor in full.GetAtomWithIdx(0).GetNeighbors():
        if neighbor.GetAtomicNum() == 1:
            neighbor.SetIsotope(isotope)
    identity = canonical(Chem.MolToSmiles(full))
    cid = next(iter(reference.relative_energies))
    overlap, counts = pharmacophore_overlap(reference.mol, full, cid, cid)
    assert counts["candidate_features"] == counts["reference_features"]
    assert overlap == pytest.approx(1.0)
    candidate = module.Ensemble(
        full, reference.energies, reference.relative_energies, reference.force_field, {}
    )
    info, aligned = module._overlay(reference, candidate, cid, SETTINGS)
    assert info["shape"] > 0.999 and info["pharmacophore"] > 0.999
    assert canonical(Chem.MolToSmiles(aligned)) == identity
    assert aligned.GetNumAtoms() == full.GetNumAtoms()
    before = full.GetConformer(cid).GetPositions()
    after = aligned.GetConformer(cid).GetPositions()
    assert (
        np.max(
            np.abs(
                np.linalg.norm(before[:, None] - before[None, :], axis=-1)
                - np.linalg.norm(after[:, None] - after[None, :], axis=-1)
            )
        )
        < 1e-4
    )
