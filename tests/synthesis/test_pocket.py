import copy

import numpy as np
import pytest
from rdkit import Chem

from chembounce_synthesis.conformers import ConformerSettings, Ensemble, reference_from_sdf, score_candidate
from chembounce_synthesis.pocket import constraints_from_config, validate_pocket_config


def atom_line(
    serial, name, element, residue, number, xyz, *, record="ATOM", chain="A", altloc="", occupancy=1.0
):
    x, y, z = xyz
    return (
        f"{record:<6}{serial:5d} {name:>4}{altloc:1}{residue:>3} {chain:1}{number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{20.0:6.2f}          {element:>2}\n"
    )


def make_pose(smiles="CO", coordinates=((0, 0, 0), (1.4, 0, 0))):
    mol = Chem.MolFromSmiles(smiles)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    conformer.Set3D(True)
    for index, position in enumerate(coordinates):
        conformer.SetAtomPosition(index, position)
    mol.AddConformer(conformer)
    return mol


@pytest.fixture
def pocket_files(tmp_path):
    sdf = tmp_path / "bound.sdf"
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(make_pose())
    pdb = tmp_path / "complex.pdb"
    pdb.write_text(
        atom_line(1, "C1", "C", "LIG", 100, (0, 0, 0), record="HETATM")
        + atom_line(2, "O1", "O", "LIG", 100, (1.4, 0, 0), record="HETATM")
        + atom_line(3, "CA", "C", "ALA", 1, (5.0, 0, 0))
        + atom_line(4, "O", "O", "HOH", 101, (1.4, 0, 0), record="HETATM")
        + atom_line(5, "C1", "C", "COF", 102, (1.4, 0, 0), record="HETATM")
    )
    config = {
        "receptor_pdb": str(pdb),
        "reference_ligand": {"chain": "A", "residue_number": 100, "residue_name": "LIG"},
        "max_clashes": 0,
        "vdw_overlap_tolerance": 0.4,
    }
    return sdf, pdb, reference_from_sdf(sdf, "CO"), config


def atom_anchor(name="oxygen"):
    return {
        "id": name,
        "kind": "atom",
        "reference_atom_index": 1,
        "candidate_smarts": "[O]",
        "candidate_atom_index": 0,
        "max_distance": 0.5,
        "required": True,
    }


def test_receptor_filtering_reference_exclusion_and_hashed_provenance(pocket_files):
    sdf, _, reference, config = pocket_files
    constraints = constraints_from_config(config, reference, sdf)
    assessment = constraints.evaluate(reference.mol, 0)
    assert assessment["passed"]
    assert assessment["clashes"]["pair_count"] == 0
    assert constraints.provenance["receptor_heavy_atoms"] == 1
    assert len(constraints.provenance["reference_sdf_sha256"]) == 64
    assert len(constraints.provenance["receptor_pdb_sha256"]) == 64
    assert constraints.provenance["frame_max_atom_distance"] == pytest.approx(0)
    assert assessment["contacts"]["residues"][0]["residue_name"] == "ALA"
    waters = constraints_from_config({**config, "include_waters": True}, reference, sdf)
    assert not waters.evaluate(reference.mol, 0)["passed"]
    assert waters.provenance["receptor_heavy_atoms"] == 2
    cofactor = constraints_from_config({**config, "include_hetero": True}, reference, sdf)
    assert cofactor.provenance["receptor_heavy_atoms"] == 2  # Excludes reference ligand and water.
    assert all(atom["residue"][3] != "LIG" for atom in cofactor.receptor)


def test_clash_geometry_and_report_cap_preserve_full_pair_count(pocket_files):
    sdf, _, reference, config = pocket_files
    constraints = constraints_from_config(
        {**config, "include_waters": True, "max_reported_pairs": 1}, reference, sdf
    )
    assessment = constraints.evaluate(reference.mol, 0)
    assert assessment["clashes"]["pair_count"] == 2
    assert len(assessment["clashes"]["pairs"]) == 1
    assert assessment["clashes"]["reported_pairs_truncated"]
    assert assessment["clashes"]["maximum_vdw_overlap"] == pytest.approx(3.1)
    assert not assessment["passed"]


def test_bound_coordinate_frame_mismatch_fails_closed(pocket_files):
    sdf, pdb, reference, config = pocket_files
    shifted = pdb.read_text().replace("   0.000   0.000   0.000", "  20.000   0.000   0.000")
    pdb.write_text(shifted)
    with pytest.raises(ValueError, match="coordinate frames"):
        constraints_from_config(config, reference, sdf)
    with pytest.raises(ValueError, match="explicit bound"):
        constraints_from_config(config, None, sdf)
    with pytest.raises(ValueError, match="explicit bound"):
        constraints_from_config(config, reference, None)


def test_in_memory_reference_coordinates_cannot_drift_after_registration(pocket_files):
    sdf, _, reference, config = pocket_files
    constraints = constraints_from_config(config, reference, sdf)
    changed = copy.deepcopy(reference)
    changed.mol.GetConformer().SetAtomPosition(0, (20, 0, 0))
    with pytest.raises(ValueError, match="fixed bound-reference"):
        constraints.validate_reference(changed)
    changed.bound_pose = False
    with pytest.raises(ValueError, match="freely generated"):
        constraints.validate_reference(changed)


def test_alternate_locations_are_coherent_and_explicitly_overridable(pocket_files):
    sdf, pdb, reference, config = pocket_files
    pdb.write_text(
        atom_line(1, "C1", "C", "LIG", 100, (0, 0, 0), record="HETATM", altloc="A", occupancy=0.6)
        + atom_line(2, "O1", "O", "LIG", 100, (1.4, 0, 0), record="HETATM", altloc="A", occupancy=0.6)
        + atom_line(3, "C1", "C", "LIG", 100, (20, 0, 0), record="HETATM", altloc="B", occupancy=0.4)
        + atom_line(4, "O1", "O", "LIG", 100, (21.4, 0, 0), record="HETATM", altloc="B", occupancy=0.4)
        + atom_line(5, "CA", "C", "ALA", 1, (5, 0, 0))
    )
    constraints = constraints_from_config(config, reference, sdf)
    assert constraints.provenance["selected_altlocs"][0]["selected_altloc"] == "A"
    with pytest.raises(ValueError, match="coordinate frames"):
        constraints_from_config({**config, "reference_altloc": "B"}, reference, sdf)
    with pytest.raises(ValueError, match="alternate location is absent"):
        constraints_from_config({**config, "reference_altloc": "C"}, reference, sdf)


def test_atom_anchors_require_unique_matching_and_report_distances(pocket_files):
    sdf, _, reference, _ = pocket_files
    single = constraints_from_config({"anchors": [atom_anchor()]}, reference, sdf)
    result = single.evaluate(reference.mol, 0)
    assert result["passed"]
    assert result["anchors"][0]["distance"] == pytest.approx(0)
    assert result["anchors"][0]["candidate_atom_indices"] == [1]
    conflict = constraints_from_config(
        {"anchors": [atom_anchor("first"), atom_anchor("second")]}, reference, sdf
    )
    result = conflict.evaluate(reference.mol, 0)
    assert not result["passed"]
    assert sum(anchor["matched"] for anchor in result["anchors"]) == 1
    assert all(anchor["nearest_distance"] == pytest.approx(0) for anchor in result["anchors"])
    optional = atom_anchor("optional")
    optional["required"] = False
    result = constraints_from_config({"anchors": [atom_anchor(), optional]}, reference, sdf).evaluate(
        reference.mol, 0
    )
    assert result["passed"]


def test_feature_anchors_match_family_and_cannot_reuse_one_feature(pocket_files):
    sdf, _, reference, _ = pocket_files
    base = constraints_from_config({"anchors": [atom_anchor()]}, reference, sdf)
    feature = next(feature for feature in base.reference_features if feature["family"] == "Acceptor")
    anchor = {
        "id": "acceptor",
        "kind": "feature",
        "reference_feature_index": feature["feature_index"],
        "max_distance": 0.5,
        "required": True,
    }
    constraints = constraints_from_config({"anchors": [anchor]}, reference, sdf)
    assert constraints.evaluate(reference.mol, 0)["passed"]
    donor_only = make_pose("CN")
    assert not constraints.evaluate(donor_only, 0)["passed"]
    duplicate = {**anchor, "id": "second-acceptor"}
    constraints = constraints_from_config({"anchors": [anchor, duplicate]}, reference, sdf)
    assert not constraints.evaluate(reference.mol, 0)["passed"]


def test_candidate_selection_rescues_feasible_conformer_before_ranking(pocket_files, monkeypatch):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    sdf, _, reference, _ = pocket_files
    constraints = constraints_from_config({"anchors": [atom_anchor()]}, reference, sdf)
    mol = make_pose(coordinates=((10, 0, 0), (11.4, 0, 0)))
    mol.AddConformer(Chem.Conformer(reference.mol.GetConformer()), assignId=True)
    candidate = Ensemble(mol, {0: 0.0, 1: 0.0}, {0: 0.0, 1: 0.0}, "fixture", {})
    monkeypatch.setattr("chembounce_synthesis.conformers.ensemble", lambda *args: candidate)

    def overlay(ref, cand, cid, settings):
        value = 0.99 if cid == 0 else 0.8
        return {"shape": value, "pharmacophore": value, "reference_adjusted_joint": value}, Chem.Mol(cand.mol)

    monkeypatch.setattr("chembounce_synthesis.conformers._overlay", overlay)
    assert score_candidate("CO", "CO", "CO", hit_pose=reference)[0]["candidate_conformer_id"] == 0
    result = score_candidate("CO", "CO", "CO", hit_pose=reference, pose_constraints=constraints)
    assert result[0]["candidate_conformer_id"] == 1
    assert result[0]["pocket"]["passed"]
    assert len(result[0]["pocket_selection"]["rejected_conformations"]) == 1
    assert constraints.evaluate(result[1], 1) == result[0]["pocket"]


def test_real_alignment_retains_bound_frame_and_audited_pose(pocket_files):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    sdf, _, reference, _ = pocket_files
    positions = reference.mol.GetConformer().GetPositions().copy()
    constraints = constraints_from_config({"anchors": [atom_anchor()]}, reference, sdf)
    result = score_candidate("CO", "CO", "CO", ConformerSettings(num_conformers=2), reference, constraints)
    assert result[0]["pocket"]["passed"]
    assert np.allclose(reference.mol.GetConformer().GetPositions(), positions)
    assert constraints.evaluate(result[1], result[0]["candidate_conformer_id"]) == result[0]["pocket"]


@pytest.mark.parametrize(
    "config",
    [
        {},
        {"receptor_pdb": "x"},
        {"anchors": [atom_anchor()], "max_clashes": True},
        {"anchors": [atom_anchor()], "coordinate_tolerance": float("nan")},
        {"anchors": [{**atom_anchor(), "candidate_smarts": "[invalid"}]},
    ],
)
def test_invalid_configuration_fails_early(config):
    with pytest.raises(ValueError):
        validate_pocket_config(config)
