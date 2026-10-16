"""Engineering fixtures, not empirical similarity threshold validation."""

import copy
import hashlib
import math
from dataclasses import replace
from pathlib import Path

import pytest
from rdkit import Chem

from chembounce_synthesis import selective_scoring as scorer
from chembounce_synthesis import similarity_constraints as constraints
from chembounce_synthesis.common import canonical, molecule
from chembounce_synthesis.conformers import Ensemble

HIT = "CC(=O)Nc1ccccc1"
DESIGN = "CC(=O)Nc1ccncc1"


def pose(smiles, count=1, positions=None):
    mol = molecule(canonical(smiles))
    for cid in range(count):
        conf = Chem.Conformer(mol.GetNumAtoms())
        conf.SetId(cid)
        conf.Set3D(True)
        for atom in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(atom, (positions or {}).get(atom, (float(atom), 0.0, 0.0)))
        mol.AddConformer(conf, assignId=False)
    return mol


def request(ident, family, atoms, distance=1.0):
    return {"id": ident, "family": family, "hit_atom_ids": atoms, "max_distance": distance}


def test_joint_assignment_rejects_separately_passing_polar_and_landmark():
    hit = "O=CC=O"
    bound = constraints.BoundSimilarity(hit, landmarks=[request("first", "Acceptor", [0])])
    reference = pose(hit, positions={0: (0, 0, 0), 3: (0.8, 0, 0)})
    candidate = pose(hit, positions={0: (0.8, 0, 0), 3: (-1.01, 0, 0)})
    result = bound.assess(reference, candidate, 0, 0)
    independent_recall = (math.exp(-(1.01**2) / 2) + 1) / 2
    assert independent_recall > 0.75
    assert result["landmarks"][0]["passed"]
    assert result["polar"]["recall"] == pytest.approx(
        (math.exp(-(0.8**2) / 2) + math.exp(-(1.81**2) / 2)) / 2
    )
    assert result["reasons"] == ["polar_reference_recall_below_minimum"]
    assigned = result["joint_assignment"]["assignments"]
    assert len({a["candidate_feature"] for a in assigned}) == len(assigned)


def test_mandatory_matching_solves_greedy_trap_and_rejects_hall_deficiency():
    hit = "O=CC=O"
    bound = constraints.BoundSimilarity(
        hit,
        polar_recall=0.01,
        landmarks=[request("flexible", "Acceptor", [0], 2), request("scarce", "Acceptor", [3], 0.1)],
    )
    reference = pose(hit, positions={0: (0, 0, 0), 3: (0, 0, 0)})
    candidate = pose(hit, positions={0: (0, 0, 0), 3: (1.5, 0, 0)})
    result = bound.assess(reference, candidate, 0, 0)
    assert result["passed"]
    assert [x["assignment"]["candidate_feature"] for x in result["landmarks"]] == [1, 0]
    impossible = constraints.BoundSimilarity(
        hit, landmarks=[request("one", "Acceptor", [0], 0.1), request("two", "Acceptor", [3], 0.1)]
    )
    result = impossible.assess(reference, candidate, 0, 0)
    assert not result["joint_assignment"]["feasible"]
    assert result["polar"]["status"] == "unassessed"
    assert result["reasons"] == ["required_landmark_joint_assignment_infeasible"]


def test_reference_denominator_and_typed_matching_no_reuse():
    bound = constraints.BoundSimilarity("NCCN")
    result = bound.assess(
        pose("NCCN", positions={0: (0, 0, 0), 3: (0, 0, 0)}), pose("CN", positions={1: (0, 0, 0)}), 0, 0
    )
    assert result["polar"]["reference_count"] == 4
    assert result["polar"]["candidate_count"] == 2
    assert result["polar"]["recall"] == 0.5
    for entry in result["joint_assignment"]["assignments"]:
        assert (
            result["reference_features"][entry["reference_feature"]]["family"]
            == result["candidate_features"][entry["candidate_feature"]]["family"]
        )


def test_zero_polar_is_explicit_not_applicable_but_global_empty_rejects():
    bound = constraints.BoundSimilarity("c1ccccc1")
    result = bound.assess(pose(bound.hit), pose(bound.hit), 0, 0)
    assert result["passed"] and result["polar"]["status"] == "not_applicable"
    assert result["polar"]["recall"] is None
    with pytest.raises(ValueError, match="nonempty original-hit features"):
        constraints.BoundSimilarity("[Na+]")
    with pytest.raises(ValueError, match="nonempty candidate features"):
        bound.preflight_candidate("[Na+]")


def test_isotope_hydrogen_does_not_reindex_user_feature_selector():
    bound = constraints.BoundSimilarity("CN[2H]", landmarks=[request("nitrogen", "Donor", [1])])
    assert bound.hit == "[2H]NC"
    result = bound.assess(pose(bound.hit), pose(bound.hit), 0, 0)
    assert result["passed"]
    donor = next(f for f in result["reference_features"] if f["family"] == "Donor")
    assert donor["atom_ids"] == [1]
    assert donor["position"] == [1, 0, 0]
    assert molecule(bound.hit).GetAtomWithIdx(0).GetIsotope() == 2


def test_canonical_input_reordering_invariant_but_wrong_pose_atom_order_fails():
    first = constraints.BoundSimilarity("CCN", landmarks=[request("N", "Donor", [2])])
    second = constraints.BoundSimilarity("NCC", landmarks=[request("N", "Donor", [2])])
    assert first.sha256 == second.sha256
    hit = pose(first.hit)
    wrong = Chem.RenumberAtoms(hit, [2, 1, 0])
    with pytest.raises(ValueError, match="atom order changed"):
        first.assess(wrong, hit, 0, 0)


def test_exact_stereo_and_chiral_core_preflight():
    hit = "C[C@H](O)CCN"
    bound = constraints.BoundSimilarity(hit, core_smiles="C[C@H](O)CC")
    bound.preflight_candidate(hit)
    with pytest.raises(ValueError, match="core/anchor"):
        bound.preflight_candidate("C[C@@H](O)CCN")
    with pytest.raises(ValueError, match="not a chiral substructure"):
        constraints.BoundSimilarity(hit, core_smiles="C[C@@H](O)CC")
    for state in ("CC(O)CCN", "FC=CF"):
        with pytest.raises(ValueError, match="fully specified"):
            constraints.BoundSimilarity(state)


@pytest.mark.parametrize(
    "landmarks",
    [
        [request("missing", "Donor", [999])],
        [request("wrong_family", "Hydrophobe", [0])],
        [request("one", "Donor", [0]), request("two", "Donor", [0])],
        [request("same", "Donor", [0]), request("same", "Donor", [3])],
        [request("bad", "Donor", [True])],
        [request("bad", "Donor", [0, 0])],
        [request("bad", "Donor", [0], float("nan"))],
        [request("bad", [], [0])],
        [request(str(i), "Donor", [0]) for i in range(33)],
    ],
)
def test_invalid_landmarks_rejected_before_any_embedding(landmarks):
    with pytest.raises(ValueError):
        constraints.BoundSimilarity("NCCN", landmarks=landmarks)


def test_feature_bound_and_ambiguous_feature_identity(monkeypatch):
    class Feature:
        def __init__(self, atom):
            self.atom = atom

        def GetFamily(self):
            return "Donor"

        def GetAtomIds(self):
            return (self.atom,)

    class Factory:
        def GetFeaturesForMol(self, mol, confId):
            return [Feature(i) for i in range(mol.GetNumAtoms())]

    monkeypatch.setattr(constraints, "FEATURES", Factory())
    assert len(constraints.feature_records(molecule("C" * 256))) == 256
    with pytest.raises(ValueError, match="feature bound exceeded"):
        constraints.BoundSimilarity("C" * 257)
    monkeypatch.setattr(Factory, "GetFeaturesForMol", lambda self, mol, confId: [Feature(0), Feature(0)])
    with pytest.raises(ValueError, match="Ambiguous duplicate"):
        constraints.BoundSimilarity("CC")


def controlled(monkeypatch, matrix, candidate_energies=None, reference_energies=None, design=DESIGN):
    """Real molecules/constraint matching, controlled scalar overlay metrics."""
    nc = max(c for c, _ in matrix) + 1
    nr = max(r for _, r in matrix) + 1
    ce, re = candidate_energies or [0] * nc, reference_energies or [0] * nr

    def value(smiles, energies, role):
        mol = pose(smiles, len(energies))
        mol.SetProp("role", role)
        return Ensemble(mol, dict(enumerate(energies)), dict(enumerate(energies)), "MMFF94s", {})

    candidate, reference, secondary = (
        value(HIT, ce, "candidate"),
        value(HIT, re, "hit"),
        value(design, [0], "design"),
    )
    pool = iter([candidate, reference, secondary])
    calls = []
    monkeypatch.setattr(scorer, "ensemble", lambda smiles, settings: next(pool))
    monkeypatch.setattr(scorer, "require_alignment", lambda: None)

    def overlay(ref, cand, cid, settings):
        rid = next(iter(ref.relative_energies))
        calls.append((ref.mol.GetProp("role"), cid, rid))
        shape, pharm = (0.2, 0.1) if ref.mol.GetProp("role") == "design" else matrix[cid, rid]
        joint = math.sqrt(shape * pharm)
        info = {
            "shape": shape,
            "pharmacophore": pharm,
            "alignment_color": 0.5,
            "joint": joint,
            "reference_adjusted_joint": joint
            * math.exp(-settings.strain_penalty * ref.relative_energies[rid]),
            "reference_conformer_id": rid,
            "reference_relative_energy": ref.relative_energies[rid],
            "reference_features": 4,
            "candidate_features": 4,
            "matched_weight": 4 * pharm,
        }
        return info, Chem.Mol(cand.mol)

    monkeypatch.setattr(scorer, "_overlay", overlay)
    return candidate, reference, calls


@pytest.mark.parametrize(
    "matrix,selected",
    [
        ({(0, 0): (1, 0.49), (0, 1): (0.70, 0.50)}, (0, 1)),
        ({(0, 0): (0.69, 1), (1, 0): (0.70, 0.50)}, (1, 0)),
    ],
)
def test_rescues_feasible_reference_or_candidate_hidden_by_scalar_best(monkeypatch, matrix, selected):
    _, _, calls = controlled(monkeypatch, matrix)
    info, hit_pose, design_pose = scorer.score_candidate(HIT, HIT, DESIGN)
    assert (info["candidate_conformer_id"], info["hit"]["reference_conformer_id"]) == selected
    assert info["pose_selection"]["hit_pairs_attempted"] == len(matrix)
    assert info["pose_selection"]["hit_pairs_feasible"] == 1
    assert calls[-1] == ("design", selected[0], 0)
    witness = info["similarity_constraints"]["witness"]
    assert witness["candidate_pose_sha256"] == constraints.pose_binding(hit_pose, selected[0])["sha256"]
    assert design_pose is not hit_pose
    assert info["hit"]["shape"] >= 0.70 and info["hit"]["pharmacophore"] >= 0.50


def test_both_energy_penalties_and_first_pair_on_tie(monkeypatch):
    controlled(monkeypatch, {(c, r): (0.9, 0.9) for c in range(2) for r in range(2)}, [2, 0], [2, 0])
    info = scorer.score_candidate(HIT, HIT, DESIGN)[0]
    assert (info["candidate_conformer_id"], info["hit"]["reference_conformer_id"]) == (1, 1)
    assert info["pose_selection"]["pairs"][0]["score"] == pytest.approx(0.9 * math.exp(-0.05 * 2) ** 2)
    controlled(monkeypatch, {(c, r): (0.9, 0.9) for c in range(2) for r in range(2)})
    info = scorer.score_candidate(HIT, HIT, DESIGN)[0]
    assert (info["candidate_conformer_id"], info["hit"]["reference_conformer_id"]) == (0, 0)


def test_all_64_rejections_abstain_without_design_embedding_or_retry(monkeypatch):
    _, _, calls = controlled(monkeypatch, {(c, r): (0.69, 1) for c in range(8) for r in range(8)})
    with pytest.raises(scorer.SimilarityAbstention) as caught:
        scorer.score_candidate(HIT, HIT, DESIGN)
    audit = caught.value.audit
    assert caught.value.reason == "no_jointly_feasible_hit_pose"
    assert audit["hit_pairs_attempted"] == audit["hit_pairs_completed"] == 64
    assert audit["hit_pairs_feasible"] == 0
    assert len(calls) == 64
    assert [e["role"] for e in audit["ensemble_requests"]] == ["candidate", "hit"]
    assert audit["status"] == "abstained"


def test_oversized_ensemble_fails_before_any_alignment(monkeypatch):
    _, _, calls = controlled(monkeypatch, {(c, 0): (0.9, 0.9) for c in range(9)})
    with pytest.raises(scorer.SimilarityAbstention, match="oversized"):
        scorer.score_candidate(HIT, HIT, DESIGN)
    assert calls == []


def test_invalid_settings_and_low_2d_fail_before_embedding(monkeypatch):
    monkeypatch.setattr(scorer, "ensemble", lambda *args: pytest.fail("unexpected embedding"))
    for settings in (
        replace(scorer.DEFAULT_SETTINGS, num_conformers=12),
        replace(scorer.DEFAULT_SETTINGS, reject_unspecified_stereo=False),
        replace(scorer.DEFAULT_SETTINGS, scoring_objective="joint-v1"),
    ):
        with pytest.raises(ValueError, match="Selective scoring requires"):
            scorer.score_candidate(HIT, HIT, DESIGN, settings)
    with pytest.raises(scorer.SimilarityAbstention, match="hit_tanimoto_below_minimum") as caught:
        scorer.score_candidate("CN", HIT, DESIGN)
    assert caught.value.audit["ensemble_requests"] == []


def test_deadline_stops_before_embedding_or_next_pair(monkeypatch):
    _, _, calls = controlled(monkeypatch, {(0, 0): (0.9, 0.9), (0, 1): (0.9, 0.9)})
    with pytest.raises(scorer.SimilarityAbstention, match="deadline"):
        scorer.score_candidate(HIT, HIT, DESIGN, deadline_monotonic=0)
    assert calls == []
    original = scorer._overlay
    clock = [0.0]
    monkeypatch.setattr(scorer.time, "monotonic", lambda: clock[0])

    def slow(*args):
        result = original(*args)
        clock[0] = 10.0
        return result

    monkeypatch.setattr(scorer, "_overlay", slow)
    with pytest.raises(scorer.SimilarityAbstention, match="deadline") as caught:
        scorer.score_candidate(HIT, HIT, DESIGN, deadline_monotonic=5)
    assert len(calls) == 1
    assert caught.value.audit["elapsed_seconds"] == 10


def test_pose_geometry_rejection_is_not_rescued_by_high_global_scores(monkeypatch):
    _, _, _ = controlled(monkeypatch, {(0, 0): (1, 1)})
    original = scorer._overlay

    def shifted(*args):
        info, aligned = original(*args)
        conf = aligned.GetConformer(0)
        for i in range(aligned.GetNumAtoms()):
            xyz = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (xyz.x + 10, xyz.y, xyz.z))
        return info, aligned

    monkeypatch.setattr(scorer, "_overlay", shifted)
    with pytest.raises(scorer.SimilarityAbstention) as caught:
        scorer.score_candidate(HIT, HIT, DESIGN)
    assert caught.value.audit["pairs"][0]["reasons"] == ["polar_reference_recall_below_minimum"]


def test_real_bounded_scoring_preserves_design_independence_and_full_stereo():
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from chembounce_synthesis.qualification import validate_result

    smiles = "C[C@H](O)C(=O)N"
    first = scorer.score_candidate(smiles, smiles, smiles)
    second = scorer.score_candidate(smiles, smiles, "CC(=O)N")
    validate_result(
        first, canonical(smiles), scorer.DEFAULT_SETTINGS, hit=canonical(smiles), design=canonical(smiles)
    )
    validate_result(
        second, canonical(smiles), scorer.DEFAULT_SETTINGS, hit=canonical(smiles), design="CC(N)=O"
    )
    assert first[0]["candidate_conformer_id"] == second[0]["candidate_conformer_id"]
    assert first[0]["hit"] == second[0]["hit"]
    assert first[0]["score"] == second[0]["score"]
    assert first[0]["similarity_constraints"] == second[0]["similarity_constraints"]
    assert canonical(Chem.MolToSmiles(first[1])) == canonical(smiles)
    assert first[0]["pose_selection"]["hit_pairs_attempted"] <= 64
    assert first[0]["similarity_constraints"]["witness"]["passed"]
    scorer.validate_result(
        first,
        canonical(smiles),
        scorer.DEFAULT_SETTINGS,
        hit=canonical(smiles),
        design=canonical(smiles),
        similarity=constraints.BoundSimilarity(smiles),
    )
    scorer.validate_result(
        second,
        canonical(smiles),
        scorer.DEFAULT_SETTINGS,
        hit=canonical(smiles),
        design="CC(N)=O",
        similarity=constraints.BoundSimilarity(smiles),
    )


def test_real_all_feasible_selection_matches_old_scalar_and_isotope_state():
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from chembounce_synthesis.conformers import score_candidate as legacy_score

    smiles = canonical("[2H]NC(C)=O")
    bound = constraints.BoundSimilarity(smiles, landmarks=[request("donor", "Donor", [1])])
    old = legacy_score(smiles, smiles, smiles, scorer.DEFAULT_SETTINGS)
    new = scorer.score_candidate(smiles, smiles, smiles, similarity=bound)
    assert new[0]["score"] == pytest.approx(old[0]["score"], abs=1e-12)
    assert new[0]["hit"] == old[0]["hit"]
    assert canonical(Chem.MolToSmiles(new[1])) == smiles
    scorer.validate_result(new, smiles, scorer.DEFAULT_SETTINGS, hit=smiles, design=smiles, similarity=bound)


@pytest.mark.parametrize(
    "mutation",
    ["recall", "reference_coordinate", "candidate_coordinate", "policy", "reference_id", "pharmacophore"],
)
def test_actual_pose_witness_validator_rejects_changed_bindings(mutation):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    smiles = canonical("CC(=O)N")
    result = scorer.score_candidate(smiles, smiles, smiles)
    info, aligned, _ = result
    witness = info["similarity_constraints"]["witness"]
    if mutation == "recall":
        witness["polar"]["recall"] = 0.99
    elif mutation == "reference_coordinate":
        witness["reference_pose"]["coordinates"][0][0] += 0.1
    elif mutation == "candidate_coordinate":
        conf = aligned.GetConformer(info["candidate_conformer_id"])
        xyz = conf.GetAtomPosition(0)
        conf.SetAtomPosition(0, (xyz.x + 0.1, xyz.y, xyz.z))
    elif mutation == "policy":
        info["similarity_constraints"]["policy"]["minimum_polar_reference_recall"] = 0.1
    elif mutation == "reference_id":
        witness["reference_pose"]["conformer_id"] += 1
    else:
        info["hit"]["pharmacophore"] -= 0.01
    with pytest.raises(ValueError):
        scorer.validate_result(
            result,
            smiles,
            scorer.DEFAULT_SETTINGS,
            hit=smiles,
            design=smiles,
            similarity=constraints.BoundSimilarity(smiles),
        )


def test_constraint_witness_can_be_revalidated_without_embedding():
    bound = constraints.BoundSimilarity("O=CC=O")
    reference, candidate = pose(bound.hit), pose(bound.hit)
    witness = bound.assess(reference, candidate, 0, 0)
    restored = bound.validate_witness(witness, candidate, 0)
    assert (
        restored.GetConformer(0).GetPositions().tolist() == reference.GetConformer(0).GetPositions().tolist()
    )
    changed = copy.deepcopy(witness)
    changed["joint_assignment"]["assignments"][0]["candidate_feature"] = 999
    with pytest.raises(ValueError, match="witness mismatch"):
        bound.validate_witness(changed, candidate, 0)


def test_v1_kernels_and_released_profile_bytes_are_unchanged():
    root = Path(__file__).resolve().parents[2] / "chembounce_synthesis"
    expected = {
        "conformers.py": "aba699db1bc04db517c4275862c08ed73e9ae20c753a72ffa042c34bd907add8",
        "common.py": "f8fabab242db6e463ff8088676c7cedd59cfe3883ad37d1d4571d9a4c3b24d07",
        "identity.py": "a07e99025182a09dbf66a9d2ac98635f8c9eba9895e127b9d25a40b2466423da",
        "data/qualified-development-v1.json": "9d3c60a895912b28c13a29aac0a04f385ddb8a0effc9fcd81d64a50a3d48bdd2",
    }
    for name, digest in expected.items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == digest
