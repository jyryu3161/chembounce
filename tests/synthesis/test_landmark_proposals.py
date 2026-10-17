"""Frozen manual geometry and controlled plumbing, not empirical efficacy."""

import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from rdkit import Chem

from chembounce_synthesis import conformers, landmark_proposals as lp, profiles, selective_scoring as scorer
from chembounce_synthesis.common import canonical, molecule
from chembounce_synthesis.refinement import RefinementSettings, refine_candidate
from chembounce_synthesis.similarity_constraints import BoundSimilarity, pose_binding, restore_pose

FIXTURES = json.loads(Path(__file__).with_name("landmark_fixtures.json").read_text())["fixtures"]
SETTINGS = scorer.DEFAULT_SETTINGS


def pose(smiles, points):
    mol = molecule(canonical(smiles))
    assert len(points) == mol.GetNumAtoms()
    c = Chem.Conformer(len(points))
    c.Set3D(True)
    c.SetId(0)
    for i, p in enumerate(points):
        c.SetAtomPosition(i, p)
    mol.AddConformer(c, assignId=False)
    return mol


def frozen(f):
    return pose(f["smiles"], f["reference_coordinates"]), pose(f["smiles"], f["candidate_coordinates"])


def ensemble(mol):
    return conformers.Ensemble(mol, {0: 0.0}, {0: 0.0}, "manual_no_energy", {})


def native_api():
    if conformers.rdShapeAlign is None or not hasattr(conformers.rdShapeAlign, "ScoreMol"):
        pytest.skip("Fixed-pose ScoreMol unavailable; explicit preflight refusal tested separately")


@pytest.fixture
def fake_api(monkeypatch):
    api = SimpleNamespace(ShapeInputOptions=lambda: SimpleNamespace(), ScoreMol=lambda *a: (0.8, 0.7))
    monkeypatch.setattr(conformers, "rdShapeAlign", api)
    return api


def bound(f, enabled=True):
    return BoundSimilarity(f["smiles"], landmarks=f["landmarks"], landmark_pose_proposals=enabled)


@pytest.mark.parametrize("fixture", FIXTURES, ids=lambda f: f["id"])
def test_known_manual_native_pair(fixture):
    native_api()
    reference, candidate = frozen(fixture)
    before = pose_binding(reference, 0), pose_binding(candidate, 0)
    b = bound(fixture, bool(fixture["landmarks"]))
    info, aligned = conformers._overlay(ensemble(reference), ensemble(candidate), 0, SETTINGS)
    witness = b.assess(reference, aligned, 0, 0)
    if not fixture["landmarks"]:
        assert witness["passed"] and info["shape"] >= 0.7 and info["pharmacophore"] >= 0.5
        assert "landmark_pose_proposals" not in b.policy
    else:
        assert not witness["passed"]
        search = lp.ProposalSearch(b.policy["landmark_pose_proposals"])
        results = search.propose(reference, aligned, 0, 0, info, witness, b, SETTINGS)
        assert search.audit["scoremol_attempted"] <= 3 and search.audit["alternatives_attempted"] <= 2
        assert any(w["passed"] and i["shape"] >= 0.7 and i["pharmacophore"] >= 0.5 for i, _, w in results)
        for _, p, w in results:
            assert pose_binding(p, 0)["sha256"] == w["candidate_pose_sha256"]
            assert pose_binding(reference, 0)["sha256"] == w["reference_pose"]["sha256"]
    assert before == (pose_binding(reference, 0), pose_binding(candidate, 0))


@pytest.mark.parametrize("isotope,n", [(2, 8), (3, 6)])
def test_native_isotope_heavy_metric(isotope, n):
    native_api()
    f = FIXTURES[2 if n == 8 else 3]
    smiles = f"[{isotope}H]N" + "C" * (n - 1)
    points = list(reversed(f["reference_coordinates"]))
    points = [[points[0][0], points[0][1], points[0][2] + 1.0]] + points
    reference = pose(smiles, points)
    candidate = pose(smiles, [[x + 4, -y + 2, -z + 1] for x, y, z in points])
    info, aligned = conformers._overlay(ensemble(reference), ensemble(candidate), 0, SETTINGS)
    b = BoundSimilarity(
        smiles,
        landmarks=[{"id": "N", "family": "Donor", "hit_atom_ids": [1], "max_distance": 0.25}],
        landmark_pose_proposals=True,
    )
    search = lp.ProposalSearch(b.policy["landmark_pose_proposals"])
    score = search._score(reference, aligned, 0, 0)
    assert score == pytest.approx((info["shape"], info["alignment_color"]), abs=1e-5)
    assert pose_binding(aligned, 0)["smiles"] == canonical(smiles)
    before = np.asarray(candidate.GetConformer().GetPositions())
    after = np.asarray(aligned.GetConformer().GetPositions())
    assert np.linalg.norm(before[:, None] - before[None, :], axis=2) == pytest.approx(
        np.linalg.norm(after[:, None] - after[None, :], axis=2), abs=1e-4
    )


@pytest.mark.parametrize("value", [None, 0, 1, "true", [], {}])
def test_strict_bool(value):
    with pytest.raises(ValueError, match="true or false"):
        BoundSimilarity("CCCCCCCN", landmark_pose_proposals=value)


def test_no_landmarks_rejected_and_disabled_no_api(monkeypatch):
    monkeypatch.setattr(conformers, "rdShapeAlign", None)
    a = BoundSimilarity("CCCCCCCN")
    b = BoundSimilarity("CCCCCCCN", landmark_pose_proposals=False)
    assert a.policy == b.policy and a.sha256 == b.sha256
    with pytest.raises(ValueError, match="nonempty landmarks"):
        BoundSimilarity("CCCCCCCN", landmark_pose_proposals=True)
    with pytest.raises(ValueError, match="ScoreMol"):
        bound(FIXTURES[2])


@pytest.mark.parametrize(
    "target",
    [[[0, 0, 0], [1, 0, 0]], [[0, 0, 0], [-1, 0, 0]], [[0, 0, 0], [0, 1, 0]], [[2, 3, 4], [2, 3, 4]]],
)
def test_proper_deterministic_degenerate_rotation(target):
    a = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
    r, t = lp.proper_transform(a, target)
    r2, t2 = lp.proper_transform(a, target)
    assert np.linalg.det(r) == pytest.approx(1)
    assert r.T @ r == pytest.approx(np.eye(3))
    assert np.array_equal(r, r2) and np.array_equal(t, t2)
    assert (a @ r.T + t)[0] == pytest.approx(target[0])


def test_injective_complete_mapping_and_cap(fake_api):
    search = lp.ProposalSearch(lp.declaration())
    ref = [{"family": "Acceptor", "atom_ids": [i], "position": [i, 0, 0]} for i in range(3)]
    requests = [{"family": "Acceptor", "hit_atom_ids": [i]} for i in range(3)]
    row = {"mapping_edge_visits": 0, "mapping_truncated": False}
    assert search._maps(ref, ref[:2], requests, row) == []
    row = {"mapping_edge_visits": 0, "mapping_truncated": False}
    maps = search._maps(ref, ref * 10, requests, row)
    assert len(maps) <= 4 and row["mapping_edge_visits"] <= 128 and row["mapping_truncated"]
    assert all(len({j for _, j in m}) == 3 for m in maps)


class Clock:
    value = 0.0

    def __call__(self):
        return self.value


def baseline_fixture(fake_api):
    f = FIXTURES[2]
    reference, candidate = frozen(f)
    b = bound(f)
    witness = b.assess(reference, candidate, 0, 0)
    info = {
        "shape": 0.8,
        "alignment_color": 0.7,
        "pharmacophore": 0.9,
        "reference_relative_energy": 0.0,
        "reference_conformer_id": 0,
    }
    return b, reference, candidate, witness, info


def test_parity_failure_no_proposals(fake_api):
    b, r, c, w, i = baseline_fixture(fake_api)
    i["shape"] = 0.9
    search = lp.ProposalSearch(b.policy["landmark_pose_proposals"])
    assert search.propose(r, c, 0, 0, i, w, b, SETTINGS) == []
    assert search.audit["pairs"][0]["status"] == "fixed_score_parity_unavailable"
    assert search.audit["scoremol_attempted"] == 1 and search.audit["alternatives_attempted"] == 0


@pytest.mark.parametrize("outer", [False, True])
def test_budget_expiry_before_native(fake_api, outer):
    b, r, c, w, i = baseline_fixture(fake_api)
    clock = Clock()
    search = lp.ProposalSearch(
        b.policy["landmark_pose_proposals"], clock=clock, deadline=1 if outer else None
    )
    if outer:
        clock.value = 2
    else:
        search.charge(1)
    assert search.propose(r, c, 0, 0, i, w, b, SETTINGS) == []
    assert search.audit["scoremol_attempted"] == 0 and search.audit["truncated"]


def test_inflight_overrun_discarded(fake_api):
    b, r, c, w, i = baseline_fixture(fake_api)
    clock = Clock()
    calls = []

    def score(*a):
        calls.append(1)
        if len(calls) == 2:
            clock.value = 1.1
        return 0.8, 0.7

    fake_api.ScoreMol = score
    search = lp.ProposalSearch(b.policy["landmark_pose_proposals"], clock=clock)
    assert search.propose(r, c, 0, 0, i, w, b, SETTINGS) == []
    assert (
        len(calls) == 2 and search.audit["cooperative_overrun"] and search.audit["alternatives_returned"] == 0
    )


@pytest.mark.parametrize("error", [KeyboardInterrupt, SystemExit])
def test_interruption_propagates_with_attempt_accounting(fake_api, error):
    b, r, c, w, i = baseline_fixture(fake_api)

    def score(*a):
        raise error()

    fake_api.ScoreMol = score
    search = lp.ProposalSearch(b.policy["landmark_pose_proposals"])
    with pytest.raises(error):
        search.propose(r, c, 0, 0, i, w, b, SETTINGS)
    assert search.audit["scoremol_attempted"] == 1 and search.audit["scoremol_returned"] == 0


def test_scoring_mutation_cannot_modify_caller(fake_api):
    b, r, c, w, i = baseline_fixture(fake_api)
    before = pose_binding(r, 0), pose_binding(c, 0)

    def score(ref, cand, *a):
        cand.GetConformer().SetAtomPosition(0, (100, 0, 0))
        return 0.8, 0.7

    fake_api.ScoreMol = score
    search = lp.ProposalSearch(b.policy["landmark_pose_proposals"])
    assert search.propose(r, c, 0, 0, i, w, b, SETTINGS) == []
    assert before == (pose_binding(r, 0), pose_binding(c, 0))
    assert search.audit["pairs"][0]["status"] == "unavailable"


def test_caps_and_full_pair_accounting(fake_api):
    b, r, c, w, i = baseline_fixture(fake_api)
    search = lp.ProposalSearch(b.policy["landmark_pose_proposals"], clock=lambda: 0)
    for _ in range(12):
        search.propose(r, c, 0, 0, i, w, b, SETTINGS)
    assert search.audit["alternatives_attempted"] <= 8 and search.audit["scoremol_attempted"] <= 16
    assert len(search.audit["pairs"]) == 12 and search.audit["truncated"]


def controlled_scorer(monkeypatch, fake_api):
    f = FIXTURES[2]
    reference, _ = frozen(f)
    b = bound(f)
    calls = []
    monkeypatch.setattr(scorer, "require_alignment", lambda: None)
    monkeypatch.setattr(scorer, "ensemble", lambda *a: ensemble(Chem.Mol(reference)))

    def overlay(ref, cand, cid, settings):
        calls.append((cid, next(iter(ref.relative_energies))))
        aligned = Chem.Mol(cand.mol)
        if len(calls) % 2:
            for a, p in enumerate(aligned.GetConformer(cid).GetPositions()):
                aligned.GetConformer(cid).SetAtomPosition(a, (p + [0, 2, 0]).tolist())
        _, counts = conformers.pharmacophore_overlap(reference, reference, 0, 0, 1)
        return {
            **counts,
            "shape": 0.95,
            "alignment_color": 0.7,
            "pharmacophore": 1.0,
            "joint": np.sqrt(0.95),
            "reference_adjusted_joint": np.sqrt(0.95),
            "reference_conformer_id": 0,
            "reference_relative_energy": 0.0,
        }, aligned

    monkeypatch.setattr(scorer, "_overlay", overlay)
    fake_api.ScoreMol = lambda ref, cand, *a: (
        (0.95, 0.7) if cand.GetConformer().GetAtomPosition(0).y > 1 else (0.75, 0.7)
    )
    return b, reference, calls


def test_feasible_lower_scalar_wins_and_origin_validates(monkeypatch, fake_api):
    b, reference, calls = controlled_scorer(monkeypatch, fake_api)
    result = scorer.score_candidate(b.hit, b.hit, b.hit, similarity=b)
    info, aligned, _ = result
    assert info["score"] < np.sqrt(0.95) and info["hit"]["shape"] == 0.75
    assert (
        info["pose_selection"]["hit_pairs_completed"] == 1
        and info["pose_selection"]["hit_pairs_feasible"] == 1
    )
    assert info["pose_selection"]["pairs"][0]["baseline"]["status"] == "rejected"
    assert len(calls) == 2
    lp.validate_selected(info, aligned, b)
    pair = info["pose_selection"]["pairs"][0]
    assert pair["reasons"] == [] and pair["score"] == info["score"]
    assert (
        pair["shape"] == info["hit"]["shape"]
        and pair["polar_recall"] == info["similarity_constraints"]["witness"]["polar"]["recall"]
    )
    assert pair["selected_pose_sha256"] == pose_binding(aligned, 0)["sha256"]
    for key in ("candidate_pose_sha256", "reference_pose_sha256", "proposal_policy_sha256"):
        changed = copy.deepcopy(info)
        changed["hit"]["landmark_proposal"][key] = "0" * 64
        with pytest.raises(ValueError):
            lp.validate_selected(changed, aligned, b)
    changed = copy.deepcopy(info)
    changed["hit"]["landmark_proposal"]["transform"]["rotation"][0][0] *= -1
    with pytest.raises(ValueError):
        lp.validate_selected(changed, aligned, b)


def test_profile_artifacts_not_rebound():
    import hashlib

    for name, digest in profiles.PROFILE_SHA256.items():
        assert hashlib.sha256(profiles._resource(name).read_bytes()).hexdigest() == digest
    config = profiles.apply_profile_defaults(
        {"recommendation_profile": profiles.DEFAULT_RECOMMENDATION_PROFILE}
    )
    assert profiles.profile_provenance(config)["status"] == "custom_unvalidated"


def test_refinement_retains_each_seed_origin(monkeypatch, fake_api):
    b, _, calls = controlled_scorer(monkeypatch, fake_api)
    result, report = refine_candidate(
        b.hit,
        b.hit,
        b.hit,
        SETTINGS,
        RefinementSettings(top_n=1, num_conformers=8, seeds=(42, 43, 44), max_total_conformers=24),
        similarity=b,
    )
    assert result is not None and report["status"] == "success"
    assert len(calls) == 6
    assert [t["seed"] for t in report["trials"]] == [42, 43, 44]
    for trial in report["trials"]:
        restored = scorer.restore_result(trial)
        seeded = conformers.ConformerSettings(**trial["settings"])
        scorer.validate_result(restored, b.hit, seeded, hit=b.hit, design=b.hit, similarity=b)
        assert restore_pose(trial["selective_pose_evidence"]["hit_candidate"]).GetNumConformers() == 1


def test_config_early_api_and_same_bound_policy(fake_api):
    config = profiles.apply_profile_defaults(
        {"recommendation_profile": profiles.DEFAULT_RECOMMENDATION_PROFILE}
    )
    config["similarity_constraints"].update(landmarks=FIXTURES[2]["landmarks"], landmark_pose_proposals=True)
    b = scorer.configured_constraints(config, FIXTURES[2]["smiles"])
    assert b.policy == bound(FIXTURES[2]).policy
    assert profiles.profile_provenance(config)["status"] == "custom_unvalidated"
    fake_api.ScoreMol = None
    with pytest.raises(ValueError, match="ScoreMol"):
        scorer.configured_constraints(config)


def test_disabled_full_result_and_call_parity(monkeypatch, fake_api):
    from .test_selective_scoring import controlled, HIT, DESIGN

    _, _, calls = controlled(monkeypatch, {(0, 0): (0.9, 0.9)})
    a = scorer.score_candidate(HIT, HIT, DESIGN)
    _, _, second_calls = controlled(monkeypatch, {(0, 0): (0.9, 0.9)})
    b = scorer.score_candidate(
        HIT, HIT, DESIGN, similarity=BoundSimilarity(HIT, landmark_pose_proposals=False)
    )
    assert calls == second_calls and len(calls) == 2
    assert a[0]["selective_pose_evidence"] == b[0]["selective_pose_evidence"]
    assert a[0]["score"] == b[0]["score"] and a[0]["similarity_constraints"] == b[0]["similarity_constraints"]
    assert a[0]["pose_selection"]["pairs"] == b[0]["pose_selection"]["pairs"]
    assert "landmark_pose_proposals" not in a[0]["pose_selection"]


@pytest.mark.parametrize(
    "field", ["score", "shape", "pharmacophore", "polar_recall", "reasons", "selected_origin"]
)
def test_selected_pair_contradiction_rejected(monkeypatch, fake_api, field):
    b, _, _ = controlled_scorer(monkeypatch, fake_api)
    info, aligned, _ = scorer.score_candidate(b.hit, b.hit, b.hit, similarity=b)
    changed = copy.deepcopy(info)
    changed["pose_selection"]["pairs"][0][field] = (
        ["contradiction"] if field == "reasons" else "baseline" if field == "selected_origin" else -0.1
    )
    with pytest.raises(ValueError):
        lp.validate_selected(changed, aligned, b)


def test_proposal_rejection_saved_in_screen_and_refinement(monkeypatch, fake_api):
    from chembounce_synthesis.rejected_pose_diagnostics import RejectedPoseDiagnostics

    b, _, _ = controlled_scorer(monkeypatch, fake_api)
    # Actual scalar values are controlled; all proposed shape values fail the
    # unchanged .70 floor. A successful landmark assignment is insufficient.
    fake_api.ScoreMol = lambda ref, cand, *a: (
        (0.95, 0.7) if cand.GetConformer().GetAtomPosition(0).y > 1 else (0.6, 0.7)
    )
    collector = RejectedPoseDiagnostics()
    with pytest.raises(scorer.SimilarityAbstention) as caught:
        scorer.score_candidate(
            b.hit,
            b.hit,
            b.hit,
            similarity=b,
            rejected_pose_diagnostics=collector,
            rejected_pose_diagnostic_phase="screening",
        )
    assert caught.value.reason == "no_jointly_feasible_hit_pose"
    assert caught.value.audit["landmark_pose_proposals"]["alternatives_returned"] >= 1
    assert caught.value.audit["pairs"][0]["proposal_outcomes"][0]["reasons"] == ["hit_shape_below_minimum"]
    # Failed hit selection never reaches design diagnostics, so every call in
    # this seed block is a shifted baseline (not alternating hit/design calls).
    b, reference, _ = controlled_scorer(monkeypatch, fake_api)
    prior_overlay = scorer._overlay

    def rejected_overlay(*args):
        info, aligned = prior_overlay(*args)
        for atom, point in enumerate(reference.GetConformer().GetPositions()):
            aligned.GetConformer().SetAtomPosition(atom, (point + [0, 2, 0]).tolist())
        return info, aligned

    monkeypatch.setattr(scorer, "_overlay", rejected_overlay)
    fake_api.ScoreMol = lambda ref, cand, *a: (
        (0.95, 0.7) if cand.GetConformer().GetAtomPosition(0).y > 1 else (0.6, 0.7)
    )
    result, report = refine_candidate(
        b.hit,
        b.hit,
        b.hit,
        SETTINGS,
        RefinementSettings(top_n=1, num_conformers=8, seeds=(42, 43, 44), max_total_conformers=24),
        similarity=b,
        rejected_pose_diagnostics=collector,
    )
    assert result is None and report["status"] != "success"
    assert all(t["similarity_failure"]["scorer_audit"]["landmark_pose_proposals"] for t in report["trials"])
    assert collector.snapshot()["records"][0]["phase"] == "screening"


def test_same_bound_reaches_workflow_screening_config(fake_api):
    from chembounce_synthesis import workflow

    config = profiles.apply_profile_defaults(
        {
            "schema_version": 1,
            "recommendation_profile": profiles.DEFAULT_RECOMMENDATION_PROFILE,
            "synthesis_evidence": {
                "predictor": {"kind": "reactiont5", "model_path": "synthetic-unused-model"}
            },
            "catalog": "synthetic.tsv",
            "rules": "synthetic.json",
            "provider": {"kind": "saved", "path": "synthetic-proposals.json"},
        }
    )
    config["similarity_constraints"].update(landmarks=FIXTURES[2]["landmarks"], landmark_pose_proposals=True)
    resolved = workflow._config(config)
    b = scorer.configured_constraints(resolved, FIXTURES[2]["smiles"])
    assert b.policy == bound(FIXTURES[2]).policy
    assert resolved["similarity_constraints"]["landmark_pose_proposals"] is True
