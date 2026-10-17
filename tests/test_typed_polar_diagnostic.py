"""Small real geometry and inert control tests; never the 72-condition study."""

import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from rdkit import Chem

from benchmarks.research import typed_polar_diagnostic as d
from chembounce_synthesis import conformers, selective_scoring as scoring
from chembounce_synthesis.common import canonical, molecule
from chembounce_synthesis.landmark_proposals import proper_transform
from chembounce_synthesis.similarity_constraints import BoundSimilarity, pose_binding

FIXTURES = json.loads((Path(__file__).parent / "synthesis/landmark_fixtures.json").read_text())["fixtures"]


def native_api():
    if conformers.rdShapeAlign is None or not hasattr(conformers.rdShapeAlign, "ScoreMol"):
        pytest.skip("Actual fixed-pose API requires modern RDKit; unsupported path tested separately")


def pose(smiles, points):
    mol = molecule(canonical(smiles))
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for i, xyz in enumerate(points):
        conf.SetAtomPosition(i, xyz)
    mol.AddConformer(conf)
    return mol


def ensemble(mol):
    return conformers.Ensemble(mol, {0: 0.0}, {0: 0.0}, "manual_no_energy", {})


@pytest.fixture
def fake_api(monkeypatch):
    api = SimpleNamespace(ShapeInputOptions=lambda: SimpleNamespace(), ScoreMol=lambda *a: (0.8, 0.7))
    monkeypatch.setattr(conformers, "rdShapeAlign", api)
    return api


def feature(kind, i):
    return {"family": kind, "atom_ids": [i], "position": [float(i), 0.0, 0.0]}


def map_row():
    return {"mapping_edge_visits": 0, "helper_edge_visits": 0, "mapping_truncated": False}


def test_typed_injective_subsets_and_bound(fake_api):
    ref = [feature("Acceptor", i) for i in range(3)] + [feature("Donor", 3)]
    cand = [feature("Acceptor", i) for i in range(12)]
    search = d.PolarSearch(clock=lambda: 0.0)
    row = map_row()
    maps = search._maps(ref, cand, [], row)
    assert len(maps) == 4 and len({tuple(m) for m in maps}) == 4
    assert all(len(m) == 3 and len({j for _, j in m}) == 3 for m in maps)
    assert all(ref[i]["family"] == cand[j]["family"] for m in maps for i, j in m)
    assert row["mapping_edge_visits"] <= 128 and row["helper_edge_visits"] <= 4096
    assert search._maps(ref, cand, [], map_row()) == maps


@pytest.mark.parametrize(
    "ref,cand,reason",
    [
        ([feature("Hydrophobe", 0)], [feature("Hydrophobe", 0)], "no_polar_reference"),
        (
            [feature("Acceptor", i) for i in range(3)],
            [feature("Acceptor", 0)],
            "insufficient_typed_polar_capacity",
        ),
        ([feature("Acceptor", 0)], [feature("Donor", 0)], "insufficient_typed_polar_capacity"),
    ],
)
def test_no_polar_or_insufficient_features(fake_api, ref, cand, reason):
    row = map_row()
    assert d.PolarSearch()._maps(ref, cand, [], row) == []
    assert row["mapping_reason"] == reason and row["mapping_edge_visits"] == 0


@pytest.mark.parametrize("target", [[[3, 2, 1], [3, 3, 1], [2, 2, 1]], [[0, 0, 0], [-1, 0, 0], [0, 1, 0]]])
def test_proper_rotation_preserves_distances_not_reflection(target):
    source = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    rotation, shift = proper_transform(source, target)
    assert np.linalg.det(rotation) == pytest.approx(1)
    transformed = source @ rotation.T + shift
    assert np.linalg.norm(transformed[:, None] - transformed[None, :], axis=2) == pytest.approx(
        np.linalg.norm(source[:, None] - source[None, :], axis=2)
    )


@pytest.mark.parametrize("fixture", FIXTURES[:2], ids=lambda f: f["id"])
def test_real_fixed_pose_score_parity_and_input_preservation(fixture):
    native_api()
    ref = pose(fixture["smiles"], fixture["reference_coordinates"])
    cand = pose(fixture["smiles"], fixture["candidate_coordinates"])
    before = pose_binding(ref, 0), pose_binding(cand, 0)
    info, aligned = conformers._overlay(ensemble(ref), ensemble(cand), 0, scoring.DEFAULT_SETTINGS)
    search = d.PolarSearch()
    assert search._score(ref, aligned, 0, 0) == pytest.approx(
        (info["shape"], info["alignment_color"]), abs=1e-5
    )
    assert (pose_binding(ref, 0), pose_binding(cand, 0)) == before


def test_real_typed_polar_rigid_proposal_is_jointly_valid():
    native_api()
    f = FIXTURES[2]
    ref = pose(f["smiles"], f["reference_coordinates"])
    cand = pose(f["smiles"], f["candidate_coordinates"])
    bound = BoundSimilarity(f["smiles"], d.FLOORS)
    # Deliberately shifted real geometry; fixed-pose ScoreMol supplies exact
    # parity input. This tests proposal transforms, not a fabricated baseline.
    search = d.PolarSearch(clock=lambda: 0.0)
    shape, color = search._score(ref, cand, 0, 0)
    pharm, counts = conformers.pharmacophore_overlap(ref, cand, 0, 0)
    info = dict(
        shape=shape,
        alignment_color=color,
        pharmacophore=pharm,
        reference_relative_energy=0.0,
        reference_conformer_id=0,
        **counts,
    )
    before = pose_binding(cand, 0), pose_binding(ref, 0)
    results = search.propose(
        ref, cand, 0, 0, info, bound.assess(ref, cand, 0, 0), bound, scoring.DEFAULT_SETTINGS
    )
    assert results and search.audit["alternatives_attempted"] <= 2
    for info, transformed, witness in results:
        assert "landmark_proposal" not in info
        bound.validate_witness(witness, transformed, 0)
        a, b = (
            np.array(cand.GetConformer().GetPositions()),
            np.array(transformed.GetConformer().GetPositions()),
        )
        assert np.linalg.norm(a[:, None] - a[None, :], axis=2) == pytest.approx(
            np.linalg.norm(b[:, None] - b[None, :], axis=2)
        )
    assert before == (pose_binding(cand, 0), pose_binding(ref, 0))


def test_score_failure_budget_and_late_result(fake_api, monkeypatch):
    f = FIXTURES[2]
    ref = pose(f["smiles"], f["reference_coordinates"])
    clock = [0.0]
    search = d.PolarSearch(clock=lambda: clock[0])

    def late(*args):
        clock[0] = 2.0
        return (0.8, 0.7)

    fake_api.ScoreMol = late
    search.active = clock[0]  # _score normally runs inside propose's active accounting scope.
    assert search._score(ref, ref, 0, 0) is None
    assert search.audit["scoremol_attempted"] == search.audit["scoremol_returned"] == 1
    assert search.audit["truncated"] and search._score(ref, ref, 0, 0) is None
    assert search.audit["scoremol_attempted"] == 1


def test_production_baseline_preserved_with_exact_ensembles(monkeypatch):
    native_api()
    f = FIXTURES[0]
    smiles = canonical(f["smiles"])
    fixed = ensemble(pose(smiles, f["reference_coordinates"]))
    provider = lambda *a: fixed
    monkeypatch.setattr(conformers, "ensemble", provider)
    monkeypatch.setattr(scoring, "ensemble", provider)
    direct = scoring.score_candidate(
        smiles, smiles, smiles, scoring.DEFAULT_SETTINGS, similarity=BoundSimilarity(smiles, d.FLOORS)
    )
    original_overlay, original_scorer = scoring._overlay, scoring.score_candidate
    row = d.paired_condition(smiles, smiles, smiles, 42)
    assert row["status"] == "baseline_feasible_preserved"
    assert row["baseline"] == row["prototype"] and row["extra_search"] is None
    assert row["baseline"]["info"]["selective_pose_evidence"] == direct[0]["selective_pose_evidence"]
    assert row["baseline"]["info"]["score"] == direct[0]["score"]
    assert scoring._overlay is original_overlay and scoring.score_candidate is original_scorer


def test_plan_mismatch_refuses_before_output(tmp_path):
    plan = tmp_path / "plan.json"
    plan.write_text("{}")
    with pytest.raises(ValueError, match="Plan byte mismatch"):
        d.main(["--plan", str(plan), "--plan-sha256", "wrong", "--output", str(tmp_path / "out")])
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("mismatch", [False, True])
def test_failed_baseline_replay_and_design_overlay_scope(monkeypatch, mismatch):
    """Controlled bad baseline pose, genuine metrics/production loop/rigid rescue.

    This is a plumbing control, not a claim the normal optimizer misses this pose.
    The separate baseline test uses the unmodified native optimizer.
    """
    native_api()
    f = FIXTURES[2]
    smiles = canonical(f["smiles"])
    fixed = ensemble(pose(smiles, f["reference_coordinates"]))
    original_overlay = conformers._overlay
    calls = []

    def bad_pose(reference, candidate, cid, settings):
        calls.append(reference)
        if len(calls) == 3:  # design display after experimental hit selection
            return original_overlay(reference, candidate, cid, settings)
        shifted = Chem.Mol(candidate.mol)
        offset = 5.0 + (1.0 if mismatch and len(calls) == 2 else 0.0)
        for i, xyz in enumerate(shifted.GetConformer(cid).GetPositions()):
            shifted.GetConformer(cid).SetAtomPosition(i, (xyz + [offset, 0.0, 0.0]).tolist())
        search = d.PolarSearch(clock=lambda: 0.0)
        shape, color = search._score(reference.mol, shifted, 0, cid)
        pharm, counts = conformers.pharmacophore_overlap(reference.mol, shifted, 0, cid)
        return dict(
            shape=shape,
            alignment_color=color,
            pharmacophore=pharm,
            joint=float(np.sqrt(shape * pharm)),
            reference_adjusted_joint=float(np.sqrt(shape * pharm)),
            reference_conformer_id=0,
            reference_relative_energy=0.0,
            **counts,
        ), shifted

    provider = lambda *a: fixed
    monkeypatch.setattr(conformers, "ensemble", provider)
    monkeypatch.setattr(scoring, "ensemble", provider)
    monkeypatch.setattr(conformers, "_overlay", bad_pose)
    monkeypatch.setattr(scoring, "_overlay", bad_pose)
    before = pose_binding(fixed.mol, 0)
    row = d.paired_condition(smiles, smiles, smiles, 42)
    assert row["baseline"]["status"] == "abstained"
    if mismatch:
        assert row["status"] == "diagnostic_error"
        assert row["prototype"]["status"] == "unknown"
        assert row["extra_search"]["alternatives_attempted"] == 0
    else:
        assert row["status"] == "rescued", row
        assert len(calls) == 3 and row["cost"]["design_overlays"] == 1
        assert len(row["extra_search"]["pairs"]) == 1
        info, audit = row["prototype"]["info"], row["extra_search"]
        assert info["similarity_constraints"]["witness"]["passed"]
        assert d.selected_alternative(info, audit) == row["selected_alternative_index"]
        # A different candidate/reference conformer or pose, score or decision
        # cannot borrow the accepted alternative's provenance.
        for path, replacement in [
            (("selective_pose_evidence", "hit_candidate", "sha256"), "other"),
            (("selective_pose_evidence", "hit_candidate", "conformer_id"), 9),
            (("similarity_constraints", "witness", "reference_pose", "sha256"), "other"),
            (("similarity_constraints", "witness", "reference_pose", "conformer_id"), 9),
            (("similarity_constraints", "witness", "passed"), False),
            (("hit", "reference_conformer_id"), 9),
            (("hit", "shape"), 0.99999),
        ]:
            altered = copy.deepcopy(info)
            node = altered
            for key in path[:-1]:
                node = node[key]
            node[path[-1]] = replacement
            with pytest.raises(d.DiagnosticError, match="not an accepted rigid alternative"):
                d.selected_alternative(altered, audit)
        empty = dict(audit, accepted_alternatives=[])
        with pytest.raises(d.DiagnosticError, match="not an accepted rigid alternative"):
            d.selected_alternative(info, empty)
    assert pose_binding(fixed.mol, 0) == before
    assert scoring._overlay is bad_pose and conformers._overlay is bad_pose


def test_all_72_statuses_survive_errors_without_running_science(tmp_path, monkeypatch):
    import hashlib

    population = [
        dict(attempt=attempt, index=i, product="inert", hit="inert", design="inert")
        for attempt in ("attempt-0000", "attempt-0001")
        for i in range(12)
    ]
    plan = {
        "version": d.VERSION,
        "source_sha256": hashlib.sha256(Path(d.__file__).read_bytes()).hexdigest(),
        "seeds": list(d.SEEDS),
        "limits": d.LIMITS,
        "minimum_scores": d.FLOORS,
        "polar_recall": 0.75,
        "population": population,
    }
    plan.update(
        runtime_versions={
            "python": d.sys.version.split()[0],
            "rdkit": d.rdBase.rdkitVersion,
            "numpy": d.np.__version__,
            "scipy": d.scipy.__version__,
        },
        python_executable_sha256=hashlib.sha256(Path(d.sys.executable).resolve().read_bytes()).hexdigest(),
        production_source_sha256={
            name: hashlib.sha256((Path(d.__file__).resolve().parents[2] / name).read_bytes()).hexdigest()
            for name in d.PRODUCTION_SOURCES
        },
    )
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    calls = []

    def controlled(*args):
        calls.append(args)
        if len(calls) == 2:
            raise d.DiagnosticError("controlled parity failure")
        return {
            "status": "unresolved_allowance_exhausted",
            "baseline": {"status": "abstained"},
            "prototype": {"status": "abstained"},
        }

    monkeypatch.setattr(d, "paired_condition", controlled)
    output = tmp_path / "output"
    args = [
        "--plan",
        str(path),
        "--plan-sha256",
        hashlib.sha256(path.read_bytes()).hexdigest(),
        "--output",
        str(output),
    ]
    assert d.main(args) == 1
    saved = json.loads((output / "summary.json").read_text())
    assert len(calls) == len(saved["rows"]) == 72 and len(saved["paired_products"]) == 24
    assert saved["rows"][1]["status"] == "diagnostic_error"
    assert not any(row["status"] == "not_started" for row in saved["rows"])
    assert not saved["paired_products"][0]["prototype"]["assessment_complete"]
    assert not saved["paired_products"][1]["prototype"]["assessment_complete"]
    assert saved["paired_products"][1]["prototype"]["bounded_calls_returned"] == 3
    assert saved["cost_coverage"]["measured_slots"] == 72
    assert saved["rows"][1]["outer_cost"]["wall_seconds"] >= 0
    assert saved["rows"][1]["outer_cost"]["cpu_seconds"] >= 0
    with pytest.raises(FileExistsError):
        d.main(args)
    assert len(calls) == 72


def test_no_polar_failed_pair_is_not_diagnostic_error(monkeypatch):
    native_api()
    smiles = "c1ccccc1"
    angles = np.arange(6) * np.pi / 3
    points = np.column_stack((1.4 * np.cos(angles), 1.4 * np.sin(angles), np.zeros(6)))
    fixed = ensemble(pose(smiles, points))
    provider = lambda *a: fixed
    monkeypatch.setattr(conformers, "ensemble", provider)
    monkeypatch.setattr(scoring, "ensemble", provider)

    def distant(reference, candidate, cid, settings):
        shifted = Chem.Mol(candidate.mol)
        for i, xyz in enumerate(shifted.GetConformer(cid).GetPositions()):
            shifted.GetConformer(cid).SetAtomPosition(i, (xyz + [6.0, 0.0, 0.0]).tolist())
        search = d.PolarSearch(clock=lambda: 0.0)
        shape, color = search._score(reference.mol, shifted, 0, cid)
        pharm, counts = conformers.pharmacophore_overlap(reference.mol, shifted, 0, cid)
        return dict(
            shape=shape,
            alignment_color=color,
            pharmacophore=pharm,
            joint=float(np.sqrt(shape * pharm)),
            reference_adjusted_joint=float(np.sqrt(shape * pharm)),
            reference_conformer_id=0,
            reference_relative_energy=0.0,
            **counts,
        ), shifted

    monkeypatch.setattr(conformers, "_overlay", distant)
    monkeypatch.setattr(scoring, "_overlay", distant)
    row = d.paired_condition(smiles, smiles, smiles, 42)
    assert row["status"] == "not_rescued", json.dumps(row, indent=2)
    assert row["baseline"]["audit"]["pairs"][0]["polar_recall"] is None
    assert row["extra_search"]["alternatives_attempted"] == 0


def test_baseline_provider_must_be_shared(monkeypatch):
    monkeypatch.setattr(scoring, "ensemble", lambda *a: pytest.fail("must refuse before ensemble work"))
    with pytest.raises(d.DiagnosticError, match="different ensemble provider"):
        d.paired_condition("CCO", "CCO", "CCO", 42)


@pytest.mark.parametrize("field,floor", [("shape", 0.7), ("pharmacophore", 0.5), ("polar", 0.75)])
@pytest.mark.parametrize("exhausted", [False, True])
def test_replay_floor_crossing_is_unknown_even_within_scalar_tolerance(monkeypatch, field, floor, exhausted):
    """Controlled scalar/witness drift, actual production pair loop; no native optimizer.

    Below -> above differs by2e-13, within scalar tolerance but changes a hard
    decision. These deliberately controlled values are not native drift evidence.
    """
    def denied_native_call(*args, **kwargs):
        pytest.fail("Controlled replay boundary must not call native geometry APIs")

    # Satisfy capability checks on minimum RDKit without permitting native work.
    monkeypatch.setattr(
        conformers,
        "rdShapeAlign",
        SimpleNamespace(
            AlignMol=denied_native_call,
            ScoreMol=denied_native_call,
            ShapeInputOptions=denied_native_call,
        ),
    )
    monkeypatch.setattr(conformers.AllChem, "EmbedMultipleConfs", denied_native_call)
    smiles = "c1ccccc1"
    angles = np.arange(6) * np.pi / 3
    points = np.column_stack((1.4 * np.cos(angles), 1.4 * np.sin(angles), np.zeros(6)))
    fixed = ensemble(pose(smiles, points))
    provider = lambda *a: fixed
    monkeypatch.setattr(conformers, "ensemble", provider)
    monkeypatch.setattr(scoring, "ensemble", provider)
    calls = []

    def controlled(reference, candidate, cid, settings):
        calls.append(cid)
        values = dict(shape=1.0, pharmacophore=1.0)
        if field != "polar":
            values[field] = floor + (-1e-13 if len(calls) == 1 else 1e-13)
        joint = float(np.sqrt(values["shape"] * values["pharmacophore"]))
        return dict(
            **values,
            alignment_color=1.0,
            joint=joint,
            reference_adjusted_joint=joint,
            reference_conformer_id=0,
            reference_relative_energy=0.0,
        ), Chem.Mol(candidate.mol)

    monkeypatch.setattr(conformers, "_overlay", controlled)
    monkeypatch.setattr(scoring, "_overlay", controlled)
    if field == "polar":
        actual_assess = BoundSimilarity.assess

        def controlled_assess(self, *args):
            witness = actual_assess(self, *args)
            below = len(calls) == 1
            witness["polar"]["recall"] = floor + (-1e-13 if below else 1e-13)
            witness["reasons"] = ["polar_reference_recall_below_minimum"] if below else []
            witness["passed"] = not below
            return witness

        monkeypatch.setattr(BoundSimilarity, "assess", controlled_assess)
    if exhausted:
        monkeypatch.setattr(d.PolarSearch, "available", lambda self: False)
    row = d.paired_condition(smiles, smiles, smiles, 42)
    assert row["status"] == "diagnostic_error"
    assert row["prototype"]["status"] == "unknown"
    assert "joint/floor decisions changed" in row["prototype"]["error"]
    assert row["extra_search"]["alternatives_attempted"] == 0
    assert row["extra_search"]["accepted_alternatives"] == []
    assert len(calls) == 2


@pytest.mark.parametrize("reference_donors,truncated", [(3, True), (2, False)])
def test_outer_subset_cap_marks_only_unexplored_maps(fake_api, reference_donors, truncated):
    # First subset has exactly4 maps. With3 donors, two more subsets remain;
    # with2 donors and full recall there is only one subset, fully enumerated.
    ref = [feature("Acceptor", i) for i in range(2)] + [
        feature("Donor", i + 2) for i in range(reference_donors)
    ]
    cand = [feature("Acceptor", i) for i in range(2)] + [feature("Donor", i + 2) for i in range(2)]
    search = d.PolarSearch(polar_recall=0.75 if truncated else 1.0)
    row = map_row()
    maps = search._maps(ref, cand, [], row)
    assert len(maps) == 4
    assert row["mapping_truncated"] is truncated
    assert search.audit["truncated"] is truncated
    assert row["helper_edge_visits"] == 1
