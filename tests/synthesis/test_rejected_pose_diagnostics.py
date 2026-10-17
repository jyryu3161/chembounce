"""Controlled diagnostic plumbing only; no empirical scoring or model inference."""

import copy
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import pytest

from chembounce_synthesis import profiles, rejected_pose_diagnostics as diagnostics, workflow
from chembounce_synthesis.refinement import RefinementSettings, refine_candidate
from .test_selective_scoring import DESIGN, HIT, constraints, controlled, pose, scorer
from .test_selective_workflow import config_for


class Clock:
    def __init__(self):
        self.now = 0.0

    def __call__(self):
        return self.now


def capture(monkeypatch, collector, matrix=None, *, settings=None, **kwargs):
    _, _, calls = controlled(monkeypatch, matrix or {(0, 0): (0.69, 0.9)})
    with pytest.raises(scorer.SimilarityAbstention) as caught:
        scorer.score_candidate(
            HIT,
            HIT,
            DESIGN,
            settings or scorer.DEFAULT_SETTINGS,
            rejected_pose_diagnostics=collector,
            **kwargs,
        )
    return caught.value, calls


def pending(collector, *, state="C[C@H](O)CCN", phase="direct"):
    bound = constraints.BoundSimilarity(state)
    mol = pose(bound.hit)
    witness = bound.assess(mol, mol, 0, 0)
    pair = {
        "candidate_conformer_id": 0,
        "reference_conformer_id": 0,
        "status": "rejected",
        "reasons": ["hit_shape_below_minimum"],
        "score": 0.4,
        "shape": 0.4,
        "pharmacophore": 0.9,
        "polar_recall": witness["polar"]["recall"],
    }
    call = collector.start_call(
        {"product": bound.hit, "hit": bound.hit, "design": bound.hit},
        asdict(scorer.DEFAULT_SETTINGS),
        bound.policy,
        None,
        phase,
    )
    call.consider(0, pair, mol, witness, {"shape": 0.4, "pharmacophore": 0.9}, 0.0, 0.0)
    return call, mol, pair, witness


def test_rejected_best_first_tie_exact_geometry_energy_and_unchanged_calls(monkeypatch):
    matrix = {(0, 0): (0.5, 0.9), (0, 1): (0.69, 0.9), (1, 0): (0.69, 0.9), (1, 1): (0.6, 0.9)}
    old, old_calls = capture(monkeypatch, None, matrix)
    collector = diagnostics.RejectedPoseDiagnostics()
    new, new_calls = capture(monkeypatch, collector, matrix, rejected_pose_diagnostic_phase="screening")
    assert old.reason == new.reason == "no_jointly_feasible_hit_pose"
    assert old_calls == new_calls and len(new_calls) == 4
    assert old.audit["pairs"] == new.audit["pairs"]
    record = collector.snapshot()["records"][0]
    assert record["pair_index"] == 1
    assert record["phase"] == "screening" and record["settings"]["seed"] == 42
    assert record["pair"]["reference_conformer_id"] == 1
    assert record["candidate_pose"]["sha256"] == record["witness"]["candidate_pose_sha256"]
    assert record["witness"]["reference_pose"]["conformer_id"] == 1
    assert record["policy_sha256"] == diagnostics.digest(record["policy"])
    assert record["global_pharmacophore_assignment"]["status"] == "not_retained"
    assert new.audit["rejected_pose_diagnostic"]["status"] == "captured_not_published"
    # No additional global feature assignment is invented from counts.
    assert "assignments" not in record["global_pharmacophore_assignment"]


def test_selection_uses_both_existing_energy_penalties(monkeypatch):
    collector = diagnostics.RejectedPoseDiagnostics()
    controlled(monkeypatch, {(c, r): (0.69, 0.9) for c in range(2) for r in range(2)}, [2, 0], [2, 0])
    with pytest.raises(scorer.SimilarityAbstention):
        scorer.score_candidate(HIT, HIT, DESIGN, rejected_pose_diagnostics=collector)
    r = collector.snapshot()["records"][0]
    assert (r["pair"]["candidate_conformer_id"], r["pair"]["reference_conformer_id"]) == (1, 1)
    assert r["candidate_relative_energy"] == r["reference_relative_energy"] == 0


def test_success_has_identical_result_and_no_diagnostic_promotion(monkeypatch):
    controlled(monkeypatch, {(0, 0): (0.9, 0.9)})
    original = scorer.score_candidate(HIT, HIT, DESIGN)
    collector = diagnostics.RejectedPoseDiagnostics()
    controlled(monkeypatch, {(0, 0): (0.9, 0.9)})
    result = scorer.score_candidate(HIT, HIT, DESIGN, rejected_pose_diagnostics=collector)
    assert result[0]["score"] == original[0]["score"]
    assert result[0]["selective_pose_evidence"] == original[0]["selective_pose_evidence"]
    assert collector.snapshot()["records"] == []
    assert collector.statistics()["not_applicable_calls"] == 1
    assert not collector.statistics()["truncated"]


@pytest.mark.parametrize(
    "reason",
    [
        "deadline_exhausted",
        "hit_tanimoto_below_minimum",
        "selected_conformer_design_diagnostic_failed",
        "selected",
    ],
)
def test_only_final_no_joint_decision_serializes(reason, monkeypatch):
    collector = diagnostics.RejectedPoseDiagnostics()
    call, *_ = pending(collector)
    monkeypatch.setattr(constraints, "pose_binding", lambda *a: pytest.fail("unexpected serialization"))
    assert call.finish(reason)["status"] == "not_applicable"
    assert not collector.snapshot()["records"]


def test_witness_snapshot_is_immutable_but_changed_candidate_pose_rejects():
    collector = diagnostics.RejectedPoseDiagnostics()
    call, mol, pair, witness = pending(collector)
    original = copy.deepcopy(witness)
    witness["polar"]["recall"] = -100
    pair["score"] = 100
    ref = call.finish("no_jointly_feasible_hit_pose")
    assert ref["status"] == "captured_not_published"
    saved = collector.snapshot()["records"][0]
    assert saved["witness"] == original and saved["pair"]["score"] == 0.4
    assert "@" in saved["candidate_pose"]["smiles"]
    saved["witness"]["polar"]["recall"] = -200
    assert collector.snapshot()["records"][0]["witness"] == original
    call, mol, *_ = pending(collector)
    mol.GetConformer().SetAtomPosition(0, (100, 0, 0))
    assert call.finish("no_jointly_feasible_hit_pose")["reason"] == "serialization_error:ValueError"
    assert collector.statistics()["captures"] == 1


def test_infeasible_joint_witness_preserved_not_fabricated():
    collector = diagnostics.RejectedPoseDiagnostics()
    bound = constraints.BoundSimilarity(
        "O=CC=O",
        landmarks=[
            {"id": "a", "family": "Acceptor", "hit_atom_ids": [0], "max_distance": 0.1},
            {"id": "b", "family": "Acceptor", "hit_atom_ids": [3], "max_distance": 0.1},
        ],
    )
    ref = pose(bound.hit, positions={0: (0, 0, 0), 3: (0, 0, 0)})
    mol = pose(bound.hit, positions={0: (0, 0, 0), 3: (1.5, 0, 0)})
    witness = bound.assess(ref, mol, 0, 0)
    call = collector.start_call(
        {"product": bound.hit, "hit": bound.hit, "design": bound.hit},
        asdict(scorer.DEFAULT_SETTINGS),
        bound.policy,
        None,
        "direct",
    )
    call.consider(
        0,
        {
            "candidate_conformer_id": 0,
            "reference_conformer_id": 0,
            "score": 0.8,
            "status": "rejected",
            "reasons": witness["reasons"],
        },
        mol,
        witness,
        {},
        0,
        0,
    )
    assert call.finish("no_jointly_feasible_hit_pose")["status"] == "captured_not_published"
    saved = collector.snapshot()["records"][0]["witness"]
    assert not saved["joint_assignment"]["feasible"]
    assert saved["polar"]["status"] == "unassessed"
    assert saved["joint_assignment"]["assignments"] == []


def test_shared_eight_capture_limit_and_one_per_call():
    collector = diagnostics.RejectedPoseDiagnostics(clock=lambda: 0.0)
    for index in range(10):
        call, *_ = pending(collector)
        ref = call.finish("no_jointly_feasible_hit_pose")
        assert ref["status"] == ("captured_not_published" if index < 8 else "unavailable")
        if index == 0:
            assert call.finish("no_jointly_feasible_hit_pose")["reason"] == "call_already_finished"
    assert collector.statistics()["captures"] == 8
    assert collector.statistics()["skips"]["record_count_exhausted"] == 2


@pytest.mark.parametrize(
    "bound,expected",
    [
        ("atoms", "atom_limit"),
        ("record_bytes", "record_bytes_exhausted"),
        ("run_bytes", "run_bytes_exhausted"),
    ],
)
def test_size_limits_no_truncated_witness(monkeypatch, bound, expected):
    collector = diagnostics.RejectedPoseDiagnostics()
    monkeypatch.setitem(
        diagnostics.LIMITS, bound, 1 if bound != "run_bytes" else diagnostics.HEADER_RESERVE + 1
    )
    call, *_ = pending(collector)
    assert call.finish("no_jointly_feasible_hit_pose")["reason"] == expected
    assert collector.statistics()["captures"] == 0


def test_time_deadline_and_failed_capture_accounting(monkeypatch):
    clock = Clock()
    collector = diagnostics.RejectedPoseDiagnostics(clock=clock)
    call, *_ = pending(collector)
    original = constraints.pose_binding

    def slow(*args):
        clock.now += 1.1
        return original(*args)

    monkeypatch.setattr(constraints, "pose_binding", slow)
    assert call.finish("no_jointly_feasible_hit_pose")["reason"] == "diagnostic_time_exhausted"
    assert collector.statistics()["time_overrun"]
    assert not collector.snapshot()["records"]
    collector = diagnostics.RejectedPoseDiagnostics(clock=Clock())
    call, *_ = pending(collector)
    call.deadline = -1
    assert call.finish("no_jointly_feasible_hit_pose")["reason"] == "workflow_deadline_exhausted"


def test_lost_later_coverage_does_not_publish_lower_scoring_best():
    clock = Clock()
    collector = diagnostics.RejectedPoseDiagnostics(clock=clock)
    call, mol, pair, witness = pending(collector)
    call.deadline = 1
    clock.now = 2
    call.consider(1, {**pair, "score": 0.9}, mol, witness, {}, 0, 0)
    call.deadline = None
    assert call.finish("no_jointly_feasible_hit_pose")["reason"] == "workflow_deadline_exhausted"


def test_full_sidecar_bound_hash_and_conditional_refs(tmp_path):
    collector = diagnostics.RejectedPoseDiagnostics()
    call, *_ = pending(collector)
    ref = call.finish("no_jointly_feasible_hit_pose")
    result = diagnostics.publish(collector, tmp_path)
    data = (tmp_path / diagnostics.ARTIFACT).read_bytes()
    assert result["status"] == "written"
    assert len(data) == result["bytes"] <= diagnostics.LIMITS["run_bytes"]
    assert hashlib.sha256(data).hexdigest() == result["sha256"]
    artifact = json.loads(data)
    assert diagnostics.digest(artifact["records"][ref["record_index"]]) == ref["record_sha256"]
    assert artifact["publication_completion_within_budget"] is None
    assert ref["status"] == "captured_not_published"
    assert artifact["provenance"]["source_sha256"]["selective_scoring.py"]
    assert diagnostics.publish(collector, tmp_path)["status"] == "unavailable"


@pytest.mark.parametrize("kind", ["write", "full_envelope", "before_deadline", "during_serialization"])
def test_optional_persistence_failure_never_implies_availability(tmp_path, monkeypatch, kind):
    clock = Clock()
    collector = diagnostics.RejectedPoseDiagnostics(clock=clock)
    call, *_ = pending(collector)
    ref = call.finish("no_jointly_feasible_hit_pose")
    deadline = None
    if kind == "write":
        monkeypatch.setattr(diagnostics, "_write_bytes", lambda *a: (_ for _ in ()).throw(OSError("disk")))
    elif kind == "full_envelope":
        monkeypatch.setitem(diagnostics.LIMITS, "run_bytes", collector.statistics()["record_bytes"] + 1)
    elif kind == "before_deadline":
        deadline = -1
    else:
        original = collector.snapshot

        def slow():
            value = original()
            clock.now += 1.1
            return value

        monkeypatch.setattr(collector, "snapshot", slow)
    result = diagnostics.publish(collector, tmp_path, deadline_monotonic=deadline)
    assert result["status"] == "unavailable"
    assert not (tmp_path / diagnostics.ARTIFACT).exists()
    assert ref["status"] == "captured_not_published"


def test_late_atomic_commit_is_observed_not_certified(tmp_path, monkeypatch):
    clock = Clock()
    collector = diagnostics.RejectedPoseDiagnostics(clock=clock)
    original = Path.replace

    def late(path, target):
        result = original(path, target)
        clock.now += 2
        return result

    monkeypatch.setattr(Path, "replace", late)
    result = diagnostics.publish(collector, tmp_path)
    assert result["status"] == "written" and result["time_overrun"]
    assert result["publication_completion_within_budget"] is None


@pytest.mark.parametrize("exc", [KeyboardInterrupt, SystemExit])
def test_interruptions_are_never_optional_failure(tmp_path, monkeypatch, exc):
    collector = diagnostics.RejectedPoseDiagnostics()
    monkeypatch.setattr(diagnostics, "_write_bytes", lambda *a: (_ for _ in ()).throw(exc()))
    with pytest.raises(exc):
        diagnostics.publish(collector, tmp_path)


def test_refinement_seed_failures_share_budget_and_preserve_references(monkeypatch):
    collector = diagnostics.RejectedPoseDiagnostics()
    original = scorer.score_candidate

    def controlled_call(*args, **kwargs):
        controlled(monkeypatch, {(0, 0): (0.69, 0.9)})
        return original(*args, **kwargs)

    monkeypatch.setattr(scorer, "score_candidate", controlled_call)
    result, report = refine_candidate(
        HIT,
        HIT,
        DESIGN,
        scorer.DEFAULT_SETTINGS,
        RefinementSettings(top_n=1, num_conformers=8, seeds=[42, 43, 44]),
        similarity=constraints.BoundSimilarity(HIT),
        rejected_pose_diagnostics=collector,
    )
    assert result is None and report["status"] == "failed"
    refs = [r["similarity_failure"]["scorer_audit"]["rejected_pose_diagnostic"] for r in report["trials"]]
    assert [r["record_index"] for r in refs] == [0, 1, 2]
    records = collector.snapshot()["records"]
    assert [r["settings"]["seed"] for r in records] == [42, 43, 44]
    assert {r["phase"] for r in records} == {"refinement"}


@pytest.mark.parametrize("value", [0, 1, None, "false", {}, []])
def test_configuration_boolean_is_strict(value):
    with pytest.raises(ValueError, match="true or false"):
        workflow._config({"schema_version": 1, "rejected_pose_diagnostics": value})


def test_configuration_requires_selective_policy():
    with pytest.raises(ValueError, match="selective"):
        workflow._config({"schema_version": 1, "rejected_pose_diagnostics": True})


@pytest.mark.parametrize("name", [profiles.SELECTIVE_PROFILE, profiles.GUARDED_SELECTIVE_PROFILE])
def test_historical_profiles_remain_honestly_unvalidated(name):
    root = Path(profiles.__file__).parent
    artifact = profiles.load_profile(name)
    config = profiles.apply_profile_defaults({"schema_version": 1, "recommendation_profile": name})
    provenance = profiles.profile_provenance(config)
    assert provenance["status"] == "custom_unvalidated"
    assert "chembounce_synthesis/selective_scoring.py" in provenance["runtime_source_mismatches"]
    assert config["minimum_scores"] == constraints.MINIMUM_SCORES
    path = root / "data" / (name + ".json")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == profiles.PROFILE_SHA256[name]
    assert artifact["status"] == "declared_engineering_policy"


def workflow_control(monkeypatch, *, fail_screening):
    """Manual conformers/controlled shape; existing matching and validation run.

    No embedding, alignment, model or production catalog. Tiny catalog/rule
    fixtures are the existing unit-test inputs, not a feasibility study.
    """
    import math
    from rdkit import Chem
    from chembounce_synthesis.conformers import Ensemble, pharmacophore_overlap

    calls = []
    phase = {"name": "screening"}
    original = scorer.score_candidate

    def ensemble(smiles, settings):
        calls.append(("ensemble", smiles, settings.seed))
        return Ensemble(pose(smiles), {0: 0.0}, {0: 0.0}, "MMFF94s", {})

    def overlay(reference, candidate, cid, settings):
        calls.append(("overlay", settings.seed))
        probe = Chem.Mol(candidate.mol)
        pharm, counts = pharmacophore_overlap(reference.mol, probe, 0, cid)
        shape = 0.1 if fail_screening or phase["name"] == "refinement" else 0.9
        joint = math.sqrt(shape * pharm)
        return {
            "shape": shape,
            "pharmacophore": pharm,
            "alignment_color": 0.5,
            "joint": joint,
            "reference_adjusted_joint": joint,
            "reference_conformer_id": 0,
            "reference_relative_energy": 0.0,
            **counts,
        }, probe

    def score(*args, **kwargs):
        phase["name"] = kwargs.get("rejected_pose_diagnostic_phase", "screening")
        return original(*args, **kwargs)

    monkeypatch.setattr(workflow, "prepare_reference", lambda *a: None)
    monkeypatch.setattr(scorer, "require_alignment", lambda: None)
    monkeypatch.setattr(scorer, "ensemble", ensemble)
    monkeypatch.setattr(scorer, "_overlay", overlay)
    monkeypatch.setattr(scorer, "score_candidate", score)
    return calls


@pytest.mark.parametrize("fail_screening", [True, False])
@pytest.mark.parametrize("persist_fails", [True, False])
def test_actual_workflow_failure_artifacts_share_one_diagnostic_table(
    tmp_path, assets, proposal, monkeypatch, fail_screening, persist_fails
):
    from .test_selective_workflow import read

    workflow_control(monkeypatch, fail_screening=fail_screening)
    config = config_for(tmp_path, assets, proposal)
    config["rejected_pose_diagnostics"] = True
    if persist_fails:
        monkeypatch.setattr(
            diagnostics, "_write_bytes", lambda *a: (_ for _ in ()).throw(OSError("synthetic disk failure"))
        )
    out = tmp_path / "out"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    assert read(out, "results.json")["results"] == []
    summary = read(out, "summary.json")
    assert summary["status"] == "no_candidates"
    stats = summary["rejected_pose_diagnostics"]
    assert stats["export"]["status"] == ("unavailable" if persist_fails else "written")
    if fail_screening:
        failures = [
            r["similarity_failure"] for r in read(out, "routes.json")["routes"] if "similarity_failure" in r
        ]
        expected_phase, seeds = "screening", [42]
    else:
        reports = read(out, "refinement.json")["candidates"]
        failures = [t["similarity_failure"] for t in reports[0]["trials"]]
        expected_phase, seeds = "refinement", [42, 43, 44]
    refs = [f["scorer_audit"]["rejected_pose_diagnostic"] for f in failures]
    assert all(r["status"] == "captured_not_published" for r in refs)
    assert all("witness" not in r and "candidate_pose" not in r for r in refs)
    assert stats["captures"] == len(seeds)
    if not persist_fails:
        artifact = read(out, diagnostics.ARTIFACT)
        assert [r["settings"]["seed"] for r in artifact["records"]] == seeds
        assert {r["phase"] for r in artifact["records"]} == {expected_phase}
        assert all(r["identities"]["product"] == DESIGN for r in artifact["records"])
        assert all(isinstance(r["hit_tanimoto"], float) for r in artifact["records"])
        for ref in refs:
            assert ref["record_sha256"] == diagnostics.digest(artifact["records"][ref["record_index"]])
    else:
        assert not (out / diagnostics.ARTIFACT).exists()


def test_default_off_workflow_does_not_construct_collector(tmp_path, assets, proposal, monkeypatch):
    workflow_control(monkeypatch, fail_screening=True)
    monkeypatch.setattr(
        diagnostics, "RejectedPoseDiagnostics", lambda **k: pytest.fail("default-off collector")
    )
    out = tmp_path / "out"
    assert workflow.run(HIT, [DESIGN], config_for(tmp_path, assets, proposal), out) == []
    assert not (out / diagnostics.ARTIFACT).exists()
    assert "rejected_pose_diagnostics" not in json.loads((out / "summary.json").read_text())


def test_init_handoff_and_publication_import_overhead_use_same_allowance(tmp_path):
    clock = Clock()
    clock.now = 1.1
    collector = diagnostics.RejectedPoseDiagnostics(clock=clock, started_monotonic=0)
    assert collector.statistics()["setup_seconds"] == 1.1
    assert diagnostics.publish(collector, tmp_path)["status"] == "unavailable"
    clock = Clock()
    collector = diagnostics.RejectedPoseDiagnostics(clock=clock)
    clock.now = 1.1
    assert diagnostics.publish(collector, tmp_path, started_monotonic=0)["status"] == "unavailable"
    assert not (tmp_path / diagnostics.ARTIFACT).exists()


def test_no_extra_feature_matching_or_ensemble_requests(monkeypatch):
    counts = []
    original = constraints.feature_records

    def measured(*args, **kwargs):
        counts.append((args[0].GetNumAtoms(), kwargs.get("cid"), len(args)))
        return original(*args, **kwargs)

    monkeypatch.setattr(constraints, "feature_records", measured)
    before, _ = capture(monkeypatch, None)
    old_counts = list(counts)
    counts.clear()
    collector = diagnostics.RejectedPoseDiagnostics()
    after, _ = capture(monkeypatch, collector)
    assert counts == old_counts
    assert (
        [r["role"] for r in before.audit["ensemble_requests"]]
        == [r["role"] for r in after.audit["ensemble_requests"]]
        == ["candidate", "hit"]
    )


def test_failed_overlay_never_substitutes_diagnostic_pose(monkeypatch):
    collector = diagnostics.RejectedPoseDiagnostics()
    controlled(monkeypatch, {(0, 0): (0.69, 0.9)})
    monkeypatch.setattr(scorer, "_overlay", lambda *a: (_ for _ in ()).throw(ValueError("overlay failed")))
    with pytest.raises(scorer.SimilarityAbstention) as caught:
        scorer.score_candidate(HIT, HIT, DESIGN, rejected_pose_diagnostics=collector)
    assert caught.value.reason == "no_jointly_feasible_hit_pose"
    assert caught.value.audit["rejected_pose_diagnostic"]["reason"] == "no_completed_rejected_pair"
    assert not collector.snapshot()["records"]


def test_success_after_rejected_pair_serializes_no_rejected_evidence(monkeypatch):
    collector = diagnostics.RejectedPoseDiagnostics()
    controlled(monkeypatch, {(0, 0): (0.69, 0.9), (0, 1): (0.9, 0.9)})
    assert (
        scorer.score_candidate(HIT, HIT, DESIGN, rejected_pose_diagnostics=collector)[0]["hit"][
            "reference_conformer_id"
        ]
        == 1
    )
    assert not collector.snapshot()["records"]
    assert collector.statistics()["rejected_pairs_considered"] == 1


@pytest.mark.parametrize("field", ["identity", "policy", "reference", "cid"])
def test_bound_capture_rejects_contradictions(field):
    collector = diagnostics.RejectedPoseDiagnostics()
    call, *_ = pending(collector)
    if field == "identity":
        call.identities["product"] = "CC"
    elif field == "policy":
        call.policy["polar_recall"] = 0.9
    elif field == "reference":
        call.best[3]["reference_pose"]["coordinates"][0][0] += 0.25
    else:
        call.best[1]["reference_conformer_id"] = 1
    assert call.finish("no_jointly_feasible_hit_pose")["reason"] == "serialization_error:ValueError"
    assert not collector.snapshot()["records"]


def test_serialization_interruption_propagates(monkeypatch):
    collector = diagnostics.RejectedPoseDiagnostics()
    call, *_ = pending(collector)
    monkeypatch.setattr(constraints, "pose_binding", lambda *a: (_ for _ in ()).throw(KeyboardInterrupt()))
    with pytest.raises(KeyboardInterrupt):
        call.finish("no_jointly_feasible_hit_pose")
    assert not collector.snapshot()["records"]


def test_already_expired_publication_is_not_retryable(tmp_path):
    collector = diagnostics.RejectedPoseDiagnostics()
    assert diagnostics.publish(collector, tmp_path, deadline_monotonic=-1)["status"] == "unavailable"
    result = diagnostics.publish(collector, tmp_path)
    assert "already_attempted" in result["reason"]
    assert not (tmp_path / diagnostics.ARTIFACT).exists()


def test_expired_setup_does_not_start_source_reads(monkeypatch):
    monkeypatch.setattr(diagnostics, "source_provenance", lambda *a: pytest.fail("expired provenance read"))
    collector = diagnostics.RejectedPoseDiagnostics(deadline_monotonic=-1)
    assert collector._reason(None) == "provenance_unavailable:TimeoutError"


def test_shared_initial_deadline_cannot_be_removed_by_later_call(tmp_path):
    clock = Clock()
    collector = diagnostics.RejectedPoseDiagnostics(clock=clock, deadline_monotonic=1)
    clock.now = 2
    result = diagnostics.publish(collector, tmp_path)
    assert result["status"] == "unavailable"
    assert "workflow_deadline_exhausted" in result["reason"]
    assert not (tmp_path / diagnostics.ARTIFACT).exists()


def test_late_commit_reports_inherited_deadline_overrun(tmp_path, monkeypatch):
    clock = Clock()
    collector = diagnostics.RejectedPoseDiagnostics(clock=clock, deadline_monotonic=1)
    original = Path.replace

    def late(path, target):
        result = original(path, target)
        clock.now = 2
        return result

    monkeypatch.setattr(Path, "replace", late)
    result = diagnostics.publish(collector, tmp_path)
    assert result["status"] == "written" and result["deadline_overrun"]
    assert result["publication_completion_within_budget"] is None


def test_interrupted_publication_keeps_complete_core_bundle(tmp_path, assets, proposal, monkeypatch):
    workflow_control(monkeypatch, fail_screening=True)
    config = config_for(tmp_path, assets, proposal)
    config["rejected_pose_diagnostics"] = True
    out = tmp_path / "out"

    def interrupted(*args, **kwargs):
        assert json.loads((out / "summary.json").read_text())["status"] == "no_candidates"
        assert json.loads((out / "results.json").read_text())["results"] == []
        assert json.loads((out / "routes.json").read_text())["routes"]
        raise KeyboardInterrupt()

    monkeypatch.setattr(diagnostics, "publish", interrupted)
    with pytest.raises(KeyboardInterrupt):
        workflow.run(HIT, [DESIGN], config, out)
    summary = json.loads((out / "summary.json").read_text())
    assert summary["rejected_pose_diagnostics"]["export"]["reason"] == "optional_publication_pending"
    assert summary["status"] == "no_candidates"


@pytest.mark.parametrize("failure", ["write", "expired", "interrupt"])
def test_optional_summary_update_preserves_original_bytes(tmp_path, monkeypatch, failure):
    initial = {"status": "no_candidates", "rejected_pose_diagnostics": {"export": {"status": "unavailable"}}}
    path = tmp_path / "summary.json"
    original = json.dumps(initial).encode()
    path.write_bytes(original)
    clock = Clock()
    collector = diagnostics.RejectedPoseDiagnostics(clock=clock)
    if failure == "expired":
        collector.charge_overhead(2)
    else:
        exc = KeyboardInterrupt if failure == "interrupt" else OSError
        monkeypatch.setattr(diagnostics, "_write_bytes", lambda *a: (_ for _ in ()).throw(exc()))
    if failure == "interrupt":
        with pytest.raises(KeyboardInterrupt):
            diagnostics.complete_summary(
                collector, tmp_path, initial, {"status": "written"}, workflow_started=0
            )
    else:
        assert (
            diagnostics.complete_summary(
                collector, tmp_path, initial, {"status": "written"}, workflow_started=0
            )["status"]
            == "unavailable"
        )
    assert path.read_bytes() == original


def test_late_summary_update_never_self_certifies_completion(tmp_path, monkeypatch):
    clock = Clock()
    collector = diagnostics.RejectedPoseDiagnostics(clock=clock)
    original = Path.replace

    def late(path, target):
        result = original(path, target)
        clock.now = 2
        return result

    monkeypatch.setattr(Path, "replace", late)
    result = diagnostics.complete_summary(
        collector, tmp_path, {"status": "no_candidates"}, {"status": "written"}, workflow_started=0
    )
    assert result["status"] == "written" and result["time_overrun"]
    assert result["completion_within_budget"] is None
    saved = json.loads((tmp_path / "summary.json").read_text())
    assert saved["rejected_pose_diagnostics"]["summary_update_completion_within_budget"] is None


def test_unexpected_optional_publication_error_preserves_core(tmp_path, assets, proposal, monkeypatch):
    workflow_control(monkeypatch, fail_screening=True)
    config = config_for(tmp_path, assets, proposal)
    config["rejected_pose_diagnostics"] = True
    monkeypatch.setattr(
        diagnostics, "publish", lambda *a, **k: (_ for _ in ()).throw(OSError("synthetic optional fault"))
    )
    out = tmp_path / "out"
    assert workflow.run(HIT, [DESIGN], config, out) == []
    summary = json.loads((out / "summary.json").read_text())
    assert summary["status"] == "no_candidates"
    assert summary["rejected_pose_diagnostics"]["export"]["status"] == "unavailable"
