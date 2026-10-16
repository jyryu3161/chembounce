"""Necessary-condition efficiency only; no empirical yield or activity claim."""

import json
from collections import Counter
from types import SimpleNamespace

import pytest

from chembounce_synthesis import polar_capacity, selective_scoring, workflow
from chembounce_synthesis.catalog import Catalog
from chembounce_synthesis.common import canonical
from chembounce_synthesis.polar_capacity import CapacityRejection, PolarCapacity
from chembounce_synthesis.route_objective import RootProductObjective, search_policy
from chembounce_synthesis.routes import SearchLimits, expand_route, resolve_proposal, verify_route
from chembounce_synthesis.similarity_constraints import BoundSimilarity, POLAR, VERSION

HIT = "CC(=O)Nc1ccncc1"
POOR = "CC(=O)Nc1ccccc1"
GOOD = "CC(=O)Nc1nccnc1"
MINIMA = dict(hit_tanimoto=0.01, hit_shape=0.01, hit_pharmacophore=0.01)


def test_real_typed_capacity_and_overlapping_atom_features():
    gate = PolarCapacity(BoundSimilarity("NCCN"))
    result = gate.assess("CN")
    # Each nitrogen contributes Donor and PosIonizable, as distinct features.
    assert result["reference_counts"] == {"Donor": 2, "PosIonizable": 2}
    assert result["candidate_counts"] == {"Donor": 1, "PosIonizable": 1}
    assert (result["numerator"], result["denominator"], result["recall_upper_bound"]) == (2, 4, 0.5)
    assert result["status"] == "capacity_rejected"
    with pytest.raises(CapacityRejection) as exc:
        gate.require("CN")
    assert exc.value.decision == result
    assert gate.assess("NCCN")["eligible"]


@pytest.mark.parametrize("threshold,eligible", [(0.1, True), (0.10000000000000002, False), (0.09, True)])
def test_float_equality_and_same_family_distinct_columns(monkeypatch, threshold, eligible):
    # Distinct same-family features may overlap atoms; count feature columns.
    reference = [{"family": "Acceptor", "atom_ids": [0, i + 1]} for i in range(10)]
    candidate = [{"family": "Acceptor", "atom_ids": [0]}]
    monkeypatch.setattr(polar_capacity, "feature_records", lambda mol: reference)
    bound = SimpleNamespace(
        hit="CO",
        sha256="fixture",
        policy={"minimum_polar_reference_recall": threshold},
        preflight_candidate=lambda product: candidate,
    )
    result = PolarCapacity(bound).assess("CO")
    assert result["numerator"] == 1 and result["denominator"] == 10
    assert result["recall_upper_bound"] == 0.1
    assert result["eligible"] is eligible


def test_no_polar_and_empty_features_are_distinct():
    gate = PolarCapacity(BoundSimilarity("c1ccccc1"))
    result = gate.assess("c1ccncc1")
    assert result["status"] == "not_applicable" and result["eligible"]
    assert result["denominator"] == 0 and result["recall_upper_bound"] is None
    rejected = gate.assess("[Na+]")
    assert rejected["status"] == "preflight_rejected"
    assert "nonempty candidate features" in rejected["reason"]
    assert gate.snapshot()["counts"]["capacity_rejections"] == 0


@pytest.mark.parametrize("smiles", ["[2H]OC", "[H]OC", "N[C@@H](C)C(=O)O", "F/C=C/F"])
def test_exact_h_isotope_and_specified_stereo_semantics(smiles):
    bound = BoundSimilarity(smiles)
    result = PolarCapacity(bound).assess(smiles)
    expected = Counter(f["family"] for f in bound.preflight_candidate(smiles) if f["family"] in POLAR)
    assert result["product"] == canonical(smiles)
    assert result["reference_counts"] == result["candidate_counts"] == dict(expected)
    assert result["eligible"]


@pytest.mark.parametrize("product,reason", [("NC(C)C(=O)O", "stereo"), ("CCO", "core")])
def test_state_and_core_errors_are_preflight_not_capacity(product, reason):
    gate = PolarCapacity(BoundSimilarity("N[C@@H](C)C(=O)O", core_smiles="N[C@@H](C)C(=O)O"))
    result = gate.assess(product)
    assert result["status"] == "preflight_rejected"
    assert reason in result["reason"].lower()
    assert "recall_upper_bound" not in result


@pytest.mark.parametrize(
    "reason", ["Similarity feature bound exceeded (256)", "Ambiguous duplicate feature identity"]
)
def test_preflight_feature_errors_never_truncated_or_reclassified(monkeypatch, reason):
    bound = BoundSimilarity(HIT)

    def fail(product):
        raise ValueError(reason)

    monkeypatch.setattr(bound, "preflight_candidate", fail)
    gate = PolarCapacity(bound)
    assert gate.assess(GOOD) == {
        "product": canonical(GOOD),
        "constraint_policy_sha256": bound.sha256,
        "eligible": False,
        "status": "preflight_rejected",
        "reason": reason,
    }
    assert gate.snapshot()["counts"]["preflight_failures"] == 1


def test_landmarks_do_not_tighten_optimistic_capacity_bound():
    bound = BoundSimilarity(
        "NCCN",
        polar_recall=0.5,
        landmarks=[
            {"id": "left", "family": "Donor", "hit_atom_ids": [0]},
            {"id": "right", "family": "Donor", "hit_atom_ids": [3]},
        ],
    )
    result = PolarCapacity(bound).assess("CN")
    assert result["eligible"] and result["recall_upper_bound"] == 0.5
    # One donor cannot fulfill both mandatory donor landmarks: final geometry
    # still decides. This deliberately does not claim sufficient acceptance.
    assert len(bound.policy["landmarks"]) == 2


def test_lru_is_bounded_per_run_and_defensive_with_honest_counts():
    bound = BoundSimilarity(HIT)
    gate = PolarCapacity(bound, cache_size=1)
    first = gate.assess(POOR)
    first["reference_counts"].clear()
    assert gate.assess(POOR)["reference_counts"]
    gate.assess(GOOD)
    gate.assess(POOR)
    audit = gate.snapshot()
    assert audit["counts"] == dict(
        requests=4,
        evaluations=3,
        cache_hits=1,
        evictions=2,
        preflight_failures=0,
        capacity_rejections=2,
        deadline_events=0,
    )
    assert audit["current_cached_products"] == len(audit["decision_sample"]) == 1
    assert audit["sample_at_capacity"]
    audit["decision_sample"][0]["reference_counts"].clear()
    assert gate.snapshot()["decision_sample"][0]["reference_counts"]
    assert PolarCapacity(bound).snapshot()["counts"]["requests"] == 0


@pytest.mark.parametrize(
    "phase", ["before_reference", "after_reference", "before_candidate", "after_candidate", "cache_hit"]
)
def test_deadline_never_turns_incomplete_perception_into_pass(monkeypatch, phase):
    state = {"expired": phase == "before_reference"}
    original = polar_capacity.feature_records

    def perceive(mol):
        records = original(mol)
        if phase == "after_reference":
            state["expired"] = True
        return records

    monkeypatch.setattr(polar_capacity, "feature_records", perceive)
    bound = BoundSimilarity(HIT)
    gate = PolarCapacity(bound, deadline_reached=lambda: state["expired"])
    original_candidate = bound.preflight_candidate

    def candidate(product):
        records = original_candidate(product)
        if phase == "after_candidate":
            state["expired"] = True
        return records

    monkeypatch.setattr(bound, "preflight_candidate", candidate)
    if phase == "cache_hit":
        assert gate.assess(GOOD)["eligible"]
    if phase in {"before_candidate", "cache_hit"}:
        state["expired"] = True
    result = gate.assess(GOOD)
    assert not result["eligible"] and result["status"] == "not_evaluated_deadline"
    assert gate.snapshot()["counts"]["deadline_events"] == 1
    assert result["preflight_completed"] is (phase == "after_candidate")
    assert result["capacity_computed"] is (phase == "after_candidate")
    if phase == "after_candidate":
        assert gate.snapshot()["counts"]["evaluations"] == 1
        assert result["recall_upper_bound"] == 1
    if phase != "cache_hit":
        assert gate.snapshot()["current_cached_products"] == 0


def objective(hit, design, gate, rounds=0):
    def check(product):
        workflow.check_constraints(product, hit, design)
        gate.require(product)

    return RootProductObjective(
        hit,
        design,
        MINIMA,
        check,
        search_policy({"recommendation": {"mode": "qualified"}, "expansion": {"rounds": rounds}}),
    )


def polar_catalog(assets):
    path = assets[0].path
    path.write_text(
        "id\tsmiles\tsource\tsnapshot\n"
        + "".join(
            f"{ident}\t{smiles}\tunit fixture\t2026-10-17\n"
            for ident, smiles in [("acetic", "CC(=O)O"), ("phenyl", "Nc1ccccc1"), ("diazine", "Nc1nccnc1")]
        )
    )
    return Catalog(path), assets[1]


def proposal_for(product=POOR, acid="CC(=O)O", amine="Nc1ccccc1"):
    return {
        "reactions": [{"reaction_template": "amide", "reactants": [acid, amine], "product": product}],
        "building_blocks": [acid, amine],
    }


def test_capacity_eligibility_precedes_complete_root_beam_pruning(assets):
    assets = polar_catalog(assets)
    limits = SearchLimits(building_block_top_k=2, beam_width=1, max_combinations=4)
    legacy, before = resolve_proposal(proposal_for(), *assets, limits)
    gate = PolarCapacity(BoundSimilarity(HIT, MINIMA))
    guided, after = resolve_proposal(
        proposal_for(), *assets, limits, root_objective=objective(HIT, POOR, gate)
    )
    assert legacy[0]["verification"]["product"] == POOR
    assert guided[0]["verification"]["product"] == canonical(GOOD)
    assert before.combinations == after.combinations <= 4
    assert gate.assess(POOR)["status"] == "capacity_rejected"
    assert verify_route(guided[0]["route"], *assets)["product"] == canonical(GOOD)


def test_capacity_poor_seed_still_expands_to_eligible_product(assets, route):
    assets = polar_catalog(assets)
    route["nodes"][1]["building_block_id"] = "phenyl"
    route["nodes"][-1]["product"] = POOR
    gate = PolarCapacity(BoundSimilarity(HIT, MINIMA))
    obj = objective(HIT, POOR, gate, rounds=1)
    initial, _ = resolve_proposal(
        proposal_for(), *assets, SearchLimits(building_block_top_k=1, beam_width=1), root_objective=obj
    )
    assert initial and not initial[0]["search_eligibility"]["eligible"]
    expanded, report = expand_route(
        route,
        *assets,
        SearchLimits(beam_width=1, max_combinations=2),
        rounds=1,
        max_routes=1,
        building_block_top_k=1,
        target_smiles=POOR,
        root_objective=obj,
    )
    assert expanded and expanded[0]["verification"]["product"] == canonical(GOOD)
    assert expanded[0]["search_eligibility"]["eligible"]
    assert report.combinations <= 2


def test_intermediate_feature_deficit_is_not_pruned(assets):
    hit = "CC(=O)NCC(=O)Nc1nccnc1"
    final = "CC(=O)NCC(=O)Nc1ncncn1"
    with assets[0].path.open("a") as stream:
        stream.write("triazine\tNc1ncncn1\ttest fixture\t2026-10-17\n")
    assets = Catalog(assets[0].path), assets[1]
    intermediate = "CC(=O)NCC(=O)O"
    gate = PolarCapacity(BoundSimilarity(hit, MINIMA, polar_recall=1.0))
    proposal = {
        "reactions": [
            {"reaction_template": "amide", "reactants": ["CC(=O)O", "NCC(=O)O"], "product": intermediate},
            {"reaction_template": "amide", "reactants": [intermediate, "Nc1ncncn1"], "product": final},
        ],
        "building_blocks": ["CC(=O)O", "NCC(=O)O", "Nc1ncncn1"],
    }
    roots, _ = resolve_proposal(
        proposal,
        *assets,
        SearchLimits(building_block_top_k=1, beam_width=1),
        root_objective=objective(hit, final, gate),
    )
    assert roots and verify_route(roots[0]["route"], *assets)["product"] == final
    assert all(d["product"] != intermediate for d in gate.snapshot()["decision_sample"])

    assert not gate.assess(intermediate)["eligible"]


def workflow_config(tmp_path, assets, proposals):
    path = tmp_path / "proposals.json"
    path.write_text(json.dumps({"schema_version": 1, "proposals": proposals}))
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "recommendation": {"mode": "qualified"},
        "minimum_scores": MINIMA,
        "similarity_constraints": {"version": VERSION},
        "conformers": {"num_conformers": 8, "scoring_objective": "original-hit-v1"},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "refinement": {"top_n": 1, "num_conformers": 8, "seeds": [42, 43, 44], "max_total_conformers": 24},
        "max_scored_products": 1,
        "max_results": 1,
    }


@pytest.mark.parametrize("search", ["original-hit-v1", "design-v1"])
def test_impossible_final_never_consumes_evidence_or_score_slot(tmp_path, assets, monkeypatch, search):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    assets = polar_catalog(assets)
    config = workflow_config(
        tmp_path,
        assets,
        [
            {"design": POOR, "response": proposal_for()},
            {"design": GOOD, "response": proposal_for(GOOD, amine="Nc1nccnc1")},
        ],
    )
    config["route_search_policy"] = search
    config["expansion"] = {"rounds": 1, "max_routes": 1, "building_block_top_k": 1}
    evidence, scoring = [], []
    original = workflow._evidence_evaluator

    def evaluator(config):
        result = original(config)
        evaluate = result.evaluate

        def track(verification, **kwargs):
            evidence.append(verification["product"])
            return evaluate(verification, **kwargs)

        result.evaluate = track
        return result

    monkeypatch.setattr(workflow, "_evidence_evaluator", evaluator)

    def score(product, *args, **kwargs):
        scoring.append(product)
        raise ValueError("Controlled stop after allocated scoring; no 3D calculation")

    monkeypatch.setattr(selective_scoring, "score_candidate", score)
    out = tmp_path / "run"
    assert workflow.run(HIT, [POOR, GOOD], config, out) == []
    summary = json.loads((out / "summary.json").read_text())
    records = json.loads((out / "routes.json").read_text())["routes"]
    poor = [r for r in records if r.get("verification", {}).get("product") == POOR]
    assert poor and all(r["polar_capacity"]["status"] == "capacity_rejected" for r in poor)
    assert all(not r["eligible_for_evidence"] and not r["eligible_for_scoring"] for r in poor)
    assert POOR not in evidence and scoring == [canonical(GOOD)]
    assert evidence and set(evidence) == {canonical(GOOD)}
    assert summary["scored_products"] == 1
    assert summary["polar_capacity"]["counts"]["capacity_rejections"] >= 1
    manifest = json.loads((out / "manifest.json").read_text())
    assert "polar_capacity.py" in manifest["code_sha256"]
    assert manifest["polar_capacity"]["version"] == "typed-polar-capacity-v1"


def test_selective_off_never_constructs_capacity_gate(tmp_path, assets, monkeypatch):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = workflow_config(tmp_path, assets, [{"design": POOR, "response": proposal_for()}])
    del config["similarity_constraints"]

    def forbidden(*args, **kwargs):
        pytest.fail("Capacity must not activate without selective constraints")

    monkeypatch.setattr(workflow, "PolarCapacity", forbidden)

    def stop(*args, **kwargs):
        raise ValueError("Controlled scoring stop")

    monkeypatch.setattr(workflow, "score_candidate", stop)
    out = tmp_path / "legacy"
    workflow.run(HIT, [POOR], config, out)
    assert "polar_capacity" not in json.loads((out / "summary.json").read_text())


def test_carried_positive_flags_never_bypass_fresh_route_and_capacity_checks(
    tmp_path, assets, route, monkeypatch
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    assets = polar_catalog(assets)
    route["nodes"][1]["building_block_id"] = "phenyl"
    route["nodes"][-1]["product"] = POOR
    source = tmp_path / "carried"
    source.mkdir()
    (source / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "hit": HIT,
                "designs": [POOR],
                "config": {},
            }
        )
    )
    (source / "routes.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "routes": [
                    {
                        "route_id": "prior",
                        "design": POOR,
                        "route": route,
                        "polar_capacity": {"eligible": True},
                        "eligible_for_scoring": True,
                    }
                ],
            }
        )
    )
    config = workflow_config(tmp_path, assets, [])
    config["carry_forward"] = [str(source)]
    verified = []
    replay = workflow.verify_route

    def verify(graph, *args, **kwargs):
        verified.append(graph)
        return replay(graph, *args, **kwargs)

    monkeypatch.setattr(workflow, "verify_route", verify)
    monkeypatch.setattr(
        selective_scoring, "score_candidate", lambda *a, **k: pytest.fail("Impossible carried final scored")
    )
    out = tmp_path / "continued"
    assert workflow.run(HIT, [POOR], config, out) == []
    records = json.loads((out / "routes.json").read_text())["routes"]
    assert len(verified) == len(records) == 1
    record = records[0]
    assert record["carried_from"]["route_id"] == "prior"
    assert record["selection_status"] == "selective_capacity_rejected"
    assert record["polar_capacity"]["recall_upper_bound"] == 2 / 3
    assert not record["eligible_for_scoring"]
    assert "synthesis_evidence" not in record


def test_completed_rejection_is_counted_but_not_applied_or_cached_after_deadline(monkeypatch):
    bound = BoundSimilarity(HIT)
    state = {"expired": False}
    gate = PolarCapacity(bound, deadline_reached=lambda: state["expired"])
    preflight = bound.preflight_candidate

    def expire(product):
        records = preflight(product)
        state["expired"] = True
        return records

    monkeypatch.setattr(bound, "preflight_candidate", expire)
    decision = gate.assess(POOR)
    assert decision["status"] == "not_evaluated_deadline"
    assert decision["preflight_completed"] and decision["capacity_computed"]
    assert decision["recall_upper_bound"] == 2 / 3
    audit = gate.snapshot()
    assert audit["counts"]["evaluations"] == audit["counts"]["capacity_rejections"] == 1
    assert audit["counts"]["deadline_events"] == 1
    assert audit["current_cached_products"] == 0 and audit["decision_sample"] == []
