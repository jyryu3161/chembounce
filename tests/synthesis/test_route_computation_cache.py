"""Pure-work reuse with exact chemistry and fresh admission checks.

Tiny synthetic routes only; call counts are not wall-time speedup measurements.
"""

import copy
import json
import sys
from collections import Counter
from types import SimpleNamespace

import pytest
from rdkit import Chem, DataStructs

from chembounce_synthesis import routes, workflow
from chembounce_synthesis.common import FP, canonical, molecule
from .test_route_objective import DESIGN, HIT, objective


@pytest.mark.parametrize(
    "left,right",
    [
        ("N[C@H](C)C(=O)O", "N[C@@H](C)C(=O)O"),
        ("CCO", "[13CH3]CO"),
        ("CN", "C[NH3+]"),
        ("CC=O", "C=CO"),
    ],
)
def test_property_cache_keeps_exact_states_separate(left, right):
    calls = []

    def calculate(mol):
        exact = Chem.MolToSmiles(mol, isomericSmiles=True)
        calls.append(exact)
        return {"exact": exact}

    cache = workflow._ProductProperties(calculate, 2)
    for _ in range(2):
        assert cache(molecule(left)) == {"exact": canonical(left)}
        assert cache(molecule(right)) == {"exact": canonical(right)}
    assert calls == [canonical(left), canonical(right)]


def test_property_cache_is_lru_scoped_and_copies_values():
    calls, retained = [], {"descriptor": [1]}

    def calculate(mol):
        calls.append(Chem.MolToSmiles(mol))
        return retained

    cache = workflow._ProductProperties(calculate, 2)
    first = cache(molecule("CCO"))
    first["descriptor"].append(2)
    retained["descriptor"].append(3)
    assert cache(molecule("OCC")) == {"descriptor": [1]}
    cache(molecule("CCN"))
    cache(molecule("CCO"))  # Touch ethanol, so ethylamine is least recent.
    cache(molecule("CCC"))
    cache(molecule("CCO"))
    cache(molecule("CCN"))
    assert calls == ["CCO", "CCN", "CCC", "CCN"]
    workflow._ProductProperties(calculate, 2)(molecule("CCO"))
    assert calls[-1] == "CCO" and len(calls) == 5


def test_property_calculation_failures_are_not_cached():
    attempts = []

    def calculate(mol):
        attempts.append(1)
        if len(attempts) == 1:
            raise ValueError("temporary test failure")
        return {"mw": 46}

    cache = workflow._ProductProperties(calculate, 2)
    with pytest.raises(ValueError, match="temporary test failure"):
        cache(molecule("CCO"))
    assert cache(molecule("CCO")) == {"mw": 46}
    assert cache(molecule("CCO")) == {"mw": 46}
    assert len(attempts) == 2


def repeated_inputs(tmp_path, assets, proposal):
    designs = ["CC(=O)Nc1ccncc1", DESIGN]
    path = tmp_path / "repeated.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "proposals": [
                    {"design": design, "response": proposal} for design in designs for _ in range(2)
                ],
            }
        )
    )
    return designs, {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "route_search_policy": "original-hit-v1",
        "route_search_cache_size": 2,
        "search": {"building_block_top_k": 1, "beam_width": 1, "max_combinations": 4},
        "refinement": {"top_n": 0},
        "max_scored_products": 2,
    }


def test_workflow_reuses_only_properties_across_routes_designs_and_not_runs(
    tmp_path, assets, proposal, monkeypatch
):
    import utils

    designs, config = repeated_inputs(tmp_path, assets, proposal)
    calculate, factory = utils._get_molecular_prop_extd_, workflow._ProductProperties
    measurements, outputs = {}, {}
    for mode in ("uncached", "cached", "next_run"):
        calls = []

        def counted(mol):
            calls.append(Chem.MolToSmiles(mol))
            return calculate(mol)

        with monkeypatch.context() as scoped:
            scoped.setattr(utils, "_get_molecular_prop_extd_", counted)
            scoped.setattr(
                workflow,
                "_ProductProperties",
                (lambda calculate, maxsize: calculate) if mode == "uncached" else factory,
            )
            out = tmp_path / mode
            rows = workflow.run(HIT, designs, config, out, scoring_mode="2d")
        assert rows
        outputs[mode] = {"rows": rows}
        for name in ("routes", "scores", "search_objectives"):
            outputs[mode][name] = json.loads((out / (name + ".json")).read_text())
        measurements[mode] = dict(Counter(calls))
    assert outputs["uncached"] == outputs["cached"] == outputs["next_run"]
    product = canonical(designs[0])
    assert measurements == {"uncached": {product: 6}, "cached": {product: 1}, "next_run": {product: 1}}
    assert len(outputs["cached"]["routes"]["routes"]) == 4
    print("PROPERTY_CALL_COUNTS " + json.dumps(measurements, sort_keys=True))


@pytest.mark.parametrize("veto", ["constraint", "property"])
def test_workflow_rechecks_dynamic_veto_after_property_cache_hit(
    tmp_path, assets, proposal, monkeypatch, veto
):
    import utils

    designs, config = repeated_inputs(tmp_path, assets, proposal)
    config["route_search_policy"] = "design-v1"
    calls, calculations = [], []
    constraints, thresholds = workflow._check_constraints, utils.determ_props
    calculate = utils._get_molecular_prop_extd_

    def check(*args, **kwargs):
        calls.append(1)
        if len(calls) == 3:
            if veto == "constraint":
                raise ValueError("fresh dynamic veto")
            return False, ["fresh dynamic veto"]
        return (constraints if veto == "constraint" else thresholds)(*args, **kwargs)

    def counted(mol):
        calculations.append(1)
        return calculate(mol)

    monkeypatch.setattr(
        workflow if veto == "constraint" else utils,
        "_check_constraints" if veto == "constraint" else "determ_props",
        check,
    )
    monkeypatch.setattr(utils, "_get_molecular_prop_extd_", counted)
    out = tmp_path / veto
    workflow.run(HIT, designs, config, out, scoring_mode="2d")
    records = json.loads((out / "routes.json").read_text())["routes"]
    assert len(calls) == 4 and len(calculations) == 1
    assert records[2]["selection_status"] == "constraint_rejected"
    assert "fresh dynamic veto" in records[2]["reason"]
    assert all(r["selection_status"] == "accepted" for r in (records[0], records[1], records[3]))


def test_workflow_deadline_still_stops_repeated_cached_candidate(tmp_path, assets, proposal, monkeypatch):
    import utils

    designs, config = repeated_inputs(tmp_path, assets, proposal)
    config["route_search_policy"] = "design-v1"
    clock, calls = SimpleNamespace(now=0.0), []
    calculate = utils._get_molecular_prop_extd_
    thresholds = utils.determ_props

    def checked(props, *args, **kwargs):
        calls.append(1)
        if len(calls) == 2:  # Second route has already hit the property cache.
            clock.now = 20.0
        return thresholds(props, *args, **kwargs)

    calculations = []

    def counted(mol):
        calculations.append(1)
        return calculate(mol)

    monkeypatch.setattr(workflow.time, "monotonic", lambda: clock.now)
    monkeypatch.setattr(utils, "determ_props", checked)
    monkeypatch.setattr(utils, "_get_molecular_prop_extd_", counted)
    out = tmp_path / "deadline"
    assert workflow.run(HIT, designs, config, out, scoring_mode="2d", deadline_monotonic=10.0) == []
    summary = json.loads((out / "summary.json").read_text())
    assert summary["deadline_reached"] and summary["scored_products"] == 0
    assert len(calls) == 2 and len(calculations) == 1


@pytest.mark.parametrize("cache_size", [1, 4096])
def test_expansion_reuses_guidance_score_without_changing_requests_or_eviction(
    assets, route, monkeypatch, cache_size
):
    route["nodes"][0]["building_block_id"] = "propionic"
    route["nodes"][-1]["product"] = DESIGN
    records, measurements = {}, {}
    for mode in ("independent_score", "reuse"):
        obj = objective(rounds=2, cache_size=cache_size)
        # The proxy delegates the identical objective but disables the reference
        # identity fast path, retaining the independent fingerprint-score oracle.
        supplied = (
            SimpleNamespace(design=HIT, assess=obj.assess, policy=obj.policy)
            if mode == "independent_score"
            else obj
        )
        calls = []

        def fingerprint(mol):
            calls.append(sys._getframe(1).f_code.co_name)
            return FP.GetFingerprint(mol)

        with monkeypatch.context() as scoped:
            scoped.setattr(routes, "FP", SimpleNamespace(GetFingerprint=fingerprint))
            result = routes.expand_route(
                copy.deepcopy(route),
                *assets,
                routes.SearchLimits(beam_width=2, max_combinations=8),
                rounds=2,
                max_routes=4,
                building_block_top_k=2,
                target_smiles="O=C(CC)Nc1ccncc1",
                root_objective=supplied,
            )
        assert result[0]
        records[mode] = (result, obj.snapshot(), list(obj.cache))
        measurements[mode] = dict(Counter(calls))
    assert records["independent_score"] == records["reuse"]
    assert measurements["independent_score"]["expand_route"] == 8
    assert "expand_route" not in measurements["reuse"]
    # Intermediate-node ranking still needs its own fingerprints; only complete
    # products with the identical design reference use the reused objective score.
    assert {
        k: n for k, n in measurements["independent_score"].items() if k != "expand_route"
    } == measurements["reuse"]
    if cache_size == 1:
        assert records["reuse"][1]["counts"]["evictions"] > 0
    print("EXPANSION_FINGERPRINT_CALL_COUNTS " + json.dumps({"cache_size": cache_size, **measurements}))


@pytest.mark.parametrize("guided", [False, True])
def test_expansion_keeps_independent_score_when_reference_does_not_match(assets, route, guided):
    obj = objective(rounds=1) if guided else None
    items, _ = routes.expand_route(
        route,
        *assets,
        routes.SearchLimits(beam_width=2, max_combinations=4),
        rounds=1,
        max_routes=4,
        building_block_top_k=2,
        target_smiles=HIT,
        root_objective=obj,
    )
    assert items
    reference = FP.GetFingerprint(molecule(HIT))
    for item in items:
        expected = DataStructs.TanimotoSimilarity(
            reference, FP.GetFingerprint(molecule(item["verification"]["product"]))
        )
        assert item["target_tanimoto"] == expected
        if guided:
            assert item["target_tanimoto"] == item["search_eligibility"]["hit_tanimoto"]
