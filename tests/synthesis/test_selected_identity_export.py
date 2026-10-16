"""Retained identity comparisons are descriptive exports, never new acceptance."""

import copy
import inspect
import json
from dataclasses import replace

import pytest

from chembounce_synthesis import identity, workflow
from chembounce_synthesis.common import canonical


HIT = "CC(=O)Nc1ccccc1"
PRODUCT = "CC(=O)Nc1ccncc1"
RECORDED_ONLY = {
    "original_hit",
    "final_product",
    "change_class",
    "atom_aware_scaffold_changed",
    "generic_ring_linker_topology_changed",
    "scaffold_presence_change",
}


def forbidden(*args, **kwargs):
    raise AssertionError("Comparison export invoked chemistry")


def retained(reference=HIT, product=PRODUCT):
    original, candidate = identity.chemical_identity(reference), identity.chemical_identity(product)
    change = identity.classify_change(reference, product)
    state = workflow._SelectedIdentityEvidence(original.exact_smiles, None)
    state.capture_hit()
    admission = {
        "original_hit_smiles": original.exact_smiles,
        "final_product_smiles": candidate.exact_smiles,
        "change_class": change,
    }
    state.retain_product(candidate.exact_smiles, candidate, admission, workflow._identity_signature())
    row = {
        "smiles": candidate.exact_smiles,
        "family_smiles": candidate.family_smiles,
        "scaffold_family_smiles": candidate.scaffold_family_smiles,
        "design": "discarded design",
        "route_id": "alternate-primary",
        "score": 0.5,
    }
    return state, row, admission


def block_export_chemistry(monkeypatch):
    identity.clear_identity_caches()
    for module, names in (
        (
            workflow,
            (
                "chemical_identity",
                "classify_change",
                "scaffold_family_key",
                "canonical",
                "molecule",
                "_identity_signature",
            ),
        ),
        (identity, ("chemical_identity", "classify_change", "family_key", "scaffold_family_key")),
    ):
        for name in names:
            monkeypatch.setattr(module, name, forbidden)


@pytest.mark.parametrize(
    "reference,product,change,atom,generic,presence",
    [
        ("Oc1ccccc1", "c1ccccc1O", "unchanged", False, False, "both_present"),
        ("O=c1cccc[nH]1", "Oc1ccccn1", "state_only", False, False, "both_present"),
        ("CC(=O)O", "CC(=O)[O-]", "state_only", False, False, "both_absent"),
        ("CC", "[13CH3]C", "state_only", False, False, "both_absent"),
        ("C[C@H](O)C(=O)O", "C[C@@H](O)C(=O)O", "state_only", False, False, "both_absent"),
        ("Cc1ccccc1", "CCc1ccccc1", "decoration", False, False, "both_present"),
        ("CC", "CCC", "decoration", False, False, "both_absent"),
        ("CC", "CC", "unchanged", False, False, "both_absent"),
        ("c1ccccc1", "c1ccncc1", "scaffold", True, False, "both_present"),
        ("C1CCCCC1", "C1CC1", "scaffold", True, True, "both_present"),
        ("CC", "C1CC1", "scaffold", True, True, "gained"),
        ("C1CC1", "CC", "scaffold", True, True, "lost"),
    ],
)
def test_retained_classification_and_graph_semantics_without_export_chemistry(
    monkeypatch, reference, product, change, atom, generic, presence
):
    state, row, _ = retained(reference, product)
    before = copy.deepcopy(row)
    block_export_chemistry(monkeypatch)
    state.attach_selected([row])
    result = row.pop("hit_identity_comparison")
    assert row == before
    assert result["status"] == "recorded"
    assert result["change_class"] == change
    assert result["atom_aware_scaffold_changed"] is atom
    assert result["generic_ring_linker_topology_changed"] is generic
    assert result["scaffold_presence_change"] == presence
    assert result["original_hit"]["exact_smiles"] == result["original_hit_smiles"] == state.hit
    assert result["final_product"]["exact_smiles"] == result["final_product_smiles"] == row["smiles"]


@pytest.mark.parametrize(
    "problem,code",
    [
        ("missing_hit", "original_hit_not_retained"),
        ("missing_product", "final_product_not_retained"),
        ("missing_class", "original_hit_classification_not_retained"),
        ("wrong_reference", "retained_identity_binding_mismatch"),
        ("wrong_product", "retained_identity_binding_mismatch"),
        ("wrong_identity_exact", "retained_identity_binding_mismatch"),
        ("wrong_family_row", "retained_identity_binding_mismatch"),
        ("missing_family_row", "retained_identity_binding_mismatch"),
        ("missing_scaffold_row", "retained_identity_binding_mismatch"),
        ("wrong_descriptor", "retained_identity_policy_mismatch"),
        ("missing_descriptor", "retained_identity_policy_mismatch"),
        ("wrong_runtime", "retained_identity_policy_mismatch"),
        ("conflicting_class", "retained_identity_conflict"),
        ("conflicting_identity", "retained_identity_conflict"),
        ("conflicting_signature", "retained_identity_conflict"),
    ],
)
def test_missing_or_conflicting_values_are_atomic_unavailable_without_normalization(
    monkeypatch, problem, code
):
    state, row, admission = retained()
    product = row["smiles"]
    value = list(state.products[product])
    if problem == "missing_hit":
        state.hit_identity = None
    elif problem == "missing_product":
        state.products.clear()
    elif problem == "missing_class":
        value[2] = None
    elif problem == "wrong_reference":
        value[0] = "wrong reference"
    elif problem == "wrong_product":
        value[1] = "wrong product"
    elif problem == "wrong_identity_exact":
        value[3] = replace(value[3], exact_smiles="wrong product")
    elif problem == "wrong_family_row":
        row["family_smiles"] = "wrong family"
    elif problem == "missing_family_row":
        row.pop("family_smiles")
    elif problem == "missing_scaffold_row":
        row.pop("scaffold_family_smiles")
    elif problem == "wrong_descriptor":
        state.descriptor["tautomer_limit"] += 1
    elif problem == "missing_descriptor":
        state.descriptor = None
    elif problem == "wrong_runtime":
        value[3] = replace(value[3], rdkit_version="other")
    else:
        changed_admission = dict(admission)
        changed_identity, changed_signature = value[3], value[4]
        if problem == "conflicting_class":
            changed_admission["change_class"] = "decoration"
        elif problem == "conflicting_identity":
            changed_identity = replace(value[3], ring_skeleton_smiles="conflict")
        else:
            changed_signature = ("different", *value[4][1:])
        state.retain_product(product, changed_identity, changed_admission, changed_signature)
        state.retain_product(product, value[3], admission, value[4])  # No last-write recovery.
    if problem not in {"missing_product"} and not problem.startswith("conflicting"):
        state.products[product] = tuple(value)
    block_export_chemistry(monkeypatch)
    result = state.comparison(row)
    assert result["status"] == "unavailable" and result["reason_code"] == code
    assert result["original_hit_smiles"] == state.hit and result["final_product_smiles"] == product
    assert not RECORDED_ONLY.intersection(result)


def test_duplicate_rows_have_independent_serialized_values(monkeypatch):
    state, row, admission = retained()
    identity_value = state.products[row["smiles"]][3]
    state.retain_product(row["smiles"], identity_value, admission, workflow._identity_signature())
    other = {**row, "design": "other design", "route_id": "other-primary"}
    block_export_chemistry(monkeypatch)
    state.attach_selected([row, other])
    assert row["hit_identity_comparison"] == other["hit_identity_comparison"]
    row["hit_identity_comparison"]["original_hit"]["exact_smiles"] = "mutated"
    assert other["hit_identity_comparison"]["original_hit"]["exact_smiles"] == state.hit
    assert len(state.products) == 1 and state.snapshot()["recorded_selected_count"] == 2


def test_public_constraint_behavior_and_disabled_retention(monkeypatch):
    assert str(inspect.signature(workflow.check_constraints)) == (
        "(product, hit, design, mode='hop', core_smiles=None, preserve_design_scaffold=False)"
    )
    monkeypatch.setattr(workflow, "chemical_identity", forbidden)
    assert workflow.check_constraints(PRODUCT, HIT, PRODUCT) is None
    assert workflow.check_constraints("CCC(=O)Nc1ccncc1", HIT, PRODUCT, mode="expand") is None
    for product, hit, kwargs, message in (
        (HIT, HIT, {"core_smiles": "invalid"}, "original hit"),
        ("CC(=O)Nc1cc[nH+]cc1", PRODUCT, {}, "only chemical state"),
        ("CCC(=O)Nc1ccccc1", HIT, {}, "did not change"),
        (PRODUCT, HIT, {"core_smiles": "S"}, "core/anchor"),
    ):
        with pytest.raises(ValueError, match=message):
            workflow.check_constraints(product, hit, PRODUCT, **kwargs)


def test_optional_capture_once_and_deadline_skip(monkeypatch):
    state = workflow._SelectedIdentityEvidence(canonical(HIT), 10)
    monkeypatch.setattr(workflow.time, "monotonic", lambda: 10)
    monkeypatch.setattr(workflow, "chemical_identity", forbidden)
    state.capture_hit()
    monkeypatch.setattr(workflow.time, "monotonic", lambda: 1)
    state.capture_hit()
    assert state.lookup_count == 0
    assert state.hit_problem["reason_code"] == "deadline_before_hit_capture"


def test_optional_capture_reports_late_completion_and_observed_miss(monkeypatch):
    identity.clear_identity_caches()
    state = workflow._SelectedIdentityEvidence(canonical(HIT), 10)
    ticks = iter((9, 9, 11))
    monkeypatch.setattr(workflow.time, "monotonic", lambda: next(ticks))
    state.capture_hit()
    state.capture_hit()
    assert state.lookup_count == 1 and state.completed_after_deadline
    assert state.capture_seconds == 2 and state.cache_miss_delta >= 1
    assert "not a measured process-wide" in state.snapshot()["post_selection_identity_resolution_scope"]


@pytest.mark.parametrize(
    "error", [ValueError("expected normalization failure"), RuntimeError("bug"), KeyboardInterrupt()]
)
def test_only_optional_value_error_is_converted_to_unavailable(monkeypatch, error):
    state = workflow._SelectedIdentityEvidence(canonical(HIT), None)

    def fail(value):
        raise error

    monkeypatch.setattr(workflow, "chemical_identity", fail)
    if isinstance(error, ValueError):
        state.capture_hit()
        result = state.comparison({"smiles": PRODUCT})
        assert result["reason_code"] == "hit_capture_failed" and result["error_type"] == "ValueError"
        assert not RECORDED_ONLY.intersection(result)
    else:
        with pytest.raises(type(error)):
            state.capture_hit()


def config_for(tmp_path, assets, proposal):
    path = tmp_path / "proposals.json"
    path.write_text(
        json.dumps({"schema_version": 1, "proposals": [{"design": PRODUCT, "response": proposal}]})
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "search": {"building_block_top_k": 2, "beam_width": 3, "max_combinations": 12},
        "expansion": {"rounds": 1, "max_routes": 4},
        "refinement": {"top_n": 0},
        "max_scored_products": 3,
    }


def test_workflow_retains_before_selection_and_exports_same_api_json(tmp_path, assets, proposal, monkeypatch):
    config = config_for(tmp_path, assets, proposal)
    original_attach = workflow._SelectedIdentityEvidence.attach_selected
    observed = []

    def attach(state, selected):
        observed.append(state)
        # Context scope leaves the rest of the existing workflow unaffected.
        with monkeypatch.context() as patch:
            block_export_chemistry(patch)
            original_attach(state, selected)

    monkeypatch.setattr(workflow._SelectedIdentityEvidence, "attach_selected", attach)
    output = tmp_path / "out"
    rows = workflow.run(HIT, [PRODUCT], config, output, scoring_mode="2d")
    assert rows
    assert json.loads((output / "results.json").read_text())["results"] == rows
    state = observed[0]
    assert state.lookup_count == 1 and state.cache_miss_delta == 0
    assert len(state.products) <= config["max_scored_products"]
    manifest = json.loads((output / "manifest.json").read_text())
    descriptor = manifest["identity_comparison_policy"]
    assert descriptor == workflow._IDENTITY_COMPARISON_POLICY
    assert descriptor["code_sha256"] == {
        key: manifest["code_sha256"][key] for key in ("identity.py", "common.py")
    }
    assert (output / "results.tsv").read_text().splitlines()[0].split("\t") == workflow.COLUMNS
    for row in rows:
        comparison = row["hit_identity_comparison"]
        assert comparison["status"] == "recorded"
        assert comparison["original_hit_smiles"] == manifest["hit"]
        assert comparison["final_product_smiles"] == row["smiles"]
        assert comparison["final_product"]["family_smiles"] == row["family_smiles"]
        assert comparison["final_product"]["scaffold_family_smiles"] == row["scaffold_family_smiles"]
        assert comparison["original_hit"]["rdkit_version"] == descriptor["rdkit_version"]


def test_capture_failure_leaves_selection_scores_and_tsv_unchanged(tmp_path, assets, proposal, monkeypatch):
    config = config_for(tmp_path, assets, proposal)
    first = workflow.run(HIT, [PRODUCT], config, tmp_path / "recorded", scoring_mode="2d")
    capture = workflow._SelectedIdentityEvidence.capture_hit

    def failed_optional(state):
        with monkeypatch.context() as patch:
            patch.setattr(
                workflow, "chemical_identity", lambda value: (_ for _ in ()).throw(ValueError("fixture"))
            )
            capture(state)

    monkeypatch.setattr(workflow._SelectedIdentityEvidence, "capture_hit", failed_optional)
    second = workflow.run(HIT, [PRODUCT], config, tmp_path / "unavailable", scoring_mode="2d")
    assert first and second
    assert [{k: v for k, v in r.items() if k != "hit_identity_comparison"} for r in first] == [
        {k: v for k, v in r.items() if k != "hit_identity_comparison"} for r in second
    ]
    assert (tmp_path / "recorded/results.tsv").read_bytes() == (
        tmp_path / "unavailable/results.tsv"
    ).read_bytes()
    assert all(row["hit_identity_comparison"]["reason_code"] == "hit_capture_failed" for row in second)


def test_expand_never_captures_original_hit_for_export(tmp_path, assets, proposal, monkeypatch):
    config = config_for(tmp_path, assets, proposal)
    config["mode"] = "expand"
    monkeypatch.setattr(workflow._SelectedIdentityEvidence, "capture_hit", forbidden)
    rows = workflow.run(HIT, [PRODUCT], config, tmp_path / "expand", scoring_mode="2d")
    assert rows
    assert all(row["hit_identity_comparison"]["reason_code"] == "original_hit_not_retained" for row in rows)


def test_mandatory_product_family_failure_remains_error(tmp_path, assets, proposal, monkeypatch):
    config = config_for(tmp_path, assets, proposal)
    real = workflow.chemical_identity

    def fail_product(value):
        if value != canonical(HIT):
            raise ValueError("mandatory product family failure")
        return real(value)

    monkeypatch.setattr(workflow, "chemical_identity", fail_product)
    with pytest.raises(ValueError, match="mandatory product family failure"):
        workflow.run(HIT, [PRODUCT], config, tmp_path / "mandatory", scoring_mode="2d")


def test_root_search_checks_do_not_capture_hit(tmp_path, assets, proposal, monkeypatch):
    config = config_for(tmp_path, assets, proposal)
    config["route_search_policy"] = "original-hit-v1"
    assess = workflow.RootProductObjective.assess
    capture = workflow._SelectedIdentityEvidence.capture_hit
    inside_search, visits, captures = False, [], []

    def tracked_assess(objective, product, report=None):
        nonlocal inside_search
        visits.append(product)
        inside_search = True
        try:
            return assess(objective, product, report)
        finally:
            inside_search = False

    def tracked_capture(state):
        assert not inside_search
        capture(state)
        captures.append(state.lookup_count)

    monkeypatch.setattr(workflow.RootProductObjective, "assess", tracked_assess)
    monkeypatch.setattr(workflow._SelectedIdentityEvidence, "capture_hit", tracked_capture)
    rows = workflow.run(HIT, [PRODUCT], config, tmp_path / "guided", scoring_mode="2d")
    assert rows and visits and captures
    assert set(captures) == {1}


def test_continuation_captures_current_run_instead_of_saved_comparison(
    tmp_path, assets, proposal, monkeypatch
):
    config = config_for(tmp_path, assets, proposal)
    first = tmp_path / "first"
    rows = workflow.run(HIT, [PRODUCT], config, first, scoring_mode="2d")
    assert rows
    artifact = json.loads((first / "results.json").read_text())
    for row in artifact["results"]:
        row["hit_identity_comparison"] = {"status": "recorded", "original_hit": "bogus saved value"}
    (first / "results.json").write_text(json.dumps(artifact))
    config["carry_forward"] = [str(first)]
    captured = []
    capture = workflow._SelectedIdentityEvidence.capture_hit

    def fresh_capture(state):
        capture(state)
        captured.append(state)

    monkeypatch.setattr(workflow._SelectedIdentityEvidence, "capture_hit", fresh_capture)
    second = workflow.run(HIT, [PRODUCT], config, tmp_path / "second", scoring_mode="2d")
    assert captured and second
    for row in second:
        result = row["hit_identity_comparison"]
        assert result["status"] == "recorded" and result["original_hit"]["exact_smiles"] == HIT


def test_empty_selection_has_descriptor_and_no_hit_capture(tmp_path, assets, monkeypatch):
    config = config_for(tmp_path, assets, {"reactions": [], "building_blocks": []})
    monkeypatch.setattr(workflow._SelectedIdentityEvidence, "capture_hit", forbidden)
    output = tmp_path / "empty"
    assert workflow.run(HIT, [PRODUCT], config, output, scoring_mode="2d") == []
    assert json.loads((output / "manifest.json").read_text())["identity_comparison_policy"]
    summary = json.loads((output / "summary.json").read_text())["identity_comparison_export"]
    assert summary["optional_hit_capture_lookup_count"] == 0
    assert summary["recorded_selected_count"] == 0
