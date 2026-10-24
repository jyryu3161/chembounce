"""Saved-only unselected route review; fixtures do not construct molecules."""

import copy
import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest

from chembounce_synthesis import report
from chembounce_synthesis import report_routes

_spec = importlib.util.spec_from_file_location("saved_report_test_fixture", Path(__file__).with_name("test_report.py"))
_fixture = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fixture)
PRODUCT, bundle, save = _fixture.PRODUCT, _fixture.bundle, _fixture.save


def bound(route):
    evidence = route["synthesis_evidence"]
    audit = evidence["route_audit"]
    evidence["predictor_provenance"] = {"model": "saved-fixture", "top_k": 3}
    parents = {}
    for block in audit["building_blocks"]:
        block.pop("building_block_binding_sha256", None)
        block["building_block_binding_sha256"] = report._binding(block)
        parents[block["node_id"]] = block["building_block_binding_sha256"]
    for audited, supported in zip(audit["steps"], evidence["steps"]):
        for slot in audited["reactant_slots"]:
            slot["parent_binding_sha256"] = parents[slot["node_id"]]
        audited.pop("step_binding_sha256", None)
        audited["step_binding_sha256"] = report._binding(audited)
        parents[audited["node_id"]] = audited["step_binding_sha256"]
        supported["step_binding_sha256"] = audited["step_binding_sha256"]
        raw = supported["forward_prediction"].get("raw")
        if raw is not None:
            raw["input"] = {
                "reactants": sorted(supported["reactants"]), "reagents": [],
                "prompt": "REACTANT:" + ".".join(sorted(supported["reactants"])) + "REAGENT:",
            }
            raw["predictor_binding_sha256"] = report._binding(evidence["predictor_provenance"])
            supported["context"] = {"status": "not_supplied", "supplied": {}}
    audit.pop("route_binding_sha256", None)
    audit["route_binding_sha256"] = report._binding(audit)
    evidence["route_binding_sha256"] = audit["route_binding_sha256"]
    return route


def other_route():
    route = copy.deepcopy(bundle()["routes.json"]["routes"][0])
    route["route_id"] = "rejected-route"
    route["route"]["schema_version"] = 1
    route["selection_status"] = "evidence_rejected"
    route.pop("qualification", None)
    evidence = route["synthesis_evidence"]
    evidence["route_audit"]["schema"] = "concrete-route-audit-v1"
    leaf = evidence["route_audit"]["building_blocks"][0]
    route["verification"]["building_blocks"] = [
        {"node_id": leaf["node_id"], "id": leaf["catalog_record"]["id"], "smiles": leaf["smiles"]}
    ]
    evidence["eligibility"] = {"passed": False, "reasons": ["bounded prediction mismatch"]}
    prediction = evidence["steps"][0]["forward_prediction"]
    prediction.update(
        status="mismatched",
        matched_rank=None,
        raw={
            "predictions": [
                {"rank": rank, "complete": True, "product": f"other-{rank}", "original_decoded_text": f"<raw-{rank}>"}
                for rank in (1, 2, 3)
            ]
        },
    )
    return bound(route)


def two_step_route():
    route = other_route()
    nodes = route["route"]["nodes"]
    nodes[-1]["reactants"] = ["ir", "bb1"]
    nodes.insert(1, {"id": "ir", "rule_id": "inner-rule", "reactants": ["bb1"], "product": "intermediate"})
    verified = route["verification"]["steps"]
    verified[-1]["reactants"] = ["intermediate"] + verified[-1]["reactants"]
    verified.insert(
        0,
        {"node_id": "ir", "rule_id": "inner-rule", "reactants": ["CC(=O)O"], "product": "intermediate"},
    )
    audit = route["synthesis_evidence"]["route_audit"]
    audit["steps"][-1]["reactant_slots"] = [
        {"slot": 0, "node_id": "ir", "smiles": "intermediate"},
        {"slot": 1, "node_id": "bb1", "smiles": "CC(=O)O"},
    ]
    audit["steps"].insert(
        0,
        {
            "node_id": "ir", "rule_id": "inner-rule", "product": "intermediate",
            "reactant_slots": [{"slot": 0, "node_id": "bb1", "smiles": "CC(=O)O"}],
            "step_binding_sha256": "inner-binding",
        },
    )
    evidence_steps = route["synthesis_evidence"]["steps"]
    evidence_steps[-1]["reactants"] = ["intermediate", "CC(=O)O"]
    evidence_steps.insert(
        0,
        {
            "node_id": "ir", "reactants": ["CC(=O)O"], "product": "intermediate",
            "step_binding_sha256": "inner-binding", "rule_replay": "matched",
            "forward_prediction": {"status": "unknown", "reason": "budget exhausted"},
            "context": {"status": "not_supplied"},
        },
    )
    return bound(route)


def test_rejected_two_step_keeps_graph_and_partial_prediction_separate():
    route = two_step_route()
    review = report_routes.build([route], [])
    item = review["routes"][0]
    assert item["identity_status"] == "consistent_graph_and_verification"
    assert item["evidence_status"] == "recorded_bound_steps"
    assert item["forward_review_status"] == "contains_recorded_nonmatch"
    assert [step["node_id"] for step in item["steps"]] == ["ir", "rxn"]
    assert item["steps"][0]["forward_status"] == "unassessed"
    assert item["steps"][1]["hypotheses"][0]["decoded_text"] == "<raw-1>"
    assert item["intermediates"] == [{"node_id": "ir", "product": "intermediate"}]


def test_missing_evidence_still_shows_verified_expected_steps():
    route = two_step_route()
    route.pop("synthesis_evidence")
    item = report_routes.build([route], [])["routes"][0]
    assert item["identity_status"] == "consistent_graph_and_verification"
    assert item["evidence_status"] == "not_recorded"
    assert [(step["node_id"], step["expected_product"]) for step in item["steps"]] == [
        ("ir", "intermediate"), ("rxn", PRODUCT),
    ]
    assert all(step["forward_status"] == "unassessed" for step in item["steps"])


@pytest.mark.parametrize("missing_index", [0, 1])
def test_missing_two_step_support_sibling_retains_bound_raw_other_step(missing_index):
    route = two_step_route()
    if missing_index == 1:
        prediction = route["synthesis_evidence"]["steps"][0]["forward_prediction"]
        prediction.update(
            status="matched", matching_policy="exact", matched_rank=1,
            raw={"predictions": [
                {"rank": 1, "complete": True, "product": "intermediate",
                 "original_decoded_text": "<raw-inner>"},
                {"rank": 2, "complete": True, "product": None, "invalid_reason": "invalid_saved_smiles"},
                {"rank": 3, "complete": True, "product": None, "invalid_reason": "invalid_saved_smiles"},
            ]},
        )
        bound(route)
    route["synthesis_evidence"]["steps"].pop(missing_index)
    item = report_routes.build([route], [])["routes"][0]
    assert item["identity_status"] == "consistent_graph_and_verification"
    assert item["evidence_status"] == "incomplete_or_inconsistent"
    assert item["forward_review_status"] == "incomplete_or_unknown"
    assert item["steps"][0]["expected_product"] == "intermediate"
    assert item["steps"][missing_index]["forward_status"] == "unassessed"
    surviving = item["steps"][1 - missing_index]
    assert surviving["hypotheses"][0]["decoded_text"] == ("<raw-1>" if missing_index == 0 else "<raw-inner>")
    assert surviving["forward_status"] == (
        "recorded_complete_nonmatch" if missing_index == 0 else "recorded_exact_match"
    )


@pytest.mark.parametrize(
    "break_record",
    [
        lambda r: r["route"]["nodes"][-1].update(product="wrong"),
        lambda r: r["verification"]["steps"][0].update(node_id="wrong"),
        lambda r: r["verification"]["steps"][0].update(reactants=["wrong"]),
        lambda r: r["synthesis_evidence"]["steps"][0].update(reactants=["wrong"]),
        lambda r: r["synthesis_evidence"]["steps"][0].update(product="wrong"),
        lambda r: r["synthesis_evidence"]["steps"][0].update(step_binding_sha256="wrong"),
        lambda r: r["synthesis_evidence"].update(route_binding_sha256="wrong"),
        lambda r: r["synthesis_evidence"]["route_audit"]["building_blocks"][0]["catalog_record"].update(id="wrong"),
    ],
)
def test_inconsistent_joins_never_look_supported(break_record):
    route = other_route()
    break_record(route)
    item = report_routes.build([route], ["route-final"])["routes"][0]
    assert item["issues"]
    assert item["forward_review_status"] in ("unassessed", "incomplete_or_unknown")


def test_incomplete_or_corrupt_raw_remains_unknown():
    route = other_route()
    prediction = route["synthesis_evidence"]["steps"][0]["forward_prediction"]
    prediction["status"] = "partially_evaluated"
    prediction["raw"]["predictions"][-1] = {"rank": 3, "complete": False, "product": None}
    item = report_routes.build([route], [])["routes"][0]
    assert item["steps"][0]["forward_status"] == "unassessed"
    prediction["raw"]["predictions"][-1]["product"] = PRODUCT
    item = report_routes.build([route], [])["routes"][0]
    assert item["forward_review_status"] == "incomplete_or_unknown"
    assert item["issues"]


def test_invalid_complete_beam_does_not_erase_other_step_detail():
    route = two_step_route()
    prediction = route["synthesis_evidence"]["steps"][-1]["forward_prediction"]
    prediction["raw"]["predictions"][1].update(product=None, invalid_reason="invalid_saved_smiles")
    item = report_routes.build([route], [])["routes"][0]
    assert [step["node_id"] for step in item["steps"]] == ["ir", "rxn"]
    assert item["steps"][1]["forward_status"] == "unassessed"
    assert item["steps"][1]["hypotheses"][1]["invalid_reason"] == "invalid_saved_smiles"


def test_rehashed_wrong_parent_binding_is_still_rejected():
    route = other_route()
    audit = route["synthesis_evidence"]["route_audit"]
    audit["steps"][0]["reactant_slots"][0]["parent_binding_sha256"] = "0" * 64
    audit["steps"][0].pop("step_binding_sha256")
    audit["steps"][0]["step_binding_sha256"] = report._binding(audit["steps"][0])
    route["synthesis_evidence"]["steps"][0]["step_binding_sha256"] = audit["steps"][0]["step_binding_sha256"]
    audit.pop("route_binding_sha256")
    audit["route_binding_sha256"] = report._binding(audit)
    route["synthesis_evidence"]["route_binding_sha256"] = audit["route_binding_sha256"]
    item = report_routes.build([route], [])["routes"][0]
    assert item["forward_review_status"] == "unassessed"
    assert any("parent bindings" in issue for issue in item["issues"])


def test_short_or_conflicting_raw_beams_cannot_claim_complete_nonmatch():
    route = other_route()
    prediction = route["synthesis_evidence"]["steps"][0]["forward_prediction"]
    prediction["raw"]["predictions"].pop()
    item = report_routes.build([route], [])["routes"][0]
    assert item["steps"][0]["forward_status"] == "unassessed"
    assert item["evidence_status"] == "incomplete_or_inconsistent"
    prediction["raw"]["predictions"].append(
        {"rank": 3, "complete": True, "product": "other-3", "invalid_reason": "decoder_invalid"}
    )
    item = report_routes.build([route], [])["routes"][0]
    assert item["steps"][0]["forward_status"] == "unassessed"
    assert item["issues"]


@pytest.mark.parametrize("change", ["matched_rank", "prompt"])
def test_raw_request_or_rank_contradiction_remains_unknown(change):
    route = other_route()
    prediction = route["synthesis_evidence"]["steps"][0]["forward_prediction"]
    if change == "matched_rank":
        prediction["matched_rank"] = 1
    else:
        prediction["raw"]["input"]["prompt"] = "REACTANT:wrongREAGENT:"
    item = report_routes.build([route], [])["routes"][0]
    assert item["steps"][0]["forward_status"] == "unassessed"
    assert item["forward_review_status"] == "incomplete_or_unknown"
    assert item["evidence_status"] == "incomplete_or_inconsistent"
    assert item["issues"]


def test_float_support_schema_is_unknown():
    route = other_route()
    route["synthesis_evidence"]["schema_version"] = 2.0
    item = report_routes.build([route], [])["routes"][0]
    assert item["identity_status"] == "consistent_graph_and_verification"
    assert item["evidence_status"] == "incomplete_or_inconsistent"


def test_wrong_support_sibling_keeps_other_bound_step_and_global_consistency(tmp_path, capsys):
    data = bundle()
    route = two_step_route()
    route["synthesis_evidence"]["steps"][0]["product"] = "wrong"
    data["routes.json"]["routes"].append(route)
    view = report.build_report(save(tmp_path, data))
    item = view["unselected_routes"]["routes"][0]
    assert item["identity_status"] == "consistent_graph_and_verification"
    assert item["steps"][0]["forward_status"] == "unassessed"
    assert item["steps"][1]["forward_status"] == "recorded_complete_nonmatch"
    assert item["steps"][1]["hypotheses"]
    assert view["artifact_consistency"] == "consistent_recorded_bundle"
    assert view["unselected_artifact_consistency"] == "incomplete_or_inconsistent"
    rendered = report.render_html(view)
    assert "Selected recommendation artifact consistency" in rendered
    assert "Unselected route artifact consistency (displayed records only)" in rendered
    assert report.main([str(tmp_path)]) == 3
    stdout = capsys.readouterr().out
    assert "Selected recommendation artifact consistency: consistent_recorded_bundle" in stdout
    assert "unselected route artifact consistency (displayed records only): incomplete_or_inconsistent" in stdout


def test_selected_dedup_cap_and_html_escaping(tmp_path):
    data = bundle()
    routes = [copy.deepcopy(data["routes.json"]["routes"][0])]
    routes += [other_route() for _ in range(18)]
    for index, route in enumerate(routes[1:]):
        route["route_id"] = f"other-{index:02d}"
    routes[1]["reason"] = "<script>alert(1)</script>"
    data["routes.json"]["routes"] = routes
    view = report.build_report(save(tmp_path, data))
    review = view["unselected_routes"]
    assert (review["total"], review["displayed"], review["omitted"]) == (18, 16, 2)
    assert view["unselected_artifact_consistency"] == "consistent_displayed_records"
    assert review["routes"][0]["route_id"] == "other-00"
    assert all(item["route_id"] != "route-final" for item in review["routes"])
    rendered = report.render_html(view)
    assert "Other evaluated and unselected routes" in rendered
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in rendered
    assert "<script>alert(1)</script>" not in rendered


def test_zero_selected_still_shows_unselected(tmp_path):
    data = bundle()
    data["results.json"]["results"] = []
    data["qualification.json"]["selected_route_ids"] = []
    data["summary.json"].update(status="no_candidates", returned_products=0, accepted_products=0, qualified_products=0)
    data["routes.json"]["routes"] = [other_route()]
    data["routes.json"]["routes"][0]["design"] = "recorded-different-design"
    data["routes.json"]["routes"][0]["tanimoto"] = {"hit_tanimoto": 0.4, "design_tanimoto": 0.7}
    data["routes.json"]["routes"][0]["polar_capacity"] = {"recall_upper_bound": 0.8, "status": "possible"}
    view = report.build_report(save(tmp_path, data))
    assert view["unselected_routes"]["total"] == 1
    rendered = report.render_html(view)
    assert "No selected candidates" in rendered
    assert "rejected-route" in rendered
    assert _fixture.HIT in rendered
    assert "recorded-different-design" in rendered
    assert "polar_recall_upper_bound" in rendered


def test_stdlib_only_import():
    code = "import sys; from chembounce_synthesis import report, report_routes; " \
        "assert 'rdkit' not in sys.modules and 'torch' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], check=True)


def test_duplicate_unselected_id_is_diagnostic():
    routes = [other_route(), copy.deepcopy(other_route())]
    review = report_routes.build(routes, [])
    assert review["total"] == 2
    assert all(item["identity_status"] == "inconsistent_or_incomplete" for item in review["routes"])

