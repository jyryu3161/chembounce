"""Bounded refill at workflow, continuation and real integrated CLI boundaries."""

import itertools
import json

import pytest
from rdkit import Chem

from chembounce_synthesis import conformers, workflow
from chembounce_synthesis.cli import main
from chembounce_synthesis.continuation import continuation_inputs

HIT, DESIGN = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"
UNSPECIFIED = "CC(O)c1ccncc1"


def config_for(tmp_path, assets, proposal):
    path = tmp_path / "saved.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "proposals": [
                    {"design": DESIGN, "response": proposal, "sampling": {"seed": 17, "sample": 0}}
                ],
            }
        )
    )
    return {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {"kind": "saved", "path": str(path)},
        "search": {"building_block_top_k": 1, "beam_width": 1},
        "conformers": {"num_conformers": 2},
        "max_designs": 1,
        "max_scored_products": 1,
        "recommendation": {"mode": "qualified"},
        "minimum_scores": {"hit_tanimoto": 0.01, "hit_shape": 0.01, "hit_pharmacophore": 0.01},
    }


def read(out, name):
    return json.loads((out / name).read_text())


def forbid_assets(monkeypatch):
    monkeypatch.setattr(conformers, "rdShapeAlign", object())
    for name in ("Catalog", "Rules", "_provider", "_evidence_evaluator", "score_candidate"):
        monkeypatch.setattr(workflow, name, lambda *a, **k: pytest.fail("No eligible design may load assets"))


@pytest.mark.parametrize("limit,expected", [(3, False), (4, True)])
def test_bad_first_designs_refill_with_complete_denominators(tmp_path, assets, proposal, limit, expected):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)
    config["max_design_attempts"] = limit
    out = tmp_path / "refill"
    rows = workflow.run(HIT, ["invalid", UNSPECIFIED, "OC(C)c1ccncc1", DESIGN], config, out)
    assert bool(rows) is expected
    preflight = read(out, "manifest.json")["design_preflight"]
    assert preflight["counts"] == {
        "raw_supplied": 4,
        "raw_consumed": limit,
        "parse_attempted": limit,
        "preflight_attempted": 2 if expected else 1,
        "invalid": 1,
        "duplicates": 1,
        "preflight_rejected": 1,
        "eligible": int(expected),
        "unattempted": 4 - limit,
    }
    assert preflight["checks"][2]["duplicate_of_raw_index"] == 1
    assert read(out, "summary.json")["scored_products"] == int(expected)
    assert preflight["policy"]["version"] == "bounded-refill-v1"


@pytest.mark.parametrize("kind", ["hostile", "infinite", "finite"])
def test_unsized_input_is_never_read_past_raw_bound(tmp_path, assets, proposal, monkeypatch, kind):
    config = config_for(tmp_path, assets, proposal)
    config["max_design_attempts"] = 3
    forbid_assets(monkeypatch)

    def hostile():
        yield from [UNSPECIFIED] * 3
        pytest.fail("Read beyond max_design_attempts")

    designs = (
        hostile()
        if kind == "hostile"
        else itertools.repeat(UNSPECIFIED)
        if kind == "infinite"
        else iter([UNSPECIFIED])
    )
    out = tmp_path / kind
    assert workflow.run(HIT, designs, config, out) == []
    selection = read(out, "manifest.json")["design_preflight"]
    consumed = 1 if kind == "finite" else 3
    assert selection["counts"]["raw_consumed"] == consumed
    assert selection["counts"]["duplicates"] == consumed - 1
    assert selection["counts"]["raw_supplied"] is None
    assert selection["counts"]["unattempted"] is None
    assert selection["input_exhausted"] is (True if kind == "finite" else None)
    assert len(selection["checks"]) == consumed


def test_eligible_quota_does_not_parse_or_consume_tail(tmp_path, assets, proposal):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)

    def designs():
        yield DESIGN
        pytest.fail("Read past eligible quota")

    out = tmp_path / "quota"
    assert workflow.run(HIT, designs(), config, out)
    selection = read(out, "manifest.json")["design_preflight"]
    assert selection["counts"]["raw_consumed"] == 1
    assert selection["stop_reason"] == "eligible_limit" and selection["input_exhausted"] is None
    assert selection["policy"]["max_raw_attempts"] == 4


def test_simultaneous_limits_and_sized_unattempted_count(tmp_path, assets, proposal):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)
    config["max_design_attempts"] = 1
    out = tmp_path / "both"
    assert workflow.run(HIT, [DESIGN, "invalid"], config, out)
    selection = read(out, "manifest.json")["design_preflight"]
    assert selection["counts"]["unattempted"] == 1
    assert selection["limits_reached"] == ["max_design_attempts", "max_designs"]
    assert selection["checks"][0]["status"] == "eligible"


@pytest.mark.parametrize("expired_initially", [False, True])
def test_deadline_stops_preflight_and_preserves_unassessed_denominator(
    tmp_path, assets, proposal, monkeypatch, expired_initially
):
    config = config_for(tmp_path, assets, proposal)
    forbid_assets(monkeypatch)
    now = [2.0 if expired_initially else 0.0]
    monkeypatch.setattr(workflow.time, "monotonic", lambda: now[0])
    original = workflow.prepare_reference
    calls = []

    def prepare(smiles, *args):
        calls.append(smiles)
        if smiles == UNSPECIFIED:
            now[0] = 2.0
        return original(smiles, *args)

    monkeypatch.setattr(workflow, "prepare_reference", prepare)
    out = tmp_path / "deadline"
    assert workflow.run(HIT, [UNSPECIFIED, DESIGN], config, out, deadline_monotonic=1.0) == []
    selection = read(out, "manifest.json")["design_preflight"]
    assert selection["counts"]["raw_consumed"] == (0 if expired_initially else 1)
    assert selection["counts"]["unattempted"] == (2 if expired_initially else 1)
    assert selection["stop_reason"] == "deadline"
    assert calls == ([] if expired_initially else [HIT, UNSPECIFIED])
    assert read(out, "summary.json")["deadline_reached"]


def test_all_parse_failures_are_recorded_without_loading_assets(tmp_path, assets, proposal, monkeypatch):
    config = config_for(tmp_path, assets, proposal)
    forbid_assets(monkeypatch)
    out = tmp_path / "invalid"
    assert workflow.run(HIT, [None, "invalid", "CC.O", "*C"], config, out) == []
    selection = read(out, "manifest.json")["design_preflight"]
    assert selection["counts"]["invalid"] == selection["counts"]["raw_consumed"] == 4
    assert selection["counts"]["preflight_attempted"] == 0
    assert all(check["reason_code"] == "invalid_design" for check in selection["checks"])
    assert read(out, "manifest.json")["provider"] is None


def test_unspecified_query_is_explicitly_unassessed_and_never_repaired(
    tmp_path, assets, proposal, monkeypatch
):
    config = config_for(tmp_path, assets, proposal)
    forbid_assets(monkeypatch)
    out = tmp_path / "query"
    with pytest.raises(ValueError, match="specified stereochemistry"):
        workflow.run(UNSPECIFIED, [DESIGN], config, out)
    failure = read(out, "summary.json")["query_preflight"]
    assert failure["input"] == failure["canonical_smiles"] == UNSPECIFIED
    assert failure["reason_code"] == "unspecified_stereochemistry"
    assert failure["qualification_status"] == "not_assessed"


def test_rejected_prefix_preserves_provider_sampling_and_scoring_budget(tmp_path, assets, proposal):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)
    direct, refill = tmp_path / "direct", tmp_path / "refill"
    assert workflow.run(HIT, [DESIGN], config, direct) == workflow.run(
        HIT, ["invalid", UNSPECIFIED, DESIGN], config, refill
    )
    assert read(direct, "proposals.json") == read(refill, "proposals.json")
    assert read(direct, "scores.json") == read(refill, "scores.json")
    assert read(refill, "summary.json")["scored_products"] == 1


def test_exact_known_stereo_states_remain_distinct_and_unchanged(tmp_path, assets, proposal):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)
    config.update(max_designs=2, max_scored_products=2)
    designs = ["CC(=O)N[C@@H](C)C(=O)O", "CC(=O)N[C@H](C)C(=O)O"]
    records = []
    for design, amine in zip(designs, ["N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"]):
        records.append(
            {
                "design": design,
                "response": {
                    "reactions": [
                        {"reaction_template": "amide", "reactants": ["CC(=O)O", amine], "product": design}
                    ],
                    "building_blocks": ["CC(=O)O", amine],
                },
            }
        )
    (tmp_path / "saved.json").write_text(json.dumps({"schema_version": 1, "proposals": records}))
    out = tmp_path / "stereo"
    hit = "C[C@H](O)c1ccccc1"
    rows = workflow.run(hit, designs, config, out)
    assert rows
    selection = read(out, "manifest.json")["design_preflight"]
    assert selection["eligible_designs"] == designs and selection["counts"]["duplicates"] == 0
    assert read(out, "manifest.json")["hit"] == hit
    assert read(out, "summary.json")["query_preflight"]["input"] == hit
    assert all(r["verification"]["product"] == r["design"] for r in read(out, "routes.json")["routes"])
    pose = next(iter(Chem.SDMolSupplier(str(out / "poses_hit.sdf"))))
    assert Chem.MolToSmiles(pose) == rows[0]["smiles"] in designs


@pytest.mark.parametrize("with_saved_other", [False, True])
def test_continuation_preserves_rejected_history_and_unselected_route_then_recovers(
    tmp_path, assets, proposal, with_saved_other
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)
    first = tmp_path / "first"
    assert workflow.run(HIT, [UNSPECIFIED, "invalid", UNSPECIFIED, DESIGN], config, first)
    hit, designs, resumed = continuation_inputs(first)
    assert designs == [UNSPECIFIED, "invalid", UNSPECIFIED, DESIGN]
    resumed["max_design_attempts"] = 2
    second = tmp_path / "second"
    assert workflow.run(hit, designs, resumed, second) == []
    records = read(second, "routes.json")["routes"]
    assert len(records) == 1 and records[0]["selection_status"] == "design_not_assessed"
    assert records[0]["route"] == read(first, "routes.json")["routes"][0]["route"]
    assert (
        read(second, "manifest.json")["carry_forward"][0]["prior_design_preflight"]["counts"]["duplicates"]
        == 1
    )
    hit, designs, resumed = continuation_inputs(second)
    assert designs == [UNSPECIFIED, "invalid", UNSPECIFIED, DESIGN]
    resumed["max_design_attempts"] = 4
    if with_saved_other:
        other = "CCC(=O)Nc1ccncc1"
        designs.append(other)
        resumed.update(max_designs=2, max_design_attempts=5)
        partial = tmp_path / "partial-proposals.json"
        partial.write_text(
            json.dumps({"schema_version": 1, "proposals": [{"design": other, "response": proposal}]})
        )
        resumed["provider"] = {"kind": "saved", "path": str(partial)}
    third = tmp_path / "third"
    rows = workflow.run(hit, designs, resumed, third)
    assert rows and rows[0]["recommendation_status"] == "qualified"
    summary = read(third, "summary.json")
    assert summary["carried_screened_products"] == summary["scored_products"] == 1
    assert summary["proposals"] == int(with_saved_other) and summary["carried_only_designs"] == [DESIGN]
    current_route = read(third, "routes.json")["routes"][0]
    assert current_route["route_id"] == rows[0]["route_id"]
    assert (
        current_route["verification"]["product"] == current_route["synthesis_evidence"]["product"] == DESIGN
    )
    assert read(third, "manifest.json")["effective_workflow"]["max_scored_products"] == 1
    with pytest.raises(ValueError, match="different original hit"):
        workflow.run("CCC(=O)Nc1ccccc1", designs, resumed, tmp_path / "wrong-hit")


@pytest.mark.parametrize("attempts", [0, -1, True, 1.5])
def test_invalid_attempt_budget_is_configuration_error(tmp_path, assets, proposal, attempts):
    config = config_for(tmp_path, assets, proposal)
    config["max_design_attempts"] = attempts
    with pytest.raises(ValueError, match="max_design_attempts"):
        workflow.run(HIT, [DESIGN], config, tmp_path / "bad-config")


def test_exploratory_refill_requires_explicit_attempt_budget(tmp_path, assets, proposal):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)
    config.pop("recommendation")
    assert workflow.run(HIT, [UNSPECIFIED, DESIGN], config, tmp_path / "legacy") == []
    config["max_design_attempts"] = 2
    rows = workflow.run(HIT, [UNSPECIFIED, DESIGN], config, tmp_path / "refill")
    assert rows and rows[0]["recommendation_status"] == "exploratory_verified"


def test_standalone_tsv_intake_is_bounded_before_reading_tail(tmp_path, assets, proposal):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    config = config_for(tmp_path, assets, proposal)
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    tsv = tmp_path / "designs.tsv"
    # A CSV field beyond the intake bound would exceed csv's default field limit.
    tsv.write_text("smiles\ninvalid\n" + UNSPECIFIED + "\n" + DESIGN + "\n" + "C" * 200000 + "\n")
    out = tmp_path / "tsv"
    assert (
        main(["--hit", HIT, "--designs-tsv", str(tsv), "--config", str(config_path), "--output", str(out)])
        == 0
    )
    assert read(out, "manifest.json")["design_preflight"]["counts"]["raw_consumed"] == 3


def test_duplicate_unselected_carried_graph_is_counted_once_across_sources(
    tmp_path, assets, proposal, route, monkeypatch
):
    config = config_for(tmp_path, assets, proposal)
    sources = [tmp_path / name for name in ("source-a", "source-b")]
    for source in sources:
        source.mkdir()
        (source / "manifest.json").write_text(
            json.dumps({"schema_version": 1, "hit": HIT, "designs": [DESIGN], "config": config})
        )
        (source / "routes.json").write_text(
            json.dumps(
                {"schema_version": 1, "routes": [{"route_id": "prior", "design": DESIGN, "route": route}]}
            )
        )
    config["carry_forward"] = [str(source) for source in sources]
    forbid_assets(monkeypatch)
    out = tmp_path / "deduplicated"
    assert workflow.run(HIT, [UNSPECIFIED], config, out) == []
    provenance = read(out, "manifest.json")["carry_forward"]
    assert [source["loaded_routes"] for source in provenance] == [1, 0]
    assert [source["duplicate_routes"] for source in provenance] == [0, 1]
    assert [source["unselected_design_routes_retained"] for source in provenance] == [1, 0]
    retained = read(out, "routes.json")["routes"]
    assert len(retained) == read(out, "summary.json")["carried_routes_loaded"] == 1
    assert retained[0]["route"] == route and retained[0]["selection_status"] == "design_not_assessed"


@pytest.mark.parametrize("explicit_cap", [None, 1])
def test_real_integrated_cli_supplies_refill_pool_and_honors_explicit_cap(
    tmp_path, assets, proposal, explicit_cap
):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    import chembounce

    config = config_for(tmp_path, assets, proposal)
    config["max_design_attempts"] = 2
    proposals = []
    for amine in ("Nc1ccncc1", "Nc1cccnc1", "Nc1ccccn1"):
        design = "CC(=O)" + amine
        proposals.append(
            {
                "design": design,
                "response": {
                    "reactions": [
                        {"reaction_template": "amide", "reactants": ["CC(=O)O", amine], "product": design}
                    ],
                    "building_blocks": ["CC(=O)O", amine],
                },
            }
        )
    (tmp_path / "saved.json").write_text(json.dumps({"schema_version": 1, "proposals": proposals}))
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    replacements = tmp_path / "replacements.smi"
    replacements.write_text("CC(O)C1CCCCC1\nc1ccncc1\n")
    out = tmp_path / "integrated"
    args = [
        "-i",
        HIT,
        "-o",
        str(out),
        "--fragments",
        "c1ccccc1",
        "--replace_scaffold_files",
        str(replacements),
        "--cand_max_n__rplc",
        "1",
        "--no-shape",
        "-t",
        "0",
        "--search-threshold",
        "0",
        "--synthesis-config",
        str(config_path),
    ]
    if explicit_cap is not None:
        args += ["--overall_max_n", str(explicit_cap)]
    assert chembounce.main(args) == (3 if explicit_cap else 0)
    supply = read(out, "run_manifest.json")["synthesis_design_supply"]
    assert supply["effective_overall_max_n"] == (explicit_cap or 2)
    assert supply["explicit_overall_cap_limits_pool"] is bool(explicit_cap)
    counts = read(out / "synthesis", "manifest.json")["design_preflight"]["counts"]
    assert read(out / "synthesis", "manifest.json")["route_search_policy"]["version"] == "original-hit-v1"
    assert counts["preflight_rejected"] == 1 and counts["eligible"] == (0 if explicit_cap else 1)
    assert read(out / "synthesis", "summary.json")["scored_products"] <= 1
