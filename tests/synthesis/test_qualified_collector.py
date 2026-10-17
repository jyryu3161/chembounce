"""Real toy chemistry and controlled raw beams; never local model/tokenizer loading."""

import json
from pathlib import Path

import pytest

# Test-only artifacts retain real loaders, policy checks and collector validation.
pytestmark = pytest.mark.usefixtures("synthetic_benchmark_profile_bindings")

from benchmarks.research import qualified_collector as collector
from benchmarks.research.protocol import digest
from chembounce_synthesis import workflow

SELECTIVE_PROFILES = (
    "qualified-selective-engineering-v2",
    "qualified-selective-engineering-v5",
    "qualified-selective-engineering-v6",
)

HIT, DESIGN = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"


def toy_request(tmp_path, assets, route, profile_name):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from rdkit.Chem import rdChemReactions
    from rdkit import Chem
    from benchmarks.research import qualified_runner as runner
    from benchmarks.research.protocol import write_json
    from chembounce_synthesis.profiles import profile_provenance
    from chembounce_synthesis.routes import verify_route

    rxn = rdChemReactions.ReactionFromSmarts("[C:1]-[O:2]>>[C:1]")
    if not rxn.RunReactants((Chem.MolFromSmiles("CO"),))[0][0].GetAtomWithIdx(0).HasProp("react_idx"):
        pytest.skip("Concrete replay metadata unavailable; strict support fails closed")
    config = {
        "schema_version": 1,
        "catalog": str(assets[0].path),
        "rules": str(assets[1].path),
        "provider": {
            "kind": "synllama",
            "model_path": "unused",
            "prompt_path": "unused",
            "samples": 2,
            "max_new_tokens": 1200,
        },
        "recommendation_profile": profile_name,
        "max_design_attempts": 12,
        "max_designs": 3,
        "max_results": 3,
        "max_scored_products": 12,
        "search": {
            "building_block_top_k": 3,
            "beam_width": 6,
            "max_combinations": 256,
            "max_steps": 12,
            "reaction_product_limit": 64,
        },
        "expansion": {"rounds": 1, "max_routes": 12, "building_block_top_k": 6, "min_tanimoto": 0.0},
        "refinement": {"top_n": 6, "max_total_conformers": 144},
        "synthesis_evidence": {
            "profile": "computational-route-support-v1",
            "predictor": {"kind": "reactiont5", "model_path": str(tmp_path / "controlled-model")},
            "max_predictor_calls": 24,
            "max_input_tokens": 12288,
            "max_output_tokens": 18432,
        },
    }
    path = tmp_path / "config.json"
    write_json(path, config)
    expected_status = (
        "engineering_policy_matched" if profile_name in SELECTIVE_PROFILES else "development_protocol_matched"
    )
    assert profile_provenance(workflow.load_config(path))["status"] == expected_status
    native, out = tmp_path / "native", tmp_path / "common"
    (native / "synthesis").mkdir(parents=True)
    out.mkdir()
    query = {"id": "fixture", "target_id": "synthetic", "smiles": HIT}
    write_json(native / "request.json", {"query": query, "seed": 97})
    write_json(
        native / "response.json", {"candidates": [{"smiles": DESIGN, "design": DESIGN, "route": route}]}
    )
    (native / "synthesis/results.tsv").write_text(
        "route_id\tdesign\tsmiles\nselected\t" + DESIGN + "\t" + DESIGN + "\n"
    )
    write_json(
        native / "synthesis/routes.json",
        {
            "routes": [
                {
                    "route_id": "selected",
                    "design": DESIGN,
                    "route": route,
                    "verification": verify_route(route, *assets),
                }
            ]
        },
    )
    request = {
        "query": query,
        "seed": 97,
        "settings": {"arm": "baseline"},
        "assets": {"synthesis_config": str(path)},
        "native": {"directory": str(native), "files_sha256": runner.file_manifest(native)},
        "output_dir": str(out),
        "deadline_monotonic": 999999999999,
        "finalization_margin_seconds": 2,
    }
    if profile_name in SELECTIVE_PROFILES:
        from chembounce_synthesis.profiles import PROFILE_SHA256

        request["benchmark_protocol"] = {
            "schema": "qualified-benchmark-protocol-"
            + dict(zip(SELECTIVE_PROFILES, ("v2", "v3", "v4")))[profile_name],
            "recommendation_profile": profile_name,
            "profile_artifact_sha256": PROFILE_SHA256[profile_name],
        }
    return request


@pytest.fixture(params=["qualified-development-v1", *SELECTIVE_PROFILES])
def completed(tmp_path, assets, route, strict_predictor, request):
    from benchmarks.research import qualified_adapter as adapter
    from benchmarks.research.protocol import write_json

    request = toy_request(tmp_path, assets, route, request.param)
    out = Path(request["output_dir"])
    write_json(out / "request.json", request)
    response = adapter.common(request)
    if "benchmark_protocol" in request:
        response["metadata"]["benchmark_protocol"] = request["benchmark_protocol"]
    assert len(response["candidates"]) == 1
    return out, response, strict_predictor


def collect(data, bound_catalog=None):
    out, response, decoder = data
    return collector.checked_products(
        out,
        response,
        {"id": "fixture", "target_id": "synthetic", "smiles": HIT},
        workflow._SOURCE_HASHES,
        decoder_factory=lambda *args: decoder,
        bound_catalog=bound_catalog,
        benchmark_protocol=response["metadata"].get("benchmark_protocol"),
    )


@pytest.fixture
def v6_full_request(tmp_path, assets, route, proposal, strict_predictor, monkeypatch):
    """Real common/collector path; controlled native supply, proposals and raw beams."""
    from benchmarks.research import qualified_runner as runner
    from benchmarks.research.protocol import write_json
    from chembounce_synthesis.synllama import SavedProposals

    request = toy_request(tmp_path, assets, route, "qualified-selective-engineering-v6")
    request["settings"] = {"arm": "full"}
    scaffold_db = tmp_path / "native-scaffolds.csv"
    scaffold_db.write_text("SMILES\nc1ccncc1\n")
    request["assets"]["scaffold_db"] = str(scaffold_db)
    native = Path(request["native"]["directory"])
    write_json(
        native / "request.json",
        {key: request[key] for key in ("query", "seed", "settings", "assets", "benchmark_protocol")},
    )
    write_json(native / "design-pool.json", {"designs": [DESIGN]})
    request["native"]["files_sha256"] = runner.file_manifest(native)
    proposals = tmp_path / "controlled-proposals.json"
    write_json(proposals, {"schema_version": 1, "proposals": [{"design": DESIGN, "response": proposal}]})
    monkeypatch.setattr(workflow, "SynLlama", lambda **kwargs: SavedProposals(proposals))
    return request, strict_predictor


@pytest.mark.parametrize("ablation", ["no_extra_refill_capacity", "design_only_route_search"])
@pytest.mark.parametrize("empty", [False, True], ids=["selected", "failed-final-seed"])
def test_v6_full_ablation_collector_roundtrip(v6_full_request, monkeypatch, ablation, empty):
    from benchmarks.research import qualified_adapter as adapter
    from benchmarks.research.protocol import write_json
    from chembounce_synthesis import selective_scoring

    request, predictor = v6_full_request
    base = adapter.common_config(request)
    request["settings"]["ablation"] = ablation
    config = adapter.common_config(request)
    key, value = (
        ("max_design_attempts", 3)
        if ablation == "no_extra_refill_capacity"
        else ("route_search_policy", "design-v1")
    )
    assert config[key] == value
    assert {k: v for k, v in config.items() if k != key} == {k: v for k, v in base.items() if k != key}
    if empty:
        original = selective_scoring.score_candidate

        def fail_seed(smiles, hit, design, settings, **kwargs):
            if settings.seed == 44:
                raise ValueError("controlled ablation final-seed abstention")
            return original(smiles, hit, design, settings, **kwargs)

        monkeypatch.setattr(selective_scoring, "score_candidate", fail_seed)
    out = Path(request["output_dir"])
    write_json(out / "request.json", request)
    response = adapter.common(request)
    response["metadata"]["benchmark_protocol"] = request["benchmark_protocol"]
    calls_before_collection = len(predictor.calls)
    rows = collect((out, response, predictor))
    assert calls_before_collection > 0 and len(predictor.calls) == calls_before_collection
    results = adapter.read(out / "synthesis/results.json")
    calibration = results["policy"]["calibration"]
    expected_status = (
        "custom_unvalidated" if ablation == "design_only_route_search" else "engineering_policy_matched"
    )
    assert calibration["status"] == expected_status
    assert calibration["material_overrides"] == (
        {"route_search_policy": {"profile": "original-hit-v1", "effective": "design-v1"}}
        if ablation == "design_only_route_search"
        else {}
    )
    if empty:
        assert rows == response["candidates"] == results["results"] == []
        assert adapter.read(out / "synthesis/summary.json")["status"] == "no_candidates"
        refined = adapter.read(out / "synthesis/refinement.json")["candidates"]
        assert refined and all(c["trials"][-1]["seed"] == 44 for c in refined)
        assert all(c["trials"][-1]["status"] == "failed" for c in refined)
    else:
        assert rows and DESIGN in {r["smiles"] for r in rows}
        assert all(r["calibration_status"] == expected_status for r in results["results"])
        assert all(set(r["selective_evidence"]["seed_witness_sha256"]) == {"42", "43", "44"} for r in rows)


@pytest.mark.parametrize("ablation", ["no_extra_refill_capacity", "design_only_route_search"])
def test_v6_ablation_rejects_unrelated_custom_base(v6_full_request, monkeypatch, ablation):
    from benchmarks.research import qualified_adapter as adapter
    from benchmarks.research.protocol import write_json

    request, predictor = v6_full_request
    request["settings"]["ablation"] = ablation
    path = Path(request["assets"]["synthesis_config"])
    config = adapter.read(path)
    config["ranking_policy"] = "joint-v1"
    write_json(path, config)
    monkeypatch.setattr(workflow, "run", lambda *a, **k: pytest.fail("custom base reached workflow"))
    with pytest.raises(ValueError, match="must match frozen engineering profile"):
        adapter.common(request)
    assert predictor.calls == []


def update(data, name, mutate):
    out, response, _ = data
    path = out / "synthesis" / name
    value = json.loads(path.read_text())
    mutate(value)
    path.write_text(json.dumps(value))
    response["metadata"]["files_sha256"]["synthesis/" + name] = digest(path)


def test_real_workflow_collector_keeps_robust_scores_separate_from_pose(completed):
    rows = collect(completed)
    assert len(rows) == 1 and rows[0]["smiles"] == DESIGN
    assert rows[0]["metrics"]["hit_tanimoto"] == pytest.approx(0.6)
    assert all(value >= 0 for value in rows[0]["margins"].values())
    assert rows[0]["bit_morgan_bits"] == sorted(set(rows[0]["bit_morgan_bits"]))
    assert len(completed[2].calls) == 1  # collector did not call predictor


@pytest.mark.parametrize("kind", ["seed", "raw", "route", "source", "row", "tsv", "sdf", "hash"])
def test_saved_artifact_tampering_rejected(completed, kind):
    if kind == "seed":
        update(completed, "refinement.json", lambda v: v["candidates"][0]["trials"][0].update(seed=97))
    elif kind == "raw":

        def change(v):
            record = next(r for r in v["routes"] if "synthesis_evidence" in r)
            record["synthesis_evidence"]["steps"][0]["forward_prediction"]["raw"]["beam_indices"][0][0] = -1

        update(completed, "routes.json", change)
    elif kind == "route":
        update(completed, "routes.json", lambda v: v["routes"][0]["route"].update(root="missing"))
    elif kind == "source":
        update(completed, "manifest.json", lambda v: v["code_sha256"].update({"workflow.py": "f" * 64}))
    elif kind == "row":
        update(completed, "results.json", lambda v: v["results"][0].update(smiles=HIT))
    else:
        out, response, _ = completed
        name = "results.tsv" if kind == "tsv" else "poses_hit.sdf"
        path = out / "synthesis" / name
        path.write_text(path.read_text().replace(DESIGN, HIT))
        if kind != "hash":
            response["metadata"]["files_sha256"]["synthesis/" + name] = digest(path)
    with pytest.raises((ValueError, KeyError, OSError)):
        collect(completed)


def test_token_replay_cannot_call_model():
    replay = object.__new__(collector.TokenReplay)
    with pytest.raises(RuntimeError, match="inference is forbidden"):
        replay.predict_raw({})


def test_bound_catalog_preserves_outputs_with_fresh_rules(completed, assets, monkeypatch):
    from chembounce_synthesis import reactions as rules

    bound = collector.BoundCatalog(
        assets[0].path, assets[0].sha256, workflow._SOURCE_HASHES, sidecar_sha256=None
    )
    original = collect(completed)
    instances = []
    constructor = rules.Rules

    def fresh(*args, **kwargs):
        instance = constructor(*args, **kwargs)
        instances.append(instance)
        return instance

    monkeypatch.setattr(rules, "Rules", fresh)
    assert collect(completed, bound) == original
    assert collect(completed, bound) == original
    assert len(instances) == 2 and instances[0] is not instances[1]
    assert bound._catalog.nearest.cache_info().currsize == 0
    assert all(v is bound._catalog.fps._PENDING for v in bound._catalog.fps._values.values())
    assert len(completed[2].calls) == 1


def is_selective(data):
    out, _, _ = data
    return json.loads((out / "request.json").read_text()).get("benchmark_protocol") is not None


def test_selective_collector_retains_all_seed_evidence_without_embedding(completed, monkeypatch):
    if not is_selective(completed):
        return
    from chembounce_synthesis import conformers, selective_scoring

    def forbidden(*args, **kwargs):
        raise AssertionError("Collector may not embed, align or score new poses")

    monkeypatch.setattr(conformers, "ensemble", forbidden)
    monkeypatch.setattr(selective_scoring, "ensemble", forbidden)
    monkeypatch.setattr(selective_scoring, "score_candidate", forbidden)
    result = collect(completed)[0]["selective_evidence"]
    assert set(result["seed_witness_sha256"]) == {"42", "43", "44"}
    assert set(result["seed_pose_evidence_sha256"]) == {"42", "43", "44"}
    assert result["all_seed_min_polar_reference_recall"] >= 0.75
    assert len(completed[2].calls) == 1


@pytest.mark.parametrize(
    "kind",
    [
        "missing_pose",
        "nonrepresentative_witness",
        "policy",
        "provenance",
        "seed_missing",
        "seed_failed",
        "row_witness",
        "pose_coordinates",
        "protocol",
    ],
)
def test_selective_collector_rejects_evidence_tamper(completed, kind):
    if not is_selective(completed):
        return
    if kind in {"missing_pose", "nonrepresentative_witness", "seed_missing", "seed_failed"}:

        def change(value):
            trials = value["candidates"][0]["trials"]
            if kind == "missing_pose":
                trials[2].pop("selective_pose_evidence")
            elif kind == "nonrepresentative_witness":
                row = json.loads((completed[0] / "synthesis/results.json").read_text())["results"][0]
                trial = next(t for t in trials if t["seed"] != row["pose_seed"])
                trial["similarity_constraints"]["witness"]["polar"]["recall"] = 1
            elif kind == "seed_missing":
                trials.pop()
            else:
                trials[2]["status"] = "failed"

        update(completed, "refinement.json", change)
    elif kind in {"policy", "provenance"}:

        def change(value):
            if kind == "policy":
                value["policy"]["similarity_constraints"]["policy"]["polar_recall"] = 0.01
            else:
                value["policy"]["calibration"]["current_numerical_runtime"]["source_sha256"].pop(
                    "chembounce_synthesis/selective_scoring.py"
                )

        update(completed, "results.json", change)
    elif kind == "row_witness":
        update(
            completed, "results.json", lambda v: v["results"][0].update(similarity_witness_sha256="0" * 64)
        )
    elif kind == "protocol":
        path = completed[0] / "request.json"
        request = json.loads(path.read_text())
        request.pop("benchmark_protocol")
        path.write_text(json.dumps(request))
    else:
        from rdkit import Chem

        path = completed[0] / "synthesis/references_hit.sdf"
        mol = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)))
        xyz = mol.GetConformer().GetAtomPosition(0)
        xyz.x += 0.1
        mol.GetConformer().SetAtomPosition(0, xyz)
        with Chem.SDWriter(str(path)) as writer:
            writer.write(mol)
        completed[1]["metadata"]["files_sha256"]["synthesis/references_hit.sdf"] = digest(path)
    with pytest.raises((ValueError, KeyError)):
        collect(completed)


def test_selective_failed_seed_exports_empty_and_collector_keeps_failure(completed, tmp_path, monkeypatch):
    if not is_selective(completed):
        return
    from benchmarks.research import qualified_adapter as adapter
    from chembounce_synthesis import selective_scoring

    original = selective_scoring.score_candidate

    def fail_final_seed(smiles, hit, design, settings, **kwargs):
        if settings.seed == 44:
            raise ValueError("controlled final seed failure")
        return original(smiles, hit, design, settings, **kwargs)

    monkeypatch.setattr(selective_scoring, "score_candidate", fail_final_seed)
    request = adapter.read(completed[0] / "request.json")
    out = tmp_path / "failed-seed-common"
    out.mkdir()
    request["output_dir"] = str(out)
    (out / "request.json").write_text(json.dumps(request))
    response = adapter.common(request)
    response["metadata"]["benchmark_protocol"] = request["benchmark_protocol"]
    assert response["candidates"] == []
    assert collect((out, response, completed[2])) == []
    trials = adapter.read(out / "synthesis/refinement.json")["candidates"][0]["trials"]
    assert [t["seed"] for t in trials] == [42, 43, 44]
    assert trials[-1]["status"] == "failed"
    assert adapter.read(out / "synthesis/summary.json")["status"] == "no_candidates"


@pytest.fixture
def empty_selective(completed, tmp_path, monkeypatch):
    from benchmarks.research import qualified_adapter as adapter
    from chembounce_synthesis import selective_scoring

    assert is_selective(completed)
    original = selective_scoring.score_candidate

    def fail_seed(smiles, hit, design, settings, **kwargs):
        if settings.seed == 44:
            raise ValueError("controlled final-seed abstention")
        return original(smiles, hit, design, settings, **kwargs)

    monkeypatch.setattr(selective_scoring, "score_candidate", fail_seed)
    request = adapter.read(completed[0] / "request.json")
    out = tmp_path / "empty-selective"
    out.mkdir()
    request["output_dir"] = str(out)
    (out / "request.json").write_text(json.dumps(request))
    response = adapter.common(request)
    response["metadata"]["benchmark_protocol"] = request["benchmark_protocol"]
    data = out, response, completed[2]
    assert collect(data) == []  # A genuine failed seed remains an auditable empty result.
    return data


@pytest.mark.parametrize("completed", SELECTIVE_PROFILES, indirect=True)
@pytest.mark.parametrize(
    "name",
    ["results.tsv", "poses_hit.sdf", "poses_design.sdf", "references_hit.sdf", "references_design.sdf"],
)
@pytest.mark.parametrize("attack", ["missing", "nonempty"])
def test_selective_empty_requires_complete_empty_exports(empty_selective, completed, name, attack):
    path = empty_selective[0] / "synthesis" / name
    if attack == "missing":
        path.unlink()
        del empty_selective[1]["metadata"]["files_sha256"]["synthesis/" + name]
    else:
        path.write_bytes((completed[0] / "synthesis" / name).read_bytes())
        empty_selective[1]["metadata"]["files_sha256"]["synthesis/" + name] = digest(path)
    with pytest.raises(ValueError):
        collect(empty_selective)


@pytest.mark.parametrize("completed", SELECTIVE_PROFILES, indirect=True)
@pytest.mark.parametrize(
    "attack", ["header", "qualified_count", "accepted_count", "routes", "refinements", "decision"]
)
def test_selective_empty_denominators_and_decisions_agree(empty_selective, attack):
    if attack == "header":
        path = empty_selective[0] / "synthesis/results.tsv"
        path.write_text("smiles\n")
        empty_selective[1]["metadata"]["files_sha256"]["synthesis/results.tsv"] = digest(path)
    elif attack in {"qualified_count", "accepted_count"}:
        key = "qualified_products" if attack == "qualified_count" else "accepted_products"
        update(empty_selective, "summary.json", lambda v: v.update({key: 1}))
    elif attack == "routes":
        update(empty_selective, "routes.json", lambda v: v.update(routes=[]))
    elif attack == "refinements":
        update(empty_selective, "refinement.json", lambda v: v.update(candidates=[]))
    else:
        update(
            empty_selective,
            "qualification.json",
            lambda v: v["decisions"][0].update(qualified=True, status="qualified"),
        )
    with pytest.raises(ValueError):
        collect(empty_selective)


@pytest.mark.parametrize("completed", SELECTIVE_PROFILES, indirect=True)
@pytest.mark.parametrize(
    "field",
    [
        "design_shape",
        "design_pharmacophore",
        "design_tanimoto",
        "relative_energy",
        "scoring_objective",
        "ranking_policy",
        "scoring_stage",
        "seed_score_min",
        "seed_score_max",
        "seed_score_median",
        "seed_score_std",
        "ranking_score",
        "refinement_conformers",
    ],
)
def test_selective_representative_fields_match_actual_retained_trial(completed, field):
    import csv

    value = "forged" if field in {"scoring_objective", "ranking_policy", "scoring_stage"} else 0.123456789
    update(completed, "results.json", lambda v: v["results"][0].update({field: value}))
    path = completed[0] / "synthesis/results.tsv"
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields, table = reader.fieldnames, list(reader)
    table[0][field] = str(value)
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(table)
    completed[1]["metadata"]["files_sha256"]["synthesis/results.tsv"] = digest(path)
    with pytest.raises(ValueError):
        collect(completed)


@pytest.mark.parametrize("completed", SELECTIVE_PROFILES, indirect=True)
@pytest.mark.parametrize(
    "name,field,value",
    [
        ("references_hit.sdf", "source_conformer_id", "999"),
        ("references_design.sdf", "source_conformer_id", "999"),
        ("references_hit.sdf", "reference_role", "design"),
        ("references_design.sdf", "scoring_objective", "joint-v1"),
        ("poses_hit.sdf", "ranking_policy", "joint-v1"),
        ("poses_design.sdf", "design_shape", ".123456789"),
        ("poses_design.sdf", "design_pharmacophore", ".123456789"),
        ("poses_hit.sdf", "seed_score_min", ".123456789"),
        ("poses_hit.sdf", "verification", "experimental_validation=true"),
    ],
)
def test_selective_sdf_properties_match_retained_trial(completed, name, field, value):
    from rdkit import Chem

    path = completed[0] / "synthesis" / name
    mol = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)))
    mol.SetProp(field, value)
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    completed[1]["metadata"]["files_sha256"]["synthesis/" + name] = digest(path)
    with pytest.raises(ValueError):
        collect(completed)


@pytest.mark.parametrize("completed", SELECTIVE_PROFILES, indirect=True)
def test_selective_false_empty_cannot_hide_qualified_routes(completed):
    update(completed, "results.json", lambda value: value.update(results=[]))
    update(completed, "qualification.json", lambda value: value.update(selected_route_ids=[]))
    update(
        completed,
        "summary.json",
        lambda value: value.update(
            status="no_candidates",
            returned_products=0,
            accepted_products=0,
            qualified_products=0,
            refinement_accepted=0,
        ),
    )
    completed[1]["candidates"] = []
    with pytest.raises(ValueError):
        collect(completed)


@pytest.mark.parametrize("completed", ["qualified-selective-engineering-v5"], indirect=True)
@pytest.mark.parametrize("kind", ["contradictory_code", "projection", "counts"])
def test_collector_refuses_new_forward_diagnostic_contradictions(completed, kind):
    if kind == "counts":
        update(completed, "summary.json", lambda v: v["forward_evidence_outcomes"].update(step_records=999))
    else:

        def mutate(document):
            evidence = next(r["synthesis_evidence"] for r in document["routes"] if "synthesis_evidence" in r)
            if kind == "contradictory_code":
                evidence["steps"][0]["forward_prediction"]["reason_code"] = "predictor_call_budget_exhausted"
                # Contradictory cause remains invalid even if copied into every diagnostic layer.
                evidence["steps"][0]["forward_outcome"]["cause_code"] = "predictor_call_budget_exhausted"
                evidence["forward_outcome"]["step_outcomes"][0]["cause_code"] = (
                    "predictor_call_budget_exhausted"
                )
            else:
                evidence["forward_outcome"]["status"] = "unassessed"

        update(completed, "routes.json", mutate)
    with pytest.raises(ValueError, match="forward"):
        collect(completed)
    assert len(completed[2].calls) == 1


@pytest.mark.parametrize("completed", ["qualified-selective-engineering-v5"], indirect=True)
def test_collector_legacy_diagnostic_absence_still_requires_raw_support(completed):
    def remove(document):
        for record in document["routes"]:
            evidence = record.get("synthesis_evidence", {})
            evidence.pop("forward_outcome", None)
            for step in evidence.get("steps", []):
                step.pop("forward_outcome", None)
                step["forward_prediction"].pop("reason_code", None)

    update(completed, "routes.json", remove)
    update(completed, "summary.json", lambda v: v.pop("forward_evidence_outcomes"))
    assert len(collect(completed)) == 1
    assert len(completed[2].calls) == 1
