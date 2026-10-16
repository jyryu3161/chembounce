"""Real toy chemistry and controlled raw beams; never local model/tokenizer loading."""

import json

import pytest

from benchmarks.research import qualified_collector as collector
from benchmarks.research.protocol import digest
from chembounce_synthesis import workflow

HIT, DESIGN = "CC(=O)Nc1ccccc1", "CC(=O)Nc1ccncc1"


@pytest.fixture
def completed(tmp_path, assets, route, strict_predictor):
    pytest.importorskip("rdkit.Chem.rdShapeAlign")
    from rdkit.Chem import rdChemReactions
    from rdkit import Chem
    from benchmarks.research import qualified_adapter as adapter, qualified_runner as runner
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
        "recommendation_profile": "qualified-development-v1",
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
    if profile_provenance(workflow.load_config(path))["status"] != "development_protocol_matched":
        pytest.skip("Runtime does not match reviewed development calibration")
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
    write_json(out / "request.json", request)
    response = adapter.common(request)
    assert len(response["candidates"]) == 1
    return out, response, strict_predictor


def collect(data):
    out, response, decoder = data
    return collector.checked_products(
        out,
        response,
        {"id": "fixture", "target_id": "synthetic", "smiles": HIT},
        workflow._SOURCE_HASHES,
        decoder_factory=lambda *args: decoder,
    )


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
