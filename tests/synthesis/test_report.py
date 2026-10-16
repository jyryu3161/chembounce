"""Pure saved-artifact fixtures; no chemistry or scientific asset constructors."""

import copy
import hashlib
import json
import os
import statistics
import subprocess
import sys

import pytest

from chembounce_synthesis import report as r


HIT = "CC(=O)Nc1ccccc1"
PRODUCT = "CC(=O)Nc1ccncc1"


def pose(smiles, cid):
    value = {
        "smiles": smiles,
        "conformer_id": cid,
        "coordinates": [[0.0, 1.0, 2.0]],
        "atom_order_sha256": "a" * 64,
    }
    return {**value, "sha256": r._binding(value)}


def bundle():
    """Schema-shaped synthetic coordinates; never claimed to be chemical proof."""
    minimum = {"hit_tanimoto": 0.35, "hit_shape": 0.7, "hit_pharmacophore": 0.5}
    bound = {
        "version": "selective-similarity-constraints-v2",
        "hit": HIT,
        "minimum_scores": minimum,
        "minimum_polar_reference_recall": 0.75,
        "core_smiles": None,
        "landmarks": [],
    }
    policy = {
        "policy_version": "qualified-v1",
        "mode": "qualified",
        "calibration": {
            "name": "qualified-selective-engineering-v2",
            "status": "engineering_policy_matched",
            "artifact_sha256": "d" * 64,
        },
        "minimum_scores": minimum,
        "evidence_profile": "computational-route-support-v1",
        "conformers": {"scoring_objective": "original-hit-v1"},
        "ranking_policy": "original-hit-v1",
        "refinement": {"seeds": [42, 43, 44], "threshold_policy": "all_seeds"},
        "constraints": {"mode": "hop", "core_smiles": None},
        "similarity_constraints": {"policy": bound, "sha256": r._binding(bound)},
    }
    policy["conformers"].update(
        num_conformers=8,
        seed=42,
        energy_window=5.0,
        strain_penalty=0.05,
        max_iterations=500,
        pharmacophore_sigma=1.0,
        reject_unspecified_stereo=True,
    )
    policy["refinement"].update(num_conformers=8, adaptive_max_conformers=None, conformer_growth_factor=2)
    trials = []
    for seed, score, shape in [(42, 0.6, 0.72), (43, 0.7, 0.88), (44, 0.8, 0.81)]:
        candidate, reference = pose(PRODUCT, seed), pose(HIT, seed + 100)
        witness = {
            "policy_sha256": r._binding(bound),
            "passed": True,
            "reasons": [],
            "reference_pose": reference,
            "candidate_pose_sha256": candidate["sha256"],
            "joint_assignment": {"feasible": True},
            "polar": {"recall": 0.8, "status": "passed", "reference_count": 5},
        }
        trials.append(
            {
                "seed": seed,
                "status": "success",
                "score": score,
                "candidate_smiles": PRODUCT,
                "hit_smiles": HIT,
                "design_smiles": PRODUCT,
                "candidate_relative_energy": 1.2,
                "candidate_conformer_id": seed,
                "scoring_objective": "original-hit-v1",
                "hit": {"shape": shape, "pharmacophore": 0.65, "reference_conformer_id": seed + 100},
                "design": {"shape": 0.9, "pharmacophore": 0.9, "reference_conformer_id": seed},
                "settings": {"seed": seed, "num_conformers": 8, "scoring_objective": "original-hit-v1"},
                "similarity_constraints": {"policy": bound, "witness": witness},
                "pose_selection_policy": "joint-feasible-hit-pose-v2",
                "pose_selection": {"status": "selected"},
                "selective_pose_evidence": {
                    "hit_candidate": candidate,
                    "design_candidate": candidate,
                    "design_reference": candidate,
                },
            }
        )
    for trial in trials:
        trial["settings"] = {**policy["conformers"], "seed": trial["seed"]}
    rep = trials[1]
    scores = [t["score"] for t in trials]
    summary_fields = {
        "pose_seed": 43,
        "seed_count": 3,
        "seed_successes": 3,
        "seed_score_median": statistics.median(scores),
        "seed_score_min": min(scores),
        "seed_score_max": max(scores),
        "seed_score_range": max(scores) - min(scores),
        "seed_score_std": statistics.pstdev(scores),
    }
    row = {
        "smiles": PRODUCT,
        "design": PRODUCT,
        "route_id": "route-final",
        "alternative_route_ids": ["route-coarse"],
        "recommendation_status": "qualified",
        "qualification_policy": "qualified-v1",
        "recommendation_profile": policy["calibration"]["name"],
        "calibration_status": policy["calibration"]["status"],
        "calibration_artifact_sha256": policy["calibration"]["artifact_sha256"],
        "evidence_profile": policy["evidence_profile"],
        "route_support_binding": "route-sha",
        "hit_tanimoto": 0.6,
        "hit_morgan": 0.6,
        "scaffold_family_smiles": "c1ccncc1",
        "refinement_conformers": 8,
        "score": 0.7,
        "relative_energy": 1.2,
        "scoring_objective": "original-hit-v1",
        "ranking_policy": "original-hit-v1",
        **summary_fields,
        **{
            f"{role}_{metric}": rep[role][metric]
            for role in ("hit", "design")
            for metric in ("shape", "pharmacophore")
        },
        "similarity_witness_sha256": r._binding(rep["similarity_constraints"]["witness"]),
        "similarity_constraint_sha256": r._binding(bound),
        "similarity_constraint_policy": bound["version"],
        "pose_selection_policy": rep["pose_selection_policy"],
        "polar_reference_recall": 0.8,
    }
    decision = {
        "policy_version": "qualified-v1",
        "similarity_constraints": policy["similarity_constraints"],
        "route_id": "route-final",
        "status": "qualified",
        "qualified": True,
        "scores": {k: row[k] for k in minimum},
        "calibration": policy["calibration"],
        "reasons": [],
    }
    origins = [
        {
            "input_index": 0,
            "record_index": 0,
            "id_origin": "input_provided",
            "original_identifier": "generated-looking-001",
            "input_path": "/outside/supplier.csv",
            "input_sha256": "e" * 64,
        },
        {
            "input_index": 1,
            "record_index": 0,
            "id_origin": "structure_generated",
            "original_identifier": None,
            "input_path": "/outside/second.csv",
            "input_sha256": "f" * 64,
        },
    ]
    leaf = {
        "node_id": "bb1",
        "smiles": "CC(=O)O",
        "catalog_record": {
            "id": "generated-looking-001",
            "source": "provided supplier label",
            "snapshot": "2020-01-01",
            "availability": "in_stock",
            "catalog_sha256": "c" * 64,
            "supplier_sku_verified": False,
            "current_stock_verified": False,
            "identifier_provenance": {
                "status": "recorded",
                "manifest_sha256": "b" * 64,
                "source_records": origins,
            },
        },
    }
    step = {
        "node_id": "rxn",
        "product": PRODUCT,
        "step_binding_sha256": "step-sha",
        "reactant_slots": [{"slot": 0, "node_id": "bb1", "smiles": "CC(=O)O"}],
        "rule_id": "amide",
        "smarts": "synthetic fixture",
        "rule_source": "fixture",
        "enumeration": {"truncated": True},
    }
    support = {
        "schema_version": 2,
        "profile": "computational-route-support-v1",
        "product": PRODUCT,
        "route_binding_sha256": "route-sha",
        "evidence_level": "conditional_computational_selection_support",
        "eligibility": {"passed": True},
        "route_audit": {
            "product": PRODUCT,
            "root": "rxn",
            "route_binding_sha256": "route-sha",
            "steps": [step],
            "building_blocks": [leaf],
            "limitations": ["Conditions and yield unassessed"],
            "enumeration_truncated": True,
        },
        "steps": [
            {
                "node_id": "rxn",
                "product": PRODUCT,
                "step_binding_sha256": "step-sha",
                "reactants": ["CC(=O)O"],
                "rule_replay": "matched",
                "forward_prediction": {"status": "matched", "matched_rank": 1, "matching_policy": "exact"},
                "context": {"status": "supplied_unvalidated", "validated": False},
                "observations": {"status": "unavailable", "source_audited": False},
            }
        ],
    }
    route = {
        "route_id": "route-final",
        "design": PRODUCT,
        "route": {"root": "rxn", "nodes": [{"id": "rxn", "product": PRODUCT}]},
        "verification": {"product": PRODUCT},
        "qualification": decision,
        "selection_status": "qualified",
        "synthesis_evidence": support,
    }
    row["reaction_steps"] = 1
    route["route"]["nodes"] = [
        {"id": "bb1", "building_block_id": "generated-looking-001"},
        {"id": "rxn", "rule_id": "amide", "reactants": ["bb1"], "product": PRODUCT},
    ]
    route["verification"]["steps"] = [
        {"node_id": "rxn", "rule_id": "amide", "product": PRODUCT, "reactants": ["CC(=O)O"]}
    ]
    support["steps"][0]["eligibility"] = {"passed": True, "unmet_requirements": []}
    return {
        "manifest.json": {
            "schema_version": 1,
            "hit": HIT,
            "effective_conformers": policy["conformers"],
            "recommendation_policy": policy,
            "recommendation_profile": policy["calibration"],
            "config": {"minimum_scores": minimum},
            "effective_workflow": {
                "minimum_scores": minimum,
                "refinement": policy["refinement"],
                "similarity_constraints": policy["similarity_constraints"],
            },
        },
        "results.json": {"schema_version": 1, "policy": policy, "results": [row]},
        "qualification.json": {
            "schema_version": 1,
            "policy": policy,
            "selected_route_ids": ["route-final"],
            "decisions": [decision],
        },
        "summary.json": {
            "schema_version": 1,
            "status": "success",
            "recommendation_profile": policy["calibration"],
            "recommendation_mode": "qualified",
            "returned_products": 1,
            "accepted_products": 4,
            "qualified_products": 4,
            "limits_reached": ["max_results"],
            "route_rejections": {"qualified": 4, "rejected": 2},
        },
        "routes.json": {"schema_version": 1, "routes": [route]},
        "refinement.json": {
            "schema_version": 1,
            "candidates": [
                {
                    "route_id": "route-coarse",
                    "product": PRODUCT,
                    "design": PRODUCT,
                    "status": "success",
                    "trials": trials,
                    "num_conformers": 8,
                    "threshold_policy": "all_seeds",
                    **summary_fields,
                }
            ],
        },
        "scores.json": {
            "schema_version": 1,
            "scores": [{"route_id": "route-final", "product": PRODUCT, "scoring_stage": "refined", **rep}],
        },
    }


def save(path, value):
    path.mkdir(exist_ok=True)
    for name, doc in value.items():
        (path / name).write_text(json.dumps(doc), encoding="utf-8")
    return path


@pytest.fixture
def data():
    return copy.deepcopy(bundle())


def test_primary_fallback_and_minima(tmp_path, data):
    run = save(tmp_path, data)
    before = {p.name: p.read_bytes() for p in run.iterdir()}
    view = r.write_report(run)
    assert view["artifact_consistency"] == "consistent_recorded_bundle"
    candidate = view["candidates"][0]
    assert candidate["refinement"]["route_id"] != candidate["row"]["route_id"]
    assert candidate["derived_minima"]["hit_shape"] == 0.72
    assert candidate["row"]["hit_shape"] == 0.88
    output = (run / "review.html").read_text()
    for text in (
        "Each seed has its own same-pose witness",
        "Derived arithmetic",
        "supplied_unvalidated",
        "input_provided",
        "structure_generated",
        "supplier_sku_verified",
        "False",
        "Generic topology comparison",
        "not recorded",
    ):
        assert text in output
    assert all((run / n).read_bytes() == raw for n, raw in before.items())
    assert all(view["artifacts"][n]["sha256"] == hashlib.sha256(raw).hexdigest() for n, raw in before.items())
    assert r.render_html(view) == r.render_html(r.build_report(run))


@pytest.mark.parametrize(
    "change",
    [
        lambda d: d["routes.json"]["routes"][0].update(design="other"),
        lambda d: d["routes.json"]["routes"][0]["verification"].update(product="other"),
        lambda d: d["refinement.json"]["candidates"].append(
            copy.deepcopy(d["refinement.json"]["candidates"][0])
        ),
        lambda d: d["refinement.json"].update(candidates=[]),
        lambda d: d["refinement.json"]["candidates"][0]["trials"].pop(),
        lambda d: d["refinement.json"]["candidates"][0]["trials"][1].update(seed=42),
        lambda d: d["refinement.json"]["candidates"][0]["trials"][1].update(status="failed"),
        lambda d: d["refinement.json"]["candidates"][0]["trials"][1]["settings"].update(seed=44),
        lambda d: d["results.json"]["results"][0].update(pose_seed=42),
        lambda d: d["results.json"]["results"][0].update(score=0.123),
        lambda d: d["results.json"]["results"][0].update(design_shape=0.123),
        lambda d: d["results.json"]["results"][0].update(relative_energy=99),
        lambda d: d["results.json"]["results"][0].update(scoring_objective="joint-v1"),
        lambda d: d["results.json"]["results"][0].update(seed_score_min=0.123),
        lambda d: d["results.json"]["results"][0].update(similarity_witness_sha256="wrong"),
        lambda d: d["results.json"]["results"][0].update(route_support_binding="wrong"),
        lambda d: d["refinement.json"]["candidates"][0]["trials"][0]["similarity_constraints"][
            "witness"
        ].update(candidate_pose_sha256="wrong"),
        lambda d: d["refinement.json"]["candidates"][0]["trials"][0]["selective_pose_evidence"][
            "hit_candidate"
        ]["coordinates"][0].__setitem__(0, 9),
        lambda d: d["scores.json"]["scores"][0].update(pose_selection={"status": "different"}),
        lambda d: d["scores.json"]["scores"][0].update(hit_smiles="other"),
        lambda d: d["qualification.json"].update(selected_route_ids=[]),
        lambda d: d["summary.json"].update(returned_products=2),
        lambda d: d["summary.json"].update(status="error"),
        lambda d: d["results.json"]["results"].append(copy.deepcopy(d["results.json"]["results"][0])),
        lambda d: d["manifest.json"].update(recommendation_profile={}),
        lambda d: d["manifest.json"].update(recommendation_policy={}),
    ],
)
def test_corruptions_are_diagnostic_not_qualified(tmp_path, data, change):
    change(data)
    view = r.write_report(save(tmp_path, data))
    assert view["artifact_consistency"] == "incomplete_or_inconsistent"
    assert "Unconfirmed recorded claim" in (tmp_path / "review.html").read_text()


def test_same_product_competing_design_is_not_joined(tmp_path, data):
    other = copy.deepcopy(data["refinement.json"]["candidates"][0])
    other.update(design="another exact design", route_id="unrelated")
    other["trials"][0]["score"] = 999
    data["refinement.json"]["candidates"].insert(0, other)
    view = r.build_report(save(tmp_path, data))
    assert view["artifact_consistency"] == "consistent_recorded_bundle"
    assert view["candidates"][0]["refinement"]["design"] == PRODUCT


def test_legacy_leaves_remain_unknown(tmp_path, data):
    record = data["routes.json"]["routes"][0]["synthesis_evidence"]["route_audit"]["building_blocks"][0][
        "catalog_record"
    ]
    for key in ("identifier_provenance", "supplier_sku_verified", "current_stock_verified"):
        del record[key]
    view = r.write_report(save(tmp_path, data))
    assert view["artifact_consistency"] == "consistent_recorded_bundle"
    assert "unknown_legacy" in (tmp_path / "review.html").read_text()


def test_no_polar_and_identical_seed_witness_allowed(tmp_path, data):
    trials = data["refinement.json"]["candidates"][0]["trials"]
    for trial in trials:
        w = trial["similarity_constraints"]["witness"]
        w["polar"] = {"recall": None, "status": "not_applicable", "reference_count": 0}
    rep = trials[1]
    data["scores.json"]["scores"][0] = {
        "route_id": "route-final",
        "product": PRODUCT,
        "scoring_stage": "refined",
        **copy.deepcopy(rep),
    }
    data["results.json"]["results"][0].update(
        polar_reference_recall=None,
        similarity_witness_sha256=r._binding(rep["similarity_constraints"]["witness"]),
    )
    view = r.write_report(save(tmp_path, data))
    assert view["artifact_consistency"] == "consistent_recorded_bundle"
    assert view["candidates"][0]["seed_rows"][0]["polar_recall"] is None
    assert "not_applicable" in (tmp_path / "review.html").read_text()
    assert "not applicable" in view["candidates"][0]["derived_minima"]["polar_recall"]


@pytest.mark.parametrize("status", ["no_candidates", "error", "running"])
def test_empty_diagnostics_preserve_units_and_status(tmp_path, data, status):
    data["results.json"]["results"] = []
    data["qualification.json"]["selected_route_ids"] = []
    data["summary.json"].update(
        status=status,
        returned_products=0,
        no_candidates_reason="Deadline reached",
        deadline_reached=True,
        failures=[{"stage": "proposal", "reason": "model failed"}],
    )
    view = r.write_report(save(tmp_path, data))
    output = (tmp_path / "review.html").read_text()
    assert view["summary"]["status"] == status
    assert "No selected candidates" in output and "Deadline reached" in output
    assert "route records, not unique molecules" in output


@pytest.mark.parametrize(
    "bad",
    [
        '{"schema_version":1,"schema_version":1}',
        '{"schema_version":true}',
        '{"schema_version":2}',
        '{"schema_version":1,"bad":NaN}',
        '{"schema_version":1,"bad":1e400}',
        "[]",
        "{",
    ],
)
def test_malformed_json(tmp_path, data, bad):
    save(tmp_path, data)
    (tmp_path / "manifest.json").write_text(bad)
    view = r.write_report(tmp_path)
    assert view["artifact_consistency"] == "incomplete_or_inconsistent"


def test_partial_directory(tmp_path):
    (tmp_path / "summary.json").write_text(
        '{"schema_version":1,"status":"error","error":"early config failure"}'
    )
    view = r.write_report(tmp_path)
    assert view["artifact_consistency"] == "incomplete_or_inconsistent"
    assert "early config failure" in (tmp_path / "review.html").read_text()


def test_html_escaping_safe_links_and_external_output(tmp_path, data):
    data["summary.json"]["error"] = '<script src="https://evil">alert(1)</script>'
    data["routes.json"]["routes"][0]["synthesis_evidence"]["route_audit"]["building_blocks"][0][
        "catalog_record"
    ]["source"] = "javascript:alert(1)"
    run = save(tmp_path / "run", data)
    view = r.write_report(run)
    output = (run / "review.html").read_text()
    assert "<script" not in output and "&lt;script" in output
    assert "href='javascript:" not in output and 'href="javascript:' not in output
    assert "href='results.json'" in output
    r.write_report(run, tmp_path / "external.html")
    assert "href=" not in (tmp_path / "external.html").read_text()
    with pytest.raises(ValueError, match="Unsafe"):
        r.render_html(view, artifact_links={"results.json": "https://evil"})


def test_no_overwrite_or_symlink(tmp_path, data):
    run = save(tmp_path / "run", data)
    with pytest.raises(ValueError, match="new file"):
        r.write_report(run, run / "results.json")
    outside = tmp_path / "secret"
    outside.write_text("DO NOT READ")
    (run / "poses_hit.sdf").symlink_to(outside)
    with pytest.raises(ValueError):
        r.write_report(run)
    assert not (run / "review.html").exists()
    assert outside.read_text() == "DO NOT READ"


def test_changed_inputs_before_and_after_write(tmp_path, data, monkeypatch):
    run = save(tmp_path, data)
    original = r.render_html

    def mutation(*args, **kwargs):
        value = original(*args, **kwargs)
        (run / "summary.json").write_text("{}")
        return value

    monkeypatch.setattr(r, "render_html", mutation)
    with pytest.raises(ValueError, match="changed"):
        r.write_report(run)
    assert not (run / "review.html").exists()
    monkeypatch.setattr(r, "render_html", original)
    save(run, data)
    original_sync = os.fsync

    def during_write(fd):
        original_sync(fd)
        (run / "summary.json").write_text("{}")

    monkeypatch.setattr(os, "fsync", during_write)
    with pytest.raises(ValueError, match="changed"):
        r.write_report(run)
    assert not (run / "review.html").exists()


def test_bounds(tmp_path, data, monkeypatch):
    save(tmp_path, data)
    monkeypatch.setattr(r, "MAX_TOTAL_BYTES", 1)
    view = r.build_report(tmp_path)
    assert any("256 MiB" in item for item in view["issues"])


@pytest.mark.parametrize("with_identity", [False, True])
def test_stdlib_import_and_file_access_sentinels(tmp_path, data, with_identity):
    run = save(tmp_path, identity_bundle() if with_identity else data)
    # Disposable child: a permanent audit hook cannot contaminate other tests.
    script = """import sys, pathlib
root = pathlib.Path(sys.argv[1]).resolve()
def guard(event, args):
    if event == 'import' and args[0].split('.')[0] in {'rdkit','numpy','torch','transformers'}:
        raise RuntimeError('scientific import forbidden')
    if event.startswith('socket.') or event == 'subprocess.Popen':
        raise RuntimeError('network/process forbidden')
    if event == 'open' and isinstance(args[0], (str, bytes)):
        p = pathlib.Path(args[0]).resolve()
        if p.suffix in {'.tsv','.sdf','.json'} and p.parent != root:
            raise RuntimeError('external artifact/asset access forbidden')
sys.addaudithook(guard)
from chembounce_synthesis.report import main
raise SystemExit(main([str(root)]))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script, str(run)], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, completed.stderr
    assert "no scientific revalidation" in completed.stdout


def test_cli_exit_codes(tmp_path, data):
    run = save(tmp_path, data)
    assert r.main([str(run)]) == 0
    with pytest.raises(SystemExit) as exc:
        r.main([str(run)])
    assert exc.value.code == 2
    (run / "manifest.json").unlink()
    assert r.main([str(run), "--output", str(run / "partial.html")]) == 3


@pytest.mark.parametrize(
    "path",
    [
        ("results.json", "policy"),
        ("results.json", "policy", "calibration"),
        ("results.json", "policy", "constraints"),
        ("results.json", "policy", "similarity_constraints"),
        ("results.json", "policy", "refinement"),
        ("routes.json", "routes", 0, "synthesis_evidence"),
        ("routes.json", "routes", 0, "synthesis_evidence", "route_audit"),
        ("routes.json", "routes", 0, "synthesis_evidence", "steps"),
        (
            "routes.json",
            "routes",
            0,
            "synthesis_evidence",
            "route_audit",
            "building_blocks",
            0,
            "catalog_record",
        ),
        ("refinement.json", "candidates", 0, "trials", 0, "similarity_constraints", "witness"),
    ],
)
@pytest.mark.parametrize("wrong_type", [[1], "invalid", None])
def test_nested_wrong_types_render_diagnostics(tmp_path, data, path, wrong_type):
    # JSON roundtrip breaks fixture aliases: mutate just this one saved artifact.
    data = json.loads(json.dumps(data))
    parent = data
    for key in path[:-1]:
        parent = parent[key]
    parent[path[-1]] = wrong_type
    view = r.write_report(save(tmp_path, data))
    assert view["artifact_consistency"] == "incomplete_or_inconsistent"
    assert (tmp_path / "review.html").is_file()


def test_failed_seed_preserved_not_silently_omitted(tmp_path, data):
    data["refinement.json"]["candidates"][0]["trials"][1] = {
        "seed": 43,
        "status": "not_evaluated",
        "reason": "deadline before seed",
    }
    view = r.write_report(save(tmp_path, data))
    assert [t["seed"] for t in view["candidates"][0]["seed_rows"]] == [42, 43, 44]
    assert view["candidates"][0]["seed_rows"][1]["status"] == "not_evaluated"
    assert view["candidates"][0]["derived_minima"] == {}


def test_identical_pose_bindings_across_seeds_not_rejected(tmp_path, data):
    trials = data["refinement.json"]["candidates"][0]["trials"]
    rep = trials[1]
    for trial in trials:
        trial["candidate_conformer_id"] = rep["candidate_conformer_id"]
        trial["hit"]["reference_conformer_id"] = rep["hit"]["reference_conformer_id"]
        trial["design"]["reference_conformer_id"] = rep["design"]["reference_conformer_id"]
        trial["similarity_constraints"] = copy.deepcopy(rep["similarity_constraints"])
        trial["selective_pose_evidence"] = copy.deepcopy(rep["selective_pose_evidence"])
    assert r.build_report(save(tmp_path, data))["artifact_consistency"] == "consistent_recorded_bundle"


def test_read_bound_and_presence_race(tmp_path, data, monkeypatch):
    save(tmp_path, data)
    original = r.render_html

    def add_artifact(*args, **kwargs):
        output = original(*args, **kwargs)
        (tmp_path / "poses_hit.sdf").write_text("new artifact")
        return output

    monkeypatch.setattr(r, "render_html", add_artifact)
    with pytest.raises(ValueError, match="presence changed"):
        r.write_report(tmp_path)
    assert not (tmp_path / "review.html").exists()
    monkeypatch.setattr(r, "render_html", original)
    monkeypatch.setattr(r, "MAX_FILE_BYTES", 3)
    with pytest.raises(ValueError):
        r.write_report(tmp_path)
    assert not (tmp_path / "review.html").exists()


def diagnostic(tmp_path, data):
    view = r.write_report(save(tmp_path, data))
    assert view["artifact_consistency"] == "incomplete_or_inconsistent"
    assert "Unconfirmed recorded claim" in (tmp_path / "review.html").read_text()
    return view


@pytest.mark.parametrize(
    "role,metric",
    [("hit", "shape"), ("hit", "pharmacophore"), ("design", "shape"), ("design", "pharmacophore")],
)
@pytest.mark.parametrize("value", [-0.01, 1.01, True])
def test_similarity_domains_all_seeds(tmp_path, data, role, metric, value):
    data["refinement.json"]["candidates"][0]["trials"][0][role][metric] = value
    diagnostic(tmp_path, data)


@pytest.mark.parametrize("field,value", [("shape", 0.699), ("pharmacophore", 0.499)])
def test_nonrepresentative_below_floor(tmp_path, data, field, value):
    data["refinement.json"]["candidates"][0]["trials"][0]["hit"][field] = value
    diagnostic(tmp_path, data)


@pytest.mark.parametrize("value", [-0.01, 0.749, 1.01, True])
def test_nonrepresentative_polar_domain_and_floor(tmp_path, data, value):
    data["refinement.json"]["candidates"][0]["trials"][0]["similarity_constraints"]["witness"]["polar"][
        "recall"
    ] = value
    diagnostic(tmp_path, data)


@pytest.mark.parametrize("value", [-0.01, 0.349, 1.01, True])
def test_selected_tanimoto_domain_and_floor(tmp_path, data, value):
    data["results.json"]["results"][0]["hit_tanimoto"] = value
    data["qualification.json"]["decisions"][0]["scores"]["hit_tanimoto"] = value
    diagnostic(tmp_path, data)


def test_floor_equality_is_permitted(tmp_path, data):
    trial = data["refinement.json"]["candidates"][0]["trials"][0]
    trial["hit"].update(shape=0.7, pharmacophore=0.5)
    trial["similarity_constraints"]["witness"]["polar"]["recall"] = 0.75
    assert r.build_report(save(tmp_path, data))["artifact_consistency"] == "consistent_recorded_bundle"


@pytest.mark.parametrize("value", [0, -0.1, 1.1, True])
def test_invalid_synchronized_floors(tmp_path, data, value):
    # Existing aliases synchronize every saved copy; refresh claimed hashes too.
    policy = data["results.json"]["policy"]
    policy["minimum_scores"]["hit_shape"] = value
    bound = policy["similarity_constraints"]["policy"]
    policy["similarity_constraints"]["sha256"] = r._binding(bound)
    for trial in data["refinement.json"]["candidates"][0]["trials"]:
        trial["similarity_constraints"]["witness"]["policy_sha256"] = r._binding(bound)
    row = data["results.json"]["results"][0]
    row["similarity_constraint_sha256"] = r._binding(bound)
    row["similarity_witness_sha256"] = r._binding(
        data["refinement.json"]["candidates"][0]["trials"][1]["similarity_constraints"]["witness"]
    )
    diagnostic(tmp_path, data)


@pytest.mark.parametrize(
    "key,value",
    [
        ("energy_window", 6.0),
        ("strain_penalty", 0.2),
        ("max_iterations", 501),
        ("pharmacophore_sigma", 2.0),
        ("reject_unspecified_stereo", False),
        ("num_conformers", 1),
        ("scoring_objective", "joint-v1"),
        ("seed", 40),
    ],
)
def test_complete_trial_settings_are_policy_bound(tmp_path, data, key, value):
    data["refinement.json"]["candidates"][0]["trials"][0]["settings"][key] = value
    diagnostic(tmp_path, data)


@pytest.mark.parametrize(
    "key,value",
    [
        ("energy_window", True),
        ("strain_penalty", -1),
        ("max_iterations", 1.5),
        ("pharmacophore_sigma", 0),
        ("reject_unspecified_stereo", 1),
        ("num_conformers", True),
        ("seed", -1),
    ],
)
def test_invalid_synchronized_settings(tmp_path, data, key, value):
    policy = data["results.json"]["policy"]
    policy["conformers"][key] = value
    for trial in data["refinement.json"]["candidates"][0]["trials"]:
        trial["settings"][key] = value
    diagnostic(tmp_path, data)


def set_final_size(data, size):
    ref = data["refinement.json"]["candidates"][0]
    ref["num_conformers"] = size
    data["results.json"]["results"][0]["refinement_conformers"] = size
    for trial in ref["trials"]:
        trial["settings"]["num_conformers"] = size


@pytest.mark.parametrize("size", [1, 12, 24])
def test_custom_fixed_size_uses_policy_not_hardcoded_eight(tmp_path, data, size):
    data["results.json"]["policy"]["refinement"]["num_conformers"] = size
    set_final_size(data, size)
    assert r.build_report(save(tmp_path, data))["artifact_consistency"] == "consistent_recorded_bundle"


@pytest.mark.parametrize("size", [1, 9, True, 8.0])
def test_final_size_must_follow_typed_policy(tmp_path, data, size):
    set_final_size(data, size)
    diagnostic(tmp_path, data)


def adaptive_bundle(data):
    ref = data["refinement.json"]["candidates"][0]
    data["results.json"]["policy"]["refinement"]["adaptive_max_conformers"] = 24
    set_final_size(data, 24)
    stages = []
    for size in (8, 16, 24):
        stage = copy.deepcopy(ref)
        stage["num_conformers"] = size
        for trial in stage["trials"]:
            trial["settings"]["num_conformers"] = size
        stage["status"] = "failed" if size < 24 else "success"
        stages.append(stage)
    ref.update(stages=stages, adaptive_stages=3)
    return data


def test_adaptive_final_has_recorded_schedule_and_complete_final_stage(tmp_path, data):
    assert (
        r.build_report(save(tmp_path, adaptive_bundle(data)))["artifact_consistency"]
        == "consistent_recorded_bundle"
    )


@pytest.mark.parametrize(
    "mutation",
    [
        lambda ref: ref.pop("stages"),
        lambda ref: ref["stages"][1].update(num_conformers=12),
        lambda ref: ref["stages"][-1]["trials"][0].update(status="failed"),
        lambda ref: ref.update(adaptive_stages=2),
    ],
)
def test_adaptive_missing_or_contradictory_stages(tmp_path, data, mutation):
    adaptive_bundle(data)
    mutation(data["refinement.json"]["candidates"][0])
    diagnostic(tmp_path, data)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda route: route["synthesis_evidence"].update(steps=[]),
        lambda route: route["synthesis_evidence"]["route_audit"].update(steps=[]),
        lambda route: route["verification"].update(steps=[]),
        lambda route: route["synthesis_evidence"]["steps"][0]["eligibility"].update(passed=False),
        lambda route: route["synthesis_evidence"]["steps"][0].update(rule_replay="mismatched"),
        lambda route: route["synthesis_evidence"]["steps"][0].update(product="CC"),
        lambda route: route["synthesis_evidence"]["route_audit"]["steps"][0].update(product="CC"),
        lambda route: route["verification"]["steps"][0].update(rule_id="other"),
        lambda route: route["route"]["nodes"][-1].update(reactants=["rxn"]),
        lambda route: route["synthesis_evidence"]["route_audit"]["steps"][0]["reactant_slots"][0].update(
            node_id="rxn"
        ),
        lambda route: route["synthesis_evidence"]["route_audit"]["steps"][0]["reactant_slots"][0].update(
            slot=True
        ),
        lambda route: route["synthesis_evidence"]["steps"].append(
            copy.deepcopy(route["synthesis_evidence"]["steps"][0])
        ),
    ],
)
def test_graph_verification_audit_and_support_correspondence(tmp_path, data, mutation):
    mutation(data["routes.json"]["routes"][0])
    diagnostic(tmp_path, data)


def test_repeated_reactant_slots_are_not_deduplicated(tmp_path, data):
    route = data["routes.json"]["routes"][0]
    route["route"]["nodes"][-1]["reactants"] = ["bb1", "bb1"]
    route["verification"]["steps"][0]["reactants"] *= 2
    route["synthesis_evidence"]["steps"][0]["reactants"] *= 2
    slots = route["synthesis_evidence"]["route_audit"]["steps"][0]["reactant_slots"]
    slots.append({**slots[0], "slot": 1})
    assert r.build_report(save(tmp_path, data))["artifact_consistency"] == "consistent_recorded_bundle"


@pytest.mark.parametrize(
    "key",
    [
        "accepted_products",
        "qualified_products",
        "returned_products",
        "verified_routes",
        "scored_products",
        "refinement_attempted",
    ],
)
@pytest.mark.parametrize("value", [-1, True, 0.5])
def test_counters_require_nonnegative_integers(tmp_path, data, key, value):
    data["summary.json"][key] = value
    diagnostic(tmp_path, data)


@pytest.mark.parametrize("key", ["accepted_products", "qualified_products"])
def test_selected_count_cannot_exceed_accepted_counts(tmp_path, data, key):
    data["summary.json"][key] = 0
    diagnostic(tmp_path, data)


def test_partial_running_counts_do_not_require_terminal_equality(tmp_path, data):
    data["results.json"]["results"] = []
    data["qualification.json"]["selected_route_ids"] = []
    data["summary.json"].update(status="running", returned_products=0, accepted_products=0)
    data["summary.json"].pop("qualified_products")
    assert r.build_report(save(tmp_path, data))["artifact_consistency"] == "consistent_recorded_bundle"


@pytest.mark.parametrize(
    "mutation",
    [
        lambda route: route["route"]["nodes"].append({"id": "unused", "building_block_id": "unused"}),
        lambda route: route["synthesis_evidence"]["route_audit"].update(strict_eligible=False),
        lambda route: route["synthesis_evidence"]["route_audit"]["steps"][0].update(strict_eligible=False),
        lambda route: route["synthesis_evidence"]["steps"][0]["eligibility"].update(
            unmet_requirements=["failed"]
        ),
    ],
)
def test_additional_direct_graph_claim_contradictions(tmp_path, data, mutation):
    mutation(data["routes.json"]["routes"][0])
    diagnostic(tmp_path, data)


def test_no_candidate_terminal_counts_are_not_running_counts(tmp_path, data):
    data["results.json"]["results"] = []
    data["qualification.json"]["selected_route_ids"] = []
    data["summary.json"].update(status="no_candidates", returned_products=0)
    view = r.build_report(save(tmp_path, data))
    assert view["artifact_consistency"] == "incomplete_or_inconsistent"


def test_stage_size_scalar_type_is_not_coerced(tmp_path, data):
    adaptive_bundle(data)
    data["refinement.json"]["candidates"][0]["stages"][0]["num_conformers"] = 8.0
    diagnostic(tmp_path, data)


def test_policy_omissions_remain_diagnostic_not_filled_from_defaults(tmp_path, data):
    data["results.json"]["policy"]["conformers"].pop("energy_window")
    diagnostic(tmp_path, data)


def nonselective_custom(data):
    policy = data["results.json"]["policy"]
    policy.pop("similarity_constraints")
    policy["calibration"].update(name=None, status="custom_unvalidated", artifact_sha256=None)
    data["manifest.json"]["effective_workflow"].pop("similarity_constraints")
    data["qualification.json"]["decisions"][0].pop("similarity_constraints")
    row = data["results.json"]["results"][0]
    row.update(
        recommendation_profile="none",
        calibration_status="custom_unvalidated",
        calibration_artifact_sha256="none",
    )
    for trial in data["refinement.json"]["candidates"][0]["trials"] + data["scores.json"]["scores"]:
        for key in (
            "similarity_constraints",
            "selective_pose_evidence",
            "pose_selection",
            "pose_selection_policy",
        ):
            trial.pop(key, None)
    for key in (
        "pose_selection_policy",
        "similarity_constraint_policy",
        "similarity_constraint_sha256",
        "similarity_witness_sha256",
        "polar_reference_recall",
    ):
        row.pop(key, None)
    return data


@pytest.mark.parametrize("key", ["design_shape", "design_pharmacophore", "design_tanimoto"])
@pytest.mark.parametrize("floor", [0, 0.5])
def test_optional_design_floors_allow_zero_and_positive(tmp_path, data, key, floor):
    nonselective_custom(data)
    row = data["results.json"]["results"][0]
    row.setdefault("design_tanimoto", 0.9)
    data["results.json"]["policy"]["minimum_scores"][key] = floor
    data["qualification.json"]["decisions"][0]["scores"][key] = row[key]
    assert r.build_report(save(tmp_path, data))["artifact_consistency"] == "consistent_recorded_bundle"


@pytest.mark.parametrize("floor", [-0.01, 1.01, True])
def test_optional_design_floor_domain_not_relaxed(tmp_path, data, floor):
    nonselective_custom(data)
    data["results.json"]["policy"]["minimum_scores"]["design_shape"] = floor
    data["qualification.json"]["decisions"][0]["scores"]["design_shape"] = 0.9
    diagnostic(tmp_path, data)


@pytest.mark.parametrize("key", ["hit_tanimoto", "hit_shape", "hit_pharmacophore"])
def test_nonselective_required_hit_floors_still_positive(tmp_path, data, key):
    nonselective_custom(data)
    data["results.json"]["policy"]["minimum_scores"][key] = 0
    diagnostic(tmp_path, data)


@pytest.mark.parametrize("identifier", ["different-id", "", None, 7, ["generated-looking-001"]])
def test_route_leaf_id_must_match_saved_catalog_record(tmp_path, data, identifier):
    data["routes.json"]["routes"][0]["route"]["nodes"][0]["building_block_id"] = identifier
    diagnostic(tmp_path, data)


def test_catalog_identifier_cannot_be_missing(tmp_path, data):
    del data["routes.json"]["routes"][0]["synthesis_evidence"]["route_audit"]["building_blocks"][0][
        "catalog_record"
    ]["id"]
    diagnostic(tmp_path, data)


def test_cross_leaf_origin_borrowing_rejected(tmp_path, data):
    route = data["routes.json"]["routes"][0]
    route["route"]["nodes"].insert(1, {"id": "bb2", "building_block_id": "second-id"})
    route["route"]["nodes"][-1]["reactants"].append("bb2")
    route["verification"]["steps"][0]["reactants"].append("CC(=O)O")
    support = route["synthesis_evidence"]
    support["steps"][0]["reactants"].append("CC(=O)O")
    audit = support["route_audit"]
    second = copy.deepcopy(audit["building_blocks"][0])
    second["node_id"] = "bb2"
    second["catalog_record"]["id"] = "second-id"
    audit["building_blocks"].append(second)
    audit["steps"][0]["reactant_slots"].append({"slot": 1, "node_id": "bb2", "smiles": "CC(=O)O"})
    valid = tmp_path / "valid"
    assert r.build_report(save(valid, data))["artifact_consistency"] == "consistent_recorded_bundle"
    first = audit["building_blocks"][0]
    first["catalog_record"], second["catalog_record"] = second["catalog_record"], first["catalog_record"]
    invalid = tmp_path / "swapped"
    invalid.mkdir()
    diagnostic(invalid, data)


def selective_design_floor(data, key, floor):
    policy = data["results.json"]["policy"]
    # Producer BoundSimilarity owns the three-hit projection, not optional design floors.
    policy["similarity_constraints"]["policy"]["minimum_scores"] = {
        name: policy["minimum_scores"][name] for name in ("hit_tanimoto", "hit_shape", "hit_pharmacophore")
    }
    policy["minimum_scores"][key] = floor
    row = data["results.json"]["results"][0]
    row.setdefault("design_tanimoto", 0.9)
    data["qualification.json"]["decisions"][0]["scores"][key] = row[key]
    return data


@pytest.mark.parametrize("key", ["design_shape", "design_pharmacophore", "design_tanimoto"])
@pytest.mark.parametrize("floor", [0, 0.5])
def test_selective_optional_design_floors_are_outside_hit_witness(tmp_path, data, key, floor):
    selective_design_floor(data, key, floor)
    assert r.build_report(save(tmp_path, data))["artifact_consistency"] == "consistent_recorded_bundle"


@pytest.mark.parametrize(
    "mutation",
    [
        lambda scores: scores.pop("hit_shape"),
        lambda scores: scores.update(design_shape=0),
        lambda scores: scores.update(hit_pharmacophore=0.6),
    ],
)
def test_witness_requires_exact_hit_projection_even_with_rebound_hashes(tmp_path, data, mutation):
    selective_design_floor(data, "design_shape", 0)
    bound = data["results.json"]["policy"]["similarity_constraints"]["policy"]
    mutation(bound["minimum_scores"])
    digest = r._binding(bound)
    data["results.json"]["policy"]["similarity_constraints"]["sha256"] = digest
    for trial in data["refinement.json"]["candidates"][0]["trials"]:
        trial["similarity_constraints"]["witness"]["policy_sha256"] = digest
    row = data["results.json"]["results"][0]
    row["similarity_constraint_sha256"] = digest
    row["similarity_witness_sha256"] = r._binding(
        data["refinement.json"]["candidates"][0]["trials"][1]["similarity_constraints"]["witness"]
    )
    view = diagnostic(tmp_path, data)
    assert "witness hit-only floors" in str(view["candidates"][0]["issues"])


@pytest.mark.parametrize("key", ["design_shape", "design_pharmacophore", "design_tanimoto"])
def test_selective_selected_design_floors_still_enforced(tmp_path, data, key):
    selective_design_floor(data, key, 0.5)
    data["results.json"]["results"][0][key] = 0.4
    data["qualification.json"]["decisions"][0]["scores"][key] = 0.4
    view = diagnostic(tmp_path, data)
    assert f"Selected {key} contradicts floor" in str(view["candidates"][0]["issues"])


@pytest.mark.parametrize("metric", ["shape", "pharmacophore"])
def test_selective_nonrepresentative_design_floor_still_enforced(tmp_path, data, metric):
    selective_design_floor(data, f"design_{metric}", 0.5)
    data["refinement.json"]["candidates"][0]["trials"][0]["design"][metric] = 0.4
    view = diagnostic(tmp_path, data)
    assert f"Required seed design_{metric} contradicts floor" in str(view["candidates"][0]["issues"])


# Producer-shaped metadata; symbolic saved strings are deliberately not parsed.
def identity_bundle():
    docs = bundle()
    manifest = docs["manifest.json"]
    descriptor = {
        "schema_version": 1,
        "policy_version": "protonation-tautomer-family-v1",
        "rdkit_version": "2025.09.6",
        "tautomer_limit": 1000,
        "transform_limit": 1000,
        "code_sha256": {"identity.py": "a" * 64, "common.py": "b" * 64},
    }
    manifest["identity_comparison_policy"] = descriptor
    manifest["code_sha256"] = copy.deepcopy(descriptor["code_sha256"])
    manifest["versions"] = {"rdkit": "2025.9.6"}  # Package spelling is not rdBase spelling.
    row = docs["results.json"]["results"][0]
    row["family_smiles"] = PRODUCT

    def identity(exact, scaffold):
        return {
            "exact_smiles": exact,
            "family_smiles": exact,
            "scaffold_family_smiles": scaffold,
            "ring_skeleton_smiles": "C1CCCCC1",
            "has_scaffold": True,
            "policy_version": descriptor["policy_version"],
            "rdkit_version": descriptor["rdkit_version"],
        }

    row["hit_identity_comparison"] = {
        "schema_version": 1,
        "reference_role": "original_hit",
        "status": "recorded",
        "original_hit_smiles": HIT,
        "final_product_smiles": PRODUCT,
        "original_hit": identity(HIT, "c1ccccc1"),
        "final_product": identity(PRODUCT, "c1ccncc1"),
        "change_class": "scaffold",
        "atom_aware_scaffold_changed": True,
        "generic_ring_linker_topology_changed": False,
        "scaffold_presence_change": "both_present",
    }
    return docs


def identity_parts(docs):
    row = docs["results.json"]["results"][0]
    return row, row["hit_identity_comparison"], docs["manifest.json"]


def identity_view(docs):
    row, _, manifest = identity_parts(docs)
    return r._identity_comparison(
        row, manifest, docs["results.json"]["policy"], r._identity_descriptor(manifest)
    )


def unavailable(docs, reason="retained_identity_policy_mismatch"):
    row, record, _ = identity_parts(docs)
    row["hit_identity_comparison"] = {
        **{
            k: record[k]
            for k in ("schema_version", "reference_role", "original_hit_smiles", "final_product_smiles")
        },
        "status": "unavailable",
        "reason_code": reason,
        "reason": "<saved & unavailable>",
    }


def test_recorded_identity_public_roundtrip_and_defensive_copy(tmp_path):
    docs = identity_bundle()
    before = copy.deepcopy(docs)
    view = r._assemble(docs, {}, [])
    candidate = view["candidates"][0]
    assert view["artifact_consistency"] == "consistent_recorded_bundle"
    assert candidate["identity_comparison"]["status"] == "recorded"
    assert candidate["refinement"]["route_id"] != candidate["row"]["route_id"]
    saved = candidate["identity_comparison"]["record"]
    saved["original_hit"]["exact_smiles"] = "mutated view"
    view["identity_descriptor"]["record"]["code_sha256"]["identity.py"] = "changed"
    assert docs == before
    run = save(tmp_path, docs)
    raw = {p.name: p.read_bytes() for p in run.iterdir()}
    result = r.write_report(run)
    assert result["artifact_consistency"] == "consistent_recorded_bundle"
    output = (run / "review.html").read_text()
    for text in ("both_present", "scaffold", "c1ccccc1", "ring/linker graph", "False"):
        assert text in output
    assert all((run / name).read_bytes() == value for name, value in raw.items())


@pytest.mark.parametrize("state", ["absent", "valid", "null", "unsupported", "conflicting"])
@pytest.mark.parametrize("kind", ["absent", "unavailable", "recorded", "null"])
def test_identity_presence_matrix(state, kind):
    docs = identity_bundle()
    row, _, manifest = identity_parts(docs)
    if state == "absent":
        del manifest["identity_comparison_policy"]
    elif state == "null":
        manifest["identity_comparison_policy"] = None
    elif state == "unsupported":
        manifest["identity_comparison_policy"]["schema_version"] = 2
    elif state == "conflicting":
        manifest["code_sha256"]["common.py"] = "c" * 64
    if kind == "absent":
        del row["hit_identity_comparison"]
    elif kind == "null":
        row["hit_identity_comparison"] = None
    elif kind == "unavailable":
        unavailable(docs)
    view = r._assemble(docs, {}, [])
    expected = (
        (kind == "absent" and state == "absent")
        or kind == "unavailable"
        or (kind == "recorded" and state == "valid")
    )
    assert (view["artifact_consistency"] == "consistent_recorded_bundle") == expected
    assert "<saved & unavailable>" not in r.render_html(view)
    if kind == "unavailable":
        assert view["candidates"][0]["identity_comparison"]["status"] == "unavailable"
        assert "&lt;saved &amp; unavailable&gt;" in r.render_html(view)


@pytest.mark.parametrize("descriptor", ["absent", None, {}, {"schema_version": 99}, "valid"])
def test_identity_zero_selected_never_invents_row_failure(descriptor):
    docs = identity_bundle()
    docs["results.json"]["results"] = []
    docs["qualification.json"]["selected_route_ids"] = []
    docs["summary.json"].update(
        status="no_candidates", returned_products=0, accepted_products=0, qualified_products=0
    )
    if descriptor == "absent":
        del docs["manifest.json"]["identity_comparison_policy"]
    elif descriptor != "valid":
        docs["manifest.json"]["identity_comparison_policy"] = descriptor
    view = r._assemble(docs, {}, [])
    assert view["artifact_consistency"] == "consistent_recorded_bundle"
    assert not view["candidates"]
    assert "No selected candidates" in r.render_html(view)


@pytest.mark.parametrize(
    "reason",
    [
        "original_hit_not_retained",
        "deadline_before_hit_capture",
        "hit_capture_failed",
        "retained_identity_conflict",
        "final_product_not_retained",
        "retained_identity_binding_mismatch",
        "original_hit_classification_not_retained",
        "retained_identity_policy_mismatch",
    ],
)
def test_all_unavailable_reasons_are_not_chemical_rejections(reason):
    docs = identity_bundle()
    unavailable(docs, reason)
    assert identity_view(docs)["status"] == "unavailable"


@pytest.mark.parametrize("field", r._IDENTITY_RECORDED)
def test_unavailable_rejects_each_recorded_field_even_null(field):
    docs = identity_bundle()
    unavailable(docs)
    identity_parts(docs)[1][field] = None
    assert identity_view(docs)["issues"]


@pytest.mark.parametrize(
    "path,value",
    [
        (("schema_version",), True),
        (("schema_version",), 2),
        (("status",), "unknown"),
        (("reference_role",), "design"),
        (("original_hit_smiles",), PRODUCT),
        (("final_product_smiles",), HIT),
        (("original_hit",), []),
        (("final_product", "exact_smiles"), HIT),
        (("final_product", "family_smiles"), "other"),
        (("final_product", "scaffold_family_smiles"), "other"),
        (("original_hit", "policy_version"), "other"),
        (("final_product", "rdkit_version"), "2025.9.6"),
        (("original_hit", "has_scaffold"), 1),
        (("original_hit", "has_scaffold"), False),
        (("original_hit", "ring_skeleton_smiles"), ""),
        (("atom_aware_scaffold_changed",), 1),
        (("atom_aware_scaffold_changed",), False),
        (("generic_ring_linker_topology_changed",), True),
        (("generic_ring_linker_topology_changed",), "false"),
        (("scaffold_presence_change",), "lost"),
        (("change_class",), None),
        (("change_class",), "unknown"),
        (("change_class",), "state_only"),
        (("change_class",), "unchanged"),
        (("change_class",), "decoration"),
    ],
)
def test_identity_record_contradictions(path, value):
    docs = identity_bundle()
    node = identity_parts(docs)[1]
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value
    result = identity_view(docs)
    assert result["status"] == "inconsistent_metadata" and result["issues"]
    assert result["record"] is None  # Never display an inferred replacement.


@pytest.mark.parametrize(
    "field,value",
    [
        ("schema_version", True),
        ("schema_version", 2),
        ("policy_version", "future"),
        ("rdkit_version", ""),
        ("tautomer_limit", True),
        ("tautomer_limit", 0),
        ("transform_limit", 1.5),
        ("code_sha256", None),
        ("code_sha256", {"identity.py": "x" * 64, "common.py": "b" * 64}),
    ],
)
def test_descriptor_malformed_values_are_diagnostic(field, value):
    docs = identity_bundle()
    docs["manifest.json"]["identity_comparison_policy"][field] = value
    view = r._assemble(docs, {}, [])
    assert view["identity_descriptor"]["warnings"]
    assert view["candidates"][0]["identity_comparison"]["issues"]
    r.render_html(view)


@pytest.mark.parametrize("label", ["unchanged", "state_only", "decoration", "scaffold"])
@pytest.mark.parametrize("acyclic", [False, True])
def test_supplied_class_implications_not_fallback_classification(label, acyclic):
    docs = identity_bundle()
    docs["results.json"]["policy"]["constraints"]["mode"] = "unconstrained"
    row, record, manifest = identity_parts(docs)
    hit, product = record["original_hit"], record["final_product"]
    if label == "unchanged":
        product.update(copy.deepcopy(hit))
    elif label == "state_only":
        product.update(
            {k: hit[k] for k in ("family_smiles", "scaffold_family_smiles", "ring_skeleton_smiles")}
        )
    elif label == "decoration":
        product["scaffold_family_smiles"] = hit["scaffold_family_smiles"]
    if acyclic:
        for identity in (hit, product):
            identity.update(has_scaffold=False, scaffold_family_smiles="", ring_skeleton_smiles="")
    row.update(
        smiles=product["exact_smiles"],
        family_smiles=product["family_smiles"],
        scaffold_family_smiles=product["scaffold_family_smiles"],
    )
    record.update(
        final_product_smiles=row["smiles"],
        change_class=label,
        atom_aware_scaffold_changed=hit["scaffold_family_smiles"] != product["scaffold_family_smiles"],
        scaffold_presence_change="both_absent" if acyclic else "both_present",
    )
    result = identity_view(docs)
    assert (result["status"] == "recorded") == (not acyclic or label != "scaffold")
    if label != "scaffold" and result["status"] == "recorded":
        docs["results.json"]["policy"]["constraints"]["mode"] = "hop"
        assert identity_view(docs)["issues"]


@pytest.mark.parametrize("presence", ["gained", "lost", "both_present"])
def test_recorded_topology_may_change_or_scaffold_presence_change(presence):
    docs = identity_bundle()
    row, record, _ = identity_parts(docs)
    if presence != "both_present":
        identity = record["original_hit" if presence == "gained" else "final_product"]
        identity.update(has_scaffold=False, scaffold_family_smiles="", ring_skeleton_smiles="")
    else:
        record["final_product"]["ring_skeleton_smiles"] = "different generic graph"
    row["scaffold_family_smiles"] = record["final_product"]["scaffold_family_smiles"]
    record.update(scaffold_presence_change=presence, generic_ring_linker_topology_changed=True)
    assert identity_view(docs)["status"] == "recorded"


@pytest.mark.parametrize("equal_part", ["exact_smiles", "family_smiles", "scaffold_family_smiles"])
def test_nested_normalization_relationships_cannot_contradict(equal_part):
    docs = identity_bundle()
    row, record, _ = identity_parts(docs)
    record["original_hit"][equal_part] = record["final_product"][equal_part]
    # Keep exact common binding, so the implication itself must reject.
    if equal_part == "exact_smiles":
        docs["manifest.json"]["hit"] = record["original_hit_smiles"] = row["smiles"]
    if equal_part == "scaffold_family_smiles":
        record["original_hit"]["ring_skeleton_smiles"] = "contradictory graph"
    assert identity_view(docs)["issues"]


@pytest.mark.parametrize("broken", ["metric", "route", "count"])
def test_identity_failure_cannot_skip_existing_checks(broken):
    docs = identity_bundle()
    identity_parts(docs)[0]["hit_identity_comparison"] = None
    if broken == "metric":
        docs["refinement.json"]["candidates"][0]["trials"][0]["hit"]["shape"] = 0.01
    elif broken == "route":
        docs["routes.json"]["routes"][0]["route"]["nodes"][-1]["product"] = "wrong"
    else:
        docs["summary.json"]["returned_products"] = 99
    view = r._assemble(docs, {}, [])
    issues = view["issues"] + view["candidates"][0]["issues"]
    assert any("identity comparison" in issue for issue in issues)
    assert any("identity comparison" not in issue for issue in issues)
    r.render_html(view)


def test_mixed_rows_do_not_borrow_identity_and_duplicates_stay_invalid():
    docs = identity_bundle()
    second = copy.deepcopy(docs["results.json"]["results"][0])
    del second["hit_identity_comparison"]
    docs["results.json"]["results"].append(second)
    view = r._assemble(docs, {}, [])
    first, second_view = view["candidates"]
    assert first["identity_comparison"]["status"] == "recorded"
    assert second_view["identity_comparison"]["issues"]
    assert view["issues"]  # Existing selected-count/duplicate validation is not bypassed.
    first["identity_comparison"]["record"]["original_hit"]["family_smiles"] = "changed"
    assert (
        docs["results.json"]["results"][0]["hit_identity_comparison"]["original_hit"]["family_smiles"] == HIT
    )


def test_descriptor_source_change_is_visible_not_an_authenticity_claim():
    docs = identity_bundle()
    docs["manifest.json"]["sources_changed_since_import"] = ["identity.py"]
    view = r._assemble(docs, {}, [])
    assert view["artifact_consistency"] == "consistent_recorded_bundle"
    assert "sources changed since import" in r.render_html(view)


@pytest.mark.parametrize(
    "field",
    ["change_class", "original_hit", "final_product", "schema_version", "atom_aware_scaffold_changed"],
)
def test_missing_recorded_field_is_never_inferred(field):
    docs = identity_bundle()
    del identity_parts(docs)[1][field]
    result = identity_view(docs)
    assert result["status"] == "inconsistent_metadata" and result["record"] is None


@pytest.mark.parametrize(
    "field,value", [("reason_code", ""), ("reason", None), ("error_type", []), ("final_product_smiles", HIT)]
)
def test_unavailable_common_contract_still_checked(field, value):
    docs = identity_bundle()
    unavailable(docs)
    identity_parts(docs)[1][field] = value
    assert identity_view(docs)["issues"]
