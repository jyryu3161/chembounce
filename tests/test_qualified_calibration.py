"""Synthetic-only checks of frozen development calibration and strict records."""

import copy
from dataclasses import asdict
from fractions import Fraction
import json
import math

import pytest
from rdkit import Chem

from benchmarks.research import qualified_calibration as calibration
from chembounce_synthesis import conformers


def row(ident, *, target="A", query="q1", scaffold="s1", label="active", value=0.7, topology=True):
    return {
        "pair_id": ident,
        "target_id": target,
        "assay_id": target + query,
        "context_id": target + query,
        "query_id": query,
        "control_scaffold_family": scaffold,
        "label": label,
        "features": [value] * 3 if value is not None else None,
        "topology_change": topology,
    }


def six_targets():
    return [
        row(f"{t}-{label}-{i}", target=t, label=label, value=0.7 if label == "active" else 0.2)
        for t in "ABCDEF"
        for label in ("active", "inactive")
        for i in range(10)
    ]


def test_weighting_and_primary_rates_use_different_frozen_hierarchies():
    rows = [row(str(i), value=0.1) for i in range(3)]
    rows += [
        row("a4", scaffold="s2", value=0.2),
        row("a5", query="q2", value=0.3),
        row("b1", target="B", value=0.4),
        row("failed", value=None),
    ]
    _, report = calibration.weighted_grid(rows)
    assert report["row_weights"] == {
        "0": "1/24",
        "1": "1/24",
        "2": "1/24",
        "a4": "1/8",
        "a5": "1/4",
        "b1": "1/2",
    }
    assert report["quantile_values"]["hit_tanimoto"]["3/20"] == 0.2
    assert report["supported_active_rows"] == 6
    assert report["excluded_active_rows"] == 1
    predictions = {r["pair_id"]: r["features"] is not None for r in rows}
    rates = calibration.rates(rows, predictions, "active")
    assert rates["targets"] == {"A": Fraction(9, 10), "B": Fraction(1)}
    assert rates["macro"] == Fraction(19, 20)
    sensitivity = calibration.rates(rows, predictions, "active", scaffold=True)
    assert sensitivity["macro"] != rates["macro"]


def test_inverse_cdf_groups_ties_and_uses_observed_endpoint_without_rounding():
    low, high = 0.1, math.nextafter(0.1, 1)
    points = [(low, Fraction(1, 4)), (low, Fraction(1, 4)), (high, Fraction(1, 2))]
    assert calibration.inverse_cdf(points, Fraction(1, 2)) == low
    assert calibration.inverse_cdf(points, Fraction(500001, 1000000)) == high
    assert calibration.inverse_cdf(points, Fraction(1)) == high
    assert calibration.accepted({"features": [low, low, low]}, [low, low, low])
    assert not calibration.accepted({"features": [low, low, low]}, [high, low, low])


def test_grid_deduplicates_and_discards_nonpositive_cutoffs():
    rows = [row("zero", value=0), row("positive", value=0.5)]
    grid, diagnostics = calibration.weighted_grid(rows)
    assert grid == []
    assert diagnostics["quantile_values"]["hit_tanimoto"]["1/5"] == 0
    assert calibration.fit_thresholds(rows)["fit_status"] == "infeasible"


def test_joint_acceptance_not_product_of_marginal_retentions(monkeypatch):
    rows = [row(f"a{i}") for i in range(10)] + [row("inactive", label="inactive", value=0.1)]
    for i in range(3):
        rows[i]["features"][i] = 0.2
    monkeypatch.setattr(calibration, "weighted_grid", lambda rows: ([(0.5, 0.5, 0.5)], {}))
    fitted = calibration.fit_thresholds(rows)
    assert fitted["fit_status"] == "infeasible"
    assert fitted["candidates"][0]["active"]["macro_rational"] == "7/10"


def test_exact_constraint_boundary_and_tie_order(monkeypatch):
    rows = [row(f"a{i}", value=0.7 if i < 8 else None) for i in range(10)]
    rows += [row("inactive", label="inactive", value=0.01)]
    monkeypatch.setattr(calibration, "weighted_grid", lambda rows: ([(0.2, 0.1, 0.3), (0.1, 0.2, 0.3)], {}))
    fitted = calibration.fit_thresholds(rows)
    assert fitted["cutoff"] == [0.1, 0.2, 0.3]
    assert all(c["feasible"] and c["active"]["macro_rational"] == "4/5" for c in fitted["candidates"])


def test_unsupported_topology_never_satisfies_constraint():
    rows = [row("active", topology=False), row("inactive", label="inactive", topology=False)]
    fitted = calibration.fit_thresholds(rows)
    assert fitted["fit_status"] == "infeasible"
    assert fitted["candidates"][0]["topology_active"]["macro"] is None
    assert fitted["candidates"][0]["topology_active"]["unsupported_targets"] == 1


def test_loto_uses_only_training_rows_and_failed_rows_remain(monkeypatch):
    rows = six_targets()
    rows[0]["features"] = None
    original = calibration.fit_thresholds
    calls = []

    def inspect(train):
        calls.append({r["target_id"] for r in train})
        return original(train)

    monkeypatch.setattr(calibration, "fit_thresholds", inspect)
    result = calibration.fit_development(rows)
    assert calls == [set("ABCDEF") - {t} for t in "ABCDEF"] + [set("ABCDEF")]
    assert result["minimum_scores"] == dict.fromkeys(calibration.METRICS, 0.7)
    assert len(result["oof"]["predictions"]) == 120
    assert result["oof"]["predictions"][rows[0]["pair_id"]] is False
    assert result["oof"]["report"]["active_retention"]["macro_rational"] == "59/60"
    assert result["oof"]["bootstrap"]["draws"] == result["final"]["bootstrap"]["draws"]


def test_infeasible_folds_abstain_without_claiming_geometric_success():
    rows = six_targets()
    for r in rows:
        r["topology_change"] = False
    result = calibration.fit_development(rows)
    assert result["minimum_scores"] is None
    assert all(f["fit_status"] == "infeasible" for f in result["folds"])
    report = result["oof"]["report"]
    assert report["inactive_operational_rejection"]["macro"] == 1
    assert report["active_retention"]["macro"] == 0
    assert report["geometric_discrimination"]["status"] == "not_assessed"
    assert report["fit_coverage"]["pairs"] == 0
    assert report["rejection_counts"]["no_fit"] == len(rows)


def test_bootstrap_fixed_predictions_uses_inverse_empirical_cdf():
    import numpy as np

    rows = [row(t, target=t) for t in "ABCDEF"]
    predictions = {t: t in "AB" for t in "ABCDEF"}
    result = calibration.bootstrap(rows, predictions)
    expected = np.random.Generator(np.random.PCG64(20261016)).integers(0, 6, size=(2000, 6)).tolist()
    assert result["draws"] == expected
    samples = sorted(Fraction(sum(i < 2 for i in draw), 6) for draw in expected)
    assert result["metrics"]["active_retention"]["interval"] == [float(samples[49]), float(samples[1949])]
    assert result["metrics"]["inactive_operational_rejection"]["valid_replicates"] == 0


def pair():
    result = {**row("synthetic-pair"), "query_smiles": "CCO", "control_smiles": "CCCO"}
    result.pop("features")
    return result


def synthetic_result(product, hit, design, settings):
    # Explicit coordinates only: no embedding, optimization, real scorer, or model.
    pose = Chem.MolFromSmiles(product)
    conf = Chem.Conformer(pose.GetNumAtoms())
    conf.Set3D(True)
    conf.SetId(7)
    for index in range(pose.GetNumAtoms()):
        conf.SetAtomPosition(index, (index * 1.123456789, index * 0.7, index * 0.3))
    pose.AddConformer(conf, assignId=False)
    joint = math.sqrt(0.6 * 0.7)
    overlay = {
        "shape": 0.6,
        "pharmacophore": 0.7,
        "reference_conformer_id": 3,
        "reference_relative_energy": 1.0,
        "joint": joint,
        "reference_adjusted_joint": joint * math.exp(-settings.strain_penalty),
    }
    info = {
        "candidate_smiles": product,
        "hit_smiles": hit,
        "design_smiles": design,
        "settings": asdict(settings),
        "scoring_objective": settings.scoring_objective,
        "candidate_conformer_id": 7,
        "candidate_energy": 10.0,
        "candidate_relative_energy": 2.0,
        "score": overlay["reference_adjusted_joint"] * math.exp(-settings.strain_penalty * 2),
        "hit": dict(overlay),
        "design": dict(overlay),
        "force_field": "synthetic_fixture",
        "ensemble_sizes": {"candidate": 2, "hit": 2, "design": 2},
    }
    return info, Chem.Mol(pose), Chem.Mol(pose)


def test_seed_records_keep_exact_metadata_both_poses_and_full_precision(monkeypatch):
    monkeypatch.setattr(conformers, "score_candidate", synthetic_result)
    scored = calibration.score_pair(pair())
    assert [t["seed"] for t in scored["trials"]] == [42, 43, 44]
    assert scored["all_3d_seeds_succeeded"]
    assert scored["features"][1:] == [0.6, 0.7]
    assert scored["convergence_status"] == "not_assessed"
    for trial in scored["trials"]:
        assert trial["info"]["ensemble_sizes"]["hit"] == 2
        for saved in trial["poses"].values():
            restored = calibration.restore_pose(saved)
            assert restored.GetConformer(7).GetAtomPosition(1).x == 1.123456789
    calibration.validate_scored_row(scored, pair())


def test_one_seed_failure_never_supplies_partial_features(monkeypatch):
    def scorer(product, hit, design, settings):
        if settings.seed == 43:
            raise ValueError("synthetic failed seed")
        return synthetic_result(product, hit, design, settings)

    monkeypatch.setattr(conformers, "score_candidate", scorer)
    scored = calibration.score_pair(pair())
    assert [t["status"] for t in scored["trials"]] == ["success", "failed", "success"]
    assert scored["features"] is None
    assert scored["status"] == "failed"
    assert scored["trials"][1]["error_type"] == "ValueError"
    assert scored["trials"][1]["settings"]["seed"] == 43
    calibration.validate_scored_row(scored, pair())


def test_all_seed_failure_is_explicit_and_preserves_denominators(monkeypatch):
    def fail(*args):
        raise ValueError("synthetic all-seed failure")

    monkeypatch.setattr(conformers, "score_candidate", fail)
    scored = calibration.score_pair(pair())
    assert scored["status"] == "failed"
    assert not scored["all_3d_seeds_succeeded"]
    assert scored["features"] is None
    assert len(scored["trials"]) == 3
    assert all(t["status"] == "failed" for t in scored["trials"])
    calibration.validate_scored_row(scored, pair())


@pytest.mark.parametrize(
    "smiles", ["C[C@H](O)[13CH2][NH3+]", "C/C=C/[C@@H](N)C(=O)[O-]", "[2H]C([2H])([2H])[C@H](O)C"]
)
def test_pose_roundtrip_preserves_real_synthetic_stereo_isotope_charge_and_coordinates(smiles):
    from rdkit.Chem import AllChem
    from chembounce_synthesis.common import canonical

    pose = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(pose, randomSeed=42) == 0
    pose = Chem.RemoveHs(pose)
    restored = calibration.restore_pose(calibration.pose_record(pose, 0))
    assert canonical(Chem.MolToSmiles(restored)) == canonical(smiles)
    assert restored.GetConformer().GetPositions().tolist() == pose.GetConformer().GetPositions().tolist()


@pytest.mark.parametrize(
    "defect", ["seed", "identity", "energy", "objective", "coordinate", "feature", "count"]
)
def test_mutated_seed_metadata_and_summaries_fail_closed(monkeypatch, defect):
    monkeypatch.setattr(conformers, "score_candidate", synthetic_result)
    scored = calibration.score_pair(pair())
    trial = scored["trials"][0]
    if defect == "seed":
        trial["seed"] = 43
    elif defect == "identity":
        trial["info"]["hit_smiles"] = "C"
    elif defect == "energy":
        trial["info"]["candidate_relative_energy"] = 9.0
    elif defect == "objective":
        trial["info"]["scoring_objective"] = "joint-v1"
    elif defect == "coordinate":
        trial["poses"]["hit"]["coordinates"][0][0] = float("nan")
    elif defect == "feature":
        scored["features"][0] = 0.987
    else:
        trial["info"]["ensemble_sizes"]["hit"] = 9
    with pytest.raises(ValueError):
        calibration.validate_scored_row(scored, pair())


def test_historical_parity_is_diagnostic_not_an_eligibility_source(monkeypatch):
    monkeypatch.setattr(conformers, "score_candidate", synthetic_result)
    scored = calibration.score_pair(pair())
    previous = {
        **pair(),
        "trials": [
            {
                "seed": t["seed"],
                "status": "success",
                "metrics": {"score": t["info"]["score"], "hit_shape": 0.6, "hit_pharmacophore": 0.7},
            }
            for t in scored["trials"]
        ],
    }
    report = calibration.historical_parity([scored], [previous])
    assert report["compared_successes"] == 3
    assert report["max_absolute_delta"] == 0
    failed = calibration.failed_row(pair(), "synthetic worker death", "worker_execution")
    report = calibration.historical_parity([failed], [previous])
    assert report["compared_successes"] == 0
    assert report["status_changes"] == 3
    assert failed["features"] is None


def test_identity_inventory_requires_exact_original_labels_and_structures():
    original = pair()
    manifest = {"expected_pairs": 1, "expected_requested_trials": 3, "pairs": [original]}
    calibration.verify_inventory([original], manifest)
    changed = {**original, "label": "inactive"}
    with pytest.raises(ValueError, match="differs"):
        calibration.verify_inventory([changed], manifest)


def test_procedure_digest_and_semantic_drift_rejected_before_data_work(tmp_path):
    path = tmp_path / "procedure.json"
    calibration.save(path, {"calibration_spec": calibration.SPEC})
    with pytest.raises(ValueError, match="SHA-256"):
        calibration.load_procedure(path, "0" * 64)
    changed = copy.deepcopy(calibration.SPEC)
    changed["quantile_levels"].append("1/2")
    calibration.save(path, {"calibration_spec": changed})
    with pytest.raises(ValueError, match="calibration_spec"):
        calibration.load_procedure(path, calibration.file_hash(path))


def test_worker_failure_ledger_contains_every_requested_seed():
    result = calibration.failed_row(pair(), "synthetic process failure", "worker_execution")
    assert [t["seed"] for t in result["trials"]] == list(calibration.SEEDS)
    assert all(t["status"] == "unreported" and t["reason"] for t in result["trials"])
    assert not result["all_3d_seeds_succeeded"]


def test_spec_command_has_no_data_access_or_fitting(monkeypatch, capsys):
    def unexpected(*args, **kwargs):
        pytest.fail("spec must not score or fit")

    monkeypatch.setattr(calibration, "rescore", unexpected)
    monkeypatch.setattr(calibration, "fit", unexpected)
    calibration.main(["spec"])
    assert json.loads(capsys.readouterr().out) == calibration.SPEC


@pytest.fixture
def synthetic_artifact_inputs(tmp_path, monkeypatch):
    from benchmarks.research import assay_matched

    pairs = [
        {**pair(), "pair_id": f"{target}-{label}", "target_id": target, "label": label}
        for target in ("A", "B")
        for label in ("active", "inactive")
    ]
    panel = tmp_path / "panel.json"
    calibration.save(panel, {"scope": "synthetic_fixture"})
    identity = tmp_path / "identities.json"
    calibration.save(identity, {"expected_pairs": 4, "expected_requested_trials": 12, "pairs": pairs})
    previous = tmp_path / "old.json"
    historical = []
    for item in pairs:
        trials = []
        for seed in calibration.SEEDS:
            settings = conformers.ConformerSettings(**calibration.SETTINGS, seed=seed)
            info = synthetic_result(
                item["control_smiles"], item["query_smiles"], item["query_smiles"], settings
            )[0]
            trials.append(
                {
                    "seed": seed,
                    "status": "success",
                    "metrics": {"score": info["score"], "hit_shape": 0.6, "hit_pharmacophore": 0.7},
                }
            )
        historical.append({**item, "trials": trials})
    calibration.save(previous, historical)
    current = calibration.runtime()
    procedure = {
        "calibration_spec": calibration.SPEC,
        "panel_sha256": calibration.file_hash(panel),
        "expected_counts": calibration.counts(pairs),
        "expected_identity_subgroups": {
            "normalized_scaffold_hops": 4,
            "topology_changes": 4,
            "same_topology_scaffold_changes": 0,
        },
        "runtime_for_scoring": {
            "interpreter": current["executable"],
            "python": current["python"],
            "numpy": current["packages"]["numpy"],
            "rdkit": current["packages"]["rdkit"],
        },
        "identity_manifest": {"path": str(identity), "sha256": calibration.file_hash(identity)},
        "historical_scores": {"path": str(previous), "sha256": calibration.file_hash(previous)},
    }
    frozen = tmp_path / "procedure.json"
    calibration.save(frozen, procedure)
    # The public CLI correctly rejects replacing the197-pair panel. Only this
    # synthetic orchestration test substitutes a smaller verified input ledger.
    monkeypatch.setattr(calibration, "load_procedure", lambda *args: procedure)
    monkeypatch.setattr(assay_matched, "validate_panel", lambda *args, **kwargs: None)
    monkeypatch.setattr(calibration, "panel_pairs", lambda *args: pairs)
    monkeypatch.setattr(conformers, "score_candidate", synthetic_result)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        monkeypatch.setenv(name, "1")
    return panel, frozen, calibration.file_hash(frozen), pairs


def test_separate_scoring_and_fit_commands_roundtrip_complete_synthetic_artifacts(
    synthetic_artifact_inputs, tmp_path
):
    panel, procedure, digest, pairs = synthetic_artifact_inputs
    scores = tmp_path / "scored"
    rows = calibration.rescore(panel, procedure, digest, scores)
    assert len(rows) == 4 and sum(len(r["trials"]) for r in rows) == 12
    assert json.loads((scores / "parity.json").read_text())["max_absolute_delta"] == 0
    assert not (scores / "fit.json").exists()
    fitted = calibration.fit(scores, procedure, digest, tmp_path / "fitted")
    assert fitted["provenance"]["counts"] == calibration.counts(pairs)
    assert (
        fitted["oof"]["report"]["geometric_discrimination"]["status"]
        == "assessed_conditional_on_fit_and_complete_computation"
    )
    candidate = json.loads((tmp_path / "fitted" / "preset-candidate.json").read_text())
    assert candidate["minimum_scores"] == fitted["minimum_scores"]
    assert candidate["status"] == "development_fit_requires_independent_review"
    assert candidate["fit_sha256"] == calibration.file_hash(tmp_path / "fitted" / "fit.json")
    with pytest.raises(FileExistsError):
        calibration.rescore(panel, procedure, digest, scores)


@pytest.mark.parametrize("failure", ["initialization", "submission", "execution"])
def test_process_failures_preserve_full_synthetic_pair_seed_ledger(
    synthetic_artifact_inputs, tmp_path, monkeypatch, failure
):
    import concurrent.futures

    class BrokenExecutor:
        def __init__(self, **kwargs):
            assert kwargs["mp_context"].get_start_method() == "spawn"
            assert kwargs["max_workers"] == 2
            if failure == "initialization":
                raise RuntimeError("synthetic initialization failure")

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def submit(self, *args):
            if failure == "submission":
                raise RuntimeError("synthetic submission failure")
            future = concurrent.futures.Future()
            future.set_exception(RuntimeError("synthetic worker failure"))
            return future

    monkeypatch.setattr(calibration.concurrent.futures, "ProcessPoolExecutor", BrokenExecutor)
    panel, procedure, digest, _ = synthetic_artifact_inputs
    output = tmp_path / "failures"
    rows = calibration.rescore(panel, procedure, digest, output, workers=2)
    assert len(rows) == 4 and sum(len(r["trials"]) for r in rows) == 12
    assert all(r["status"] == "error" and r["features"] is None for r in rows)
    assert all(t["stage"] == "worker_" + failure for r in rows for t in r["trials"])
    assert json.loads((output / "progress.json").read_text())["complete"]


def test_mutated_artifact_and_runtime_fail_before_fit(synthetic_artifact_inputs, tmp_path, monkeypatch):
    panel, procedure, digest, _ = synthetic_artifact_inputs
    output = tmp_path / "scores"
    calibration.rescore(panel, procedure, digest, output)
    with (output / "scores.json").open("a") as stream:
        stream.write(" ")
    with pytest.raises(ValueError, match="integrity"):
        calibration.fit(output, procedure, digest, tmp_path / "fit")
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    with pytest.raises(ValueError, match="CUDA_VISIBLE_DEVICES"):
        calibration.rescore(panel, procedure, digest, tmp_path / "not-created")
    assert not (tmp_path / "not-created").exists()


def test_changed_pharmacophore_feature_asset_blocks_fitting(synthetic_artifact_inputs, tmp_path, monkeypatch):
    panel, procedure, digest, _ = synthetic_artifact_inputs
    output = tmp_path / "scores"
    calibration.rescore(panel, procedure, digest, output)
    changed = copy.deepcopy(calibration.runtime())
    changed["pharmacophore_feature_definition"]["sha256"] = "0" * 64
    monkeypatch.setattr(calibration, "runtime", lambda: changed)
    with pytest.raises(ValueError, match="scientific runtime"):
        calibration.fit(output, procedure, digest, tmp_path / "fit")
