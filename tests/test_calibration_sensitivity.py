"""Synthetic-only invariants for the frozen reporting addendum."""

import copy
import json
from pathlib import Path

import pytest

from benchmarks.research import calibration_sensitivity as sensitivity
from benchmarks.research import qualified_calibration as calibration


def row(
    ident,
    *,
    target="A",
    label="active",
    bit=0.6,
    shape=(0.4, 0.6, 0.8),
    pharm=(0.8, 0.4, 0.6),
    failed=None,
    query="q",
    scaffold="s",
    topology=True,
):
    trials = [
        {
            "seed": seed,
            "status": "failed" if seed == failed else "success",
            "info": {"hit": {"shape": shape[i], "pharmacophore": pharm[i]}},
        }
        for i, seed in enumerate(calibration.SEEDS)
    ]
    result = {
        "pair_id": ident,
        "target_id": target,
        "assay_id": target + query,
        "context_id": target + query,
        "query_id": query,
        "control_scaffold_family": scaffold,
        "label": label,
        "topology_change": topology,
        "trials": trials,
        "hit_tanimoto": bit,
        "all_3d_seeds_succeeded": failed is None,
        "status": "success" if failed is None else "failed",
    }
    result["features"] = calibration.features(result) if failed is None else None
    return result


def locked_fit(rows, *, final=None, folds=None):
    """Hand-specified locked rules; deliberately never invokes a fitter."""
    if final is None:
        final = [0.5, 0.5, 0.5]
    targets = sorted({r["target_id"] for r in rows})
    folds = dict.fromkeys(targets, final) if folds is None else folds
    cutoff = lambda c: {"cutoff": c, "fit_status": "fitted" if c is not None else "infeasible"}
    blocks, oof = [], {}
    for target in targets:
        test = [r for r in rows if r["target_id"] == target]
        predicted = {r["pair_id"]: calibration.accepted(r, folds[target]) for r in test}
        blocks.append(
            {
                "test_target": target,
                "training_targets": sorted(set(targets) - {target}),
                **cutoff(folds[target]),
                "test_predictions": predicted,
                "test_report": calibration.prediction_report(
                    test, predicted, fitted_targets=[target] if folds[target] is not None else []
                ),
            }
        )
        oof.update(predicted)
    predicted = {r["pair_id"]: calibration.accepted(r, final) for r in rows}
    return {
        "scope": calibration.SPEC["scope"],
        "folds": blocks,
        "final": {
            **cutoff(final),
            "predictions": predicted,
            "report": calibration.prediction_report(rows, predicted, fitted_targets=targets),
        },
        "oof": {
            "predictions": oof,
            "report": calibration.prediction_report(
                rows, oof, fitted_targets=[t for t, c in folds.items() if c is not None]
            ),
        },
        "minimum_scores": dict(zip(calibration.METRICS, final)),
    }


def variant(application, name):
    return next(v for v in application["variants"] if v["id"] == name)


def test_variants_never_fit_score_or_mutate_inputs(monkeypatch):
    rows = [row("a"), row("b", target="B", label="inactive")]
    fit = locked_fit(rows)
    before = copy.deepcopy((rows, fit))

    def prohibited(*args, **kwargs):
        pytest.fail("Reporting attempted fitting/scoring")

    for name in ("fit", "fit_thresholds", "fit_development", "rescore", "score_pair", "bootstrap"):
        monkeypatch.setattr(calibration, name, prohibited)
    report = sensitivity.build_report(rows, fit)
    assert len(report["final"]["variants"]) == 16
    assert report["refitting_performed"] is False
    assert (rows, fit) == before
    assert report["final"]["baseline"]["accepted_count"] == 0
    assert variant(report["final"], "three_seed_median")["accepted_pair_ids"] == ["a", "b"]


def test_omit_seed_changes_both_3d_metrics_and_preserves_failed_seed_eligibility():
    rows = [
        row("pass", shape=(0.4, 0.7, 0.8), pharm=(0.3, 0.6, 0.9)),
        row("failed", shape=(0.4, 0.7, 0.8), pharm=(0.3, 0.6, 0.9), failed=42),
    ]
    result = sensitivity.build_report(rows, locked_fit(rows))
    v = variant(result["final"], "omit_seed42")
    assert v["accepted_pair_ids"] == ["pass"]
    assert v["rejected_to_accepted_pair_ids"] == ["pass"]
    assert v["report"]["active_retention"]["macro_rational"] == "1/2"
    for v in result["final"]["variants"]:
        assert not v["predictions"]["failed"]
        assert v["report"]["computational_failure_pairs"] == 1
    assert result["computational_failure_pair_ids"] == ["failed"]


def test_omitting_one_seed_does_not_select_best_shape_and_pharm_independently():
    rows = [row("crossed", shape=(0.3, 0.8, 0.8), pharm=(0.8, 0.3, 0.8))]
    result = sensitivity.build_report(rows, locked_fit(rows))["final"]
    assert variant(result, "three_seed_median")["accepted_count"] == 1
    for seed in calibration.SEEDS:
        assert variant(result, f"omit_seed{seed}")["accepted_count"] == 0


@pytest.mark.parametrize("index,metric", enumerate(calibration.METRICS))
def test_cutoff_perturbs_exactly_one_coordinate_and_uses_original_minima(index, metric):
    rows = [
        row("base", shape=(0.5,) * 3, pharm=(0.5,) * 3, bit=0.5),
        row("near", shape=(0.46,) * 3, pharm=(0.46,) * 3, bit=0.46),
    ]
    result = sensitivity.build_report(rows, locked_fit(rows))["final"]
    for factor in sensitivity.FACTORS:
        case = variant(result, f"{metric}_times_{factor}")
        expected = [0.5] * 3
        expected[index] *= factor
        assert case["cutoffs_by_target"]["A"] == expected
        assert case["cutoff_provenance"]["A"]["requested"] == 0.5 * factor
        assert case["predictions"]["near"] is False
        assert case["predictions"]["base"] is (factor < 1)
        if factor > 1:
            assert case["accepted_to_rejected_pair_ids"] == ["base"]


def test_boundary_gain_full_float_precision_and_fixed_binary_feature():
    rows = [row("gain", bit=0.5 * 0.9, shape=(0.8,) * 3, pharm=(0.8,) * 3)]
    result = sensitivity.build_report(rows, locked_fit(rows))["final"]
    case = variant(result, "hit_tanimoto_times_0.9")
    assert case["rejected_to_accepted_pair_ids"] == ["gain"]
    assert case["cutoff_provenance"]["A"]["requested_hex"] == (0.5 * 0.9).hex()
    assert variant(result, "three_seed_median")["accepted_count"] == 0


def test_clipping_and_duplicate_cases_are_retained_and_marked():
    rows = [row("a", bit=1, shape=(1,) * 3, pharm=(1,) * 3)]
    report = sensitivity.build_report(rows, locked_fit(rows, final=[0.99, 1, 1]))["final"]
    first, second = (variant(report, f"hit_tanimoto_times_{f}") for f in (1.05, 1.1))
    assert first["cutoff_provenance"]["A"]["clipped"] is True
    assert first["cutoffs_by_target"]["A"] == [1, 1, 1]
    assert second["duplicate_of"] == first["id"]
    unchanged = variant(report, "hit_shape_times_1.1")
    assert unchanged["duplicate_of"] == "baseline"
    assert unchanged["cutoff_provenance"]["A"]["unchanged"] is True
    assert len(report["variants"]) == 16


def test_locked_oof_rules_and_no_fit_folds_do_not_use_final_rule():
    rows = [
        row("a", shape=(0.8,) * 3, pharm=(0.8,) * 3),
        row("b", target="B", label="inactive", shape=(0.8,) * 3, pharm=(0.8,) * 3),
    ]
    fit = locked_fit(rows, folds={"A": [0.8, 0.5, 0.5], "B": None})
    result = sensitivity.build_report(rows, fit)
    assert result["final"]["baseline"]["accepted_count"] == 2
    for case in result["locked_oof"]["variants"]:
        assert case["accepted_count"] == 0
        assert case["report"]["fit_coverage"]["no_fit_pairs"] == 1
    for case in result["folds"][1]["variants"]:
        assert case["report"]["geometric_discrimination"]["status"] == "not_assessed"
        assert case["report"]["inactive_operational_rejection"]["macro_rational"] == "1"


def test_final_no_fit_remains_no_fit():
    rows = [row("a")]
    fit = locked_fit(rows, folds={"A": None})
    fit["final"].update(cutoff=None, fit_status="infeasible")
    fit["final"]["report"] = calibration.prediction_report(rows, {"a": False}, fitted_targets=[])
    fit["minimum_scores"] = None
    result = sensitivity.build_report(rows, fit)
    for case in result["final"]["variants"]:
        assert case["accepted_count"] == 0
        assert case["report"]["geometric_discrimination"]["status"] == "not_assessed"


def test_rates_keep_full_query_target_topology_and_scaffold_denominators():
    rows = [
        row("a1", shape=(0.4, 0.8, 0.9), pharm=(0.4, 0.8, 0.9)),
        row("a2", failed=42),
        row("a3", query="q2", topology=False, shape=(0.8,) * 3, pharm=(0.8,) * 3),
        row("b", target="B", shape=(0.8,) * 3, pharm=(0.8,) * 3),
    ]
    case = variant(sensitivity.build_report(rows, locked_fit(rows))["final"], "three_seed_median")
    assert case["report"]["active_retention"]["macro_rational"] == "7/8"
    assert case["report"]["topology_active_retention"]["macro_rational"] == "3/4"
    assert case["report"]["scaffold_sensitivity"]["active"]["macro_rational"] == "7/8"


@pytest.mark.parametrize(
    "mutation",
    [
        lambda rows, fit: rows[0]["trials"].pop(),
        lambda rows, fit: rows[0]["trials"].reverse(),
        lambda rows, fit: rows[0]["features"].__setitem__(1, 0.9),
        lambda rows, fit: rows[0]["trials"][0]["info"]["hit"].__setitem__("shape", float("nan")),
        lambda rows, fit: fit["final"]["predictions"].__setitem__("a", True),
        lambda rows, fit: fit["oof"]["predictions"].__setitem__("a", True),
        lambda rows, fit: fit["folds"][0].__setitem__("training_targets", ["A"]),
        lambda rows, fit: fit["folds"].append(copy.deepcopy(fit["folds"][0])),
        lambda rows, fit: fit["final"].__setitem__("cutoff", [0, 0.5, 0.5]),
        lambda rows, fit: fit["final"]["report"]["fit_coverage"].__setitem__("no_fit_pairs", 10),
    ],
)
def test_semantic_mutations_rejected(mutation):
    rows = [row("a")]
    fit = locked_fit(rows)
    mutation(rows, fit)
    with pytest.raises(ValueError):
        sensitivity.build_report(rows, fit)


@pytest.fixture
def bound_inputs(tmp_path, monkeypatch):
    rows = [row("a"), row("b", target="B", label="inactive")]
    pairs = [
        {
            k: r[k]
            for k in (
                "pair_id",
                *calibration.QUERY_FIELDS,
                "label",
                "control_scaffold_family",
                "topology_change",
            )
        }
        for r in rows
    ]
    source, fitted, code = [tmp_path / name for name in ("scores", "fit", "code")]
    for path in (source, fitted, code):
        path.mkdir()
    for name in calibration.CODE_FILES:
        path = code / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(
            Path(calibration.__file__).read_bytes()
            if name.endswith("qualified_calibration.py")
            else b"synthetic source"
        )
    hashes = {name: calibration.file_hash(code / name) for name in calibration.CODE_FILES}
    panel = {"synthetic": True}
    calibration.save(source / "panel.json", panel)
    calibration.save(source / "pairs.json", pairs)
    identity = {"pairs": pairs, "expected_pairs": len(rows), "expected_requested_trials": len(rows) * 3}
    calibration.save(source / "identity_manifest.json", identity)
    procedure = {
        "panel_sha256": calibration.file_hash(source / "panel.json"),
        "expected_counts": calibration.counts(rows),
    }
    calibration.save(source / "procedure.json", procedure)
    procedure_sha = calibration.file_hash(source / "procedure.json")
    runtime = {
        "python": "synthetic",
        "packages": {},
        "pharmacophore_feature_definition": {"sha256": "synthetic"},
    }
    protocol = {
        "calibration_spec": calibration.SPEC,
        "procedure_sha256": procedure_sha,
        "pairs_sha256": calibration.file_hash(source / "pairs.json"),
        "source_sha256": hashes,
        "runtime": runtime,
    }
    for name, value in (
        ("protocol", protocol),
        ("scores", rows),
        ("parity", {}),
        ("progress", {"complete": True, "expected_pairs": 2, "reported_pairs": 2, "expected_trials": 6}),
    ):
        calibration.save(source / f"{name}.json", value)
    artifacts = {name: calibration.file_hash(source / name) for name in sensitivity.SCORING_FILES}
    calibration.save(source / "artifacts.json", artifacts)
    fit = locked_fit(rows)
    fit["provenance"] = {
        "procedure_sha256": procedure_sha,
        "scoring_artifacts": artifacts,
        "runtime": runtime,
        "source_sha256": hashes,
        "counts": calibration.counts(rows),
    }
    calibration.save(fitted / "fit.json", fit)
    fit_sha = calibration.file_hash(fitted / "fit.json")
    calibration.save(fitted / "artifacts.json", {"fit.json": fit_sha})
    addendum = tmp_path / "addendum.json"
    calibration.save(addendum, {"schema": sensitivity.SCHEMA, "primary_procedure": {"sha256": procedure_sha}})
    addendum_sha = calibration.file_hash(addendum)
    monkeypatch.setattr(sensitivity, "ADDENDUM_SHA256", addendum_sha)

    def load(path, digest):
        assert calibration.file_hash(path) == digest
        return procedure

    monkeypatch.setattr(calibration, "load_procedure", load)
    monkeypatch.setattr(calibration, "verified_reference", lambda *args: source / "identity_manifest.json")
    return dict(
        scores_dir=source,
        fit_dir=fitted,
        fit_sha256=fit_sha,
        procedure_path=source / "procedure.json",
        procedure_sha256=procedure_sha,
        addendum_path=addendum,
        addendum_sha256=addendum_sha,
        scoring_source_root=code,
    )


def test_cli_integrity_chain_and_new_output_only(bound_inputs, tmp_path):
    kwargs = {**bound_inputs, "output": tmp_path / "report"}
    result = sensitivity.report(**kwargs)
    assert result["provenance"]["fit_sha256"] == bound_inputs["fit_sha256"]
    assert json.loads((kwargs["output"] / "sensitivity.json").read_text()) == result
    assert set(json.loads((kwargs["output"] / "artifacts.json").read_text())) == {"sensitivity.json"}
    with pytest.raises(FileExistsError):
        sensitivity.report(**kwargs)


@pytest.mark.parametrize("filename", sensitivity.SCORING_FILES)
def test_each_scoring_artifact_is_hash_bound(bound_inputs, filename):
    path = bound_inputs["scores_dir"] / filename
    path.write_text(path.read_text() + " ")
    with pytest.raises((ValueError, AssertionError)):
        sensitivity.load_inputs(**bound_inputs)


@pytest.mark.parametrize("key", ["fit_sha256", "addendum_sha256", "procedure_sha256"])
def test_external_digests_cannot_be_replaced_with_self_consistent_local_manifest(bound_inputs, key):
    bound_inputs[key] = "0" * 64
    with pytest.raises((ValueError, AssertionError)):
        sensitivity.load_inputs(**bound_inputs)


def test_original_source_cannot_change_while_reporter_is_added(bound_inputs):
    (bound_inputs["scoring_source_root"] / "chembounce_synthesis/conformers.py").write_text("changed")
    with pytest.raises(ValueError, match="Immutable scoring source"):
        sensitivity.load_inputs(**bound_inputs)


def test_fit_file_hash_cannot_be_updated_only_in_local_manifest(bound_inputs):
    path = bound_inputs["fit_dir"] / "fit.json"
    path.write_text(path.read_text() + " ")
    calibration.save(path.parent / "artifacts.json", {"fit.json": calibration.file_hash(path)})
    with pytest.raises(ValueError, match="externally approved fit"):
        sensitivity.load_inputs(**bound_inputs)
