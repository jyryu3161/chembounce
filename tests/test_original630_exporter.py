"""Synthetic original630 presentation fixtures only; no real results or chemistry."""

import csv
import importlib.util
from pathlib import Path
import shutil

import pytest

from benchmarks.diagnostics.generation_replay import export_figures as exporter

replay = exporter.replay
FIXTURE = importlib.util.spec_from_file_location(
    "synthetic_replay_fixture", Path(__file__).with_name("test_generation_statistical_replay.py")
)
fixture = importlib.util.module_from_spec(FIXTURE)
FIXTURE.loader.exec_module(fixture)


@pytest.fixture(scope="module")
def completed(tmp_path_factory):
    base = tmp_path_factory.mktemp("SYNTHETIC-original630")
    blobs = fixture.make_fixture()
    # Give one synthetic target no successful 3D scores, preserving cache joins.
    # Recompute/replay synthetic records only, so null CSV values have a valid origin.
    helper = replay.historical(blobs[replay.HELPER])
    code = replay.parse(blobs[replay.PROVENANCE])["code_sha256"]
    for split in replay.SPLITS:
        names = replay.paths(split)
        scored = replay.parse(blobs[names["common_generated_scores"]])
        for row in scored:
            if row["query_id"] in {"q0", "q1"}:
                row["all_3d_seeds_succeeded"] = False
                row["scores"] = {"bit_tanimoto": 0.2, "count_tanimoto": 0.25}
                row["trials"] = [
                    {"seed": seed, "status": "failed", "reason": "SYNTHETIC no 3D result", "seconds": 0.1}
                    for seed in replay.SEEDS
                ]
        blobs[names["common_generated_scores"]] = replay.encoded(scored)
        progress = replay.parse(blobs[names["common_progress"]])
        progress["failed_3d_occurrences"] = sum(not r["all_3d_seeds_succeeded"] for r in scored)
        blobs[names["common_progress"]] = replay.encoded(progress)
        runs, scores, targets, caps, metadata = replay.validate_split(split, blobs, code)
        blobs[names["comparison"]] = replay.encoded(
            {
                **metadata,
                **helper.analyze_runs(runs, scores, targets, caps),
                **replay.report_text(blobs[replay.HELPER]),
            }
        )
    fixture.refresh_review(blobs)
    bundle = base / "bundle"
    fixture.write_bundle(bundle, blobs)
    directory = base / "replay"
    replay.replay(bundle, replay.sha((bundle / "manifest.json").read_bytes()), directory)
    return directory


def anchor(directory):
    return replay.sha((directory / "validation.json").read_bytes())


def alter(directory, name, edit):
    path = directory / name
    value = replay.parse(path.read_bytes())
    edit(value)
    path.write_bytes(replay.encoded(value))
    if name != "validation.json":
        validation = replay.parse((directory / "validation.json").read_bytes())
        validation["recomputed_sha256"][name.split("-")[0]] = replay.sha(path.read_bytes())
        (directory / "validation.json").write_bytes(replay.encoded(validation))


@pytest.fixture
def mutable(completed, tmp_path):
    directory = tmp_path / "SYNTHETIC-replay"
    shutil.copytree(completed, directory)
    return directory


@pytest.fixture(scope="module")
def rendered(completed):
    matplotlib = pytest.importorskip("matplotlib")

    matplotlib.use("Agg", force=True)
    from matplotlib.figure import Figure

    output = completed.parent / "SYNTHETIC-presentation"
    captured = {}
    original = Figure.savefig

    def save(fig, path, **kwargs):
        captured[Path(path).name] = (
            [t.get_text() for t in fig.texts]
            + [
                t.get_text()
                for ax in fig.axes
                for t in [ax.title, ax.xaxis.label, *ax.texts, *ax.get_yticklabels()]
            ]
            + [t.get_text() for legend in fig.legends for t in legend.get_texts()]
        )
        return original(fig, path, **kwargs)

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(Figure, "savefig", save)
        manifest = exporter.export(completed, anchor(completed), output, synthetic_fixture=True)
    return output, manifest, captured


def rows(output, name):
    with (output / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def test_synthetic_figures_scope_labels_and_manifest(completed, rendered):
    output, manifest, captured = rendered
    assert manifest["synthetic_fixture"] is True
    assert manifest["original_source_commit"] == replay.COMMIT and manifest["original_scientific_files"] == 39
    assert len(manifest["files_sha256"]) == 10
    assert set(manifest["runtime"]) == {"python", "matplotlib", "numpy"}
    for name, expected in manifest["files_sha256"].items():
        assert replay.sha((output / name).read_bytes()) == expected
    for split in replay.SPLITS:
        report = replay.parse((completed / f"{split}-comparison.json").read_bytes())
        assert manifest["method_order"][split] == sorted(report["methods"])
        text = "\n".join(captured[f"{split}.png"])
        assert "SYNTHETIC FIXTURE" in text and "Original frozen generation study" in text
        assert "target-bootstrap 95% CI" in text and "NOT a CI" in text
        assert "Rule replay ≠ experimental synthesis or common-forward audit" in text
        assert f"{exporter.TARGETS[split]} independent targets" in text
        assert all(f"A={x['attempts']}" in text for x in report["methods"].values())
        assert (output / f"{split}.png").read_bytes().startswith(b"\x89PNG")
        assert (output / f"{split}.pdf").read_bytes().startswith(b"%PDF")


def assert_value(text, value):
    assert text == ("" if value is None else str(value))


def test_all_recorded_csv_metrics_and_paired_scope_are_preserved(completed, rendered):
    output, manifest, _ = rendered
    reports = {s: replay.parse((completed / f"{s}-comparison.json").read_bytes()) for s in replay.SPLITS}
    metrics = rows(output, "metrics.csv")
    missing = 0
    for row in metrics:
        report = reports[row["split"]]
        group = (
            report["methods"]
            if row["population"] == "realized_all_attempts"
            else report["complete_case_secondary"]["methods"]
        )
        item = group[row["method"]][row["level"]]
        if row["level"] != "target_macro":
            item = item[row["query_id"] or row["target_id"]]
        value = item[row["metric"]]
        if row["level"] == "target_macro":
            assert_value(row["ci95_low"], value["ci95"][0] if value["ci95"] else None)
            assert_value(row["ci95_high"], value["ci95"][1] if value["ci95"] else None)
            assert int(row["independent_targets"]) == value["independent_targets"]
            value = value["mean"]
        assert_value(row["value"], value)
        missing += value is None
    assert missing > 0
    assert len(metrics) == manifest["table_rows"]["metrics"]
    for split, report in reports.items():
        for method, item in report["methods"].items():
            actual = {
                r["metric"]
                for r in metrics
                if r["split"] == split
                and r["method"] == method
                and r["population"] == "realized_all_attempts"
            }
            assert actual == set(item["target_macro"])
    for row in rows(output, "paired_contrasts.csv"):
        report = reports[row["split"]]
        if row["population"] == "recorded_target_contrast":
            stat = report["paired_pipeline_minus_comparator"][row["comparator"]][row["metric"]]
            assert row["contrast_scope"] == exporter.PRIMARY_SCOPE
            assert "not like-for-like" in row["contrast_scope"]
        else:
            item = report["complete_case_secondary"]["paired_pipeline_minus_comparator_bounded_metrics"][
                row["comparator"]
            ][row["level"]]
            stat = (
                item[row["metric"]]
                if row["level"] == "target_macro"
                else {"mean": item[row["query_id"] or row["target_id"]][row["metric"]], "ci95": None}
            )
        assert_value(row["pipeline_minus_comparator"], stat["mean"])
        assert_value(row["ci95_low"], stat["ci95"][0] if stat["ci95"] else None)
        assert_value(row["ci95_high"], stat["ci95"][1] if stat["ci95"] else None)
    for row in rows(output, "infrastructure_bounds.csv"):
        key = (
            "methods" if row["population"] == "method_missing_outcome" else "paired_pipeline_minus_comparator"
        )
        item = reports[row["split"]]["infrastructure_missing_outcome_sensitivity"][key][
            row["method_or_comparator"]
        ][row["level"]]
        if row["level"] != "target_macro":
            item = item[row["query_id"] or row["target_id"]]
        assert_value(row["lower"], item[row["metric"]]["lower"])
        assert_value(row["upper"], item[row["metric"]]["upper"])
        assert "NOT confidence" in row["interval_type"]
    availability_rows = rows(output, "availability.csv")
    assert any("runtime_missing_attempts" in r["recorded_field"] for r in availability_rows)
    for row in availability_rows:
        report = reports[row["split"]]
        groups = {
            "realized_all_attempts": report["methods"],
            "complete_case_secondary": report["complete_case_secondary"]["methods"],
            "complete_case_paired_secondary": report["complete_case_secondary"][
                "paired_pipeline_minus_comparator_bounded_metrics"
            ],
            "method_missing_outcome": report["infrastructure_missing_outcome_sensitivity"]["methods"],
            "paired_pipeline_minus_comparator": report["infrastructure_missing_outcome_sensitivity"][
                "paired_pipeline_minus_comparator"
            ],
        }
        value = groups[row["population"]][row["method"]]
        for key in row["recorded_field"].split("."):
            value = value[key]
        assert_value(row["count"], value)
    projected = exporter.table_rows(reports)
    for name, values in projected.items():
        assert (output / (name + ".csv")).read_bytes() == exporter.csv_bytes(values)


@pytest.mark.parametrize(
    "failure",
    [
        "anchor",
        "missing_split",
        "split_hash",
        "helper",
        "tolerance",
        "incomplete",
        "query_count",
        "source_hash",
        "ci",
        "boolean_metric",
        "unsafe_provenance",
        "extra_file",
    ],
)
def test_gates_refuse_before_output_or_plot(mutable, tmp_path, monkeypatch, failure):
    expected = anchor(mutable)
    if failure == "anchor":
        expected = "0" * 64
    elif failure == "missing_split":
        (mutable / "target-comparison.json").unlink()
    elif failure == "extra_file":
        (mutable / "unexpected.txt").write_text("not permitted")
    elif failure == "split_hash":
        alter(mutable, "validation.json", lambda v: v["recomputed_sha256"].update(target="0" * 64))
    elif failure == "helper":
        alter(mutable, "validation.json", lambda v: v.update(historical_helper_sha256="0" * 64))
    elif failure == "tolerance":
        alter(mutable, "validation.json", lambda v: v.update(absolute_tolerance=1e-6))
    elif failure == "unsafe_provenance":
        alter(
            mutable,
            "validation.json",
            lambda v: v["input_before_after_sha256"].update({"parent-review-evidence/../escape": "0" * 64}),
        )
    else:

        def change(v):
            if failure == "incomplete":
                v["grid_audit"]["observed_runs"] -= 1
            elif failure == "query_count":
                v["methods"]["pipeline"]["per_query"].pop(next(iter(v["methods"]["pipeline"]["per_query"])))
            elif failure == "source_hash":
                v["source_sha256"][next(iter(v["source_sha256"]))] = "0" * 64
            elif failure == "ci":
                v["methods"]["pipeline"]["target_macro"][exporter.METRICS[0]]["ci95"] = [0.8, 0.2]
            else:
                v["methods"]["pipeline"]["target_macro"][exporter.METRICS[0]]["mean"] = True

        alter(mutable, "target-comparison.json", change)
    if failure not in {"anchor", "missing_split", "extra_file"}:
        expected = anchor(mutable)
    monkeypatch.setattr(exporter, "draw", lambda *a, **k: pytest.fail("Must reject before plotting"))
    output = tmp_path / "must-not-exist"
    with pytest.raises((ValueError, OSError)):
        exporter.export(mutable, expected, output)
    assert not output.exists()


@pytest.mark.parametrize("raw", [b'{"mean":NaN}', b'{"mean":1e999}', b'{"x":1,"x":2}', b"[]"])
def test_nonfinite_duplicate_and_malformed_json(mutable, tmp_path, raw):
    (mutable / "target-comparison.json").write_bytes(raw)
    alter(mutable, "validation.json", lambda v: v["recomputed_sha256"].update(target=replay.sha(raw)))
    with pytest.raises(ValueError):
        exporter.export(mutable, anchor(mutable), tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_symlinks_existing_output_and_input_change_preserve_evidence(mutable, tmp_path, monkeypatch):
    output = tmp_path / "existing"
    output.mkdir()
    (output / "keep").write_text("preserve")
    with pytest.raises(ValueError):
        exporter.export(mutable, anchor(mutable), output)
    assert (output / "keep").read_text() == "preserve"
    linked = tmp_path / "linked"
    linked.symlink_to(mutable, target_is_directory=True)
    with pytest.raises(ValueError):
        exporter.export(linked, anchor(mutable), tmp_path / "unused")
    dangling = tmp_path / "dangling-output"
    dangling.symlink_to(tmp_path / "absent")
    with pytest.raises(ValueError):
        exporter.export(mutable, anchor(mutable), dangling)
    assert dangling.is_symlink() and not (tmp_path / "absent").exists()
    original = mutable / "target-comparison.json"
    preserved = tmp_path / "preserved-report.json"
    original.rename(preserved)
    original.symlink_to(preserved)
    with pytest.raises(ValueError):
        exporter.export(mutable, anchor(mutable), tmp_path / "unused-file-link")
    original.unlink()
    preserved.rename(original)

    def mutate(*args, **kwargs):
        with (mutable / "time-comparison.json").open("ab") as stream:
            stream.write(b" ")

    monkeypatch.setattr(exporter, "draw", mutate)
    partial = tmp_path / "partial"
    with pytest.raises(ValueError):
        exporter.export(mutable, anchor(mutable), partial)
    assert partial.exists() and not (partial / "manifest.json").exists()


def test_missing_plot_value_remains_missing(completed, tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    from matplotlib.figure import Figure

    report = replay.parse((completed / "target-comparison.json").read_bytes())
    report["methods"]["pipeline"]["target_macro"][exporter.METRICS[0]] = {
        "mean": None,
        "ci95": None,
        "independent_targets": 0,
    }
    texts = []

    def capture(fig, *args, **kwargs):
        texts.extend(t.get_text() for ax in fig.axes for t in ax.texts)

    monkeypatch.setattr(Figure, "savefig", capture)
    exporter.draw("target", report, tmp_path, synthetic_fixture=True)
    assert "missing" in texts and "n=0" in texts
