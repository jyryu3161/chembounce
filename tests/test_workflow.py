import json
from pathlib import Path
import pytest
from rdkit import Chem
import chembounce as cb
import utils
from conftest import ADP

LIB = Path(__file__).parent / "fixtures/scaffolds.smi"


def run(tmp_path, **kwargs):
    return cb.chembounce(ADP, output_dir=tmp_path, scaffold_db_file=LIB, calculate_shape=False, **kwargs)[0]


def test_adp_end_to_end_with_diagnostics(tmp_path):
    result = run(tmp_path, candidate_thresholds={}, overall_max_n=4, frag_max_n=3)
    assert 0 < len(result) <= 4
    assert not result["Final structure"].duplicated().any()
    assert result["Tanimoto Similarity"].ge(0.5).all()
    assert result["Standardized final structure"].notna().all()
    report = json.loads((tmp_path / "run_summary.json").read_text())
    assert report["counts"]["accepted"] == len(result)
    assert report["status"] == "success"
    for fragment_file in tmp_path.glob("fragment_result_*.txt"):
        assert len(fragment_file.read_text().splitlines()) <= 4
    manifest = json.loads((tmp_path / "run_manifest.json").read_text())
    assert manifest["sources"]
    assert manifest["parameters"]["effective_limits"]["overall"] == 4


def test_property_exclusion_is_explicit(tmp_path):
    result = run(tmp_path, candidate_thresholds={"mw": (0, 1)}, overall_max_n=4)
    assert result.empty
    report = json.loads((tmp_path / "run_summary.json").read_text())
    assert report["status"] == "no_candidates"
    assert report["counts"]["property_rejected"] > 0
    assert "filters" in report["message"]


def test_core_is_enforced_on_candidates(tmp_path):
    core = "Nc1ncnc2c1nc[nH]2"
    result = run(tmp_path, core_smiles=core, candidate_thresholds={}, tanimoto_threshold=0.2, frag_max_n=6)
    assert len(result)
    assert all(
        Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmiles(core), useChirality=True)
        for s in result["Final structure"]
    )
    report = json.loads((tmp_path / "run_summary.json").read_text())
    assert report["counts"]["core_rejected"] > 0


@pytest.mark.parametrize(
    "options",
    [
        {"core_smiles": "broken"},
        {"core_smiles": "ClCl"},
        {"tanimoto_threshold": float("nan")},
        {"frag_max_n": 0},
        {"overall_max_n": -1},
        {"candidate_thresholds": {"mw": (500, 100)}},
        {"fragments": ["CCCl"]},
        {"scaffold_db_file": "/does/not/exist"},
    ],
)
def test_invalid_requests_fail_with_report(tmp_path, options):
    with pytest.raises((ValueError, OSError, utils.SmilesInputError)):
        cb.chembounce(ADP, output_dir=tmp_path, calculate_shape=False, **options)
    assert json.loads((tmp_path / "run_summary.json").read_text())["status"] == "error"


def test_output_requires_explicit_overwrite(tmp_path):
    run(tmp_path, candidate_thresholds={}, overall_max_n=1)
    before = (tmp_path / "overall_result.txt").read_bytes()
    with pytest.raises(FileExistsError):
        run(tmp_path)
    assert (tmp_path / "overall_result.txt").read_bytes() == before
    run(tmp_path, candidate_thresholds={}, overall_max_n=1, overwrite=True)


def test_acyclic_no_database_needed(tmp_path):
    results, _ = cb.chembounce("CCO", output_dir=tmp_path, calculate_shape=False)
    assert results.empty
    assert "No ring" in json.loads((tmp_path / "run_summary.json").read_text())["message"]


def test_shared_replacement_file_and_api_no_filters(tmp_path):
    library = tmp_path / "replacements.smi"
    library.write_text("c1ccncc1\nC1CCNC1\n")
    result, _ = cb.chembounce(
        ADP,
        output_dir=tmp_path / "out",
        fragments=["C1CCOC1", "c1ncc2nc[nH]c2n1"],
        replace_scaffold_files=[library],
        candidate_thresholds={},
        calculate_shape=False,
        overall_max_n=4,
        tanimoto_threshold=0.2,
    )
    assert len(result) == 4


def test_default_filters_adp_warning(tmp_path):
    run(tmp_path, overall_max_n=2)
    report = json.loads((tmp_path / "run_summary.json").read_text())
    assert any("n_hdonor" in reason for reason in report["input_filter_failures"])
    assert any("n_hacceptor" in reason for reason in report["input_filter_failures"])


def test_two_fingerprint_databases_do_not_share_stale_cache(tmp_path):
    from generate_fingerprints import generate_fingerprints

    first, second = tmp_path / "first.smi", tmp_path / "second.smi"
    first.write_text("c1ccncc1\n")
    second.write_text("C1CCCCC1\n")
    for source, name in ((first, "first"), (second, "second")):
        database = tmp_path / (name + ".npz")
        generate_fingerprints(source, database)
        cb.chembounce(
            "CCc1ccccc1",
            output_dir=tmp_path / name,
            fingerprint_db_file=database,
            _search_scf_thr_=0,
            tanimoto_threshold=0,
            calculate_shape=False,
            candidate_thresholds={},
        )
    one = (tmp_path / "first/replace_scaffold_list.fragment_0001.tsv").read_text()
    two = (tmp_path / "second/replace_scaffold_list.fragment_0001.tsv").read_text()
    assert "c1ccncc1" in one and "C1CCCCC1" not in one
    assert "C1CCCCC1" in two and "c1ccncc1" not in two


def test_mismatched_database_source_is_rejected(tmp_path):
    from generate_fingerprints import generate_fingerprints

    source = tmp_path / "source.smi"
    source.write_text("c1ccncc1\n")
    database = tmp_path / "fp.npz"
    generate_fingerprints(source, database)
    source.write_text("C1CCCCC1\n")
    with pytest.raises(ValueError, match="source hash"):
        cb.chembounce(
            "CCc1ccccc1", output_dir=tmp_path / "out", scaffold_db_file=source, fingerprint_db_file=database
        )


def test_shape_failure_keeps_candidate_and_explains_nan(tmp_path, monkeypatch):
    def fail(*args):
        raise ValueError("embedding failed for test")

    monkeypatch.setattr(utils, "calc_electron_shape", fail)
    result, _ = cb.chembounce(
        ADP, output_dir=tmp_path, scaffold_db_file=LIB, candidate_thresholds={}, overall_max_n=1
    )
    assert len(result) == 1
    assert result["Electron shape Similarity"].isna().all()
    report = json.loads((tmp_path / "run_summary.json").read_text())
    assert report["counts"]["shape_failed"] == 1
    assert "embedding failed" in report["shape_errors"][0]["reason"]


@pytest.mark.parametrize("low_mem", [False, True])
def test_text_library_skips_and_reports_unusable_records(tmp_path, low_mem):
    library = tmp_path / "mixed.smi"
    library.write_text("CC*\nCC.CC\nbad\nc1ccncc1\n")
    result, _ = cb.chembounce(
        "CCc1ccccc1",
        scaffold_db_file=library,
        output_dir=tmp_path / "out",
        use_fingerprint_db=False,
        low_mem=low_mem,
        candidate_thresholds={},
        calculate_shape=False,
        _search_scf_thr_=0,
        tanimoto_threshold=0,
        overall_max_n=2,
    )
    assert len(result) > 0
    report = json.loads((tmp_path / "out/run_summary.json").read_text())
    assert report["text_library"]["rejected_records"] == 3
    assert report["text_library"]["valid_records"] == 1
