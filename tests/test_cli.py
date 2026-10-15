import json
from pathlib import Path
import subprocess
import sys
import pytest

ROOT = Path(__file__).resolve().parents[1]


def invoke(tmp_path, *args):
    return subprocess.run(
        [sys.executable, str(ROOT / "chembounce.py"), "-o", str(tmp_path), *args],
        capture_output=True,
        text=True,
    )


def test_no_candidates_exit_code_and_header(tmp_path):
    result = invoke(tmp_path, "-i", "CCO", "--no-shape")
    assert result.returncode == 3
    assert len((tmp_path / "overall_result.txt").read_text().splitlines()) == 1
    assert (tmp_path / "resource_cost.json").exists()


@pytest.mark.parametrize("smiles", ["CC.O", "not-smiles"])
def test_invalid_input_is_nonzero_without_traceback(tmp_path, smiles):
    result = invoke(tmp_path, "-i", smiles)
    assert result.returncode == 2
    assert "Traceback" not in result.stderr
    assert json.loads((tmp_path / "run_summary.json").read_text())["status"] == "error"


def test_no_fingerprint_db_uses_custom_text_library(tmp_path):
    result = invoke(
        tmp_path,
        "-i",
        "CCc1ccncc1",
        "--scaffold-db",
        str(ROOT / "tests/fixtures/scaffolds.smi"),
        "--no-fingerprint-db",
        "--no-shape",
        "--search-threshold",
        "0",
        "-t",
        "0",
        "--overall_max_n",
        "2",
    )
    assert result.returncode == 0, result.stderr
    assert len((tmp_path / "overall_result.txt").read_text().splitlines()) == 3


def test_explicit_synllama_off_skips_even_an_unavailable_configuration(tmp_path):
    result = invoke(
        tmp_path,
        "-i",
        "CCc1ccncc1",
        "--scaffold-db",
        str(ROOT / "tests/fixtures/scaffolds.smi"),
        "--no-fingerprint-db",
        "--no-shape",
        "--search-threshold",
        "0",
        "-t",
        "0",
        "--overall_max_n",
        "2",
        "--no-synllama",
        "--synthesis-config",
        "missing-model-config.json",
    )
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "synthesis").exists()
    assert len((tmp_path / "overall_result.txt").read_text().splitlines()) == 3


def test_explicit_synllama_requires_config_before_hopping(tmp_path):
    result = invoke(tmp_path, "-i", "CCO", "--synllama")
    assert result.returncode == 2
    assert "requires --synthesis-config" in result.stderr
    assert not (tmp_path / "overall_result.txt").exists()


@pytest.mark.parametrize("args", [("--refine-top-n", "10"), ("--include-run", "prior-run")])
def test_synthesis_overrides_require_config_before_hopping(tmp_path, args):
    result = invoke(tmp_path, "-i", "CCO", *args)
    assert result.returncode == 2
    assert "require --synthesis-config" in result.stderr
    assert not (tmp_path / "overall_result.txt").exists()


def test_existing_synthesis_outputs_are_protected_before_hopping_overwrite(tmp_path):
    old = tmp_path / "overall_result.txt"
    old.write_text("previous hopping results\n")
    synthesis = tmp_path / "synthesis"
    synthesis.mkdir()
    (synthesis / "results.tsv").write_text("previous synthesis results\n")
    result = invoke(
        tmp_path,
        "-i",
        "CCO",
        "--overwrite",
        "--synthesis-config",
        str(ROOT / "examples/synthesis/demo-config.json"),
    )
    assert result.returncode == 2
    assert "fresh output directory" in result.stderr
    assert old.read_text() == "previous hopping results\n"
    assert (synthesis / "results.tsv").read_text() == "previous synthesis results\n"


@pytest.mark.parametrize(
    "setting",
    [
        {"search": {"max_steps": 0}},
        {"conformers": {"num_conformers": 0}},
        {"conformers": {"energy_window": "5"}},
        {"conformers": {"seed": True}},
    ],
)
def test_invalid_synthesis_settings_fail_before_hopping(tmp_path, setting):
    config = json.loads((ROOT / "examples/synthesis/demo-config.json").read_text())
    config.update(setting)
    path = tmp_path / "invalid.json"
    path.write_text(json.dumps(config))
    result = invoke(tmp_path, "-i", "CCO", "--synthesis-config", str(path))
    assert result.returncode == 2
    assert "Traceback" not in result.stderr
    assert not (tmp_path / "overall_result.txt").exists()
