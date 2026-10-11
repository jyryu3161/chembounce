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
