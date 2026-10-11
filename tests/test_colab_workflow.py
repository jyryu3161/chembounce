"""Exercise notebook orchestration with real child processes and archive contents."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import colab_support as colab  # noqa: E402
from scripts import download_data  # noqa: E402


@pytest.fixture
def checkout(tmp_path):
    repo = tmp_path / "checkout with spaces"
    repo.mkdir()
    (repo / "data").mkdir()
    (repo / "data/scaffold_fingerprints_mw250.npz").write_bytes(b"placeholder")
    (repo / "chembounce.py").write_text("""
import json, pathlib, sys
args = sys.argv[1:]
out = pathlib.Path(args[args.index('-o') + 1])
mode = args[args.index('-i') + 1]
if mode == 'import-error':
    import intentionally_missing_chembounce_dependency
if mode == 'no-output':
    sys.exit(0)
count = 0 if mode == 'zero-candidates' else 1
(out / 'overall_result.txt').write_text('Final structure\\n' + ('CCO\\n' if count else ''))
(out / 'run_summary.json').write_text(json.dumps({
    'status': 'success' if count else 'no_candidates',
    'counts': {'accepted': count}, 'message': 'Result diagnostics'}))
(out / 'received_args.json').write_text(json.dumps(args))
print('CLI completed')
if mode == 'failed-after-output':
    sys.exit(7)
""")
    return repo


@pytest.mark.parametrize("mode", ["import-error", "no-output", "failed-after-output"])
def test_failed_runs_never_become_downloads(checkout, mode):
    with pytest.raises((RuntimeError, FileNotFoundError)):
        colab.run_and_archive(
            sys.executable, checkout, "result", {"input_smiles": mode}
        )
    assert not (checkout / "result.tar.gz").exists()
    report = json.loads((checkout / "result/notebook_run.json").read_text())
    assert report["status"] == "failed"
    log = (checkout / "result/run.log").read_text()
    if mode == "import-error":
        assert "ModuleNotFoundError" in log
    elif mode == "failed-after-output":
        assert report["returncode"] == 7


def test_missing_interpreter_is_logged_without_archive(checkout):
    with pytest.raises(FileNotFoundError):
        colab.run_and_archive(
            checkout / "missing/python", checkout, "result", {"input_smiles": "CCO"}
        )
    assert not (checkout / "result.tar.gz").exists()
    assert "missing/python" in (checkout / "result/run.log").read_text()


@pytest.mark.parametrize("mode,count", [("CCO", 1), ("zero-candidates", 0)])
def test_valid_results_archive_includes_data_and_diagnostics(checkout, mode, count):
    archive = colab.run_and_archive(
        sys.executable, checkout, "result", {"input_smiles": mode}
    )
    with tarfile.open(archive) as bundle:
        names = bundle.getnames()
        assert all(
            f"result/{name}" in names
            for name in (
                "overall_result.txt",
                "run_summary.json",
                "run.log",
                "notebook_run.json",
            )
        )
        report = json.load(bundle.extractfile("result/notebook_run.json"))
    assert report["candidate_count"] == count
    assert report["status"] == ("success" if count else "no_candidates")


def test_rerun_preserves_previous_results_even_when_new_run_fails(checkout):
    original = colab.run_and_archive(
        sys.executable, checkout, "result", {"input_smiles": "CCO"}
    )
    old_bytes = original.read_bytes()
    with pytest.raises(RuntimeError):
        colab.run_and_archive(
            sys.executable, checkout, "result", {"input_smiles": "import-error"}
        )
    assert original.read_bytes() == old_bytes
    assert len(list(checkout.glob("*.tar.gz"))) == 1
    assert len(list(checkout.glob("*/notebook_run.json"))) == 2


def test_literal_paths_core_and_fragments_reach_child_process(checkout):
    library = checkout / "custom library $(touch UNEXPECTED).txt.npz"
    library.write_bytes(b"placeholder")
    colab.run_and_archive(
        sys.executable,
        checkout,
        "output $(touch UNEXPECTED)",
        {
            "input_smiles": "CCO",
            "core_smiles": "CO",
            "fragments": ["CO", "CC"],
            "use_custom_library": True,
            "custom_lib_file_name": str(library)[:-4],
        },
    )
    args = json.loads(
        (checkout / "output $(touch UNEXPECTED)/received_args.json").read_text()
    )
    assert args[args.index("--fingerprint-db") + 1] == str(library)
    assert args[args.index("--core_smiles") + 1] == "CO"
    assert args.count("--fragments") == 2
    assert not (checkout / "UNEXPECTED").exists()


def test_missing_custom_library_does_not_silently_use_default(checkout):
    with pytest.raises(FileNotFoundError, match="Fingerprint database"):
        colab.run_and_archive(
            sys.executable,
            checkout,
            "result",
            {
                "input_smiles": "CCO",
                "use_custom_library": True,
                "custom_lib_file_name": "absent.txt",
            },
        )
    assert not (checkout / "result.tar.gz").exists()


def test_inconsistent_result_count_is_rejected(tmp_path):
    (tmp_path / "overall_result.txt").write_text("Final structure\nCCO\n")
    (tmp_path / "run_summary.json").write_text(
        json.dumps({"status": "success", "counts": {"accepted": 2}})
    )
    with pytest.raises(ValueError, match="disagree"):
        colab.validate_results(tmp_path)


def test_notebook_runs_with_all_optional_cells_skipped(checkout):
    notebook = json.loads((ROOT / "chembounce_colab.ipynb").read_text())
    sources = {c["metadata"]["id"]: "".join(c["source"]) for c in notebook["cells"]}
    downloaded = []

    class Files:
        download = staticmethod(downloaded.append)

    scope = {
        "DEFAULTS": colab.DEFAULTS,
        "checkout": checkout,
        "chembounce_python": sys.executable,
        "run_and_archive": colab.run_and_archive,
        "files": Files,
        "os": __import__("os"),
    }
    # Execute the same default initialization as the setup cell, without network installation.
    defaults = sources["stable-setup"].split("# Initialize every optional setting")[1]
    exec(
        compile(
            "# Initialize every optional setting" + defaults, "notebook-setup", "exec"
        ),
        scope,
    )
    exec(compile(sources["stable-input"], "notebook-input", "exec"), scope)
    exec(compile(sources["stable-run"], "notebook-run", "exec"), scope)
    assert len(downloaded) == 1 and Path(downloaded[0]).is_file()
    scope["input_smiles"] = "import-error"
    with pytest.raises(RuntimeError):
        exec(compile(sources["stable-run"], "notebook-run", "exec"), scope)
    assert len(downloaded) == 1  # Failed rerun cannot download the previous archive.
    assert scope["archive"] is None


def test_reference_data_check_does_not_require_obsolete_pickle(checkout, monkeypatch):
    (checkout / "data/Scaffolds_processed.txt").write_text("c1ccccc1\n")
    with zipfile.ZipFile(
        checkout / "data/scaffold_fingerprints_mw250.npz", "w"
    ) as archive:
        archive.writestr("fingerprints.npy", b"fixture")

    def unexpected(*args, **kwargs):
        pytest.fail("Complete reference data must not be downloaded again")

    monkeypatch.setattr(subprocess, "run", unexpected)
    colab.ensure_data(checkout, sys.executable)


@pytest.mark.parametrize("unsafe", [False, True])
def test_data_installer_checks_checksum_and_archive_paths(
    tmp_path, monkeypatch, unsafe
):
    archive = tmp_path / "reference.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.writestr(
            "../outside.txt" if unsafe else "data/Scaffolds_processed.txt", "CCO\n"
        )
    monkeypatch.setattr(download_data, "URL", archive.as_uri())
    monkeypatch.setattr(
        download_data, "MD5", hashlib.md5(archive.read_bytes()).hexdigest()
    )
    destination = tmp_path / "installed"
    if unsafe:
        with pytest.raises(ValueError, match="Unsafe"):
            download_data.install(destination)
        assert not destination.exists() and not (tmp_path / "outside.txt").exists()
    else:
        download_data.install(destination)
        assert (destination / "Scaffolds_processed.txt").read_text() == "CCO\n"
        assert (destination / "download_manifest.json").is_file()
    monkeypatch.setattr(download_data, "MD5", "0" * 32)
    with pytest.raises(ValueError, match="checksum"):
        download_data.install(tmp_path / "bad-download")
    assert not (tmp_path / "bad-download").exists()
