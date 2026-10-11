"""Checked setup and execution for the stable notebook (host standard library only)."""

import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess
import tarfile
import time
import uuid
import zipfile


COUNT_OPTIONS = ("overall_max_n", "scaffold_top_n", "cand_max_n__rplc")
THRESHOLD_OPTIONS = tuple(
    f"{name}_{bound}"
    for name in ("qed", "sa", "logp", "mw", "h_donor", "h_acceptor")
    for bound in ("min", "max")
)
DEFAULTS = {
    "core_smiles": "",
    "frag_max_n": 100,
    "tanimoto_threshold": 0.5,
    "lipinski": True,
    "use_custom_library": False,
    "custom_lib_file_name": "",
    "fragments": [],
    "replace_scaffold_files": [],
    **dict.fromkeys(COUNT_OPTIONS + THRESHOLD_OPTIONS),
    "cand_max_n__rplc": 10,
}


def setup_environment(checkout, environment, uv):
    """Repair incomplete installs on rerun; never alter the notebook's interpreter."""
    checkout, environment = Path(checkout).resolve(), Path(environment).resolve()
    python = environment / "bin/python"
    if not python.is_file():
        subprocess.run(
            [str(uv), "venv", "--python", "3.10", "--managed-python", str(environment)],
            check=True,
        )
    subprocess.run(
        [
            str(uv),
            "--config-file",
            str(checkout / "uv.toml"),
            "pip",
            "install",
            "--python",
            str(python),
            "-r",
            str(checkout / "requirements-colab.txt"),
        ],
        check=True,
    )
    subprocess.run([str(uv), "pip", "check", "--python", str(python)], check=True)
    probe_environment(python, checkout)
    return python


def probe_environment(python, checkout):
    """Import the actual CLI dependencies before reporting setup success."""
    probe = (
        "import json, sys, importlib.metadata as m; import chembounce, utils; "
        "assert sys.version_info[:2] == (3, 10), 'ChemBounce 1.x needs Python 3.10'; "
        "print(json.dumps({'python': sys.version, 'packages': "
        "{n: m.version(n) for n in ('rdkit', 'numpy', 'pandas', 'scipy', 'oddt', "
        "'MolVS', 'scaffoldgraph')}}))"
    )
    result = subprocess.run(
        [str(python), "-c", probe],
        cwd=checkout,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode:
        raise RuntimeError(f"ChemBounce environment check failed:\n{result.stdout}")
    # Some legacy dependencies print import notices before the JSON line.
    runtime = json.loads(result.stdout.strip().splitlines()[-1])
    print(f"ChemBounce runtime ready: {runtime}")
    return runtime


def ensure_data(checkout, python):
    """Check files used by 1.x; preserve an incomplete download before retrying."""
    checkout = Path(checkout).resolve()
    data = checkout / "data"

    def valid():
        text = data / "Scaffolds_processed.txt"
        fingerprints = data / "scaffold_fingerprints_mw250.npz"
        return (
            text.is_file()
            and text.stat().st_size > 0
            and fingerprints.is_file()
            and zipfile.is_zipfile(fingerprints)
        )

    if valid():
        return
    if data.exists():
        backup = data.with_name(f"data.incomplete-{uuid.uuid4().hex[:8]}")
        data.rename(backup)
        print(f"Preserved incomplete reference data at {backup}")
    subprocess.run(
        [
            str(python),
            str(checkout / "scripts/download_data.py"),
            "--output",
            str(data),
        ],
        check=True,
    )
    if not valid():
        raise RuntimeError(
            f"Reference data are incomplete: {data}. Rerun setup to retry."
        )


def build_command(python, checkout, output, settings):
    checkout = Path(checkout).resolve()
    options = {**DEFAULTS, **settings}
    command = [
        str(python),
        "-u",
        str(checkout / "chembounce.py"),
        "-o",
        str(output),
        "-i",
        options["input_smiles"],
        "--tanimoto_threshold",
        str(options["tanimoto_threshold"]),
        "--frag_max_n",
        str(options["frag_max_n"]),
    ]
    if options["core_smiles"] and options["core_smiles"] != "None":
        command += ["--core_smiles", options["core_smiles"]]
    if not options["lipinski"]:
        command.append("--wo_lipinski")
    for name in COUNT_OPTIONS + THRESHOLD_OPTIONS:
        value = options[name]
        if value is not None:
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise ValueError(f"{name} must be a finite number or None")
            if name in COUNT_OPTIONS:
                if value < 1 or int(value) != value:
                    raise ValueError(f"{name} must be a positive integer or None")
                value = int(value)
            command += [f"--{name}", str(value)]
    database = (
        Path(str(options["custom_lib_file_name"]) + ".npz")
        if options["use_custom_library"]
        else Path("data/scaffold_fingerprints_mw250.npz")
    )
    database = (checkout / database).resolve()
    if not database.is_file():
        raise FileNotFoundError(
            f"Fingerprint database not found: {database}. Run setup/custom-library cell."
        )
    command += ["--fingerprint-db", str(database)]
    fragments = list(options["fragments"])
    files = list(options["replace_scaffold_files"])
    if files and len(files) != len(fragments):
        raise ValueError(
            "Supply one replacement-file entry per fragment (empty entries use the default library)."
        )
    for fragment in fragments:
        command += ["--fragments", fragment]
    for filename in files:
        path = (checkout / filename).resolve() if filename else None
        if path is not None and not path.is_file():
            raise FileNotFoundError(f"Replacement scaffold file not found: {path}")
        command += ["--replace_scaffold_files", str(path) if path else ""]
    return command


def validate_results(output):
    """Distinguish genuine zero-candidate results from a failed or partial run."""
    output = Path(output)
    summary = json.loads((output / "run_summary.json").read_text())
    with (output / "overall_result.txt").open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or "Final structure" not in reader.fieldnames:
            raise ValueError("Missing result-table header")
        count = 0
        for row in reader:
            if not row.get("Final structure"):
                raise ValueError("Incomplete row in result table")
            count += 1
    expected_status = "success" if count else "no_candidates"
    if (
        summary.get("status") != expected_status
        or summary.get("counts", {}).get("accepted") != count
    ):
        raise ValueError("Result table and run summary disagree")
    return count, summary


def run_and_archive(python, checkout, output, settings):
    """Execute in a fresh directory; archive only checked results, with diagnostics."""
    checkout = Path(checkout).resolve()
    output = (checkout / output).resolve()
    archive = output.with_name(output.name + ".tar.gz")
    if archive.exists() or (
        output.exists() and (not output.is_dir() or any(output.iterdir()))
    ):
        suffix = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        output = output.with_name(f"{output.name}-{suffix}-{uuid.uuid4().hex[:8]}")
        archive = output.with_name(output.name + ".tar.gz")
    output.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output}")
    started = time.monotonic()
    manifest = {
        "status": "running",
        "output": str(output),
        "cwd": str(checkout),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path = output / "notebook_run.json"

    def save_manifest():
        manifest["elapsed_seconds"] = round(time.monotonic() - started, 3)
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    save_manifest()
    try:
        command = build_command(python, checkout, output, settings)
        manifest["command"] = command
        save_manifest()
        with (output / "run.log").open("w") as log:
            with subprocess.Popen(
                command,
                cwd=checkout,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=1,
            ) as process:
                try:
                    for line in process.stdout:
                        log.write(line)
                        log.flush()
                        print(line, end="", flush=True)
                    manifest["returncode"] = process.wait()
                except BaseException:
                    process.kill()
                    process.wait()
                    raise
        if manifest["returncode"] != 0:
            raise RuntimeError(
                f"ChemBounce exited with status {manifest['returncode']}"
            )
        count, summary = validate_results(output)
        manifest.update(status=summary["status"], candidate_count=count)
        save_manifest()
        # Write to a temporary name so an interrupted compression is not a valid download.
        partial = archive.with_name(archive.name + ".partial")
        try:
            with tarfile.open(partial, "w:gz") as bundle:
                bundle.add(output, arcname=output.name)
            partial.replace(archive)
        finally:
            partial.unlink(missing_ok=True)
    except BaseException as error:
        manifest.update(status="failed", error=str(error))
        save_manifest()
        with (output / "run.log").open("a") as log:
            log.write(f"\nNotebook run failed: {error}\n")
        print(
            f"Run failed. No result archive created. See {output / 'run.log'} and {manifest_path}."
        )
        raise
    print(f"Completed: {count} candidates. {summary['message']}")
    print(f"Checked result archive: {archive}")
    return archive
