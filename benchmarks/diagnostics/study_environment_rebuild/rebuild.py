"""Bounded same-host installation diagnostic. Running this downloads and installs."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PINS_SHA256 = "2bec9937d4df0195c57dc308471929c1773a4cf1fa066ae8b9c82e646ea13de6"
CANCELLED = False
INSTALL_ORDER = {
    "research": ["research-pypi.txt"],
    "forward-cpu": ["forward-cpu-pypi.txt", "forward-cpu-pytorch-cpu.txt"],
    "crem": ["crem-pypi.txt"],
    "mol2mol": ["mol2mol-pypi.txt", "mol2mol-pytorch-cu121.txt", "mol2mol-nvidia.txt", "mol2mol-openeye.txt"],
}

METADATA_CHECK = r'''
import importlib.metadata as metadata, json, pathlib, re, sys
from urllib.parse import unquote, urlparse
profile = json.load(open(sys.argv[1]))['profiles'][sys.argv[2]]
normalize = lambda name: re.sub(r'[-_.]+', '-', name).lower()
rows = list(metadata.distributions())
actual = {normalize(d.metadata['Name']): d.version for d in rows}
assert sys.version_info[:3] == (3, 10, 12)
assert actual == profile['installed_packages'], {k:(profile['installed_packages'].get(k),actual.get(k)) for k in profile['installed_packages'].keys() | actual.keys() if profile['installed_packages'].get(k)!=actual.get(k)}
assert len(rows) == len(actual), 'Unexpected duplicated installed metadata'
prefix = pathlib.Path(sys.prefix).resolve()
assert all(pathlib.Path(d._path).resolve().is_relative_to(prefix) for d in rows)
source = profile['source_installation']
direct = {}
for d in rows:
    raw = d.read_text('direct_url.json')
    if raw:
        direct[normalize(d.metadata['Name'])] = json.loads(raw)
assert set(direct) == {source['distribution']}
entry = direct[source['distribution']]
assert bool(entry.get('dir_info',{}).get('editable',False)) == source['editable']
url = urlparse(entry['url'])
assert url.scheme == 'file' and url.netloc in ('','localhost')
assert pathlib.Path(unquote(url.path)).resolve() == pathlib.Path(sys.argv[3]).resolve()
assert 'include-system-site-packages = false' in (prefix/'pyvenv.cfg').read_text()
print(json.dumps({'status':'passed','profile':sys.argv[2],'installed_count':len(actual),'versions':actual,'direct_source_metadata':direct,'python':sys.version,'prefix':str(prefix),'scope':'Metadata only; no scientific modules imported'},sort_keys=True))
'''

ORIGIN_CHECK = r'''
from importlib.machinery import PathFinder
import json,pathlib,sys
root = pathlib.Path(sys.argv[1]).resolve()
result = {}
for name in sys.argv[2:]:
    locations = [str(root)]
    spec = None
    for part in name.split('.'):
        spec = PathFinder.find_spec(part, locations)
        assert spec is not None, name
        locations = list(spec.submodule_search_locations or [])
    assert spec.origin and pathlib.Path(spec.origin).resolve().is_relative_to(root), name
    result[name] = spec.origin
print(json.dumps({'status':'passed','root':str(root),'origins':result,'scope':'PathFinder restricted to declared checkout; no module execution or fallback to editable/site-packages'},sort_keys=True))
'''


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def write(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def cancelled(signum, frame):
    global CANCELLED
    CANCELLED = True


def stop_group(process):
    actions = []
    for sig, grace in ((signal.SIGTERM, 2), (signal.SIGKILL, 5)):
        try:
            os.killpg(process.pid, sig)
            actions.append(sig.name)
        except ProcessLookupError:
            break
        try:
            process.wait(timeout=grace)
        except subprocess.TimeoutExpired:
            pass
    try:
        process.wait(timeout=1)
        reaped = True
    except subprocess.TimeoutExpired:
        reaped = False
    try:
        os.killpg(process.pid, 0)
        group_absent = False  # Includes unreaped orphan zombies; conservatively stop.
    except ProcessLookupError:
        group_absent = True
    return {"signals": actions, "leader_reaped": reaped, "process_group_absent": group_absent}


class Diagnostic:
    def __init__(self, root, environment, pins, uv, python, started):
        self.root, self.environment, self.pins = root, environment, pins
        self.uv, self.python = str(uv), str(python)
        self.deadline = started + pins["limits"]["total_seconds"]
        self.records = []

    def save(self, status, **extra):
        write(self.root / "status.json", {"status": status, "updated_utc": now(),
              "pid": os.getpid(), "scope": self.pins["scope"], "steps": self.records, **extra})

    def check_time(self):
        if CANCELLED or time.monotonic() >= self.deadline:
            raise RuntimeError("Cancelled or overall 3600-second budget exhausted; no retry")

    def run(self, name, argv, *, cwd=None, stdin=None, environment=None):
        self.check_time()
        directory = self.root / "steps" / f"{len(self.records):02d}-{name}"
        directory.mkdir()
        env = dict(self.environment if environment is None else environment)
        argv = ["/usr/bin/nice", "-n", "19", *map(str, argv)]
        record = {"name": name, "argv": argv, "cwd": str(cwd or self.root),
                  "environment": env, "started_utc": now(), "status": "running",
                  "stdout": str(directory / "stdout.log"), "stderr": str(directory / "stderr.log")}
        if stdin is not None:
            (directory / "stdin.txt").write_text(stdin)
            record["stdin_sha256"] = digest(directory / "stdin.txt")
        self.records.append(record)
        self.save("running")
        started = time.monotonic()
        cap = min(self.pins["limits"]["stage_seconds"], self.deadline - started)
        reason = None
        process = None
        with (directory / "stdout.log").open("x") as stdout, (directory / "stderr.log").open("x") as stderr:
            input_file = (directory / "stdin.txt").open() if stdin is not None else None
            try:
                process = subprocess.Popen(argv, cwd=cwd or self.root, env=env,
                    stdin=input_file or subprocess.DEVNULL, stdout=stdout, stderr=stderr,
                    start_new_session=True)
                record["child_pid"] = process.pid
                self.save("running")
                while process.poll() is None:
                    if CANCELLED or time.monotonic() - started >= cap:
                        reason = "cancelled" if CANCELLED else "timeout"
                        break
                    time.sleep(0.2)
                if reason is None and (CANCELLED or time.monotonic() - started >= cap):
                    reason = "cancelled" if CANCELLED else "timeout"
                if reason:
                    record["cleanup"] = stop_group(process)
                else:
                    try:
                        os.killpg(process.pid, 0)
                    except ProcessLookupError:
                        pass
                    else:
                        reason = "remaining_process_group"
                        record["cleanup"] = stop_group(process)
                record["returncode"] = process.returncode
            except BaseException as exc:
                reason = "launch_or_monitor_error"
                record["error"] = type(exc).__name__ + ": " + str(exc)
                if process is not None:
                    record["cleanup"] = stop_group(process)
                record["returncode"] = process.returncode if process is not None else None
            finally:
                if input_file:
                    input_file.close()
        record.update(ended_utc=now(), wall_seconds=time.monotonic() - started,
                      stage_cap_seconds=cap, status=reason or ("passed" if record["returncode"] == 0 else "failed"))
        record["stdout_sha256"] = digest(directory / "stdout.log")
        record["stderr_sha256"] = digest(directory / "stderr.log")
        write(directory / "record.json", record)
        self.save("running")
        if record["status"] != "passed":
            raise RuntimeError(f"Step {name}: {record['status']}; first failure preserved, no retry")
        self.check_time()
        return (directory / "stdout.log").read_text()

    def source(self, name, spec):
        target = self.root / "sources" / name
        self.run(name + "-clone", ["/usr/bin/git", "clone", "--filter=blob:none", "--no-checkout", spec["repository"], target])
        self.run(name + "-sparse", ["/usr/bin/git", "sparse-checkout", "set", "--no-cone", "--stdin"],
                 cwd=target, stdin="".join("/" + path + "\n" for path in spec["sparse_files"]))
        self.run(name + "-checkout", ["/usr/bin/git", "checkout", "--detach", spec["commit"]], cwd=target)
        head = self.run(name + "-head", ["/usr/bin/git", "rev-parse", "HEAD"], cwd=target).strip()
        if head != spec["commit"]:
            raise ValueError("Wrong checkout commit: " + name)
        self.verify_source(target, spec)
        return target

    def verify_source(self, target, spec, phase="checkout"):
        self.check_time()
        for relative, expected in spec["expected_files_sha256"].items():
            if digest(target / relative) != expected:
                raise ValueError("Source hash mismatch: " + relative)
        for relative, info in spec["additional_files"].items():
            if digest(target / relative) != info["sha256"] or (target / relative).stat().st_size != info["bytes"]:
                raise ValueError("Build/package-data hash mismatch: " + relative)
        write(self.root / (target.name + "-source-validation-" + phase + ".json"), {
            "commit": spec["commit"], "scientific_source_hashes": len(spec["expected_files_sha256"]),
            "additional_files": spec["additional_files"], "sparse_files": spec["sparse_files"], "passed": True})

    def execute(self, profiles):
        version = self.run("uv-version", [self.uv, "--version"]).strip()
        if version != self.pins["uv"]["version"]:
            raise ValueError("Wrong uv version")
        self.run("python-version", [self.python, "-I", "-B", "-c",
            "import sys; assert sys.version_info[:3] == (3,10,12); print(sys.version)"])
        sources = {name: self.source(name, spec) for name, spec in self.pins["sources"].items()}
        for name in ("research", "forward-cpu", "crem", "mol2mol"):
            profile = profiles["profiles"][name]
            prefix = self.root / "venvs" / name
            python = prefix / "bin/python"
            self.run(name + "-venv", [self.uv, "--no-config", "venv", "--python", self.python, prefix])
            if set(profile["requirements"]) != set(INSTALL_ORDER[name]):
                raise ValueError("Unexpected requirement group membership")
            for filename in INSTALL_ORDER[name]:
                channel = profile["requirements"][filename]
                self.run(name + "-" + filename.removesuffix(".txt"), [self.uv, "--no-config", "pip", "install",
                    "--python", python, "--no-deps", "--index-url", channel["index_url"],
                    "-r", self.root / "bundle" / filename])
            install = profile["source_installation"]
            argv = [self.uv, "--no-config", "pip", "install", "--python", python, "--no-deps", "--index-url", "https://pypi.org/simple"]
            if install["editable"]:
                argv.append("--editable")
            argv.append(sources[install["source_id"]])
            self.run(name + "-source-install", argv)
            self.run(name + "-metadata", [python, "-I", "-B", "-c", METADATA_CHECK,
                self.root / "bundle/profiles.json", name, sources[install["source_id"]]])
            self.run(name + "-dependency-check", [self.uv, "--no-config", "pip", "check", "--python", python])
        probes = {
            "chembounce_frozen": ("research", ["chembounce_synthesis.identity", "fingerprint_db"]),
            "chembounce_current": ("research", ["chembounce_synthesis.identity", "fingerprint_db", "benchmarks.research.runner", "benchmarks.research.controls"]),
            "crem_adapter": ("research", ["benchmarks.external.crem_baseline"]),
            "mol2mol_adapter": ("research", ["benchmarks.external.mol2mol.adapter"]),
        }
        for name, (profile, modules) in probes.items():
            env = dict(self.environment, PYTHONPATH=str(sources[name]))
            self.run(name + "-source-origin", [self.root / "venvs" / profile / "bin/python",
                "-I", "-B", "-c", ORIGIN_CHECK, sources[name], *modules], cwd=sources[name], environment=env)
        for name, spec in self.pins["sources"].items():
            self.verify_source(sources[name], spec, phase="after-installation")
        self.save("passed_metadata_installation_diagnostic", cold_installation_performed=True,
                  functional_scientific_checks_performed=False, model_or_gpu_calls_performed=False,
                  numerical_or_independent_os_equivalence_claim=False)


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets-root", type=Path, required=True)
    parser.add_argument("--uv", type=Path, required=True)
    parser.add_argument("--python", type=Path, required=True)
    args = parser.parse_args()
    pins_file = HERE / "pins.json"
    if digest(pins_file) != PINS_SHA256:
        raise ValueError("Changed diagnostic pins")
    pins = json.loads(pins_file.read_text())
    bundle = REPO / "benchmarks/environments/research-validation-20261016"
    if {p.name for p in bundle.iterdir() if p.is_file()} != set(pins["bundle_files_sha256"]):
        raise ValueError("Bundle file membership changed")
    for name, expected in pins["bundle_files_sha256"].items():
        if digest(bundle / name) != expected:
            raise ValueError("Bundle changed: " + name)
    uv, python, root = args.uv.resolve(), args.python.resolve(), args.assets_root.resolve()
    if digest(uv) != pins["uv"]["sha256"] or digest(python) != pins["python"]["sha256"]:
        raise ValueError("Installer/base interpreter hash mismatch")
    if root.exists() or args.assets_root.is_symlink() or root.is_relative_to(REPO):
        raise ValueError("Assets must be a fresh directory outside the diagnostic checkout; no overwrite/resume")
    profiles = json.loads((bundle / "profiles.json").read_text())
    if any(pins["sources"][name]["commit"] != source["commit"] for name, source in profiles["sources"].items()):
        raise ValueError("Source commits differ from bundle")
    root.mkdir(parents=True, exist_ok=False)
    for name in ("sources", "venvs", "steps", "tmp", "uv-cache", "config", "data", "python-cache"):
        (root / name).mkdir()
    shutil.copytree(bundle, root / "bundle")
    shutil.copy2(pins_file, root / "driver-pins.json")
    shutil.copy2(Path(__file__), root / "driver.py")
    environment = {
        "PATH": "/usr/bin:/bin", "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8",
        "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
        "PYTHONNOUSERSITE": "1", "PYTHONDONTWRITEBYTECODE": "1", "PYTHONHASHSEED": "0",
        "TMPDIR": str(root / "tmp"), "UV_CACHE_DIR": str(root / "uv-cache"),
        "UV_PYTHON_INSTALL_DIR": str(root / "python-cache"), "UV_PYTHON_DOWNLOADS": "never",
        "UV_HTTP_RETRIES": "0", "UV_HTTP_TIMEOUT": "120", "UV_NO_PROGRESS": "1",
        "UV_CONCURRENT_DOWNLOADS": "1", "UV_CONCURRENT_BUILDS": "1",
        "XDG_CACHE_HOME": str(root / "uv-cache"), "XDG_CONFIG_HOME": str(root / "config"),
        "XDG_DATA_HOME": str(root / "data"), "GIT_TERMINAL_PROMPT": "0",
        "GIT_CONFIG_GLOBAL": "/dev/null", "GIT_CONFIG_NOSYSTEM": "1", "GIT_LFS_SKIP_SMUDGE": "1",
        "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1", "HF_HUB_DISABLE_TELEMETRY": "1",
        "TOKENIZERS_PARALLELISM": "false", "USE_TF": "0",
    }
    diagnostic = Diagnostic(root, environment, pins, uv, python, started)
    write(root / "launch.json", {"started_utc": now(), "argv": sys.argv, "cwd": str(Path.cwd()),
          "environment_whitelist": environment, "uv_sha256": digest(uv), "python_sha256": digest(python),
          "driver_sha256": digest(Path(__file__)), "pins_sha256": digest(pins_file),
          "limits": pins["limits"], "scope": pins["scope"]})
    signal.signal(signal.SIGINT, cancelled)
    signal.signal(signal.SIGTERM, cancelled)
    try:
        diagnostic.execute(profiles)
    except Exception as exc:
        diagnostic.save("failed_no_retry", error=type(exc).__name__ + ": " + str(exc))
        raise


if __name__ == "__main__":
    main()
