# Research validation environment specifications

These files record the **four observed environments** used by the supplementary comparison and forward evaluation. They are a portable version specification, not a hashed lockfile, container, or independently tested cold rebuild. They do not replace the original frozen generation campaign's other runtime. The repository's general `requirements-lock.txt` pins RDKit 2026.9.1 and does **not** describe these environments.

[profiles.json](profiles.json) contains complete installed package maps, source commits/file hashes, interpreter/config hashes, and links by SHA-256 to retained preparation evidence. [validation.json](validation.json) records metadata and source checks. [installation-origins.json](installation-origins.json) preserves public download URLs and archive hashes **only where already present in installation evidence**. Those hashes are evidence; the requirements files do not enforce a complete wheel lock. Local evidence records, models, DBs, datasets and private paths are not bundled.

| Profile | Installed distributions | RDKit / NumPy | Torch / Triton | Direct source installation |
|---|---:|---|---|---|
| `research` | 64 | 2025.9.6 / 1.26.4 | 2.8.0 / 3.4.0 | editable chembounce 2.0.0.dev0 |
| `forward-cpu` | 49 | 2025.9.6 / 1.26.4 | 2.8.0+cpu / absent | editable chembounce 2.0.0.dev0 |
| `crem` | 25 | 2025.9.6 / 2.2.6 | absent / absent | editable crem 0.3.3.dev1+gb84003cd8 |
| `mol2mol` | 108 | 2023.9.5 / 1.26.4 | 2.3.1+cu121 / 2.3.1 | non-editable local-source reinvent 4.4.22 |

The observed host was Linux x86_64, kernel 5.15.0-139-generic, glibc 2.31, Python **3.10.12** (GCC 11.2.0). All four virtual environments excluded system site-packages. `research` uses the PyPI distribution version `torch==2.8.0`; its retained runtime reports the CUDA 12.8 build. Do not replace that pin with a different CUDA channel. The profile includes the actual NVIDIA dependencies and setuptools 84.0.0, which was installed during editable project construction after the original 62 dependency pins were recorded. Python distributions alone do not reproduce the OS, driver, compiler, CPU/GPU hardware or temporary PEP 517 build environment.

CReM's isolated `importlib.metadata` view contains 25 installed distributions. Importing `crem` also exposes setuptools' vendored metadata: 41 raw entries, 37 unique names, including 12 additional names and four duplicated names. These vendored copies are **not separate installation requirements**. The prepared execution plan's 37-name map agrees after name normalization and duplicate handling; the 25-name installed map is the rebuild target. An import-based metadata count alone can therefore be misleading.

**Rebuild in new directories only.** The following is an unexecuted recipe. Choose an available Python 3.10.12 executable and a separately installed `uv`; the packaging assessment observed uv 0.11.7, without establishing which installer version produced every historical environment. Use `--no-deps` because the files enumerate the observed installed closure. Build isolation may resolve temporary tools; those tools and wheel/platform compatibility still require a fresh-build assessment. Do not reuse or modify any running study environment.

```bash
PYTHON310=python3.10
"$PYTHON310" -c 'import sys; assert sys.version_info[:3] == (3, 10, 12)'
BUNDLE="$PWD/benchmarks/environments/research-validation-20261016"
STUDY_ENVS="$PWD/rebuilt-study-environments"
STUDY_SOURCES="$PWD/rebuilt-study-sources"
mkdir -p "$STUDY_SOURCES"
for profile in research forward-cpu crem mol2mol; do
  uv venv --python "$PYTHON310" "$STUDY_ENVS/$profile"
done

uv pip install --python "$STUDY_ENVS/research/bin/python" --no-deps --index-url https://pypi.org/simple -r "$BUNDLE/research-pypi.txt"
uv pip install --python "$STUDY_ENVS/forward-cpu/bin/python" --no-deps --index-url https://pypi.org/simple -r "$BUNDLE/forward-cpu-pypi.txt"
uv pip install --python "$STUDY_ENVS/forward-cpu/bin/python" --no-deps --index-url https://download.pytorch.org/whl/cpu -r "$BUNDLE/forward-cpu-pytorch-cpu.txt"
uv pip install --python "$STUDY_ENVS/crem/bin/python" --no-deps --index-url https://pypi.org/simple -r "$BUNDLE/crem-pypi.txt"
uv pip install --python "$STUDY_ENVS/mol2mol/bin/python" --no-deps --index-url https://pypi.org/simple -r "$BUNDLE/mol2mol-pypi.txt"
uv pip install --python "$STUDY_ENVS/mol2mol/bin/python" --no-deps --index-url https://download.pytorch.org/whl/cu121 -r "$BUNDLE/mol2mol-pytorch-cu121.txt"
uv pip install --python "$STUDY_ENVS/mol2mol/bin/python" --no-deps --index-url https://pypi.nvidia.com -r "$BUNDLE/mol2mol-nvidia.txt"
uv pip install --python "$STUDY_ENVS/mol2mol/bin/python" --no-deps --index-url https://pypi.anaconda.org/OpenEye/simple -r "$BUNDLE/mol2mol-openeye.txt"

git clone https://github.com/jyryu3161/chembounce "$STUDY_SOURCES/chembounce-frozen"
git -C "$STUDY_SOURCES/chembounce-frozen" checkout --detach d972124d587e647ef48455cf1a038e6be96e5b04
for profile in research forward-cpu; do
  uv pip install --python "$STUDY_ENVS/$profile/bin/python" --no-deps --editable "$STUDY_SOURCES/chembounce-frozen"
done
git clone https://github.com/DrrDom/crem "$STUDY_SOURCES/crem"
git -C "$STUDY_SOURCES/crem" checkout --detach b84003cd8c84bc2afe98b5473570b7328e703ead
uv pip install --python "$STUDY_ENVS/crem/bin/python" --no-deps --editable "$STUDY_SOURCES/crem"
git clone https://github.com/MolecularAI/REINVENT4 "$STUDY_SOURCES/reinvent"
git -C "$STUDY_SOURCES/reinvent" checkout --detach eb61443a1a93a8cbd06510548b25cc77cfd9d3a6
uv pip install --python "$STUDY_ENVS/mol2mol/bin/python" --no-deps "$STUDY_SOURCES/reinvent"
```

Mol2Mol's original installation report records torch/torchvision/triton from PyTorch, twelve NVIDIA distributions from `pypi.nvidia.com`, and OpenEye-toolkits 2025.2.3 from the OpenEye index. Its pinned upstream [`pyproject.toml`](https://github.com/MolecularAI/REINVENT4/blob/eb61443a1a93a8cbd06510548b25cc77cfd9d3a6/pyproject.toml#L50) explicitly declares that dependency and index. A plain PyPI requirement is insufficient evidence for reproducing this installation. Installing the dependency does not establish a licensed OpenEye scoring capability; this comparison does not use licensed OpenEye scoring. This record neither supplies a license nor guarantees future package access. Do not silently omit that distribution while claiming the complete 108-package profile.

**Installation source and execution source differ.** The two editable chembounce distributions point at the original frozen checkout; `direct_url.json` encodes a local directory, not a commit. The frozen scientific file provenance establishes d972124d, while later documentation commits are immaterial to those scientific hashes. The supplementary runner explicitly uses a separate checkout at **2328588d0969369e694d9293a1dccc89392a8203** through its working directory/PYTHONPATH. CReM's adapter is source at **107c906ef98bd8345a7ddfcf0eeff7f4e1821835**, and Mol2Mol's adapter is source at **be6204e43060d6d168e7151177c0694f78a389c2**; neither adapter is an installed distribution. Prepare separate checkouts at these commits and compare the relative source hashes in `profiles.json` before building new plans. The historical forward audit uses the frozen source; supplementary evaluation specifies the current source separately. Do not substitute the latest develop source or rely solely on whichever editable import resolves first.

**Verify before any scientific execution.** This metadata-only check imports neither Torch nor a model. Run it for each new environment, then use `uv pip check --python "$STUDY_ENVS/$profile/bin/python"`. Research has no installed pip distribution, so `python -m pip` is not a universal verification command. Extra bootstrap packages must be explained or removed before claiming exact inventory parity.

```bash
profile=research  # repeat with forward-cpu, crem, mol2mol
"$STUDY_ENVS/$profile/bin/python" -I -B - "$BUNDLE/profiles.json" "$profile" <<'PY'
import importlib.metadata as metadata, json, re, sys
normalize = lambda name: re.sub(r"[-_.]+", "-", name).lower()
expected = json.load(open(sys.argv[1]))["profiles"][sys.argv[2]]["installed_packages"]
actual = {normalize(d.metadata["Name"]): d.version for d in metadata.distributions()}
assert sys.version_info[:3] == (3, 10, 12)
assert actual == expected, {k: (expected.get(k), actual.get(k)) for k in expected.keys() | actual.keys() if expected.get(k) != actual.get(k)}
print("Installed version inventory matches", len(actual))
PY
```

During preparation, all four maps matched the retained execution plan; 351 active dependency constraints and Python requirements passed metadata checks. The source inventories and recorded interpreter/config hashes also matched. No package was installed, downloaded, or exercised on a GPU, and no model/DB inference or cold rebuild was run for this bundle. New environments still require module-origin, source/asset pin, functional and resource checks before use. Matching versions alone does not establish numerical equivalence, reproducibility of timings, activity prediction, synthesis success, or journal readiness.
