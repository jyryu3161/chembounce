# Same-host study environment reinstall diagnostic

`rebuild.py` prepares four fresh environments from the [version specification](../../environments/research-validation-20261016/README.md). **Running it downloads source/packages and installs them.** The first installation diagnostic stopped at its 900-second command limit on 2026-10-16 at 04:36:30 UTC. The amended fresh installation subsequently completed; its installation review and separate fixed CPU comparison are recorded in [execution-result.json](execution-result.json).

The parent review confirmed all 54 installation/verification steps exited 0, the installed package counts were 64/49/25/108, and 658 source/build files, 13 original bundle files and 180 logs/records/source-validation files matched. Each fresh profile's selected CPU cases matched its own original environment: exact identity/fingerprint/reaction/arithmetic/API fields, 232 conformer scalars with maximum absolute difference 0 (declared tolerance `1e-10`), and 32 module-origin checks. The 232 scalars comprise one energy and 57 coordinates (19 atoms) per profile; profiles with different RDKit versions were not compared against one another. CPU probes explicitly selected the 41-file study source at `2328588d0969369e694d9293a1dccc89392a8203`; the frozen 39-file source remains the editable installation source. A separate installed REINVENT check verified 273 Python files plus one SAScore data file against pinned upstream/original bytes, RECORD hashes and new-environment origins.

These are same-host, selected-case results. No checkpoint or chemical DB was loaded, and no learned inference or GPU execution was performed. They do not establish an independent OS rebuild, fully hashed build, general numerical or whole-pipeline equivalence. Detailed local parent reports are hash-bound in the compact result and are not bundled.

[budget-amendment.json](budget-amendment.json) preserves the original failure's timestamps, resource/cleanup outcome and retained evidence hashes. The first PyPI stage resolved 63 packages in 3.82 seconds, downloaded Torch and cuDNN, and was downloading cuBLAS when the driver terminated it after 900.0648 seconds; no resolver conflict was reported. Its exit code was -15, its leader was reaped and its process group was absent. That failure alone established budget exhaustion, not an installation incompatibility. The amendment changes only the command/total budgets to 3,600/14,400 seconds and makes the timeout message use the configured budget. Package/source pins, channels, order, concurrency and cache policy remain unchanged. The original failed output is retained unchanged; the amended execution used a separate fresh directory and cache. This is a new environment engineering diagnostic, not a rerun or alteration of scientific campaign attempts.

The driver uses Python 3.10.12 and uv 0.11.7, with exact executable hashes in [pins.json](pins.json). It verifies all 13 original bundle files, including their manifest, before creating a fresh output directory. The eight requirement groups retain their original PyPI, PyTorch CPU/cu121, NVIDIA and OpenEye channels and install with `--no-deps`. Source packages retain the specified editable/non-editable modes. There is no dependency substitution or automatic retry.

The original installed metadata targets are 64 (`research`), 49 (`forward-cpu`), 25 (`crem`) and 108 (`mol2mol`) distributions. CReM's import-augmented 37-name view includes vendored setuptools metadata and is not a separate installation target. The driver imports neither CReM nor scientific/model packages when inspecting metadata. OpenEye is an installed dependency only; no licensed scoring is invoked and package access or licensing is not guaranteed.

## Source and package-data scope

All repositories are cloned from the URLs in `pins.json`, using `--filter=blob:none --no-checkout`, a recorded exact sparse-file list, detached fixed commits and `GIT_LFS_SKIP_SMUDGE=1`. This differs from the bundle's general full-checkout recipe. Git history remains available for CReM's setuptools-scm version calculation. Scientific file hashes and additional build-file hashes are checked before and after installation.

| Checkout | Fixed commit | Sparse files | Previously prepared source hashes |
|---|---|---:|---:|
| Frozen chembounce installation | `d972124d587e647ef48455cf1a038e6be96e5b04` | 85 | 39 |
| Current study execution | `2328588d0969369e694d9293a1dccc89392a8203` | 105 | 41 |
| CReM adapter | `107c906ef98bd8345a7ddfcf0eeff7f4e1821835` | 107 | 11 |
| Mol2Mol adapter | `be6204e43060d6d168e7151177c0694f78a389c2` | 115 | 10 |
| Official CReM | `b84003cd8c84bc2afe98b5473570b7328e703ead` | 80 | 80 |
| Official REINVENT4 | `eb61443a1a93a8cbd06510548b25cc77cfd9d3a6` | 461 | 455 |

The chembounce checkouts include all tracked Python files plus their build metadata and README. CReM includes all 80 tracked files, including the declared scripts, README and license. REINVENT4 includes its 455 prepared Python/metadata files plus six additional files needed to retain the package tree and build inputs. These include README, repository metadata, package CSV/TOML files and the declared `*.pkl.gz` package data: `reinvent_plugins/components/SAScore/fpscores.pkl.gz` (3,848,394 bytes; SHA-256 `10dcef9340c873e7b987924461b0af5365eb8dd96be607203debe8ddf80c1e73`). The scorer is not executed. Model priors/checkpoints, training molecules, chemical databases and study inference assets are excluded. Sparse membership and additional file sizes/hashes are explicit in the pins; the metadata checks alone do not establish complete installed-file parity.

Research and forward environments install the frozen chembounce source. The current execution source and both adapters are separate checkouts, not replacement editable installations. The final source-origin checks use `PathFinder` restricted to each declared checkout, without importing packages or falling back to the frozen editable checkout. They establish source availability; the parent evaluator must explicitly select the intended working directory/import path for later functional checks.

## Reviewed execution only

After independent review, invoke the driver from this checkout with the hash-matching executables and a **nonexistent directory outside this checkout**:

```bash
PYTHON310=/absolute/path/to/pinned/python
UV_EXECUTABLE=/absolute/path/to/pinned/uv
REBUILD_ASSETS=/absolute/path/to/new-study-rebuild-assets
"$PYTHON310" -I -B benchmarks/diagnostics/study_environment_rebuild/rebuild.py \
  --python "$PYTHON310" --uv "$UV_EXECUTABLE" --assets-root "$REBUILD_ASSETS"
```

The driver creates six source checkouts, four virtual environments, a private uv cache and temporary/configuration directories beneath that fresh root. Child processes receive an explicit environment whitelist, `CUDA_VISIBLE_DEVICES=''`, CPU thread limits of one, `nice -n 19`, disabled model-network access and no uv download retries. Existing source checkouts, environments, models, databases and campaigns are not reused or changed. No GPU/device discovery, model loading, inference or scientific scoring is performed.

Each command has a 3,600-second limit; the diagnostic has a 14,400-second total budget. The first failed command, metadata/source mismatch, cancellation or timeout ends the diagnostic. Each child receives its own process group. Timeout/cancellation and exceptions after process creation trigger bounded TERM/KILL/reap cleanup; cleanup may take up to eight additional seconds. An unreaped process or remaining group is recorded as a failure and blocks subsequent commands. OS-level uninterruptible processes cannot be guaranteed to terminate. There is no resume, overwrite, hidden retry or background queue.

`launch.json` binds executable/driver/pin hashes, arguments and the environment. Each `steps/` directory records argv, cwd, environment, timestamps, exit code, log hashes, and cleanup details where needed. `status.json` is an atomic progress/final snapshot. Source checks preserve both pre-installation and post-installation reports. The copied bundle, source pins, driver and complete stdout/stderr logs remain with the result, including on failure.

After installation, isolated stdlib metadata checks require the exact installed maps, direct source URL/editable state and environment-local distribution metadata. `uv pip check` verifies the declared installed dependency closure. Temporary isolated PEP 517 build tools may still be resolved from public indexes; they are not a hashed build lock. Matching installed versions and source hashes alone does not prove byte-identical wheels, functional or numerical equivalence, timing reproducibility, an independent OS rebuild, or journal readiness. Functional CPU checks belong to a separately recorded parent evaluation.
