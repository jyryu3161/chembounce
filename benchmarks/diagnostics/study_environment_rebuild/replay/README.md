# Replay the selected fixed CPU cases

This directory makes the selected CPU cases underlying [execution-result.json](../execution-result.json) executable and checkable. It is not an archive of every local parent report. The four [expected snapshots](expected/) preserve the original scientific fields and module SHA-256 values, with absolute module paths removed. [provenance.json](provenance.json) binds each original report, each projected snapshot, the unchanged [criteria](criteria.json), source pins and helper bytes.

`probe_scientific.py` and `audit_reinvent_installed.py` are **byte-identical copies** of the retained parent helpers. The first runs 14 fixed chemical identity/fingerprint cases, invalid-input checks, one small reaction, a seeded 19-atom conformer, CPU arithmetic and selected API imports. Each profile is compared only to its own original environment; Mol2Mol's RDKit version differs from the other profiles. The second independently checks installed non-editable REINVENT membership, all 273 Python files plus SAScore package data, wheel RECORD hashes and installed API origins. It hashes SAScore data without unpickling or running it. Neither helper loads a checkpoint or chemical DB, constructs a learned model, generates study candidates, or calls GPU APIs. They use assertions: **do not run them with `-O`/`-OO`**.

The local Ruff configuration preserves these retained helper bytes: only those two files are excluded from formatting, and only `probe_scientific.py` has an `F401` exception for its original unused import. Other lint rules and the comparator/test formatting follow the repository configuration.

`compare.py` uses only the standard library and **validates an existing report**; it does not run the probe or prove when/how that report was produced. It verifies the pinned artifacts, exact typed scientific fields, the per-profile Python/RDKit/API values, exactly 58 finite conformer scalars (one energy plus 19 × 3 coordinates), absolute tolerance `1e-10` with zero relative tolerance, mandatory module membership, actual module bytes and their declared locations. All 41 scientific files in the explicit `chembounce_current` checkout must match the source pins. Project origins must resolve there, dependency origins inside the matching new venv, and CReM/REINVENT source-probe origins inside their matching new upstream. Old editable/site-package fallback and symlinks escaping those roots fail even when bytes match. The separate installed REINVENT helper validates import from the installed package without upstream-path fallback.

Missing/extra fields, duplicate JSON keys, changed exact types or values, NaN/infinity, changed artifacts/files and wrong origins fail closed. Comparator output uses exclusive creation: existing files and dangling symlinks are refused. A failed comparison is recorded in a new result file and exits nonzero. Inputs and expected snapshots are never edited; do not loosen criteria to conceal a difference. The original helpers likewise require unused output paths; use a fresh result directory as below.

## Re-execute after preparing the pinned environments

These commands **run the small CPU cases**; they do not install/download anything. First prepare and validate the four environments and six source checkouts with the [reinstall instructions](../README.md). The comparator expects that driver's directory layout: `venvs/{profile}`, `sources/chembounce_current`, `sources/crem_upstream` and `sources/mol2mol_upstream`. `REBUILD_ASSETS` must name the new installation, and `REPLAY_OUTPUT` must not exist. Source files/model assets are not included here. The CLI resolves module bytes in those local checkouts/environments, so a saved report alone is insufficient.

From the repository root, in Bash:

```bash
set -e
REPLAY="$PWD/benchmarks/diagnostics/study_environment_rebuild/replay"
REBUILD_ASSETS=/absolute/path/to/completed-fresh-study-rebuild
REPLAY_OUTPUT=/absolute/path/to/new-fixed-cpu-replay
mkdir "$REPLAY_OUTPUT"
export CUDA_VISIBLE_DEVICES='' OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONHASHSEED=0
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 TOKENIZERS_PARALLELISM=false USE_TF=0

for profile in research forward-cpu crem mol2mol; do
  upstream_args=()
  case "$profile" in
    crem|mol2mol) upstream_args=(--upstream "$REBUILD_ASSETS/sources/${profile}_upstream") ;;
  esac
  "$REBUILD_ASSETS/venvs/$profile/bin/python" -I -B "$REPLAY/probe_scientific.py" \
    --profile "$profile" --source "$REBUILD_ASSETS/sources/chembounce_current" \
    "${upstream_args[@]}" --output "$REPLAY_OUTPUT/$profile.json"
  "$REBUILD_ASSETS/venvs/$profile/bin/python" -I -B "$REPLAY/compare.py" \
    --profile "$profile" --assets-root "$REBUILD_ASSETS" \
    --report "$REPLAY_OUTPUT/$profile.json" --output "$REPLAY_OUTPUT/$profile-comparison.json"
done

"$REBUILD_ASSETS/venvs/mol2mol/bin/python" -I -B "$REPLAY/audit_reinvent_installed.py" \
  --upstream "$REBUILD_ASSETS/sources/mol2mol_upstream" --pins "$REPLAY/../pins.json" \
  --output "$REPLAY_OUTPUT/mol2mol-installed-package.json"
```

To validate an already saved probe report, run only its `compare.py` command with a new output filename. Successful report comparison is separate from evidence of executing the probe. Keep command logs, source/environment provenance and raw outputs when recording a new execution.

The portable artifacts were packaged without rerunning scientific probes. Fixture tests use synthetic files and stored values, not actual molecular calculations:

```bash
python3.10 -I -B -m unittest discover \
  -s benchmarks/diagnostics/study_environment_rebuild/replay -p test_replay.py -v
```

This remains a same-host, selected-case assessment. It does not reproduce model inference, whole-pipeline outputs, GPU behavior, an independent OS, every dependency's binary bytes or a fully hashed build. Matching these cases does not establish general numerical equivalence, timing, activity/synthesis accuracy or scientific superiority. Original environment specifications and earlier result records remain historical, unchanged evidence.
