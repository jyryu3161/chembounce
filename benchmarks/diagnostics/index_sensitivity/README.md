# Frozen-generator index sensitivity helper archive

This standalone diagnostic leaves the original 630-run campaign, its code, assets,
and outcomes unchanged. It makes **126 paired comparisons / 252 generator calls**:
90 diverse conditions (30 distinct public input strings × seeds 42/43/44) and 36
legacy conditions (12 main-panel inputs × the same seeds). These represent 270
planned pipeline attempts after deduplicating their identical generator inputs.
It does not call SynLlama, ReactionT5, synthesis reconstruction, or 3D scoring.

`sensitivity.py` is a **byte-identical archive of the actual execution helper**,
whose SHA-256 is
`db04fac3f57fd82b978442477bc8eef9627e62a6934461d8c7484040a9e726a5`.
The original was retained at the local path
`benchmark-results/research-validation/fingerprint-index-sensitivity-tool/sensitivity.py`.
This packaging does not replace that running file or modify its plan, processes,
assets, or original campaign outcomes.

At initial implementation, only fixture tests and read-only preflight checks had
been performed. **At packaging on 2026-10-16, the parent reviewer has started the separately
reviewed 252-call diagnostic; the comparison is incomplete.** No counterfactual
result or performance superiority is reported here. A fixture pass is engineering
evidence, not an activity or synthesis validation.

This directory contains only helper source, tests, and documentation. Local
indices, source libraries, panels/protocols, exhaustive audit reports, plans,
archived reference snapshots, model weights, and generated molecule records are
not bundled and must be prepared separately. No sensitive source data is added.
Packaging the code grants no third-party dataset/model redistribution rights or
new license to those assets; obtain them under their respective upstream terms.
The local-path examples below are preparation requirements, not downloadable
public artifacts or a claim of self-contained public replay.

## Preconditions

- Use the **original runtime and original frozen scientific checkout**, associated
  with implementation commit `d972124d587e647ef48455cf1a038e6be96e5b04`.
  Current develop contains later repairs and is not a substitute for either arm.
  The original interpreter was `.venv-synllama-min/bin/python`; its generator
  environment used Python 3.10.12, RDKit distribution 2025.9.6 (runtime 2025.09.6),
  NumPy 1.26.4, pandas 2.3.3, SciPy 1.13.1, ScaffoldGraph 1.1.2, ODDT 0.7,
  and psutil 7.2.2.
  Matching package names alone does not establish runtime equivalence. The helper
  checks all 39 source hashes from `frozen-campaign-code-provenance.json`, verifies
  the generator source-file inventory, and checks the runtime in every worker.
  It does not import the repaired develop generator into either arm.
- Supply a **separate final index path** and the actual exhaustive parent audit
  `validation.json`. The provisional current-runtime index is not acceptable.
  `validate_index_report` requires the existing audit fields `complete`,
  `input_asset_sha256`, `fingerprint_parameters`, `rows`, `fingerprints_checked`,
  `runtime_canonical_unique_strings`, `raw_strings_strictly_sorted_unique`,
  `database_rdkit_version`, and `runtime_rdkit_version`. Invalid rows, bit
  mismatches, canonical changes, and canonical duplicates must all be empty
  lists. The actual report bytes and hash are retained. A
  `passed_fingerprint_compatibility=true` flag alone does not suffice.
- The index metadata must match the scaffold source SHA and current RDKit,
  radius 2 / 2048 bits / no chirality / big bit order, and the finalized
  `rdkit-roundtrip-fixed-point-v1` policy with an 8-roundtrip bound. The policy
  marker is a builder declaration; the separate exhaustive report supplies the
  validation evidence. Full metadata, including exclusion counts, is retained.
- Existing protocol/panel paths are explicit inputs. The extractor checks the
  fixed seed, thread, output, wall, generator, and downstream design budgets.
  Only public query SMILES enter a worker. It does not read models or send
  assay labels or held-out structures into generation.
  Secondary protocols omit a direct thread declaration; their primary-protocol
  hash linkage is verified and this absence is recorded, rather than invented.
  Both diagnostic arms use the explicitly recorded primary CPU thread settings.

## Commands (review before execution)

Set `FROZEN_ROOT` to the separately prepared original checkout and `ARCHIVE_ROOT`
to the checkout containing this directory, then run from the frozen repository
root. Set `FINAL_INDEX` and `FINAL_INDEX_REPORT`
to the final independently validated files; the commands intentionally supply no
provisional default. `--workers 1` is the default and `2` is the only alternative,
to be selected before planning after a resource check.
Process niceness is also recorded and must match between planning, controller,
and workers. The examples use `nice -n 19` for both commands.
**Set `PYTHONHASHSEED=0` before Python starts in both commands.** Setting it inside
an already-running Python process does not fix that process's hash seed. The
helper propagates and verifies the environment value in its fresh workers.
The original references did not record a fixed hash seed; old-arm replay remains
an explicit check rather than an assumption of identical archived output.

```bash
FROZEN_ROOT=/path/to/original-frozen-checkout
ARCHIVE_ROOT=/path/to/checkout-containing-this-archive
cd "$FROZEN_ROOT"
B=benchmark-results/research-validation
HELPER="$ARCHIVE_ROOT/benchmarks/diagnostics/index_sensitivity/sensitivity.py"
PYTHONHASHSEED=0 nice -n 19 .venv-synllama-min/bin/python "$HELPER" plan \
  --source-root . \
  --frozen-manifest "$B/frozen-campaign-code-provenance.json" \
  --python .venv-synllama-min/bin/python \
  --new-index "$FINAL_INDEX" --index-validation "$FINAL_INDEX_REPORT" \
  --campaign "$B/panels/direct-human-scaffold/panel.json" "$B/main-protocol-gpu0.json" "$B/comparison-gpu0" \
  --campaign "$B/panels/direct-human-scaffold/panel.json" "$B/main-protocol-gpu1.json" "$B/comparison-gpu1" \
  --campaign "$B/panels/direct-human-time/panel.json" "$B/secondary-time-protocol-gpu0.json" "$B/comparison-time-gpu0" \
  --campaign "$B/panels/direct-human-time/panel.json" "$B/secondary-time-protocol-gpu1.json" "$B/comparison-time-gpu1" \
  --campaign "$B/panels/direct-human-target/panel.json" "$B/secondary-target-protocol-gpu0.json" "$B/comparison-target-gpu0" \
  --campaign "$B/panels/direct-human-target/panel.json" "$B/secondary-target-protocol-gpu1.json" "$B/comparison-target-gpu1" \
  --workers 1 --output "$B/reproduction-index-sensitivity-plan-hashseed0"
```

This first command creates a plan and an immutable, hash-checked reference
snapshot **at planning time**. It does not generate designs. After reviewing the
plan, its exact SHA, references, and index audit:

```bash
PYTHONHASHSEED=0 nice -n 19 .venv-synllama-min/bin/python "$HELPER" run \
  --plan "$B/reproduction-index-sensitivity-plan-hashseed0/plan.json" \
  --plan-sha256 "$REVIEWED_PLAN_SHA256" \
  --output "$B/reproduction-index-sensitivity-run-hashseed0"
```

Existing output directories are refused. There is deliberately no automatic
resume or overwrite. Outputs under a campaign, index, or plan snapshot are
refused. Input/code/index/snapshot hashes are rechecked before and after execution.
Plans bind the helper's local path and hash. A code archive is not a replacement
for an existing live plan; prepare a separate reviewed plan for any reproduction.

## Scope and output interpretation

Both arms use the frozen generator and identical explicit options. Only the
fingerprint-index path changes, besides the output directory. The indexed search
uses `research-v1`, scaffold top 8, fragment cap 10, candidate-request cap 512,
exploration fraction 0.25, and final cap 3. Shape calculation is false and
candidate thresholds are empty. Legacy retains its frozen implementation,
including any internal defaults; its three requested seeds are not collapsed.

The downstream projection uses exactly the frozen workflow expression:

```python
all_designs = list(dict.fromkeys(canonical(s) for s in designs))
selected = all_designs[:2]
```

It imports only the frozen `common.canonical`, not the synthesis workflow.
Canonicalization errors remain unassessed downstream comparisons, even if the
generator itself completed. An empty completed generator result is retained and
is distinguishable from a failed generator call.

- `plan.json` binds public inputs, full options, source/runtime/asset/report hashes,
  all original grid members, and reference availability. Complete generator
  references are snapshotted even when their original pipeline timed out; their
  original status stays recorded and **is never promoted to campaign success**.
  Missing, incomplete, and failed generator references remain explicit.
- Per-arm `request.json`, `worker.log`, original `designs/` artifacts,
  `worker-result.json` when available, and `report.json` preserve outcomes.
  `runs.json`, `summary.json`, and `progress.json` are atomic snapshots.
- Raw ordered output SMILES, canonical ordered deduplicated designs, and the
  actual first two are compared separately. Full table cells are lossless strings,
  preserving numeric text without rounding; raw file bytes are another comparison.
  Fragment/replacement/decision files and summaries are separate comparisons.
  Only the explicit index path is removed from structured summary comparisons.
  Provenance manifests retain their legitimate path/time/index differences.
- Every available complete reference must match the old-arm replay in table bytes,
  selected first two, intermediate file hashes, and normalized summary. A mismatch
  blocks an asset-only attribution for that condition. Missing references do not
  count as successful replay validation.
- The fixed denominator is 126 conditions and 252 arm calls, including timeout,
  worker error, malformed/missing artifact, interruption, and unattempted calls.
  Partial-run summaries do not relabel unattempted calls as chemical failures.

Each arm starts a fresh CPU-only subprocess with the same thread/environment
settings. A process-group timeout enforces 360 seconds; TERM then KILL also removes
descendants that ignore TERM. RSS is sampled across live descendants and checked
against 24,000 MB per arm. At concurrency 2 the aggregate process limit may reach
48,000 MB; this is why concurrency is predeclared. SIGTERM/keyboard interruption
requests cancellation and retains completed records. A SIGKILL cannot finalize
the current record; the last atomic progress/summary and remaining per-arm files
require parent review. No automatic restart is attempted.
Reaping is bounded (0.5 seconds after TERM, then 1.0 after KILL). An unreaped
process, such as one stuck in uninterruptible kernel I/O, is explicitly recorded;
dispatch stops and parent review is required. No finite user-space timeout can
guarantee termination of such a process.

Timing includes subprocess startup, imports, provenance, generation, and helper
finalization; original adapter timing and downstream remaining budget have different
scopes. Concurrent work, filesystem cache, and shared-host load confound timing.
**Even identical first-two designs do not establish whole-pipeline performance
equivalence under 360 seconds or corrected-method activity/synthesis superiority.**
The 994 historical canonical spelling changes can affect hash sampling and tie
ordering; the contrast is not restricted to the two mismatched fingerprint rows.

## Separate fragment-boundary evidence

The index-only contrast does not by itself test the repaired query canonicalization.
The parent independently audited the original fragment Mol objects and corrected
fixed-point handling under the exact campaign interpreter and thread environment
in `fingerprint-fragment-boundary-campaign-runtime-parent/validation.json`.
The earlier `fingerprint-fragment-boundary-parent/validation.json` and
`snapshot-completeness-check.json` remain historical records. That separate evidence is not treated as an
extra generation arm or silently included in this helper's tests/results.

Relevant frozen call boundaries for independent reconstruction are
`chembounce.py:563` (fragment construction), `:578` (`fragment_info.tsv` writes),
`:633` / `:657` (diverse/legacy fast-search calls), and `fingerprint_db.py:145`
(fingerprinting the original Mol and exclusion identity). The fragment TSV has
`Fragment_no` and `SMILES` columns in order. Reconstruct with `utils.init_mol` then
`scaffoldgraph_fragmenter.get_all_murcko_fragments(..., break_fused_rings=False,
iteration_round=1)`; check saved order. Re-parsing the saved strings alone does not
recover the original-Mol fingerprint boundary. Any unresolved fragment or query-bit
difference requires review before claiming repaired-generation equivalence.

## Tests

The 39 original fixture cases are preserved. `HERE.parents[2]` still resolves to
the repository root at this packaged location. The only substantive test
adaptation is the small real-generator smoke: it now requires the explicit
`CHEMBOUNCE_FROZEN_SOURCE_ROOT` environment variable, verifies the original commit
identity and source manifest, and uses that checkout read-only. It skips with a
preparation explanation when the variable is unset. A supplied path lacking the
manifest, or containing changed source, fails clearly. It never substitutes
current develop for the frozen generator. All fixtures load **this public copy**
of `sensitivity.py`.

From `ARCHIVE_ROOT`, using the original runtime:

```bash
cd "$ARCHIVE_ROOT"
PYTHONHASHSEED=0 PYTHONPATH="$PWD" \
  CHEMBOUNCE_FROZEN_SOURCE_ROOT="$FROZEN_ROOT" \
  "$FROZEN_ROOT/.venv-synllama-min/bin/python" -m pytest \
  benchmarks/diagnostics/index_sensitivity/test_sensitivity.py -q
```

To exercise the default scope without supplying the local frozen checkout:

```bash
env -u CHEMBOUNCE_FROZEN_SOURCE_ROOT PYTHONHASHSEED=0 PYTHONPATH="$PWD" \
  "$FROZEN_ROOT/.venv-synllama-min/bin/python" -m pytest \
  benchmarks/diagnostics/index_sensitivity/test_sensitivity.py -q -rs
```

The explicit-root run exercises all 39 cases. The default run exercises 38 and
skips only the real frozen-generator smoke. These commands require the original
runtime dependencies; the skip does not remove other package requirements.
Packaging validation of this public copy passed all 39 cases in 9.24 seconds
with the explicit read-only frozen checkout. The default run passed 38 cases
and skipped that one smoke in 6.62 seconds. Ruff check and format check passed.

Tests use synthetic campaigns, mock orchestration, and one small synthetic index
through the real frozen generator. They check the complete grid, exact selection
semantics, private-label exclusion, report/source/asset rejection paths, snapshot
immutability, process-group timeout, memory limit, bounded concurrency, failure
denominators, and output protection. The mocked 252-task orchestration test does
not execute 252 generators or use public campaign chemistry.
