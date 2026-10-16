# Stored common-forward statistical replay

This tool replays statistics from independently reviewed **normalized records** for
the original 630 generation attempts (scaffold 360, time 180, target 90). It uses
Python and NumPy only. It does not load RDKit, Torch, model weights, reaction rules,
or a building-block catalog; perform inference; reinterpret raw predicted products;
or revalidate chemical routes. No actual study results or public bundle are included
with this implementation. The supplementary 270-attempt study is outside its scope.

The forward analysis differs from the earlier generation comparison: it forms
differences at **exact query/seed pairs before target aggregation**, including
conditional metrics only when both paired observations are eligible. It retains
all planned attempts for coverage, infrastructure identification bounds, and
separate selection-biased complete-case summaries. The original target bootstrap,
10,000 draws, seed 20261016, eligibility accounting, and dependence diagnostics are
unchanged. No pooling of targets across splits or independent-route/seed uncertainty
is introduced. Few targets and shared reactant/context inputs across targets limit
the descriptive intervals.

`no_route_supplied` remains an absence of supplied evidence, not a chemical synthesis
failure. The common forward predictor is also a selection-gate model, so enrichment
is not independent validation. Missing reaction conditions and empirical context
remain limitations. Agreement does not establish activity, experimental synthesis,
yield, synthesis probability, or superiority of the complete pipeline.

## Frozen code boundary

The original helper is archived **byte-identically** as
`analyze_forward_historical.py`:

```text
01b4c9468d23bb2870083c0b8f66886fafbf0262b429e08e5e102e683e20d83e
```

The unchanged analysis plan is `forward-analysis-plan.json`:

```text
385be03fee3d2df6cfc04975cb63b42717088eabe62077bdbd3c6db9ca5db529
```

`statistics_extracted.py` copies only the original top-level statistical imports,
constants, error/require helpers, `run_key`, `mean`, `bootstrap`, `summarize`,
`paired_rows`, and `aggregate_split`. Its SHA256 is:

```text
0bf96138f0a137a6048a700d78bc064f6b278766d8b09da87f3f54ca45739cd1
```

`extraction.json` records each original line boundary and node-byte SHA. Whole
selected-node bytes, including their terminal newline, are joined in original
source order by one additional newline. Every verification reconstructs these
bytes from the archived helper using an explicit AST allowlist. The full original
helper is **not imported or executed** by the tool: its `sys.path` changes,
`exact_smiles`, `validate_step`, `validate_cache`, `build_replay`,
`validate_occurrences`, `analyze`, and original command are outside the execution
boundary. The extraction manifest SHA is
`d9c53088f3f5cdf04fdbeecc4851a55a6517dc3088311740db09992178ca2e27`.

Filesystem/JSON/hash primitives reuse the unchanged sibling
[`generation_replay/replay.py`](../generation_replay/replay.py), pinned at
`9751209528300cdc38250d68c3a6a8748ed873362af1d1dabfa14ecbb2234113`.
Use this tool from its repository checkout. That dependency is copied into bundles
as `transport_helpers.py`; it introduces no chemistry imports. The original and
extracted statistical files have precise local Ruff formatting exceptions to retain
their hash-bound bytes. Lint and normal formatting still cover new implementation.

## Inputs and normalization contract

Seven input files are preserved byte-for-byte under these bundle names:

| Bundle name | Source supplied to `pack` |
| --- | --- |
| `attempts.json` | `--inputs/attempts.json` |
| `normalized-occurrences.json` | `--inputs/normalized-occurrences.json` |
| `integrity-validation.json` | `--inputs/integrity-validation.json` |
| `recorded-analysis.json` | `--analysis/analysis.json` |
| `recorded-validation.json` | `--analysis/validation.json` |
| `audit-protocol.json` | `--audit/protocol.json` |
| `audit-progress.json` | `--audit/progress.json` |

An attempt is exactly `{split, method_id, query_id, target_id, seed, status}`.
The plan supplies exact methods, query/target mappings, seeds, and output caps.
All 630 unique terminal slots must be present, including failures and infrastructure
interruptions. No automatic preparation or chemistry-dependent normalization occurs.

An occurrence preserves its original fields and adds `matched` and `evaluated`
booleans, `input_keys` (ordered step/context SHA256 strings, with repeats retained),
and `route_key` (SHA256 or null). Required original fields are `split`, `method_id`,
`query_id`, `target_id`, integer `generation_seed`/`rank`, stored `smiles`, `status`,
and `forward_all_steps_matched`. The replay checks:

- Occurrences belong to successful attempts, with the planned target and contiguous
  valid-output ordinals `1..n`, at most three per attempt. These are valid-output
  ordinals, not raw response ranks.
- `assessed` records have nonempty input keys, `route_cache_key == route_key`, and
  original matched/evaluated flags equal the normalized booleans. Matched implies
  evaluated. Reused route keys retain the same outcomes, candidate identity, and
  ordered input keys.
- `no_route_supplied` has false normalized flags, empty input keys, null route key,
  no cache key, and null/absent original evaluated flags. In particular, the actual
  no-route schema can omit `all_steps_evaluated`.
- `audit_failed` has false normalized flags, empty input keys, a supplied route key,
  null original outcome flags, and an explicit error.

The independent normalization assessment must bind both normalized-file hashes,
630 attempts, occurrence and unique-route counts. Raw run and audit-input hashes
must agree between that assessment, the original validation, and progress. The
assessment's original `forward-cpu-runtime/final-analysis-20261016/analysis.json`
and `validation.json` hashes must also match the two bundled recorded files, using
one unambiguous historical path-suffix match per file. Verification diagnostics
record both checked hashes. The original validation must report the fixed
helper/plan, frozen 39-source inventory
and commit, completed grids, and matching model provenance. No recorded path is
dereferenced: absolute paths remain historical strings.

The **parent integrity assessment**, not this tool, establishes the raw accepted
output identity/ordinal join, route chemistry, prediction interpretation, and
normalization completeness. This replay checks their stored normalized form and
hash bindings. It cannot replace that assessment. The original validation metadata
is preserved and compared exactly; it is not regenerated by running its scientific
checks. All statistical split sections are recomputed, and the entire reconstructed
analysis is compared to the original with exact JSON types/structure and absolute
float tolerance `1e-12`, without relative tolerance.

## Final parent gate

The parent independently reviews integrity **and statistics** before preparing a
completion projection with these required fields:

```json
{
  "schema_version": 1,
  "status": "passed",
  "total_attempts": 630,
  "expected_attempts": {"scaffold": 360, "time": 180, "target": 90},
  "analysis_helper_sha256": "01b4c9468d23bb2870083c0b8f66886fafbf0262b429e08e5e102e683e20d83e",
  "analysis_plan_sha256": "385be03fee3d2df6cfc04975cb63b42717088eabe62077bdbd3c6db9ca5db529",
  "extracted_statistics_sha256": "0bf96138f0a137a6048a700d78bc064f6b278766d8b09da87f3f54ca45739cd1",
  "original_source_commit": "d972124d587e647ef48455cf1a038e6be96e5b04",
  "normalization_source_validation_sha256": "<original validation.json SHA256>",
  "input_sha256": {"<each of the seven bundle input names>": "<SHA256>"},
  "assessment_evidence": [{"artifact": "final-parent-assessment.json", "sha256": "<SHA256>"}]
}
```

At least one actual final assessment JSON is required. Its path is relative to
the projection's directory and it is copied under `parent-review-evidence/`.
Additional descriptive projection fields are allowed. A hash or a `passed` string
does not automatically establish scientific independence: evidence text must state
which checks the parent actually performed. The integrity report alone does not
substitute for the separately required final review.

## Commands and trust anchor

Python 3.10+ and NumPy suffice; the original analysis used NumPy 1.26.4. There are no
installs/downloads in the commands. `pack` requires an externally retained parent
review SHA. `verify` checks the complete bundle and normalization/provenance gates;
`replay` additionally executes the pinned statistical extraction.

```sh
python -B benchmarks/diagnostics/forward_replay/replay.py pack \
  --inputs /path/to/parent-normalized-records \
  --analysis /path/to/original-final-analysis \
  --audit /path/to/completed-common-forward-audit \
  --parent-review /path/to/final-review/completion.json \
  --parent-review-sha256 "$REVIEW_SHA256" \
  --output /path/to/new-forward-bundle

python -B benchmarks/diagnostics/forward_replay/replay.py verify \
  --bundle /relocated/forward-bundle \
  --expected-manifest-sha256 "$PUBLISHED_MANIFEST_SHA256" \
  --output /path/to/new-verification.json

python -B benchmarks/diagnostics/forward_replay/replay.py replay \
  --bundle /relocated/forward-bundle \
  --expected-manifest-sha256 "$PUBLISHED_MANIFEST_SHA256" \
  --output /path/to/new-replay-directory
```

Retain the manifest SHA printed by `pack` **outside the bundle** in a reviewed
release record. Computing and trusting a hash from the same untrusted download
does not establish authenticity. All output destinations must be new; replay and
verification outputs must be outside the input bundle. There is no resume/overwrite
mode. Duplicate JSON keys, nonfinite numbers, unsafe artifact paths, symlinks,
extra/missing bundle files, changed inputs, and incomplete grids are rejected.
Input bytes are rehashed after copying/analysis. The bundle manifest and passing
replay validation are written last; incomplete output is not a passing result.

The bundle also preserves the original frozen-code provenance and existing ChEMBL
37 source attribution from the sibling generation replay. ChEMBL's recorded
CC BY-SA 3.0 attribution does not assign licenses to the tool or models. Raw model
weights, package binaries, catalog/database files, and per-route raw prediction
artifacts are not copied. Publisher review remains necessary for historical path
strings and any supplied parent assessment before public release.

## Synthetic tests

```sh
python -B -m pytest -q tests/test_forward_statistical_replay.py
```

Fixtures construct 630 synthetic attempts with absent routes, matched/mismatched
or incomplete routes, audit failures, timeouts, infrastructure interruptions, and
reused inputs across targets. They compare the extracted function against the
original pure function, test relocation and full-report equality, and check failure,
tamper, type, path, and overwrite boundaries. Import guards reject Torch/RDKit/model
imports in the pack/replay fixture. The full original helper is loaded only for
synthetic function-equivalence testing, with its temporary `sys.path` change restored.

To avoid publishing real protocol inputs as fake fixtures, tests replace **only
the fixture plan hash/path trust anchor** with a synthetic plan's hash. Statistical
code, original/extracted helper pins, completion rules, and validators are unchanged;
the production CLI exposes no such override. Separate tests reject substitutions
under the real default plan/helper pins. These are engineering tests, not real
forward evidence or a scientific result.
