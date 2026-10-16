# Original 630-attempt statistical replay

This tool packages and recomputes **stored statistics**, without running generation,
models, chemical identity calculations, conformer generation, or 3D scoring. It is
restricted to the original scaffold/time/target study: 360 + 180 + 90 attempts,
12 + 12 + 6 queries, seeds 42/43/44, and the original method groups and caps.
It excludes the supplementary 270-attempt study and the separate common forward
model re-inference audit. No real study results or bundle are supplied here.

The [figure and table exporter](EXPORT.md) accepts only a completed, externally
hash-pinned replay. It presents recorded statistics for the original frozen
39-file study, with separate infrastructure bounds and explicit conditional-score
limitations; its tests use synthetic data only.

`analyze_generation_historical.py` is a byte-identical archive of the original
`analyze-generation.py`, SHA256
`cb6992cda28895fa42791ad3969df5cbafd48a71b4ab8362ecc180cfeb4f2cc0`.
Only its pure `analyze_runs(runs, scored, targets, caps)` function is invoked;
its original filesystem-dependent command is not run. The local Ruff formatting
exception preserves this hash-bound file; lint remains enabled. New code is formatted.

The estimator is unchanged: average seeds within queries, queries within targets,
then give targets equal weight; 10,000 target bootstrap draws use seed 20261016.
Timeouts and other method failures retain their attempted denominators. An
infrastructure interruption retains its missing-outcome interpretation and bounds;
lost wall time stays null. Common similarity means remain success-conditional.
The original `paired_pipeline_minus_comparator` contrasts pair **target means**:
method-specific success-conditional populations can differ. These contrasts do not
establish like-for-like conditional-score improvements. Full-grid hop, route, and
output-slot coverage metrics retain the same attempted grid. Repeated targets
across splits are not new independent cohorts. Recovery/similarity are proxies,
not experimental activity, synthetic feasibility, or yield evidence.

## Completed-input contract

Both packaging and replay refuse incomplete grids, unfinished common scoring, or
missing comparisons. All three splits must pass before any statistics are computed.
`verify` checks integrity, completion, targets, caps, provenance, and score joins;
`replay` additionally recomputes every report field and compares it to the archived
comparison, with exact JSON types/structure and absolute float tolerance `1e-12`
(no relative tolerance). Integer/boolean differences, NaN, infinity, duplicate JSON
keys, unsafe artifact paths, symlinks, and unexpected bundle files are rejected.

The 44 required input files are copied without changing their bytes:

| Role | Original relative path |
| --- | --- |
| Three panels | `panels/direct-human-{scaffold,time,target}/panel.json` |
| Three primary protocols | `main-protocol-{baselines,gpu0,gpu1}.json` |
| Six secondary protocols | `secondary-{time,target}-protocol-{baselines,gpu0,gpu1}.json` |
| Nine run lists and nine manifests | `comparison-{baselines,gpu0,gpu1}/{runs,manifest}.json` and `comparison-{time,target}-{baselines,gpu0,gpu1}/{runs,manifest}.json` |
| Three common-score triplets | `common-generated-scores-8c/{protocol,progress,generated_scores}.json` and `common-generated-{time,target}-scores-8c/{protocol,progress,generated_scores}.json` |
| Three original comparisons | `analysis-generation-{scaffold,time,target}/comparison.json` |
| Frozen implementation provenance | `frozen-campaign-code-provenance.json` |
| ChEMBL source attribution | `chembl-source-provenance.json` (copied from the archive beside this tool) |

For the primary study, the method groups contain respectively 3, 3, and 4 methods;
the secondary groups contain 3, 1, and 1. Every method must have its exact
query/seed grid. Query targets come from the panel, each recorded run manifest must
bind the raw protocol and panel, and all nine manifests must match the original
39-source hash inventory and retain package versions/Python/platform metadata.
The run manifest's `dataset.source_artifacts` contains source ID → original
absolute filename strings, **not hashes**. Their IDs and lexical path construction
are checked against the panel's `sources[].artifact`; the separate source SHA/URL
fields are bound to the archived ChEMBL attribution. Original path strings are
never opened or resolved on the replay host. The original panels use absolute
POSIX paths; a relative source artifact without its historical panel directory
cannot be reconstructed by this compact schema and is rejected.
The original 360-second/24,000-MB/three-output outer budget is checked, including
method-specific cap overrides. Loader-added protocol fields are projected out only
when absent in the raw protocol; explicit original fields remain checked.

Common-score progress must bind its actual protocol and all three run files.
Score joins use `enumerate(evaluation.candidates, 1)`—the **valid-output ordinal**,
which can differ from a candidate's raw-response `rank` after invalid outputs are
excluded. The tool checks each target, seed, identity, cache key, failure status,
and complete occurrence count. It does not re-establish the underlying chemistry.
For three successful seeds, their stored metrics must reproduce the median and
all stability fields (`n`, mean, median, min, max, range, sample SD) within `1e-12`.
Any failed seed prohibits a 3D summary while preserving available 2D scores and
all trial records. This check is stored-number arithmetic, not 3D rescoring.

The byte-identical frozen provenance record binds commit
`d972124d587e647ef48455cf1a038e6be96e5b04` and 39 source hashes. Its timestamped
workspace checks are historical observations. Original run manifests supply package
versions; this bundle does not recreate package binaries, model files, or execution.

## Parent review and trust anchor

`pack` requires a separately prepared completion projection and its externally
retained SHA256. The parent evaluator creates this only after auditing the actual
completed study. Required fields are:

```json
{
  "schema_version": 1,
  "status": "passed",
  "total_runs": 630,
  "split_counts": {"scaffold": 360, "time": 180, "target": 90},
  "historical_helper_sha256": "cb6992cda28895fa42791ad3969df5cbafd48a71b4ab8362ecc180cfeb4f2cc0",
  "original_source_commit": "d972124d587e647ef48455cf1a038e6be96e5b04",
  "source39_sha256": {"<each original source path>": "<its SHA256>"},
  "input_sha256": {"<each of the 44 paths above>": "<its SHA256>"},
  "assessment_evidence": [{"artifact": "final-audit.json", "sha256": "<its SHA256>"}]
}
```

At least one actual audit JSON is required. Evidence paths must be safe relative
paths under the review file's directory. Evidence is copied under
`parent-review-evidence/`; the projection is copied as `parent-review.json`.
Additional descriptive projection fields are permitted. The audit's own text must
state what was checked: a hash or a `passed` string cannot itself prove independence,
scientific validity, or execution. A partial 576-attempt audit cannot satisfy this
contract. The compact bundle omits individual request/evaluation/response files;
those checks belong to the separately retained parent assessment, not this replay.

The attribution input uses the bundled file's SHA256
`5d680d5dad8a273989d3548035d1e3352453e4939bd6fa22f041c7ad0da5f94b`.
The frozen provenance input must match
`8c8d9a883fe29acc8ffc1d6c3826a0056def13a93fa852508617069f3d48075e`.

## Commands

Python 3.10+ and NumPy are required; the original analysis environment used
NumPy 1.26.4. Record the replay's actual version (written to `validation.json`).
No project installation, RDKit, GPU, network, or model assets are needed.

After independent completion review, with all paths supplied by the reviewer:

```sh
python -B benchmarks/diagnostics/generation_replay/replay.py pack \
  --source /path/to/completed/original-study \
  --parent-review /path/to/review/completion.json \
  --parent-review-sha256 "$REVIEW_SHA256" \
  --output /path/to/new-bundle

python -B benchmarks/diagnostics/generation_replay/replay.py verify \
  --bundle /relocated/new-bundle \
  --expected-manifest-sha256 "$PUBLISHED_MANIFEST_SHA256" \
  --output /path/to/new-verification.json

python -B benchmarks/diagnostics/generation_replay/replay.py replay \
  --bundle /relocated/new-bundle \
  --expected-manifest-sha256 "$PUBLISHED_MANIFEST_SHA256" \
  --output /path/to/new-replay-directory
```

Retain the manifest SHA printed by `pack` **outside** the bundle, for example in a
reviewed release record. Supply that independently obtained value when verifying;
computing a digest from an untrusted bundle and immediately trusting it is not an
external integrity check. Bundle paths are relative; recorded original absolute
paths inside protocols/manifests are preserved as historical strings and are never
followed. The publisher must review their disclosure before release. No credentials
or new private study records are bundled with this tool.

All outputs must be new; no overwrite/resume/partial mode exists. Verification and
replay outputs must be outside the input bundle. Input hashes are checked again
after copying or analysis. A pack writes its manifest last, and a replay writes its
passing validation last, so interrupted partial output directories cannot represent
a completed result. Retain a failed directory for diagnosis and use a new destination.

The provenance archive preserves [ChEMBL](https://www.ebi.ac.uk/chembl/) 37 source
URLs, bounded-sample scope, and upstream CC BY-SA 3.0 attribution. It does not contain
raw ChEMBL records; the panels contain the selected data. The bounded source sample
is not the complete target inventory. This data attribution does not assign a new
license to code, models, or other artifacts; consult their existing terms.

## Synthetic regression checks

```sh
python -B -m pytest -q tests/test_generation_statistical_replay.py
```

The fixture constructs 630 explicitly synthetic attempts, including timeouts,
infrastructure loss, failed 3D seeds, reused pairs, and noncontiguous raw ranks.
It checks byte-preserving packing, relocation after deleting the fixture source,
recomputation, input tampering, incomplete grids, numeric/types/fields, path/symlink
boundaries, output preservation, and before/after hashes. It does not contain or
evaluate real molecules, real study estimates, or model outputs. These tests run
under the repository's normal `tests/` CI discovery.
