# Portable stored-statistics replay: supplementary 270 attempts

This implementation replays CPU statistics from **completed, independently reviewed
stored records**. It contains no supplementary study results. Actual normalization,
common scoring and final analysis must finish and receive parent review before
packing. The original 630-attempt study remains separate.

Python 3.10+ and NumPy are sufficient; the original statistics used NumPy 1.26.4.
No chemistry package, model, database, catalog, route replay, generation, or 3D
scoring is imported or invoked. Historical absolute paths are inert provenance
strings, including the original interpreter and scientific-source paths.

The replay reproduces every field of `comparison.json`, including attempt/status
counts, all metrics, observed/missing coverage, seed→query→target hierarchy, paired
query/seed contrasts, infrastructure identification bounds and complete-case
sensitivities. It uses the original 10,000-draw target bootstrap with seed 20261016
and no interval for one target. Comparison uses exact JSON types and structure,
with **absolute** floating tolerance `1e-12` and no relative tolerance.

Raw-accounting, unique-identity, and scaffold-retention summaries are preserved
byte-for-byte and bound to the final parent assessment. They are not re-derived by
this tool. The stored per-attempt metrics are reconstructed from the reviewed
candidate/annotation/score records; this establishes statistical consistency, not
independent chemical validity. Raw accepted-response joins, normalization, graph
annotations and route flags remain the parent's responsibility.

## Frozen execution boundary

`assess_historical.py` is the byte-identical reviewed assessor from commit
`37a1f2b838da0dc0af1577d059ba4fb633d0feee`, SHA256:

```
b9e24075c3940121c0dea6f0a6323a826824f372ca00129fe576194494f4a295
```

`evaluation_historical.py` preserves the evaluator from the tested 41-file source
commit `2328588d0969369e694d9293a1dccc89392a8203`. Neither full file is executed.
An explicit top-level AST allowlist extracts:

- Assessor: collections/statistics imports, `INFRA`, `RETENTION_FIELDS`, `BOUNDED`,
  `key`, `occurrence_key`, `mean`, `attempt_metrics`, `hierarchical`, `comparisons`.
- Evaluator: math/NumPy imports, `finite_scores`, `score_stability`, `cluster_interval`.

Every check reconstructs `statistics_extracted.py` from whole original node bytes,
including the final newline and original order. `extraction.json` records source
file, line boundaries, and individual node hashes. Only the extracted bytes execute.
`preparation` is injected as an inert namespace of fixed split/method constants;
the generator is never imported. No filesystem-dependent helper function is in
the allowlist. The extracted source SHA256 is
`737c13a0ccf8014459e7eb524c6fcbe12a6f13e48f72e1e2ec54fb1238c3d293`.

`pins.json` binds the source archives, extraction and immutable prepared/generation
metadata. The known preparation manifest is
`db76ab6355526ce833c06b443ccf9bd9b1222cbb26e253441f2f23ee2d4d9004`;
the generation-snapshot manifest is
`7d01e374c6e4e2a983f0c322b3036afc9258e3e0c7d7da334bf0403ab8f4d2f2`.
Frozen metadata is copied, not interpreted as permission to execute commands.

Strict JSON, path, hashing and immutable-write primitives reuse the reviewed
sibling `generation_replay/replay.py`, SHA256
`9751209528300cdc38250d68c3a6a8748ed873362af1d1dabfa14ecbb2234113`.
Its historical replay functions are never called. The transport copy and existing
ChEMBL attribution are included in the bundle. Source and extracted files have
local formatting exceptions to retain their reviewed bytes.

## Required input contract

`pack` takes four explicit local inputs. A completed parent review and its external
SHA are checked **before reading result archives**. There is no automatic preparation
or normalization step and no production override of the fixed study/code pins.

| Argument | Required contents |
| --- | --- |
| `--inputs` | `normalized-occurrences.json`, `integrity-validation.json` |
| `--analysis` | Complete original final analysis directory and its `manifest.json`; every original file is copied. |
| `--common` | Complete original sequential common-3D directory and its manifest, including every pair cache and failure record. |
| `--parent-review` | Independently produced final completion projection JSON, with evidence files relative to its directory. |

A normalized occurrence is the assessor's original full occurrence object:

```
split, method_id, query_id, target_id, seed, rank, smiles, query_exact,
candidate, route_supplied, serialized_route_sha256, scaffold_retention
```

`candidate` must equal the stored final recomputed candidate object. `rank` is the
**raw accepted-response rank**, an integer 1..3, increasing and unique within an
attempt. Gaps such as `[1,3]` are valid and are retained. Do not replace these with
contiguous valid-output ordinals from the original-630 forward replay.

The tool verifies exact 270-slot grid membership **and order**, target/seed/method
identity, terminal status, output cap, full occurrence order, candidate object,
original common-occurrence/cache join, and stored scaffold-retention occurrence
join. Failed attempts cannot acquire accepted output. Null annotation values remain
unassessed. Route-supplied flags are booleans; absent routes have null serialized
route hashes. The parent must establish their correspondence to original responses.

All unique common caches must be present exactly once in archive membership. Each
payload retains exact query/candidate strings, all 41 source hashes, conformer
configuration, seeds 42/43/44, identity policy and RDKit version. Every seed outcome
and reason is retained; only all-seed success permits a median and stability
summary. Recorded 2D values are compared with stored evaluated candidates, never
recomputed. Original common and final recomputed-attempt bytes must be identical.

The parent normalization integrity record requires these fields (extra explanatory
fields are allowed):

```json
{
  "passed": true,
  "total_attempts": 270,
  "candidate_occurrences": "<integer N>",
  "raw_candidate_joins_verified": true,
  "route_flags_verified": true,
  "scaffold_annotations_verified": true,
  "source_sha256": {
    "analysis/<each original file, including manifest.json>": "<SHA256>",
    "common/<each original file, including manifest.json and cache/*.json>": "<SHA256>"
  },
  "output_sha256": {"normalized-occurrences.json": "<SHA256>"}
}
```

This record must describe checks actually performed; a boolean does not establish
scientific independence. Its `source_sha256` map must contain every analysis/common
file and no others, using the displayed bundle-relative prefixes.

## Final parent completion projection

The projection requires the following fields. Angle-bracket examples are a schema
illustration, not executable input or a claim that final results exist.

```json
{
  "schema_version": 1,
  "status": "passed",
  "complete": true,
  "stage": "final_supplementary_270_completion",
  "total_attempts": 270,
  "split_attempts": {"scaffold": 108, "time": 108, "target": 54},
  "analysis_helper_sha256": "b9e24075c3940121c0dea6f0a6323a826824f372ca00129fe576194494f4a295",
  "prepared_manifest_sha256": "db76ab6355526ce833c06b443ccf9bd9b1222cbb26e253441f2f23ee2d4d9004",
  "generation_manifest_sha256": "7d01e374c6e4e2a983f0c322b3036afc9258e3e0c7d7da334bf0403ab8f4d2f2",
  "generation_parent_gate_sha256": "b9d1bc72af859c40cc2b8df67b3a156f71b21ebb4c3a9309eb9dc7b122f4c776",
  "source_commit": "2328588d0969369e694d9293a1dccc89392a8203",
  "input_sha256": {
    "analysis/<every file>": "<SHA256>",
    "common/<every file>": "<SHA256>",
    "normalized-occurrences.json": "<SHA256>",
    "integrity-validation.json": "<SHA256>"
  },
  "assessment_evidence": [
    {"role": "generation_parent_gate", "artifact": "generation-gate.json", "sha256": "b9d1bc72af859c40cc2b8df67b3a156f71b21ebb4c3a9309eb9dc7b122f4c776"},
    {"role": "final_integrity_and_statistics", "artifact": "final-review.json", "sha256": "<future actual final assessment SHA256>"}
  ]
}
```

The final evidence must be a JSON object with `passed: true`. Additional distinct
roles can retain normalization, independent graph checks, route assessment, or
other actual evidence. Artifact names must be safe relative paths; duplicate names
or roles are rejected. All evidence is copied under `parent-review-evidence/`.
The external projection SHA, not a hash self-calculated from an untrusted download,
is the entry trust anchor. Review text and recorded checks remain necessary.

## Commands and portability

Run the reviewed repository CLI; the resulting data bundle can move to another
machine without the original paths, campaign directories, chemistry libraries or
model assets. No installation or download is performed by these commands.

```sh
python -B benchmarks/diagnostics/supplement_replay/replay.py pack \
  --inputs /path/to/reviewed-normalization --analysis /path/to/final-analysis \
  --common /path/to/completed-common --parent-review /path/to/review/completion.json \
  --parent-review-sha256 "$PARENT_REVIEW_SHA" --output /path/to/new-bundle

python -B benchmarks/diagnostics/supplement_replay/replay.py replay \
  --bundle /relocated/bundle --expected-manifest-sha256 "$PUBLISHED_MANIFEST_SHA" \
  --output /path/to/new-replay
```

`verify` uses the same bundle arguments and a new JSON file as `--output`; it checks
stored joins/provenance/per-attempt statistics without recomputing the full
comparison. `replay` writes the recomputed comparison and a validation record last.
`pack` also reproduces the entire comparison before publishing a bundle manifest.

Keep the printed manifest SHA externally in the parent-reviewed release record.
Extra/missing archive files, symlinks including parent-directory symlinks, unsafe
paths, duplicate JSON keys, nonfinite numbers, incomplete grids, and overwrite or
resume are rejected. Original input bytes are checked again after analysis/copying.
A partial directory without its final manifest/validation is not a passing result.
Historical paths and source attribution remain visible metadata; publisher review
can assess those records before release. No raw model weights or executable model
assets are bundled.

## Interpretation and tests

Engineering delivery and retrospective similarity/recovery remain proxies.
No-route means absence of supplied evidence, not synthesis failure. Core expansion
and query-core non-containment are descriptive graph annotations, not proof of
activity retention or experimental synthesis. Forward evidence shares the pipeline
selection predictor. Reused pairs/seeds and targets repeated across splits are not
independent replication. Intervals are descriptive, unadjusted and limited by few
targets; infrastructure bounds are not confidence intervals and complete cases can
be selection biased.

`tests/test_supplement_statistical_replay.py` uses synthetic records with invalid
raw-rank gaps, scoring failures, no outputs, infrastructure interruption and
conditional paired missingness. Hand arithmetic checks complement full-report
reproduction, tamper tests, relocation, and import guards. Only synthetic metadata
trust anchors are substituted in fixtures; frozen helper/extraction bytes remain
unchanged and the production CLI has no fixture override. Tests neither claim real
study outcomes nor perform chemistry/model/3D calculations.
