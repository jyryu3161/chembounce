# Original 630-attempt figures and recorded tables

`export_figures.py` presents an already completed, independently hash-pinned
[statistical replay](README.md). It does not run the archived estimator, generation,
models, chemistry, conformers, or common 3D scoring. Its implementation tests use
synthetic fixtures only; they are not a completed study analysis or parent review.

The scope is **Original frozen generation study**, source commit
`d972124d587e647ef48455cf1a038e6be96e5b04` (39 scientific files). It is not current
develop's 41-file pipeline or the supplementary 270-attempt comparison. The three
splits remain separate: scaffold 360 attempts / 6 targets, time 180 / 6, and target
90 / 3. Repeated targets across splits are not independent cohorts. Original method
identifiers are retained and sorted alphabetically, never by observed results.

## Completed-input gate

Supply the directory containing the replay's `validation.json` and three
`{scaffold,time,target}-comparison.json` files, plus an independently retained
SHA256 for `validation.json`. The exporter checks the passed 630-attempt contract,
the fixed historical helper hash, absolute tolerance `1e-12`, recorded scopes,
input provenance joins, all three comparison hashes, exact split/method grids,
recorded denominators, numeric types, finite values, and interval structure.
Missing or extra files, malformed JSON, duplicate keys, symlink inputs or ancestors,
incomplete grids, and stale hashes are rejected before plotting or output creation.
The output directory must be new and outside the input directory.

This gate verifies the recorded replay contract and externally supplied trust
anchor. It does not independently reassess the chemical results or prove that an
externally supplied hash came from an independent reviewer. All input bytes and
directory membership are checked again before the output manifest is written.

## Figures

Each split produces one PDF and PNG with these fixed panels:

| Recorded metric | Interpretation |
| --- | --- |
| `attempt_with_hop` | Attempt-level observed scaffold-hop output indicator |
| `attempt_with_verified_route` | Attempt-level output with a verified rule route |
| `common_3d_output_slot_coverage` | Available common-3D output slots divided by the cap of 3 |

Points are recorded target-macro means. Solid intervals are the recorded
descriptive, unadjusted target-bootstrap 95% CIs; the archived estimator used
10,000 draws and seed 20261016. The exporter does not draw bootstrap samples or
recalculate estimates. Orange dotted intervals and diamonds show infrastructure
missing-outcome **identification bounds, not confidence intervals**. Attempt,
infrastructure-interruption, and contributing independent-target counts are shown.
Missing means remain visibly missing.

Only 6, 6, and 3 independent targets support the respective splits. There are no
significance stars, rankings, tests of CI overlap, or superiority claims. This is
post hoc descriptive presentation, not preregistration. Common-3D slot coverage is
an output availability measure, not a conditional quality score. Verified rule
routes are bounded computational replay evidence, not experimental synthesis,
activity, yield, or probability of success. The panels do not include the separate
common-forward model re-inference audit.

## Recorded CSVs

All recorded metrics are exported, including metrics not plotted. Sorting is
deterministic. `None` becomes an empty field; integers and floats use their Python
string representation without decimal rounding. No pooling, imputation, or new
statistical calculation is performed.

- `metrics.csv`: recorded per-query, per-target, and target-macro values, CIs and
  contributing target counts, with primary or secondary population labels.
- `paired_contrasts.csv`: original pipeline-minus-comparator target-macro
  contrasts for all metrics, plus the recorded secondary complete-case bounded
  contrasts and their available query/target summaries. Each row names its scope.
- `infrastructure_bounds.csv`: recorded method and paired-contrast lower/upper
  bounds, separately labelled as identification bounds rather than CIs.
- `availability.csv`: recorded attempt, status, infrastructure, runtime-observation,
  query, target, and secondary complete-case availability counts.

The original paired contrasts pair **target means**, not matched candidate
populations. In particular, success-conditional common similarity means can involve
different successful attempts and candidates for each method; their contrasts
cannot establish like-for-like improvements. Secondary complete-case bounded
contrasts exclude infrastructure-missing pairs and can be selection-biased.
Timeouts remain failures. An attempt denominator is not a count of observed
conditional scores. Missing runtime observations are not zero runtimes.
Unbounded metrics have no finite identification upper bound supplied by this
presentation; only the recorded bounded metrics have infrastructure bounds.

## Invocation and artifacts

Use the existing plotting environment; the exporter adds no dependency or runtime
installation step:

```bash
python -B benchmarks/diagnostics/generation_replay/export_figures.py \
  --replay /path/to/completed/replay \
  --expected-validation-sha256 "$INDEPENDENTLY_RETAINED_VALIDATION_SHA256" \
  --output /path/to/new/original630-presentation
```

Matplotlib uses the noninteractive Agg backend. `manifest.json` is written last
and records input hashes before/after, exporter and replay-contract source hashes,
the historical helper hash, actual Python/Matplotlib/NumPy versions, all ten
figure/table hashes, table row counts, method order, and scope. An interrupted or
failed output directory is retained without a completion manifest and is never
overwritten or resumed by this command.

`--synthetic-fixture` visibly watermarks every figure and labels the manifest; it
does not change values, waive validation, or certify that a fixture is synthetic.
Use it for test artifacts, not study outputs. Focused tests are:

```bash
python -B -m pytest -q tests/test_original630_exporter.py
```

The synthetic 630-attempt fixture includes a target with unavailable conditional
3D scores, preserved failures and an infrastructure interruption. Tests cover
external trust anchors, full-grid and integrity rejection, null and numeric CSV
fidelity, labels, deterministic tables, input immutability, and existing-output
preservation. No real-study estimates or unfinished method comparisons are
produced by these tests.
