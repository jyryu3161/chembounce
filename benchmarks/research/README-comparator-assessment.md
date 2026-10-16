# Assessment of the separate 270-attempt comparison

`assess_comparator_followup.py` consumes completed output from the separately
declared scaffold/time/target comparison. It provides no generation, scheduling,
model loading, CReM database API, automatic retry, or backfill. The original 630
attempts cannot substitute for the exact supplementary method/query/seed grid.
Implementing this tool does not authorize its execution on study data.

The primary endpoint remains engineering delivery of at least one output meeting
the existing scaffold-hop definition. This is not evidence of biological activity
or a calibrated probability of successful synthesis. All splits remain separate.

## Three immutable stages

Invoke the file directly, with the exact Python recorded in the prepared execution
plan. For `score` and `analyze`, set `CUDA_VISIBLE_DEVICES=''`,
`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`.
Do not invoke it as a package module from another checkout: the evaluator must
import only the tested 41-file scientific source and its pinned runtime.

Every action requires `--prepared PATH --prepared-sha256 PARENT_RETAINED_SHA`
and a fresh `--output PATH`. Each stage prints its manifest hash for independent
review. A hash read from an unreviewed input manifest is not an independent trust
anchor. Existing directories cannot be overwritten or resumed.

| Action | Additional required inputs | Output and completion gate |
| --- | --- | --- |
| `freeze` | None | Snapshot of all files in nine completed campaigns, requiring exactly 270 unique planned slots, valid protocol/panel/source/runtime provenance, and matching response hashes and per-attempt evaluation records. |
| `score` | `--generation PATH --generation-sha256 PARENT_RETAINED_SHA` | Fresh chemistry evaluation, all valid candidate occurrences, exact-pair common-3D caches and coverage. Requires the reviewed generation snapshot. |
| `analyze` | Generation arguments, plus `--scores PATH --scores-sha256 PARENT_RETAINED_SHA` | Fresh chemistry evaluation again, exact occurrence/cache joins, all-attempt metrics, split-specific target-macro summaries and paired contrasts. Requires complete reviewed common scoring. |

`score` is sequential and CPU only, outside generation budgets. It may be slow
for hundreds of unique pairs. After the complete generation gate and chemistry
recomputation, it creates `context.json` with `complete=false` and immediately
persists each completed, validated pair. An interruption leaves diagnostic
artifacts without a final manifest; partial work cannot produce comparison
estimates. There is no automatic resume or retry. Only full coverage plus repeated
input/source checks permit the final complete manifest.

## Scientific and accounting rules

The adapter reuses `evaluate_candidates`, `similarity_scores`, `score_stability`,
and `cluster_interval` from the tested source. It recomputes identity, Murcko and
generic scaffolds, bit/count Morgan similarity, held-out structural proxies, and
catalog/rule route replay. Saved evaluation fields must agree with recomputation.
Route replay reads the pinned catalog/rules and uses RDKit; it does not invoke a
learned forward model. Declared forward support is reported as such, including
the selection bias from sharing the pipeline's predictor.

Every accepted valid candidate occurrence retains its original raw rank, even
when invalid raw rows create gaps. Shared exact query/product pairs reuse common
3D computation, while occurrence denominators remain intact. Cache identity binds
the exact pair, tested source, identity policy, RDKit version, conformer settings,
and all three seeds (42/43/44). Each seed outcome and failure reason is retained.
Only three-seed success permits median 3D scores and stability summaries.

The primary denominator contains every attempt, including timeouts, invalid
responses, no-candidate outcomes, and infrastructure interruption. Terminal
failures deliver zero accepted candidates. Missing conditional similarities stay
null, with coverage counts; no supplied route is unassessed, not synthesis failure.
Raw CReM products and Mol2Mol sequences are reconciled with selection metadata and
returned outputs. Failed/partial work retains observed lower bounds or unknown
counts. Native work units are method-specific and cannot be compared as universal
evaluation counts. Sampled process-tree resources retain cold-start scope; these
shared-host measurements are not a controlled speed benchmark.

Estimates use seed means per query, query means per target, and equal target
weighting. Paired contrasts first subtract at the identical query/seed slot;
conditional contrasts require both values observed. Target bootstrap intervals
use 10,000 draws and seed 20261016, with no interval for a single target. They are
descriptive and unadjusted. Infrastructure sensitivity separately reports bounded
[0,1] counterfactual endpoints and paired complete cases; endpoints are not
confidence intervals and complete cases can be selection biased.

## Secondary scaffold-retention annotation

The generator additionally pins the parent addendum
`comparator-scaffold-retention-parent-addendum-20261016.json`, SHA-256
`38d0505cd67eedb5f19cae68925794406f84e6d815e12887f5bedebb2a4878d2`, in input files
and embeds its definitions in the analysis plan. This was specified after one
fixed nonbenchmark diagnostic and before supplementary study outputs; it is not
preregistered. Original decisions, budgets, ranking, and the hop endpoint remain
unchanged.

The annotation records atom/bond-aware query-scaffold containment in the candidate
scaffold, generic-skeleton containment in the same direction, core-retaining
scaffold expansion, and hop outputs where the query scaffold is not contained.
Containment uses actual graphs without chirality, not SMILES substrings or MCS.
Scaffold atom counts, occurrence counts, unique qualifying families, assessable
coverage, and all-attempt indicators are reported. Acyclic or unavailable
scaffolds yield explicit null annotation values. Neither graph expansion nor
non-containment proves a synthetic transformation, activity retention, or a
universal definition of successful scaffold hopping.

The focused fixtures use actual CPU identity/2D evaluation with hand-set 3D
outcomes. They cover all 270 slots, invalid-rank gaps, tampering, conditional and
paired missingness, infrastructure bounds, interrupted cache persistence, and
hand-derived extension/fusion/heteroatom/ring-size/acyclic graph cases. They do not
constitute a live study, model preflight, or experimental validation.
