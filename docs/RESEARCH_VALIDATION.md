# Research validation protocol

This protocol supports a computational methods study. Passing software tests does not establish improved scaffold hopping, preserved activity, or experimental synthesizability. Report actual comparisons, including negative results; the runner does not create results for unexecuted methods.

## Claims and corresponding evidence

| Question | Required evidence | Executable support |
| --- | --- | --- |
| Was the scaffold changed, rather than only protonation, tautomerism, isotope or stereo? | Exact state plus a separate standardized chemical family, Murcko family and generic ring skeleton. Preserve the input chemical state for scoring. | Shared `chembounce_synthesis.identity`; benchmark deduplication and scaffold metrics use the same versioned policy. |
| Can early product deduplication discard a valid design? | Failure-injection regression with more than one design for a product; surviving alternate must be attempted. | Production workflow regression tests; benchmark counts a failed query instead of silently removing it. |
| Does the method generalize? | Pinned external assay sources, chemical grouping, held-out reference scaffolds, development/test separation, and independently audited model/prompt inventories. | `protocol.py`, `import_chembl.py`; unknown or partial training inventories remain explicitly unknown/partial. |
| Is the synthetic route credible? | Catalog/rule replay, explicit missing atom/reagent/context information, separately reported external model/literature/experimental support. | Independent replay of returned graphs in `runner.py`; original external evidence records remain available. Replay alone is not experimental validation. |
| Are similarity results reproducible? | Common evaluation settings, all scoring failures, multiple seeds, score ranges, input-pose provenance, positive and negative controls. | `controls.py`; success-conditional ROC/AP and all-control missing-score bounds. Similarities are not calibrated probabilities. |
| Are more candidates also more different scaffolds? | Unique chemical families, standardized Murcko scaffolds and generic ring skeletons; active scaffold recovery as an explicitly labelled proxy. | `evaluation.py`; no activity inference from a scaffold match. |
| Does a 2D prefilter suppress known hops? | Recall over every held-out active scaffold before applying 3D scoring, with the tested threshold recorded. | `controls.py` reports full-population prefilter recall separately from capped 3D controls. |
| Is the comparison fair under finite compute? | Equal outer time/candidate caps, declared objective-call units, all timeouts/errors/overruns retained, process-tree resource records. | `runner.py` enforces wall/memory/output caps; objective units are explicit and not assumed equivalent. |

## External activity panel

`import_chembl.py` consumes retained ChEMBL API artifacts and a manifest with the exact request URLs and SHA-256 values. It does not silently re-fetch a newer database release. The API snapshot's version and bounded sampling rule should accompany a publication. Re-running a live API URL is not a guarantee of identical content.

```bash
python -m benchmarks.research.import_chembl \
  --manifest path/to/retained/chembl/manifest-with-assays.json \
  --output benchmark-results/panel-scaffold \
  --split scaffold --seed 42 --queries-per-target 2
```

The default panel first joins retained assay metadata and requires direct single-protein assignment (`confidence_score=9`, `relationship_type=D`), `assay_organism=Homo sapiens`, an exact target identifier match, binding assay type, and an explicitly null variant-sequence annotation. Human *target* annotation alone is insufficient. These are structured-metadata criteria: a null variant field does not prove a wild-type construct, and confidence 9 does not establish comparable assay conditions. Missing metadata fails closed; `--allow-unverified-assay-targets` produces an explicitly diagnostic panel without a direct-human-assay claim.

The default labels then use binding measurements in nM: active at or below 1,000 nM, inactive at or above 10,000 nM. These are declared benchmark definitions, not universal medicinal chemistry thresholds. Censored relations are respected: `<100 nM` supports the active definition, while `>100 nM` does not establish inactivity. Nonstandard, flagged, variant-target, conflicting-family, ambiguous and out-of-scope records are logged as exclusions. Assay identifiers, endpoint, relation, value, document publication year and supporting activity IDs are retained. Assay heterogeneity still matters; target-level concordance cannot replace a matched experimental assay.

Use separate `--split scaffold`, `--split target`, and `--split time --cutoff-year 2020` panels to test different claims. Family overlaps across development/test partitions are rejected. Scaffold and target partitions enforce their additional group constraints. Time partitions use the recorded **publication year**, not an invented assay date, and do not establish unseen pretrained-model data. Historical and post-cutoff measurements are curated separately: a future contradictory observation cannot retroactively remove or relabel a training compound. Families already observed in development are excluded only from the later test set.

A dataset consists of `schema_version`, `dataset_id`, `scope`, `sources`, `compounds`, `queries`, and `split`. Each source has `id`, `url`, `artifact`, and `sha256`. Compound records contain `id`, `smiles`, `target_id`, `label` (`active`, `inactive`, `decoy`, or `unlabelled`), `split` (`train`, `validation`, `test`, or `query`), and `source_id`; active/inactive labels also require `assay_id` and `label_basis`. A query points to its public `compound_id`; it may be a known development ligand. Other test-labelled actives on that target with different scaffold families form the private recovery set. A `query` partition permits a public input on a target whose reference compounds are held out. Synthetic fixtures must declare `scope: synthetic_fixture` and are never scientific evidence.

An optional query `reference_sdf: {"artifact": "...", "sha256": "..."}` supplies only that query's public bound pose. Its structure is validated by the scorer. Held-out reference ligand poses must not be passed to a generator as designs or input structures.


Retain an auditable offline snapshot with the companion source tool:

```bash
python -m benchmarks.research.sources package \
  --manifest path/to/retained/chembl/manifest-with-assays.json \
  --output benchmark-results/source-snapshot --archive benchmark-results/source-snapshot.zip
python -m benchmarks.research.sources verify benchmark-results/source-snapshot
```

The packager verifies raw bytes, preserves original manifests and declared license metadata, and emits portable rewritten manifests. `fetch` retrieves the exact declared URLs and requires their pinned hashes; a changed live API response is rejected rather than silently replacing the snapshot.

## Baselines and ablations

The built-in methods are:

- `nearest_catalog`: 2D Morgan retrieval from a declared library, with a seeded bounded subset when the score budget is smaller than the library.
- `nearest_3d`: fixed-protocol 3D retrieval from that library; embedding failures count as attempted evaluations.
- `reaction_enumeration`: one-step forward enumeration from allowed building blocks and rules, with a global tuple budget and round-robin reaction allocation. This is an explicit simple baseline, not an exhaustive synthesis planner.
- `synllama_only`: direct model projection of the public query, exact catalog/rule reconstruction, exclusion of the unchanged input family, and Morgan ranking. It does not receive a known held-out active as a design.
- `pipeline`: actual ChemBounce generation from the query, then the synthesis workflow. Saved-proposal configurations are reproducible diagnostic replay; they omit inference cost and must not be presented as a live-model wall-time comparison.
- `external`: an explicit command consuming `{request}` and producing `{response}`. This supports a real [RuSH](https://github.com/molML/RUSH) execution when its dependencies and ROCS license are available. A missing installation/license remains an unexecuted comparison; another method is not relabelled as RuSH.

Use end-product retrieval as well as building-block retrieval to avoid a molecular-complexity mismatch. Export a development-only library without assay labels or test/query families:

```bash
python -m benchmarks.research.training_library \
  --dataset benchmark-results/panel-scaffold/panel.json \
  --output benchmark-results/development-library.tsv

python -m benchmarks.research.build_protocol \
  --output benchmark-results/methods.json \
  --catalog path/to/building-blocks.tsv --rules path/to/rules.json \
  --retrieval-library benchmark-results/development-library.tsv \
  --methods nearest_catalog nearest_3d reaction_enumeration \
  --seeds 42 43 44 --wall-seconds 300 --max-candidates 20 --max-evaluations 1000000
```

The exported library contains experimental development structures, not supplier stock records. Exact held-out molecules are absent by construction. In a strict scaffold split, a train-only retrieval baseline also cannot return an unseen test scaffold; it is a useful floor and similarity/complexity control, not the sole novelty comparator. Report generated-method and published-method comparisons too.

For `pipeline`, add `--synthesis-config`, `--scaffold-db`, optionally `--fingerprint-db`, and an explicit generator-options JSON. For `synllama_only`, add a local `--model-directory`, `--prompt`, `--device`, and optionally a complete chemical `--prompt-inventory`. All model files and assets are hashed; generated protocol paths are relative to the protocol location. This includes provider/prompt/evidence-model assets indirectly referenced by the synthesis config and directory-based fingerprint databases. Changing nested model contents after protocol creation fails validation. Training inventories are explicitly unknown by default. A lack of overlap in a sampled or partial inventory does not prove a clean holdout. Unsupported normalization and invalid inventory entries are retained as unassessed records with original indices, input hashes and reasons; they downgrade audit coverage and cannot establish a clean holdout.

The [recorded training-inventory audit](../benchmarks/results/training-inventory-audit-20261016/README.md) found no complete inventory tied to either exact model checkpoint. A post hoc comparison with SynLlama's upstream train-labelled 1,000-molecule evaluation sample found zero exact matches and zero matches to its 928 assessable families; 72 sample families remained unassessed. The 14-source snapshot covers 19,387 appearances (5,618 exact structures, 5,616 families), including only the 137 main-campaign runs recorded at 2026-10-15 18:27:40–41 UTC. Independent numerical recomputation passed. Exact checkpoint membership of the sample remains unverified, later candidates are outside this snapshot, and these negative sample results do not establish clean pretrained-model holdout. An author-published ZINC inventory exists, but its contents and exact ReactionT5 base-checkpoint association were not verified here.

`benchmarks/research/configs/single_factor_ablations.json` specifies legacy search, no 3D evaluation, no refinement, and encounter-order allocation. Supply it through `--variants` with a pipeline parent. `disable_3d` calls the same production workflow with `scoring_mode="2d"`: route collection, property/identity constraints, independent evidence queues, candidate allocation and screening caps are shared. Set `generator_options.calculate_shape=false` in both the parent and no-3D variant (the builder default); the no-3D adapter rejects generator-side shape work. Only the complete synthesis 3D evaluation component is removed, including its dependent refinement; ranking uses 2D query/design similarity, unavailable geometric values are null, and no poses are exported. Pocket constraints or 3D score thresholds must be removed explicitly or the run rejects the incompatible configuration. Under binding wall budgets, different scorer costs can still change completed-work counts. Compare variants to a parent that explicitly enables the removed feature. The separate `configs/compound_baselines.json` contains `design_only_2d`, which removes synthesis, evidence, 3D scoring and refinement together and cannot support a single-factor conclusion.

## Execution, outputs and failure accounting

```bash
python -m benchmarks.research.runner \
  --dataset benchmark-results/panel-scaffold/panel.json \
  --protocol benchmark-results/methods.json \
  --output benchmark-results/comparison
```

Output directories must be new. Each method/query/seed has public request, raw response, stdout/stderr, evaluation and resource files. Ground-truth references and labels are excluded from generator requests. Methods are trusted local programs, not sandboxed adversaries; this contract does not prevent a deliberately malicious external method from reading shared filesystem data.

The external response schema is:

```json
{
  "candidates": [{"smiles": "...", "score_samples": [0.7, 0.8]}],
  "evaluations": 2,
  "evaluation_unit": "3D candidate screens",
  "metadata": {}
}
```

Candidates may include a full `route` graph. The evaluator replays it against pinned catalog/rule assets; a provider's claimed evidence label cannot promote an invalid route. External synthesis evidence is retained separately and does not become experimentally validated by being reported. Quantitative evidence level claims should distinguish rule replay, forward-model support, supplied observations and actual audited experiments.

The primary matched compute budget is wall time, including imports and model loading. A fingerprint comparison, reaction tuple, 3D screen and language-model token are not equivalent operations. Their count units and internal counters remain explicit; the pipeline also reports generator/refinement work. The common outer evaluation ceiling is a safety bound; method-specific `budget_overrides.max_evaluations` may be lower, but wall time and output count cannot be overridden. Builder defaults are a full-library-cap 2D retrieval, 256 3D/pipeline evaluations, and 4,096 reaction/SynLlama reconstruction evaluations beneath the declared common maximum; override these explicitly with `--method-evaluation-caps` when needed. This permits full inexpensive retrieval without forcing every 3D method to attempt the same arbitrary count.

Adapters receive the same monotonic deadline. Retrieval/enumeration stop with a serialization margin; 3D retrieval reserves time based on completed screen durations. Reaction enumeration lazily prepares each rule's pools on first visit and immediately evaluates one tuple before the next rule; pool scans check the deadline every 256 structures. Prepared/unvisited rules, incomplete pools and preparation costs remain explicit. SynLlama preserves its original sampling modes and seed offsets while checking sample and token boundaries. The pipeline retains fully screened/refined partial results when it stops at a work boundary. Incomplete refinement seed blocks never become successful refined results. Internal stopping reasons are recorded; exceeding the outer hard deadline is still a failed run. A single unpredictable long operation can cross that deadline, so inspect actual failure rates. Compare warm-model throughput separately if it is part of the research question. RSS and CPU are sampled for the process tree; short peaks and GPU memory are not claimed measured. Setup/evaluation costs outside generation are separate.

Timeouts, launch errors, malformed responses, memory limits and output/evaluation-budget violations remain in the attempted-query denominator. Partial files from failed runs never count as successful final candidates. Invalid structures, repeated chemical states and unchanged input families cannot inflate recovery. Summary statistics average seeds and queries within a target and then weight targets equally; confidence intervals resample targets. A single-target panel cannot identify between-target uncertainty.

SIGTERM and SIGINT stop the campaign without starting another attempt. When a worker is interrupted, its `status` is `interrupted`, accepted `response` is empty, and previously sampled resource values remain in its evaluation record and `runs.json`. Its raw response file and hash remain available for diagnosis, including a complete-looking response written before interruption. The runner signals only the worker's own process group: SIGTERM, at most one second for the leader to exit, then SIGKILL to the group even if the leader has already exited. It waits at most one further second to reap the leader; `cleanup` records the signals sent and whether the leader was reaped. No GPU-wide or job-wide process search is used. Repeated termination signals during this bounded cleanup and checkpointing retain the first signal. Uncatchable SIGKILL and a process stuck in the kernel cannot provide the same graceful-finalization guarantee.

Each attempt's evaluation, the accumulated runs and summary are replaced atomically as complete JSON files. `campaign.json` separately records `completed` or `interrupted`, the attempted-run count and the termination signal/phase. A validated worker result remains valid if termination arrives during subsequent evaluation or checkpointing; campaign termination is not a reason to relabel that finished result. Setup cancellation records zero attempts. These records describe attempted work only; unstarted entries are not fabricated as failures, and this is not an automatic resume mechanism. The CLI exits 143 for SIGTERM and 130 for SIGINT. The library APIs `run_benchmark` and `bounded_method` require the main thread and restore previous SIGTERM/SIGINT handlers on exit. `run_benchmark` raises `BenchmarkInterrupted` (a `BaseException` with `runs`, `termination` and `exit_code`) after checkpointing; `bounded_method` returns its report with termination metadata. Signal delivery during setup or evaluation is observed at the next Python work boundary; no asynchronous exception interrupts file replacement or process creation. Unrelated exceptions still propagate.

## Common scoring, negative controls and filter recall

```bash
python -m benchmarks.research.controls \
  --dataset benchmark-results/panel-scaffold/panel.json \
  --output benchmark-results/common-score-audit \
  --max-per-label 20 --seeds 42 43 44 \
  --count-tanimoto-threshold 0.5 \
  --run-results benchmark-results/comparison/runs.json
```

This command performs 2D evaluation by default. Add `--conformers fixed-conformer-settings.json` for common multi-seed 3D scoring of controls and generated candidates, independent of each method's ranking score. Selection of the capped control sample is deterministic; full-population prefilter recall is evaluated before that cap. Bit Morgan similarity matches the synthesis fingerprint convention; sparse-count Morgan similarity matches the generator filter.

Inactive experimental observations and decoys are separate. ROC/AP are reported only for active/inactive labels, with the selected population and failed-score counts. Success-conditional estimates are accompanied by pessimistic and optimistic all-control bounds, so failed scores cannot disappear from the interpretation. Similarity scores are not probabilities; probability calibration requires a separately fitted and held-out calibrated model, which this framework does not pretend to provide.

## Recorded controls: 2026-10-16

The [recorded control evaluation](../benchmarks/results/research-validation-20261016/README.md) contains 240 query–control pairs from 120 compound IDs (119 exact structures), 12 queries and six targets in the strict direct-human-assay ChEMBL scaffold panel. Full-population held-out active-scaffold recall at count Tanimoto ≥0.5 was **3.16%**, with a target-bootstrap 95% interval of **1.15–5.18%** (10,000 resamples). This describes loss under that filter, not the recall of a completed generation method.

On the same complete-case controls, target-macro ROC AUC was 0.639 for bit Morgan, 0.637 for count Morgan, 0.542 for shape, 0.564 for pharmacophore similarity and 0.550 for the joint 3D score. These results do not show improved active/inactive discrimination from 3D scoring. All 42 failed pairs remain recorded; each failed because of the specified-stereochemistry policy. Joint-score ROC bounds including missing scores were 0.386–0.644, which are missing-data bounds rather than confidence intervals. One CA2 query had no successful 3D controls and remains in failure accounting and these bounds.

See the [figure](../benchmarks/results/research-validation-20261016/controls.pdf), [independent analysis](../benchmarks/results/research-validation-20261016/control_analysis.json) and reproduction instructions for capped balanced sampling, complete-case selection and six-target uncertainty limits. A [bundled ChEMBL source snapshot](../benchmarks/results/research-validation-20261016/chembl-source-snapshot.zip) supports control reconstruction; model-training inventories remain unknown. At this record, the 360-run main generation comparison was in progress and 270 secondary runs were queued; no method-superiority or activity-probability claim follows from these controls.

A subsequent [assay-context audit](../benchmarks/results/research-validation-20261016/assay_heterogeneity.json) found that 1,103/5,567 records have a representative assay annotated as cell-based, and 1,229 have representative descriptions saying the protein origin is unknown. Counting every retained supporting assay raises the latter to 1,250 records. Structured human assignment and these descriptions therefore leave unresolved provenance ambiguities. This is not a uniform biochemical assay panel. Requiring a shared retained assay identifier between query and control leaves only six pairs, four with successful 3D scores, and only one target with both labels. That subset cannot establish general activity discrimination. These post hoc sensitivity checks preserve the frozen panel and expose data limitations; they are not threshold selection or evidence of improved performance. This motivated a separately selected same-assay evaluation; the frozen panel and original estimates remain retained.

The frozen [publication-year and held-out-target controls](../benchmarks/results/research-validation-20261016/secondary_controls.json) also completed. Time controls contain 240 pairs across six targets, with 21 stereo-policy failures; target controls contain 120 pairs across three targets, with 29 such failures. Paired bit/joint-3D ROC AUCs are 0.564/0.510 and 0.438/0.340, respectively. Count-Tanimoto ≥0.5 retains only 0.298% and 2.005% of eligible active scaffolds on target-macro averaging. These alternative partitions preserve the same assay limitations and do not establish a 3D discrimination benefit. Independent metric recomputation and reconstruction from the published source archive matched both records; see the result README for figures, uncertainty and reproduction commands.

## Completed original frozen generation study

The original study used commit `d972124d587e647ef48455cf1a038e6be96e5b04`
and its pinned 39-file scientific source inventory. It completed all 630 planned
method/query/seed attempts: scaffold 360, publication-year 180 and held-out-target
90. These are the original frozen outputs, including their failures. They are not
measurements of the later 41-file current pipeline or the separate supplementary
270-attempt comparator study. Earlier partial-status records remain historical
snapshots rather than being rewritten as completed experiments.

Metrics were averaged over seeds 42, 43 and 44 within each query, queries within
each target, and target clusters equally. The splits contain respectively 12, 12 and 6 queries
and 6, 6 and 3 target clusters. Reused targets and reaction inputs do not establish
statistical independence, and the splits must not be pooled as independent
cohorts. The fixed outer budget was 360 seconds including cold start, 24,000 MB
sampled process-tree RSS, and at most three returned candidates per attempt;
method-specific inner limits remain recorded in the protocols.

### All-method output and failure accounting

The complete grid contains 434 successful attempts, 140 with no candidates,
55 timeouts and one infrastructure interruption. The interruption is an unknown
outcome, not a demonstrated chemistry failure. Conservative realized-output
estimates retain its attempted slot; separately reported identification bounds
show its possible contribution to bounded metrics. These bounds are not confidence
intervals. Its lost runtime remains missing, and observed-runtime summaries do not
substitute zero. Timeouts retain zero accepted outputs. Complete-case analyses
excluding only infrastructure-missing outcomes are secondary and selection-biased.

The tables below retain every method and split in alphabetical method order.
Hop is the attempt-level scaffold-hop output indicator; route is the attempt-level
output with a verified rule route; 3D coverage is available successful common-3D
output slots divided by the output cap of three. Coverage measures available
output, not success-conditional score quality. Zero route coverage for retrieval
methods means no verified route was supplied, not that synthesis was attempted
and failed. Displayed values are target-macro means and descriptive, unadjusted
target-bootstrap 95% CIs (10,000 draws, seed
20261016). Few target clusters limit inference; CI overlap is not a statistical
test, and no significance or superiority claim is made.

**Scaffold split: 360 attempts, 6 target clusters.** [PNG](../benchmarks/results/research-validation-20261016/original-generation-study/presentation/scaffold.png) / [PDF](../benchmarks/results/research-validation-20261016/original-generation-study/presentation/scaffold.pdf).

| Method | Attempts | Hop, mean [95% CI] | Route, mean [95% CI] | 3D coverage, mean [95% CI] |
| --- | ---: | ---: | ---: | ---: |
| `nearest_3d` | 36 | 0.694 [0.611, 0.778] | 0.000 [0.000, 0.000] | 0.722 [0.611, 0.861] |
| `nearest_catalog` | 36 | 0.667 [0.500, 0.833] | 0.000 [0.000, 0.000] | 0.889 [0.722, 1.000] |
| `pipeline` | 36 | 0.167 [0.028, 0.333] | 0.167 [0.028, 0.333] | 0.139 [0.028, 0.250] |
| `pipeline_encounter_allocation` | 36 | 0.417 [0.278, 0.556] | 0.417 [0.278, 0.556] | 0.315 [0.213, 0.444] |
| `pipeline_legacy_search` | 36 | 0.361 [0.139, 0.639] | 0.361 [0.139, 0.639] | 0.343 [0.139, 0.583] |
| `pipeline_no_3d` | 36 | 0.806 [0.667, 0.944] | 0.806 [0.667, 0.944] | 0.546 [0.444, 0.648] |
| `pipeline_no_refinement` | 36 | 0.444 [0.306, 0.556] | 0.444 [0.306, 0.556] | 0.398 [0.278, 0.528] |
| `pipeline_rule_only` | 36 | 0.306 [0.139, 0.500] | 0.306 [0.139, 0.500] | 0.269 [0.111, 0.444] |
| `reaction_enumeration` | 36 | 0.833 [0.611, 1.000] | 0.833 [0.611, 1.000] | 0.620 [0.407, 0.833] |
| `synllama_only` | 36 | 0.722 [0.611, 0.833] | 0.833 [0.694, 0.944] | 0.667 [0.472, 0.861] |

**Publication-year split: 180 attempts, 6 target clusters.** [PNG](../benchmarks/results/research-validation-20261016/original-generation-study/presentation/time.png) / [PDF](../benchmarks/results/research-validation-20261016/original-generation-study/presentation/time.pdf).

| Method | Attempts | Hop, mean [95% CI] | Route, mean [95% CI] | 3D coverage, mean [95% CI] |
| --- | ---: | ---: | ---: | ---: |
| `nearest_3d` | 36 | 0.972 [0.917, 1.000] | 0.000 [0.000, 0.000] | 1.000 [1.000, 1.000] |
| `nearest_catalog` | 36 | 0.667 [0.333, 0.917] | 0.000 [0.000, 0.000] | 0.944 [0.889, 1.000] |
| `pipeline` | 36 | 0.444 [0.250, 0.639] | 0.444 [0.250, 0.639] | 0.380 [0.194, 0.565] |
| `reaction_enumeration` | 36 | 1.000 [1.000, 1.000] | 1.000 [1.000, 1.000] | 0.843 [0.796, 0.898] |
| `synllama_only` | 36 | 0.972 [0.917, 1.000] | 0.972 [0.917, 1.000] | 0.870 [0.713, 0.981] |

**Held-out-target split: 90 attempts, 3 target clusters.** [PNG](../benchmarks/results/research-validation-20261016/original-generation-study/presentation/target.png) / [PDF](../benchmarks/results/research-validation-20261016/original-generation-study/presentation/target.pdf).

| Method | Attempts | Hop, mean [95% CI] | Route, mean [95% CI] | 3D coverage, mean [95% CI] |
| --- | ---: | ---: | ---: | ---: |
| `nearest_3d` | 18 | 0.833 [0.500, 1.000] | 0.000 [0.000, 0.000] | 0.833 [0.500, 1.000] |
| `nearest_catalog` | 18 | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 0.778 [0.500, 1.000] |
| `pipeline` | 18 | 0.056 [0.000, 0.167] | 0.056 [0.000, 0.167] | 0.037 [0.000, 0.111] |
| `reaction_enumeration` | 18 | 0.944 [0.833, 1.000] | 0.944 [0.833, 1.000] | 0.537 [0.444, 0.667] |
| `synllama_only` | 18 | 0.500 [0.000, 1.000] | 0.667 [0.500, 1.000] | 0.444 [0.333, 0.500] |

These display values are rounded to three decimals; the linked recorded JSON and
[metrics CSV](../benchmarks/results/research-validation-20261016/original-generation-study/presentation/metrics.csv) preserve full numerical precision.
The figures separately show infrastructure identification bounds; the complete
[bounds](../benchmarks/results/research-validation-20261016/original-generation-study/presentation/infrastructure_bounds.csv) and
[availability](../benchmarks/results/research-validation-20261016/original-generation-study/presentation/availability.csv) tables retain that distinction.

The original `pipeline` produced successful outputs in only **6/36 scaffold
attempts** (19 timeouts, 10 no-candidate outcomes and one infrastructure-missing
outcome) and **1/18 held-out-target attempts** (17 no-candidate outcomes). These
results do not support pipeline superiority. They document limited realized
output under this frozen implementation and finite budget; later repairs and
different successful-score populations cannot turn these failures into successes.

### What was checked, and what remains unmeasured

The completed common-scoring record covers 1,264 returned candidate occurrences;
repeated molecules, queries and seeds are not independent new compounds.
Independent direct catalog/rule replay matched all 763 supplied-route occurrences;
501 candidate occurrences supplied no route and remain separate. Replay used the
same RDKit chemistry kernel, a direct enumeration cap of 4,096 per step and the
pinned catalog/rule records. It is bounded route consistency, not proof of all
reagent sources, conditions, chemical feasibility, experimental synthesis or yield.
It does not imply chemical validation of every catalog molecule or registered rule.
The separate common-forward audit is not part of the generation-statistics bundle.

Common 3D means and the original pipeline-minus-comparator conditional-score
contrasts can involve different successful attempts and candidates for each method.
The historical paired contrasts pair target means, not identical observed
query/seed or candidate populations. They cannot establish like-for-like 3D
improvements. The full-grid hop, route and output-slot metrics retain the common
attempt grid; successful-score subsets and runtime availability are separately
recorded. Model agreement and similarity are not calibrated probabilities,
preserved biological activity, synthesis success or reaction yield.

This is a post hoc descriptive evaluation. The bounded ChEMBL panels have the
assay-context limitations recorded above, and unknown complete checkpoint training
inventories prevent a clean pretrained holdout claim. The separate partial training
sample audit cannot resolve that absence. Reusing a forward predictor for selection
and evaluation can enrich agreement through the selection rule itself; it is not
independent experimental evidence. Later engineering repairs, index sensitivity
checks and current-pipeline studies do not retrospectively change the original
failed attempts or establish their full-pipeline equivalence.

### Reproduce the recorded statistics and presentation

The [publication validation record](../benchmarks/results/research-validation-20261016/original-generation-study/validation.json),
[SHA256 checksums](../benchmarks/results/research-validation-20261016/original-generation-study/SHA256SUMS),
[portable statistical bundle](../benchmarks/results/research-validation-20261016/original-generation-study/statistical-bundle.zip),
[recorded replay validation](../benchmarks/results/research-validation-20261016/original-generation-study/replay/validation.json), and
[presentation manifest](../benchmarks/results/research-validation-20261016/original-generation-study/presentation/manifest.json) bind the completed
inputs and outputs. The 55-file ZIP contains the original panels, protocols, run
lists, common-score records, historical analysis helper and hash-bound completion
reviews. It does not contain model weights or a runnable generation environment.
Original absolute paths remain historical provenance strings and are not opened
by statistical replay. Existing source attribution is retained without assigning
a new project or model license.

The parent replayed all three splits with Python 3.10.12 / NumPy 1.26.4,
Python 3.10.12 / NumPy 1.24.4 and Python 3.12.13 / NumPy 2.5.4; the three comparison
JSON files were byte-identical across these environments. This establishes stored
statistical replay under those runtimes, not chemical or neural-runtime equivalence.
An independent arithmetic implementation also checked the original per-query,
per-target, bootstrap, contrast and infrastructure records. The actual presentation
review checked all 21,440 metric rows, 1,405 contrast rows, 1,866 bounds rows and
2,565 availability rows, including preservation of 756 null numeric rows.

From the repository root, use fresh destinations:

```bash
STUDY_ARTIFACTS=benchmarks/results/research-validation-20261016/original-generation-study
(cd "$STUDY_ARTIFACTS" && sha256sum -c SHA256SUMS)
mkdir /tmp/original630-reproduction
python -m zipfile -e "$STUDY_ARTIFACTS/statistical-bundle.zip" /tmp/original630-reproduction

python -B benchmarks/diagnostics/generation_replay/replay.py verify \
  --bundle /tmp/original630-reproduction/original630-bundle \
  --expected-manifest-sha256 147e2b874974e2a99586570315bf72f9a0e97aac1647285e8a97533d3c9b300a \
  --output /tmp/original630-reproduction/verification.json

python -B benchmarks/diagnostics/generation_replay/replay.py replay \
  --bundle /tmp/original630-reproduction/original630-bundle \
  --expected-manifest-sha256 147e2b874974e2a99586570315bf72f9a0e97aac1647285e8a97533d3c9b300a \
  --output /tmp/original630-reproduction/recomputed

python -B benchmarks/diagnostics/generation_replay/export_figures.py \
  --replay "$STUDY_ARTIFACTS/replay" \
  --expected-validation-sha256 87a0210ea63fedd89d480850e5b205e26d1313aa6e66af258af21ef9dec1db94 \
  --output /tmp/original630-reproduction/presentation
```

The final command deliberately renders the published, externally anchored replay.
A new replay records its own runtime in `validation.json`; its validation-file hash
can therefore differ even when the three comparison JSONs match. Review and retain
that new validation hash before using it for a new presentation. The published
figures used Matplotlib 3.10.8; rendering bytes can vary across plotting runtimes.
The [replay contract](../benchmarks/diagnostics/generation_replay/README.md) and
[exporter contract](../benchmarks/diagnostics/generation_replay/EXPORT.md) describe
the trust anchors, fixed estimator and interpretation limits.

The portable statistical replay verifies the completed input grids, original
source/protocol identities, source and output hashes, and common-score joins
before recomputing recorded aggregates. It requires Python and NumPy; it performs
no generation, model inference, chemistry or common-3D recomputation. The separate
figure exporter accepts a passed replay and an externally retained validation SHA,
then exports its recorded values without new statistical calculations. Matplotlib
is needed only for presentation. All destinations must be fresh, and the external
hashes must come from the reviewed publication record rather than being trusted
solely because an unverified archive contains them.

The complete-precision CSVs preserve null values as blanks and include every
recorded metric, per-query and per-target summaries, target-macro estimates,
conditional-population contrast scope, and separate infrastructure bounds and
availability. The publication archive preserves existing source attribution; it
does not assign a new license to project code, models or third-party artifacts.

## Completed original common-forward audit

This separate audit covers the original 39-file study at
`d972124d587e647ef48455cf1a038e6be96e5b04`: all 630 planned attempts and 1,264
accepted candidate occurrences. It does not measure the later 41-file pipeline or
the supplementary 270-attempt study. The generation statuses remain 434 success, 140
no-candidate outcomes, 55 timeouts and one infrastructure-missing outcome. A
failed generation attempt has no accepted returns; the audit cannot promote it to
success because a diagnostic or partial file exists.

The common ReactionT5v2-forward predictor is also the synthesis selection-gate
model. Agreement after selection is therefore not independent experimental
validation. Missing reagents, conditions and empirical context remain limitations.
No-route controls show absent supplied evidence, not failed chemical synthesis.
This audit establishes neither activity, yield, synthesis-success probability nor
superiority of the pipeline or 3D scoring.

### All-method forward coverage and eligible fractions

The primary metric, `attempt_with_common_forward_matched_route`, is an indicator
that an accepted candidate supplies a nonempty exact replay-verified route with
all steps evaluated and matched by the common forward model. Its denominator
retains every planned attempt. The conditional metric,
`matched_fraction_given_evaluable_routes`, first divides matched route occurrences
by fully evaluated route occurrences **within each eligible attempt**, then averages
by query and target. It is a secondary selected-subset summary, not the pooled
fraction across routes and not a synthesis-success probability.

Each table keeps every method in alphabetical order. Primary and conditional
values are recorded target-macro means with descriptive, unadjusted 95%
target-cluster bootstrap CIs, rounded here to three decimals. S/E/N reports exact
candidate-occurrence counts with a supplied route / a fully evaluated route / no
supplied route. A/Q/T reports eligible attempts / queries / target clusters for
the conditional fraction. Planned query counts are 12, 12 and 6 for the scaffold,
publication-year and held-out-target splits. `NA` means no eligible observations; a single target
cluster has no CI. Complete-precision records remain in the CSVs and analysis JSON.

**Scaffold split: 360 attempts, 6 target clusters.**

| Method | Attempts | Primary matched-route indicator [95% CI] | S/E/N occurrences | Conditional fraction [95% CI] | Eligible A/Q/T |
| --- | ---: | ---: | ---: | ---: | ---: |
| `nearest_3d` | 36 | 0.000 [0.000, 0.000] | 0/0/78 | NA | 0/0/0 |
| `nearest_catalog` | 36 | 0.000 [0.000, 0.000] | 0/0/108 | NA | 0/0/0 |
| `pipeline` | 36 | 0.167 [0.028, 0.333] | 15/15/0 | 1.000 [1.000, 1.000] | 6/4/3 |
| `pipeline_encounter_allocation` | 36 | 0.417 [0.278, 0.556] | 34/34/0 | 1.000 [1.000, 1.000] | 15/7/6 |
| `pipeline_legacy_search` | 36 | 0.361 [0.139, 0.639] | 37/37/0 | 1.000 [1.000, 1.000] | 13/6/5 |
| `pipeline_no_3d` | 36 | 0.806 [0.667, 0.944] | 82/82/0 | 1.000 [1.000, 1.000] | 29/11/6 |
| `pipeline_no_refinement` | 36 | 0.444 [0.306, 0.556] | 43/43/0 | 1.000 [1.000, 1.000] | 16/7/6 |
| `pipeline_rule_only` | 36 | 0.306 [0.139, 0.500] | 29/29/0 | 0.817 [0.683, 0.950] | 11/6/5 |
| `reaction_enumeration` | 36 | 0.556 [0.417, 0.694] | 90/90/0 | 0.375 [0.241, 0.528] | 30/11/6 |
| `synllama_only` | 36 | 0.611 [0.333, 0.861] | 90/90/0 | 0.653 [0.375, 0.875] | 30/12/6 |

**Publication-year split: 180 attempts, 6 target clusters.**

| Method | Attempts | Primary matched-route indicator [95% CI] | S/E/N occurrences | Conditional fraction [95% CI] | Eligible A/Q/T |
| --- | ---: | ---: | ---: | ---: | ---: |
| `nearest_3d` | 36 | 0.000 [0.000, 0.000] | 0/0/108 | NA | 0/0/0 |
| `nearest_catalog` | 36 | 0.000 [0.000, 0.000] | 0/0/108 | NA | 0/0/0 |
| `pipeline` | 36 | 0.444 [0.250, 0.639] | 41/41/0 | 1.000 [1.000, 1.000] | 16/7/5 |
| `reaction_enumeration` | 36 | 0.722 [0.556, 0.861] | 108/108/0 | 0.417 [0.306, 0.500] | 36/12/6 |
| `synllama_only` | 36 | 0.667 [0.444, 0.861] | 105/105/0 | 0.593 [0.343, 0.843] | 35/12/6 |

**Held-out-target split: 90 attempts, 3 target clusters.**

| Method | Attempts | Primary matched-route indicator [95% CI] | S/E/N occurrences | Conditional fraction [95% CI] | Eligible A/Q/T |
| --- | ---: | ---: | ---: | ---: | ---: |
| `nearest_3d` | 18 | 0.000 [0.000, 0.000] | 0/0/45 | NA | 0/0/0 |
| `nearest_catalog` | 18 | 0.000 [0.000, 0.000] | 0/0/54 | NA | 0/0/0 |
| `pipeline` | 18 | 0.056 [0.000, 0.167] | 2/2/0 | 1.000 [CI unavailable] | 1/1/1 |
| `reaction_enumeration` | 18 | 0.667 [0.500, 0.833] | 51/51/0 | 0.380 [0.194, 0.500] | 17/6/3 |
| `synllama_only` | 18 | 0.333 [0.167, 0.500] | 36/36/0 | 0.444 [0.111, 0.667] | 12/4/3 |

For the original `pipeline`, conditional agreement is 1.000, but the eligible
population is only 6/36 scaffold attempts and 1/18 held-out-target attempts.
The target conditional value comes from one target cluster and has no CI.
These selected populations do not support a general synthesis-success claim.
Retrieval methods have no supplied routes: their primary coverage is zero and
their conditional fraction is missing, even when generation returned molecules.

[Method metrics](../benchmarks/results/research-validation-20261016/original-common-forward-study/method_metrics.csv),
[coverage/status counts](../benchmarks/results/research-validation-20261016/original-common-forward-study/route_coverage.csv),
[paired metrics](../benchmarks/results/research-validation-20261016/original-common-forward-study/paired_metrics.csv) and
[infrastructure bounds](../benchmarks/results/research-validation-20261016/original-common-forward-study/infrastructure_bounds.csv) are direct projections
of recorded values, with nulls left blank and no additional statistical calculation.
Method/paired CSVs include all 19 target-macro metrics, their eligibility and the
secondary complete-case population. Full per-query and per-target records remain
in [analysis.json](../benchmarks/results/research-validation-20261016/original-common-forward-study/analysis.json).

Unlike the original generation target-mean contrasts above, the forward audit
forms pipeline-minus-comparator differences at **identical query/seed slots before
aggregation**. A conditional contrast requires both observations to be eligible;
its pair eligibility is preserved. This pairing does not make selected conditional
populations representative or establish an independent chemical experiment.
The fixed hierarchy averages seeds within queries, queries within targets and
targets equally, with 10,000 bootstrap draws and seed 20261016. Only 6/6/3 target
clusters support the three splits. Reused routes, reactant/context inputs and
targets create within- and cross-target dependence; repeated splits are not new
independent cohorts. Intervals are descriptive and unadjusted, not significance
or superiority tests.

The infrastructure interruption retains an unknown counterfactual outcome.
Recorded lower/upper identification bounds are **not CIs**. Secondary complete
cases exclude only infrastructure-missing slots, retain other terminal failures,
and are selection-biased; paired complete cases require both outcomes observed.
The post hoc analysis plan was fixed before comparative analysis while generation
and the audit were unfinished. It is not preregistered confirmatory evidence.
Unknown complete model-training inventories and the original panels' assay/source
limitations remain unresolved; model selection/evaluation reuse cannot establish
clean holdout or independent experimental validation.

### Forward replay artifacts and verification scope

The [18-file portable ZIP](../benchmarks/results/research-validation-20261016/original-common-forward-study/statistical-bundle.zip),
[checksums](../benchmarks/results/research-validation-20261016/original-common-forward-study/SHA256SUMS),
[publication validation](../benchmarks/results/research-validation-20261016/original-common-forward-study/validation.json) and
[copied replay validation](../benchmarks/results/research-validation-20261016/original-common-forward-study/replay-validation.json) preserve the reviewed
normalized attempts, occurrences, original analysis and provenance. Existing source
attribution is retained; no new project or model license is assigned. Historical
absolute path strings remain metadata and are never opened during portable replay.
Model weights, raw per-route predictions, catalogs and generation environments are
not included. The original generation-statistics bundle remains separate.

The parent independently checked accepted-output joins, saved raw prediction
interpretation and route records using the same RDKit kernel, without rerunning
the neural model. Its separate arithmetic implementation checked all **52,514
numeric cells**, including exact query/seed conditional pairing, eligibility,
bounds and dependence diagnostics. Portable replay then reproduced the original
`analysis.json` byte-for-byte under Python 3.10.12 / NumPy 1.26.4,
Python 3.10.12 / NumPy 1.24.4 and Python 3.12.13 / NumPy 2.5.4. This validates
stored statistical replay; it does not reproduce neural inference or substitute
for the raw integrity assessment.

From the repository root, use fresh destinations:

```bash
FORWARD_ARTIFACTS=benchmarks/results/research-validation-20261016/original-common-forward-study
(cd "$FORWARD_ARTIFACTS" && sha256sum -c SHA256SUMS)
mkdir /tmp/original630-forward-reproduction
python -m zipfile -e "$FORWARD_ARTIFACTS/statistical-bundle.zip" /tmp/original630-forward-reproduction

python -B benchmarks/diagnostics/forward_replay/replay.py verify \
  --bundle /tmp/original630-forward-reproduction/original630-forward-bundle \
  --expected-manifest-sha256 9650f71cdd4a7e2de9f8b90a839be9477153fef94a938325e033b924572b05c8 \
  --output /tmp/original630-forward-reproduction/verification.json

python -B benchmarks/diagnostics/forward_replay/replay.py replay \
  --bundle /tmp/original630-forward-reproduction/original630-forward-bundle \
  --expected-manifest-sha256 9650f71cdd4a7e2de9f8b90a839be9477153fef94a938325e033b924572b05c8 \
  --output /tmp/original630-forward-reproduction/recomputed

cmp "$FORWARD_ARTIFACTS/analysis.json" /tmp/original630-forward-reproduction/recomputed/analysis.json
```

Python and NumPy suffice: no chemistry/model imports, network, GPU or new inference
are required. The [replay contract](../benchmarks/diagnostics/forward_replay/README.md)
describes the pinned statistical extraction and strict full-grid, normalization
and hash gates. Retain the external manifest SHA from this reviewed publication;
a hash obtained only from an untrusted archive does not authenticate that archive.
Runtime metadata can change a new `replay-validation.json` hash even when its
`analysis.json` matches the published bytes.

## Publication interpretation

A separately selected [same-assay follow-up](../benchmarks/results/assay-matched-controls-20261016/README.md) now contains 197 pairs, 12 queries and 12 assays across those same six targets. Selection/scoring criteria were fixed before its chemical scores were calculated, after the initial heterogeneity findings. Three pairs failed the stereo policy; all failures remain in the denominators and bounds. Paired complete-case target-macro ROC AUC was 0.7392 for bit Tanimoto, 0.7370 for shape and 0.7131 for joint 3D. The paired joint-minus-bit difference was −0.0261 (target-bootstrap 95% interval −0.1017 to +0.0402), so this evaluation does not establish additional 3D discrimination. All hierarchical ROC/AP estimates and missing-score bounds matched independent scikit-learn calculations within `1e-12`. Reconstruction from the published raw ZIP reproduced all chemical records and selected pairs. Same assay ID narrows context without resolving every metadata ambiguity or establishing new-molecule activity; changed cohorts prevent attributing differences from earlier panels to assay matching alone.

[SynLlama](https://pubs.acs.org/doi/10.1021/acscentsci.5c01285) and [RuSH](https://pubs.acs.org/doi/10.1021/acs.jcim.5c00029) already address synthesis-aware generation and/or combined 2D/3D scaffold hopping. A new paper needs a defensible improvement under a stated setting, supported by independent comparisons and ablations. Known active scaffold recovery remains a structural proxy; activity of a new candidate and chemical reaction yield remain unmeasured without additional evidence. Experimental collaboration can strengthen an application claim, but the computational-method evaluation should stand on its own with declared limits and reproducible artifacts.

For a separate same-assay evaluation, use the [assay-matched control protocol](ASSAY_MATCHED_CONTROLS.md). It audits all retained activities before score-blind selection and preserves assay-local labels; it does not replace the frozen generation benchmark.
