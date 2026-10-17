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

## Completed supplementary comparator study

This separate **270-attempt study** uses the fixed 41-file scientific source at
`2328588d0969369e694d9293a1dccc89392a8203`. It is not the original 630-attempt
study at `d972124d587e647ef48455cf1a038e6be96e5b04`; their results are neither pooled
nor a causal old-versus-new improvement estimate. Each of three methods used the
same 12 scaffold, 12 publication-year and 6 held-out-target public query slots,
with seeds 42/43/44. The three splits have 6/6/3 target clusters; repeated targets,
inputs, routes and cached exact pairs are not independent replications.

**The current pipeline was not superior in this comparison.** It delivered hops
less often and had lower common-3D output coverage and lower conditional common
shape, pharmacophore and joint means than both comparators in each split. All six
paired joint-score point estimates were negative; both publication-year intervals
include zero. Selected conditional populations and few target clusters preclude
sweeping superiority or equivalence conclusions in either direction.

There were **227 success, 37 no-candidate and 6 worker-error outcomes**, with no
retries or backfill. Pipeline success was 47/90: 19/36 scaffold, 24/36 publication-year
and 4/18 held-out-target attempts. Its six errors were the strict unspecified-stereo
reference preflight rejections for `CHEMBL205-query-2` (scaffold) and
`CHEMBL220-query-2` (held-out-target), each at all three seeds. They remain failures
in the planned denominator. CReM and Mol2Mol each returned valid outputs in all
90 attempts, with 270 valid candidate occurrences per method.

A separate [post hoc native-stage diagnostic](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/pipeline-failure-diagnostic.json)
joins all 90 pipeline slots to their saved requests and 84 available synthesis
summaries. Of the 37 no-candidate outcomes, 12 had all capped design references
rejected by stereochemistry preflight, 17 had an eligible design but no verified
route, and 8 had a verified route but no final candidate. The six query-preflight
worker errors are separate. Native no-candidate rejection counters record 1,070
`evidence_rejected`, 2 `constraint_rejected` and 20 `3d_failed` **search events**,
which can repeat products or routes. They are not unique molecules, causal effects,
synthesis probabilities or comparable method-level denominators. This diagnostic
changes none of the 69 metrics, original stage records or planned failures.

### All-method delivery, routes and common-3D coverage

Methods are ordered alphabetically, not by their results. S/N/E are exact success /
no-candidate / worker-error counts. All nine combinations have zero timeout,
budget-exceeded and infrastructure-interruption counts. Hop and rule-route
indicators retain every planned attempt. A hop here requires a ring scaffold,
a different chemical family and a different normalized Murcko scaffold family;
it does not require generic ring topology replacement or measured activity.
Rule-route coverage means at least one supplied catalog/rule-verified route,
not experimental synthesis. Common-3D slot coverage is the number of returned
candidates with all three common-scoring seeds successful divided by **cap 3**
within each attempt; it is not conditional similarity quality.

Values below are recorded target-macro means [descriptive 95% target-cluster
bootstrap CI], displayed to three decimals. The hierarchy averages seeds within
queries, queries within targets, then targets equally. There are 10,000 bootstrap
draws with seed 20261016. Historical JSON `independent_targets` fields count target
clusters; they do not demonstrate statistical independence. Intervals are
unadjusted descriptions, not significance tests.

| Split | Method | Attempts | S/N/E | Hop indicator [95% CI] | Rule-route indicator [95% CI] | Common-3D slot coverage [95% CI] |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Scaffold | `crem_train_partial_exo` | 36 | 36/0/0 | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 0.852 [0.685, 0.991] |
| Scaffold | `mol2mol_pretrained` | 36 | 36/0/0 | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 0.833 [0.667, 0.972] |
| Scaffold | `pipeline_current` | 36 | 19/14/3 | 0.528 [0.333, 0.750] | 0.528 [0.333, 0.750] | 0.426 [0.241, 0.676] |
| Publication-year | `crem_train_partial_exo` | 36 | 36/0/0 | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 1.000 [1.000, 1.000] |
| Publication-year | `mol2mol_pretrained` | 36 | 36/0/0 | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 1.000 [1.000, 1.000] |
| Publication-year | `pipeline_current` | 36 | 24/12/0 | 0.667 [0.333, 0.917] | 0.667 [0.333, 0.917] | 0.519 [0.278, 0.731] |
| Held-out-target | `crem_train_partial_exo` | 18 | 18/0/0 | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 0.796 [0.500, 1.000] |
| Held-out-target | `mol2mol_pretrained` | 18 | 18/0/0 | 1.000 [1.000, 1.000] | 0.000 [0.000, 0.000] | 0.667 [0.500, 1.000] |
| Held-out-target | `pipeline_current` | 18 | 4/11/3 | 0.222 [0.000, 0.500] | 0.222 [0.000, 0.500] | 0.204 [0.000, 0.500] |

Of 653 accepted valid candidate occurrences, 113 supplied routes and 540 supplied
no route. All 113 supplied routes matched direct parent catalog/rule replay:
58 distinct serialized route/product inputs, 137 step occurrences and 250 catalog
leaf occurrences within the assessed 4,096-product step cap.
Those are not 58 independent syntheses or exhaustive validation of every catalog
molecule or registered rule. CReM graph edits are not synthetic routes; absent
routes in either comparator are unassessed evidence, not chemical failure.
Declared pipeline forward evidence uses its own selection-gate predictor.
**No fresh uniform forward re-inference was performed for these 270 attempts.**
The separate original-630 common-forward audit does not supply that missing check.
Neither route replay nor forward agreement establishes synthesis probability,
yield, reagent/condition completeness or experimental success.

### Conditional similarity and exact paired comparisons

The common scorer evaluated 252 unique exact query/candidate pairs, reused across
653 candidate occurrences. All three seeds succeeded for 590 occurrences; 63
failed or unassessed occurrences retain their reasons and missing scores. A
successful pair contributes its recorded median over seeds. The attempt then
averages eligible candidate occurrences before the same query/target hierarchy.
Seed spreads describe algorithmic stability, not calibrated uncertainty.

Binary Tanimoto uses Morgan radius 2, 2,048 bits and chirality; sparse-count
Tanimoto uses radius 2 without chirality. These are distinct representations,
not interchangeable thresholds. The common ligand-only joint score also retains
the pinned energy/strain term; it is neither pure shape similarity nor a calibrated
probability. No universal activity cutoff is inferred. Structural active-scaffold
recovery uses scaffold intersection with eligible references, without a Tanimoto
0.5 gate; it differs from the earlier prefilter-recall endpoint.

A/Q/T below reports eligible attempts / queries / target clusters separately for
2D and 3D. Means are rounded only for display; the complete-precision CSV includes
all recorded CIs and eligibility. Missing values are `NA`, never similarity zero.

| Split | Method | Bit T | Count T | Eligible 2D A/Q/T | Shape | Pharmacophore | Joint [95% CI] | Eligible 3D A/Q/T |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Scaffold | `crem_train_partial_exo` | 0.705 | 0.740 | 36/12/6 | 0.809 | 0.554 | 0.624 [0.538, 0.703] | 33/11/6 |
| Scaffold | `mol2mol_pretrained` | 0.560 | 0.569 | 36/12/6 | 0.794 | 0.579 | 0.640 [0.521, 0.751] | 33/11/6 |
| Scaffold | `pipeline_current` | 0.347 | 0.392 | 19/7/6 | 0.646 | 0.330 | 0.427 [0.322, 0.520] | 19/7/6 |
| Publication-year | `crem_train_partial_exo` | 0.828 | 0.839 | 36/12/6 | 0.761 | 0.530 | 0.585 [0.499, 0.676] | 36/12/6 |
| Publication-year | `mol2mol_pretrained` | 0.559 | 0.580 | 36/12/6 | 0.735 | 0.549 | 0.585 [0.495, 0.676] | 36/12/6 |
| Publication-year | `pipeline_current` | 0.463 | 0.472 | 24/8/5 | 0.671 | 0.372 | 0.484 [0.394, 0.613] | 24/8/5 |
| Held-out-target | `crem_train_partial_exo` | 0.725 | 0.745 | 18/6/3 | 0.761 | 0.496 | 0.568 [0.394, 0.750] | 15/5/3 |
| Held-out-target | `mol2mol_pretrained` | 0.511 | 0.530 | 18/6/3 | 0.717 | 0.467 | 0.560 [0.395, 0.745] | 12/4/3 |
| Held-out-target | `pipeline_current` | 0.469 | 0.466 | 4/2/2 | 0.663 | 0.317 | 0.424 [0.391, 0.456] | 4/2/2 |

Method conditional means use different selected populations. The paired summaries
instead subtract comparator from `pipeline_current` at **identical query/seed
slots with both values observed**, before seed/query/target aggregation. They are
not subtraction of the displayed method means. All-slot coverage contrasts retain
failures; conditional joint contrasts exclude unavailable scores in either arm.
This differs from the original-630 generation target-mean contrast and does not
remove selection bias.

| Split | Comparator | Paired hop difference [95% CI] | Paired common-3D coverage difference [95% CI] | Paired joint difference [95% CI] | Joint eligible A/Q/T |
| --- | --- | ---: | ---: | ---: | ---: |
| Scaffold | `crem_train_partial_exo` | -0.472 [-0.667, -0.250] | -0.426 [-0.630, -0.204] | -0.175 [-0.246, -0.117] | 19/7/6 |
| Scaffold | `mol2mol_pretrained` | -0.472 [-0.667, -0.250] | -0.407 [-0.611, -0.204] | -0.217 [-0.420, -0.076] | 19/7/6 |
| Publication-year | `crem_train_partial_exo` | -0.333 [-0.667, -0.083] | -0.481 [-0.722, -0.269] | -0.048 [-0.138, 0.046] | 24/8/5 |
| Publication-year | `mol2mol_pretrained` | -0.333 [-0.667, -0.083] | -0.481 [-0.722, -0.269] | -0.064 [-0.130, 0.003] | 24/8/5 |
| Held-out-target | `crem_train_partial_exo` | -0.778 [-1.000, -0.500] | -0.593 [-1.000, 0.000] | -0.068 [-0.134, -0.003] | 4/2/2 |
| Held-out-target | `mol2mol_pretrained` | -0.778 [-1.000, -0.500] | -0.463 [-0.889, 0.000] | -0.034 [-0.064, -0.004] | 4/2/2 |

### Scaffold retention is a separate structural diagnostic

A Murcko change can retain the query core. The secondary **core-retaining scaffold
expansion** annotation requires the query normalized Murcko scaffold to be an
atom/bond-aware subgraph of a larger candidate scaffold, in addition to the
existing hop definition. **Query scaffold not contained** requires that directed
query-within-candidate test to fail for a hop. The match is asymmetric and uses
`useChirality=False`; generic ring-skeleton containment is separately recorded.
Neither annotation establishes pharmacophore retention, target activity or a
universally accepted scaffold-hop class.

These attempt indicators can **both be true for the same attempt**, because
different returned candidates can meet the two definitions. They are not mutually
exclusive attempt categories and must not be summed as a partition. The annotation
was specified after a fixed-example diagnostic, before supplementary comparative
output, without changing primary eligibility, ranking or budgets.
All 653 retained occurrences were assessable for these annotations; missing
annotations would remain null and separately counted.

| Split | Method | Attempt with core-retaining expansion [95% CI] | Attempt with query scaffold not contained [95% CI] |
| --- | --- | ---: | ---: |
| Scaffold | `crem_train_partial_exo` | 0.333 [0.167, 0.500] | 0.944 [0.833, 1.000] |
| Scaffold | `mol2mol_pretrained` | 0.167 [0.000, 0.333] | 0.917 [0.750, 1.000] |
| Scaffold | `pipeline_current` | 0.000 [0.000, 0.000] | 0.528 [0.333, 0.750] |
| Publication-year | `crem_train_partial_exo` | 0.083 [0.000, 0.250] | 1.000 [1.000, 1.000] |
| Publication-year | `mol2mol_pretrained` | 0.417 [0.083, 0.750] | 1.000 [1.000, 1.000] |
| Publication-year | `pipeline_current` | 0.083 [0.000, 0.250] | 0.667 [0.333, 0.917] |
| Held-out-target | `crem_train_partial_exo` | 0.500 [0.000, 1.000] | 0.667 [0.500, 1.000] |
| Held-out-target | `mol2mol_pretrained` | 0.333 [0.000, 0.500] | 0.667 [0.500, 1.000] |
| Held-out-target | `pipeline_current` | 0.000 [0.000, 0.000] | 0.222 [0.000, 0.500] |

### Scope, budgets and interpretation

The [CReM baseline](../benchmarks/external/README-crem.md) uses its declared train-only
fragment library and partial/exocyclic mutation settings. The
[official REINVENT4 Mol2Mol baseline](../benchmarks/external/mol2mol/README.md) is plain
pretrained sampling with declared selection, not full reinforcement-learning
REINVENT or RuSH reproduction. Complete Mol2Mol and other model pretraining
inventories remain unknown. A public split or train-only CReM construction does
not prove a clean holdout across all models and sources.

The native work caps differ: up to 12 initial synthesis product/design 3D screens
for the pipeline, with refinement separately bounded and reported, versus 256
CReM API-yielded products or 256 Mol2Mol sampled sequences. The pinned preparation
plan retains its historical "proposal evaluations" label; the actual adapter
counts initial 3D screens, not generated proposals. Each attempt has the declared 360-second
cold-start wall limit, 24,000 MB memory limit and output cap 3. These are resource
constraints, not equal amounts of scientific search or universal compute equality.
Shared-host process-tree RSS/CPU measurements and nonuniform GPU-memory observation
are not a controlled hardware performance comparison. Full native/raw accounting
is retained; unavailable pipeline raw counts are not invented.

This is a post hoc computational-method study. The recorded analysis/annotation
plans were fixed before their comparative calculation; that is not preregistration.
There are only 6/6/3 target clusters, with shared inputs and cross-split dependence.
Infrastructure identification bounds are **not CIs**, and infrastructure-only
complete-case sensitivity can be selection-biased. No infrastructure outcome is
missing in this supplementary grid, so these stored sensitivities coincide with
their observed counterparts. Other failures remain in full-attempt denominators;
conditional scores and unbounded resource metrics are never zero-imputed.

The negative delivery and similarity findings remain part of the conclusion.
Structural recovery and 3D proxies do not establish new-compound activity,
experimental synthesizability, 3D superiority or journal readiness. Experimental
collaboration remains undecided; a computational-method claim must stand on the
stated evidence and limits. The original frozen39, fixed current41 and separate
original common-forward records remain distinct.

### Supplementary artifacts and reproduction

The [291-file portable ZIP](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/statistical-bundle.zip),
[checksums](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/SHA256SUMS), [publication validation](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/validation.json)
and [copied replay validation](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/replay-validation.json) bind the completed
study. The [final comparison](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/analysis/comparison.json) and all eight
companion analysis files are byte-preserving copies. [Method metrics](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/method_metrics.csv),
[paired metrics](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/paired_metrics.csv), [availability](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/availability.csv),
[separate infrastructure bounds](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/infrastructure_bounds.csv) and
[attempt/status counts](../benchmarks/results/research-validation-20261016/supplementary-comparator-study/attempt_status.csv) directly project stored values.
The first three CSVs include all 69 target-macro metrics and both primary and
secondary complete-case populations. No new estimates were calculated for this
publication. Numeric precision is preserved and nulls are blank; complete
per-query/per-target and per-query seed availability remain in the JSON.

Parent assessment checked all **92,398 statistical numeric cells**, full normalized
joins for 653 candidates, 393 independent graph checks including controls, and
1,306 directly calculated bit/count fingerprint values with maximum error zero.
The direct route/fingerprint checks and graph checks retain their shared RDKit
chemistry kernel and pinned identity-policy limitations. Common-score cache
integrity and saved-seed arithmetic were checked separately; these checks do not
turn cached values into fresh chemical or experimental validation.

The stored comparison replay was then byte-identical in Python 3.10.12 / NumPy
1.26.4, Python 3.10.12 / NumPy 1.24.4 and Python 3.12.13 / NumPy 2.5.4. The bundle
retains normalized occurrences, common-score caches and failures, original plans,
extracted statistics and parent evidence. Historical absolute paths are inert
provenance; model weights, live runtimes and full catalog assets are not bundled.
Existing source attribution is preserved; no new license is assigned.

From the repository root, use fresh destinations:

```bash
SUPPLEMENT_ARTIFACTS=benchmarks/results/research-validation-20261016/supplementary-comparator-study
(cd "$SUPPLEMENT_ARTIFACTS" && sha256sum -c SHA256SUMS)
mkdir /tmp/supplementary270-reproduction
python -m zipfile -e "$SUPPLEMENT_ARTIFACTS/statistical-bundle.zip" /tmp/supplementary270-reproduction

python -B benchmarks/diagnostics/supplement_replay/replay.py verify \
  --bundle /tmp/supplementary270-reproduction/supplementary270-bundle \
  --expected-manifest-sha256 9ed59273170348c0b910957e6c43989313a64143acdb6e8b0a24fe66a0332b69 \
  --output /tmp/supplementary270-reproduction/verification.json

python -B benchmarks/diagnostics/supplement_replay/replay.py replay \
  --bundle /tmp/supplementary270-reproduction/supplementary270-bundle \
  --expected-manifest-sha256 9ed59273170348c0b910957e6c43989313a64143acdb6e8b0a24fe66a0332b69 \
  --output /tmp/supplementary270-reproduction/recomputed

cmp "$SUPPLEMENT_ARTIFACTS/analysis/comparison.json" /tmp/supplementary270-reproduction/recomputed/comparison.json
```

The [replay contract](../benchmarks/diagnostics/supplement_replay/README.md) requires
Python and NumPy only and rejects incomplete grids, stale hashes and unsafe inputs.
It performs stored statistical replay, not new model/chemistry/3D calculation.
Keep the externally reviewed manifest SHA as the trust anchor. Runtime metadata
can change a new validation-file hash even when comparison bytes match.

## Publication interpretation

A separately selected [same-assay follow-up](../benchmarks/results/assay-matched-controls-20261016/README.md) now contains 197 pairs, 12 queries and 12 assays across those same six targets. Selection/scoring criteria were fixed before its chemical scores were calculated, after the initial heterogeneity findings. Three pairs failed the stereo policy; all failures remain in the denominators and bounds. Paired complete-case target-macro ROC AUC was 0.7392 for bit Tanimoto, 0.7370 for shape and 0.7131 for joint 3D. The paired joint-minus-bit difference was −0.0261 (target-bootstrap 95% interval −0.1017 to +0.0402), so this evaluation does not establish additional 3D discrimination. All hierarchical ROC/AP estimates and missing-score bounds matched independent scikit-learn calculations within `1e-12`. Reconstruction from the published raw ZIP reproduced all chemical records and selected pairs. Same assay ID narrows context without resolving every metadata ambiguity or establishing new-molecule activity; changed cohorts prevent attributing differences from earlier panels to assay matching alone.

[SynLlama](https://pubs.acs.org/doi/10.1021/acscentsci.5c01285) and [RuSH](https://pubs.acs.org/doi/10.1021/acs.jcim.5c00029) already address synthesis-aware generation and/or combined 2D/3D scaffold hopping. A new paper needs a defensible improvement under a stated setting, supported by independent comparisons and ablations. Known active scaffold recovery remains a structural proxy; activity of a new candidate and chemical reaction yield remain unmeasured without additional evidence. Experimental collaboration can strengthen an application claim, but the computational-method evaluation should stand on its own with declared limits and reproducible artifacts.

For a separate same-assay evaluation, use the [assay-matched control protocol](ASSAY_MATCHED_CONTROLS.md). It audits all retained activities before score-blind selection and preserves assay-local labels; it does not replace the frozen generation benchmark.

### Opt-in allocation method binding (v6 / protocol v4)

`qualified-benchmark-protocol-v4` admits only the exact packaged
`qualified-selective-engineering-v6` artifact hash. Profile v6 combines ranked
original-hit screening with design/scaffold evidence rotation; protocol v4 is the
benchmark envelope version, not a numerical profile version. The artifact cites
prior implementation review and the accepted two-input development comparison;
those records are not self-certification of this release or evidence of general
superiority. The combination did not increase qualified yield over ranked
screening alone on those inputs.

The adapter still requires `engineering_policy_matched`; there is no generic
`custom_unvalidated` bypass. Protocol v4 preserves protocol v3's explicit native
DB selection and cross-phase request binding, including rejection of unbound
replacement libraries. The runner still pins all required source and assets and
uses the existing fixed budgets. Selected and empty outcomes traverse the same
selective collector validation, including all three seed witnesses and JSON,
TSV and all four SDF exports. Legacy None/v1, protocol v2 -> profile v2 and protocol
v3 -> profile v5 bindings remain unchanged; the recommendation default remains v5.
