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

## Publication interpretation

A separately selected [same-assay follow-up](../benchmarks/results/assay-matched-controls-20261016/README.md) now contains 197 pairs, 12 queries and 12 assays across those same six targets. Selection/scoring criteria were fixed before its chemical scores were calculated, after the initial heterogeneity findings. Three pairs failed the stereo policy; all failures remain in the denominators and bounds. Paired complete-case target-macro ROC AUC was 0.7392 for bit Tanimoto, 0.7370 for shape and 0.7131 for joint 3D. The paired joint-minus-bit difference was −0.0261 (target-bootstrap 95% interval −0.1017 to +0.0402), so this evaluation does not establish additional 3D discrimination. All hierarchical ROC/AP estimates and missing-score bounds matched independent scikit-learn calculations within `1e-12`. Reconstruction from the published raw ZIP reproduced all chemical records and selected pairs. Same assay ID narrows context without resolving every metadata ambiguity or establishing new-molecule activity; changed cohorts prevent attributing differences from earlier panels to assay matching alone.

[SynLlama](https://pubs.acs.org/doi/10.1021/acscentsci.5c01285) and [RuSH](https://pubs.acs.org/doi/10.1021/acs.jcim.5c00029) already address synthesis-aware generation and/or combined 2D/3D scaffold hopping. A new paper needs a defensible improvement under a stated setting, supported by independent comparisons and ablations. Known active scaffold recovery remains a structural proxy; activity of a new candidate and chemical reaction yield remain unmeasured without additional evidence. Experimental collaboration can strengthen an application claim, but the computational-method evaluation should stand on its own with declared limits and reproducible artifacts.
