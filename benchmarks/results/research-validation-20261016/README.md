# Recorded control evaluation — 2026-10-16

These records assess similarity discrimination and prefilter loss on retrospective labelled controls. They do not establish improved generation, biological activity of new molecules, or experimental synthesizability. At this record, the 360-run main generation comparison was in progress and 270 secondary runs were queued.

## Scope and provenance

The strict scaffold panel contains 5,567 target–chemical-family records, with two public queries for each of EGFR (CHEMBL203), CA2 (CHEMBL205), PIM1 (CHEMBL2147), AChE (CHEMBL220), A2a (CHEMBL251) and JNK3 (CHEMBL2637). Its SHA-256 is `e87a3b539b8c41dde4887a7d52ec67b29a7a63f688584b579272398ef6522204`. Scaffold assignment uses seed 42 and the versioned chemical-family policy in the [protocol](../../../docs/RESEARCH_VALIDATION.md).

Source data are from EMBL-EBI [ChEMBL](https://www.ebi.ac.uk/chembl/), release 37, retained on 2026-10-15. The bounded acquisition selected at most the first 1,000 activity IDs per target and publication-year band (≤2010, 2011–2020, >2020); it is not a complete target inventory or a random population sample. ChEMBL-derived records and adaptations retain [CC BY-SA 3.0 attribution](https://creativecommons.org/licenses/by-sa/3.0/); this data license does not change the repository's code license. Changes comprise assay filtering, chemical-family grouping, label assignment, partitioning and computed similarity summaries. See the [ChEMBL database citation](https://doi.org/10.1093/nar/gkad1004) and [41-source URL/hash manifest](source_provenance.json).

Included measurements require direct human binding-assay metadata: confidence 9, relationship D, exact target assignment, *Homo sapiens* assay organism and an explicitly null variant-sequence annotation. These metadata checks do not prove wild-type constructs or comparable experimental conditions. IC50 is used for EGFR/PIM1/AChE/JNK3; Ki for CA2/A2a. Active means ≤1,000 nM and inactive ≥10,000 nM, respecting censoring; conflicting and ambiguous labels are excluded. Assays are still heterogeneous. Inputs are connected molecules with at most 60 heavy atoms.

The control sample contains **240 query–control pairs, representing 120 compound IDs and 119 distinct exact structures**: 10 active and 10 inactive controls per query, deterministically selected with seed 42 from held-out families with a different scaffold from the query. The same compounds recur across a target's two queries. There are no decoys in this record. Balanced capped sampling does not reproduce deployment prevalence.

## Results and analysis

Full-population recall counts a held-out active scaffold as surviving when at least one member passes count Tanimoto ≥0.5. This uses every eligible active reference before control subsampling. Mean recall across queries within each target, then equally across six targets, was **0.0316201**; the independent target-bootstrap 95% interval was **0.0114856–0.0518026**.

| Score | ROC AUC on the same complete cases |
| --- | ---: |
| Bit Morgan Tanimoto | 0.6389815 |
| Count Morgan Tanimoto | 0.6365741 |
| 3D shape | 0.5415741 |
| 3D pharmacophore | 0.5639815 |
| Joint 3D | 0.5500000 |

The 3D protocol requests eight conformers per seed (42, 43, 44), 500 optimization iterations, energy window 5, pharmacophore sigma 1 and strain penalty 0.05; unspecified stereochemistry is rejected. Query conformers are generated; no bound pose or receptor is supplied. Scores are seed medians only when all three seeds succeed. Shape, pharmacophore and joint scores therefore share the same complete cases.

**42/240 pairs failed**, all with the specified-stereochemistry policy reason, producing 126 failed seed trials. Complete-case comparisons use the remaining 198 pairs across 11 queries. CA2 query 2 has no complete cases; its cage-stereo flag is retained in [query_scope.json](query_scope.json), and its 20 failed pairs remain in failure totals and all-control bounds. Within-target complete-case averages use available query estimates, so the CA2 estimate uses one query. Complete-case selection can bias discrimination estimates.

ROC AUC and average precision are recomputed per query with scikit-learn, then averaged within target and equally across targets. All recomputed ROC/AP estimates and missing-score bounds matched the framework within `1e-12`. Bootstrap intervals use 10,000 resamples of the six target means, seed 20261016; compounds and seeds are not treated as independent target replicates. The interval above comes from [control_analysis.json](control_analysis.json); the original controls summary uses a different bootstrap draw configuration and has slightly different endpoints.

Missing-score bounds assign failed active controls below every observed score and failed inactive controls above them, then reverse those assignments for the upper bound. The joint 3D all-control ROC interval is **0.3858333–0.6441667**. These are identification bounds, not confidence intervals. The record does not support a claim that 3D scoring improves discrimination; none of these similarities is a calibrated activity probability.

![Labelled control evaluation](controls.png)

Figure: (A) Count-Tanimoto ≥0.5 active-scaffold recall over all eligible held-out actives; bars are target means and points are the two queries. (B) Target-macro ROC AUC on paired complete cases with 95% target-bootstrap intervals; the dashed line marks 0.5. (C) Pairs with at least one failed 3D seed, out of 40 per target. All recorded failures are stereo-policy rejections. [PDF figure](controls.pdf).

## Reproduction and files

### Assay-context sensitivity

The later [assay-context audit](assay_heterogeneity.json) is explicitly post hoc. Of 5,567 records, 1,103 have a representative assay annotated as cell-based and 1,229 have representative descriptions containing “unknown origin” or “unspecified origin.” Across all retained supporting assays, 1,250 records have that origin ambiguity and 79 span multiple assay-format annotations. The text flags do not prove a non-human experiment; they expose unresolved ambiguity despite the structured human assignment. A null variant annotation is not experimental proof of wild-type sequence. ChEMBL's [confidence score concerns assay-to-target assignment](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions), not equivalent experimental conditions.

Requiring query and control to share a retained assay ID leaves six pairs, four with successful 3D scoring, and only one target with both active and inactive observations. That is insufficient for a general activity-discrimination claim. The audit also reports format-matched and single-protein-format sensitivity subsets, their missing scores and actual target counts. Those subsets change the evaluated population; their estimates cannot be treated as improvement over the original cohort or as confirmatory findings. Frozen panels, thresholds and generation settings were not changed in response to these results.

### Recompute the records

From the repository root, use the recorded Python 3.10 environment, project synthesis dependencies and [pinned control requirements](control-requirements.txt). Reproduce the independent metric check, bootstrap and figures directly from the adjacent compact records:

```bash
python benchmarks/results/research-validation-20261016/analyze_controls.py \
  --output benchmark-results/recomputed-controls
```

To reconstruct the panel and recompute chemical scores, extract the [bundled raw ChEMBL snapshot](chembl-source-snapshot.zip), verify its hashes, and use new output directories:

```bash
python -m zipfile -e \
  benchmarks/results/research-validation-20261016/chembl-source-snapshot.zip \
  benchmark-results/reproduced-chembl-sources
python -m benchmarks.research.sources verify benchmark-results/reproduced-chembl-sources
python -m benchmarks.research.import_chembl \
  --manifest benchmark-results/reproduced-chembl-sources/00-chembl-manifest-with-assays/manifest.json \
  --output benchmark-results/reproduced-control-panel \
  --split scaffold --seed 42 --queries-per-target 2 \
  --active-nm 1000 --inactive-nm 10000 --max-heavy-atoms 60
python -m benchmarks.research.controls \
  --dataset benchmark-results/reproduced-control-panel/panel.json \
  --output benchmark-results/recomputed-control-scores \
  --conformers benchmarks/results/research-validation-20261016/common-conformers.json \
  --max-per-label 10 --seeds 42 43 44 --count-tanimoto-threshold 0.5
python benchmarks/results/research-validation-20261016/analyze_controls.py \
  --scores benchmark-results/recomputed-control-scores/control_scores.json \
  --summary benchmark-results/recomputed-control-scores/summary.json \
  --output benchmark-results/recomputed-chemical-analysis
```

The [conformer JSON](common-conformers.json) is exactly `control_summary.json["protocol"]["conformers"]`. Reconstructed panel paths differ from the original local paths, so its byte-level JSON hash can differ despite equivalent chemical records and partitions. [control_summary.json](control_summary.json) records the original dataset hash, local source paths, scoring protocol and full prefilter population; [control_scores.json](control_scores.json) preserves pair scores and seed failures. Original source URLs and per-source hashes are in the source manifests, not in the compact summary. [artifacts.json](artifacts.json) identifies artifact hashes; live source URLs do not guarantee the same bytes. [engineering.json](engineering.json) separately records software checks and environment limitations. The archive covers ChEMBL controls; it does not package all generation-model assets or establish an unseen pretraining inventory.

The [archive relocation check](snapshot_roundtrip.json) verified all 43 packaged artifacts and rebuilt the strict panel in a fresh location. All 5,567 biological records, 12 queries, split definitions and 41 original source URL/hash entries matched. The portable analysis script also reproduced the recorded numerical estimates exactly.

After the source extraction and panel reconstruction above, reproduce the assay audit with the adjacent [audit script](audit_assay_heterogeneity.py):

```bash
python benchmarks/results/research-validation-20261016/audit_assay_heterogeneity.py \
  --panel benchmark-results/reproduced-control-panel/panel.json \
  --manifest benchmark-results/reproduced-chembl-sources/00-chembl-manifest-with-assays/manifest.json \
  --output benchmark-results/recomputed-assay-heterogeneity.json
```

The default control scores are the bundled primary controls; use `--controls` to audit an explicitly recomputed score file. Path and script hashes can differ after relocation; assay counts, retained pair IDs and numerical summaries are the reproducibility targets.
