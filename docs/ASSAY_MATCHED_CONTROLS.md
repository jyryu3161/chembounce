# Same-assay similarity controls

`benchmarks.research.assay_matched` builds a separate retrospective control evaluation. It does not replace the existing generation benchmark or its train/test partitions. The earlier capped 240-pair control sample retained only six pairs after a shared-assay filter; that did not establish that the complete retained source sample lacks suitable controls.

A context is the exact tuple **assay ID, target ID, endpoint, units**. The importer reuses `activity_label`, `assess_direct_human_assay` and `chemical_identity`, verifies every retained source hash, and preserves the original observations. Active/inactive disagreements remove only that family inside that context. Opposite labels for the same family in another assay or endpoint remain separate observations. Duplicate activity IDs with different payloads fail rather than being silently merged.

Structured human assignment, confidence 9 and a null variant annotation are the existing inclusion policy. A null annotation does not prove wild type. Assay descriptions, structured organism metadata and BAO annotations are retained, including unresolved or contradictory information. Prose patterns are not used to infer biological truth or to remove selected experiments. A shared assay ID narrows context but does not establish identical experimental conditions for every observation.

## Audit before selecting or scoring

From the repository root in the recorded Python/RDKit environment:

```bash
python -m benchmarks.research.assay_matched audit \
  --manifest path/to/retained/chembl/manifest-with-assays.json \
  --output benchmark-results/assay-matched-audit \
  --active-nm 1000 --inactive-nm 10000 --max-heavy-atoms 60
```

The existing [ChEMBL snapshot](../benchmarks/results/research-validation-20261016/chembl-source-snapshot.zip) can also be extracted and its portable `00-chembl-manifest-with-assays/manifest.json` supplied. The audit emits `pool.json`, full per-context `availability.json`, `exclusions.json` and artifact hashes. Availability includes all encountered contexts, active/inactive family counts, active scaffolds and the available references for every active query, including contexts with no qualifying records. An audit reports availability in the bounded retained source sample, not in all ChEMBL.

Review availability and fix minima/caps before looking at chemical scores. The first policy specified for this follow-up is:

- Seed 42; at most two contexts per target and one active query per context.
- At least five active and five inactive control families, each different from the query family **and** query scaffold.
- At most ten controls per label. At least one active reference scaffold; there is no additional requirement that all active controls have different scaffolds from one another.
- Never lower these minima to fill a missing target. Record insufficient context/target groups and cap exclusions.

```bash
python -m benchmarks.research.assay_matched select \
  --pool benchmark-results/assay-matched-audit/pool.json \
  --pool-sha256 REVIEWED_POOL_SHA256 \
  --output benchmark-results/assay-matched-selected \
  --seed 42 --min-active-controls 5 --min-inactive-controls 5 \
  --min-active-scaffolds 1 --max-controls-per-label 10 \
  --max-queries-per-context 1 --max-contexts-per-target 2
```

Replace `REVIEWED_POOL_SHA256` with `artifacts.json["pool_sha256"]` from the reviewed audit. Additional assay/context caps are explicit in `panel.json["selection"]`. Selection ranks SHA-256 hashes of the seed, stage namespace and stable identifiers, independent of input ordering. Assays are ordered within target, contexts within assay and query families within context. Query scaffolds do not repeat within a context. Control selection reserves the required number of active scaffolds, then fills the label cap. Neither potency magnitude nor similarity is used as a ranking variable; potency is used only to determine the declared label.

`panel.json` uses the independent `assay-matched-controls-v1` schema. It intentionally cannot be passed to the target-pooled `protocol.load_dataset` or `query_references` APIs. `validate_panel(..., verify_sources=True)` rechecks exact retained activity and assay payloads, labels, identities, representative input states, source hashes, context separation and declared caps.

## Scoring interface and aggregation

`control_pairs.json` contains evaluator-only rows with stable pair ID, query/control SMILES, label, family/scaffold, context, assay, target, endpoint, units and source IDs. `public_queries.json` exposes only query ID, SMILES, target and context; it omits labels and control structures. Context/family repeats are retained rather than presented as new independent compounds.

The curator and selector calculate **no similarity or activity-performance scores**. The separate [scoring CLI](../benchmarks/research/assay_control_scoring.py) validates the selected panel and every retained source payload before calling the existing `controls.similarity_scores` unchanged. It does not select new pairs or tune thresholds:

```bash
python -m benchmarks.research.assay_control_scoring \
  --panel benchmark-results/assay-matched-selected/panel.json \
  --output benchmark-results/assay-matched-scores \
  --conformers benchmarks/results/research-validation-20261016/common-conformers.json \
  --seeds 42 43 44 --workers 3 \
  --bootstrap-seed 20261016 --resamples 10000
```

The referenced conformer file specifies eight conformers and the previously fixed common scoring settings. Other assessments may supply another `ConformerSettings` JSON object, declared before scoring. Requested seeds override the JSON's optional `seed` field. The default worker count is three; workers use separate spawned processes. An existing output path is rejected, including incomplete earlier runs; there is no restart/overwrite mode. Source artifacts referenced by the panel must still be accessible. Recreate the audit/selection from the portable snapshot above when relocating retained sources.

The output preserves the input panel and conformer file, regenerated `control_pairs.json` and `public_queries.json`, and these evaluation records:

- `protocol.json`: exact input, selected-pair, source, scoring-code and effective-settings hashes, runtime versions, and status/aggregation meanings. This records the executed configuration; it does not establish that an external analysis plan was preregistered.
- `scores.json`: every pair ID with its label/context/source fields, per-seed trials, metric aggregates and duration. Atomic snapshots are written after each result in the panel's original pair order.
- `progress.json`: completed pair records, unsuccessful 3D pairs, unexpected errors and stage. `complete: true` means all pairs and both summaries were recorded, even if some or all scoring failed.
- `framework-summary.json`: five metric-specific summaries (bit/count Tanimoto, shape, pharmacophore and joint score). Valid 2D scores remain available after a 3D seed failure.
- `framework-paired-summary.json`: the same metrics restricted to pairs where every requested 3D seed succeeded. Success-conditional estimates compare the same complete cases; missing-score bounds still include **every selected pair**.

Pair `status: success` means the common scorer returned. It does **not** imply that 3D scoring succeeded: inspect `all_3d_seeds_succeeded` and each trial. If any requested 3D seed fails, no 3D aggregate is exported, while valid 2D scores are kept. Unexpected scorer/worker exceptions produce `status: error`; requested seeds without returned results are marked `unreported`, because their execution/outcome is unknown. Such pairs remain in all denominators and bounds. No score is an activity probability, and this control assessment does not replace generation performance evaluation.

An external evaluator can alternatively return one object per `pair_id`:

```json
{"pair_id": "pair-...", "scores": {"score": 0.4}, "status": "success"}
```

The number above is a schema example, not a result. Supply each metric's completed seed aggregate, using the same scoring settings for all pairs. A missing metric or missing pair is retained as unscored. `failed`, `timeout`, `error` and `not_evaluated` statuses suppress stale numeric scores. Duplicate/unknown pair IDs and nonfinite scores are rejected. Keep per-seed trials and scoring-settings provenance alongside score records; the curator does not infer them.

```bash
python -m benchmarks.research.assay_matched summarize \
  --panel benchmark-results/assay-matched-selected/panel.json \
  --scores path/to/independently-computed-pair-scores.json \
  --metric score --bootstrap-seed 20261016 --resamples 10000 \
  --output benchmark-results/assay-matched-summary.json
```

ROC AUC and average precision are computed within each query's context. Summaries average queries within context, contexts within assay, assays within target, then targets equally. Multiple endpoints cannot give an assay extra weight. Every selected pair enters missing-score bounds, including entirely failed query groups. Success-conditional estimates keep unestimable groups as null and report their counts; bounds are not confidence intervals. Target bootstrap resampling never treats assay repeats, compounds or conformer seeds as independent targets. Shared chemistry and related target biology can still correlate across target clusters, especially with only six or fewer targets.

## Recorded availability, without chemical scores

The [2026-10-16 availability record](../benchmarks/results/assay-matched-availability-20261016/availability.json) used all 15,533 retained measurements and yielded 5,890 context-family records. The first policy has the following feasible choices **before caps**:

| Target | Contexts containing both labels | Eligible contexts | Eligible active queries |
| --- | ---: | ---: | ---: |
| EGFR | 24 | 3 | 20 |
| CA2 | 22 | 4 | 28 |
| PIM1 | 16 | 3 | 30 |
| AChE | 48 | 3 | 59 |
| A2a | 23 | 4 | 40 |
| JNK3 | 17 | 5 | 164 |

The fixed hash selection contains 12 queries and 197 query-control pairs across all six targets; caps are maxima, so pairs are not padded to 240. Of the curated target-family groups, 197 recur across contexts and 20 have different labels in different contexts. Those labels are preserved separately. Four within-context conflicting families were excluded. These are availability/selection findings only; no new chemical scores, discrimination estimates or generation-performance claims accompany this record.
