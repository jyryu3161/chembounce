# Qualified recommendations: implementation and evidence gates

The complete objective is to recommend final, route-regenerated molecules that
remain similar to the **original input hit**, change its normalized scaffold,
preserve configured anchors, and have concrete synthetic-route evidence. A
qualification filter alone does not establish useful yield, search quality,
biological activity, or experimental synthesizability.

## Six-stage delivery tracker

| Stage | Required outcome and evidence gate | Status |
| --- | --- | --- |
| 1. Unified final qualification | One versioned conjunction at the real API/CLI export boundary; final-product/original-hit similarity, normalized hop, anchor, exact selected route, requested evidence and complete scoring/seed checks. Adversarial workflow tests, alternative survival, continuation, and matching JSON/TSV/SDF selection. | Engineering gate passed at `e47c38a1734bab4f154fc47b4f4fb7fb5493aa43`; explicit engineering thresholds only. |
| 2. Stereo/design eligibility and bounded refill | Explicit stereo/design eligibility with an auditable, bounded replacement supply when designs fail; preserve caps and complete rejection denominators. | Engineering gate passed at `8461d9a492d181165df53f031f5c12f6a3dde167`. No scientific yield claim. |
| 3. Similarity-aware route/search | Search and route construction retain original-hit similarity while finding changed scaffolds and viable route products; compare yield and compute under fixed budgets. | Engineering gate passed at `6adc2cedd73c14cf5a0ab88fbcd61004f3267ef6`. Controlled toy diagnostics only; scientific yield/compute benchmark remains Stage 6. |
| 4. Calibrated acceptance and hit-prioritized ranking | Calibrate thresholds using development-only controls, prioritize original-hit acceptance/ranking, and ship a usable recommendation preset and CLI default flow through this contract. Previously inspected data cannot be described as held-out evidence. Multi-seed verification and threshold provenance must support the intended profile. | Full gate passed at `f6ab40e0d1d536636015d94c75f1b815898dc907`; reviewed development artifact and both CLI/API paths verified. Explicit opt-in, permissive development thresholds—not a strong-similarity guarantee. |
| 5. Independent synthesis evidence | Establish independent, route-specific evidence and provenance beyond catalog/rule replay; define treatment of unknowns, missing context, atom sources, model disagreement and observation quality. Never present a heuristic as a synthesis probability. | Engineering and known-fixture execution gate passed at `4d47b30c925152198d2fb8e68cdad3b502bacc7d`. Conditional computational support; no experimental or post-selection independence claim. |
| 6. Fresh fair benchmarks and ablations | Fresh artifact directories; comparable proposals, budgets, splits and resources; full input/proposal/route/scoring/seed/qualification failure denominators; fair ablations and confidence/uncertainty reporting. | Ten-attempt inspected development pilot completed; fresh evaluation pending. No scientific performance improvement claimed. |

Each stage requires implementation, focused validation and independent review
before the next stage. Historical published results and source inventories are
not revised to imply these stages were present.

Stage 1 root validation: research runtime **878 passed, 11 skipped**; minimum
runtime **750 passed, 47 skipped, 8 warnings**; Python 3.12 **868 passed,
20 skipped**. Ruff lint/format and 20 portable replay tests passed. The minimum
RDKit runtime skips modern 3D alignment checks, so its pass count is not evidence
that it can execute the qualified 3D path. Independent review passed the
reachable contract; synthetic post-refiner malformed-report persistence remains
a nonblocking robustness limitation. Root validation record:
`qualified-development-20261016/stage1-e47c38a-root/validation.json`, SHA-256
`a357511436d0b3c6005636caf252855cbcb89c68747f0155a8dec4dc53aca122`.

Stage 2 root validation: parent `bc12d92` full research suite **901 passed,
11 skipped**; focused minimum runtime **63 passed, 22 skipped, 8 warnings**;
focused Python 3.12 **146 passed**. The minimum-runtime 3D limitation above
still applies. The final narrow unique-route-counter correction passed **45**
root tests and **17** unchanged independent tests. Root gate record:
`qualified-development-20261016/stage2-gate.json`, SHA-256
`2b8be505ad9ba2cb5f67714d2520e9e52ba41db7352f50c9a51054b98c613df7`.

Stage 3 root validation: full research suite **917 passed, 11 skipped**;
focused Python 3.12 **142 passed**; final minimum runtime **59 passed,
22 skipped**; independent review **14 passed**. Modern 3D integration tests
remain skipped on the minimum RDKit runtime. Gate record:
`qualified-development-20261016/stage3-gate.json`, SHA-256
`827ad5b9f5adcbb3d6b33fe9fbad7a2d2ae12f0dccb428b44c9b6843de9881ac`.

Stage 4a scorer/ranking gate passed after the SDF provenance correction at
`7394d30c90ba1d54ffc5933cb94270cf9d9ed18a`. Root parent validation:
**933 passed, 11 skipped**; focused Python 3.12 **117 passed**; minimum
runtime **25 passed, 17 skipped**. The correction passed **16** objective/export
tests and **19** unchanged independent tests. Gate record SHA-256:
`792edba1c666d61a8f4bc7c6173325ea3ef980ed03d3cd39eab1acdf33ad3365`.

Stage 4 full root gate at `f6ab40e`: **1028 passed, 11 skipped**; focused
minimum runtime **35 passed, 7 skipped**; Python 3.12 **102 passed, 1 skipped**;
Ruff lint/format passed for 140 files. Independent release review: **37 passed**
on research, **30 passed, 7 skipped** on minimum, **35 passed, 2 skipped** on
Python 3.12. Calibration, reporting-only sensitivity and release artifact chains
were independently reviewed. Root gate record:
`qualified-development-20261016/stage4-gate.json`, SHA-256
`9b419280f4b19c996a0cff0290cf2eecaf5513b346376320832cb50daaf131a2`.

## Stage 1 interface

The canonical API remains `chembounce_synthesis.workflow.run(hit, designs,
config, output_dir)`. Set:

```json
"recommendation": {"mode": "qualified", "policy_version": "qualified-v1"}
```

Also explicitly supply `minimum_scores.hit_tanimoto`, `hit_shape`, and
`hit_pharmacophore`, each finite and in `(0, 1]`, using your declared engineering
policy. There are deliberately no invented empirical defaults. Other configured
minimum scores remain mandatory. `hit_tanimoto` is binary Morgan Tanimoto,
radius 2, 2048 bits, with chirality; it is **not** the sparse-count Morgan value
used by the hopping CLI's `-t` option. `hit_morgan` is the existing alias.

Both CLIs accept `--recommendation-mode qualified` with a synthesis configuration
containing those minima:

```sh
chembounce-synthesis --hit "$HIT_SMILES" --design "$DESIGN_SMILES" \
  --config engineering-policy.json --recommendation-mode qualified --output fresh-run
chembounce -i "$HIT_SMILES" -o fresh-integrated \
  --synthesis-config engineering-policy.json --recommendation-mode qualified
```

Qualified mode requires hop mode, strict specified stereochemistry, 3D scoring,
and `all_seeds` refinement thresholds. It cannot be combined with
`--no-synllama`. A single screening seed qualifies only that requested protocol;
it does not establish stability. If refinement is configured, every requested
seed at the final ensemble size must succeed and pass all its pose minima. The
exported representative pose, score components, identity and seed must agree.
Configured convergence tolerances must also be assessed and satisfied. Budgets
and deadlines cannot convert incomplete evaluations into recommendations.

Catalog membership and reaction-rule replay are always required for the exact
final product. `synthesis_evidence.policy` controls additional existing
requirements. Missing or unknown **required** evidence fails; unrequested
independent evidence can remain unknown and is reported as such. Passing the
baseline policy does not claim independent or experimental synthesis validation.
The named development profile additionally requires the Stage 5 computational-support contract below; this legacy engineering policy is labelled `legacy-conditional-v1`.

## Artifacts and compatibility

`manifest.json` and `qualification.json` save the effective versioned policy.
Every route attempt receives a qualification decision or an explicit reason it
was not qualified/assessed; `routes.json`, `scores.json`, `refinement.json`,
`proposals.json` and `summary.json` retain earlier failures and unassessed work.
`qualification.json.selected_route_ids` and `results.json.results` identify the
same final selection as `results.tsv` and the pose SDFs. Rows and poses explicitly
state `qualified` or `exploratory_verified`. Screening tables are intermediate
diagnostics and are not the final recommendation list.

Qualification occurs before family representative selection and the result cap,
and inside refined design fallback. Failed route/design alternatives are not
labelled as qualified alternate routes. A fresh run with no qualified products
has `status=no_candidates`, an explanatory reason, empty results and header-only
TSV. The CLI returns 3 for an empty result and 2 for configuration/runtime errors.

Omitting `recommendation` or specifying `mode=exploratory` preserves historical
search/scoring settings and outputs explicitly unqualified exploratory records.
Continuation carries routes and saved proposals, then replays verification and
recomputes scoring/qualification under the new effective policy; historical
scores or qualification claims are never authorization to pass a new run.

## Stage 2 bounded design supply

Qualified runs use selection policy `bounded-refill-v1`. `max_designs` is the
eligible unique-design quota; `max_design_attempts` bounds raw intake and defaults
to four times that quota. This multiplier is an engineering resource limit,
not a calibrated scientific setting. Invalid SMILES, exact-state duplicates,
and stereo/reference failures each consume a raw slot. Intake stops at the
eligible quota, raw bound, exhausted input, or deadline. It never enumerates or
assigns unspecified stereoisomers, and distinct specified enantiomers remain
distinct designs. Provider sample seeds and per-design proposal order are
unchanged; rejection adds no model, scoring, or refinement retries.

Exploratory configurations without `max_design_attempts` retain
`legacy-capped-v1` behavior. Explicitly setting that field opts exploration into
bounded refill. Both CLIs expose `--synthesis-max-designs` and
`--synthesis-max-design-attempts`. The integrated CLI supplies at most the raw
attempt bound from hopping; an explicit `--overall_max_n` remains a hard upper
bound and can prevent refill. `run_manifest.json.synthesis_design_supply` records
the requested and effective upstream limits. Product scoring, route search,
refinement, and final result limits remain separate.

`manifest.json.design_preflight` records each consumed input, exact canonical
identity where valid, reason codes, duplicate origin, eligibility, counts and
all simultaneous limits reached. Sized inputs have exact supplied/unattempted
counts without reading the tail. Streaming API/TSV inputs are consumed only
within the bound; their unread totals remain explicitly unknown. There is no
extra read merely to establish exhaustion. Query stereo failures retain the
original input and a structured `query_preflight` reason with qualification
unassessed; the existing API error remains an error, not evidence that no useful
chemistry exists. Runs with no eligible designs load no synthesis assets.

Continuation preserves prior preflight lineage and all carried graphs,
including designs currently rejected or unassessed. Already saved raw prefixes
and known carried design identities remain available to a later larger budget;
unread tails cannot be reconstructed. Replayed inputs consume the new attempt
budget, and carried routes do not count as newly scored work unless selected
and actually scored. Saved scores never authorize eligibility or qualification.
When a bounded continuation has a carried graph but its latest saved-proposal
archive has no response for that design, the run records a carry-only design:
it rechecks the graph without fabricating a proposal or calling the provider.
Other saved responses are replayed normally; malformed provider data still fails.

Historical benchmark adapters and frozen studies remain unchanged. Their
upstream candidate caps and query-failure conventions must be addressed in the
new Stage 6 benchmark profile before attributing any benchmark result to refill.

## Stage 3 original-hit route search

Qualified runs default to `route_search_policy: original-hit-v1`; exploratory
runs retain `design-v1`. Both CLIs accept `--route-search-policy` for an explicit
override, including `design-v1` as a qualified-mode search ablation. The final
qualified-v1 acceptance contract and configured thresholds remain unchanged.

Projection and expansion assess each complete regenerated root product before
streaming and final beam/frontier truncation. Cheap eligibility reuses the
configured hop, core, property and 2D checks, with strict specified stereo in
qualified mode. Eligible roots rank first by original-hit binary Morgan, then
actual workflow-design similarity and a deterministic route key. Distinct exact
products precede alternate routes within each eligibility tier; spare slots
retain route alternatives. Intermediate reactions and role-based catalog leaf
shortlists keep their own targets, without whole-hit thresholds on fragments.

When expansion is enabled, ineligible roots may fill remaining existing beam
slots as explicit seed-only states. Earlier expansion rounds may retain them
within the same frontier capacity; they do not consume cumulative eligible
output slots. The final round discards them. They never receive evidence/3D
scoring or recommendation status. This permits a bounded sequence of catalog
substitutions to recover a qualifying hop. Combination limits, frontier widths,
round limits and model calls are unchanged. Actual substitutions can increase
relative to legacy early termination, while remaining within those limits.
Existing workflow deadline boundaries remain in force; a single route call is
bounded by its search limits rather than interrupted midway by the deadline.

`manifest.json.route_search_policy` records the effective policy.
`search_objectives.json` records each design's exact hit/design references,
2D minima, LRU cache counts and bounded decision sample. The positive integer
`route_search_cache_size` defaults to 4096 entries per design; the first
`min(size, 256)` distinct product decisions are sampled separately. Evaluation
counts are cache misses, including repeat evaluations after eviction, not an
unbounded global unique-product count. Search reports count generated root
states, repeated pruning visits, ineligible/filtered visits, capacity truncation,
retained seed visits and per-round eligible/frontier counts. Samples are not
complete failure inventories; aggregate counters retain all measured events.
The cache and seed-only distinction are inactive under `design-v1`.

The beam-one amide fixture tests original-hit guidance against legacy search
with equal actual combination counts. Its explicit 0.5 threshold is an
engineering diagnostic, not a production cutoff. Separate two-round seed
recovery tests use equal configured caps. These tests establish bounded search
behavior, not global recall, general scientific improvement or calibrated yield.

## Stage 4a scorer and ranking interface

`conformers.scoring_objective` defaults to historical `joint-v1`.
`original-hit-v1` selects each candidate conformation using only its original-hit
overlay: `sqrt(hit_shape * hit_pharmacophore)` multiplied by separate exponential
penalties for the selected hit-reference and candidate relative energies. It
retains the legacy strict first-on-tie choice and equals the historical scorer
with the hit supplied as both references. The actual design is overlaid only
after selecting the candidate conformation, with the same internal geometry.
Both real-design ensemble generation and the secondary overlay must still
succeed; otherwise existing scoring/refinement failure counts include that
candidate. This availability requirement is not a hit-similarity failure.

The explicit `ranking_policy` defaults to the selected scoring objective.
`original-hit-v1` ranking requires the corresponding scorer and uses only
original-hit bitMorgan and minimum-over-final-seed hit shape/pharmacophore for
Pareto dominance. Within a front, minimum seed hit-objective scalar precedes
observed score spread, design similarity, route steps and deterministic identity
ties. Screening uses its one assessed seed. Robust ranking metrics remain
separate from the exported representative seed's pose metrics. The same policy
orders alternative designs before family selection. Design/route advantages
cannot rescue a hit-dominated candidate. Dispersion is not a convergence or
confidence claim.

Both CLIs expose `--scoring-objective` and `--ranking-policy`; API configurations
use the fields above. Saved settings and scoring results record the objective;
strict validation checks seed identity, realized ensemble counts, finite pose
metrics and consistency of both energy penalties. Refinement preserves the
objective across every requested seed and final-ensemble stage.

This interface alone is **not Stage 4 completion**. The reviewed development
profile below adds fitted thresholds and provenance. Its protocol uses fixed eight
conformers with seeds 42/43/44, matching eight-conformer seed-42 screening, no
adaptive enlargement and no arbitrary range tolerance. Previously inspected
controls are development data only. The scorer change itself introduces no
empirical cutoffs.

## Stage 4c recommendation profile (gate passed)

The API field `recommendation_profile: "qualified-development-v1"` and both
CLIs' `--recommendation-profile qualified-development-v1` select one named
development profile through **explicit opt-in**. Catalog, rules and
proposal-provider configuration remain explicit; no manual cutoffs are needed. The current recommendation flow also requires an explicit pinned local forward-model bundle for Stage 5 support.
The default exploratory configuration and explicit engineering-qualified mode
remain available. Missing or corrupt profile assets fail clearly rather than
substituting values or falling back.

```sh
chembounce-synthesis --hit "$HIT_SMILES" --design "$DESIGN_SMILES" \
  --config synthesis-assets.json --recommendation-profile qualified-development-v1 \
  --forward-model-path "$REACTIONT5_LOCAL_BUNDLE" --output fresh-development-run
chembounce -i "$HIT_SMILES" -o fresh-development-integrated \
  --synthesis-config synthesis-assets.json \
  --recommendation-profile qualified-development-v1 \
  --forward-model-path "$REACTIONT5_LOCAL_BUNDLE"
```

The [packaged reviewed profile](../chembounce_synthesis/data/qualified-development-v1.json)
sets these conjunctive minima, with both 3D criteria required for **every** final
seed:

| Original-hit metric | Minimum |
| --- | --- |
| Binary Morgan Tanimoto | 0.08695652173913043 |
| Shape overlap | 0.5055221213607054 |
| Pharmacophore overlap | 0.11206335132501284 |

These are **permissive development thresholds**, especially the 2D and
pharmacophore minima. Passing does not establish strong structural similarity,
biological activity retention, or superiority of 3D scoring. CLI output states
this scope explicitly. The fit is preserved unchanged; raising a threshold is
an explicit custom policy, with no inherited calibration-retention claim.

The frozen panel contains **197 pairs** (111 labelled active, 86 inactive),
12 queries and **six previously inspected targets**. All three computationally
failed pairs remain in performance denominators. Results are target-macro rates:

| Development diagnostic | Active retention | Topology-changing active retention | Inactive rejection |
| --- | --- | --- | --- |
| Final fit on all development targets | 85.00% | 84.17% | 24.44% |
| Locked leave-one-target-out rules | 80.83% | 80.00% | 32.78% |

The six leave-one-target-out fits were feasible, but active retention on the
excluded target was only **50% for CHEMBL2147** and **55% for CHEMBL220**. These
retrospective development diagnostics are not unseen validation; operational
inactive rejection includes computational failures and is distinct from purely
geometric discrimination. No full-pipeline 80% retention promise follows.

Profile defaults apply to raw configuration **before** generic protocol defaults
and validation. CLI overrides are applied before this resolution, so an existing
`research-v1` declaration cannot inject adaptive stages or arbitrary range gates
into the fixed protocol. Explicit values survive and are audited. API callers
should specify the profile in the raw configuration or before `load_config`;
settings already materialized by a prior load are effective overrides.

The complete profile supplies eight-conformer seed-42 screening, final eight
conformers at seeds 42/43/44, all-seed acceptance, minimum-seed hit ranking,
original-hit route search and strict hop/stereo handling. Engineering limits
default to 20 eligible designs, 80 raw attempts, 256 scored products, 10 refined
families, 240 final candidate-conformer requests and 10 final results. Explicit
resource limits remain authoritative and are recorded separately from numerical
protocol changes. Additional anchors, property bounds and evidence requirements
remain explicit and may reduce coverage.

The release artifact has its own reviewed-profile schema/status, binds the
untouched fit candidate and root/independent review record digests, and is pinned
by its complete-file SHA-256. Its numerical values, procedure and provenance
must match the reviewed candidate; packaging never fits or adjusts cutoffs.
Package data and wheel-import tests cover the JSON resource.

Release SHA-256:
`9d3c60a895912b28c13a29aac0a04f385ddb8a0effc9fcd81d64a50a3d48bdd2`.
The release binds candidate
`56cc18999701745b89eca83a5798d3a5f67d1a2789e80d8b6c87aff3b6e7fad7`,
fit `2294a4d86c38e2e3c5fede92ab5095f6f563693ee1fd76f04c7a19ee74aae77d`,
root fit gate `264cfaaae02039680f7f12931eef6c9b2e4c8b948c1ac9cc27a5c4a7ad46e39b`
and independent fit audit
`532faeb3c6e14929633f56063cdd579d7ed9dda7243b3d8835d1ceb0284f7a67`.

`development_protocol_matched` requires qualified mode, the complete fixed
three-seed acceptance settings, matching minima, original-hit scoring/ranking,
and the recorded numerical runtime. Material overrides—including stricter
thresholds, seeds, ensembles, energy settings, bound references, pocket
conditioning, exploratory mode or adaptive/range gates—are labelled
`custom_unvalidated` (or rejected when internally inconsistent). Exact Python,
RDKit, NumPy and SciPy versions, `BaseFeatures.fdef` SHA-256, and the frozen
`conformers.py`, `common.py`, `identity.py` kernel hashes are compared. Full
calibration and current source/runtime snapshots remain auditable; reviewed
orchestration/packaging changes are distinct from numerical-kernel drift.

Manifest and qualification policy record fit/procedure/source/runtime/review
digests, effective overrides and scope. JSON, TSV, hit/design pose SDFs and
reference SDFs carry profile identity, match status and artifact digest.
Continuation retains prior provenance but re-resolves the current artifact,
settings, runtime and every route/seed; historic profile claims never authorize
a new result.

A matched score protocol is **not a pipeline retention promise**. These controls
are previously inspected development data, and real-design diagnostic
availability, anchors/properties, finite search/refinement budgets and synthesis
evidence introduce additional vetoes. No biological activity, future-seed
robustness, experimental synthesis or held-out validation claim follows.

## Stage 5 computational route support (gate passed)

The named `qualified-development-v1` recommendation flow now requires
`synthesis_evidence.profile="computational-route-support-v1"`. The API supplies
`synthesis_evidence.predictor.model_path`; both CLIs accept
`--forward-model-path`. Missing configuration/assets or an attempted evidence
profile downgrade fails clearly. To reproduce historical exploratory or explicit
engineering-qualified runs, omit the named recommendation profile and use the
separately labelled `legacy-conditional-v1` evidence policy. Supplied observations
or attribution prose in that legacy policy do not become verified evidence.

The stronger profile adds selection constraints to the unchanged Stage 4
numerical artifact and kernels. Its success rate does **not** inherit the
calibration panel's retention figures. It requires every step of the exact
selected route to have fresh catalog/rule replay, concrete atom/stereo audit,
and at least one completed exact full-product hypothesis among three pinned
ReactionT5 beams. A different model from SynLlama is useful selection support,
but training independence is unknown, and using ReactionT5 during selection is
not independent post-selection validation. No output is an experimental success
claim or a synthesis probability.

`route_audit.py` binds route topology, ordered reactant slots, parent/intermediate
bindings, catalog/rule content digests, rule SMARTS and exact regenerated products.
All retained heavy atoms and explicit isotopic hydrogens require concrete
reactant-slot/atom origins; a template-created atom without a donor fails even
when free-text attribution is supplied. Ordinary hydrogen/proton transfer,
byproducts and reaction conditions are not mass-balanced. All observed matching
replay witnesses must have consistent defined atom/stereo correspondence.
Symmetry and bounded enumeration/truncation are explicit; this is not a claim
of exhaustive reaction outcomes. Fully specified new/inverted/template-asserted
stereo additionally requires exact forward-product agreement and is not a
selectivity guarantee. Unspecified or unsupported route-state stereo fails.
Older RDKit without concrete `react_idx` origin metadata fails closed.

Only sorted canonical molecular reactant/reagent **multisets**, preserving
duplicates, enter the predictor's documented `REACTANT:...REAGENT:...` payload.
Ordered route slots remain separately recorded. Expected product, rule/template,
original hit, design, source prose and physical conditions never enter the model.
Supplied reagents must parse as molecular SMILES. Temperature/time remain
unvalidated physical context; this model does not consume them. Verified physical
context and empirical-success policies are unavailable in this profile, and
supplied success records cannot enable them.

The strict adapter pins the complete numerical checkpoint/tokenizer/config
bundle and declared revision `933114058cb2604dc1bf536dbebdfcefbe83d4fc`, including
weights SHA-256 `69c23ca364e63e0e2be22d114652aad9f3df4e161d43936eee0dc3ce74a87fb7`.
It requires Python 3.10.12, Torch 2.8.0, Transformers 4.51.3, tokenizers 0.21.4,
safetensors 0.6.2 and RDKit 2025.9.6. It uses local files only, no remote code,
CPU with one Torch thread, deterministic top-three decoding, at most 256 new
tokens per beam and 512 input tokens. Bundle, adapter, runtime, input-token,
thread and decoding provenance is retained. Standard local Hugging Face download
metadata is separately recorded; it cannot replace any pinned numerical file.
The adapter also pins the Transformers generation implementation, SHA-256
`dc9044d10c00abe7d341f2e777fb623925ff59dda2c1d814ecc7d6a685ad2ee1`.
Default global evidence limits are
128 calls, 32768 input tokens and 32768 output tokens; explicit lower limits and
the workflow deadline remain authoritative. A failed call charges its reserved
output allowance and is not retried.

The normalization policy was frozen before model smoke execution: remove ASCII
tokenization spaces, canonicalize isomeric molecular SMILES and remove atom-map
labels. Reject unknown/special content, other whitespace, empty components,
names/CXSMILES and terminal dots; preserve every molecular component. The local
model card documents trailing-dot cleanup, so this deliberately stricter policy
may lose coverage. It is not a claim that upstream preprocessing is a bug.
Original decoded text, full/generated/content token IDs, actual `beam_indices`,
EOS completion and every invalid or unfinished beam remain visible. Padding
policy `transformers-4.51.3-beam-ancestry-v1` requires a nonempty contiguous
prefix of beam indices 0–2, followed only by -1. Rows must align with the returned
sequences, excluding the decoder start. Only the witnessed prefix contains
generated tokens; its unique terminal EOS establishes completion. The -1 suffix
must contain the exact EOS1 fill used by this pinned runtime when PAD0 is falsey.
EOS in fill cannot complete an unfinished prefix. Generated PAD/unknown/special
content is never erased into a valid molecule. Every witnessed generated position
is charged; malformed metadata fails closed, preserves the returned tensors and
charges the reserved allowance. Missing ancestry has no permissive fallback. One valid exact full-product beam can support
a step; absent a match, invalid/unfinished beams mean partial evaluation, not
complete model disagreement. There is no component-only match or punctuation
repair.

The cache stores actual molecular inputs, model/settings bindings and raw
predictions, not expected-product decisions. Each route comparison is recomputed.
Final qualification performs another fresh route audit and re-decodes/rechecks
the input-bound cached token evidence without another generation call. Supplied
matched flags, hashes, attribution claims or a different route's success cannot
substitute for this check. Alternative routes qualify separately; continuation
replays and predicts under the current policy rather than trusting prior reports.
JSON/TSV and selected pose/reference SDFs identify the effective evidence profile
and selected route-support binding, with complete evidence in `routes.json`.

The initial source-review failure for undercounted unfinished PAD tokens was
corrected before execution. The metadata-only load failure and first two-call
adapter exercise are preserved. The latter returned zero recommendations under the old decoder and
exposed its incorrect assumption that unused suffixes contain PAD0. That run
discarded beam ancestry and cannot be retroactively certified. The correction
changes interpretation of witnessed tensor padding, not SMILES normalization,
chemical thresholds or model settings. A separately frozen correction regression
then repeated only the same two known inputs with actual returned ancestry.

Stage 5 passed at `4d47b30`: root research **1113 passed, 11 skipped**, focused
minimum **102 passed, 36 skipped**, Python 3.12 **181 passed, 1 skipped**.
Independent source tests passed **65** on research and **65** on Python 3.12;
minimum **41 passed, 24 skipped**. Minimum-runtime skips do not establish support
for the modern concrete-origin or qualified 3D path. Independent actual-artifact
review passed **321 checks**. Root gate:
`qualified-development-20261016/stage5-gate.json`, SHA-256
`3e980bca062871f177e253c80cf3abe06474d9e7f475b1275f26f3889f5aba8b`;
independent actual audit SHA-256
`2653824c74f26afe40a9258f1fc2fd15bf6490f52c1bfac383b48689e429c244`.

The corrected CPU regression used **two new calls**, 75 input tokens and 120
output tokens in 16.64 seconds; cumulative actual calls across the preserved
runs were **four**. The known model-card example matched at rank 1. The existing
saved-proposal amide fixture produced **one qualified recommendation** through
the actual named-profile workflow, eight conformers at each of seeds 42/43/44,
with matching JSON/TSV and all four SDF bindings. This benzene-to-pyridine change
is an atom-aware scaffold hop with unchanged generic topology. It is a known
fixture using a saved proposal, not live end-to-end SynLlama generation, fresh
benchmark performance, experimental synthesis or activity validation. The earlier
zero-recommendation run remains unchanged. Fresh comparisons and ablations remain
Stage 6, whose execution requires a separate gate.

## Stage 6 source path (inspected pilot completed; fresh evaluation pending)

The new `benchmarks.research.qualified_runner` consumes an immutable
`qualified-benchmark-protocol-v1` JSON with public `queries` (`id`, `target_id`,
`smiles` only), generation `seeds`, ordered `methods`, pinned file/directory
`assets`, complete `source_sha256`, explicit stage budgets, and
`collector_catalog_policy: "runner-preload-catalog-v1"`. This policy requires a
new reviewed protocol freeze; older pilot protocols are not silently upgraded.
Each method has
an `arm` (`full`, `baseline`, `crem`, `mol2mol`, or `synllama`), common
`synthesis_config`/catalog/rules assets, and settings. External native methods
also declare their exact command, native asset aliases, settings and
`request_kind`; the frozen baseline receives `pipeline`, while the external
comparators receive `external`. Commands explicitly consume `{request}` and
produce `{response}`. This is a new path; historical adapters/results are not
relabelled or changed.

```sh
python -m benchmarks.research.qualified_runner \
  --protocol frozen-qualified-protocol.json --output fresh-qualified-run
```

The complete query/seed/method grid is saved before work. Method order rotates
cyclically by query/seed block, with the full order bound into the plan digest.
Missing, failed, unsupported, interrupted and timed-out attempts remain in the
planned denominator. No automatic retry, partial native rescue, cross-method
prediction cache or fallback model is available. Source tests use synthetic
inputs/controlled beams only; this command requires a separate execution gate.

The prospective engineering pilot uses native **360 seconds**, common **600
seconds**, outer **960 seconds**, one sequential attempt and a 24000 MB
limit on each sampled child process tree. This limit does not enforce combined
parent-plus-child host memory.
The earlier planning-only 300/600/900 allocation was revised before any pilot
inference because unchanged native Mol2Mol validation requires a literal
360-second budget. Its validator and sampling protocol are preserved. Native
and common subprocess groups use the existing bounded runner, including cold
imports, asset checks and model loading. Unused native time cannot enlarge the
common 600-second allocation. The outer clock starts before native setup and
includes both cleanup phases, collector and accounting; late outputs are
rejected. Parent collector/hash/cleanup checks are cooperative and may finish
after that clock, so this is not an exact wall-clock kill guarantee for parent
work. Shared administrative setup measures preflight asset validation plus one
full eager catalog construction per unique pinned path/checksum. The runner owns
these read-only catalog snapshots for this campaign only. Each method must pin
`catalog_manifest` to the adjacent `<catalog>.manifest.json` path and SHA256, or
set `catalog_manifest_absent: true` when that sidecar is absent. Presence, path,
TSV/sidecar hashes, source and runtime are checked before reuse; per-route catalog
digest checks remain. Earlier protocols that pinned only the TSV have a narrower
provenance scope and are not retroactively relabelled.

Only collector record lookup is shared. Every nonempty collector call constructs fresh
Rules, and native/common children remain independent cold invocations, including
their own catalog costs. No reaction-result, fingerprint-search, prediction,
score, route or artifact cache is shared. Setup failure or cancellation starts no
attempts and preserves the full planned denominator. Cancellation during eager
construction is cooperative and is checked after that validation completes.

`setup.json`, collection and analysis record setup wall/parent CPU time, parent
RSS point observations and its lifetime high-water mark. Setup cost is reported
separately and allocated equally over **all planned attempts** for accounting;
it is not free compute and cannot enlarge any attempt budget. Each `phases.json`
likewise records parent CPU/RSS separately from sampled child-tree resources.
These distinct RSS maxima are never summed into a purported simultaneous peak;
combined parent-plus-child peak memory remains unknown. Total campaign wall time
through checkpoint includes setup and attempts.

Current generation supplies at most twelve raw designs for three eligible
inputs. CReM and Mol2Mol retain their complete fresh, in-budget raw pools (at
most 256) and select twelve by their own count-Morgan or NLL ordering. The
baseline preserves old native algorithms and diagnostics, then passes **only
its selected primary three routes/products** through an empty SavedProposals
carry bridge with zero expansion. No unselected baseline route, new proposal or
current reconstruction can rescue it. Direct SynLlama receives the original
hit with six proposal samples; other common generative arms use two samples
per eligible design, at most six. Common screening/refinement remain eight
conformers at fixed scoring seeds 42/43/44 regardless of generation seed.
Common budgets are twelve screens, six refinement candidates/144 requested
candidate conformers and 24 strict forward calls/12288 input/18432 reserved
output tokens. The old baseline retains its separate 504-conformer native
refinement allowance; these are bounded system comparisons, not equal actual
work or isolated causal effects.

The declared full-method ablations change only raw-attempt capacity (three
instead of twelve, retaining the same supplied pool) or route-search policy
(`design-v1`). The final numerical/evidence contract is unchanged; the latter
is disclosed as a custom search configuration. Joint-objective ablation is
not implemented because its screening/final-rescoring protocol needs a separate
prospective design.

The collector replays the selected DAG/atom/stereo audit, validates saved raw
forward tensors with the pinned tokenizer only, checks source/request/config
bindings, every final seed and JSON/TSV/SDF selection consistency, and computes
bit-Morgan fingerprints and per-metric all-seed minima itself. It never loads
model weights, calls a predictor, reconstructs an evaluator cache or promotes
supplied accepted flags into live evidence. It attests integrity of trusted,
pinned workflow artifacts; arbitrary invented token tensors cannot prove a
model call. Robust seed-minimum metrics are reported separately from exported
representative-pose scores.

The analysis primary endpoint is at least one jointly qualified output per
planned attempt, averaged over generation seeds within queries and then equally
over targets. Secondary unique-family yield is divided by three, including
unfilled slots. Structural diversity with insufficient outputs is unavailable,
not zero. Paired target-cluster bootstrap uses shared draws; single-target
intervals are unavailable. Sampled process-tree CPU is a **lower bound** and
RSS is the sampled aggregate peak, not exact CPU time or memory. GPU time and
memory remain null when unmeasured. Native/common work units and partial-file
accounting failures remain explicit; missing work counts are not invented as
zero. All retained costs include failed attempts. Parent CPU/RSS observations
include collector work and retained shared catalogs, with the separate scope
described above; subprocess metrics remain unchanged. The bounded runner and
collector passed source review at `1841ed8`; the two-input, ten-attempt inspected
development pilot used that immutable source. All ten planned attempts completed
with at least one qualified output, yielding 27 output occurrences across methods
and queries; this is not a count of globally distinct compounds. The final
collection SHA-256 is
`9ccee0460f2c9ad38661c43bda8ad207335428f8c00812e4781ac3a94c189caf`.
The final ledger passed 926 independent checks, and each of the ten selected-output
detail audits passed (3,059 checks in total). The root result gate is
`qualified-development-20261016/stage6-inspected-pilot-result-gate-1841ed8-v3.json`,
SHA-256 `647b4859859dece15a5f529dfa745691e2d902215cdd8e68446f26852102f762`.
These results attest engineering behavior under the existing
permissive development profile, not high similarity or scientific superiority.

| Method | CHEMBL203 output count | CHEMBL205 output count | Mean attempt wall seconds |
| --- | --- | --- | --- |
| Full method | 3 | 3 | 386.09 |
| Frozen develop, common assessment | 2 | 1 | 325.24 |
| CReM, common backend | 3 | 3 | 368.09 |
| Mol2Mol, common backend | 3 | 3 | 382.42 |
| Direct SynLlama, common assessment | 3 | 3 | 368.62 |

Both inputs were previously inspected, and only generation seed 42 was used.
The full method, CReM and direct SynLlama returned the same three products for
CHEMBL203. Thus this pilot establishes neither a unique-output advantage nor a
speed advantage for the full method. Some selected products lose input functional
groups: preservation is guaranteed only for explicitly configured anchors. The
full method's CHEMBL205 candidates span all-seed-minimum pharmacophore scores
0.191–0.314, illustrating why the current acceptance profile cannot support a
strong similarity claim. Later source changes do not modify this completed pilot.
Revised protocol activation, fresh input selection, frozen final grids, execution
and scientific interpretation remain pending.

The first shared-catalog diagnostic reproduced all 18 ordered qualified records
from seven saved attempts exactly, but its final guard assertion failed with one
blocked operation. That failed report remains unchanged. A separate import-only
probe reproduced one blocked `socket.bind` from urllib3's IPv6 loopback capability
check while importing Transformers. The operation was denied before binding, and
the library caught the guard exception. This is consistent with the old counter
failure; the old run did not retain a call stack.

The corrected diagnostic retained every operation block and recorded each event
before raising the guard exception. It passed with exactly one expected, pinned
IPv6 import probe blocked and no unexpected blocked operations. One shared
read-only catalog, seven fresh rule sets/replay decoders, and 18 fresh route/raw
prediction checks reproduced all 18 ordered outputs exactly; all 1,443 pinned
input hashes were unchanged. Independent actual-result review passed 1,511 checks.
The root result gate is
`qualified-development-20261016/stage6-catalog-equivalence-result-gate-v2.json`,
SHA-256 `6febbafb4f257f29d814fb5719f37c408007c542fbec713dd90d663f0f4143c1`.
This establishes saved collector-output equivalence only. Its measured setup
and collector costs do not establish an end-to-end speedup or retime the pilot.

Fresh-panel preparation is implemented separately in
`benchmarks.research.qualified_fresh_panel`; actual acquisition and selection
remain gated. Its `acquire --protocol ... --output ...` command is the only
network path. `select --protocol ... --snapshot ... --output ...` verifies and
replays retained bytes offline. A `qualified-fresh-panel-protocol-v1` binds the
literal selection policy, source hashes, Python/RDKit versions, exclusion
manifest and official target-mapping source hashes, plus explicit request,
response-byte, total-byte, timeout and wall limits. There are no automatic
retries or changed-release fallbacks.

Independent review of `7b073f2` found that a same-origin redirect could change
the frozen activity query. Correction `f4ca7cd` rejects redirects before any
follow-up request and binds the endpoint and full query multimap in both live
acquisition and offline replay. Harmless parameter-order/encoding differences
are allowed; extra filters, duplicate/blank changes and fragments reject. Root
focused tests passed **73** in each of research, minimum and Python 3.12;
independent review passed **40**, including the original 29 regression
assertions. These were synthetic/local tests, not fresh-panel acquisition.

Catalog-only setup commit `d7b24f5` passed **94** root research tests and **75
passed, 19 skipped** in each alternate runtime, plus **42** independent research
checks. Its integration changes neither the live pilot nor child workflows.
Retain `setup.json` or runner-produced `analysis.json` with cost reporting:
regenerating the standalone pure-analysis report does not propagate optional
administrative setup metadata. Saved-output equivalence has passed as described
above; exact future protocol activation remains pending.

The finalized pilot's retained structured artifacts have also been added to the
future-panel exclusion inventory. The reviewed union contains 17,563 exact
chemical states, 17,021 normalized families and 7,240 nonempty normalized
scaffolds, including 3,400 newly added exact states. Its immutable manifest has
SHA-256 `1d9e5b99a13a217fb0ec004abf66f1cf7d7d3e9a115d494666ceb98eb4fd817c`;
the later root gate is
`qualified-development-20261016/stage6-pilot-augmented-inventory-root-gate-v1.json`,
SHA-256 `b63ad7f86918d299186c361ce75a25fd8a4c4e7aa5520514e769c36474ee8053`.
The inventory covers the declared retained artifacts, not every historical or
model-training structure. Forty nonstructured artifacts were not newly parsed;
unreturned/internal-filtered structures, unresolved normalization inputs and
checkpoint-training overlap remain limitations. No fresh panel has been acquired
or selected under this augmented inventory.

The policy retains the complete ordered human single-protein target catalog
(cap 25000), ranks eligible targets by the frozen SHA-256 rule, examines at
most sixty, and retains only the first thousand ordered activities per target.
Every referenced assay ID is joined through lexicographic batches of fifty.
Both release checks must report `ChEMBL_37` and `UP`. Exact raw bytes, URLs,
parameters, response status, pagination and failures remain in the acquisition
ledger. The shared pure selector controls the preregistered stop after six
qualifying targets; offline replay must reproduce that decision ledger exactly.

Activity labeling and direct-human-assay assignment reuse the existing strict
helpers. Query structures require specified stereo, a connected sanitized
structure and nonempty normalized atom/bond scaffold. Family label/endpoint
conflicts reject the family; the earliest surviving publication-year/activity-ID
record supplies its exact state. Deterministic query ranks select two globally
distinct scaffolds per target, twelve across the panel. A target lacking two
available scaffolds is skipped without reserving its tentative scaffold. An
infeasible panel exports no executable public queries, while its partial
selection and complete rejection ledger remain private. No model vocabulary,
3D score, route feasibility or method outcome enters selection.

Known excluded IDs and every declared target alias must first resolve against
hash-bound official target records. Matching accession/component identities
exclude renamed target IDs; missing or conflicting mapping fails clearly.
These endpoint records establish accession/component mapping, not sequence-byte
or protein-family independence. Exact, normalized-family and atom/bond-scaffold
exclusions use the named preparation inventory and reviewed identity kernel;
known unresolved normalization strings and the remaining inventory gaps are
retained explicitly. Generic topology overlap is a reported stratum, not a
rejection. The claim remains absence from the named resolved inventory and
explicitly target-labeled studies, never exhaustive unseen inspection or model
training independence. Public query output contains only IDs, target IDs and
exact SMILES, separate from private activity/assay/selection evidence.

## Remaining delivery and research decisions

1. The inspected pilot root result gate is closed with the final independent
   ledger and all ten selected-output audits. Preserve the negative comparisons
   above. Finish a revised catalog diagnostic with exact blocked-operation
   telemetry; the original failed run must not be relabelled as passed.
2. Review the new, source-bound pilot exclusion inventory for actual coverage,
   including retained structures in previously unclassified fields.
   Freeze the reviewed selection source, exclusion mappings, models, assets,
   runtime, methods, budgets and analysis before fresh acquisition. Resolve
   known target aliases from official records and disclose normalization,
   incomplete inspection and checkpoint-training gaps.
3. Evaluate the proposed six-target/two-query/three-seed panel with five methods
   and two ablations (252 planned attempts), only after the complete panel and
   protocol pass review. Measure joint recommendation success, unique output
   yield, per-target behavior, scaffold changes, similarity distributions and
   full costs. Failures and missing attempts stay in the denominator; six target
   clusters support descriptive uncertainty, not a strong significance claim.
4. Before fresh acquisition, establish and freeze a separately versioned stricter
   structural acceptance policy. The existing development profile is deliberately
   permissive and does not guarantee high
   similarity or retained activity. A stricter profile needs a separately
   declared acceptance objective and frozen external evaluation;
   it must not be fitted to the fresh benchmark outcomes or retroactively
   relabel the existing profile. Active/inactive labels alone do not establish
   structural similarity: inactive compounds can be close structural analogues.
   If structural reference labels are unavailable, explicitly identify any
   stricter thresholds as engineering constraints and report development
   tradeoffs without claiming a calibrated similarity classifier. Keep 2D, 3D
   shape and pharmacophore evidence separate, and report the tradeoff with changed
   scaffolds and output yield.
5. Keep practical synthesis questions explicit: model/rule agreement does not
   supply reaction conditions, yield, selectivity or current purchasing stock.
   Computational route support is the current delivered evidence level.
   Experimental synthesis/activity validation and an additional independent
   synthesis assessment are separate research extensions, with collaboration
   still undecided.

The publication claim must follow these results. If generator ablations and
comparators do not show an advantage, report that outcome and reconsider which
component is the methodological contribution. Source correctness and successful
known examples alone do not establish a top-journal contribution.
