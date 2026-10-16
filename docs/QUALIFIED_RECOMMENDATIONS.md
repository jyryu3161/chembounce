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
| 4. Calibrated acceptance and hit-prioritized ranking | Calibrate thresholds using development-only controls, prioritize original-hit acceptance/ranking, and ship a usable recommendation preset and CLI default flow through this contract. Previously inspected data cannot be described as held-out evidence. Multi-seed verification and threshold provenance must support the intended profile. | Scorer/ranking and frozen development fit gates passed. Reviewed artifact packaged; release integration gate pending. Explicit opt-in, permissive development thresholds—not a strong-similarity guarantee. |
| 5. Independent synthesis evidence | Establish independent, route-specific evidence and provenance beyond catalog/rule replay; define treatment of unknowns, missing context, atom sources, model disagreement and observation quality. Never present a heuristic as a synthesis probability. | Pending. Stage 1 enforces the existing explicitly requested evidence policy. |
| 6. Fresh fair benchmarks and ablations | Fresh artifact directories; comparable proposals, budgets, splits and resources; full input/proposal/route/scoring/seed/qualification failure denominators; fair ablations and confidence/uncertainty reporting. | Pending. No scientific campaign or performance improvement claimed. |

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
Stage 5 must establish the production evidence profile.

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

## Stage 4c recommendation profile (release integration gate pending)

The API field `recommendation_profile: "qualified-development-v1"` and both
CLIs' `--recommendation-profile qualified-development-v1` select one named
development profile through **explicit opt-in**. Catalog, rules and
proposal-provider configuration remain explicit; no manual cutoffs are needed.
The default exploratory configuration and explicit engineering-qualified mode
remain available. Missing or corrupt profile assets fail clearly rather than
substituting values or falling back.

```sh
chembounce-synthesis --hit "$HIT_SMILES" --design "$DESIGN_SMILES" \
  --config synthesis-assets.json --recommendation-profile qualified-development-v1 \
  --output fresh-development-run
chembounce -i "$HIT_SMILES" -o fresh-development-integrated \
  --synthesis-config synthesis-assets.json \
  --recommendation-profile qualified-development-v1
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
