# Qualified recommendations: implementation and evidence gates

The complete objective is to recommend final, route-regenerated molecules that
remain similar to the **original input hit**, change its normalized scaffold,
preserve configured anchors, and have concrete synthetic-route evidence. A
qualification filter alone does not establish useful yield, search quality,
biological activity, or experimental synthesizability.

## Six-stage delivery tracker

| Stage | Required outcome and evidence gate | Status |
| --- | --- | --- |
| 1. Unified final qualification | One versioned conjunction at the real API/CLI export boundary; final-product/original-hit similarity, normalized hop, anchor, exact selected route, requested evidence and complete scoring/seed checks. Adversarial workflow tests, alternative survival, continuation, and matching JSON/TSV/SDF selection. | Implemented; review pending. Explicit engineering thresholds only. |
| 2. Stereo/design eligibility and bounded refill | Explicit stereo/design eligibility with an auditable, bounded replacement supply when designs fail; preserve caps and complete rejection denominators. | Pending. Existing design preflight does not refill. |
| 3. Similarity-aware route/search | Search and route construction retain original-hit similarity while finding changed scaffolds and viable route products; compare yield and compute under fixed budgets. | Pending. Stage 1 gates outputs but does not improve search. |
| 4. Calibrated acceptance and hit-prioritized ranking | Calibrate thresholds using development-only controls, prioritize original-hit acceptance/ranking, and ship a usable recommendation preset and CLI default flow through this contract. Previously inspected data cannot be described as held-out evidence. Multi-seed verification and threshold provenance must support the intended profile. | Pending. No calibrated cutoffs or production recommendation preset is supplied by Stage 1. |
| 5. Independent synthesis evidence | Establish independent, route-specific evidence and provenance beyond catalog/rule replay; define treatment of unknowns, missing context, atom sources, model disagreement and observation quality. Never present a heuristic as a synthesis probability. | Pending. Stage 1 enforces the existing explicitly requested evidence policy. |
| 6. Fresh fair benchmarks and ablations | Fresh artifact directories; comparable proposals, budgets, splits and resources; full input/proposal/route/scoring/seed/qualification failure denominators; fair ablations and confidence/uncertainty reporting. | Pending. No scientific campaign or performance improvement claimed. |

Each stage requires implementation, focused validation and independent review
before the next stage. Historical published results and source inventories are
not revised to imply these stages were present.

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
