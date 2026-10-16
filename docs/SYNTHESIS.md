# Synthesis-aware scaffold hopping

## Request qualified recommendations

Use `--recommend` on either CLI or `from chembounce_synthesis import recommend`
for the version-fixed `qualified-selective-engineering-v2` contract. The
[asset template](../examples/synthesis/strict-recommendation-config.json) declares
the same profile, so selection is idempotent. Replace its placeholder paths with
existing local files/directories; the example is not a ready-made vendor catalog,
model download or promise that its proposals cover a new input.

```python
from chembounce_synthesis import recommend
from chembounce_synthesis.workflow import load_config

config = load_config("my-strict-assets.json")
rows = recommend(hit, designs, config, "NEW_RUN", deadline_monotonic=None)
```

The API accepts supplied designs, including a lazy iterable. It validates the
request before consuming them and delegates once to the existing final workflow;
it does not generate designs or retry an empty result. JSON paths are relative to
the JSON parent; direct API paths are relative to the current working directory.
`--forward-model-path` is relative to the CLI working directory. Existing generic
`run(...)`, hopping defaults and exploratory examples retain their behavior.

The template uses saved proposals for known designs. To generate proposals for
new hopping designs, replace only its `provider` section with an explicit local
configuration, for example:

```json
{
  "kind": "synllama",
  "model_path": "/path/to/local/SynLlama",
  "prompt_path": "/path/to/local/vars.py",
  "device": "cpu",
  "samples": 2,
  "max_new_tokens": 1600,
  "seed": 42
}
```

Both provider choices still require the separate pinned local ReactionT5 bundle
for strict route support. No model or dataset is acquired implicitly. Early
preflight checks configuration, hit/core/landmarks, alignment capability and local
path types without constructing catalogs, models, reference ensembles or consuming
designs. Asset contents/hashes are still validated by their authoritative loaders;
path presence alone is not certification.

The default original-hit floors are Tanimoto **0.35**, shape **0.70**,
pharmacophore **0.50** and polar reference recall **0.75**. Screening uses eight
conformers at seed 42; final acceptance requires eight at each of 42/43/44. Final
products must change the normalized atom-aware scaffold, which need not change
generic ring topology, and have concrete conditional computational route support.
The actual effective floors and matched/custom status appear in the CLI summary
and complete effective policy remains in the output artifacts. They are engineering
preferences, not activity predictions or experimental synthesis evidence.

`--recommend` rejects a different profile declared in either the original JSON or
CLI, even if ordinary CLI precedence would otherwise hide that declaration.
Unsupported effective modes, seeds, stereo policies or evidence settings fail
clearly. Supported explicit positive numerical/ranking/route-search overrides are
preserved as `custom_unvalidated`; valid resource budgets and extra core/landmark/
property vetoes retain their existing provenance categories. The shortcut never
silently resets conflicting explicit 12/24-conformer settings in old exploratory
examples. Use the new template or explicitly migrate those settings.

`--no-synllama` conflicts with `--recommend`, including with an unavailable config.
`--no-shape` only skips an upstream hopping descriptor; the final geometry gate
still runs and requires an RDKit build with `rdShapeAlign`. The integrated CLI
validates the exact intersection of upstream and synthesis property bounds before
generation and passes that same configuration to final evaluation. It cannot
silently relax an empty property interval.

Continuation with `--recommend` preserves inherited explicit settings and rejects
incompatible ones. An explicit replacement `--config` retains existing migration
and carry-forward semantics. Saved routes are freshly verified/scored; saved
scores, statuses or prior profile labels are never acceptance authority. Carry-only
saved archives remain supported. Deadlines and budgets are unchanged, with no new
top-ups or retries. CLI exit 0 requires final qualified rows, exit 3 means no qualified
rows within the budget, and exit 2 reports an invalid/failed request. Raw hopping
files remain exploratory; final recommendations are in `synthesis/results.*` for
the integrated CLI or `results.*` for the standalone workflow.


The optional synthesis workflow turns scaffold-hopping designs into **catalog and
reaction verified computational candidates**. SynLlama supplies route proposals.
The [SynLlama publication](https://doi.org/10.1021/acscentsci.5c01285) describes
the external model. ChemBounce resolves building blocks, executes and replays registered reactions,
checks the final structure, and ranks regenerated products against both the
original hit and the design using multiple 3D conformers.
Optional independent forward-model evidence and bound-pose geometry checks are
reported separately from this catalog/rule consistency check.

This verifies a specific, reproducible computational claim. It does not establish
experimental yield, selectivity, protection requirements, current stock, biological
activity, or that the selected compound is globally closest in synthesizable space.

## Install and run the small offline example

```bash
python -m pip install -e '.[synthesis]'
python -m chembounce_synthesis.cli \
  --hit 'CC(=O)Nc1ccccc1' \
  --design 'CC(=O)Nc1ccncc1' \
  --config examples/synthesis/demo-config.json \
  --output output_synthesis_demo
```

The demo contains **handwritten proposals and fixture building blocks**, explicitly
labelled as such. It exercises replay, expansion, constraints and scoring; it is
not evidence of model accuracy or commercial availability. Use a new output
directory for each run. Existing nonempty synthesis directories are protected.
The integrated CLI checks this before starting hopping, including when
`--overwrite` is supplied, so a failed synthesis rerun cannot replace the hopping
results underneath an older synthesis result.

The base ChemBounce installation retains support for RDKit 2023.9. The optional
3D alignment workflow needs a build with `rdShapeAlign`; `.[synthesis]` selects
RDKit >=2025.9. Route verification and asset conversion do not need that module.

## Use actual SynLlama inference and a catalog snapshot

```bash
python -m pip install -e '.[synthesis,synllama]'
mkdir -p .chembounce-data-synthesis
git clone https://github.com/THGLab/SynLlama.git .chembounce-data-synthesis/SynLlama
```

Download the trained model/reconstruction archive from the authors'
[Figshare record](https://figshare.com/s/39a37d31cea2c190498d), following their
[inference guide](https://github.com/THGLab/SynLlama/blob/main/assets/docs/inference_guide.md).
Extract `synllama-data` under `.chembounce-data-synthesis`. Model weights and
catalogs remain external assets and are not bundled in ChemBounce. The upstream
[software license](https://github.com/THGLab/SynLlama/blob/main/LICENSE) grants
educational, research and not-for-profit use; model/base-model and supplier data
terms apply separately. The adapter reads the prompt's `TEMPLATE` as an AST literal;
it does not import or execute upstream Python code.

Convert the reaction set matching the chosen model:

```bash
python -m chembounce_synthesis.assets rules \
  --input .chembounce-data-synthesis/SynLlama/data/91_rxns/91_rxn_templates.txt \
  --output .chembounce-data-synthesis/rules.json \
  --source 'THGLab/SynLlama 91-reaction set; record the source commit'
```

The importer retains source line numbers in rule IDs and records duplicate source
lines. For example, the published 115-line set contains 114 unique SMARTS.
Unregistered model templates are rejected; string-similarity substitution is not
used. A reaction registry is an explicit allowlist, not an experimental validation
of each substrate/condition combination.

Supply your inventory or supplier SDF/CSV/TSV with real IDs and dated provenance:

```bash
python -m chembounce_synthesis.assets catalog \
  --input supplier-snapshot.csv \
  --smiles-column smiles --id-column supplier_id \
  --source 'supplier and export description' --snapshot 'snapshot identifier/date' \
  --output .chembounce-data-synthesis/catalog.tsv
```

SMILES line files (`.smi`, `.smiles`, `.txt`) are also supported, with an optional
whitespace-separated identifier after each SMILES.
The catalog importer also accepts the authors' reconstruction `fpindex.pkl`,
`fpindex_<n>.pkl` or `reaction_matrix.pkl`. Repeat `--input` to union per-reaction
indices. Every global referenced by a pickle is replaced with an inert local
placeholder; only known `Molecule` SMILES states are read. No upstream class,
function, model, sklearn object, or pickle reduction is executed. Other schemas
are rejected. This conversion can need substantial memory for large indices.

If an upstream snapshot lacks supplier SKUs, the converter assigns
`structure-sha256:<digest>` identifiers and explicitly records that these are
chemical identities, **not supplier SKUs**. A historical Enamine-derived snapshot
is evidence of membership in that snapshot; use a current supplier export when
purchasability or delivery time matters. No stock-check API is implied.

New conversions write a v2 `<catalog>.manifest.json` sidecar. Its `record_origins`
maps every normalized catalog ID to all accepted source occurrences, each with
`input_index`, `record_index`, and `identifier_origin` (`input_provided` or
`generated_structure_sha256`). `input_index` references the existing `inputs`
path/SHA256 list; repeated input arguments have distinct indices. IDs are stripped
of surrounding whitespace before deduplication, and missing or blank IDs receive
a structure digest. A supplied ID with a `structure-sha256:` prefix remains
`input_provided`. Identical IDs and structures retain all origins, including mixed
supplied/generated occurrences; conflicting structures for the same normalized
ID fail before publication. Source bytes and record coordinates retain the link
to the original ID spelling. `records` counts unique TSV IDs, not occurrences.

The `zero_based_logical_records` policy counts parsed CSV/TSV data records after
the header (multiline records count once), SDF supplier records including invalid
entries, nonblank/noncomment SMILES entries, or inert-pickle molecule positions.
Rejected records consume indices and include `input_index` alongside the existing
`file`, `record`, and `reason`. The sidecar also records `smiles_column` and
`id_column`. Input hashes are checked before and after parsing; a detected change
fails conversion before publishing the TSV or sidecar.

`Catalog.evidence` adds `identifier_provenance` with a `recorded` status, the
sidecar SHA256, and `source_records` expanded to include source paths and hashes.
Route verification and concrete route audits carry these fields for each leaf.
Absent, legacy v1 or unversioned checksum-only sidecars produce `unknown_legacy`
with an empty source-record list; identifier prefixes never fill missing provenance. Malformed v2 data and
unsupported explicit schema versions fail closed. Every v2 ID, including an
unavailable entry, requires origins; accepted and rejected source coordinates
must be unique. Audits detect sidecar changes, removal or addition since loading.
Present directories, other nonregular sidecars and dangling symlinks are rejected
during loading and freshness checks; a symlink to a regular file uses content-hash
semantics. The v2 fields `record_origins` and `record_index_policy` require explicit
`schema_version: 2`, so deleting or downgrading the version cannot silently discard
lineage. Genuine old v1 and unversioned metadata without those fields remain valid.

Both `supplier_sku_verified` and `current_stock_verified` remain false, including
for supplied identifiers and `in_stock` entries. Hashes bind source bytes and
record traceability; a self-consistent rewritten sidecar is not independently
authentic without an external trusted pin. These fields add no experimental
yield, reaction-condition or synthesis-feasibility evidence. Existing TSV files
and v1 sidecars remain readable; newly exported evidence includes additive fields
and therefore receives new audit bindings. Historical artifacts are unchanged.

The normalized TSV requires `id`, `smiles`, `source`, `snapshot`; optional
`availability` is `catalogued` (default), `in_stock`, or `unavailable`. Unavailable
entries are excluded. Even `in_stock` refers to the dated input evidence, not a
live query. Invalid/disconnected molecules are reported during conversion, never
silently desalted. Multiple SKUs for one isomeric structure remain in the catalog;
shortlists use one deterministic representative.

The Python `Catalog` API validates every row and creates its canonical molecular
snapshot during construction. `catalog.fps` is a mutable mapping with the complete
eligible structure key set; length, key iteration and membership checks do not
calculate fingerprints. Each value is computed on first lookup using the
fingerprint generator captured at construction, then cached within that catalog
instance. Reading all values or copying the mapping computes any remaining
fingerprints. Route replay without similarity search does not require them.
There is no disk cache. Callers should treat `catalog.mols` as a fixed snapshot
and use the mapping API for `fps`, without relying on it being a concrete `dict`.

Edit [the model configuration](../examples/synthesis/synllama-config.json) for
your asset paths, model/reaction set, device, and limits. Paths inside a JSON
configuration are relative to that file. Run:

```bash
chembounce-synthesis \
  --hit 'CC(=O)Nc1ccccc1' --design 'CC(=O)Nc1ccncc1' \
  --config examples/synthesis/synllama-config.json \
  --output output_synllama
```

Local Hugging Face weights are required. No implicit model download or remote
model code execution occurs. Safetensors and the authors' PyTorch `.bin` weights
are supported; `.bin` loading uses `weights_only=True` with Torch >=2.6. The
adapter records weight/tokenizer/config hashes, prompt hash, device, library
versions, sampling parameters, seeds, raw outputs and EOS/truncation status.
Four samples cycle through the authors' frugal temperature/top-p settings. More
samples repeat that cycle with distinct seeds. Reproducibility across GPU types
or library versions is not guaranteed by a seed alone.

## Connect to ChemBounce and expand hits

Existing ChemBounce results can be fed directly into the synthesis stage:

```bash
chembounce-synthesis --hit 'CC(=O)Nc1ccccc1' \
  --designs-tsv output_hopping/overall_result.txt \
  --config synthesis-config.json --output output_verified
```

Or append `--synthesis-config synthesis-config.json` to the normal `chembounce`
command. Artifacts then appear in `OUTPUT/synthesis/`. The final CLI exit code
reflects the synthesis stage. The original property bounds and required core
are applied again after reconstruction, intersected with synthesis-specific
bounds. Existing hopping results remain available for auditing.

The main script also accepts an explicit switch:

```bash
# Hopping alone; no synthesis/model/catalog configuration is loaded.
chembounce -i 'CC(=O)Nc1ccccc1' -o output_hopping --no-synllama

# Enable the actual local model, then expand compatible building blocks.
chembounce -i 'CC(=O)Nc1ccccc1' -o output_with_synthesis \
  --synllama --synthesis-config examples/synthesis/synllama-config.json \
  --expansion-top-k 16 --expansion-rounds 1 --expansion-max-routes 128
```

`--synllama` requires a config with `provider.kind: "synllama"`; conflicting or
missing configuration is rejected before hopping. `--no-synllama` skips the entire
synthesis stage even if a config is supplied, which allows one shell script to
toggle the stage. With neither flag, no config means hopping alone; supplying a
config runs its configured provider, preserving the previous behavior. A
`provider.kind: "saved"` config replays stored model proposals without inference;
omit `--synllama` in that case. The three expansion override options also work with
`chembounce-synthesis`.

### Research preset and bounded 2D exploration

Both CLIs accept `--protocol research-v1`. After preparing the local assets above,
edit [research-config.json](../examples/synthesis/research-config.json) for your
model/device and run the integrated protocol:

```bash
chembounce -i 'CC(=O)Nc1ccccc1' -o output_research \
  --protocol research-v1 \
  --scaffold-db tests/fixtures/scaffolds.smi \
  --overall_max_n 20 \
  --synllama --synthesis-config examples/synthesis/research-config.json
```

The fixture scaffold library is for a small operational check. Replace it with
your pinned research library for evaluation. To exercise the synthesis preset
without a model or commercial catalog, use the handwritten fixture:

```bash
chembounce-synthesis \
  --hit 'CC(=O)Nc1ccccc1' --design 'CC(=O)Nc1ccncc1' \
  --protocol research-v1 --config examples/synthesis/demo-config.json \
  --output output_research_fixture
```

The preset makes these choices explicit:

| Stage | `research-v1` defaults |
| --- | --- |
| Hopping | `search_policy="diverse"`, `exploration_fraction=0.25`, `search_candidate_budget=4096` |
| Synthesis allocation | `allocation_policy="balanced"` across scaffold families and designs |
| Refinement | Up to 10 chemical families; seeds 42, 43, 44; 8 conformers initially, then 16 and 32 when escalation is needed |
| Acceptance and ranking | `threshold_policy="all_seeds"`, `ranking_statistic="minimum"` |
| Sampling stability | `max_score_range=0.10`, `max_metric_range=0.15`, `borderline_margin=0.02` |
| Refinement budget | `max_total_conformers=1680` candidate-conformer requests across all stages and fallback designs |
| Similarity acceptance | `minimum_scores={}`; no assumed activity or pharmacophore cutoff |

These spread tolerances are engineering choices for observed sampling stability,
not calibrated biological thresholds. Explicit configuration and CLI overrides
take precedence over preset defaults. In particular, the older
`synllama-config.json` explicitly sets 24 refinement conformers; use
`research-config.json` to obtain the preset's initial 8 without that override.
Setting `protocol` in synthesis JSON configures the synthesis stage only; use the
main CLI's `--protocol research-v1` or Python `protocol="research-v1"` to select
the corresponding upstream hopping policy. Effective settings and protocol
provenance are saved in each manifest.

The default named `legacy` search policy retains mandatory scaffold-bit
`--search-threshold` (0.3) and final count-Tanimoto `-t` (0.5) cutoffs. The `diverse`
policy uses those values to prioritize inexpensive 2D matches while reserving a
bounded fraction of scaffold, intermediate attachment and final candidate slots
for diverse alternatives, including candidates below those cutoffs. It does not
simply set both thresholds to zero. The reservoir and tie-breaking are seeded;
fragment allocation is round-robin. `--exploration-fraction` controls the reserved
fraction; zero disables lower-2D admission. `--search-candidate-budget` caps
requested candidate slots, including requests that produce no valid molecule.
Existing beam, per-replacement, output, core and property limits also apply.
`search_decisions.json` records admissions, exclusions and budget truncation.
This is bounded exploration, not exhaustive retrieval of all possible hops.

### Chemical identity and scaffold constraints

Three identities serve different purposes:

- **Exact state:** canonical isomeric SMILES used for reaction replay, scoring,
  core checks and route-to-product attribution. Input protonation, tautomer,
  stereo and isotope choices are retained.
- **Chemical family:** the versioned normalization groups reversible protonation,
  tautomers, stereo and isotopic variants for novelty and final deduplication.
  Permanent charges are retained. Normalized structures never replace the exact
  state in scoring or reaction replay.
- **Scaffold family:** Murcko scaffolds of normalized families distinguish scaffold
  changes from decoration changes. The generic ring skeleton is an additional
  diagnostic; a heteroatom substitution can change the scaffold without changing
  that generic topology.

Set `mode: "hop"` to require a changed scaffold family versus the original hit.
State-only and decoration-only changes do not qualify. Normalization failures
are explicit: an invalid query fails, and an invalid candidate is excluded with
diagnostics rather than counted as a novel scaffold.
Set `mode: "expand"` to preserve each selected design's Murcko scaffold while
varying blocks; the unchanged design family is excluded. This preservation check
uses the exact isomeric scaffold. `preserve_design_scaffold: true` optionally
keeps that scaffold in hop mode too. `core_smiles` protects a required substructure with specified
stereochemistry. These checks apply to **regenerated final molecules**, not just
LLM proposals. Acyclic scaffolds are empty Murcko scaffolds: hop mode cannot
certify an acyclic-to-acyclic scaffold change; use explicit cores and expansion
mode for that case.

`expansion` accepts `rounds` (0 disables it), `max_routes`,
`building_block_top_k` (alternative blocks per mutable leaf), `min_tanimoto`
(optional block-to-block cutoff, default 0), and
`locked_building_block_ids` (catalog IDs; equivalent structures under other SKUs
are protected too). For a single known route, the Python `expand_route` API also
accepts `locked_nodes` with route-local building-block node IDs. These local IDs
must not be reused as global selectors across different proposals.
For example, merge these keys into an existing synthesis configuration:

```json
{
  "protocol": "research-v1",
  "expansion": {
    "building_block_top_k": 16,
    "rounds": 1,
    "max_routes": 128,
    "min_tanimoto": 0.0,
    "locked_building_block_ids": []
  },
  "minimum_scores": {}
}
```

Each round substitutes one leaf
of each parent route using catalog blocks matching its consuming reaction roles.
All steps are replayed, including changed intermediates. Multiple rounds permit
combinations of substitutions. Neighbours satisfy all consuming reaction roles
before shortlisting. Search uses Morgan-bit Tanimoto with a diverse shortlist tail,
and final core/scaffold checks still apply. Each retained substitution records the
old/new catalog IDs and structures, its Tanimoto and its round. Search/result
capacity is shared across rounds so the first round does not consume the whole
budget. Parent routes and mutable leaves are visited in turn; each frontier prefers
distinct products and higher Tanimoto to the design. Increasing these limits can
increase candidates but cannot guarantee better 3D matches. Explicit
beam/combination/product caps prevent combinatorial growth; reports identify
truncation per round. Repeated neighbour searches and identical starting-route
expansions are cached. `max_routes` and combination bounds apply per starting route;
`max_scored_products` bounds attempted initial exact-product/design comparisons.
Refinement has a separate global candidate-conformer request budget, including
additional stages and alternative designs, as described below.

The Python interfaces are `workflow.run(hit, designs, config, output_dir)`,
`routes.resolve_proposal(...)`, `routes.verify_route(...)`, and
`routes.expand_route(...)`. The original ChemBounce Python API remains compatible.
`workflow.run(..., deadline_monotonic=...)` optionally accepts an absolute
`time.monotonic()` deadline. It stops at work boundaries and exports fully
evaluated candidates plus explicit unassessed evidence. An incomplete refinement
seed block cannot become a successful refined result. A single long operation
can exceed the deadline; use an outer process timeout for a hard limit.
The research benchmark's `scoring_mode="2d"` ablation shares route collection,
constraints, evidence and allocation with the default `"3d"` mode, but uses only
hit/design Morgan Tanimotos and emits no pose SDFs. Its 3D/energy/seed values are
null and refinement is disabled explicitly. Combining this mode with pocket
constraints or 3D `minimum_scores` is rejected rather than silently ignored.

## Refinement and continuation

The older [model configuration](../examples/synthesis/synllama-config.json)
uses 12 screening conformers and 24 refinement conformers at three seeds. The
[research configuration](../examples/synthesis/research-config.json) uses the
adaptive preset above. Configurations without a `refinement` object or a research
preset retain screening-only behavior (`top_n=0`); configurations without
`expansion` disable expansion.

Screening retains every eligible exact-product/design comparison, including
multiple designs or explicit chemical states in one family. With
`allocation_policy="balanced"` (the default), scoring and refinement shortlists
visit distinct scaffold families and then designs in turn, preserving ranking
within those groups. `"encounter"` is an explicit allocation ablation. Carried
routes retain priority over newly projected routes. When a family's preferred
design fails refinement, stability or threshold acceptance, eligible alternatives
are attempted deterministically within the same budget. Final deduplication
selects one accepted exact state/design per chemical family; it does not erase
alternatives before they can be tested.

Each refinement stage evaluates all requested seeds at one ensemble size. Every
seed must succeed before the stage can supply a result. Candidate and reference
ensembles are generated independently per seed; a supplied hit `reference_sdf`
remains fixed. Raising conformer counts does not produce nested ensembles.
The **median-score seed** supplies one actual candidate conformation and its two
reference comparisons. With an even seed count, the lower middle trial supplies
the pose. Exported `score`, shape, pharmacophore and energy components belong to
that trial. Its hit and design poses differ by a rigid alignment, not by candidate
conformation. `seed_score_median` can differ from `score` for an even seed count.

`seed_score_min`, `seed_score_max`, `seed_score_range`, `seed_score_std` and
per-component `metric_ranges` expose sampling sensitivity. With
`ranking_statistic="median"`, Pareto ranking uses the representative pose's
components and the scalar tie-breaker uses the seed median. With `"minimum"`,
Pareto ranking uses the per-component minima across final-stage seeds and maximum
relative energy; `ranking_score` is the minimum seed score. These conservative
ranking diagnostics are stored separately in `ranking_metrics`; they are not
presented as the components of an invented pose.

For explicitly configured 3D `minimum_scores`, `threshold_policy="all_seeds"`
requires every final-stage seed to meet each cutoff. `"representative"` requests
only the representative pose check. Optional `max_score_range` and
`max_metric_range` bound the observed across-seed score and four shape/pharmacophore
component ranges. `convergence_status="converged"` means those engineering range
criteria passed at the final ensemble size. It does **not** establish asymptotic
convergence, a confidence interval, activity probability or biological validity.
Without range criteria, stability is reported as unassessed.

`adaptive_max_conformers` allows larger stages when seeds fail, ranges are unstable
or a configured cutoff is within `borderline_margin` of the sampled component
range. `conformer_growth_factor` controls the increase. Each stage restarts all
seeds; reports retain every trial and the change between stages. Only the final
complete, eligible stage can supply refined results. Adaptive refinement requires
an explicit `max_total_conformers`. A candidate that remains unstable or cannot
complete its necessary stage is not silently replaced by its earlier screening
score.

The budget counts **requested candidate conformers**, not successful conformers,
physical CPU work or reference ensembles. In `research-v1`, 1,680 requests cover
10 families × 3 seeds × (8 + 16 + 32), if each family uses all stages without
fallbacks. Alternative designs share this same budget. Separate counters record
score calls and actual ensemble outcomes. Before escalation or fallback, the
workflow reserves a complete initial seed block for each remaining shortlisted
family. If the total budget cannot cover that initial allocation, a deterministic
balanced subset is evaluated and omitted families receive explicit
`not_evaluated`/budget-exhausted records.

When refinement is enabled, `results.tsv` and the four SDF files contain only
successfully refined shortlisted families. `screening_alternatives.tsv` preserves
all accepted initial exact-product/design scores; `screening.tsv` contains one
initial representative per family. `refinement.json` records attempted alternatives,
stages, failures, thresholds, stability and budgets. `routes.json` also preserves
routes beyond scoring/output budgets. The selected row's `alternative_route_ids`
refer only to the **same exact product**; `alternative_family_candidates` names
other state/design candidates separately, so a stereoisomer's route cannot be
mistaken for a route to the selected state.

Optional 3D cutoffs are checked after refinement; one poor screening pose cannot
prematurely eliminate a candidate. Explicit 2D, core and property constraints still
apply before scoring. No similarity cutoff is enabled by default. Shortlisting
can miss candidates; increasing `--refine-top-n` should be accompanied by a review
of the global conformer budget.

Continue an existing synthesis run without loading SynLlama again:

```bash
chembounce-synthesis --continue-from output_with_synthesis/synthesis \
  --output output_expanded \
  --expansion-top-k 32 --expansion-rounds 2 --expansion-max-routes 256 \
  --synthesis-max-scored 512 --synthesis-max-results 100 \
  --refine-top-n 30 --refine-conformers 24 --refine-seeds 42 43 44
```

`--continue-from` is a standalone synthesis CLI option. It inherits the original
hit, designs, absolute asset paths and settings from `manifest.json`, reuses
`proposals.json` through the saved provider, and carries all previous routes
forward. `--expansion-rounds 2` means two total rounds from projected proposals,
not two additional rounds on every carried route. Save a new output directory;
the earlier run is never overwritten. `--config FILE` can explicitly replace
inherited settings, including the provider; a model provider in that replacement
config will perform new inference. For multiple earlier runs, add
`--include-run DIR` repeatedly (available on both CLIs), or set the config's
`carry_forward` list. Config-relative paths resolve beside the config file.

Increasing neighbour counts or search limits alone does **not** guarantee a
superset, because diverse-tail selection and pruning can change earlier choices.
Carrying routes forward retains that evidence independently of the new search.
Every carried route is replayed against the current catalog and registered rules,
then rechecked for core, scaffold and property constraints. Old scores and
verification claims are never trusted. Ineligible routes are retained with
explicit rejection reasons; eligible products are rescored with current settings.
Carried routes from **all designs** receive scoring priority before any newly
projected routes consume the global scoring budget. Finite scoring and
shortlist budgets can still exclude products from final results; increase those
budgets to compare more of the accumulated pool. A changed hit is rejected;
routes for designs omitted from the current run are counted as excluded in the
manifest. `summary.json` distinguishes carried and new eligible products, initial
screening counts, refinement counts and active limits.

The same refinement and budget flags work in the integrated `chembounce` command.
Use `--refine-top-n 0` for a screening-only run, or increase conformer counts and
seeds for a more expensive sampling check. Users can also provide their own
catalog snapshots and registered reaction SMARTS through the asset formats above;
route verification and core constraints apply to those extensions too.

## Verification and search semantics

1. Parse strict JSON and validate a connected, acyclic reaction dependency graph.
   Reject undeclared leaves, disconnected steps, ambiguous intermediates and
   unregistered reactions. Resolve dependencies by canonical isomeric structure,
   not the order in which the model listed reactions.
2. For each terminal reactant, search only catalog entries matching the required
   reaction role. Morgan fingerprints use radius 2, 2,048 bits and chirality;
   the bounded shortlist reserves a diverse tail. This
   coarse filter is efficient but can miss 3D-similar blocks; increase the limits
   or supply a focused inventory for broader coverage.
3. Enumerate reactant-slot permutations, execute RDKit reaction rules, sanitize
   each principal product, and retain bounded intermediate states. Shared nodes
   must retain a consistent block/intermediate assignment across branches.
   Chirality matching applies during reaction execution as well as catalog
   filtering; a molecule with one matching stereocenter cannot authorize reaction
   at another site that violates the template's stereochemistry.
4. Replay each selected route in a separate verification pass using the registered
   rules. Every leaf must resolve to an eligible catalog ID, every rule to the
   allowlist, and every recorded intermediate/final
   product must match the replay product's **canonical isomeric SMILES exactly**.
   Fingerprint similarity of 1.0 never substitutes for identity.
5. Reapply scaffold/core/property constraints and any configured independent
   synthesis-evidence requirements, then score actual regenerated molecules. Keep
   model-claimed and regenerated structures separate in reports.

Rule replay supports one principal product template and 1–3 reactants.
Replay alone does not validate reaction conditions, stoichiometry,
chemoselectivity, byproducts, kinetics or protecting-group strategy. Product
template heavy atoms without explicit reactant atom-map sources are listed per
step as `unmapped_product_template_atoms`, with a route-level
`has_implicit_atom_sources` flag. Some published rules implicitly supply atoms
from fixed reagents; this flag exposes that missing reagent provenance rather
than claiming every atom originated in a catalogued leaf. Such reagents and
conditions require separate review (`reagents_and_conditions_validated: false`).
Product
enumeration is bounded and may retain multiple regioisomers. An unspecified
stereocenter is not assigned an experimentally justified configuration.

## Independent synthesis evidence

Catalog/rule replay establishes consistency with the allowed reaction templates.
A forward predictor provides a separate check: it receives exact reactants and
supplied reagents, without the claimed product or the generating reaction template.
Its agreement is model support, not a synthesis probability. The research preset
does not automatically load this second model or require an agreement.

The optional local ReactionT5 adapter has its own dependency extra:

```bash
python -m pip install -e '.[synthesis,reactiont5]'
# Add synllama as well when running both models in this environment:
python -m pip install -e '.[synthesis,synllama,reactiont5]'
```

ReactionT5 requires Transformers >=4.51,<5 and Torch >=2.6. The adapter rejects
Transformers 5 because its changed tied-weight loading produced incompatible
predictions for this checkpoint. Supply a compatible local ReactionT5 model
directory containing safetensors, configuration and tokenizer assets; the adapter
does not download weights or execute remote model code. File hashes and runtime
versions are recorded. An optional `model_revision` is a declared attribution,
not a substitute for the recorded file hashes.

Merge this example into your synthesis configuration, replacing the model path:

```json
{
  "synthesis_evidence": {
    "predictor": {
      "kind": "reactiont5",
      "model_path": "models/reactiont5",
      "device": "cpu",
      "top_k": 5,
      "max_new_tokens": 128,
      "max_input_tokens": 512
    },
    "product_matching": "exact",
    "policy": {
      "require_forward_agreement": true,
      "require_complete_context": false,
      "require_explicit_atom_sources": false,
      "require_observed_success": false
    },
    "max_predictor_calls": 128,
    "max_input_tokens": 32768,
    "max_output_tokens": 32768,
    "reaction_contexts": [],
    "observations": []
  }
}
```

Here the predictor's `max_input_tokens` is a per-request limit; the outer token
and call limits are shared across the workflow. Deterministic beam search returns
up to `top_k` predictions. `product_matching="exact"` requires the complete
canonical isomeric predicted product to match. Explicit `"component"` matching
allows the expected product among predicted disconnected components, preserving
all other components in the report. Truncated, missing or budget-exhausted
predictions remain `unknown`; a complete nonmatching prediction is `mismatched`.
Runtime failures remain explicit. None counts as agreement.

Each evidence dimension is reported separately:

| Dimension | What is recorded |
| --- | --- |
| Rule replay | Exact intermediate/product agreement with registered templates and catalog leaves |
| Forward prediction | `matched`, `mismatched`, `unknown` or `failed`, prediction ranks and model/input provenance |
| Reaction context | Explicit reagents and conditions; `supplied_unvalidated` when present, otherwise `incomplete` |
| Atom sources | Template mappings or supplied reagent-atom attributions; missing sources remain `unknown` |
| Observations | Exact reaction/context matches to supplied literature or experiment records: `success`, `failure`, `conflicting` or `unknown` |

All four `require_*` flags default to false; enabling one requires it to pass for
**every route step**. `require_complete_context` checks that a reagents list
(which may be empty) and nonempty conditions object were supplied. It does not
validate those conditions. ReactionT5 consumes reactants and reagents but does
not use temperature/time conditions; reports state this limitation explicitly.
Omit `predictor` to record context/observations without forward inference.

`reaction_contexts` entries bind an exact `reactants` SMILES multiset and `product`
SMILES to optional `reagents` (SMILES list), `conditions` (object), `source` and
`atom_sources`. Reusing a route-local node name such as `r0` would not safely
identify the same chemistry across expanded routes, so config entries use exact
reaction selectors. Each atom-source attribution requires
`template_atom_index`, `reagent_index`, `reagent_atom_index`, `source` and
`record_id`. Indices are zero-based and must refer to the reported unmapped
product-template atom and the supplied reagent SMILES representation; element
matching is checked. Supplied attributions remain unaudited source claims.

Each `observations` entry requires `kind` (`literature` or `experiment`), `source`,
`record_id`, exact `reactants` and `product`, `reagents`, nonempty `conditions`, and
`outcome` (`success` or `failure`). These records are matched to the same exact
reaction/context, never propagated by similarity. Their provenance remains
`source_audited=false`: users must verify the underlying sources. Conflicting
records cannot satisfy `require_observed_success`.

Evidence evaluation visits candidate/route queues in rounds before spending work
on additional routes for one candidate. Cached predictions share the global
budget. When a preferred route fails a requirement, other routes to that exact
product can still qualify; all unsuccessful or unassessed routes remain in
`routes.json`. Exhausted budgets are unknown evidence, not negative chemistry
labels, and cannot silently satisfy a required agreement.

## 3D and pharmacophore methods

Scoring uses explicit input states. Family normalization is used only for
identity/novelty checks; no pH-based protonation, tautomer or stereoisomer ensemble
is generated for scoring. Prepare
and evaluate plausible states explicitly. By default, 3D scoring rejects
unspecified potential stereochemistry; the failure is recorded. Setting
`reject_unspecified_stereo: false` is exploratory and does not validate a
particular stereoisomer.

Before loading the catalog, model or carried routes for a 3D run, the workflow
checks alignment availability and the fixed hit's existing stereochemistry
policy. A supplied `reference_sdf` instead passes the existing exact-state,
stereochemistry and coordinate checks and is reused unchanged. The research
pipeline performs this same input check before generating designs. A failed
precondition is an explicit error, not a successful empty search. This check
neither embeds molecules nor assigns missing stereo, including cage stereo;
`reject_unspecified_stereo: false` and the complete 2D ablation retain their
existing behavior.

The workflow then checks each capped design reference against the same configured
stereochemistry policy, before loading the catalog or calling synthesis/evidence
models. Designs always use generated references, even when the hit has a bound
SDF. A rejected design is skipped individually; valid alternatives keep their
original order and provider sampling seeds. Canonical deduplication and
`max_designs` apply **before** this check, with no backfill from later designs.
Product stereochemistry is still checked during scoring. Passing reference
preflight does not guarantee later conformer embedding, optimization or alignment.

`manifest.json` retains the original capped `designs` for continuation.
`design_preflight` records the requested inputs, canonical unique inputs, capped
and eligible lists, and per-design checks. Rejections appear in `summary.failures`
with `role: "design"` and `stage: "design_reference_preflight"`;
`designs_used` remains the capped count, while `designs_eligible_for_scoring` and
`designs_preflight_rejected` report eligibility. When every capped design is
rejected, the run returns `no_candidates` with explicit reasons and empty result
artifacts. Catalog, rules, provider, evidence model and optional pocket assets are
not initialized or validated in that case; `assets_not_loaded` lists them and
unavailable asset hashes are null. This outcome is not a successful candidate.

Carried route graphs for rejected designs are preserved with
`selection_status: "design_preflight_rejected"` and
`verification_status: "not_evaluated"`; prior verification and scores are not
trusted. The carry source's schema and original hit are still checked. An explicit
later policy change can reconsider those routes. A skipped design has no newly
generated proposal, so continuation that newly enables it must also supply a
provider or saved proposal source containing that design. The complete 2D ablation
does not perform these 3D reference checks; it records `not_applicable_2d`.

For each candidate, original hit and design, generate seeded, single-threaded
ETKDGv3 ensembles. Optimize with MMFF94s when all parameters exist, otherwise
UFF. Discard unconverged/nonfinite conformers and those outside the configured
energy window (default 5 kcal/mol above the lowest sampled energy of that
**same molecule/state**). Energies are not compared as absolute energies across
molecules. They are approximate force-field energies, not binding free energies.

For every retained candidate conformer, align separately to hit and design
reference conformers with [RDKit `rdShapeAlign`](https://www.rdkit.org/docs/source/rdkit.Chem.rdShapeAlign.html)
and `opt_param=0.5` (joint shape
and color optimization). On each aligned pose, compute typed pharmacophore
features on a heavy-atom view; explicit ordinary or isotopic hydrogens are excluded
from shape/color and pharmacophore perception. The rigid transform is transferred
to the full candidate for SDF export, preserving isotope-labelled chemical identity.
This avoids extra features caused solely by explicit D/T atom representation;
isotope-specific kinetic or quantum effects are not modelled. The features use
RDKit `BaseFeatures.fdef`: donor, acceptor, aromatic, positive/
negative ionizable and hydrophobic groups. One-to-one feature assignment uses
Gaussian position weights with configurable sigma (default 1 Å). Aromatic plane
normal agreement multiplies that weight; its sign is ignored. Grouped hydrophobes
replace overlapping atom features. Donor/acceptor direction vectors are not part of this similarity score. Optional
receptor geometry is assessed separately by the pocket constraints below.

If feature-match weight is `M`, pharmacophore overlap is
`M / (N_reference + N_candidate - M)`. Two featureless structures have overlap 1
(the pharmacophore adds no discrimination); one featureless structure has 0.
The per-pose joint score is `sqrt(shape * pharmacophore)`. Shape and pharmacophore
always come from the **same alignment**. Scores from different poses are never
combined by taking independent maxima.

A single candidate conformer must support both references. The selection score
is the geometric mean of the best hit and design joint scores, adjusted by
`exp(-strain_penalty * relative_energy)` (default penalty 0.05 per kcal/mol).
Each selected reference conformer's relative energy is also penalized in its
joint score. A supplied `reference_sdf` replaces the hit ensemble with one exact
3D pose of the same chemical state; its coordinates are preserved without
minimization or a force-field strain estimate.

Final results report `hit_tanimoto` and `design_tanimoto`: Tanimoto similarity of
binary Morgan fingerprints (radius 2, 2048 bits, specified chirality included).
`hit_morgan` remains an alias for `hit_tanimoto` for compatibility. These are
different from the base hopping output's legacy Morgan-count similarity; do not
interpret them as interchangeable scores. The six similarity components are
exported in both the table and candidate SDFs, with 2D scores also in `scores.json`.

Final results receive Pareto fronts over both Tanimotos, hit/design shape and
pharmacophore, relative conformer energy and reaction count. Refinement can use
representative or conservative across-seed components as described above. Within
a front, `ranking_score`, route length and canonical SMILES determine order. These defaults
are transparent ranking choices, not experimentally calibrated thresholds. Use
`minimum_scores` for explicit `hit_tanimoto`, `design_tanimoto`, `hit_shape`,
`design_shape`, `hit_pharmacophore`, `design_pharmacophore` cutoffs and
`property_thresholds` for the existing `mw`, `logp`, `qed`, `sa`, `n_hdonor`,
`n_hacceptor` bounds. SA is an optional property filter, not route verification.
The default `minimum_scores: {}` applies **no hard similarity cutoff**. Earlier
illustrative shape 0.70/pharmacophore 0.50 comparisons are not defaults or validated
quality thresholds. Core retention and catalog/rule replay remain mandatory when
specified regardless of similarity. Optional Tanimoto cutoffs are applied before
3D scoring; 3D cutoffs follow the configured final-stage seed acceptance policy (or the
selected joint pose in a screening-only run). The scalar `score` remains the
3D/energy heuristic above; Tanimoto contributes to Pareto ranking rather than an
arbitrary weighted scalar sum.

## Optional bound-pose pocket constraints

A `pocket` object adds geometric assessment of candidate poses aligned to a supplied
bound hit. It always requires an explicit `reference_sdf` containing one molecule
of the exact hit state. To assess receptor clashes, provide the full PDB containing
that ligand and identify it by chain, residue number and residue name. The SDF
coordinates must already be in that PDB frame: heavy-atom count, elements and
matched coordinates are checked. A ligand-deleted receptor or an independently
embedded reference cannot establish that frame.

Merge this illustrative configuration into an existing synthesis config and
replace the paths and residue selector with those from your own structure:

```json
{
  "reference_sdf": "structures/bound-hit.sdf",
  "pocket": {
    "receptor_pdb": "structures/complex.pdb",
    "reference_ligand": {
      "chain": "A",
      "residue_number": 401,
      "residue_name": "LIG"
    },
    "model": 1,
    "coordinate_tolerance": 0.25,
    "include_waters": false,
    "include_hetero": false,
    "vdw_overlap_tolerance": 0.4,
    "max_clashes": 0,
    "contact_distance": 4.0,
    "max_reported_pairs": 1000,
    "anchors": []
  }
}
```

Distances and overlap tolerances are in Å. `max_clashes: 0` is an explicit example
gate, not a default or an affinity cutoff. Omitting `max_clashes` reports clashes
without excluding poses on that basis. Clash counts are heavy-atom pairs with
van der Waals overlap exceeding the specified tolerance; contact distances are
reported separately. `max_reported_pairs` limits the detailed list, not the count
used for acceptance.

The selected reference ligand is automatically excluded from receptor atoms.
Protein residues are included by default; `protein_chains` can restrict their
chains. Waters and other hetero residues are opt-in. `include_residues` and
`exclude_residues` accept lists of selectors with the same fields as
`reference_ligand`, optionally including `insertion_code`. `reference_altloc`
selects an explicit ligand alternate location; the selected model, alternate
locations and receptor filters are recorded in provenance.

`anchors` can require a candidate atom or pharmacophore feature near a specified
reference position. An atom-anchor entry has exactly these fields:

```json
{
  "id": "reference-carbon",
  "kind": "atom",
  "reference_atom_index": 0,
  "candidate_smarts": "[#6]",
  "candidate_atom_index": 0,
  "max_distance": 1.5,
  "required": true
}
```

This is a schema example, not a recommended chemical anchor. Replace the indices,
SMARTS and distance after inspecting the actual bound ligand. All indices are
zero-based: `reference_atom_index` indexes the scoring reference after ordinary
hydrogen removal and must select a heavy atom; `candidate_atom_index` indexes
the SMARTS query, not the entire candidate. A
feature anchor instead requires exactly `id`, `kind: "feature"`,
`reference_feature_index`, `max_distance` and `required`. The reference feature
index, family, atoms and position are exported in pocket provenance. An optional
anchor (`required: false`) reports its match without gating acceptance. Anchors
alone may be configured without a receptor PDB, but still require a bound
`reference_sdf` and make no receptor-clash claim.

Constraints are evaluated on each sampled, hit-aligned candidate pose before
selecting the best eligible pose, during both screening and refinement. If all
poses fail, the candidate has an explicit scoring failure. Accepted pose geometry
and selection diagnostics are included in `scores.json`. These checks do not
perform receptor relaxation, exhaustive docking or affinity prediction; a failed
alignment does not prove that no feasible binding pose exists.

## Artifacts and evaluation

- `results.tsv`: one selected exact product/design per chemical family, pose
  components, 2D similarities, route IDs, family/scaffold identities, ranking
  diagnostics, representative seed and observed seed spread.
- `screening.tsv`: one representative per initially scored family, before refinement
  and `max_results` truncation.
- `screening_alternatives.tsv`: all accepted initial exact-product/design scores,
  including different chemical states and fallback designs.
- `refinement.json`: attempted alternatives, complete and incomplete stages, seed
  trials, failures, acceptance/stability criteria and candidate-conformer budgets.
- `routes.json`: actual blocks, catalog hash, rule IDs/SMARTS/source, all actual
  reactants/intermediates/products, projection or expansion provenance, final
  constraint/scoring acceptance or rejection reasons, and separately reported
  forward-model/context/observation evidence for evaluated routes.
- `poses_hit.sdf`, `poses_design.sdf`: the selected candidate conformation aligned
  to each reference, linked to its route ID.
- `references_hit.sdf`, `references_design.sdf`: the exact reference conformers
  used for those scores, in matching result order, enabling coordinate-level checks.
- `scores.json`: complete score components, feature counts, conformer IDs,
  energies, force field, settings and ensemble failures; configured pocket
  assessments and candidate-pose selection diagnostics. Every 2D screening,
  3D screening and refined record includes `hit_smiles` and `design_smiles`,
  the canonical workflow input references used for that assessment. Use these
  fields to identify references, including when one product is scored against
  multiple designs; `product` identifies the candidate and `route_id` links its
  route. For compatibility, legacy `hit`/`design` fields retain the 3D metric
  dictionaries; in 2D records the legacy `design` field remains a SMILES string
  and `hit` is absent. Earlier outputs lack the explicit reference fields.
- `proposals.json`: raw model responses and sampling metadata; reusable by the
  `saved` provider without Torch or Transformers.
- `manifest.json`: effective settings, exact chemical-state policy, input/model/
  catalog/rule/code hashes and software/Git versions, protocol/identity policies,
  optional forward-model provenance and bound-pose/receptor frame provenance.
- `summary.json`: status, counts, failed proposals, rejected routes, search
  truncations, unscored products, evidence/refinement usage, active limits and
  elapsed time. Cooperative deadline stops are explicit when requested through
  the Python API.

Exit codes are 0 (results), 3 (valid run without accepted results), or 2 (invalid
request/data or model runtime failure). Model runtime errors are logged and retained
in workflow diagnostics; other unexpected programming failures remain visible.
An empty successful search writes a TSV header and diagnostic files.

Tests cover catalog/rule replay, independent forward-model evidence, multistep
dependencies, unknown rules, missing blocks, false products, stereochemical
mismatches, expansion locks,
post-reconstruction constraints, shared 3D pose scoring, deterministic ensembles,
bound-pose preservation, provenance, empty results and bounded search. Model
accuracy and experimental synthesis success require separate held-out and
laboratory evaluation. For publication, compare projection-only, 2D, shape-only
and joint ranking on the same catalog/rule/model budgets, report failed as well
as successful targets, and use scaffold/time separated evaluation sets. See the executable
[research validation protocol](RESEARCH_VALIDATION.md) for held-out panels,
baselines, ablations, negative controls, compute accounting and failure denominators.

An earlier local-model/full-catalog integration run and exported-coordinate audit
are recorded in [the validation artifacts](../benchmarks/results/synthesis/README.md).
The small case produced five scaffold-hopping products and four expansion products.
The later [expansion comparison](../benchmarks/results/synthesis-expansion/README.md)
withholds the exact precursor, protects a specified core, and increases returned
analogues from 4 to 35 by widening the neighbour and round budgets under the updated
algorithm. It reports both Tanimotos and all 3D components, including the tradeoff
between design and original-hit pharmacophore similarity.
These are integration results, not an estimate of general success rate.

The earlier [multi-seed refinement and continuation validation](../benchmarks/results/synthesis-refinement/README.md)
uses the same withheld-precursor case with 12-conformer screening and 24-conformer
reevaluation at three seeds. It retains all four earlier candidates, screens 14
products in the first stage and 36 after continuation, and refines the top 10 in
each stage. All 14 first-stage products remain in the continued pool. Seed ranges
and independently checked exported coordinates accompany the results; increasing
conformer sampling changes earlier low-budget 3D conclusions. These historical
artifacts describe their pinned versions and settings, not an evaluation of the
new `research-v1` preset.

### Offline candidate review

Create a static, human-readable report from a saved synthesis run:

```bash
python -m chembounce_synthesis.report path/to/synthesis-run
# Or write a new HTML file elsewhere:
python -m chembounce_synthesis.report path/to/synthesis-run --output review.html
```

The default output is a new `review.html` inside the run directory. Existing files
are never overwritten. The Python API is `chembounce_synthesis.report.write_report(
run_dir, output=None)`, returning a structured review view. `build_report(run_dir)`
returns that view without writing; `render_html(view)` renders it without links.
No chemistry library, model, catalog, route replay, scoring or network access is
needed. This opt-in command leaves generation, acceptance and existing exports
unchanged.

Candidate cards show the exact original hit and regenerated product, recorded
scaffold information, effective floors and profile status, representative pose
metrics, each final seed, concrete route steps and catalog identifier origins.
Each seed has its own common-pose witness; the report does not claim that different
seeds share one pose. Minima are labelled arithmetic derived from the saved final
seed values, separately from the representative exported pose. Missing original-hit
scaffold or generic-topology comparisons stay **not recorded**; no chemical
identities or scores are recomputed. Legacy catalog provenance remains unknown,
and recorded availability or supplier labels do not verify SKU or current stock.

The report distinguishes **workflow-recorded qualification** from **artifact
consistency**. It checks supported JSON schemas, exact product/design joins,
primary route bindings, saved representative/seed values and witness references.
It does not scientifically validate geometry, chemical identity, reactions or
forward predictions. SDF and TSV files are linked but not parsed or verified.
Hashes identify the same bytes used for the report; they do not authenticate a
self-consistent bundle. An inconsistent or partial bundle displays diagnostics and
unconfirmed recorded claims. Zero-result reports retain error, rejection, budget
and deadline distinctions without inferring new chemical causes.

Only the seven named JSON artifacts and the optional `results.tsv` and four pose/
reference SDF files are read. Each file is limited to 64 MiB and the bundle to
256 MiB. Unsafe symlinks, changes during rendering and output collisions are
refused. Local links are emitted only when the output is directly in the run
directory; an external output shows filenames and hashes instead. Technical
provenance is collapsible, with no scripts, external fonts or other resources.
Exit 0 means the recorded bundle passed the report's consistency checks, **not**
that the workflow succeeded or chemistry was revalidated; exit 3 means a diagnostic
report was written for incomplete/inconsistent artifacts; exit 2 means the report
could not safely be written. These report exits are separate from workflow exits.

Saved scalar consistency includes `[0, 1]` similarity domains, configured floors,
complete conformer settings, and recorded final-stage size/seed correspondence.
Custom fixed sizes are compared with their policy; adaptive sizes require the
recorded growth sequence and final-stage evidence. Missing or unsupported settings
remain diagnostic rather than being filled from current defaults. Step coverage,
products, rules and reactant slots are compared across the saved route graph,
verification, audit and support records. These are record comparisons, not chemical
replay. Final accepted/qualified counts may exceed returned counts under output
caps, but cannot contradict them; partial running/error counters remain distinct
from terminal no-candidate counts.
Qualified reports require positive original-hit Tanimoto, shape and pharmacophore
floors; optional design floors can be zero, as in the workflow policy. Each displayed
catalog origin is joined through the saved route leaf's exact building-block ID to
its catalog record ID. This prevents borrowing another leaf's recorded origins;
it does not authenticate the identifier or infer supplier provenance from its text.
