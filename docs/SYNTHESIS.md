# Synthesis-aware scaffold hopping

The optional synthesis workflow turns scaffold-hopping designs into **catalog and
reaction verified computational candidates**. SynLlama supplies route proposals.
The [SynLlama publication](https://doi.org/10.1021/acscentsci.5c01285) describes
the external model. ChemBounce independently resolves building blocks, executes registered reactions,
checks the final structure, and ranks regenerated products against both the
original hit and the design using multiple 3D conformers.

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

The normalized TSV requires `id`, `smiles`, `source`, `snapshot`; optional
`availability` is `catalogued` (default), `in_stock`, or `unavailable`. Unavailable
entries are excluded. Even `in_stock` refers to the dated input evidence, not a
live query. Invalid/disconnected molecules are reported during conversion, never
silently desalted. Multiple SKUs for one isomeric structure remain in the catalog;
shortlists use one deterministic representative.

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
`chembounce-synthesis`. Existing base hopping 2D cutoffs remain separate.

The normal hopping search still uses its configured 2D cutoffs. For exploration
of substantially different cores, choose appropriate `-t` and
`--search-threshold` explicitly; a downstream 3D method cannot recover designs
already discarded upstream. No cutoff is automatically relaxed.

Set `mode: "hop"` to require a changed Murcko scaffold versus the original hit.
Set `mode: "expand"` to preserve each selected design's Murcko scaffold while
varying blocks; the unchanged design itself is excluded. `preserve_design_scaffold: true` optionally keeps that scaffold
in hop mode too. `core_smiles` protects a required substructure with specified
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
For example, an exploratory configuration is:

```json
"expansion": {
  "building_block_top_k": 16,
  "rounds": 1,
  "max_routes": 128,
  "min_tanimoto": 0.0,
  "locked_building_block_ids": []
},
"refinement": {"top_n": 10, "num_conformers": 24, "seeds": [42, 43, 44]},
"minimum_scores": {}
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
`max_scored_products` bounds initial product/design comparisons. Refinement has
its own budget: at most `top_n` products times the number of seeds.

The Python interfaces are `workflow.run(hit, designs, config, output_dir)`,
`routes.resolve_proposal(...)`, `routes.verify_route(...)`, and
`routes.expand_route(...)`. The original ChemBounce Python API remains compatible.

## Refinement and continuation

The supplied **model configuration** starts with 16 alternative blocks per leaf,
one substitution round and 12 conformers at seed 42 for screening. The top 10
unique products are then reevaluated with 24 conformers at seeds 42, 43 and 44.
This is an adjustable starting protocol, not a universal optimum. Historical
configurations without a `refinement` object retain their previous single-seed
behavior (`top_n=0`); configurations without `expansion` still disable expansion.

Screening uses the existing Pareto ranking. Refinement reevaluates each shortlisted
product against its selected design and the original hit. Every requested seed
must produce a valid score; numerical failures are recorded and that candidate
is excluded from refined results. Each seed independently builds candidate and
reference ensembles; a supplied reference SDF remains fixed. Scoring and core
rules are unchanged. Raising conformer counts does not produce nested ensembles.

The **median-score seed** supplies one representative candidate conformer and both
reference comparisons. With an even seed count, the lower middle trial supplies
the pose. The six similarity components and `score` always belong to that actual
trial; they are never assembled from different seeds. `seed_score_median`,
`seed_score_min` and `seed_score_max` expose sampling sensitivity. Refined products
are ranked by Pareto fronts of their representative components, then median seed
score, route length and SMILES. The median is a heuristic summary, not a confidence
interval or activity probability. At the default odd seed count, `score` equals
`seed_score_median`; they can differ for even counts.

`results.tsv` and the four SDF files contain **only successfully refined shortlisted
products** when refinement is enabled. All scored products remain in
`screening.tsv`; all route evidence, including candidates beyond scoring/output
budgets, remains in `routes.json`. This avoids mixing screening and refinement
scores in one final ranking. `refinement.json` contains every seed's components,
settings, diagnostics and failures. The reference SDFs use the representative
seed's exact conformers. Optional 3D cutoffs apply after refinement, so a poor
initial pose alone cannot trigger premature rejection. 2D, core and property
constraints still apply before scoring. No similarity cutoff is enabled by default.
Shortlisting can miss candidates; increase `--refine-top-n` to reevaluate more.

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
Within each design, carried routes receive scoring priority. Finite scoring and
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
4. Independently replay each selected route. Every leaf must resolve to an eligible
   catalog ID, every rule to the allowlist, and every recorded intermediate/final
   product must match the replay product's **canonical isomeric SMILES exactly**.
   Fingerprint similarity of 1.0 never substitutes for identity.
5. Reapply scaffold/core/property constraints, then score actual regenerated
   molecules. Keep model-claimed and regenerated structures separate in reports.

The initial implementation supports rules with one principal product template
and 1–3 reactants. It does not validate reaction conditions, stoichiometry,
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

## 3D and pharmacophore methods

Chemical states are explicit input states. No automatic neutralization,
pH-based protonation, tautomer or stereoisomer enumeration is performed. Prepare
and evaluate plausible states explicitly. By default, 3D scoring rejects
unspecified potential stereochemistry; the failure is recorded. Setting
`reject_unspecified_stereo: false` is exploratory and does not validate a
particular stereoisomer.

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
replace overlapping atom features. Donor/acceptor direction vectors and receptor
contacts are **not** included in this version.

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
pharmacophore, relative conformer energy and reaction count. Within a front they are ordered by
the recorded heuristic score, route length and canonical SMILES. These defaults
are transparent ranking choices, not experimentally calibrated thresholds. Use
`minimum_scores` for explicit `hit_tanimoto`, `design_tanimoto`, `hit_shape`,
`design_shape`, `hit_pharmacophore`, `design_pharmacophore` cutoffs and
`property_thresholds` for the existing `mw`, `logp`, `qed`, `sa`, `n_hdonor`,
`n_hacceptor` bounds. SA is an optional property filter, not route verification.
The default `minimum_scores: {}` applies **no hard similarity cutoff**. Earlier
illustrative shape 0.70/pharmacophore 0.50 comparisons are not defaults or validated
quality thresholds. Core retention and catalog/rule replay remain mandatory when
specified regardless of similarity. Optional Tanimoto cutoffs are applied before
3D scoring; 3D cutoffs use the selected joint pose. The scalar `score` remains the
3D/energy heuristic above; Tanimoto contributes to Pareto ranking rather than an
arbitrary weighted scalar sum.

## Artifacts and evaluation

- `results.tsv`: unique selected products, 3D components, 2D diagnostic similarity,
  relative energy, route IDs, scoring stage, representative seed and seed-score range.
- `screening.tsv`: all accepted initial scores before refinement and `max_results` truncation.
- `refinement.json`: every reevaluated candidate and all requested seed trials.
- `routes.json`: actual blocks, catalog hash, rule IDs/SMARTS/source, all actual
  reactants/intermediates/products, projection or expansion provenance, final
  constraint/scoring acceptance or rejection reasons.
- `poses_hit.sdf`, `poses_design.sdf`: the selected candidate conformation aligned
  to each reference, linked to its route ID.
- `references_hit.sdf`, `references_design.sdf`: the exact reference conformers
  used for those scores, in matching result order, enabling coordinate-level checks.
- `scores.json`: complete score components, feature counts, conformer IDs,
  energies, force field, settings and ensemble failures.
- `proposals.json`: raw model responses and sampling metadata; reusable by the
  `saved` provider without Torch or Transformers.
- `manifest.json`: effective settings, exact chemical-state policy, input/model/
  catalog/rule/code hashes and software/Git versions.
- `summary.json`: status, counts, failed proposals, rejected routes, search
  truncations, unscored products and elapsed time.

Exit codes are 0 (results), 3 (valid run without accepted results), or 2 (invalid
request/data). Unexpected programming/model-runtime failures remain visible.
An empty successful search writes a TSV header and diagnostic files.

Tests cover independent forward replay, multistep dependencies, unknown rules,
missing blocks, false products, stereochemical mismatches, expansion locks,
post-reconstruction constraints, shared 3D pose scoring, deterministic ensembles,
bound-pose preservation, provenance, empty results and bounded search. Model
accuracy and experimental synthesis success require separate held-out and
laboratory evaluation. For publication, compare projection-only, 2D, shape-only
and joint ranking on the same catalog/rule/model budgets, report failed as well
as successful targets, and use scaffold/time separated evaluation sets.

An actual local-model/full-catalog integration run and exported-coordinate audit
are recorded in [the validation artifacts](../benchmarks/results/synthesis/README.md).
The small case produced five scaffold-hopping products and four expansion products.
The later [expansion comparison](../benchmarks/results/synthesis-expansion/README.md)
withholds the exact precursor, protects a specified core, and increases returned
analogues from 4 to 35 by widening the neighbour and round budgets under the updated
algorithm. It reports both Tanimotos and all 3D components, including the tradeoff
between design and original-hit pharmacophore similarity.
These are integration results, not an estimate of general success rate.

The [multi-seed refinement and continuation validation](../benchmarks/results/synthesis-refinement/README.md)
uses the same withheld-precursor case with 12-conformer screening and 24-conformer
reevaluation at three seeds. It retains all four earlier candidates, screens 14
products in the first stage and 36 after continuation, and refines the top 10 in
each stage. All 14 first-stage products remain in the continued pool. Seed ranges
and independently checked exported coordinates accompany the results; increasing
conformer sampling changes earlier low-budget 3D conclusions.
