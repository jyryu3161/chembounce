# Synthesis-aware scaffold hopping

The latest fixed 12-input development run completed with zero qualified
products despite saved computational route support. See the
[synthesis-route priority plan](SYNTHESIS_ROUTE_PLAN.md) for the baseline,
bounded route-generation work, and checks planned before another large run.

## Experimental route-first generation (explicit opt-in)

The `catalog-rule-route-first-v2` option generates concrete one- and two-step
catalog/rule routes before the existing recommendation intake. It does not run
native scaffold hopping or SynLlama proposals. Use it only with `--recommend`
and a synthesis configuration that declares `route_first`; omit `provider` and
route `expansion`. For example, adapt the [strict configuration](../examples/synthesis/strict-recommendation-config.json)
with local catalog, rules, and pinned forward-model assets, remove `provider`,
and add:

```json
"route_first": {"version": "catalog-rule-route-first-v2", "max_steps": 2}
```

```bash
python chembounce.py -i 'HIT_SMILES' -o NEW_RUN --recommend \
  --route-first --synthesis-config route-first.json
python -m chembounce_synthesis.cli --hit 'HIT_SMILES' --recommend \
  --route-first --config route-first.json --output NEW_SYNTHESIS_RUN
```

The Python entry point is `recommend_routes(hit, config, output_dir)` from
`chembounce_synthesis`. The policy has explicit global caps on rule scans,
catalog index rows, charged bucket/role visits and fingerprints, combinations,
products and route-audit attempts,
beam width, two-step depth and the inherited wall deadline. Reverse-template
products are retrieval hints only; every retained route is forward replayed
and audited, and its evidence and 3D gates remain separate. An exact generated
product is its own design reference, so design-to-product similarity 1 is
tautological; original-hit constraints remain the meaningful comparison.

Specified hit tetrahedral and E/Z centers compile unique, self-matching local
anchors from the original hit before any product is considered. The radius
grows only within fixed atom/match bounds until source stereo is representable;
the queries are then frozen. The comparison includes complete rings touching
the context and may conservatively reject changed or ambiguous anchors,
including some legitimate peripheral replacements. It does
not infer a global hit/product atom mapping. `route_first_candidates.json`
records reviewable generated routes and work counters; `results.json` contains
only products that pass the existing final qualification. This material
generation override is reported as `custom_unvalidated`; a development-panel
comparison does not establish independent calibration or experimental synthesis.
V1 produced no supported products on the reused 12-input development panel.
The seventh v2 development probe retained coverage at three of 12 inputs before
model and 3D checks. Distinct input/product pairs rose from 10 to 30, including
22 eligible two-step routes, all on one input. Full eligible route graphs and
near rejections are retained in the diagnostic record. More than five
pre-evidence inputs is required before a model run, so no forward/3D comparison
was performed for this revision.
Adoption requires more than five of the 12 inputs to have complete, audited
changed-scaffold routes with all-step forward support. The same
original-hit similarity, property, polar-capacity, route, forward-support and
3D gates still apply. Reverse outcomes are checked for a complete, disjoint
partition of the target atoms. A local copy of each reverse template can make
one role-conforming retrieval hint by replacing an originless terminal query
placeholder with a simple allowed atom. When that bounded construction fails,
the previous safely normalized approximate hint remains available. Typed and
approximate hints are both search guidance; neither establishes a reactant or
a viable reaction. Explicit fallback hints are reported separately.

Typed hint construction tries neutral OH, F, Cl, Br, I, then O− in that fixed
order. Each precursor-part attempt is either one representative construction
or one check for a part without a placeholder. At most 16 attempts are charged
across all parts of a raw reverse outcome and 1,024 across one search,
including rejected attempts.
Each role check enumerates at most 65 matches; reaching match 65 means the
enumeration may be truncated and rejects typed construction. These bounds do
not raise the existing reverse-outcome, catalog, product, fingerprint or route
caps.

An eighth v2 development probe that changed the second catalog-block choice
reached only two of 12 inputs before model and 3D checks, so that change was
reverted. The retained revision instead allocates the existing eight
outer search plans across distinct original-hit cut bonds recovered from the
same capped reverse outcomes. It prefers complete role-conforming hints while
retaining approximate, unary and explicit fallback plans. Each original-parent
origin atom or bond inspected for cut metadata is charged to a 32,768-check
search cap; unavailable metadata remains labelled unknown. No extra reverse
outcomes or catalog rows are read for this decision. A catalog-free control
retains both central typed aryl–amine disconnections that the prior plan order
omitted on one development input. That control is a plan-selection result, not
evidence that the catalog search or final gates produce a route. This revision was subsequently evaluated together with the topology screen
in the ninth fixed-panel development probe described below.

The next bounded revision screens inner rule matches by an explicit graph
condition. For a supported two-reactant template, it removes all mapped
cross-reactant product bonds from each matched target graph and rejects the
rule only when every completely enumerated match still connects mapped atoms
from different reactant slots. Unclear maps, retained-atom bond deletions,
truncated matches, unsupported reactions and exhausted work budgets keep the
existing rule eligibility. Disappearing leaving-group atoms do not count as
retained bond deletions. At most 128 extra match enumerations per search use
64 matches plus one truncation sentinel each, with 262,144 shared rule-template
and target-graph visits. The existing 512 rule inspections and four inner
rules remain capped as before.

For a statically supported template with exactly one such cross-reactant bond,
one local reverse-query copy can require that bond to be nonring. This removes
ring cuts that cannot partition the original target while keeping the same
single reverse call and four raw outcomes. A zero-output restricted call cannot
create an implicit whole-hit fallback; explicitly selected fallback rules
remain available. Original rules, catalog molecules, forward replay and audit
are untouched. Catalog-free original-hit controls found an already selected
coupling rule with interaromatic raw cuts and a fourth-ranked coupling rule
with typed interaromatic hints after this screen. They establish search-center
reachability under the local bounds, not catalog availability or a complete
supported route. The ninth v2 development probe using this policy reached
three of 12 inputs before model and 3D checks, with 38 eligible products and
routes (16 one-step and 22 two-step). It did not meet the fixed entry condition
of more than five pre-evidence inputs, so no model comparison was run.
Before the topology screen, a separate catalog-free check found that the
needed inner coupling rule was outside the four-rule selection on those two
guides. With the original reverse query and unchanged four-outcome cap, one
guide yielded no partition-valid outcome and the other only a peripheral
methyl-bond cut. The revised local nonring query changes that bounded center
enumeration, as shown above, without demonstrating recovery of the saved
two-step route.

The ninth probe also exposed shared-budget starvation in the two-step search.
On one input, a central outer plan selected four inner rules, including the
intended coupling rule. Those four rules generated nine guides in total, yet
the slot's 18 catalog queries received no role checks
after the preceding plan spent most available catalog-row visits. The next
revision prepares only the already bounded inner rules and reverse guides,
orders typed guides first within each outer plan, and executes one guide per
plan per sweep. It records preparation and actual per-plan catalog work
separately. A preparation limit flag records that a shared cap was observed;
it does not claim a count of unenumerated guides. Deferred-item counts include
only guides already prepared but not attempted. The original global rule,
reverse, catalog, fingerprint, product, combination, deadline and qualification
limits remain in force. Preparation
can still consume reverse-normalization and rule-inspection budgets before
execution, so this scheduler does not guarantee that every plan gets a full
catalog query. The tenth development probe completed on the fixed panel with
pre-evidence coverage of four of 12 inputs, 35 distinct input/product pairs
and 40 eligible routes (16 one-step and 24 two-step). The fixed entry condition
of more than five pre-evidence inputs remains unmet, so no model comparison
was run and P1 remains incomplete.

Strict route audit requires every concrete source reactant to have fully
specified stereochemistry, even when its unresolved center disappears from
the principal product. The source stereo precheck applies this necessary condition to
already visited, role-compatible catalog blocks before choosing the nearest
two blocks for any route step. A block proven to fail that source-stereo requirement is
skipped and the same ranked pool supplies the next option. Perception failure,
deadline or budget exhaustion leaves the block at its original rank; no
stereoisomer is assigned. At most 2,048 source-block stereo perceptions and
131,072 charged atom/bond units are allowed per search, with at most 2,048
cached completed results and 16 recorded exclusion examples. This adds no
catalog scan, role match enumeration or reaction call and does not change the
strict audit or other final gates. The eleventh fixed-panel development probe
completed with the same four eligible inputs, 35 exact input/product pairs
and 40 routes as the tenth probe; the eligible product identities were
unchanged. Across the 12 inputs, 672 stereo perceptions produced 183 exclusion
decisions, including repeated blocks across queries. Unspecified-product
first-failure records fell from 50 to 14, but this did not improve eligible
coverage. No model or 3D comparison was run; P1 remains incomplete.

The twelfth route-first revision changed only bounded catalog retrieval.
It selects the two most frequent exact reactant roles among the already chosen
outer plans, with normalized SMARTS as the deterministic tie-break. Atom maps
are removed from those role queries; recursive, stereo and other SMARTS
conditions remain. The existing catalog index pass checks those two roles at
most 500,000 additional times in total. If either index is incomplete or fails,
retrieval uses the previous full-catalog pool. A matching-role subset is ranked
first by the existing bucket distance and stable hash order, then ordinary rows
fill any remaining places without duplication. The total pool stays at 8,192.

For each reverse guide, a bounded role-matched atom window keeps the matched
atoms and immediate neighbors, closes ring systems sharing atoms, preserves
specified stereo centers and bond controls, then adds one final neighbor step.
This atom set supplies only four coarse retrieval counts. The policy takes the
componentwise arithmetic midpoint, rounded up, between those raw counts and
the full normalized guide's counts before applying the existing bin sizes.
The original full guide still supplies fingerprint ranking and reactant-role
checks. An indexed role can also be reused for a different caller role only
when complete bounded matches on the full guide have an identical atom-index
set at a site; this is a retrieval heuristic, not proof that the roles entail
each other. Exact role equality takes priority. Incomplete or failed metadata
keeps the legacy retrieval choice. Limits are 128 guide/role pairs, 128 window
and 384 site match requests, 64 matches plus a truncation sentinel per request,
256 guide atoms, 512 bonds and 262,144 combined graph visits. These add work
within the inherited deadline; catalog-row, fingerprint, reaction, product and
route limits and all qualification gates are unchanged.

Saved-metadata prototypes located the four catalog leaves from two previously
supported routes outside the old pool for their needed queries. The fixed
midpoint plus exact-first same-site calculation put all three leaves of one
saved route inside at least one relevant capped query, including a guide reused
from the first stage in a later inner stage. It did not run candidate fingerprint
ranking, forward reactions or route audit, and does not establish that one
work item selects those three leaves or recovers the route. The current code
was then evaluated on the fixed panel. The twelfth probe reached five of 12
inputs, 41 distinct input/product pairs and 41 eligible routes (20 one-step,
21 two-step). It gained four two-step products on one previously uncovered
input, while another input lost five eligible products. The latter generated
120 exact products but reached the 65,536 total fingerprint budget before any
product admission was observed. No model or 3D comparison was run because the
fixed entry condition of more than five inputs was not met; P1 remains
incomplete.

The next untested revision reserves fingerprint work for product admission:
`min(2 × max_products, max_fingerprint_calls // 2)` calls, leaving a separate
generation ceiling. For two-step searches, the first-stage ceiling is the smaller of half the total
and that generation ceiling; a one-step search uses the full generation ceiling. The initial hit fingerprint counts as generation
setup work. Admission can also use generation allowance that was not spent,
without exceeding the same global cap. Once generation reaches its ceiling, already pending products can
still be ranked and admitted. A deterministic subset considers cached product
fingerprints first, then canonical SMILES; each selected product reserves its
zero-or-one new rank fingerprint and at most one admission-callback fingerprint.
Unranked products are recorded as fingerprint-unassessed, not chemical
rejections. Ranked, admitted products still undergo the unchanged concrete
forward replay and strict route audit. The total fingerprint cap, deadline,
product and route caps, and all qualification gates remain unchanged. Separate
generation, admission-ranking and admission-callback counters make budget use
reviewable. This revision has not yet been run on the fixed panel.

## Optional catalog byte binding

A configuration can pin both the catalog TSV and its import manifest:

```json
"catalog_binding": {
  "catalog_sha256": "<64 lowercase hexadecimal characters>",
  "manifest_sha256": "<64 lowercase hexadecimal characters>"
}
```

When present, this object requires exactly both hashes; null, missing fields and
unknown fields are rejected during configuration validation. The loaded catalog
checks both pins before constructing molecules, then parses the exact captured
TSV and manifest bytes that matched them. Missing or changed sidecars fail instead
of falling back to legacy provenance. Existing route audits still check on-disk
freshness. Unpinned configurations retain the existing streaming/legacy behavior.

These pins establish byte identity, not supplier SKU, current stock, source
completeness or experimental synthesis validation. An explicitly pinned legacy
manifest still reports unknown legacy record provenance. The declaration is
saved in `manifest.config` and inherited by continuation; changing only the
catalog path retains the pins. Providing a separate replacement configuration is
an explicit new configuration, which must declare its own binding when wanted.

Recommendation preflight validates the declaration without reading catalog
contents. Actual pin enforcement occurs at catalog loading before synthesis
provider/forward-predictor work; it does not precede the integrated upstream
design generator or a workflow's reference preparation. If all designs are
rejected before catalog loading, the run records `assets_not_loaded` and makes no
claim that the pin was checked.

## Request qualified recommendations

Use `--recommend` on either CLI or `from chembounce_synthesis import recommend`
for the version-fixed `qualified-selective-engineering-v7` contract. The
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

The local provider parses model responses as strict JSON by default. A new run
can explicitly opt into `smiles-bond-backslash-v1` by supplying a **JSON** prompt
template at `prompt_path` with the same `instruction` and `input` strings as its
chosen prompt plus this metadata key:

```json
"response_json_dialect": "smiles-bond-backslash-v1"
```

This policy interprets a literal backslash before a supported SMILES bond
following atom only inside `reactions[*].product`,
`reactions[*].reactants[*]`, or `building_blocks[*]` when the whole model
response is otherwise invalid JSON. It preserves the raw response and records
the exact character and UTF-8 byte offsets, derived JSON, hashes, and policy
source in `proposals.json`. Ambiguous escapes, unrelated syntax errors, and
unrecognized bond following tokens remain explicit failures. Valid JSON keeps
its existing meaning. Saved replay rederives the interpretation and rejects a
missing or changed audit. The prompt hash changes when this metadata is added;
new proposal acceptance under this policy needs separate evaluation and does
not change any chemical, route, or qualification gate.

The opt-in decoder checks the 250,000-character input and 128-level nesting
limits before strict parsing. It records a rejected raw response when nesting
exceeds the limit. For oversized generated text it retains at most 250,001
characters, records the original character count and any omitted count, and
marks the response rejected without scanning or hashing the oversized text.
Direct decoder calls reject oversized input before traversal. Rejected and
deadline-exceeded rows have no route-parsable `response`; their original text
remains in `raw_response` within the capture bound. Completed opt-in rows
require a timing seal on saved replay.
When a sample has already reached its deadline before decoding, replay keeps
that bounded, unparsed text as terminal evidence without claiming a nesting
check the live decoder never performed.
Sampling elapsed time includes generation, decoding, and audit work through the
final timing seal. Its audit hashes fixed sampling fields and records the final
elapsed time and deadline flag separately, avoiding a self-referential hash.
Saved replay checks those values and the original provider metadata against
the raw response. It accepts at most 32 consistent provenance levels and
requires unique sample indices within each design's declared sample count.
These checks establish internal consistency, not authentication of external
clock readings or model provenance.

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

Selected API rows and `results.json` include `hit_identity_comparison` (schema
version 1). Its reference is always `manifest.hit`, and its product is the selected
row's exact `smiles`, even after design or route fallback. A `recorded` comparison
contains:

- `reference_role: "original_hit"`, `original_hit_smiles`, `final_product_smiles`;
- `original_hit` and `final_product`: retained identities with `exact_smiles`,
  `family_smiles`, `scaffold_family_smiles`, `ring_skeleton_smiles`, `has_scaffold`,
  `policy_version` and `rdkit_version`;
- `change_class`: the existing classifier's actual result, with precedence
  `unchanged`, then `state_only`, then `decoration`, then `scaffold`;
- `atom_aware_scaffold_changed`: inequality of the normalized scaffold strings;
- `generic_ring_linker_topology_changed`: inequality of the generic Murcko graph
  strings, including ring atoms and their connecting linkers;
- `scaffold_presence_change`: `both_present`, `both_absent`, `gained` or `lost`.

The generic graph flag is descriptive; it adds no acceptance requirement. Both
acyclic identities have empty scaffold strings and `both_absent`; different
acyclic families are decoration changes, subject to the unchanged/state-only
precedence. Ordinary accepted hop rows therefore record `scaffold`, while an
atom substitution can leave the generic graph unchanged. These fields do not
measure experimental feasibility, activity or synthesis yield.

The producer retains the hit identity once next to an existing successful hop
constraint check at candidate admission, and retains each scored exact product's
identity at the existing family-field calculation. Root-search callbacks do not
perform the optional capture. The product map is bounded by the scored-product
allocation, and only final selected comparisons are serialized. After selection,
export uses retained values only, even if identity caches have been evicted.
Exploratory expansion does not normalize the original hit just for this export;
carried routes receive current-run admission evidence.

When retention, exact bindings, classification or policy agreement is missing or
conflicting, the comparison is atomically `unavailable`. It retains the schema,
reference role and exact SMILES fields, plus `reason_code`, `reason` and optionally
`error_type`; it omits both identities, classification, graph flags and presence
label. Missing legacy fields mean **not exported**, not a negative comparison.
An unavailable comparison can be displayed without a usable policy descriptor;
a recorded comparison requires the saved descriptor and exact bindings. A
manifest descriptor with zero selected rows is valid.

`manifest.identity_comparison_policy` records `schema_version: 1`, `policy_version`,
the actual `rdkit_version`, `tautomer_limit`, `transform_limit`, and `code_sha256`
for `identity.py` and `common.py`, captured at workflow import. Those hashes match
the corresponding saved manifest source snapshot. Offline reports must use the
saved identities and descriptor without reconstructing chemistry in another
runtime. They must not substitute installed package-version spelling for the
recorded RDKit runtime version.

`summary.identity_comparison_export` records retained/selected counts, unavailable
reasons, optional hit-capture lookup count and elapsed time, observed identity-cache
miss delta, late completion, and comparison serialization time. The once-only
optional lookup checks the existing deadline before starting and is not retried;
a cache miss can still normalize and finish after that deadline. Its cache delta
is a process observation interval, not attribution across concurrent threads or a
wall-time guarantee. `post_selection_identity_resolution_calls: 0` states the
export implementation invariant, not a measured process-wide call count. Added
overhead can affect availability in wall-limited runs. Mandatory chemistry errors,
selection/ranking rules, scores and TSV columns retain their existing behavior.

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

The optional `evidence_allocation_policy` controls only independent evidence
ordering, before 3D screening. Its default, `"inherit-v1"`, preserves existing
balanced or encounter ordering. Set `"design-scaffold-round-robin-v1"` to try
rotation across designs, then across scaffold families within each design. Both
levels use first encounter order from the actual admitted entries. Carried-route
priority remains first; a priority tier containing one design keeps its exact
inherited order. The separate screening sequence is never reordered by this option.

Each product/design entry still tries one route before its alternatives, with
shorter routes first. A multistep route can consume several fresh predictor calls;
cache hits consume none. This policy promises neither equal calls per design nor
equal runtime. The configured evaluator cache and product-matching policy are unchanged. Strict
`computational-route-support-v1` shares prepared-input/predictor-bound raw predictions
within the run and checks exact-product agreement independently for every route.
The legacy evaluator retains its existing matching and cache policies. Budgets, deadlines, unknown outcomes and final qualification
are unchanged. The option applies to both 2D and 3D workflows and is an explicit
`custom_unvalidated` material override for v1–v5 and v7 profiles. V6 explicitly
declares this policy together with ranked screening (see below). Summary and
manifest record it separately from screening and route allocation. It is a coverage
baseline, not a claim of improved yield or default superiority.

The optional `screening_allocation_policy` separates 3D screening order from
`allocation_policy`. Its default, `"inherit-v1"`, preserves the existing balanced
or encounter order and carried-route priority exactly. To try the experimental
original-hit similarity baseline, set:

```json
{"screening_allocation_policy": "original-hit-tanimoto-distinct-first-v1"}
```

This option reorders only evidence-passing product/design groups **after** the
unchanged independent-evidence queue. It ranks by original-hit Morgan Tanimoto,
descending, with ties resolved by first eligible route encounter and then exact
product. One group per distinct exact product is visited first; remaining design
alternatives follow in the same deterministic rank order. Every eligible route
association remains recorded. Each attempted design is scored and validated
independently; unattempted alternatives remain unscored. Carried priority still
applies to evidence allocation, but the experimental screening order is global.
It does not change refinement, acceptance floors, or final route-support checks.

`max_scored_products` retains its existing product/design screening-attempt unit,
including preflight failures; it is not a distinct-molecule or equal-time quota.
Summary counters distinguish attempts, distinct products, deferred design attempts
and actual scorer invocations. The summary and manifest name the screening policy
separately from `allocation_policy`. This option requires 3D scoring and is a
material `custom_unvalidated` override for v1–v5 and v7 profiles.
V2–V5 artifact bytes are unchanged. Greater 2D similarity can reduce scaffold
diversity; this baseline has no default-superiority or efficacy claim. V6 below
explicitly declares it together with evidence rotation; fresh-generation and
unseen evaluation remain separate work.

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

### Recorded forward outcomes and operational causes

Strict computational support retains the existing `forward_prediction.status` and
adds a `reason_code`, a per-step `forward_outcome`, and a route-level outcome.
A completed bounded top-three nonmatch differs from incomplete/invalid hypotheses,
unassessed work, and execution or validation errors. Budget causes distinguish
predictor calls, input tokens and output tokens; deadline, unavailable predictor,
and preparation errors retain their own codes and original error text. When
several budgets are exhausted, the first failing check in the existing
call/input/output order supplies the cause. No extra model work is performed.

These fields describe the available computation. They do not describe physical
synthesis success, failure or infeasibility. A route with different step outcomes
is marked `mixed` and retains every step cause. Exact completed product support
can still pass when other beams are incomplete; the existing strict agreement
rule, fresh raw-evidence validation and acceptance Boolean are unchanged.
Unsupported routes keep `selection_status: evidence_rejected` and remain
unqualified. Their qualification reasons now preserve the original plural route
reasons rather than collapsing to that selection status.

`summary.json` includes `forward_evidence_outcomes`, displayed in the existing
review diagnostics. Route-record and step-record totals have separate denominators;
alternative routes repeat records, so these counts are neither unique products
nor fresh predictor calls. Actual call/token/cache usage remains under
`synthesis_evidence_usage`. Missing historical cause codes remain `not_recorded`;
the reader does not infer them from prose. A legacy `mismatched` status without
a cause code is `nonmatch_recorded`, not proof of a completed strict top-three
evaluation. Reports and collectors reject newly
present diagnostic fields that contradict their saved statuses or record counts,
without treating diagnostic labels as proof of chemical support. These checks
cover every retained route, including unselected routes and empty-result runs,
and any optional forward outcome copied into qualification decisions.

This diagnostic correction does not change historical profile artifacts, score
minima, budgets or ranking. Historical v6 pins the workflow source, so a changed
current workflow is honestly reported as a source mismatch; the frozen campaign
and its source snapshot remain unchanged.

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
- `results.json`: selected API rows, including the retained original-hit identity
  comparison described above; the TSV schema is unchanged.
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

When strict forward hypotheses are recorded, their completed exact product and
lowest matching rank must agree with the step's saved match claim. Duplicate ranks,
contradictory raw/reported hypothesis lists and malformed completion flags are
diagnostic. This compares saved strings only; it does not decode tokens, normalize
SMILES or rerun a predictor. If only the raw record retains hypotheses, the report
validates and displays that list with its source labelled. Explicit malformed raw
records or hypothesis lists remain diagnostic. Only when neither location contains
hypothesis details does the report display that the match is a saved claim only,
rather than inventing missing evidence.

The offline report also reads selected `hit_identity_comparison` records and the
saved manifest `identity_comparison_policy`. It displays the producer-recorded
change class, atom-aware normalized scaffold change, generic ring/linker graph
change, and scaffold presence separately. An atom-aware scaffold hop can retain
the same generic graph. These labels do not add a similarity or synthesis test.

Old bundles without either field show `not_exported`. Atomic `unavailable`
records show their saved reason, including when the saved descriptor is missing
or unusable; this alone does not negate recorded qualification. Explicit null,
contradictory recorded values, or a declared descriptor with a missing selected
comparison produce artifact diagnostics. An empty selected array does not require
comparison records. The report checks saved strings, flags, labels and descriptor
hash copies without importing chemistry code, parsing SMILES, reconstructing a
class, or consulting the installed runtime/source. Matching saved copies do not
authenticate scientific provenance. Existing metric, seed, route and count checks
run independently, even if identity metadata is malformed. Each seed still has
its own common-pose witness; the comparison does not assert a shared pose across
seeds.

Nonabsolute enhanced stereochemistry is unsupported. OR groups describe unknown
absolute configurations; AND groups describe mixtures. The synthesis and core
input boundaries reject these groups (including mixed groups) before ordinary
SMILES serialization can discard them. The program does not choose an enantiomer
or infer experimental selectivity. Ordinary defined stereochemistry and explicit
ABS groups retain their existing absolute-state meaning. Direct raw-molecule
route diagnostics also report unsupported enhanced groups rather than calling
them fully specified.

Catalog conversion checks raw SDF molecules before serialization. In `.smi`,
`.smiles` and `.txt` files a complete CXSMILES annotation belongs to the molecular
input; an identifier may follow its closing `|`. Unsupported molecular records
retain rejection ordinals/reasons while later valid records survive. They are
never converted into identifiers or supplier claims. Previously serialized
ordinary SMILES cannot reveal an enhanced group already lost upstream; this fix
does not recover or certify old catalogs.

This safety change updates profile-bound source files without changing score
floors or packaged profile artifacts. Existing source-bound profiles therefore
report `custom_unvalidated` with explicit source mismatches. Any later engineering
profile rebinding requires separate review; the calibrated development profile
has not been silently rebound.


Exact-state text accepts ordinary SMILES or one complete `|a:...|` CX block whose
unique atom references are retained as specified tetrahedral ABS centers. Other
CX annotations, including coordinates/labels, axial `wU`/`wD`, extra blocks and
trailing names are explicitly unsupported formats; they are not declarations
that the molecule is chemically invalid. Optional identifiers in catalog SMILES
files are parsed separately. CX annotations on reaction SMARTS are unsupported;
ordinary stereochemical SMARTS remain available.

Explicit atropisomeric bond states are also rejected before serialization.
Older SDF readers may lose axial interpretation while retaining a bond wedge;
wedging not accounted for by a specified tetrahedral origin is rejected as
unsupported/unresolved, without inferring an axial configuration. This conservative
rule can reject decorative or ambiguous wedges. It cannot recover stereo metadata
already discarded by an upstream reader or an old serialized catalog, and it does
not classify potential atropisomerism in an unannotated structure.

### Engineering profile revisions and explicit migration

The v3 release introduced `qualified-selective-engineering-v3` as the then-current
recommendation default (superseded by v4 below). V3 has exactly the v2 floors, eight-conformer
screening, final seeds 42/43/44, budgets and joint-pose constraint policy. It binds
the reviewed exact-state input guards to new artifact bytes; it is not a new fit
or an empirically calibrated similarity/activity policy. Ordinary defined stereo
and supported tetrahedral ABS remain supported; unsupported representations are
rejected before lossy serialization (see the input-format contract above).

The original `qualified-development-v1` and `qualified-selective-engineering-v2`
artifacts and names remain available. On guarded source their historical kernel
bindings differ, so their provenance honestly reports `custom_unvalidated`.
Explicit historical config/profile requests are never silently migrated. Combining
`--recommend` with an explicit historical profile in either the original config
or CLI is an error before model/catalog work, even when the other declaration is
the current default. V3 is now also historical for this conflict check. Without the shortcut, `--recommendation-profile qualified-selective-engineering-v2`
continues to request v2.

For a continuation, inherited v2 without the shortcut remains v2. To migrate,
supply an explicit v3 replacement config, or use the profile-only
`--recommendation-profile qualified-selective-engineering-v3` option without
`--recommend`. Inherited v2 plus `--recommend` still conflicts before a v3 CLI
override can conceal it. Prior provenance remains retained, and carried routes
still undergo the existing fresh qualification gates; migration does not reuse
old scores or reset budgets/deadlines.

V3's `provenance.source_revision` names its preserved v2 predecessor, the reviewed
guarded commit, and two prior engineering review records. Release review verifies
those record bytes; the installed loader checks the portable packaged declaration
and pinned artifact digest without opening local review files. These records do
not certify calibration, biological activity, experimental synthesis, route yield
or a fresh comparison. `engineering_policy_matched` compares the five recorded
numerical modules, scientific runtime and pharmacophore feature content; it does
not certify every current pipeline source file. Supported material overrides
remain `custom_unvalidated`; resource limits and additional selection retain their
separate provenance fields.

Historical benchmark protocol-v2 and the saved-response development assessor stay
v2-only. Current campaign and prior capacity outputs must not be relabelled v3.
Any future v3 performance comparison needs an explicitly reviewed benchmark
protocol/adapter revision; this profile release launches no evaluation.

### Bounded reconstruction-index lineage migration

`chembounce_synthesis.catalog_migration` converts **one pinned processed index**
into a new catalog and independently validates complete record lineage. It does
not change the live catalog, switch a configuration, recover supplier IDs, or
approve adoption. A reviewed, externally SHA256-pinned JSON plan is required;
there is no automatic input discovery, retry or budget expansion.

The plan contains `schema_version: 1`, a resolved `input` path/hash/byte count,
`source` checkout/commit/tree and exact dependency hashes, `runtime` interpreter
path/bytes/version plus RDKit/NumPy versions and distribution-metadata digests,
a pinned GNU `watchdog` executable,
a resolved `parent` path/device/inode, an absent `output_dir` directly beneath
that parent, exact converter `labels`, a frozen legacy `baseline`, and `limits`.
The source dependency list is `__init__.py`, `assets.py`, `catalog.py`, `common.py`,
`reactions.py` and `catalog_migration.py`. Dirty or mismatched sources fail closed.
The baseline includes output SHA256, unique-row count and exact rejected-record
coordinates/reasons. It is supplied evidence, never repaired from new output.

The Linux coordinator creates the output directory exclusively with mode 0700,
then launches four sequential children: `startup`, `converter`, `validator` and
`bookkeeping`. Every child is restricted to one allowed CPU, an empty CUDA device
list, offline flags and one thread per numerical library, with `RLIMIT_AS` applied
before interpreter startup. The chemistry processes release their separate memory
graphs when they exit. The initial reviewed ceilings are:

| Limit | Ceiling |
| --- | ---: |
| Entire external invocation | 1,800 seconds |
| Shared normal-work deadline | 1,740 seconds |
| Startup / converter / validator / bookkeeping phase | 120 / 600 / 900 / 120 seconds |
| Owned TERM grace / KILL and reap / final receipt reserve | 15 / 15 / 30 seconds |
| Per-child address space | 16 GiB |

All phase deadlines are clamped to the same work deadline. Smaller fixture limits
are permitted; larger limits are rejected. An address-space cap is not a promise
that the full inert object graph, origin map or eager `Catalog` molecules fit.
The exact reviewed invocation must also have an **external watchdog**, for example:

```bash
/usr/bin/timeout --signal=TERM --kill-after=60s 1740s \
  /absolute/pinned/python -m chembounce_synthesis.catalog_migration run \
  --plan /absolute/reviewed-plan.json --plan-sha256 REVIEWED_PLAN_SHA256
```

The exact executable/hash, environment, working directory and source/runtime pins
must be frozen for the actual attempt. The external 1,740-second TERM and 60-second
KILL grace cover interpreter startup and final writes, with a total 1,800-second
ceiling. Each child additionally has an autonomous watchdog bounded by remaining
normal work, so interruption of the coordinator does not leave an unguarded child.
The coordinator uses Linux subreaping and cleans only its owned process groups.
It records failure when a group remains alive or cannot be reaped; no universal
termination guarantee is made for uninterruptible native/kernel states.
After safe owned cleanup it restores the process's previous subreaper flag,
including when a phase fails or terminal persistence later fails. Unsafe cleanup
retains subreaper ownership; a failed restoration cannot return success. The
terminal snapshot records the prior/observed flag and restoration status. This
process-wide flag is not isolation from concurrent unrelated child activity;
the real invocation runs in its own coordinator process. No foreign children are
reaped or signaled to repair process state.

The converter is invoked once, unchanged. The separate validator rejects extra
pickle payloads or trailing bytes and replaces every pickle global with an inert
placeholder. It independently checks each raw SMILES through the existing pinned
canonicalization policy, regenerates its structure ID, and compares the exact
TSV structure and all v2 manifest fields. Accepted occurrences **A** and rejected
coordinates **R** must partition all **N** logical source positions; **U** is the
separate unique TSV count. Duplicate accepted occurrences remain origins, not
rejections. Unsupported pickle/molecule schema makes the source incomplete.
Ordinary canonicalization failures retain exact coordinates and typed reasons.
Every loaded ID's `Catalog.evidence` is checked, including source hashes and both
SKU/stock verification flags being false, without fingerprints or similarity search.

Input/source bindings are checked before and after phases. Output checksums are
bound across conversion, validation and final bookkeeping. The whole output must
match the frozen legacy TSV baseline; any difference remains a failed attempt
for review. Do not change normalization or omit rejections to force equality.
`owner.json`, phase logs/receipts, failure receipts and `terminal.json` preserve
partial output, nonzero exits, timeouts and unavailable cost counters. Records
are written exclusively and fsynced. Missing counters are unavailable, not zero;
child self-usage and pre-receipt coordinator elapsed time have explicit scopes.
The terminal receipt is a snapshot taken before its own serialization and fsync.
It is never rewritten to certify the completion time of that same write. A
successful migration requires a `validated` snapshot, coordinator exit 0 after
the write strictly before the original work and total deadlines without an
interruption, and external watchdog exit 0 within the complete ceiling. A late
write can therefore leave a `validated` snapshot with a nonzero exit; the snapshot
alone is insufficient. A killed writer or unavailable terminal receipt never
proves completion. Independent result review and adoption remain separate gates.
Existing and dangling/symlink destinations are refused.

Processed-index lineage cannot certify original supplier-state completeness.
The retained upstream preprocessing converts supplier molecules to ordinary
SMILES before serialization and may deduplicate using nonisomeric canonical
SMILES. Original enhanced-stereo annotations or separate supplier stereo records
cannot be recovered or certified from that processed index. This is a source-
processing limitation, not a finding that particular input records contain such
states. A change to ingress/state policy requires a new source pin and baseline
review before migration; newly rejected states must remain visible. Obtaining
original supplier exports and authenticating SKU/current-stock claims are
separate dependencies.

### Offline structure and saved-pose review

New successful synthesis runs optionally export `review_visuals.json`. Use the
existing offline report command, `python -m chembounce_synthesis.report RUN_DIR`,
to see actual original-hit and final-product 2D PNG depictions and an interactive
saved-pose overlay. The report itself remains standard-library-only. It does not
parse chemistry, load models, align molecules, or calculate new scores. Old
bundles remain readable with `not_exported` display status; they are not enriched.

The seed selector defaults to the recorded representative seed. Each option uses
that seed's own saved hit reference and aligned candidate, with one shared camera.
Rotate by dragging/touch or arrow keys, zoom with the wheel, buttons or +/−, pan
with shift-drag, and reset with 0 or Reset camera. An unavailable representative is
not replaced automatically. Across-seed minima in the evidence tables remain
separate from the selected seed's metrics. The 3D lines show connectivity; they
are not a stereochemical or experimental selectivity certificate.

The producer verifies exact isomeric structure identity around 2D preparation,
exports indexed graph metadata, and reuses saved coordinates without embedding,
optimization, alignment, feature scoring or route replay. The reader checks
saved row/trial/pose bindings and the actual indexed topology projection. That
projection binds elements, isotopes, charges and bond connectivity/types, **not
stereo**. Stereo display labels have separate integrity bindings; PNG semantics
remain the producer's responsibility. Fully rewritten self-consistent files are
not independently authenticated by hashes. Display consistency is reported
separately from the workflow-recorded recommendation and route support.

Display is optional and cannot veto a recommendation. Its whole cooperative
allowance is **5 seconds**, further capped by the original workflow deadline,
even when no workflow deadline was supplied. Limits are 32 selected candidates,
32 retained seeds per candidate, 512 atoms and 1,024 bonds per graph, 1 MiB per
static RGB/RGBA PNG, 1,024 pixels per dimension and 1,048,576 total pixels, and
16 MiB serialized display JSON. Coordinate magnitudes are bounded at 1,000,000;
the shared camera clamps zoom to 0.25–8 and uses a finite shared scale for
coincident coordinates. Limits restrict optional display, never scientific
selection. Missing Cairo/PNG capability, missing saved geometry or other display
failures remain explicit unavailable evidence. Published artifacts must cover
every selected row and retained seed; no partial available artifact is accepted.

Core results are written and flushed before optional work. The manifest starts
with an unavailable display descriptor; only complete optional output is
advertised through an atomic manifest extension. Known deadline exhaustion
prevents starting further work or a commit. Native work and an already-started
filesystem commit cannot be interrupted by this cooperative mechanism. A valid
commit that returns late is content-available, with
`publication_completion_within_budget: null` and
`timing_scope: final_persistence_unverified`, never certified as within budget.
`generation_elapsed_seconds` excludes publication; the workflow logger
`chembounce_synthesis.workflow` reports observed post-return total display cost,
and an actual overrun also emits a warning. Logging is best effort, not a durable
completion receipt. Existing scientific summary timing is not relabeled to
include this separate display phase. A failed metadata write can leave only the
original unavailable descriptor or an unadvertised orphan, which the report
ignores. No recursive completion receipt, rollback or scientific retry is used.

The HTML contains only bounded data PNGs and fixed bundled viewer code with a
restrictive content security policy; it fetches no external resources. Legacy
reports without display evidence remain script-free. Conditions, activity,
experimental feasibility, supplier SKU and current stock remain unverified
unless separate evidence explicitly records them.

### Optional rejected-pose diagnostics (source-review candidate)

Set `"rejected_pose_diagnostics": true` in a selective synthesis configuration to
retain a bounded diagnostic for future runs that finish with
`no_jointly_feasible_hit_pose`. The default is `false`; the option must be a JSON
boolean and requires `similarity_constraints`. This does not rescue or qualify a
candidate, change a floor, or add an embedding, alignment, feature match, route
search, or model call. Opting in adds overhead and can affect work completed under
a fixed wall deadline.

For each eligible failed scoring call, the collector chooses the highest existing
scalar among completed rejected pairs, keeping the first pair on a tie. It saves
that pair's exact original-hit reference coordinates and aligned candidate
coordinates, identities, screening/refinement phase, seed/settings, both relative
energies, scores/reasons, and the already-computed polar/landmark witness. Infeasible
joint assignments remain infeasible or unassessed. The **global pharmacophore
assignment is not retained**: its existing score/counts cannot identify a unique
feature correspondence. This is one diagnostic pose, not an exhaustive
conformational explanation or biological/synthesis evidence. Failed overlays,
preflight failures, deadline abstentions, successful selections and later design
failures do not acquire substitute rejected poses.

One shared collector spans screening and refinement: at most one record per call,
eight per workflow, 512 atoms in each pose, 256 KiB per record and 2 MiB for the
**complete** `rejected_pose_diagnostics.json` including envelope/provenance. Excess
records are skipped, never truncated into incomplete witnesses or retried. The
one-second cumulative cooperative allowance includes import/setup handoffs,
retention bookkeeping, failed captures, hashing, serialization and publication;
checks are also capped by the existing workflow deadline. Already-entered native
or filesystem operations can overrun. The artifact's final publication timing is
explicitly unverified; the summary reports observed costs and overruns separately.
The updated workflow elapsed time includes diagnostic sidecar publication when
enabled; the final summary write endpoint remains explicitly unverified. A
best-effort log records post-return publication and summary-update overhead.

Existing failure audits in `routes.json`, `scores.json` and `refinement.json`
carry only `captured_not_published` references (index and canonical record SHA),
so repeated routes do not duplicate complete coordinates. Resolve a reference
only against this run's successful `summary.json` diagnostic export receipt,
matching the **whole-file SHA** and then the record SHA. An unavailable export
leaves references unresolved, including when a stale sidecar happens to exist.
The core summary is saved first with an unavailable/pending descriptor. One
budgeted atomic status update may then resolve publication; if it fails or the
allowance expires, the initial summary survives and references stay unresolved.
Ordinary optional persistence errors preserve the original rejection and core
outputs; interruption/termination is not swallowed. The offline recommendation
report does not present these rejected records as selected evidence. Private
serialized records yield defensive copies; source files/runtime labels are bound
once per collector, without model/asset discovery or per-record package scans.

The historical source-only diagnostic stages (`40324`/`612596`) changed the
hash-bound selective scorer while retaining the then-current v3 default. At that
stage v2/v3 correctly reported `custom_unvalidated`; it was not a deployable
matched-profile update. The v4 integration below adds a separately reviewed
declaration and changes the default together, while preserving the old artifacts.
Nothing reconstructs diagnostics for older runs that did not retain these poses.


### Engineering v4 source revision (historical)

The v4 release introduced `qualified-selective-engineering-v4` as its
`--recommend`/`recommend()` default; v5 and then v7 superseded it. It preserves the exact v3 numerical
specification, runtime declaration, floors, seeds, budgets and strict route
support. Its new artifact binds the accepted rejected-pose diagnostic scorer
revision and retains the prior exact-state guard provenance. Source review is
not calibration, efficacy, yield or a fresh scientific evaluation.

`rejected_pose_diagnostics` defaults to `false`. When enabled, its separate
profile provenance records version, declared limits and the diagnostic/scorer/
refinement/workflow source bindings. A mismatch in those enabled bindings makes
the profile custom. It does not change numerical floors or rescue a rejection.
Profile-declaration/source hashing is setup work outside collector statistics;
it consumes the enclosing workflow deadline when supplied. Collector costs, skips
and overruns are in `summary.rejected_pose_diagnostics`;
final write completion requires the existing post-return observation. Diagnostic
work consumes the existing deadline and can change deadline-limited coverage.

Historical v1/v2/v3 artifact bytes, names and benchmark authorities remain
unchanged. Explicit v3 still selects v3 and honestly reports a scorer source
mismatch on this checkout; `--recommend` conflicts with an explicit v3 rather
than silently upgrading it. Continue old runs under their saved declaration.
A new local `qualified-selective-engineering-v4-lineage-v1.json` configuration
must be separately reviewed and installed after release; the installed v3
configuration is not modified by this source change. Generic runs remain
exploratory unless qualified behavior is explicitly requested.

### Experimental landmark-aware rigid proposals

`similarity_constraints.landmark_pose_proposals: true` explicitly enables
`landmark-anchored-rigid-v1` for **nonempty user-supplied landmarks**. Omitting the
option or setting it to `false` retains the existing policy and alignment path;
it does not require the fixed-pose `ScoreMol` API. Enabling it without landmarks,
with a non-boolean value, or without that API fails during configuration preflight.
The initial prototype changed selective kernel hashes: historical v2–v4
profiles remain unchanged and are honestly `custom_unvalidated`. V5 below binds
the reviewed source with this option still disabled; enabling the option is
not a matched policy, calibration, or default recommendation method.

Only a completed baseline alignment with an infeasible mandatory-landmark
assignment triggers proposals. The helper enumerates family-correct, one-to-one
mandatory correspondences, then applies a proper rigid transformation to **all**
original atoms, including isotopic hydrogen. A single anchor proposes translation
and, when a distinct typed helper feature exists, a deterministic shortest
rotation. Multi-anchor fits correct the determinant to prohibit reflection.
Internal geometry, conformer IDs and stereo are preserved. Proposal
correspondences are never presented as the final polar/landmark assignment.

Fixed-pose shape/color scoring uses the same heavy-atom view as `_overlay`, with
normalization disabled. Its baseline shape/color must agree with the actual
original alignment within `1e-5`, otherwise the proposal path is explicitly
unavailable. Every alternative gets the unchanged full pharmacophore overlap and
joint polar/landmark assessment on that exact pose. Selection considers only
feasible poses with the existing scalar, **both** relative-energy penalties and
first-on-tie behavior; a lower-scoring feasible pose can replace a rejected
baseline. The ordinary design diagnostic and all final-seed/route gates remain.

The fixed limits are two alternatives per pair, eight per scoring call, sixteen
`ScoreMol` calls including parity checks, 128 mandatory correspondence edge visits
and four complete maps per pair, and 4096 helper-anchor edge visits per pair.
Added work has a one-second **cooperative** allowance per scoring call, capped by
the existing deadline. Native calls may return late; late proposals are discarded
and costs/overruns retained. Search truncation is distinct from geometric
infeasibility. No budget is transferred to conformer, generation or route search.
Opting in consumes time and can change deadline-limited outcomes.

`pose_selection.landmark_pose_proposals` records separate pair/proposal counts,
parity, mapping truncation, transforms and timing. Each pair preserves its baseline
result; its top-level decision describes the actually selected pose. Selected
`hit.landmark_proposal` binds the exact baseline and transformed coordinates,
reference pose, helper policy/source, and corresponding audit entry. Saved-seed
validation checks this association without embedding, aligning or searching.
Optional rejected-pose diagnostics retain the actual proposal origin when that
pose supplies the diagnostic. The workflow manifest also includes the helper's
source hash. Config/policy preparation and final saved-result validation are
setup/verification costs outside the per-scoring-call proposal timer and still
consume any enclosing workflow deadline.

The motivating manual fixture demonstrated an omission only for an explicitly
requested 0.25 Å donor landmark. Default fixtures already passed. This feature is
not evidence of improved ordinary recommendation coverage, recovery of historical
failed controls, biological activity, or experimental synthesis feasibility.


### Engineering v5 and optional landmark pose proposals

V5 introduced `qualified-selective-engineering-v5` as the recommendation default;
the current v7 transition is described below. Its declared engineering floors remain
original-hit Tanimoto/shape/pharmacophore **0.35/0.70/0.50**, polar recall **0.75**,
eight conformers at seed 42 and final seeds 42/43/44. Route support and all final
qualification checks are unchanged. V5 binds the reviewed scorer and the bounded
landmark helper, and inherits the prior guarded-stereo and diagnostic source
review declarations. These are portable software-review records authenticated by
the packaged artifact and release review, not calibration, route-performance,
activity or experimental synthesis certification. Loading does not reopen the
historical record paths; the declarations contain no current self-hash cycle.

Landmark pose proposals remain **off by default**. Explicit configuration:

```json
"similarity_constraints": {
  "version": "selective-similarity-constraints-v2",
  "polar_recall": 0.75,
  "landmarks": [
    {"id": "requested-acceptor", "family": "Acceptor", "hit_atom_ids": [2], "max_distance": 0.25}
  ],
  "landmark_pose_proposals": true
}
```

Atom IDs refer to the exact canonical original hit; the example index must be
checked against that hit. `--landmark-pose-proposals` overrides this boolean on
either CLI, and `--no-landmark-pose-proposals` explicitly disables a carried
setting. Enabling without nonempty landmarks or fixed-pose `ScoreMol` fails
preflight. Disabled requests require no `ScoreMol` capability and do not import
the helper for profile resolution. The option is a material search override:
**enabled requests remain `custom_unvalidated`**, even with a v5 source match.
Landmarks alone remain additional acceptance constraints. They are not evidence
of biological essentiality.

Only already generated conformers are used, with at most two rigid alternatives
per pair, eight per scoring call, and 16 fixed-pose score calls including parity
checks. Typed injective correspondences propose proper rotations/translations;
the original same-pose evaluator makes the final decision. Mapping and helper
search caps, partial returns and a one-second cooperative added-work allowance
are recorded and capped by the existing absolute deadline. Native work can
return late; exhausted search is truncation, not proof of geometric
infeasibility. Each selected proposal retains its exact pose, transform and
witness association through screening, refinement and saved-result validation.
Proposal correspondences are distinct from final feature assignments.

Historical v1–v4 artifacts and the benchmark's v2 authority remain unchanged.
Explicit historical requests are never silently upgraded; combining one with
`--recommend` is a conflict. Their source mismatches remain visible. Continuation
retains its explicit profile and option unless replaced by an explicit user
setting, and revalidates/rescores the carried routes. The offline report displays
the saved profile/custom status, policy and seed witnesses; it performs no new
alignment or scientific validation. This integration's controlled fixtures are
software plumbing tests, not recommendation-yield evidence; no historical
campaign is relabelled as v5 evaluation.

### Explicit v6 allocation method

`qualified-selective-engineering-v6` is an opt-in engineering profile. It combines
`screening_allocation_policy: original-hit-tanimoto-distinct-first-v1` with
`evidence_allocation_policy: design-scaffold-round-robin-v1`. Select it explicitly
with `--recommendation-profile qualified-selective-engineering-v6` (alongside the
usual synthesis configuration), or set `recommendation_profile` in that config.
Do not combine this choice with `--recommend`: `--recommend` and `recommend()`
request default v7 and reject a conflicting explicit profile.

V6 retains the .35/.70/.50 original-hit floors, .75 polar recall, exact-product
conditional route support, scaffold-change checks and three-seed qualification.
Screening prioritizes distinct exact products by original-hit Tanimoto, then
retains deferred design alternatives without copying their scores or poses.
Evidence rotation changes request ordering within the existing carried tiers;
cache sharing, per-route multistep costs, unknown status and deadlines remain
unchanged. Product/design screening attempts and predictor caps are not promises
of equal end-to-end computation.

On the two inspected development inputs with saved generator proposals, the
combination retained two qualified products on the first input and zero on the
second, identical to ranked screening alone. It broadened bounded evidence
coverage but did **not** demonstrate better qualified yield. This is not a
calibrated method, a default recommendation change, or unseen-input performance.
Subsequent fresh native/SynLlama generation and full workflow integration on the
same two seen inputs returned two/zero qualified products and passed independent
saved-result review. This is a separate integration result, not a new comparative
or held-out evaluation. See [current evidence and next work](QUALIFIED_RECOMMENDATIONS.md#explicit-v6-and-current-evidence).

V6 binds `workflow.py` as primary method source in addition to the six v5
scoring/identity/landmark modules. Its optional diagnostics and landmark proposal
semantics are retained; enabling landmark proposals is still a material custom
override. Any allocation change from the v6 declaration remains
`custom_unvalidated`. Old profile artifacts and their source scopes are unchanged.
A matching profile is not full-source or asset certification.

### Engineering v7: v5 policy with optimized validation

`--recommend` and `recommend()` explicitly select
`qualified-selective-engineering-v7`. V7 preserves all v5 effective settings:
the .35/.70/.50 original-hit floors, .75 polar recall, exact regenerated-product
route support, normalized scaffold change, seeds, budgets, search and ranking.
Both screening and evidence allocation retain `inherit-v1`; v6 allocation is
not selected. Landmark proposals and rejected-pose diagnostics remain off.

V7 changes the primary `common.py` source binding to the reviewed indexed graph
traversal implementation. The remaining five primary source pins, runtime and
optional feature declarations are inherited from v5. The retained optional
diagnostic declaration still reports its historical `workflow.py` mismatch;
enabling those diagnostics therefore remains `custom_unvalidated`.

The v1–v6 artifact bytes, names and exact source checks are unchanged. Their
`common.py` mismatch remains visible on this implementation. Combining an
explicit historical profile, including v5 or v6, with either recommendation
shortcut rejects before work starts; use an explicit profile without the
shortcut to retain that historical request. Continuation retains its saved
profile unless explicitly replaced. No historical benchmark protocol follows
the new default or accepts v7 implicitly.

One sequential full-catalog pair observed lower constructor CPU time with
matching complete molecular/catalog/lineage digests and 60 sampled fingerprints.
This is not a robust or end-to-end speedup result, a changed similarity contract,
or evidence of improved qualified yield. See the
[measurement and compatibility note](CATALOG_LOAD_OPTIMIZATION.md).
