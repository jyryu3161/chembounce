# CReM train-only fragment baseline

This is a thin adapter around the official published CReM implementation, with
separate preparation and validation of a split-specific train-only fragment DB.
It is a proposed post-hoc follow-up (90 attempts: 30 public queries × 3 seeds),
not a change to the original 630-run study, a preregistration, a matched-training-
data comparison, or a completed performance result. No full training DB or
benchmark query generation was run during implementation. Small synthetic
engineering fixtures are not activity/scaffold-recovery evidence.

The observed-before-generation settings and primary source hashes are in
[crem-proposal.json](crem-proposal.json). CReM is fixed to upstream commit
[`b84003cd8c84bc2afe98b5473570b7328e703ead`](https://github.com/DrrDom/crem/tree/b84003cd8c84bc2afe98b5473570b7328e703ead),
which includes the generated-product stereo fix. Code is BSD-3-Clause; source
data retain their own attribution/license obligations. This bundle contains no
private molecules, large database, models, or upstream source redistribution.
See the [CReM paper](https://doi.org/10.1186/s13321-020-00431-w).

## Fixed method and counting

One official `mutate_mol` call uses `partial_exo`, radius 3, core sizes 1–15,
relative core sizes 0–1, size increment −2…+2, frequency ≥1, set `train`, one
core, and `max_replacements=256`. Seeds are 42/43/44. The official builder uses
`both_optimal`, radius 3, maximum core size 15, `keep_stereo=True`, and one core.
No iterative mutation, expansion, retry, beam search, or target-label selection
is added. Radius/size-increment defaults retain local context; `partial_exo`
adds the documented narrower ring-arc edit path, compatible with `both_optimal`.
See [official API, lines 1893–2035](https://github.com/DrrDom/crem/blob/b84003cd8c84bc2afe98b5473570b7328e703ead/crem/crem.py#L1893)
and [DB compatibility](https://github.com/DrrDom/crem/blob/b84003cd8c84bc2afe98b5473570b7328e703ead/docs/fragment-databases/build-v2.md).

The shared cap is 360 seconds including imports, hashing, DB initialization,
mutation, normalization and ranking; process-tree RSS cap is 24,000 MB; at most
three products are returned. CPU threads are one; GPU visibility is empty.
The declared count unit is **API-yielded products**, capped at 256. Upstream also
caps sampled replacements and suppresses invalid/duplicate/geometry-rejected
products internally. Internal fragment matching, rejected products and SQL work
are **not** counted as 256 comparable objective evaluations. Their unavailable
counts remain null, not zero, and their work is bounded by the outer wall cap.
API yields, postfilter reasons, progress, exceptions and budget stops are retained.
No partial candidates are accepted after an internal or outer timeout/error.
The outer runner may discard a failed response but retains raw/progress files.

Official `discard_ring_geometry=True` is kept. Its rigid aromatic system helper
uses ETKDGv3 with fixed seed 42 and `maxIterations=4`, followed by a 15° torsion
scan for one hinge; it is separate from mutation seed 42/43/44. Geometry errors
fail open, so this is not a 3D validity guarantee. No upstream patch or new
force-field optimization is introduced. Same-seed checks demonstrate only the
tested environment/fixture behavior; no cross-RDKit or universal determinism claim.
[Embedding and geometry handling](https://github.com/DrrDom/crem/blob/b84003cd8c84bc2afe98b5473570b7328e703ead/crem/crem.py#L1018)

## Structures, partition roles and limits

Only connected valid ring-containing products are eligible. The common versioned
family identity excludes the query family and unchanged normalized Murcko
scaffold. Ranking is descending **unfolded count Morgan Tanimoto**, radius 2,
`includeChirality=False`; ties use exact canonical isomeric SMILES. The best
representative per family is chosen before the output cap. There is no additional
similarity threshold. Scores are neither activity nor synthesis probabilities.

Atom maps are removed as identifiers; exact charge, specified stereochemistry
and isotope labels are retained. Bounded fixed-point canonicalization fails
explicitly on cycles/nonconvergence; family normalization failures are recorded.
There is no protonation/tautomer/stereo enumeration of scoring structures.
Unspecified potential stereo is recorded, **not excluded from 2D generation or
selection**. It can fail the common later 3D scorer and must remain in coverage
denominators. No arbitrary organic-element or total-heavy-atom cap is imposed.
The upstream ordinary fragmenter has `maxCutBonds=30`; excess matching bonds are
recorded, not treated as a blanket rejection of the partial-ring branch.

The official DB builder strips isotope labels. Isotopically labelled **training**
records are explicitly excluded with their original exact/family identities;
they are not silently rewritten. By contrast the official generator backs up
and restores isotope labels on retained query atoms, so isotopic queries/products
are supported. This distinction is directly visible in
[query property handling](https://github.com/DrrDom/crem/blob/b84003cd8c84bc2afe98b5473570b7328e703ead/crem/crem.py#L29)
and [builder normalization](https://github.com/DrrDom/crem/blob/b84003cd8c84bc2afe98b5473570b7328e703ead/crem/scripts/cremdb_create.py#L422).

`prepare` verifies original source hashes through the existing dataset validator,
TSV/provenance/dataset hashes, exact train membership, and absence of every
nontrain partition family (test, validation and query partitions). Public queries
in the **train partition remain permitted** in the DB; their IDs, families and
hashes are explicitly counted in the private preparation inventory. The existing
scaffold/time panels use public train queries; target queries have query roles.
No broad “query contamination absent” or pretrained clean-holdout claim follows.
The adapter never opens the panel, labelled references or preparation inventory.

Known fragments can recombine into new full Murcko scaffolds. This does not
guarantee invention of unseen ring systems. Radius/core-size/context/DB coverage
limit reachability. Ordinary acyclic cuts can already replace a small whole ring
attached by an acyclic bond; `partial_exo` additionally enables ring-arc changes.
`forced`, which bypasses cyclic core size limits, is not used. DB preparation cost
is measured separately: hashes, inventories/exclusions, official builder counters,
DB tables/ring rows/bytes, CPU/wall and sampled process-tree peak RSS. Sampled
resource measures are labelled as such. The public builder does not expose a
per-molecule zero-fragment inventory; this remains unknown rather than zero.

## Environment and preparation

Use a new CPU venv and exact official checkout; do not alter campaign environments.
The tested interpreter is Python 3.10.12 / RDKit 2025.09.6. The committed lock
records the resolved fixture environment, including test tools. `prepare`/`sample`
verify import origins, isolated-venv configuration, dependency metadata, tracked
upstream file hashes, generated version source and project implementation hashes.
Discovered metadata can include setuptools-vendored packages; it is not a count
of independently installed distributions. No tool installs or downloads anything.

```bash
python3.10 -m venv /path/to/new-crem-venv
/path/to/new-crem-venv/bin/python -m pip install -r benchmarks/external/requirements-crem.txt
git clone https://github.com/DrrDom/crem.git /path/to/crem-upstream
git -C /path/to/crem-upstream checkout --detach b84003cd8c84bc2afe98b5473570b7328e703ead
/path/to/new-crem-venv/bin/python -m pip install --no-deps --no-build-isolation -e /path/to/crem-upstream

# Set hash seed before starting Python. Apply the same environment to both
# preparation/build and the benchmark runner launching the adapter.
export CUDA_VISIBLE_DEVICES=''
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED=0

nice -n 19 /path/to/new-crem-venv/bin/python benchmarks/external/crem_baseline.py prepare \
  --dataset /path/to/panels/direct-human-scaffold/panel.json \
  --train-tsv /path/to/train-direct-human-scaffold.tsv \
  --provenance /path/to/train-direct-human-scaffold.tsv.provenance.json \
  --upstream /path/to/crem-upstream --output /path/to/fresh-preparation
```

Repeat preparation separately for `time` and `target`. Review and pin each
`plan.json` before building. The normalized inventory is private; hashes and
counts can be published separately without distributing chemical lists. Existing
output paths are refused. After independent review, the bounded **offline** build
uses a separate output directory (the default 3600-second build cap is not the
360-second per-query generation cap):

```bash
nice -n 19 /path/to/new-crem-venv/bin/python benchmarks/external/crem_baseline.py build \
  --plan /path/to/fresh-preparation/plan.json --output /path/to/fresh-build \
  --wall-seconds 3600 --max-memory-mb 24000
```

A failed/incomplete build retains diagnostic artifacts and never emits a usable
`database-manifest.json`. Reuse/append/resume is deliberately unsupported: the
official builder's additive mode cannot silently merge prior data. The successful
public manifest contains DB/code/runtime hashes and counts, not private panel paths,
labels or molecule inventories. `build-report.json` and preparation files stay
separate from sampler assets. A full build has not yet been authorized/executed
as part of this implementation.

## Existing runner integration

Use `kind: external`; do not modify the runner. Its command should be an absolute
path to the dedicated interpreter and this exact helper, followed by `sample
--request {request} --response {response}`. Declare only these two assets:

```json
{
  "id": "crem_train_partial_exo",
  "kind": "external",
  "command": ["/path/to/new-crem-venv/bin/python", "/path/to/checkout/benchmarks/external/crem_baseline.py",
              "sample", "--request", "{request}", "--response", "{response}"],
  "assets": {
    "crem_db": {"artifact": "/path/to/fresh-build/fragments.db", "sha256": "ACTUAL_DB_SHA256"},
    "crem_manifest": {"artifact": "/path/to/fresh-build/database-manifest.json", "sha256": "ACTUAL_MANIFEST_SHA256"}
  },
  "settings": {"policy": "crem-train-only-partial-exo-v1", "manifest_sha256": "ACTUAL_MANIFEST_SHA256"},
  "budget_overrides": {"max_evaluations": 256}
}
```

This is a schema example, not an executable protocol. Use the parent-approved
actual hashes, seeds, matched wall/memory/output caps, and source provenance.
`public_request` can include a public query pose; the adapter deliberately does
not open it, use target IDs to select a DB subset, or accept added labels/reference
fields. It uses one DB per split, never a query-specific DB. The sampler requires
the runner's monotonic deadline; standalone untimed sampling is not supported.

`response.raw.jsonl` retains every yielded string; `response.progress.json` keeps
the latest accounting if killed. The final response records decisions, unknown
internal rejection counts and `no_route_supplied`. Unsupported structures produce
zero candidates with the reason retained. Errors/timeouts exit unsuccessfully and
accept none of the partial products. Common evaluator failure denominators remain
unchanged. The existing evaluator may label unverified evidence while retaining
the declared `no_route_supplied`; zero verified routes denotes coverage, **not** a
chemical synthesis failure. CReM `return_rxn` is disabled because an edit is not
a validated synthesis route.

After generation, the evaluator can compute the same family/scaffold diversity
and structural active-reference recovery definitions, and a separate common
`similarity_scores` pass with the same conformers/seeds/pose policy used for other
methods. Preserve all 3D failures, including unspecified stereo, in denominators.
Basic `evaluate_candidates` does not automatically compute 3D. Do not use held-out
activity labels or 3D results to retune this generator. Aggregate by target with
nested seeds/queries; the three panels do not create 30 independent targets.
Synthetic feasibility comparison would require the same separately declared
retrosynthesis planner for every method's final products. No planner is added here.

## Engineering checks

The dedicated CPU fixture run passed **44 tests in 9.68 seconds**. Without the
explicit upstream opt-in, 34 passed and 10 optional integrations skipped (0.99 s).
[Validation record](crem-engineering-validation.json) binds source versions,
file/log hashes and scope; independent parent assessment is pending.

```bash
CHEMBOUNCE_CREM_SOURCE=/path/to/crem-upstream \
  /path/to/new-crem-venv/bin/python -m pytest -q tests/test_crem_baseline.py
```

Tests cover real official DB/API roundtrip, ring-row compatibility, extra
partial-ring reachability, same-seed fresh-process order and output, immutable DB
bytes, direct count-fingerprint ranking, state/stereo preservation, family
deduplication, public train-query roles, nontrain-family and exact-membership
guards, hashes/source/manifest tampering, public-input schema, output protection,
count/deadline limits and retention of failed partial raw outputs. Tests requiring
the optional upstream package explicitly skip unless `CHEMBOUNCE_CREM_SOURCE` is
set; generic guard tests still run in a compatible RDKit environment. These small
fixtures do not measure full-DB build cost or benchmark-method performance.
