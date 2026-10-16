# ChemBounce 2.0 (development)

ChemBounce generates scaffold-replacement candidates for a single molecular
SMILES input. Candidates are computational hypotheses; similarity and property
scores do not establish biological activity or synthetic feasibility.

Development takes place on **`develop`**; the current Python package version is
**`2.0.0.dev0`**. The **`main`** branch contains the stable 1.x version, with
maintenance fixes carried on **`release/1.x`**.
Release tags follow [Semantic Versioning](https://semver.org/):

- **`v1.0.0`**: original legacy baseline and environment (`5780007`).
- **`v1.0.1`**: scaffold-assembly fixes and updated legacy environment.

Use `main` for stable use, `develop` to evaluate the latest improvements, and
`v1.0.0` with the original library for legacy reproduction. There is no stable 2.0 release yet. Published `vX.Y.Z` tags
are immutable; subsequent fixes get a new patch version. Short-lived `feat/...`
branches merge into `develop`; stable fixes use `fix/...` branches and are
incorporated into both the stable and development lines. Python development versions use the
[PEP 440 format](https://packaging.python.org/en/latest/specifications/version-specifiers/).
See [changes](CHANGELOG.md) and [methods](docs/METHODS.md).

For **qualified recommendations**, use the explicit `--recommend` shortcut:

```bash
chembounce -i "$HIT" -o NEW_RUN --recommend --synthesis-config my-strict-assets.json
chembounce-synthesis --hit "$HIT" --design "$DESIGN" --recommend \
  --config my-strict-assets.json --output NEW_SYNTHESIS_RUN
```

Start from the [strict asset template](examples/synthesis/strict-recommendation-config.json)
and configure existing local catalog/rules, proposal provider and pinned ReactionT5
assets. Saved proposals must cover the supplied designs; for newly generated designs,
configure a local SynLlama provider as explained in [setup](docs/SYNTHESIS.md#request-qualified-recommendations).
No assets are downloaded implicitly. The shortcut selects
`qualified-selective-engineering-v2`: final original-hit Tanimoto/shape/pharmacophore
floors **0.35/0.70/0.50**, polar recall **0.75**, changed normalized scaffold,
complete three-seed geometry and concrete per-step computational route support.
Supported explicit overrides survive and are labelled `custom_unvalidated`; the
summary prints the actual thresholds. These engineering constraints do not guarantee
activity or experimental synthesis. An empty final result is bounded abstention,
with no exploratory fallback. Generic calls without a recommendation request retain
the existing exploratory behavior. See the [Python API and conflict rules](docs/SYNTHESIS.md#request-qualified-recommendations).

The optional [synthesis workflow](docs/SYNTHESIS.md) connects SynLlama proposals
to a dated building-block catalog and an explicit reaction registry. It
replays routes, expands compatible building blocks, reapplies scaffold/core constraints,
and ranks regenerated molecules using Morgan-bit Tanimoto and joint multi-conformer
3D shape/pharmacophore similarity to both the original hit and the proposed design.
Use `--synllama --synthesis-config FILE` to enable inference, or `--no-synllama`
to run hopping alone. Expansion neighbour counts and rounds are configurable;
no minimum synthesis-stage similarity is imposed unless configured. Run the offline fixture with
`python -m chembounce_synthesis.cli`; use `.[synthesis,synllama]` and external model/catalog
assets for actual SynLlama inference. Verified reaction replay is computational evidence,
not an experimental synthesis guarantee.

The explicit `--protocol research-v1` preset reserves bounded search capacity for
diverse, lower-2D-similarity scaffolds and balances synthesis work across scaffold
families. It checks three refinement seeds, increases conformer effort for unstable
or borderline candidates, and ranks using conservative seed summaries while
exporting one actual representative pose. Similarity cutoffs remain opt-in;
stability tolerances are engineering criteria, not activity probabilities.
Exact reaction products remain separate from chemical families used to avoid
counting tautomer, protonation, isotope or stereo variants as novel scaffolds.
Optional independent ReactionT5 evidence and bound-pose receptor/anchor checks
have explicit configuration and provenance; they are not enabled merely by
selecting the preset. See [research setup](docs/SYNTHESIS.md#research-preset-and-bounded-2d-exploration)
and the [research validation protocol](docs/RESEARCH_VALIDATION.md).

Use `chembounce-synthesis --continue-from PREVIOUS_RUN --output NEW_RUN` to
reuse saved model proposals and retain previous routes while expanding the search.
See [refinement and continuation](docs/SYNTHESIS.md#refinement-and-continuation)
for adjustable budgets and output interpretation.

## Install

Python 3.10 or newer:

```bash
git clone --branch develop https://github.com/jyryu3161/chembounce.git
cd chembounce
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

Alternatively, `conda env create -f environment.yml` creates a Python 3.11
environment and installs this checkout. `environment_simple.yml` is an alias for
this same environment. `requirements-lock.txt` records the tested clean Linux /
Python 3.10 environment; install it with `pip install -r requirements-lock.txt`
before `pip install -e . --no-deps` for the recorded dependency versions.

[Google Colab notebook](https://colab.research.google.com/github/jyryu3161/chembounce/blob/develop/chembounce_colab.ipynb)
uses a small demo library by default and bundles diagnostics with the download.

## Quick start

Run a small example using the included scaffold library:

```bash
python chembounce.py \
  -i 'CC(=O)NCCc1ccccc1' \
  -o output_example \
  --scaffold-db tests/fixtures/scaffolds.smi \
  --overall_max_n 20
```

The included library verifies operation without downloading external data; it is
**not the publication library**. Results depend on library, thresholds, search
limits, and RDKit version. The default `legacy` policy applies mandatory 2D
cutoffs. To explicitly select bounded diversity exploration without synthesis:

```bash
chembounce -i 'CC(=O)NCCc1ccccc1' -o output_research_hopping \
  --scaffold-db tests/fixtures/scaffolds.smi --overall_max_n 20 \
  --protocol research-v1 --no-synllama
```

In this preset, 2D cutoffs prioritize one search lane; a reserved diversity lane
can admit candidates below them. Core/property constraints remain mandatory.
The selected policy, effective budgets and exclusions are recorded, so a change
of search policy is explicit rather than an automatic response to an empty run.
For integrated hopping, model projection and adaptive refinement, follow the
[local-asset example](docs/SYNTHESIS.md#research-preset-and-bounded-2d-exploration).

## Reference data and fast search

```bash
bash install.sh
python generate_fingerprints.py \
  --input data/Scaffolds_processed.txt \
  --output data/fingerprints_v2 --verify
```

The installer downloads [Zenodo record 16741967](https://zenodo.org/records/16741967)
and checks its published archive digest before extraction. It refuses to overwrite
an existing data directory. Fingerprint generation canonicalizes and deduplicates
connected scaffolds, reports rejected lines, and records source hash and RDKit
version. Existing output databases are also protected from overwrite.
Generation checks that each serialized isomeric SMILES is unchanged after
reparsing, with cycle detection and at most eight roundtrips. Fingerprint bits
come from that exact reparsed molecule. Unresolved inputs are excluded as
`canonicalization_unstable`; metadata records cycle/limit counts and original
input examples. Stereo is not removed to force convergence. New assets declare
`canonicalization_policy=rdkit-roundtrip-fixed-point-v1`; only this policy with a
matching RDKit version enables the canonical-string fast path for selected plain
strings without CX delimiters or whitespace. Annotated/name-bearing selected rows
always use guarded normalization and refill, even if metadata claims canonical
storage. This prevents annotation bypass; it does not certify arbitrary forged
metadata for ordinary strings.

Fingerprint input uses the same exact-state stereo guards as synthesis: ordinary
SMILES and supported tetrahedral ABS-only CX are accepted; OR/AND, axial CX,
other CX metadata and trailing names are unsupported. Builders check text before
parsing; raw-Mol query/serialization helpers check retained stereo before use.
Unsupported rows count separately as `unsupported_molecular_representation`,
with original one-based line numbers and reasons in the first ten rejection
examples. All-invalid input publishes no database. Existing bits, roundtrip
limits and ranking budgets are unchanged. Already-erased upstream annotations
cannot be recovered from stored ordinary SMILES; sampled verification also checks
raw text but certifies only its selected sample.

Fingerprint bits and canonical SMILES can differ between RDKit releases. Search
rejects a database with a known different RDKit version before ranking. Build a
separate index with the runtime you will use; for RDKit 2025.09.6 and the mw250
library:

```bash
python generate_fingerprints.py \
  --input data/Scaffolds_processed_mw250.txt \
  --output data/fingerprints_mw250_rdkit_2025.09.6 --verify
# Explicit use; this also takes precedence over automatic database selection.
chembounce -i 'CC(=O)Nc1ccccc1' -o output_runtime_db \
  --fingerprint-db data/fingerprints_mw250_rdkit_2025.09.6
```

Use your actual RDKit version in the output directory name. `--verify` recomputes
sampled fingerprint bits; a passing sample does not establish whole-database
compatibility or enable cross-version ranking. Neither search nor verification
rewrites an existing asset. Missing-version legacy databases remain readable and
searchable, with fingerprint compatibility unverified. Assets without the new
canonicalization policy have selected hits checked for a stable identity; exact
query matches, duplicates and unresolved hits are excluded, and the shortlist is
refilled. An unresolved query fails explicitly. Specified stereo, isotope and
charge differences remain distinct.

Use `--fingerprint-db data/fingerprints_v2` when running ChemBounce. A directory
stores packed bits in memory-mappable NPY files; a `.npz` output is compressed but
must be decompressed into memory. Packed bits use one eighth of the legacy uint8
bit-matrix storage. `--scaffold-db FILE` selects a text library and never silently
uses unrelated bundled fingerprints. `--no-fingerprint-db` forces text search;
`--low_mem` streams a text library per fragment. Unusable text-library records
are excluded with counts and examples in `run_summary.json`.

Old NPZ files can contain Python pickle objects. Regenerate from SMILES where
possible. To read a trusted legacy file explicitly, use
`--fingerprint-db PATH --trust-legacy-db`. This authorizes pickle loading only; it
does not override RDKit compatibility checks. Never enable it for an untrusted
file. Automatic selection first prefers
`data/fingerprints_mw250_rdkit_<runtime-version>`, then `data/fingerprints_mw250_v2`,
then the legacy filename `data/scaffold_fingerprints_mw250.npz`. An explicit
`--fingerprint-db` always takes precedence.

## Options

Run `python chembounce.py --help` for the full interface. Key settings:

| Option | Meaning |
| --- | --- |
| `-i`, `-o` | Input SMILES and output directory |
| `--protocol legacy\|research-v1` | Named defaults for hopping and synthesis; explicit overrides take precedence |
| `--search-policy legacy\|diverse` | Mandatory 2D gates or bounded 2D-priority/diversity lanes |
| `-t` | Final Morgan **count** Tanimoto cutoff in legacy mode, priority boundary in diverse mode; default 0.5 |
| `--search-threshold` | Scaffold Morgan **bit** Tanimoto cutoff in legacy mode, priority boundary in diverse mode; default 0.3 |
| `--exploration-fraction` | Reserved diversity fraction in diverse mode; default 0.25; zero disables lower-2D admission |
| `--search-candidate-budget` | Global requested candidate slots in diverse mode, including unsuccessful assembly requests; default 4096 |
| `--core_smiles` | Substructure required in both input and every accepted candidate, including specified stereo |
| `--overall_max_n` | Hard cap on unique results across fragments |
| `-n`, `--frag_max_n` | Hard cap per fragment; default 1000 |
| `--scaffold_top_n` | Maximum replacement scaffolds tested per fragment |
| `--cand_max_n__rplc` | Maximum generated candidates per replacement; default 10 |
| `--beam-width` | Maximum intermediate attachment states; default 100 |
| `--fragment-depth` | Maximum ring-system removal depth; default 1 |
| `--fragments SMILES` | Explicit fragment(s); repeat for multiple fragments |
| `--replace_scaffold_files FILE` | SMILES lines or SMILES/priority TSV; one shared file or one per fragment |
| `--wo_lipinski` | Disable default MW ≤500, logP ≤5, donors ≤5, acceptors ≤10 |
| `--qed_min/max`, `--sa_min/max`, `--logp_min/max`, `--mw_min/max` | Explicit property bounds |
| `--h_donor_min/max`, `--h_acceptor_min/max` | Explicit donor/acceptor bounds |
| `--seed` | Deterministic conformer seed; default 42 |
| `--no-shape` | Skip the base 3D descriptor; the optional synthesis stage uses its own conformer ensemble |
| `--synthesis-config FILE` | Project hopping designs onto catalog reaction space and rank with joint 3D/pharmacophore scoring; see [setup](docs/SYNTHESIS.md) |
| `--synllama` / `--no-synllama` | Explicitly enable model inference with its config, or skip the synthesis stage; without either flag, supplying a config retains the existing opt-in behavior |
| `--expansion-top-k N` | Number of alternative building blocks per mutable reaction leaf |
| `--expansion-rounds N` | Building-block substitution rounds; 0 disables expansion |
| `--expansion-max-routes N` | Maximum retained expansion routes per starting route, shared across rounds |
| `--refine-top-n N` | Reevaluate up to N screened chemical families under the configured allocation policy; 0 disables refinement |
| `--refine-conformers N` / `--refine-seeds S ...` | Conformer count per seed and seeds for reevaluation |
| `--synthesis-max-scored N` / `--synthesis-max-results N` | Synthesis screening and final output budgets |
| `--include-run DIR` | Reverify and retain previous synthesis routes; repeat to combine runs |
| `--overwrite` | Explicitly replace prior run artifacts in an output directory |

Single connected molecules are required. Salts/mixtures, invalid SMILES, dummy
atoms, invalid cores, contradictory bounds, and nonpositive limits produce errors.
Acyclic inputs require explicit `--fragments` because Murcko scaffolds contain rings.
Explicit fragments must not cut through aromatic bonds within a fused ring system.
The first deterministic substructure match is replaced. The bounded beam is a
heuristic and does not enumerate all substitution isomers or all stereoisomers.

## Results and diagnostics

- `overall_result.txt`: unique accepted candidates in tab-separated columns.
- `fragment_info.tsv`, `replace_scaffold_list.fragment_XXXX.tsv`,
  `fragment_result_XXXX.txt`: fragment-level search evidence.
- `run_summary.json`: status, input filter violations, generation/acceptance counts,
  rejected properties, and optional 3D failures.
- `search_decisions.json`: diverse-policy candidate admissions, exclusions and
  bounded attachment/search decisions, including identity failures.
- `run_manifest.json`: inputs, effective limits, seeds, dependency versions, source
  code hashes, Git revision/dirty state, and actual library checksums.
- `resource_cost.json`: CLI wall time, process CPU time, and sampled process RSS.

The legacy result column names are retained, including `Electron shape Similarity`.
A shape failure is recorded as NaN and explained in the summary; it does not discard
an otherwise valid candidate. New 3D scores use a seeded conformer and can differ
from 1.0. `Standardized final structure` uses RDKit Cleanup, without tautomer or
protonation-state enumeration. Details and limitations are in [methods](docs/METHODS.md).

Synthesis outputs are separate under `OUTPUT/synthesis/`. `results.tsv` selects
one accepted exact product/design per chemical family; `screening_alternatives.tsv`
retains scored alternatives so a failed preferred design need not erase the family.
`alternative_route_ids` always refers to the selected exact product; different
state/design alternatives are labelled separately. Seed trials, adaptive budgets,
independent synthesis evidence and optional pocket geometry are preserved in the
[synthesis artifacts](docs/SYNTHESIS.md#artifacts-and-evaluation).

Exit codes: **0** candidates found; **3** valid run with no passing candidates;
**2** invalid request/data. Unexpected programming failures remain visible.
A zero-result run still writes a header and diagnostic report.

Python API:

```python
from chembounce import chembounce

results, costs = chembounce(
    "CC(=O)NCCc1ccccc1",
    output_dir="output_api",
    scaffold_db_file="tests/fixtures/scaffolds.smi",
    candidate_thresholds={},  # explicitly disable property filters
    calculate_shape=False,
)
```

API and CLI both apply Lipinski filters by default. To repeat an existing output
location, explicitly pass `overwrite=True` / `--overwrite`.

## Validation and performance

```bash
pytest
ruff check .
python benchmarks/benchmark_search.py --rows 100000 --repeats 5
```

Regression tests cover scaffold attachments, retained stereochemistry, chemical
families, core preservation, limits, diagnostics, deterministic 3D scores,
CATS descriptors, packed-search equivalence to RDKit, route alternatives,
multi-seed/adaptive refinement, independent reaction evidence and bound-pose
constraints. Historical environment checks and their scope are recorded in
[validation](docs/VALIDATION.md); optional 3D and model runtimes have additional
requirements described in [synthesis setup](docs/SYNTHESIS.md).

The recorded [100,000-row kernel benchmark](benchmarks/search_100k.json) found
identical top-100 scores, **2.24×** lower median search time, **8×** smaller fingerprint
arrays, and sampled process RSS of about **162 vs 528 MiB** on the recorded machine.
This is a synthetic search-kernel measurement, not an end-to-end speedup or a claim
of improved molecular quality. End-to-end runs on both cleaned public libraries
(522,959 and 3,267,123 unique scaffolds) succeeded on six ring-containing inputs; [results and provenance are included](benchmarks/results/README.md). See [validation](docs/VALIDATION.md) for evidence and
remaining publication work.

An additional [versioned panel](benchmarks/publication_molecules.json) contains the
15 drugs named in the original publication plus a reported regression case. Controlled attachment tests
compare the original 1.0 snapshot (`v1.0.0`) with 2.0 and vary beam width.
A separate 60-second-per-input workflow
evaluation on the smaller library returned candidates for 11 inputs, no candidates
for tirzepatide, and timed out on four complex inputs. All outcomes, including
incomplete runs, are [recorded](benchmarks/results/README.md). These engineering
checks do not establish retained activity or superiority over other design methods.

## Attribution and contact

Original publication: Jang et al., [ChemBounce: a computational framework for
scaffold hopping in drug discovery](https://doi.org/10.1093/bioinformatics/btaf501),
*Bioinformatics* 41(9), 2025. See [CITATION.cff](CITATION.cff). For 2.0 experiments,
also report the commit, library hashes and settings from the run manifest;
2.0 changes the implementation and does not reproduce every published result.

ChemBounce uses RDKit and ScaffoldGraph. CATS2D originated from
[iwatobipen/CATS2D](https://github.com/iwatobipen/CATS2D); local bug fixes and method
conventions are documented in [methods](docs/METHODS.md). ElectroShape follows
[Armstrong et al. (2010)](https://doi.org/10.1007/s10822-010-9374-0).
Contact: [CSB Lab](https://www.csb-lab.net/).
