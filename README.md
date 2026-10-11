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

The optional [synthesis workflow](docs/SYNTHESIS.md) connects SynLlama proposals
to a dated building-block catalog and an explicit reaction registry. It independently
replays routes, expands compatible building blocks, reapplies scaffold/core constraints,
and ranks regenerated molecules by joint multi-conformer 3D shape/pharmacophore
similarity to the original hit and the proposed design. Run the offline fixture with
`python -m chembounce_synthesis.cli`; use `.[synthesis,synllama]` and external model/catalog
assets for actual SynLlama inference. Verified reaction replay is computational evidence,
not an experimental synthesis guarantee.

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
limits, and RDKit version. Filters are never relaxed automatically to produce rows.

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

Use `--fingerprint-db data/fingerprints_v2` when running ChemBounce. A directory
stores packed bits in memory-mappable NPY files; a `.npz` output is compressed but
must be decompressed into memory. Packed bits use one eighth of the legacy uint8
bit-matrix storage. `--scaffold-db FILE` selects a text library and never silently
uses unrelated bundled fingerprints. `--no-fingerprint-db` forces text search;
`--low_mem` streams a text library per fragment. Unusable text-library records
are excluded with counts and examples in `run_summary.json`.

Old NPZ files can contain Python pickle objects. Regenerate from SMILES where
possible. To read a trusted legacy file explicitly, use
`--fingerprint-db PATH --trust-legacy-db`. Never enable this for an untrusted file.
Automatic database selection prefers `data/fingerprints_mw250_v2` when present;
otherwise the legacy filename is `data/scaffold_fingerprints_mw250.npz`.

## Options

Run `python chembounce.py --help` for the full interface. Key settings:

| Option | Meaning |
| --- | --- |
| `-i`, `-o` | Input SMILES and output directory |
| `-t` | Minimum final Morgan **count** Tanimoto; default 0.5 |
| `--search-threshold` | Minimum scaffold Morgan **bit** Tanimoto; default 0.3 |
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
- `run_manifest.json`: inputs, effective limits, seeds, dependency versions, source
  code hashes, Git revision/dirty state, and actual library checksums.
- `resource_cost.json`: CLI wall time, process CPU time, and sampled process RSS.

The legacy result column names are retained, including `Electron shape Similarity`.
A shape failure is recorded as NaN and explained in the summary; it does not discard
an otherwise valid candidate. New 3D scores use a seeded conformer and can differ
from 1.0. `Standardized final structure` uses RDKit Cleanup, without tautomer or
protonation-state enumeration. Details and limitations are in [methods](docs/METHODS.md).

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

Regression tests cover scaffold attachments, retained stereochemistry, identity
reconstruction, core preservation, limits, diagnostics, deterministic 3D scores,
CATS descriptors, and packed-search equivalence to RDKit. Clean environments on
Python 3.10, 3.11 and 3.12 each passed 81 tests; a separate environment also passed
the optional ODDT comparison (82 total). The minimum
RDKit 2023.09 dependency combination passed as well. After integration, 88 tests
passed with the optional comparison installed; the minimum environment passed 87
and skipped that optional check. CI mirrors these checks.

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
