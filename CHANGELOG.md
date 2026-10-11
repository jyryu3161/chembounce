# Changes

## 2.0.0.dev0

- Integrate the 2.0 development work and 1.x maintenance history into `main`.
  Retain the six backport regressions with current database and error-reporting
  behavior; use general examples in README and Colab.
- Replace the aromatic-atom removal/recombination code with explicit boundary
  bonds and stereo-aware molzip assembly. Bound intermediate search with a
  deterministic beam instead of factorial remnant permutations.
- Reject provably over-capacity neutral C/N/O replacements before enumeration;
  preserve normal assembly for charged scaffolds, other elements and unusual bonds.
- Fix core preservation checks, ignored global limits, off-by-one fragment limits,
  shared replacement-file indexing, stale global database caching, disabled text
  fallback, and zero-fragment division.
- Report zero-candidate outcomes, property exclusions, 3D failures and actual
  effective settings; fail visibly on invalid inputs and unexpected errors.
- Use Morgan bits consistently for both scaffold search paths. Preserve Morgan
  sparse counts for final candidate similarity. Stop using SMILES string length
  as a molecular size proxy and stop treating fingerprint equality as identity.
- Add packed, chunked top-k search, a memory-mapped database format, deterministic
  tie ordering, canonical deduplication, and executable database verification.
- Seed ETKDG/UFF conformers, cache target descriptors, calculate shape only for
  accepted property-filtered candidates, and remove per-scaffold garbage collection.
- Fix CATS feature matching (all alternatives), canonical pair indexing and
  occupancy scaling; use RDKit shortest-path distances.
- Remove import-time environment-variable requirements, global warning suppression,
  and mandatory MolVS, ODDT, PubChemPy dependencies. Standardize with RDKit Cleanup.
- Add installation metadata, a tested dependency lock, CI, regression tests,
  benchmark artifacts, safer data installation, and a revised Colab notebook.
- Add the original publication citation, a versioned PubChem drug panel, isolated
  1.0/2.0 attachment comparisons, beam ablations, negative stereo-check controls,
  and per-input evaluation time budgets that retain failures in reported results.

These changes deliberately alter results; 2.0 is not a drop-in reproduction of
1.0 rankings, CATS values, 3D scores, or candidate enumeration.

## 1.0.1

Maintenance fixes on `release/1.x`, tagged `v1.0.1`:

- Backport boundary-bond scaffold assembly and preserve retained stereochemistry.
- Preserve input stereo, remove the mandatory `RDBASE` variable, and explain
  property-filter rejection and assembly errors in `run_summary.json`.
- Update the legacy environment to Python 3.10 / RDKit 2023.09–2025 and point
  the Colab notebook at the patched tag. Keep the legacy descriptor dependencies.
- Add six ADP regression tests, including real CLI output and property rejection.

This does not include the broader 2.0 search, descriptor and global-limit changes.
The original benchmark artifacts retain their historical labels and source hashes;
their baseline is now explicitly `v1.0.0`.

## 1.0.0

Legacy implementation from `eaadfa7`, with only a README version notice added in
`5780007`. Preserved unchanged under annotated tag `v1.0.0`.

## Version-label migration

Releases use immutable `vMAJOR.MINOR.PATCH` tags. The original `1.0-original` and
`1.0` names remain as compatibility aliases for `v1.0.0` and `v1.0.1`, respectively;
new documentation and reproduction scripts use the semantic tags. The numeric
`2.0` branch is integrated into `main` and retired. Python package metadata uses
`2.0.0.dev0` until a versioned prerelease or stable release is cut.
