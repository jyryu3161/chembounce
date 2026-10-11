# Changes

## 2.0.0.dev0

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

## 1.0

Legacy implementation from `eaadfa7`, with only a README version notice added in
`5780007`. Annotated tag `1.0`; preserved maintenance branch `release/1.x`.
