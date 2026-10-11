# Recorded validation artifacts

These are observed local results, not claims of biological activity or a
statistically representative chemical benchmark. See [validation](../../docs/VALIDATION.md)
for commands, conditions, and remaining manuscript work.

- `full/`: final full-library panel executed at clean commit `619de8a`. Each input
  has its manifest, diagnostics, tables and resource metrics. `summary.tsv` reports
  runtime/counts; `candidate_statistics.tsv` reports descriptive molecular metrics.
- `mw250/`: corresponding smaller-library panel at clean commit `b5040e9`, before
  the empty-vector compatibility patch; the generated molecular fingerprints are
  nonempty. Its source hashes identify the exact implementation.
- `adp_shape/`: ADP with ten accepted candidates, seeded 3D scoring and filters off.
- `adp_default_filters/`: ADP with 43 property-rejected candidates and explicit
  `no_candidates` status. This demonstrates the filter issue independently of the
  repaired attachment code.
- `text_lowmem/`: successful streaming search of the public MW≤250 text library.
- `test-matrix.json` and test logs: five tested dependency combinations. The ODDT
  oracle test is skipped when the optional legacy library is absent.
- `coverage-final.txt`: application-module coverage on the final implementation;
  CLI subprocess tests execute but are not included in this in-process coverage.
- `library_preparation.json`, `reference_download.json`: actual input hashes,
  preparation counts and verified download digest. Both legacy and generated DBs
  also have 1,000-row recomputation verification logs.
- `wheel-build.txt`: successful installation artifact build.

Individual manifests preserve the paths used on the test machine, full parameters,
seed, software versions, code hashes, Git state and library checksums. Summary-table
manifest links are relative to their containing directory for portability.
Large databases are excluded from Git; retrieve the original archive and regenerate
using the recorded RDKit and fingerprint settings. Later generator versions may
add metadata fields without changing fingerprint arrays. Recompute and inspect
actual file hashes before combining results from separately prepared libraries.

The stress panel uses `--no-shape`, explicitly disables property filters, and caps
results at 20 overall / 10 per fragment, with at most 20 scaffolds per fragment and
a beam width of 50. Refer to each manifest for all actual settings. Candidate counts
are bounded observations, not an enumeration of all possible molecules. Comparing
the cleaned and legacy libraries directly confounds representation, deduplication,
invalid-record handling, library composition and sometimes tie order; use the
isolated kernel benchmark for the recorded search-kernel speed comparison.
