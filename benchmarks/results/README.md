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
- `test-matrix.json` and original test logs: five tested dependency combinations
  through `619de8a`. `test-matrix-capacity.json` and `tests-capacity-*.txt` record
  the application code in `d7d4d61` (81 core tests, plus an optional ODDT comparison).
  The ODDT
  oracle test is skipped when the optional legacy library is absent.
- `coverage-final.txt`: application-module coverage through `619de8a`;
  CLI subprocess tests execute but are not included in this in-process coverage.
- `library_preparation.json`, `reference_download.json`: actual input hashes,
  preparation counts and verified download digest. Both legacy and generated DBs
  also have 1,000-row recomputation verification logs.
- `wheel-build.txt`: successful installation artifact build.
- `wheel-capacity.txt`: successful isolated wheel build including capacity pruning.
- `pubchem/`: raw structure-property responses for the fifteen publication drugs;
  the versioned panel in `../publication_molecules.json` records retrieval URLs,
  CIDs, hashes and the tirzepatide name-resolution decision.
- `attachment-beam-ablation.json` and `attachment-beam-summary.json`: 256 controlled
  calls comparing tagged 1.0 and pre-capacity 2.0 at three beam widths. Full input,
  fragment, replacement and candidate lists are retained, including empty/error/
  timeout outcomes. These are attachment tests without workflow filtering.
- `attachment-final.json` and `attachment-final-summary.json`: 64 calls repeating
  beam 50 after capacity pruning at clean commit `d7d4d61`. All mutually completed
  candidate lists match the earlier beam-50 run exactly; 560 entries passed checks.
- `capacity-comparison.json`: the same impossible tirzepatide replacement at
  `0f14a8f` and `d7d4d61`, with a 120-second budget, in the same RDKit environment.
- `publication-workflow/`: all sixteen inputs on the smaller cleaned library,
  clean commit `d7d4d61`, and a 60-second wall budget per input. Eleven successes,
  one completed zero-result input, and four timeouts are retained. Per-input
  `evaluation_status.json` is authoritative for completion: files in timed-out
  directories are incomplete evidence, not finalized output. Their candidate
  counts are unknown and excluded from molecular statistics. Root protocol and
  summary files retain their place in the evaluation denominator.

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
