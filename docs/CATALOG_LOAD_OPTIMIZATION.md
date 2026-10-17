# Catalog loading traversal candidate

The isolated traversal candidate was **not release-ready**. It changes `common.py`, which is bound
by the existing numerical profile source hashes. The packaged v1–v6 profiles,
including former default v5, correctly report `custom_unvalidated` on this source.
Their artifacts and matching rules remain unchanged. The prospective
`qualified-selective-engineering-v7` profile and explicit default transition
bind this implementation to v5's unchanged effective policy. That separate
release candidate requires independent review; implementation parity alone
does not authorize release. See the [v7 transition](SYNTHESIS.md#engineering-v7-v5-policy-with-optimized-validation).
The frozen `e97036e` v6 campaign is unaffected.

The change reduces Python dispatch during eager molecular validation. The stereo
guard visits bonds by index. Connected-molecule validation fetches each atom
once, checks for dummy atoms in order, and retains its handle until every atom
has passed; only then are atom maps cleared. Parsing, stereo checks, canonical
serialization, record lineage, byte binding and lazy fingerprints are unchanged.
The public stereo guard still reads the current graph on every call.

On CPU 2 with one thread and RDKit 2025.09.6, five alternating synthetic pairs
gave median CPU times of 0.5812 → 0.4919 seconds for 2,800 molecule validations
and 0.7670 → 0.6730 seconds for a 2,800-row synthetic catalog. These are small
synthetic diagnostics, not full-catalog speedup or interactive-latency claims.
The initial full-catalog baseline profile was 116.93 CPU / 116.98 wall seconds;
profiling overhead prevents comparison with those unprofiled synthetic timings.

One subsequent unprofiled full-catalog pair, in sequential fresh processes on the
same CPU, measured baseline `62e1bdf` at 94.3769 CPU / 94.4292 wall seconds and
candidate `84efc546` at 79.3623 CPU / 79.4008 wall seconds. This single observation
is 15.91% lower constructor CPU time; it does not establish a robust speedup.
Peak process RSS was 3,505,388 / 3,504,748 KiB. The separate post-construction
audits took 11.21 / 11.31 CPU seconds. All 238,542 ordered molecule keys and full
property binaries, ordered building blocks and structure indexes, manifest
lineage, and 60 predetermined direct fingerprint samples had identical digests.
Every catalog fingerprint remained pending and both search caches remained empty.

Parity checks compare valid and rejected inputs, exact exception types/messages,
canonical and noncanonical serialization, complete RDKit molecule binaries with
properties, ABS/explicit-H handling, maps, fingerprints, ordered catalog records,
lineage evidence and nearest results. The focused regression stage passed 263
tests. It also checks mutable `Mol`/`RWMol` guard calls, validation before map
writes, and rejection before later atoms or bonds are accessed.

Local development evidence is under
`benchmark-results/qualified-development-20261016/catalog-load-optimization-20261018-author-v1/`:
`baseline.json`, `baseline-profile.txt`, `synthetic-comparison.json`, and
`profile-compatibility.json`. `full-pair-v1/comparison.json` records the subsequent
pair and hashes its helper, exact invocations, per-process outputs and reports.
Each process independently verified clean source heads and both pinned asset
hashes before and after its load. This was exactly one pair, with no retries.
