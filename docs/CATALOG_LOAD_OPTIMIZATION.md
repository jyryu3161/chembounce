# Catalog loading traversal candidate

This candidate is **not release-ready**. It changes `common.py`, which is bound
by the existing numerical profile source hashes. The packaged v1–v6 profiles,
including default v5, correctly report `custom_unvalidated` on this source.
Their artifacts, defaults and matching rules remain unchanged. Release requires
independent review and a separate prospective profile compatibility decision.
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

Parity checks compare valid and rejected inputs, exact exception types/messages,
canonical and noncanonical serialization, complete RDKit molecule binaries with
properties, ABS/explicit-H handling, maps, fingerprints, ordered catalog records,
lineage evidence and nearest results. The focused regression stage passed 263
tests. It also checks mutable `Mol`/`RWMol` guard calls, validation before map
writes, and rejection before later atoms or bonds are accessed.

Local development evidence is under
`benchmark-results/qualified-development-20261016/catalog-load-optimization-20261018-author-v1/`:
`baseline.json`, `baseline-profile.txt`, `synthetic-comparison.json`, and
`profile-compatibility.json`. No additional full-catalog construction has been
performed for this candidate.
