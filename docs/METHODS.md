# Methods and reproducibility

This document describes the base scaffold-hopping workflow. The optional
[synthesis workflow](SYNTHESIS.md) independently verifies catalog reaction routes
and adds multi-conformer joint shape/pharmacophore ranking. Its scores and
verification status are separate from the base single-conformer descriptors.

## Scope and representation

ChemBounce proposes scaffold replacements. The input is sanitized as a single
RDKit molecule; isotopes and specified stereochemistry are retained, atom-map
annotations are cleared, and salts/mixtures are rejected. Canonical isomeric
SMILES represent identity. Tautomers and protonation states are not enumerated.
RDKit Cleanup provides the separately reported standardized structure; ranking
and deduplication use the generated structure before Cleanup.

## Fragmentation and assembly

ScaffoldGraph produces a Murcko scaffold and ring-system parents. Enumeration is
sorted by canonical SMILES and bounded by graph depth, with a visited set. The
1.0 code counted siblings instead of recursion depth. Acyclic inputs have no
automatic scaffold; users may supply explicit fragments.

For each fragment, the first matching atom-index tuple is selected. Every bond
crossing its boundary is cut with RDKit FragmentOnBonds. Outside components retain
isotope-labelled attachment dummies. Each proposed scaffold attachment is zipped
with its corresponding outside dummy, sanitized, scored and deduplicated.
RDKit molzip preserves stereo on retained atoms. Explicit fragments that cut an aromatic bond inside a fused ring system are
rejected with a request to select the whole aromatic system. Invalid valences are
expected attachment rejections. Unexpected exceptions propagate into the run error report.

A conservative capacity precheck skips enumeration when retained bond-order demand
exceeds the total available valence of a neutral C/N/O replacement, after allowing
all removable hydrogens to be substituted. The upper bounds are C=4, N=5 and O=2;
N=5 deliberately allows sanitization of hypervalent nitro-like forms. This check
applies only to retained C/N/O/H components and ordinary single/double/triple
boundary bonds. Charged replacements, other elements and unusual bonds fall back
to normal assembly. It rejects some provably impossible cases; passing the bound
does not establish that attachment is possible. See RDKit's
[sanitization and valence rules](https://www.rdkit.org/docs/RDKit_Book.html#molecular-sanitization).

At each attachment step only `beam_width` highest-scoring intermediate states are
kept. Scores use the entire partial molecule, including not-yet-connected outside
components and remaining dummies. Ties use canonical SMILES. This is an explicit
heuristic: it is order-dependent and not an exhaustive search. The final list is
capped per replacement, per fragment, and globally. Fragment order is deterministic;
a global cap can stop before later fragments are screened, so counts are not an
unbiased sample across fragments. Total and screened fragment counts are recorded. Generated unspecified stereo
is not automatically enumerated. A fragment matched in multiple input locations
is replaced at the first match only. Core protection checks every accepted product
with `useChirality=True`; other existing stereo outside the replacement is preserved.

## Similarity and filters

Scaffold search uses Morgan radius-2 **binary** fingerprints (2048 bits by default,
chirality disabled). Database metadata determines radius and size when a custom
fingerprint database is provided. Threshold comparisons are inclusive. Identical
canonical SMILES are excluded, but a different structure with Tanimoto=1 remains
eligible because fingerprints can collide. Generated databases are canonical,
sorted and unique under their recorded RDKit release. Known build/runtime version
differences are rejected before ranking because fingerprint bits may also differ;
regenerate into a new asset with the current runtime. The loader remains available
for read-only cross-version audits. With missing version metadata, legacy search
remains available with fingerprint compatibility unverified: it canonicalizes
selected hits, removes exact-query matches and equivalent representations, then
refills top-k as necessary. The string-only identity fast path requires matching
version metadata, the canonical/sorted/unique claim and the
`rdkit-roundtrip-fixed-point-v1` policy marker. This policy requires serialized
isomeric SMILES to be unchanged after reparsing; generation calculates bits from
that exact reparsed molecule. It detects cycles and permits at most eight
roundtrips. Unresolved source lines are rejected with original-input examples and
separate cycle/limit counts. Queries use the same bounded policy and fail if
unresolved; uncertified selected hits are checked and skipped with refill if
unresolved. No stereo removal, tautomer enumeration or charge normalization is
used to force convergence. Search does
not canonicalize the entire library before each bounded lookup. Atom maps are
removed from a query copy, preserving the caller's molecule; specified stereo,
isotopes and charges remain part of exact identity. Ties retain stored database
row order, including when canonical spellings change between versions.
Empty-vector similarity follows the installed RDKit convention
(1 for two empty vectors in 2023.09; 0 in the tested 2025.09/2026.09 releases).
Legacy databases retain their stored row order at score ties.

Before any potentially lossy parse/write, fingerprint building and selected legacy
rows reuse the exact-state stereo parser; raw query and fixed-point helpers reuse
the raw-Mol guard. Supported tetrahedral ABS is normalized with ordinary hydrogen
handling; OR/AND, unsupported axial states and other CX/name metadata are rejected
without choosing an absolute state. A selected row containing `|` or whitespace
cannot use the string-only shortcut, regardless of its metadata claims; guarded
normalization, exact-query/duplicate exclusion and bounded refill apply instead.
Ordinary canonical selected strings retain their existing no-reparse work path.
No full-library molecular parsing or new search cache is added. This is a check
against annotated-state/name bypass, not universal authentication of arbitrary
ordinary-SMILES metadata. Representation rejection counts/reasons are separate
from canonicalization cycle/limit diagnostics. Neither these checks nor sampled
bit verification can restore source annotations already erased upstream.
The text path uses canonical SMILES tie ordering.

Final candidate Tanimoto uses unfolded Morgan radius-2 **counts**, retaining the
legacy final-ranking definition. Binary and count similarities are different
metrics and their thresholds must not be interpreted as interchangeable. See the
[RDKit fingerprint API](https://www.rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html).

The default candidate limits are MW≤500, logP≤5, HBD≤5, HBA≤10. They are heuristics
for a particular chemical domain, not universal requirements. ADP violates donor
and acceptor limits. The input is not itself rejected for this; its violations are
reported and candidate filters remain unchanged unless the user changes them.
QED, SA score, logP and molecular weight come from RDKit and its SA_Score contrib.
Core and similarity checks precede property and optional 3D calculation.

## 3D descriptor

A hydrogen-complete molecule is embedded once with ETKDGv3, one thread and a
recorded integer seed. UFF minimizes it for at most 1000 iterations. Missing force
field parameters, failed embedding, nonconvergence, non-finite charges, or
geometrically degenerate descriptor reference points are reported as failures.
Explicit hydrogens are removed for the final descriptor. Gasteiger charges and
heavy-atom coordinates define ElectroShape, following
[Armstrong et al.](https://doi.org/10.1007/s10822-010-9374-0). Charge is multiplied by
25; five 4D reference points contribute mean, standard deviation and cube-root
third central moment. Similarity is `1/(1 + mean(abs(d1-d2)))`.

The NumPy implementation is compared numerically with the existing ODDT reference
on fixed coordinates. It eliminates the legacy package's installation coupling.
This is a single-conformer descriptor, not an ensemble model, docking calculation,
free-energy estimate or binding prediction. Rotatable nucleotides can have
substantial conformational variability. The seed assures repeatability within a
fixed software/hardware setup, not bitwise identity across RDKit releases.
A descriptor failure produces NaN and an explanation without removing the row.

## CATS2D

The retained SMARTS feature definitions and 21 pair types originate from
[iwatobipen/CATS2D](https://github.com/iwatobipen/CATS2D). All SMARTS alternatives
for a feature are evaluated and labels per atom are deduplicated. Pair lookup is
canonicalized, distances use RDKit graph shortest paths, and bins 0–9 are counted.
The inherited convention counts ordered atom pairs, including distance-zero
self-pairs. Raw output has 210 components and reported distance is Euclidean.
This convention must be specified when comparing with other CATS implementations.
The optional `num` and `occ` scalings divide by atom count and feature occurrence,
respectively. The 1.0 occurrence scaler reused its first divisor for all rows.

## Provenance and performance

Each manifest stores input/settings, actual library checksums, effective limits,
seed, dependency versions, Git state and Python source hashes captured at import.
Changes to source files after import are explicitly flagged in the manifest and summary. Each summary records
where candidates were excluded. CPU measurements refer to this process and use
100% per fully occupied core; RSS is sampled at 10 ms and may miss shorter peaks.
Results depend on the library and its composition; a smaller library is a distinct
experimental condition, not a drop-in optimization of the full library.

Packed bits reduce array storage by eight relative to uint8 bit matrices. Search
processes chunks and retains top-k scores without making whole-database AND/OR
matrices. NPY directory databases support mmap; NPZ loading decompresses the data.
The schema records the generator version, input digest, radius, width and bit order.
Verification recomputes sampled fingerprint rows against the runtime RDKit and
raises on bit mismatch. Passing a sample does not establish compatibility of
unsampled bits or canonical SMILES across releases. Search and verification leave
the stored fingerprint bits and metadata unchanged.
