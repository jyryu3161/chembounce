# Official Mol2Mol pretrained sampling follow-up

This is a **preparation-stage, separate proposed study**. No benchmark sampling,
checkpoint inference or GPU preflight has been executed by this change. It does
not alter the frozen 630-attempt study. Parent review of the prepared artifacts
and a separate fixed-example preflight must precede inference. The prior decision
was made after observing earlier controls and partial original-campaign results;
this is post hoc study design, not preregistration.

The comparator is plain **pretrained Mol2Mol sampling with declared selection**,
not full REINVENT reinforcement learning, fine-tuning, or a RuSH reproduction.
The full chemical training inventory is unknown. Public training descriptions or
a public subset cannot establish clean holdout, and its training data are not
matched to this project's retrieval catalog.

## Official source and runtime

Use REINVENT4 [v4.4.22 commit
`eb61443a1a93a8cbd06510548b25cc77cfd9d3a6`](https://github.com/MolecularAI/REINVENT4/tree/eb61443a1a93a8cbd06510548b25cc77cfd9d3a6).
The exact Git asset is `priors/mol2mol_medium_similarity.prior`, 80,198,341 bytes:

- SHA256: `a0f40a2fb4baedaf51ae199bb38a61a8c8d26f03062977eaf06e34079e9f22ca`
- Official SHA512: `6c7a1153aecfeef7fafe0774a66f305f1e1a7bebd4ef8781cd3bd1ef113ff5629ec502c23384e24b28b4796b4e3449cf90df0135dec18048250cf9c8b1bb7c41`

`support.py` pins that checksum and the original `priors/sha512.sum`, LICENSE,
pyproject, and Linux lock. The selected code and checkpoint are distributed in
that Git tree under Apache-2.0. The [registered Zenodo prior release
metadata](https://api.datacite.org/dois/10.5281/zenodo.15641297) also specifies
Apache-2.0 and Copyright 2025 AstraZeneca; its concept DOI is
`10.5281/zenodo.15641296`. **Byte equivalence between the selected Git checkpoint
and a Zenodo file has not been established.** Git is the actual download origin.
Only this size- and SHA256/SHA512-verified official pickle may be loaded.

There is an upstream dependency inconsistency: this commit's pyproject requires
Torch `2.3.1+cu121` and Torchvision `0.18.1+cu121`, whereas its Linux lock lists
`2.2.1+cu121` and `0.17.1+cu121`. Preserve both original files. Resolve the
**pyproject**, with the compatible RDKit `2023.9.5` and NumPy `1.26.4` constraints
in this directory, into a new dedicated Python 3.10 environment. Preserve the
pip installation report, new resolved freeze, and `pip check`. This is not an
exact reconstruction of the contradictory upstream lock. All declared pyproject
requirements are checked, including OpenEye; no upstream dependency or code is
rewritten. Mol2Mol sampling does not invoke OpenEye's licensed scoring features.

The parent/common-filter process uses the separately pinned RDKit `2025.09.6`
environment and this repository's unmodified `chemical_identity` policy. The
sampler and common filter are separate interpreters under one outer process
group/deadline. Do not assume cross-version RDKit identity equivalence. All
actually imported `reinvent.*` Python files must resolve inside the pinned Git
checkout and match its tracked blob bytes. Third-party distribution versions,
Python binary hashes, pip checks and module origins are recorded; these records
are not a complete hardware or operating-system reproducibility guarantee.

## Preparation only

Use a **new** asset directory and venv; never install into an existing research
runtime. `COMMON_PYTHON` below denotes an existing read-only RDKit 2025.09.6
interpreter; `BASELINE_ROOT` is this isolated checkout. These shell variables are
illustrative and must be set explicitly to the intended separate paths.

```bash
mkdir "$MOL2MOL_ASSETS"
git init "$MOL2MOL_ASSETS/upstream"
git -C "$MOL2MOL_ASSETS/upstream" remote add origin https://github.com/MolecularAI/REINVENT4.git
git -C "$MOL2MOL_ASSETS/upstream" config extensions.partialClone origin
git -C "$MOL2MOL_ASSETS/upstream" config remote.origin.promisor true
git -C "$MOL2MOL_ASSETS/upstream" config remote.origin.partialclonefilter blob:none
git -C "$MOL2MOL_ASSETS/upstream" fetch --depth=1 --filter=blob:none origin eb61443a1a93a8cbd06510548b25cc77cfd9d3a6
git -C "$MOL2MOL_ASSETS/upstream" sparse-checkout init --no-cone
git -C "$MOL2MOL_ASSETS/upstream" sparse-checkout set --no-cone '/*' '!/priors/*' '/priors/sha512.sum' '/priors/mol2mol_medium_similarity.prior'
git -C "$MOL2MOL_ASSETS/upstream" checkout --detach FETCH_HEAD
"$COMMON_PYTHON" -m venv "$MOL2MOL_ASSETS/venv"
mkdir "$MOL2MOL_ASSETS/preparation"
CUDA_VISIBLE_DEVICES='' PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 nice -n 19 \
  "$MOL2MOL_ASSETS/venv/bin/python" -m pip install \
  --report "$MOL2MOL_ASSETS/preparation/pip-install-report.json" \
  -c "$BASELINE_ROOT/benchmarks/external/mol2mol/resolution-constraints.txt" \
  --extra-index-url https://download.pytorch.org/whl/cu121 \
  --extra-index-url https://pypi.anaconda.org/OpenEye/simple "$MOL2MOL_ASSETS/upstream"
CUDA_VISIBLE_DEVICES='' "$COMMON_PYTHON" \
  "$BASELINE_ROOT/benchmarks/external/mol2mol/prepare.py" \
  --upstream "$MOL2MOL_ASSETS/upstream" \
  --sampler-python "$MOL2MOL_ASSETS/venv/bin/python" \
  --model "$MOL2MOL_ASSETS/upstream/priors/mol2mol_medium_similarity.prior" \
  --installation-report "$MOL2MOL_ASSETS/preparation/pip-install-report.json" \
  --output "$MOL2MOL_ASSETS/prepared-v1"
```

This official commit contains actual blobs for two LFS-attributed files. If a
local LFS clean filter reports them as changed despite exact Git blob equality,
set repository-local `filter.lfs.process=''`, `filter.lfs.clean=cat`,
`filter.lfs.smudge=cat`, `filter.lfs.required=false`; do not rewrite file bytes or
change global Git configuration. Preparation compares tracked bytes directly.

`prepare.py` performs import-only inspection with CUDA hidden. It never unpickles
the model or calls sampling. It emits an immutable preparation bundle, actual
import origins, new resolved environment, full pyproject requirement checks and
a runner method fragment. A new output directory is required on each preparation.
The fragment gives the existing runner separate hash-pinned assets for every
adapter Python file, common identity source and checkpoint, as well as the
bundle. Nested bundle inventories also bind upstream Python files, interpreter hashes,
distribution versions, and the official model. The adapter revalidates these
before sampling; the worker repeats validation before importing upstream code.
Every adapter path in the command is therefore independently content-checked
by the runner before execution, rather than relying only on adapter self-checks. A final executable protocol must pin
this fragment/bundle and the dataset panel; no original protocol is edited.

## Fixed proposed sampling and selection

The proposed panel is the existing 30 public queries (12 scaffold, 12 time,
6 target), each with seeds 42, 43, 44: **90 future attempts**. The only molecular
input is the public query's canonical isomeric SMILES. Public query/target IDs are
bookkeeping. Labels, reference structures, split membership, target SMILES files
and `check_nll` are rejected; exact settings are whitelisted.

The outer `bounded_method` budget is 360 seconds including Python imports,
provenance validation, cold checkpoint load and selection; 24,000 MB sampled
aggregate host RSS; at most 3 returned candidates; at most 256 raw sequences.
Sampler batches are 32 with multinomial decoding, temperature 1.0, no input
randomization, one CPU helper thread and no CUDA-to-CPU fallback. Activation must expose exactly
one authorized GPU; the worker refuses a multi-GPU visible environment. The sampler is
seeded once per query; it is never reseeded between batches. Native NLL is the
unrounded official sampling output, not a recalculated held-out likelihood.
It is a within-query ranking statistic, not activity or synthesis probability.

The exact API is upstream `Mol2MolModel.load_from_file`, official
`Dataset`/`SMILESTokenizer`, then `Mol2MolModel.sample(..., "multinomial")`.
The high-level `Mol2MolSampler` applies MolVS `charge_parent` standardization;
this baseline deliberately does **not** invoke that extra input transformation.
Before the official tokenizer/Dataset encode operation, unsupported tokens are
rejected explicitly because upstream encoding otherwise silently drops them or
maps them to UNK. Input BOS/EOS-inclusive length must fit the official encoder's positional
buffer capacity (5000 in this source). The checkpoint's `max_sequence_length`
controls generated output and does not become an extra input restriction. Both
actual loaded values remain to be observed during the separately authorized
checkpoint preflight.

Upstream sampling starts with BOS and makes at most `max_sequence_length - 1`
decoding steps. Rows that reach EOS receive zero padding while other rows finish;
NLL includes sampled EOS and ignores padding. A transparent observer copies
integer arrays passed to the original vocabulary `decode`, preserving its return
object, exception and RNG. It binds the N input rows and N output rows in their
original order, validates token/text correspondence and padding, and labels
no-EOS length-limit rows as truncated. All raw strings, token IDs, NLLs, duplicate,
invalid and truncated records remain in local artifacts. Truncation is ineligible.
Fixture tests establish the observer contract, **not actual checkpoint ON/OFF
numerical equivalence**; the latter is a pending independent preflight.

Before starting another batch, reserve 15 seconds for filtering plus 1.25 times
the longest completed batch, in addition to the runner's 2-second finalization
margin. This fixed anytime rule can stop below 256 raw samples; it is not a
promise that a variable-duration model call will finish in time. Any outer
timeout, memory limit, nonzero exit or interruption accepts **no partial
candidates**, even if local raw/response artifacts exist. The adapter and sampler
inherit the outer runner's owned process group, including descendants. Do not
launch the adapter outside `bounded_method`; its group-leader guard rejects that
usage. The runner owns hard group-wide cleanup after success or failure. Before any
attempt artifact is written, existing owned filenames, dangling links and stale
atomic temporary outputs are rejected; prior evidence is never reused or overwritten.

In the common RDKit runtime, require valid connected molecules with rings;
exclude the input chemical family and unchanged state-normalized Murcko scaffold.
Normalization failures are recorded and excluded. Sort eligible rows by native
NLL ascending, then exact canonical isomeric SMILES, then raw index; deduplicate
chemical families and return the first three. Exact stereochemistry is preserved
in each returned structure, while the existing family policy deliberately
collapses stereo/isotope/protonation/tautomer variants for novelty counting.
A returned candidate has `route_evidence.level="none"`: no synthetic route was
supplied. This is unassessed route coverage, not a failed synthesis experiment.

Local artifacts include raw rows/batches/import origins, per-row selection
reasons, response SHA256 linkage and source/model/environment pins. PyTorch's
peak CUDA allocated/reserved bytes are reported separately from the outer
host-RSS metric; allocator figures do not measure total device memory or other
processes. Public aggregate results should omit molecular rows. A separate
common evaluator is needed for 3D and reference comparisons; this generator does
not secretly optimize those evaluator outputs.

## Model-free review and remaining validation

An [independent parent review](implementation_review.json) of implementation
commit `6e70ad9` passed 59 model-free tests and a separate ten-row literal
selection oracle. The parent also checked the real checkpoint checksums,
before/after source hashes and preservation of the original frozen 39 files.
The checkpoint was not loaded. This engineering result leaves runtime readiness
and actual checkpoint/GPU preflight pending; it does not validate scientific
accuracy or performance.

Run the model-free fixtures with the common runtime:

```bash
"$COMMON_PYTHON" -m pytest -q tests/test_mol2mol_baseline.py
```

After parent review, a separate fixed, non-benchmark input must establish official
checkpoint observer ON/OFF raw strings, token arrays, NLL and RNG-state equality,
actual max sequence semantics, CUDA placement/memory evidence and deadline/group
cleanup. Only then may the explicitly pinned new 90-attempt protocol be activated.
No claim of benchmark accuracy, synthesis success, GPU performance, or completed
REINVENT/RuSH reproduction follows from preparation and synthetic fixtures.
