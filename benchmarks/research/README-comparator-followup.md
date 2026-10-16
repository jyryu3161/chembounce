# Separate comparator follow-up preparation

`prepare_comparator_followup.py` prepares the separately declared 270-attempt
study: 30 existing public queries, three seeds, and three methods. It has **no
execution command**. It never builds a fragment DB, loads a model, launches a
benchmark, retries a failure, or resumes an existing run.

The two parent decisions are immutable inputs, checked against their SHA-256
values in the tool. This is a post hoc supplementary study, not preregistration.
The original 630 attempts retain their frozen source and remain a separate
analysis. The new `pipeline_current` arm uses the tested detached source
`2328588d0969369e694d9293a1dccc89392a8203`, the isolated `.venv-research`
runtime, and the independently audited repaired fingerprint index. Its settings,
synthesis configuration, and effective 12-evaluation cap are copied unchanged
from the three original pipeline protocols.

## Inputs and preparation

Preparation requires three **completed and verified** CReM DB method fragments
and one finalized Mol2Mol method fragment. The fragment index is a JSON object
with exactly these keys:

- `crem_train_partial_exo`: an object containing `scaffold`, `time`, and `target`.
- `mol2mol_pretrained`: one artifact specification shared by all three splits.

Each leaf is `{"artifact": "/absolute/path/to/final-method-fragment.json",
"sha256": "the externally retained SHA-256"}`. Relative artifact paths, if used,
are resolved against the fragment index. Missing or changed assets fail before
the output directory is created. Do not create placeholder fragments to bypass
this requirement.

Run from the implementation checkout with the isolated research Python:

```text
python -m benchmarks.research.prepare_comparator_followup prepare \
  --decision /path/to/comparator-followup-parent-decision.json \
  --amendment /path/to/comparator-followup-current-pipeline-amendment.json \
  --fragments /path/to/final-fragment-index.json \
  --output /path/to/new-prepared-plan
```

The tool writes nine single-method protocols, the complete 270-slot grid,
public-query projections, separate execution/analysis/common-3D plans, and a
hashed manifest. The printed manifest SHA-256 must be retained outside the
prepared directory. Every file uses exclusive creation; existing output or
planned run directories are rejected. The planned run directory is a sibling
named `new-prepared-plan-runs`, so results cannot modify the plan manifest.

Verify with the externally retained anchor:

```text
python -m benchmarks.research.prepare_comparator_followup verify \
  --output /path/to/new-prepared-plan --manifest-sha256 RETAINED_SHA256
```

Verification rechecks source commits, input/output hashes, source-directory
membership and runtime distribution inventories. It does not grant execution
permission. The parent must separately assess the finalized plan, CReM build
evidence, and Mol2Mol observer/CUDA/process-group preflight, then record activation.

## Resource, provenance and information boundaries

All nine planned runner commands are sequential. GPU methods receive only the
physical GPU0 UUID from the amendment. CReM receives an empty
`CUDA_VISIBLE_DEVICES` and one CPU thread. The existing runner owns the wall-time,
memory and process-group limits; this preparation tool does not implement a new
scheduler. Every method has 360 seconds including cold method setup, 24000 MB
and at most three returned candidates. The external methods cap their native
work at 256; pipeline keeps 12 proposal evaluations. These are different units,
not matched universal objective-call counts.

The pipeline's 41 scientific files are checked against the exact tested Git
commit. External fragments bind their complete declared project-source
allowlists, upstream inventories, model/DB bytes, interpreter/venv configuration
and declared package versions. The Mol2Mol sampler explicitly uses RDKit
2023.9.5; its common postfilter and every runner use 2025.09.6. CReM metadata is
collected after a pinned **import-only** `crem` import, because its version helper
adds setuptools-vendored distribution metadata to `sys.path`. No model or DB
generation API is called. Version maps and origin files are not a hash of every
installed dependency binary; that limitation remains recorded.

Output paths resolve parent symlinks before checking protected trees. Inside the
original campaign checkout, writes are permitted only under its ignored
`benchmark-results/research-validation` subtree. The tested detached source and
the explicitly protected worktree are never output destinations. GPU2/GPU3 are
never placed in a command environment.

The runner is the trusted evaluator and can read labelled panels. Its existing
`public_request` whitelist passes only query ID, SMILES and target ID to methods;
no held-out labels, reference structures or panel paths enter method requests.
External assets have exact source/model/DB allowlists. This is information-flow
discipline, not OS filesystem isolation.

## Interpretation and remaining evaluation work

The primary endpoint is engineering delivery of at least one valid scaffold hop.
It does not establish useful biological activity or synthetic feasibility. The
analysis plan records exact evaluator metric names and reference populations.
Active-scaffold recovery has **no Tanimoto >= 0.5 gate**; the earlier control
prefilter-recall analysis used a different definition.

Common 3D assessment is an explicit separate pass with fixed conformer settings
and seeds. The basic benchmark runner does not compute these common 3D scores.
Every attempted outcome remains in coverage denominators, while missing scores
remain null in conditional similarity means. No supplied route means unassessed
route coverage. CReM edit reactions are not synthetic routes, and reuse of the
pipeline's forward predictor carries selection bias. Splits remain separate;
paired query/seed slots are aggregated by query and then equally by target.

Preparation does not duplicate scientific metric implementations. Final scoring
and paired analysis must reuse the pinned evaluator functions and undergo
independent assessment before comparative claims. The synthetic tests exercise
the plan contract without running real panels, DB construction or model inference.
