# Conditional saved-supply capacity ablation

This tooling measures the capacity gate's downstream behavior on the entire
`qualified_full` slice of an already inspected development pilot. It does not
measure fresh performance, activity, experimental synthesis, or end-to-end
proposal generation. The 197 similarity controls motivate the gate but provide
no route-yield evidence.

The reviewed design is `polar-capacity-development-ablation-design-v2`, SHA-256
`352a62807ae320a04d943721641061e12b86ce2804154fe7d07606c5cc7d5cde`.
Its root design gate authorizes implementation only. A clean source review,
fully bound literal, separate launch authorization, and scheduling clearance
are still required. **No actual run accompanies this source change.**

## Fixed comparison

The baseline is `2443720dc74b740a9a1761fc30b53e28b667b3f4`; the capacity arm is
`2556b14db3e7c239352e37168c857ccb185cc5a1`. Both workers use the same existing
`workflow.run` path and unchanged numerical kernels/profile artifacts. The
order is pilot slot 0000 baseline/capacity, then slot 0009 capacity/baseline.
All four slots exist durably before the first child starts.

Each input retains all 12 ordered raw designs and all six archived response
records, including malformed responses. Responses cover only zero-based raw
indices `[1,1,2,2,3,3]`, with sampling seeds `[42,43,42,43,42,43]`.
Nine raw positions per input have no response. A newly requested uncovered
design produces `archive_coverage_failure`; no response is fabricated or
borrowed, and no generation, retry, rescue, or top-up is allowed. The results
remain conditional on this historically selected response availability.

Both arms require the declared selective engineering policy: Tanimoto/shape/
pharmacophore floors 0.35/0.70/0.50, polar recall 0.75, eight conformers with
screening seed 42 and final seeds 42/43/44. No core, pocket, reference SDF,
landmarks or carried routes are admitted. The actual config key
`reference_sdf` must be absent or null; unknown `reference_pose` is rejected.
Search/refinement/evidence caps are checked explicitly. Each worker may make
at most 24 new strict forward calls; upstream generator calls remain zero.
No cap or configuration changes are allowed after outcomes.

## Authority and entrypoints

`benchmarks/research/polar_capacity_development.py` provides `prepare`, `run`,
`worker`, `assess` and `admin`. Importing it performs no scientific work.
A future literal supplies absolute `--plan`, `--plan-sha256`, and `--output`.
`run` additionally requires a separately reviewed `--authorization` and its
`--authorization-sha256`; an approved design alone is not launch authority.
Workers and assessors are internal child entrypoints: the same authorization,
exact request, deadline, and live PID/start-bound parent launch receipt are
required. Direct scientific dispatch and output-directory reuse are rejected.

The plan schema is `polar-capacity-development-plan-v2`. It binds:

- The exact reviewed design and root design gate, ordered old inputs and query
  records, complete resolved scientific config, immutable arm commits/trees,
  full clean tracked-file inventories, and one fixed checker source checkout.
- Both profile artifacts, the catalog TSV and present provenance sidecar,
  rules, full local model bundle, feature definitions, interpreter, and runtime
  distribution files (excluding mutable bytecode caches). Before/after hashing
  is real measured work, including all pinned model bytes; it is not inference.
- Per-arm and checker runtime metadata in their actual cwd/PYTHONPATH contexts.
  Editable distribution metadata is not assumed equal to a root-cwd inventory.
- The observer's ordered hook inventory and original callable code/wrapper-chain
  identity, source/module hashes, aliases and signatures. Native callables have
  no Python code object: their bound binary/module/type limitation is explicit.
- The root-declared known workload identities and exact command markers, fresh
  output path, and administrative execution context.

A `polar-capacity-development-request-v1` is derived exactly from the plan.
The development wrapper authenticates it before calling the shared
`checked_workflow_artifacts` body. The existing benchmark `checked_products`
continues to authenticate its own benchmark request and expected config first;
no benchmark envelope is spoofed and no trusted boolean bypass is provided.

The single pinned checker runs separately from either arm. It freshly checks
route/atom/stereo evidence, raw token/beam ancestry, all three saved pose/witness
records, representative scores/energy/ranking, and JSON/TSV/all four SDF exports.
Empty exports must also be consistent. It permits tokenizer-only decoding and
saved-coordinate feature recomputation, while denying model construction,
weights tensor loading, generation, embedding, alignment and workflow execution.
Attestations retain exact route IDs and support bindings alongside molecules.
A changed route to the same product therefore remains visible.

## Resource and failure contract

No worker or assessor may overlap the active frozen campaign or another known
ChemBounce scientific job. Root must establish terminal/absence before launch.
Known PID/start identities and command markers are checked before children and
at least every second during supervision. Unexpected overlap stops this study;
other jobs are never signaled. This is not a shared-host reservation.

All children use CPU threads 1, empty CUDA visibility and offline local assets.
Workers have 600 seconds including startup and durable completion/owned cleanup;
TERM at 630 and KILL within 30 more seconds are failure cleanup ceilings only.
Assessors similarly have 300 seconds, with 330 plus at most 30 cleanup seconds.
A successful child returning at 601 or 301 seconds is still a failed outcome.
Setup and final administrative work each have 300 seconds, including parent
bookkeeping. The parent has an absolute 4,680-second outer watchdog and does
not launch work that cannot retain its full ceiling plus finalization reserve.
Uninterruptible kernel operations can prevent cleanup; this is an explicit
unsafe failure, never evidence of completed bounded work.

The Linux x86_64 supervisor uses pidfds (explicit syscalls for the pinned Python
build), identity rechecks, a subreaper and anchored session/descendant ownership.
It never falls back to an unverified `kill(pid)` or unconditional `killpg`.
Preexisting children prevent launch. Owned children must be gone before another
phase; unsafe setup cleanup also prevents finalization.

Owned-tree RSS is sampled every 0.2 seconds with a 16 GiB limit. An observed
exceedance or inaccessible live observation fails the phase. This is not an
instantaneous memory cap; short excursions may be missed and shared pages may
be counted twice. Sample gaps and peak scope are retained. Linux per-reaped-child
`wait4` CPU and lifetime maximum RSS are recorded separately; independent maxima
are never summed as a simultaneous peak. Parent memory is reported separately.
There is no `RLIMIT_AS`, and old RSS does not establish address-space feasibility.

## Observation and reporting

The observer never clears or replaces the original ensemble LRU. It records
started/returned/failed calls durably, preserving nonreturned work as unknown.
Aliases imported by selective scoring, refinement and workflow are inventoried.
Physical embedding/overlay keys include exact enclosing molecule/settings/seed
context; empty coordinate hashes or atom counts alone are not comparable keys.
Absolute deadlines and runtime object IDs are recorded separately, not used as
scientific matching keys. Logical requests and physical model/embedding calls
are separate layers, with actual cache statistics and reference requests.

Primary success requires authenticated, in-budget worker and assessor receipts,
complete cleanup and zero archive coverage failures, plus at least one validated
recommendation. Labels alone do not establish success. Each arm's denominator
is two; unique-family yield among top three uses denominator six. Missing,
not-started, failed, interrupted, late and unvalidated attempts remain zeros.
Diversity is the number of normalized/generic scaffold identities and is null
with a reason if fewer than two qualified families exist. Costs include failures
and startup; unknowns remain null. Physical call multisets retain multiplicity;
capacity rejections are not equated to avoided calls.

With two targets, report paired outcomes and arithmetic differences only. No
significance, winner-selection, causal wall-speed or GPU claim is supported.
Historical generator cost is shared sunk provenance, not new physical work.
There is no follow-up input selection if both arms fail or the gate never triggers.
