# Administrative capacity-study preparation

`benchmarks/research/polar_capacity_prepare.py` constructs a candidate literal
for the already reviewed, two-input saved-supply development ablation. It does
not run a chemistry workflow or authorize the study. Actual preparation itself
requires a separate source/literal review, root authorization and current
no-overlap evidence. This source change includes synthetic tests only.

The public command is `prepare`. `coordinate`, `probe` and `handoff` are owned
internal children; they require a live PID/start-bound parent receipt, the exact
command/deadline and hash-bound spec/authorization. Run the script by its absolute
path with the checker as cwd/PYTHONPATH. Workers later import their unchanged
2443720/2556b14 arm; the checker may be the newly reviewed preparation checkout.
No source checkout or historical result is modified.

## Inputs and authority

The strict `polar-capacity-preparation-spec-v1` has these keys:

- `design`, `gate`, `reviews`, `clearance`: existing path/SHA-256 evidence records.
- `sources`, `checker`: exact path, commit, tree and full tracked-file SHA map.
- `contexts`: expected interpreter/version/package metadata and actual
  cwd/PYTHONPATH for baseline, capacity and checker, using the existing runtime
  schema. These are expected metadata, not invented runtime file hashes.
- `pins`: already authoritative design/review/source/historical/model/catalog/
  sidecar/rules/feature/interpreter/generation hashes. All reviewed historical
  bindings are required; model and scientific asset hashes are never discovered
  as replacements for approved hashes.
- `python`, `feature_definition`, `generation_source`: explicit interpreter
  invocation and path/hash records for the numerical feature asset and exact
  generation implementation. A venv interpreter symlink is preserved as the
  invocation; the actual executable bytes are pinned at its resolved path.
- `output`, `run_output`: distinct absolute fresh directories outside the source
  checkouts. Preparation never creates `run_output`.
- `workloads`: root-maintained PID/start identities and exact command markers
  covering all known scientific jobs. The historical campaign PID alone is not
  current clearance; root must provide terminal/absence evidence before launch.

`polar-capacity-preparation-authorization-v1` separately binds `spec_sha256`,
`source_review`, `literal_review`, `clearance`, `allow_asset_hashing: true`,
`allow_import_discovery: true` and `external_owner`. The latter must bind pinned
`source` and `review` path/hash records (the review must be in `reviews`), literal
`success_seconds: 300`, `watchdog_seconds: 330`, `cleanup_seconds: 30`, and scope
`owned-parent-and-descendants-through-cleanup-and-final-writes`. Missing or
different enforcement contracts reject preparation authority. This declaration
records required reviewed enforcement; it does not prove that owner is running.
The later literal gate must verify the actual owner invocation and terminal
receipt, including its PID/start ownership and all parent/descendant cleanup. Unknown fields and study-launch schemas are
rejected. Both spec and authorization file hashes are command-line inputs.
This tool never emits `polar-capacity-development-launch-v1` authorization.

Full tracked source maps are inputs because committed source authority already
exists. Runtime file inventories are constructed inside the bounded preparation:
the declared distributions, interpreter, actual generation source, urllib3
metadata/source, and files of observed transitive imports are hashed. Newly
measured runtime hashes are prospective candidate evidence, requiring independent
review. They cannot override any existing expected pin. Discovery and subsequent
verification independently derive the actual arm inventory in fresh processes.

## Fixed preparation sequence

The coordinator verifies approved bytes and reads only the bound historical
slot 0000/0009 requests, twelve ordered designs and six saved responses per slot.
It explicitly assembles the approved v2 policy; it never imports historical v1
thresholds or outcomes as new policy. Both arm probes resolve `_config` and
`_resolve_paths`, require idempotence and `engineering_policy_matched`, and must
produce equal scientific configurations. The old archive's limited response
coverage is preserved; no responses, molecules or routes are generated.

Each arm has a discovery probe and a separate verification probe. The existing
Observer installs/seals/closes its exact temporary wrappers. Discovery alone
locally seals its captured inventory to record instrumentation, not to approve
it. Verification supplies that inventory as an expectation and independently
reconstructs it. Original aliases and the ensemble callable/cache statistics
must be restored unchanged. Every scientific `started` event is forbidden.
Guard wrappers never replace inventoried original callables.

Probes use a fixed import/config/introspection program. Python profiling denies
scientific constructors, workflow, model loading/generation, scoring and related
entrypoints before invocation; observed wrappers also deny calls before their
original functions. Forbidden calls terminate the disposable probe after a
fsynced event so a swallowed exception cannot disable the profiler and continue.
Native/Boost.Python calls are not all visible to profiling: this is defense in
depth for reviewed fixed code, **not a general native-code sandbox**. No arbitrary
user Python or backend callback is accepted by the CLI.

Catalog/rules/sidecar opens are denied using both declared lexical and initially
resolved identities, including subsequent symlink retargeting. All model-root
descendants, including symlink aliases, are denied throughout probe runtime hashing/imports. Authorized model
and asset byte hashing occurs explicitly in the coordinator/handoff, without
constructing assets. Feature-definition reading and its import-time RDKit feature
factory are allowed; feature matching, molecule scoring and embedding are not.
Covered CPython socket audit events (connect, bind, sendto/sendmsg and address/name
resolution) are denied; socket creation remains allowed. This is not interception
of arbitrary native networking outside Python's audit hooks. The known urllib3 IPv6 capability bind may be
classified only with its exact pinned source/callsite/socket/address and one
blocked event; it still raises before binding. Attempted network events remain
in the evidence and are never reported as zero attempted operations.

The coordinator assembles the existing plan-v2 schema, runs its real validator,
and hands it to the unchanged development `prepare` command in a separate fresh
subdirectory. The output must contain exactly four `not_started` requests and
`prepared_not_authorized`. No fake scientific envelope or trusted accepted flag
is created. The candidate plus discovery logs go to independent literal review.

## Resource and terminal contract

This administrative operation has its own 300-second success allocation,
including process startup, all hashing/probes/handoff, bookkeeping, durable writes
and owned cleanup. It does not consume or enlarge the future study's 4,680-second
allocation. Process boottime age is converted conservatively into the monotonic
clock; clocks with different suspend offsets are never directly subtracted.
Every nested child receives only the shared deadline's remaining time.

The internal supervisor has a 330-second watchdog plus at most 30 seconds of
cleanup. The separately reviewed external owner must enforce this ceiling over
the entire preparation tree, including the parent and persistence/cleanup. These are failure ceilings, never successful extra preparation time.
Owned-tree RSS is sampled with the existing 16 GiB limit; gaps and per-child
usage remain explicit, and RSS is not an instantaneous memory guarantee. All
processes use CPU threads 1, empty CUDA and offline libraries. Unexpected known
workload overlap fails preparation; no unrelated job is signaled.

Parent final receipt persistence is checked after writing. Crossing 300 seconds
or a controlled interrupted write makes the candidate ineligible even if files
exist. Failure correction writes grant no new work allowance. The parent's one-shot
alarm is cooperative defense, **not proof of hard process termination**: native
calls or supervising cleanup may catch/delay it. Its sticky expiry state prevents
new parent reads or receipt/correction writes after expiry, even when a lower
layer catches the alarm. Deadline checks also apply before/after persistence.
Expiry leaves existing running/partial evidence and returns nonzero; no new
terminal receipt is required after the ceiling. An interrupted write may already
have persisted bytes, so a file alone never establishes eligibility. External
owned-tree enforcement and its terminal receipt remain mandatory. An owner
cannot guarantee immediate termination/reaping of uninterruptible kernel work:
unsafe or unreapable cleanup must remain an explicit failure, with retained
identities and unknown terminal costs; never report that a ceiling held when it
did not. A successful
preparation requires its external process exit 0 and intact receipt/candidate,
not a standalone file left by hard termination. Output reuse, missing pins,
config disagreement, changed hooks, unsafe cleanup or timeout stops preparation
without retries, budget expansion or automatic approval.

## Concrete preparation owner

`benchmarks/research/polar_capacity_preparation_owner.py run` is the fixed
external owner. It accepts the same `--spec`, `--spec-sha256`, `--authorization`,
`--authorization-sha256` and preparation `--output` as the preparation module.
Use the pinned interpreter and absolute script path, with the checker as cwd and
PYTHONPATH. It derives a sibling `<preparation-output>.owner` directory and an
`adapter` subdirectory; all must be fresh and disjoint from source/study outputs.
There is no backend command option, study launcher or automatic retry.

The internal `adapter` command additionally receives the owner's process-start
anchor. It authenticates the live parent and its own PID/start identities,
actual `/proc` commands/cwd, the exact supervisor receipt, interpreter, source
pins, complete environment digest and controlled CPU/offline fields. An owner
receipt and a matching command alone do not replace these live checks. The
owner and adapter use the same supervisor module object. No preparation begins
when the command, context, source, current workload clearance or shared lease
fails. The actual owner invocation and external exit must still be reviewed and
recorded at the separate literal/launch gate.

The first call to the supervisor clock returns the owner's conservative startup
anchor; every later call returns actual monotonic time. This charges supervisor
mkdir/pidfd setup and earlier owner work to the original 300/330/360 cutoffs.
The receipt labels this **owner elapsed scope**, not child lifetime. A final
pre-launch workload/lease check prevents launching after setup has consumed the
success allocation. The adapter calls the unchanged preparation `_run` with the
same anchor and never receives a second 300-second lease.

Within this owner/adapter process only, a scoped wrapper around `supervisor.write`
checks the sticky/absolute 360-second ceiling before delegating unchanged
arguments to the captured writer. It restores that exact callable in `finally`.
This covers the existing supervisor's caught-alarm final-write path without
changing supervisor source or signaling/reaping semantics. It prevents **new
durable writes**, not all I/O: necessary best-effort supervisor cleanup/proc and
terminal metadata inspection can still occur after expiry. New owner
administrative work stops after expiry. Already-started writes are not cancelled
retroactively, and uninterruptible/unsafe cleanup is never represented as a
proven hard bound or successful completion.

The adapter validates actual preparation and coordinator receipts, the candidate
plan using the existing full validator, the handoff terminal, and the actual
existing-prepare terminal and four-row `not_started` ledger. The latter terminal's
recorded lifetime RSS is retained as metadata, not interpreted as an instantaneous
or additive peak. These eight artifact files are hash-bound in the adapter
attestation. After child cleanup, the owner independently rereads the artifact
bindings and semantics; it does not repeat bulk model/runtime hashing or import
scientific packages. No boolean supplied by a generator is used as evidence.

Owner final artifact reads, safe cleanup and durable receipt writes remain inside
the 300-second success allocation. A write crossing that boundary is ineligible;
interrupted persistence or sticky ceiling expiry returns nonzero and can leave
running/partial receipts. A standalone file cannot prove terminal success. Child
RSS/CPU and sampled observation gaps come from the existing supervisor; the
owner's own memory is reported separately and peaks are not summed.

## Acyclic literal review

The fixed prospective review file is `literal-review.json` beside the declared
authorization file. Before it exists, `review_payload(spec_core, auth_core,
spec_path, authorization_path)` accepts records with that review entry/pin and
`auth.spec_sha256`/`auth.literal_review` absent. No placeholder digest is needed.
The payload binds those two absolute invocation paths plus every remaining
spec/auth value using sorted-key, finite canonical JSON; list order is retained.

After independent review, the review record has exactly `schema`, `status` and
`payload_sha256`, with values `polar-capacity-preparation-literal-review-v1`,
`approved_authority_payload` and the reviewed payload digest. Add its path/hash
record and pin to the spec, persist the final spec, then add its byte hash and
literal-review record to authorization. Source review, external-owner review,
clearance and every other field were already in the reviewed payload.

`authority_payload` reconstructs that payload by removing exactly the one
designated review entry/pin and those two authorization fields. Review-path
collisions with source, executable, protected authority or owner-review records
are rejected. No other mutation is permitted by the final payload comparison.
A final root gate, stored outside these inputs, rehashes final spec/auth bytes
and checks equality to the reviewed payload. Source pins, the review and the
spec therefore never need to contain their own eventual hashes. This source
implements the contract; it does not supply an actual literal or authorize a run.
