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
`source_review`, `literal_review`, `clearance`, `allow_asset_hashing: true` and
`allow_import_discovery: true`. Unknown fields and study-launch schemas are
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

Catalog/rules/sidecar opens and all model-root descendants, including symlink
aliases, are denied throughout probe runtime hashing/imports. Authorized model
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

The outer supervisor has a 330-second watchdog plus at most 30 seconds of
cleanup. These are failure ceilings, never successful extra preparation time.
Owned-tree RSS is sampled with the existing 16 GiB limit; gaps and per-child
usage remain explicit, and RSS is not an instantaneous memory guarantee. All
processes use CPU threads 1, empty CUDA and offline libraries. Unexpected known
workload overlap fails preparation; no unrelated job is signaled.

Parent final receipt persistence is checked after writing. Crossing 300 seconds
or a controlled interrupted write makes the candidate ineligible even if files
exist. Failure correction writes grant no new work allowance. A successful
preparation requires its external process exit 0 and intact receipt/candidate,
not a standalone file left by hard termination. Output reuse, missing pins,
config disagreement, changed hooks, unsafe cleanup or timeout stops preparation
without retries, budget expansion or automatic approval.
