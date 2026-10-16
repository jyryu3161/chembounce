# Supplementary comparator campaign supervisor

This supervisor implements orchestration for the separately prepared nine
campaigns and 270 attempts. It does not run a study during preparation or testing.
Use it only with the independently reviewed immutable preparation directory and
its separately retained manifest SHA-256. The original 630-attempt study, source,
models and prepared protocols are not changed.

The supervisor calls `prepare_comparator_followup.verify_prepared` before launch
and before each campaign, then checks the exact ordered public query/seed grid,
nine commands and their GPU visibility. It executes each plan `argv` directly,
with its exact `cwd` and environment overrides applied to the inherited
environment; there is no shell evaluation. CReM sees an empty CUDA device list;
pipeline and Mol2Mol see only the declared physical GPU0 UUID. The supervisor
does not call `nvidia-smi`, inspect GPUs2/3, load models or alter runtime packages.

After a separate execution decision, launch from the supervisor checkout using
a new runtime directory outside both the immutable preparation and planned runs:

```bash
"$COMMON_PYTHON" benchmarks/diagnostics/supervise_comparator_followup.py \
  --plan "$PREPARED/execution-plan.json" \
  --manifest-sha256 "$INDEPENDENTLY_RETAINED_MANIFEST_SHA256" \
  --runtime "$NEW_SUPERVISOR_RUNTIME" --detach
```

Omit `--detach` for foreground operation. These commands are usage examples,
not evidence that a study has been activated. The detached launcher records the
worker PID, creation time, Linux start ticks and its own process group. Its zero
exit means the worker was launched, not that the study completed. Monitor
`progress.json`; only `status="completed"`, nine completed campaigns and 270
validated attempts establish supervisor completion. Comparative analysis still
requires the separate strict assessment gate.

A read-only descriptor of the immutable manifest carries an exclusive advisory
`flock` through the detached worker. This blocks another cooperating supervisor
even if it chooses a different runtime directory, without writing into the
prepared directory. Existing runtime paths, dangling links and any existing
planned run tree are refused. There is no retry, resume, backfill or repair mode.

Campaigns run sequentially in separate owned process groups. The existing runner
enforces each attempt's 360-second and 24,000 MB limits. Exit 1 can represent
retained terminal failures, so the supervisor continues only after verifying
the completed campaign marker, exact planned rows, per-attempt evaluation joins,
response hashes and source/protocol/dataset provenance. Missing, duplicated,
reordered, interrupted or unplanned outcomes, a contradictory exit status, an
unexpected exit code or changed immutable inputs stop the study. A timeout
inside a completed campaign remains a failed attempt; it is never retried.
Failed evaluations must contain no accepted candidates. Their raw/report
responses may retain unaccepted candidates (for example `budget_exceeded`);
these remain failure evidence and do not stop an otherwise complete campaign.

SIGTERM/SIGINT record interruption, stop further scheduling and send SIGTERM
only to the active child runner group after checking its PID, birth identity
and group ID. After 15 seconds a still-verified group may receive SIGKILL.
The runner owns cleanup of its own separately grouped method workers; the
supervisor never enumerates or kills unrelated groups. A killed/unresponsive
runner cannot provide the same graceful cleanup guarantee. SIGKILL or host loss
cannot be checkpointed, and a stale `running` record must not be promoted to
completion. An interrupted runtime is retained and cannot be resumed.

The runtime contains immutable launch/start records, detached-worker identity,
per-campaign stdout/stderr logs and atomically replaced progress. Progress binds
the launch, exact command hashes, child identities, exit codes, log hashes,
validated completion-file hashes and validated attempt counts. An active
campaign's observed row count is explicitly unvalidated until the completion
checks pass. No scientific estimates or method comparisons are produced here.

Model-free fixtures exercise sequential full-grid execution, a completed campaign
with a retained timeout, incomplete/unexpected termination, exact-grid/GPU/order
refusal, immutable paths, duplicate locks, birth-identity refusal and signal
handling without touching an unrelated process. A real detached subprocess
fixture uses synthetic runner files and the actual preparation validator; it
loads no scientific model or benchmark panel.
