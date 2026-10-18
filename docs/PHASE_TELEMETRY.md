# Optional bounded phase telemetry

`chembounce_synthesis.phase_telemetry.PhaseTelemetry(phase)` observes existing
work in one process. It is optional instrumentation, not a scientific method or
an execution owner. The exact accepted `bounded-phase-telemetry-v1` module has
39 fixed buckets (12 native, 9 supply, 18 common), maximum depth four, aggregate
counters, and no event list. Enabled construction uses Linux `/proc` identity.

Pass an explicit `telemetry=` keyword to `chembounce`, `make_scaffold_hopping`,
`benchmarks.research.qualified_adapter.native_designs`, or synthesis `workflow.run`.
Existing positional arguments retain their meanings. There is no CLI flag or
JSON configuration option; the adapter rejects telemetry in generator options.
No context is created automatically, and no telemetry files or result fields are
published automatically. For example, an existing caller may write:

```python
from chembounce_synthesis.phase_telemetry import PhaseTelemetry

telemetry = PhaseTelemetry("common")
with telemetry.span("common.worker"):
    with telemetry.span("common.workflow"):
        rows = run(hit, designs, config, output_dir, telemetry=telemetry)
snapshot = telemetry.snapshot()
```

The caller must retain its existing budget/termination authority. A snapshot
only describes timing coverage; it does not authenticate saved execution, source,
inputs, results or chemistry. This port does not include a comparison harness,
worker, saved-envelope interpreter, provider wrapper or execution release.

| Phase | Production scopes | Coverage |
| --- | --- | --- |
| Native | `setup` | Existing parameter/provenance and fragment/library setup; repeated intervals |
| Native | `generation`, `allocation` | Existing diverse generation and unchanged allocator |
| Native | `output` | Existing diverse writes and public CSV read; repeated intervals |
| Common | `reference_setup`, `reference_admission` | Existing reference setup and design intake |
| Common | `catalog`, `rules`, `evidence_prepare` | Existing asset/evaluator construction |
| Common | `route_search`, `forward_evidence` | Existing search and evidence queue |
| Common | `screening`, `refinement`, `selection` | Existing scoring/refinement/final selection |
| Common | `workflow_finalization` | Existing finally publication/display |

Native legacy generation/allocation detail is explicitly `unavailable`; its
public CSV read is timed, while legacy generation writes remain in the caller's
API residual. Native reference-allocation is not included: its reserved
`native.reference_preflight` bucket stays `not_reached` in both paths, with null
durations. Supply and caller-owned buckets may also remain unused. The unchanged
registry is a vocabulary, not a claim that every operation exists or ran. Disabled
refinement is `skipped_by_policy`; stages never entered remain `not_reached`.
No CATS reserve or native reference-admission intervention is installed.

Counters sum inclusive/exclusive active wall and self-process CPU nanoseconds.
CPU includes this process's threads, excludes descendants, and is not added to
external owner CPU. Never sum nested inclusive measurements. Repeated scope
first/last timestamps are envelopes, not active durations. Uninstrumented work
remains parent-exclusive or outside spans. PID/start ticks/boot/thread and clock
names accompany snapshots; contexts cannot cross process/thread or stack order.
No CUDA synchronization or extra scientific work is added.

Exceptions propagate and mark interrupted scopes. Cooperative deadline markers
use an already-recorded workflow flag; they do not read a new deadline or change
stop conditions. A completed root can contain an interrupted child. `closed`
means spans closed, not successful execution. Hard termination can leave no
snapshot; missing timing is unavailable, never zero or molecular infeasibility.

Close all root spans before snapshot creation and serialization. Snapshot/export,
final writes, process overhead and work outside spans remain externally charged
and uncovered by the closed snapshot. This module grants no extra budget and
cannot recursively time its own final export.

`telemetry=None` performs no telemetry clock/proc reads, records or writes. Import,
function and null-context plumbing still has overhead. Enabled instrumentation
consumes the existing budget and can alter deadline-limited coverage; structural
parity is not a speedup or equal-coverage claim.

Source hashes change honestly. Numerical kernels, profiles and thresholds are
unchanged. v7 primary numerical binding excludes workflow; optional rejected-pose
diagnostics still expose workflow mismatch when enabled. v6 workflow binding and
historical source drift are preserved, not rebound. The new disabled-path parity
test pins exact develop `94b5c014` bytes, ignores only known empty declaration
`type_params` fields added in Python3.12, and needs no Git history in CI.
