# Mol2Mol fixed-example preflight plan

This diagnostic **passed on 2026-10-16 after separate parent review and explicit
activation**. The [aggregate result and evidence hashes](../external/mol2mol/actual_preflight_review.json)
preserve its scope: one fixed nonbenchmark input, observer comparison and one
actual adapter attempt. The historical preparation records remain unchanged.
Execution did not modify the [reviewed baseline](../external/mol2mol/README.md),
its official weights, or any original study protocol. Comparative execution is
a separate gate.

The predeclared nonbenchmark input is `CCOc1ccccc1`, seed 42. Model-free inspection
of the three SHA-pinned original panels found no exact canonical match among
30 public queries. This is not training-set or held-out membership proof. If the
check later fails, the tool stops instead of selecting a different example.

The diagnostic has two stages:

1. Load the verified official checkpoint once on logical `cuda:0`, exposing only
   physical GPU0 UUID `GPU-ea3de395-f78f-9f79-e363-cfac7415b6fd`. Seed Python, NumPy
   and Torch before checkpoint loading exactly as the baseline worker does. Use
   the official Dataset/tokenizer, 32 identical inputs and multinomial temperature
   1.0. Save all RNG states after model/Dataset setup. Sample once with the decode
   observer OFF; a scoped Python profiler copies the exact unwrapped official
   sample function's `src` and `ys` only on its return event. Restore the pre-sample
   RNG state, then sample with the existing observer ON. Compare input/output
   strings, full token arrays, native NLL dtype/shape/bytes, and full post-call
   Python/NumPy/Torch CPU/every-visible-CUDA RNG states. Check decode and profiler
   restoration and the observer's EOS/padding/truncation binding. Tensor copies
   can synchronize CUDA: timings are diagnostic observations, not performance
   comparisons. Existing profilers are refused, not overwritten.
2. Run one real `bounded_method` adapter attempt on the same fixed input, unchanged
   360-second cold-start budget, 24,000 MB host-RSS cap, up to 256 raw samples and
   three returned candidates. Preserve the real outcome, including failure;
   `no_candidates` is an acceptable lifecycle outcome if the full adapter succeeds.
   Check response hash linkage and independently observed sampler PID/GPU identity.
   This is an engineering diagnostic, not an additional benchmark query.

The comparison worker records every network parameter/buffer and input tensor
on `cuda:0`, actual input positional capacity and output length cap, its own PID
and `nvidia-smi` match, driver/device identity and Torch allocator peaks. Before
sampling, its separate controller must observe that PID alive on the approved
GPU and write an acknowledgement. The controller also polls live compute-apps
while the real adapter runs. No inference or process-detail collection on GPUs2/3 is permitted.
All SMI queries are restricted to the approved GPU. A busy GPU at either stage
entry causes refusal without killing any foreign process. Controller and worker
PID creation times/start ticks are recorded to distinguish process reuse. Actual
upstream module origins are checked again after checkpoint loading and sampling.
Host limits, owned process-group termination and failure artifacts are explicit.
Existing diagnostic outputs are refused. Raw strings and RNG details stay in
local diagnostic artifacts; public reports should use aggregate checks/hashes.

Prepare the immutable plan after `prepare.py` has produced its validated bundle
and method fragment; this command imports no neural runtime or checkpoint:

```bash
"$COMMON_PYTHON" "$BASELINE_ROOT/benchmarks/diagnostics/mol2mol_preflight.py" \
  --prepare-plan \
  --bundle "$MOL2MOL_ASSETS/prepared-v1/bundle.json" \
  --method-fragment "$MOL2MOL_ASSETS/prepared-v1/method-fragment.json" \
  --panels "$PANEL_SCAFFOLD" "$PANEL_TIME" "$PANEL_TARGET" \
  --output "$MOL2MOL_ASSETS/fixed-preflight-plan.json"
```

Review the resulting plan, its helper/source/model/runtime pins and fixed GPU
identity before any execution. The explicit execution command below is provided
as a record of the separately authorized step. Existing result paths must never
be reused or overwritten:

```bash
"$COMMON_PYTHON" "$BASELINE_ROOT/benchmarks/diagnostics/mol2mol_preflight.py" \
  --execute-plan "$MOL2MOL_ASSETS/fixed-preflight-plan.json" \
  --expected-plan-sha256 "$REVIEWED_PLAN_SHA256" \
  --output "$MOL2MOL_ASSETS/fixed-preflight-result"
```

Model-free helper fixtures exercise scoped return-frame observation, original
exception preservation, RNG neutrality, bitwise NLL comparison, complete RNG
stream comparison, plan tampering, exact-query exclusion and immutable output
refusal. The later actual execution passed all nine observer checks, confirmed
260 parameter tensors, two buffers and the input/mask on `cuda:0`, and observed
the official input capacity of 5000 and output sequence cap of 128. The actual
adapter returned three candidates from 256 EOS-terminated raw rows. The parent
independently checked exact strings/NLL bytes/tokens/RNG, selection, response
provenance and 544 unchanged pins without repeating model inference.

Live GPU evidence is the reviewed controller/worker recording with PID and
creation-identity joins. The parent's additional manual SMI snapshot was taken
after completion and observed an idle device; it is not another live observation.
Comparison and adapter wall times were 4.4099 and 6.4980 seconds, respectively;
sampled adapter host RSS was 979.68 MB. These are single shared-host observations,
not performance comparisons. CUDA allocator peaks are separately scoped to
PyTorch; actual Triton kernel execution was not assessed.

All three outputs retained the query's normalized Murcko core despite having
different normalized Murcko scaffolds. The hash-bound secondary annotation
addendum in the aggregate report records this distinction without changing method
eligibility. The diagnostic establishes this one example's engineering checks,
not held-out accuracy, general observer or CPU/GPU equivalence, 3D preservation,
activity, route replay, synthesis success, or comparative runtime performance.
