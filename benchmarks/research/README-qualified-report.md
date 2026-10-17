# Qualified benchmark snapshot review

`qualified_report.py` renders a captured plan and externally pinned collector
collection as one offline HTML file. It calls `qualified_analysis.analyze` to
validate bindings and recompute outcomes; it does not use saved analysis JSON,
read sidecars, inspect running campaigns, or run chemistry/models.

```bash
python -m benchmarks.research.qualified_report \
  --plan captured/plan.json \
  --collection captured/collection.json \
  --collection-sha256 "$TRUSTED_COLLECTION_SHA256" \
  --output new-review.html
```

Obtain the collection **byte** SHA-256 from a trusted external capture record.
Calculating a digest from an untrusted input does not establish collector trust.
The output must not exist; parents must already exist. Invalid inputs fail before
writing. Open the resulting file locally; it needs no network, JavaScript, or
external assets. Every artifact string is escaped.

The small Python API has the same trust contract:

```python
from pathlib import Path
from benchmarks.research.qualified_report import render_report, write_report

html = render_report(plan_bytes, collection_bytes, collection_sha256=trusted_digest)
write_report(
    Path("plan.json"),
    Path("collection.json"),
    collection_sha256=trusted_digest,
    output=Path("new-review.html"),
)
```

The report retains all planned denominators and shows exact hierarchical fractions
for all three endpoints. Seeds receive equal weight within query, queries within
target, then targets. Family deduplication is within attempt, with the first-ranked
representative retained. Missing records are distinguished from explicit
`not_started`; `no_candidates` is separate from operational statuses. Full record
coverage does not establish operational completion. Per-attempt rows locate zero
outputs and their recorded failure reasons.

Failure records can overlap and are not a sequential funnel. Stage tables show
both records and affected attempts; reason tables count records. Costs preserve
unknowns, include failed-attempt observations, and show per-field coverage. Time
sums and memory maxima describe observed attempt measurements only. Separate
administrative setup, its allocation shares, and campaign checkpoint time are
omitted because the analyzer does not validate those auxiliary fields. Attempt
cost sums are neither total campaign costs nor campaign duration. Conditional
similarity/margin summaries include all selected outputs, including duplicate
families, and are not planned-attempt or equal-target distributions.

Numeric field coverage means record availability, not complete resource
measurement. The captured `measurement_policy` appears beside the cost table as
escaped recorded provenance, distinct from analyzer-validated arithmetic. Missing,
empty, or non-text policies are unavailable; no scope is inferred for them. A
separately labeled current-runner reference explains its sampled child-process-tree
CPU lower bounds and sampled child-tree RSS, excluding parent/collector CPU and
RSS. That reference is not automatically applied to historical collections. The
combined parent+child memory peak remains unknown; separate maxima are never
summed. No sidecars or new measurements are used.

Input byte digests, canonical JSON bindings, and scope are in collapsible details.
These are collector-consistent computational outcomes, not fresh chemical
validation, activity, successful synthesis, or experimental evidence. The view
makes no efficacy, infeasibility, significance, or winner claims. It intentionally
omits intervals and paired differences. A saved single-query smoke render is only
a compatibility check, not comparative scientific evidence.
