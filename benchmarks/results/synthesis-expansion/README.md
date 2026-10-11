# Building-block expansion comparison

Recorded 2026-10-11, application revision `f9fdd2d`. This is one retrospective
integration case with fixed model proposals, not a held-out accuracy benchmark.

The four stored predictions come from the actual SynLlama run documented in
[the initial validation](../synthesis/README.md). The original reference is the
same hypothetical oxazole variant, and the design is the same published thiazole
example. The model was **not resampled** for this comparison. Both runs use the
same 238,542-entry catalog with the exact aminothiazole precursor and its
deuterated counterpart marked unavailable (238,540 eligible entries), the same
91 rules, and the same six-conformer ensemble and seed 42.

The required core `c1cccnc1N1CCOCC1` and the design Murcko scaffold are preserved.
Minimum similarity cutoffs are disabled. Projection uses top-k 3 and beam 6;
the per-starting-route combination bound is 256, and the workflow's 3D/result
budget is 64. Only expansion settings differ between cases.

| Setting/result | Small expansion | Broader expansion |
| --- | ---: | ---: |
| Alternative blocks per leaf | 3 | 16 |
| Substitution rounds | 1 | 2 |
| Retained routes per starting route | 12 | 128 |
| Unique returned products | 4 | 35 |
| Unique structures after isotope removal | 4 | 35 |
| Products preserving the required core | 4/4 | 35/35 |
| Products sent to 3D scoring | 4 | 38 |
| Elapsed seconds, including catalog load | 78.46 | 129.53 |

All returned products differ from the design after isotope removal. Three broader
search products could not be scored because their stereochemistry was unspecified;
they were explicitly rejected. Search caps remain active; the 64-product 3D/result
budget was not reached. The runtime comparison is descriptive, not a controlled
hardware performance benchmark.

For the **single best joint-3D candidate in each run**, the metrics are:

| Metric | Small expansion | Broader expansion |
| --- | ---: | ---: |
| Tanimoto vs original hit | 0.554 | 0.576 |
| Tanimoto vs design | 0.729 | 0.774 |
| Shape vs original hit | 0.885 | 0.712 |
| Pharmacophore vs original hit | 0.758 | 0.357 |
| Shape vs design | 0.631 | 0.996 |
| Pharmacophore vs design | 0.354 | 1.000 |
| Joint 3D/energy score | 0.622 | 0.709 |

These are all components from each candidate's selected conformation, not
independently selected maxima. The broader run's leading product replaces the
design's methyl substituent with bromine. Its high design alignment scores do
not establish conserved receptor interactions, activity or synthesis success.
The original-hit pharmacophore component is lower: more candidates improve the
available choices, but do not guarantee improvement in every objective. All six
components remain visible in the result table and contribute to Pareto ranking.

Independent export checks reproduced every Tanimoto exactly and pharmacophore
scores from exported coordinates within 0.0001. The two candidate poses represent
the same conformation within 0.001 Å pair-distance error. All selected routes were
replayed independently using their blocks joined back to the actual catalog, and
all final structures matched. See [the audit](export-audit.json) and
[machine-readable comparison](comparison.json). The selected best candidate and
reference coordinates are included; complete local records are under
`benchmark-results/expansion-comparison` and their hashes are recorded.

Tests: 141 passed in the primary environment; 138 passed/2 optional modules skipped
with modern RDKit; 115 passed/4 optional modules skipped with minimum dependencies.
The minimum environment continues to support hopping alone without synthesis 3D
or model dependencies. Logs are included alongside this report.

## Reproduction

Prepare the assets and actual model predictions as described in the initial
validation. Create a separate diagnostic catalog and a saved-proposal config:

```python
import csv
import json
from pathlib import Path
from chembounce_synthesis.workflow import load_config

work = Path("benchmark-results/expansion-input")
work.mkdir(parents=True, exist_ok=False)
config = load_config("benchmarks/synthesis_smoke_config.json")
withheld = {"Cc1csc(N)n1", "[2H]C([2H])([2H])c1csc(N)n1"}
found = set()
catalog = work / "catalog.tsv"
with open(config["catalog"], newline="") as source, catalog.open("w", newline="") as target:
    reader = csv.DictReader(source, delimiter="\t")
    writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t")
    writer.writeheader()
    for row in reader:
        if row["smiles"] in withheld:
            row["availability"] = "unavailable"
            found.add(row["smiles"])
        writer.writerow(row)
assert found == withheld
config["catalog"] = str(catalog.resolve())
config["provider"] = {
    "kind": "saved",
    "path": str(Path("benchmark-results/synthesis-live-hop-final/proposals.json").resolve()),
}
config["core_smiles"] = "c1cccnc1N1CCOCC1"
(work / "config.json").write_text(json.dumps(config, indent=2))
```

Then run:

```bash
python benchmarks/compare_synthesis_expansion.py \
  --config benchmark-results/expansion-input/config.json \
  --hit 'Cc1coc(NC(=O)c2cccnc2N2CCOCC2)n1' \
  --design 'Cc1csc(NC(=O)c2cccnc2N2CCOCC2)n1' \
  --output benchmark-results/expansion-comparison-reproduced --max-scored 64
```

This compares two budgets under the updated algorithm. It is not a direct
before/after count comparison with the earlier three-product diagnostic, which
used the previous expansion algorithm. No experimental synthesis or biological
assay has been performed as part of this validation.
