# Partial training-sample audit — 2026-10-16

No complete, exact-checkpoint training inventory was established for the campaign's SynLlama or ReactionT5 models. An author-published ZINC dataset exists, but its molecular contents and exact forward-checkpoint association were not verified here. [EVIDENCE.md](EVIDENCE.md) separates these findings and links the primary releases, papers and code revisions.

The post hoc audit compares the upstream **train-labelled 1,000-molecule evaluation sample** with 14 current evaluation-input artifacts. At the recorded read time, 19,387 repeated appearances represent 5,618 distinct exact structures and 5,616 families. There are **zero exact matches to the 1,000 sample structures** and **zero family matches to 928 successfully normalized sample families**. The other **72 sample structures remain unassessed for family matching** because tautomer enumeration reached the existing transform limit. Exact identities remain available for all 1,000; no evaluation-input identity failed.

This establishes neither that the sample belonged to the exact 2M checkpoint's fine-tuning rows nor that evaluation molecules were absent from its full training corpus. It does not assess broader scaffold similarity. The snapshot covers 137 recorded main-campaign runs at **2026-10-15 18:27:40–41 UTC**, while campaigns were still running. Later candidates and empty secondary-generation outputs are not silently treated as audited outcomes.

The [aggregate](aggregate.json) retains input and code hashes, read times, role-specific denominators, duplicates and normalization coverage. The [inventory evidence](inventory-evidence.json) distinguishes complete/partial/unverified provenance and code/model/underlying-data licenses. The [fetch manifest](fetch-manifest.json) contains source URLs and hashes only. Full papers, the source molecular sample, molecular row snapshots and detailed match/failure records are **not included**. Synthetic SMILES in the unit tests are software fixtures, not source data.

## Reproduction and limits

Use the repository environment with the recorded RDKit **2025.09.6** and `protonation-tautomer-family-v1` implementation; the aggregate pins the source-file hashes. From the repository root, test the bundled helper:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. python \
  benchmarks/results/training-inventory-audit-20261016/test_sample_overlap_audit.py
```

The helper verifies the full SHA256 of an already obtained upstream archive, reads only the named regular member, and refuses existing output files. The original archive is 7.77 GB; it is not bundled or automatically downloaded. Its [immutable release](https://doi.org/10.6084/m9.figshare.28719125.v1) and licensing distinctions are documented in the evidence report.

```bash
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=.
bundle=benchmarks/results/training-inventory-audit-20261016
audit_dir=benchmark-results/training-inventory-audit-local
python "$bundle/sample_overlap_audit.py" extract \
  --archive /path/to/synllama-data.tar.gz \
  --archive-sha256 d46ceb49c3e2949c6e65386187d07c8866852de93f20d710daebb57597acebed \
  --member ./synllama-data/inference/smiles/syn-planning/1k_train_91rxns.smi \
  --source-url https://doi.org/10.6084/m9.figshare.28719125.v1 \
  --output "$audit_dir/1k_train_91rxns.private.smi"
python "$bundle/sample_overlap_audit.py" audit \
  --sample "$audit_dir/1k_train_91rxns.private.smi" \
  --sample-sha256 fa00401f509a4b77d121f22824400aae63ecba48c665f9c46d4e9fbe9b7940ac \
  --spec "$bundle/audit-sources.local-example.json" \
  --output "$audit_dir/recomputed"
```

The [example source specification](audit-sources.local-example.json) uses checkout-relative paths for locally reconstructed panels and archived campaign files. It performs no download; adjust paths to available inputs. The [source curation protocol](../../../docs/RESEARCH_VALIDATION.md) and [assay-matched report](../assay-matched-controls-20261016/README.md) describe the public ChEMBL sources. These records are attributed to [EMBL-EBI ChEMBL](https://www.ebi.ac.uk/chembl/) under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/); this does not change the repository code license or establish rights to redistribute the SynLlama sample.

For an exact replay of the reported live-campaign read, use `--snapshot /path/to/input-rows.private.json` instead of `--spec`. That private row snapshot has SHA256 `0d555f7b71b2041c7e4baae9b0357ac965ae9a9cfd5214ed21fd6835f56fa30a`; it is not publicly bundled. Consequently this bundle provides the procedure and hash-bound aggregate, **not a self-contained public replay of the complete live snapshot**. A fresh `--spec` run after campaigns advance is a new audit and must retain its new hashes/time. Every output directory must be new.

Five helper tests passed, including an independent parent rerun. The [independent projection audit](projection-validation.json) checked all 19,387 input appearances across 14 sources, verifying immutable source bytes and the retained attempts/identities in growing files. It did not reconstruct the complete original bytes of those growing files. The [independent numerical audit](numerical-validation.json) passed at **2026-10-15 19:03:19 UTC**: 6,618 normalization inputs and 100 group-level checks reproduced the identities, overlaps, sample-line indices and aggregate consistency. Exact canonicalization and set/count aggregation were implemented independently of the helper; family normalization reused the project's versioned identity policy, so this is not independent validation of that policy's chemistry. [validation.json](validation.json) records these checks and verification of 33 retained source documents and three audit/identity implementations. The frozen campaign's 39 source hashes remained unchanged.
