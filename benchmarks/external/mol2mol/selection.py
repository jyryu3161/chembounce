"""Common-runtime selection; no neural or route model is called here."""

import math


def select(query, rows, maximum=3, *, deadline=None):
    import time
    from chembounce_synthesis.identity import chemical_identity

    reference = chemical_identity(query)
    annotated, eligible = [], []
    for index, raw in enumerate(rows):
        if deadline is not None and time.monotonic() >= deadline:
            raise TimeoutError("Common selection exceeded outer deadline")
        row = dict(raw, raw_index=index)
        reason = None
        if row["truncated"] or not row["terminated"]:
            reason = "unterminated_at_model_limit"
        elif not isinstance(row["native_nll"], (int, float)) or not math.isfinite(row["native_nll"]):
            reason = "nonfinite_native_nll"
        else:
            try:
                identity = chemical_identity(row["smiles"])
                row.update(
                    canonical_smiles=identity.exact_smiles,
                    family_smiles=identity.family_smiles,
                    scaffold_family_smiles=identity.scaffold_family_smiles,
                )
                if not identity.has_scaffold:
                    reason = "acyclic"
                elif identity.family_smiles == reference.family_smiles:
                    reason = "input_family"
                elif identity.scaffold_family_smiles == reference.scaffold_family_smiles:
                    reason = "unchanged_normalized_murcko"
            except ValueError as exc:
                reason = "invalid_disconnected_or_identity_failure"
                row["detail"] = str(exc)
        row["selection_status"] = reason or "eligible"
        annotated.append(row)
        if reason is None:
            eligible.append(row)
    eligible.sort(key=lambda row: (row["native_nll"], row["canonical_smiles"], row["raw_index"]))
    selected, seen = [], set()
    for row in eligible:
        family = row["family_smiles"]
        if family in seen:
            row["selection_status"] = "duplicate_family"
        else:
            seen.add(family)
            if len(selected) == maximum:
                row["selection_status"] = "output_cap"
            else:
                row["selection_status"] = "selected"
                selected.append(
                    {
                        "smiles": row["canonical_smiles"],
                        "native_nll": row["native_nll"],
                        "raw_index": row["raw_index"],
                        "route_evidence": {
                            "level": "none",
                            "reason": "No synthetic route supplied by sampler",
                        },
                    }
                )
    return selected, annotated
