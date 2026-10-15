"""Export an unlabelled end-product retrieval library from development compounds.

No test/query-partition chemical family is admitted. This is a benchmark search
library, not a supplier catalog and not an assertion of chemical availability.
"""

import argparse
import csv
from pathlib import Path

from benchmarks.research.protocol import digest, load_dataset, write_json


def export_library(dataset_path, output):
    dataset, audit = load_dataset(dataset_path)
    output = Path(output)
    if output.exists():
        raise ValueError("Library exists")
    output.parent.mkdir(parents=True, exist_ok=True)
    private = {
        dataset.identities[r["id"]].family_smiles
        for r in dataset["compounds"]
        if r["split"] in {"test", "query"}
    }
    selected = {}
    for row in dataset["compounds"]:
        if row["split"] != "train":
            continue
        identity = dataset.identities[row["id"]]
        if identity.family_smiles not in private:
            selected.setdefault(identity.family_smiles, identity.exact_smiles)
    if not selected:
        raise ValueError("No development compounds remain outside held-out families")
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["id", "smiles", "source", "snapshot"], delimiter="\t")
        writer.writeheader()
        for index, smiles in enumerate(sorted(selected.values())):
            writer.writerow(
                {
                    "id": f"train-{index:06d}",
                    "smiles": smiles,
                    "source": "Benchmark development compound; not supplier availability",
                    "snapshot": audit["dataset_sha256"],
                }
            )
    report = {
        "dataset_sha256": digest(dataset_path),
        "library_sha256": digest(output),
        "records": len(selected),
        "partition": "train only",
        "labels_exposed": False,
        "test_or_query_family_overlap": 0,
        "availability_claim": False,
        "known_active_recovery_limit": "Exact test compounds are absent by construction; novel scaffold recovery requires generation, not nearest retrieval",
    }
    write_json(output.with_suffix(output.suffix + ".provenance.json"), report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    export_library(args.dataset, args.output)


if __name__ == "__main__":
    main()
