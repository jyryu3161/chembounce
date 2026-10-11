"""CLI for existing designs or ChemBounce's tab-separated output."""

import argparse
import csv
import logging

from .workflow import load_config, run
from .options import add_expansion_options, apply_expansion_options


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hit", required=True, help="Original hit SMILES")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--design", action="append", help="Design SMILES; repeat for several designs")
    source.add_argument("--designs-tsv", help="ChemBounce overall_result.txt or TSV with a smiles column")
    parser.add_argument("--config", required=True, help="Synthesis workflow JSON configuration")
    parser.add_argument("--output", required=True, help="Fresh output directory")
    add_expansion_options(parser)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        designs = args.design
        if args.designs_tsv:
            with open(args.designs_tsv, newline="") as stream:
                reader = csv.DictReader(stream, delimiter="\t")
                column = next(
                    (c for c in ("Final structure", "smiles") if c in (reader.fieldnames or [])), None
                )
                if column is None:
                    raise ValueError("Design TSV requires a 'Final structure' or 'smiles' column")
                designs = [row[column] for row in reader]
        rows = run(args.hit, designs, apply_expansion_options(load_config(args.config), args), args.output)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        logging.error("%s", exc)
        return 2
    print(f"Verified and ranked products: {len(rows)}; results: {args.output}/results.tsv")
    return 0 if rows else 3


if __name__ == "__main__":
    raise SystemExit(main())
