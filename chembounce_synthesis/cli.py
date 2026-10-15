"""CLI for existing designs or ChemBounce's tab-separated output."""

import argparse
import csv
import logging

from .workflow import load_config, run
from .options import add_expansion_options, apply_expansion_options
from .continuation import continuation_inputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hit", help="Original hit SMILES (inherited with --continue-from)")
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--design", action="append", help="Design SMILES; repeat for several designs")
    source.add_argument("--designs-tsv", help="ChemBounce overall_result.txt or TSV with a smiles column")
    parser.add_argument("--config", help="Synthesis workflow JSON configuration")
    parser.add_argument(
        "--continue-from", help="Expand a previous run using its inputs, routes and saved model proposals"
    )
    parser.add_argument("--output", required=True, help="Fresh output directory")
    add_expansion_options(parser)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        inherited_hit, inherited_designs, config = None, None, None
        if args.continue_from:
            inherited_hit, inherited_designs, config = continuation_inputs(args.continue_from)
        if args.config:
            config = load_config(args.config)
            if args.continue_from:
                from pathlib import Path

                config.setdefault("carry_forward", []).append(str(Path(args.continue_from).resolve()))
        if config is None or (not args.hit and not inherited_hit):
            raise ValueError("Provide --hit, --design/--designs-tsv and --config, or --continue-from")
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
        if designs is None:
            designs = inherited_designs
        if designs is None:
            raise ValueError("Provide --design/--designs-tsv or --continue-from")
        rows = run(args.hit or inherited_hit, designs, apply_expansion_options(config, args), args.output)
    except (ValueError, OSError, KeyError, TypeError, RuntimeError) as exc:
        logging.error("%s", exc)
        return 2
    print(f"Verified and ranked products: {len(rows)}; results: {args.output}/results.tsv")
    return 0 if rows else 3


if __name__ == "__main__":
    raise SystemExit(main())
