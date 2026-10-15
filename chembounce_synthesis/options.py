"""Lightweight CLI controls; importing these never loads the model or 3D runtime."""

import copy


def add_expansion_options(parser):
    parser.add_argument(
        "--protocol",
        choices=("legacy", "research-v1"),
        help="Versioned computational defaults; explicit settings override them. research-v1 enables diverse search and uncertainty-aware bounded refinement",
    )
    group = parser.add_argument_group("building-block expansion")
    group.add_argument(
        "--expansion-top-k", type=int, help="Alternative blocks per mutable leaf (overrides config)"
    )
    group.add_argument("--expansion-rounds", type=int, help="Substitution rounds; 0 disables expansion")
    group.add_argument(
        "--expansion-max-routes", type=int, help="Retained expansion routes per starting route"
    )
    group = parser.add_argument_group("synthesis reevaluation and accumulated candidates")
    group.add_argument(
        "--refine-top-n", type=int, help="Reevaluate the top N screened products; 0 disables refinement"
    )
    group.add_argument("--refine-conformers", type=int, help="Conformers per refinement seed (default 24)")
    group.add_argument(
        "--refine-seeds", type=int, nargs="+", help="Distinct refinement seeds (default 42 43 44)"
    )
    group.add_argument("--synthesis-max-scored", type=int, help="Maximum screened product/design pairs")
    group.add_argument("--synthesis-max-results", type=int, help="Maximum final products written")
    group.add_argument(
        "--include-run",
        action="append",
        help="Retain and reverify routes from a previous synthesis directory; repeatable",
    )


def expansion_overrides(args):
    return {
        key: value
        for attr, key in (
            ("expansion_top_k", "building_block_top_k"),
            ("expansion_rounds", "rounds"),
            ("expansion_max_routes", "max_routes"),
        )
        if (value := getattr(args, attr, None)) is not None
    }


def apply_expansion_options(config, args):
    config = copy.deepcopy(config)
    if getattr(args, "protocol", None) is not None:
        config["protocol"] = args.protocol
    overrides = expansion_overrides(args)
    if overrides:
        config.setdefault("expansion", {}).update(overrides)
    for attr, key in (
        ("refine_top_n", "top_n"),
        ("refine_conformers", "num_conformers"),
        ("refine_seeds", "seeds"),
    ):
        if (value := getattr(args, attr, None)) is not None:
            config.setdefault("refinement", {})[key] = value
    for attr, key in (
        ("synthesis_max_scored", "max_scored_products"),
        ("synthesis_max_results", "max_results"),
    ):
        if (value := getattr(args, attr, None)) is not None:
            config[key] = value
    if getattr(args, "include_run", None):
        from pathlib import Path

        config.setdefault("carry_forward", []).extend(str(Path(p).resolve()) for p in args.include_run)
    return config


def has_synthesis_overrides(args):
    return bool({key: value for key, value in apply_expansion_options({}, args).items() if key != "protocol"})
