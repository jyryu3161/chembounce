"""Lightweight CLI controls; importing these never loads the model or 3D runtime."""

import copy


def add_expansion_options(parser):
    group = parser.add_argument_group("building-block expansion")
    group.add_argument(
        "--expansion-top-k", type=int, help="Alternative blocks per mutable leaf (overrides config)"
    )
    group.add_argument("--expansion-rounds", type=int, help="Substitution rounds; 0 disables expansion")
    group.add_argument(
        "--expansion-max-routes", type=int, help="Retained expansion routes per starting route"
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
    overrides = expansion_overrides(args)
    if overrides:
        config.setdefault("expansion", {}).update(overrides)
    return config
