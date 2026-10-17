"""Lightweight CLI controls; importing these never loads the model or 3D runtime."""

import copy
from argparse import BooleanOptionalAction


def add_expansion_options(parser):
    parser.add_argument(
        "--recommend",
        action="store_true",
        help="Request qualified-selective-engineering-v7: original-hit Tanimoto/shape/pharmacophore floors .35/.70/.50, polar recall .75, changed scaffold and computational route support; supported explicit overrides remain custom",
    )
    parser.add_argument(
        "--forward-model-path",
        help="Explicit local pinned ReactionT5 bundle required by the named recommendation profiles",
    )
    parser.add_argument(
        "--recommendation-profile",
        choices=(
            "qualified-development-v1",
            "qualified-selective-engineering-v2",
            "qualified-selective-engineering-v3",
            "qualified-selective-engineering-v4",
            "qualified-selective-engineering-v5",
            "qualified-selective-engineering-v6",
            "qualified-selective-engineering-v7",
        ),
        help="Explicit opt-in profile: permissive development fit or selective engineering constraints (not fitted); fixed three-seed scoring and computational route support, no activity guarantee",
    )
    parser.add_argument(
        "--landmark-pose-proposals",
        action=BooleanOptionalAction,
        default=None,
        help="Explicit bounded rigid proposals for configured mandatory landmarks; requires fixed-pose ScoreMol, changes search and remains custom_unvalidated. Defaults off; --no-landmark-pose-proposals disables a retained setting",
    )
    parser.add_argument(
        "--scoring-objective",
        choices=("joint-v1", "original-hit-v1"),
        help="Conformer selection objective; original-hit-v1 selects against the hit before the secondary design overlay",
    )
    parser.add_argument(
        "--ranking-policy",
        choices=("joint-v1", "original-hit-v1"),
        help="Original-hit ranking uses hit-only dominance and robust seed minima; defaults to the scoring objective",
    )
    parser.add_argument(
        "--recommendation-mode",
        choices=("exploratory", "qualified"),
        help="qualified enforces qualified-v1 on final products using explicit positive hit similarity minima; exploratory results are not recommendations",
    )
    parser.add_argument(
        "--route-search-policy",
        choices=("original-hit-v1", "design-v1"),
        help="Root-product search objective; qualified defaults to original-hit-v1; design-v1 is an explicit search ablation with unchanged final acceptance",
    )
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
        "--synthesis-max-designs", type=int, help="Maximum eligible designs in bounded refill mode"
    )
    group.add_argument(
        "--synthesis-max-design-attempts",
        type=int,
        help="Bound raw design attempts including invalid/duplicate entries; enables refill for exploratory runs",
    )
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
    if getattr(args, "recommend", False):
        from .profiles import request_recommendations

        # Inspect both declarations before ordinary CLI precedence can hide one.
        config = request_recommendations(config, getattr(args, "recommendation_profile", None))
    if getattr(args, "forward_model_path", None) is not None:
        from pathlib import Path

        config.setdefault("synthesis_evidence", {}).setdefault("predictor", {}).update(
            kind="reactiont5", model_path=str(Path(args.forward_model_path).resolve())
        )
    if getattr(args, "recommendation_profile", None) is not None:
        config["recommendation_profile"] = args.recommendation_profile
    if getattr(args, "landmark_pose_proposals", None) is not None:
        config.setdefault("similarity_constraints", {})["landmark_pose_proposals"] = (
            args.landmark_pose_proposals
        )
    if getattr(args, "scoring_objective", None) is not None:
        config.setdefault("conformers", {})["scoring_objective"] = args.scoring_objective
    if getattr(args, "ranking_policy", None) is not None:
        config["ranking_policy"] = args.ranking_policy
    if getattr(args, "recommendation_mode", None) is not None:
        config.setdefault("recommendation", {}).update(mode=args.recommendation_mode)
    if getattr(args, "route_search_policy", None) is not None:
        config["route_search_policy"] = args.route_search_policy
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
        ("synthesis_max_designs", "max_designs"),
        ("synthesis_max_design_attempts", "max_design_attempts"),
    ):
        if (value := getattr(args, attr, None)) is not None:
            config[key] = value
    if getattr(args, "include_run", None):
        from pathlib import Path

        config.setdefault("carry_forward", []).extend(str(Path(p).resolve()) for p in args.include_run)
    return config


def has_synthesis_overrides(args):
    return bool({key: value for key, value in apply_expansion_options({}, args).items() if key != "protocol"})
