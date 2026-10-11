"""Resolve LLM proposals to catalog routes, then verify by independent replay."""

import copy
import itertools
import json
from dataclasses import dataclass, field

from rdkit import Chem, DataStructs

from .common import FP, canonical, molecule, positive


@dataclass
class SearchLimits:
    building_block_top_k: int = 8
    beam_width: int = 24
    max_combinations: int = 4096
    reaction_product_limit: int = 64
    max_steps: int = 12

    def __post_init__(self):
        for key, value in vars(self).items():
            positive(value, key)


@dataclass
class SearchReport:
    combinations: int = 0
    reaction_truncations: int = 0
    combination_truncations: int = 0
    beam_truncations: int = 0
    failures: list = field(default_factory=list)


def _untag(value, tag):
    if not isinstance(value, str):
        raise ValueError("Reaction and building-block fields must be strings")
    return value.strip().removeprefix(f"<{tag}>").removesuffix(f"</{tag}>").strip()


def parse_response(text):
    """Accept one JSON object (optionally a Markdown fence); never eval model text."""
    if isinstance(text, dict):
        data = copy.deepcopy(text)
    else:
        if not isinstance(text, str) or len(text) > 250_000:
            raise ValueError("Invalid or oversized model response")
        text = text.strip()
        if text.startswith("```json") and text.endswith("```"):
            text = text[7:-3].strip()
        elif text.startswith("```") and text.endswith("```"):
            text = text[3:-3].strip()
        data = json.loads(text)
    if (
        not isinstance(data, dict)
        or not isinstance(data.get("reactions"), list)
        or not isinstance(data.get("building_blocks"), list)
    ):
        raise ValueError("Model response requires reactions and building_blocks lists")
    return data


def _topology(nodes, root):
    if not isinstance(root, str) or not root:
        raise ValueError("Route root must be a nonempty node ID")
    if not isinstance(nodes, list) or not nodes or len(nodes) > 256:
        raise ValueError("Route requires 1–256 nodes")
    lookup = {}
    for node in nodes:
        if not isinstance(node, dict) or not isinstance(node.get("id"), str) or not node["id"]:
            raise ValueError("Each route node requires a nonempty string ID")
        if node["id"] in lookup:
            raise ValueError("Duplicate route node ID")
        lookup[node["id"]] = node
    visiting, visited, order = set(), set(), []

    def visit(ident):
        if ident in visiting:
            raise ValueError("Cyclic reaction route")
        if ident in visited:
            return
        if ident not in lookup:
            raise ValueError(f"Missing route node: {ident}")
        visiting.add(ident)
        node = lookup[ident]
        is_bb = "building_block_id" in node
        if is_bb == ("rule_id" in node):
            raise ValueError("A route node must be exactly one building block or reaction")
        if is_bb:
            if set(node) != {"id", "building_block_id"}:
                raise ValueError("Unexpected building-block node fields")
        else:
            if set(node) != {"id", "rule_id", "reactants", "product"}:
                raise ValueError("Unexpected reaction node fields")
            if not isinstance(node["reactants"], list) or not all(
                isinstance(x, str) for x in node["reactants"]
            ):
                raise ValueError("Reactants must be a list of node IDs")
            for ref in node["reactants"]:
                visit(ref)
        visiting.remove(ident)
        visited.add(ident)
        order.append(node)

    visit(root)
    if visited != set(lookup):
        raise ValueError("Route contains disconnected or unused nodes")
    return order


def verify_route(route, catalog, rules, product_limit=64):
    """No model-derived structure is accepted without exact isomeric replay identity."""
    if route.get("schema_version") != 1:
        raise ValueError("Route requires schema_version=1")
    order = _topology(route.get("nodes"), route.get("root"))
    if "rule_id" not in order[-1]:
        raise ValueError("A validated synthetic route must contain a reaction")
    products, evidence, steps, truncated = {}, [], [], False
    for node in order:
        ident = node["id"]
        if "building_block_id" in node:
            block = catalog.get(node["building_block_id"])
            products[ident] = block.smiles
            evidence.append({"node_id": ident, **catalog.evidence(block.id)})
        else:
            reactants = tuple(products[ref] for ref in node["reactants"])
            choices, capped = rules.apply(node["rule_id"], reactants, product_limit)
            target = canonical(node["product"])
            if target not in choices:
                raise ValueError(
                    f"Reaction {ident} does not regenerate its claimed product (including stereochemistry)"
                )
            products[ident] = target
            rule = rules.get(node["rule_id"])
            steps.append(
                {
                    "node_id": ident,
                    "rule_id": rule.id,
                    "smarts": rule.smarts,
                    "source": rule.source,
                    "unmapped_product_template_atoms": rule.unmapped_product_template_atoms,
                    "template_validation_warnings": rule.validation_warnings,
                    "reactants": list(reactants),
                    "product": target,
                    "enumeration_truncated": capped,
                }
            )
            truncated |= capped
    return {
        "status": "catalog_and_rule_verified",
        "product": products[route["root"]],
        "building_blocks": evidence,
        "steps": steps,
        "catalog_sha256": catalog.sha256,
        "rules_sha256": rules.sha256,
        "enumeration_truncated": truncated,
        "experimental_validation": False,
        "reagents_and_conditions_validated": False,
        "has_implicit_atom_sources": any(s["unmapped_product_template_atoms"] for s in steps),
    }


def _proposal_graph(proposal, rules, limits):
    data = parse_response(proposal)
    if not 1 <= len(data["reactions"]) <= limits.max_steps:
        raise ValueError("Model reaction count is outside configured limits")
    bbs = {canonical(_untag(s, "bb")) for s in data["building_blocks"]}
    products = {}
    for index, row in enumerate(data["reactions"]):
        if not isinstance(row, dict):
            raise ValueError("Each proposed reaction must be an object")
        product = canonical(row["product"])
        if product in products or product in bbs:
            raise ValueError("Ambiguous intermediate/building-block identity")
        rule_id = rules.resolve(row["reaction_template"])
        reactants = tuple(canonical(_untag(s, "bb")) for s in row["reactants"])
        if len(reactants) != len(rules.get(rule_id).roles):
            raise ValueError("Proposed reactant count does not match rule")
        products[product] = {"id": f"r{index}", "rule_id": rule_id, "reactants": reactants}
    used = {s for p in products.values() for s in p["reactants"]}
    roots = set(products) - used
    if len(roots) != 1:
        raise ValueError("Proposal must have exactly one terminal product")
    root = roots.pop()
    visiting, visited = set(), set()

    def visit(s):
        if s in visiting:
            raise ValueError("Cyclic model proposal")
        if s in visited:
            return
        visiting.add(s)
        if s in products:
            for reactant in products[s]["reactants"]:
                visit(reactant)
        elif s not in bbs:
            raise ValueError("Proposed leaf is not a declared building block")
        visiting.remove(s)
        visited.add(s)

    visit(root)
    if visited != set(products) | bbs:
        raise ValueError("Proposal contains unused building blocks or reactions")
    return products, bbs, root


def _prune(states, target, width, report):
    """Bound intermediate search; preserve different routes to the same molecule."""
    fp = FP.GetFingerprint(molecule(target))
    unique = {}
    for state in states:
        key = json.dumps(state, sort_keys=True)
        unique.setdefault(key, state)
    ranked = sorted(
        unique.values(),
        key=lambda s: (
            -DataStructs.TanimotoSimilarity(fp, FP.GetFingerprint(molecule(s["product"]))),
            json.dumps(s, sort_keys=True),
        ),
    )
    if len(ranked) > width:
        report.beam_truncations += 1
    return ranked[:width]


def resolve_proposal(proposal, catalog, rules, limits=None):
    """Project the proposed DAG onto real blocks and rules; keep regenerated products.

    Invalid rules/SMILES are rejected. Reactant permutations are enumerated rather
    than trusting LLM slot order. Shared intermediates/blocks must remain consistent.
    """
    limits = limits or SearchLimits()
    report = SearchReport()
    graph, bbs, root = _proposal_graph(proposal, rules, limits)
    bb_ids = {s: f"b{i}" for i, s in enumerate(sorted(bbs))}
    memo = {}

    def solve(smiles):
        if smiles in memo:
            return memo[smiles]
        step = graph[smiles]
        rule = rules.get(step["rule_id"])
        states = []
        count = 0
        stop = False
        for reactants in sorted(set(itertools.permutations(step["reactants"]))):
            options = []
            for reactant, role in zip(reactants, rule.roles):
                if reactant in graph:
                    options.append(solve(reactant))
                else:
                    options.append(
                        [
                            {
                                "root": bb_ids[reactant],
                                "product": bb.smiles,
                                "nodes": {
                                    bb_ids[reactant]: {"id": bb_ids[reactant], "building_block_id": bb.id}
                                },
                            }
                            for bb in catalog.nearest(reactant, role, limits.building_block_top_k)
                        ]
                    )
            for combination in itertools.product(*options):
                if count >= limits.max_combinations:
                    report.combination_truncations += 1
                    stop = True
                    break
                count += 1
                report.combinations += 1
                nodes, consistent = {}, True
                for item in combination:
                    for key, node in item["nodes"].items():
                        if key in nodes and node != nodes[key]:
                            consistent = False
                        nodes[key] = node
                if not consistent:
                    continue
                generated, capped = rules.apply(
                    rule.id, tuple(c["product"] for c in combination), limits.reaction_product_limit
                )
                report.reaction_truncations += int(capped)
                for product in generated:
                    new_nodes = nodes | {
                        step["id"]: {
                            "id": step["id"],
                            "rule_id": rule.id,
                            "reactants": [c["root"] for c in combination],
                            "product": product,
                        }
                    }
                    states.append({"root": step["id"], "product": product, "nodes": new_nodes})
                if len(states) >= 4 * limits.beam_width:
                    states = _prune(states, smiles, limits.beam_width, report)
            if stop:
                break
        memo[smiles] = _prune(states, smiles, limits.beam_width, report)
        if not memo[smiles]:
            report.failures.append(
                {
                    "model_intermediate": smiles,
                    "rule_id": rule.id,
                    "reason": "No valid regenerated products within catalog and search limits",
                }
            )
        return memo[smiles]

    routes = []
    for state in solve(root):
        route = {"schema_version": 1, "root": state["root"], "nodes": list(state["nodes"].values())}
        evidence = verify_route(route, catalog, rules, limits.reaction_product_limit)
        routes.append(
            {
                "route": route,
                "verification": evidence,
                "model_claimed_product": root,
                "exact_model_product": evidence["product"] == root,
            }
        )
    return routes, report


def _rebuild(route, catalog, rules, target, limits, report):
    order = _topology(route["nodes"], route["root"])
    states = [{"nodes": {}, "values": {}}]
    for node in order:
        expanded = []
        for state in states:
            if "building_block_id" in node:
                products = (catalog.get(node["building_block_id"]).smiles,)
            else:
                products, capped = rules.apply(
                    node["rule_id"],
                    tuple(state["values"][r] for r in node["reactants"]),
                    limits.reaction_product_limit,
                )
                report.reaction_truncations += int(capped)
            for product in products:
                new_node = node if "building_block_id" in node else node | {"product": product}
                expanded.append(
                    {
                        "nodes": state["nodes"] | {node["id"]: new_node},
                        "values": state["values"] | {node["id"]: product},
                        "product": product,
                    }
                )
        # Intermediate products are compared to their own former product, not the final hit.
        comparison = node.get("product", target)
        states = _prune(expanded, comparison, limits.beam_width, report)
        if not states:
            break
    return [{"schema_version": 1, "root": route["root"], "nodes": list(s["nodes"].values())} for s in states]


def expand_route(
    route, catalog, rules, limits=None, rounds=1, max_routes=64, locked_nodes=(), locked_building_block_ids=()
):
    """Replace compatible catalog blocks and replay the entire affected DAG.

    Each round changes one leaf per parent, allowing bounded combinatorial expansion.
    All returned pathways pass the independent verifier; protected leaves stay fixed.
    """
    limits = limits or SearchLimits()
    positive(rounds, "expansion_rounds")
    positive(max_routes, "expansion_max_routes")
    initial = verify_route(route, catalog, rules, limits.reaction_product_limit)
    leaves = {n["id"] for n in route["nodes"] if "building_block_id" in n}
    if not set(locked_nodes) <= leaves:
        raise ValueError("Locked expansion nodes must be building-block node IDs")
    locked_smiles = {catalog.get(ident).smiles for ident in locked_building_block_ids}
    report, results = SearchReport(), []
    frontier = [route]
    seen = {json.dumps(route, sort_keys=True)}
    for _ in range(rounds):
        next_frontier = []
        for parent in frontier:
            nodes = {n["id"]: n for n in parent["nodes"]}
            for node in parent["nodes"]:
                if "building_block_id" not in node or node["id"] in locked_nodes:
                    continue
                requirements = []
                for step in parent["nodes"]:
                    if "rule_id" in step:
                        for slot, ref in enumerate(step["reactants"]):
                            if ref == node["id"]:
                                requirements.append(rules.get(step["rule_id"]).roles[slot])
                old = catalog.get(node["building_block_id"])
                if old.smiles in locked_smiles:
                    continue
                role_queries = [Chem.MolFromSmarts(role) for role in requirements]
                for bb in catalog.nearest(old.smiles, requirements[0], limits.building_block_top_k):
                    if bb.smiles == old.smiles or not all(
                        catalog.mols[bb.smiles].HasSubstructMatch(query, useChirality=True)
                        for query in role_queries
                    ):
                        continue
                    if report.combinations >= limits.max_combinations:
                        report.combination_truncations += 1
                        return results, report
                    report.combinations += 1
                    substituted = {
                        "schema_version": 1,
                        "root": parent["root"],
                        "nodes": list((nodes | {node["id"]: node | {"building_block_id": bb.id}}).values()),
                    }
                    for rebuilt in _rebuild(substituted, catalog, rules, initial["product"], limits, report):
                        key = json.dumps(rebuilt, sort_keys=True)
                        if key in seen:
                            continue
                        seen.add(key)
                        evidence = verify_route(rebuilt, catalog, rules, limits.reaction_product_limit)
                        results.append(
                            {"route": rebuilt, "verification": evidence, "expanded_from": initial["product"]}
                        )
                        next_frontier.append(rebuilt)
                        if len(results) >= max_routes:
                            report.beam_truncations += 1
                            return results, report
        frontier = next_frontier
        if not frontier:
            break
    return results, report
