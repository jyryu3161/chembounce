"""Bounded, saved-byte review of routes outside the selected recommendation list.

This module deliberately has no chemistry imports. Its checks establish record
identity and correspondence, not reaction feasibility or route eligibility.
"""

from collections import Counter
import hashlib
import json


LIMIT = 16


def _need(condition, reason):
    if not condition:
        raise ValueError(reason)


def _list(value, name):
    _need(isinstance(value, list) and all(isinstance(x, dict) for x in value), f"invalid {name}")
    return value


def _same(actual, expected, name):
    _need(actual == expected, f"{name} differs")


def _binding(record, key, name):
    value = record.get(key)
    _need(
        isinstance(value, str) and len(value) == 64 and all(ch in "0123456789abcdef" for ch in value),
        f"{name} missing or malformed",
    )
    body = {field: item for field, item in record.items() if field != key}
    digest = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()
    _same(value, digest, name)
    return value


def _by_node(records, expected, name):
    ids = [item.get("node_id") for item in _list(records, name)]
    _need(all(isinstance(x, str) and x for x in ids), f"invalid {name} node ID")
    _need(len(ids) == len(set(ids)) and set(ids) == expected, f"{name} node coverage differs")
    return dict(zip(ids, records))


def _graph(record):
    graph = record.get("route")
    _need(isinstance(graph, dict), "route graph missing")
    _need(type(graph.get("schema_version")) is int and graph["schema_version"] == 1, "route graph schema unsupported")
    nodes = _list(graph.get("nodes"), "graph nodes")
    _need(len(nodes) <= 4096, "graph node review bound exceeded")
    ids = [node.get("id") for node in nodes]
    _need(all(isinstance(x, str) and x for x in ids) and len(ids) == len(set(ids)), "invalid graph node IDs")
    by_id = dict(zip(ids, nodes))
    reactions = {ident: node for ident, node in by_id.items() if "rule_id" in node}
    leaves = {ident: node for ident, node in by_id.items() if "building_block_id" in node}
    _need(reactions and graph.get("root") in reactions, "reaction root missing")
    _need(set(by_id) == set(reactions) | set(leaves) and not (set(reactions) & set(leaves)), "invalid node kind")
    _need(all(isinstance(node.get("rule_id"), str) and node["rule_id"] for node in reactions.values()),
          "reaction rule ID missing")
    _need(all(isinstance(node.get("building_block_id"), str) and node["building_block_id"] for node in leaves.values()),
          "leaf catalog ID missing")
    seen, visiting, pending = set(), set(), [(graph["root"], False)]
    while pending:
        ident, done = pending.pop()
        _need(isinstance(ident, str) and ident in by_id, "unknown graph dependency")
        if done:
            visiting.remove(ident)
            seen.add(ident)
            continue
        if ident in seen:
            continue
        _need(ident not in visiting, "cyclic route graph")
        visiting.add(ident)
        pending.append((ident, True))
        if ident in reactions:
            refs = reactions[ident].get("reactants")
            _need(isinstance(refs, list) and refs and all(isinstance(x, str) for x in refs), "invalid reactants")
            pending.extend((ref, False) for ref in reversed(refs))
    _same(seen, set(by_id), "root-reachable graph")
    verification = record.get("verification")
    _need(isinstance(verification, dict), "verification missing")
    root_product = reactions[graph["root"]].get("product")
    _need(isinstance(root_product, str) and root_product, "root product missing")
    _same(verification.get("product"), root_product, "verification product")
    verification_blocks = _by_node(verification.get("building_blocks"), set(leaves), "verification leaves")
    products = {}
    for ident, node in leaves.items():
        block = verification_blocks[ident]
        _same(block.get("id"), node.get("building_block_id"), "leaf catalog ID")
        _need(isinstance(block.get("smiles"), str) and block["smiles"], "leaf structure missing")
        products[ident] = block["smiles"]
    verified = _by_node(verification.get("steps"), set(reactions), "verification steps")
    for step in verification["steps"]:
        ident = step["node_id"]
        node = reactions[ident]
        refs = node.get("reactants", [])
        _need(all(ref in products for ref in refs), "verification dependency order differs")
        _same(step.get("rule_id"), node.get("rule_id"), "verification rule")
        _same(step.get("reactants"), [products[ref] for ref in refs], "verification reactants")
        _same(step.get("product"), node.get("product"), "verification step product")
        _need(isinstance(node.get("product"), str) and node["product"], "reaction product missing")
        products[ident] = node["product"]
    _same(products.get(graph["root"]), root_product, "root product")
    return graph, reactions, leaves, verification, verification_blocks, products, verified


def _hypotheses(prediction, expected, reactants, context, provenance):
    """Return display fields and a conservative status from retained raw beams."""
    status = prediction.get("status", "not_recorded")
    raw = prediction.get("raw")
    if raw is None:
        return [], "unassessed", "raw forward hypotheses not recorded"
    _need(isinstance(raw, dict), "raw forward evidence malformed")
    beams = raw.get("predictions")
    _need(isinstance(beams, list) and beams, "raw forward hypotheses missing")
    if "predictions" in prediction:
        _same(prediction["predictions"], beams, "raw and reported hypotheses")
    displayed, complete_exact, complete_other = [], [], []
    raw_input = raw.get("input")
    binding_issue = None
    if not isinstance(provenance, dict) or type(provenance.get("top_k")) is not int or provenance["top_k"] != 3:
        binding_issue = "bounded top-three predictor provenance not recorded"
    elif len(beams) != 3:
        binding_issue = "retained raw top-three beam count differs"
    elif not isinstance(raw_input, dict):
        binding_issue = "raw input or predictor provenance not recorded"
    else:
        supplied = context.get("supplied") if isinstance(context, dict) else None
        if not isinstance(supplied, dict):
            binding_issue = "step context not recorded"
        else:
            reagents = supplied.get("reagents", [])
            if (raw_input.get("reactants") != sorted(reactants)
                    or not isinstance(reagents, list)
                    or raw_input.get("reagents") != sorted(reagents)):
                binding_issue = "raw input reactants or reagents differ"
            elif raw_input.get("prompt") != (
                "REACTANT:" + ".".join(sorted(reactants)) + "REAGENT:" + ".".join(sorted(reagents))
            ):
                binding_issue = "raw prompt differs from recorded reactants and reagents"
            elif raw.get("predictor_binding_sha256") != hashlib.sha256(
                json.dumps(provenance, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
            ).hexdigest():
                binding_issue = "raw predictor provenance binding differs"
    for index, beam in enumerate(beams, 1):
        _need(isinstance(beam, dict) and type(beam.get("rank")) is int and beam["rank"] == index,
              "raw hypothesis rank/order differs")
        _need(type(beam.get("complete")) is bool, "raw hypothesis completion malformed")
        product = beam.get("product")
        _need(product is None or isinstance(product, str), "raw hypothesis product malformed")
        invalid_reason = beam.get("invalid_reason")
        _need(invalid_reason is None or isinstance(invalid_reason, str), "raw hypothesis invalid reason malformed")
        if invalid_reason and product is not None:
            binding_issue = "invalid raw hypothesis also claims a product"
        if beam["complete"] and isinstance(product, str) and product and not invalid_reason:
            (complete_exact if product == expected else complete_other).append(index)
        elif not beam["complete"]:
            _need(product is None, "incomplete hypothesis has a product")
        displayed.append(
            {
                "rank": index,
                "complete": beam["complete"],
                "valid_product": beam["complete"] and isinstance(product, str) and bool(product) and not invalid_reason,
                "expected_product_match": beam["complete"] and product == expected and not invalid_reason,
                "product": product,
                "raw_product": beam.get("raw_product"),
                "decoded_text": beam.get("original_decoded_text"),
                "invalid_reason": invalid_reason,
            }
        )
    if binding_issue:
        return displayed, "unassessed", binding_issue
    if (status == "matched" and prediction.get("matching_policy") == "exact"
            and complete_exact and type(prediction.get("matched_rank")) is int
            and prediction["matched_rank"] == min(complete_exact)):
        return displayed, "recorded_exact_match", None
    if (status == "mismatched" and prediction.get("matching_policy") == "exact"
            and prediction.get("matched_rank") is None
            and not complete_exact and len(complete_other) == len(beams)):
        return displayed, "recorded_complete_nonmatch", None
    if status in ("unknown", "partially_evaluated", "failed"):
        return displayed, "unassessed", prediction.get("reason") or status
    return displayed, "unassessed", "forward status contradicts retained raw hypotheses"


def _evidence(record, graph_data, item):
    graph, reactions, leaves, verification, _, products, _ = graph_data
    evidence = record.get("synthesis_evidence")
    if evidence is None:
        item["evidence_status"] = "not_recorded"
        item["reason"] = record.get("reason") or record.get("reasons") or "forward evidence not recorded"
        return
    _need(isinstance(evidence, dict), "route evidence malformed")
    _need(type(evidence.get("schema_version")) is int and evidence["schema_version"] == 2,
          "support schema unsupported")
    _same(evidence.get("profile"), "computational-route-support-v1", "support profile")
    audit = evidence.get("route_audit")
    _need(isinstance(audit, dict), "route audit missing")
    _same(audit.get("schema"), "concrete-route-audit-v1", "route audit schema")
    root_product = products[graph["root"]]
    for name, value in (("evidence product", evidence.get("product")), ("audit product", audit.get("product"))):
        _same(value, root_product, name)
    _same(audit.get("root"), graph["root"], "audit root")
    binding = _binding(audit, "route_binding_sha256", "audit route binding")
    _same(evidence.get("route_binding_sha256"), binding, "audit/support route binding")
    audit_blocks = _by_node(audit.get("building_blocks"), set(leaves), "audit leaves")
    for ident, block in audit_blocks.items():
        parent = _binding(block, "building_block_binding_sha256", "audit building-block binding")
        catalog = block.get("catalog_record")
        _need(isinstance(catalog, dict), "audit catalog record missing")
        _same(catalog.get("id"), leaves[ident].get("building_block_id"), "audit catalog ID")
        _same(block.get("smiles"), products[ident], "audit leaf structure")
        _same(catalog.get("smiles", products[ident]), products[ident], "catalog structure")
        item["sources"] = [source for source in item["sources"] if source["node_id"] != ident]
        item["sources"].append({"node_id": ident, "smiles": products[ident], "catalog_record": catalog,
                                "building_block_binding_sha256": parent})
    audits = _by_node(audit.get("steps"), set(reactions), "audit steps")
    support_records = _list(evidence.get("steps"), "support steps")
    support_ids = [step.get("node_id") for step in support_records]
    _need(all(isinstance(ident, str) and ident in reactions for ident in support_ids), "support step node invalid")
    _need(len(support_ids) == len(set(support_ids)), "support step node duplicated")
    audit_ids = [step["node_id"] for step in audit["steps"]]
    _same(audit_ids, [step["node_id"] for step in verification["steps"]], "verification/audit step order")
    _same(support_ids, [ident for ident in audit_ids if ident in support_ids], "audit/support relative order")
    complete_coverage = set(support_ids) == set(reactions)
    if not complete_coverage:
        item["issues"].append("support step coverage incomplete; missing steps remain unassessed")
    step_views = item["steps"]
    step_by_id = {step["node_id"]: step for step in step_views}
    parent_bindings = {ident: block["building_block_binding_sha256"] for ident, block in audit_blocks.items()}
    for audited in audit["steps"]:
        ident = audited["node_id"]
        node, audited = reactions[ident], audits[ident]
        expected = [products[ref] for ref in node["reactants"]]
        slots = _list(audited.get("reactant_slots"), "audit reactant slots")
        _same([s.get("node_id") for s in slots], node["reactants"], "audit reactant node order")
        _same([s.get("slot") for s in slots], list(range(len(expected))), "audit slot order")
        _same([s.get("smiles") for s in slots], expected, "audit ordered reactants")
        _same([s.get("parent_binding_sha256") for s in slots],
              [parent_bindings[ref] for ref in node["reactants"]], "audit parent bindings")
        _same(audited.get("rule_id"), node["rule_id"], "audit rule")
        _same(audited.get("product"), node["product"], "audit product")
        _binding(audited, "step_binding_sha256", "audit step binding")
        parent_bindings[ident] = audited["step_binding_sha256"]
    for step in support_records:
        ident = step["node_id"]
        node, audited = reactions[ident], audits[ident]
        expected = [products[ref] for ref in node["reactants"]]
        step_view = step_by_id[ident]
        try:
            _same(step.get("reactants"), expected, "support reactants")
            _same(step.get("product"), node["product"], "support product")
            _same(step.get("step_binding_sha256"), audited.get("step_binding_sha256"), "support step binding")
            prediction = step.get("forward_prediction")
            _need(isinstance(prediction, dict), "forward prediction missing")
        except (KeyError, TypeError, ValueError, AttributeError) as exc:
            item["issues"].append(f"Step {ident} support correspondence incomplete or inconsistent: {exc}")
            step_view["reason"] = "saved support correspondence incomplete"
            continue
        step_view.update({
                "node_id": ident,
                "rule_id": node["rule_id"],
                "reactants": expected,
                "expected_product": node["product"],
                "step_binding_sha256": step["step_binding_sha256"],
                "rule_replay": step.get("rule_replay"),
                "audit_status": {k: audited.get(k) for k in ("strict_eligible", "atom_origin_status", "stereo_status", "enumeration")},
                "forward_status": "unassessed",
                "prediction_status": prediction.get("status"),
                "reason_code": prediction.get("reason_code"),
                "reason": prediction.get("reason"),
                "hypotheses": [],
                "context": step.get("context"),
                "observations": step.get("observations"),
                "atom_sources": step.get("atom_sources"),
            })
        try:
            hypotheses, forward_status, reason = _hypotheses(
                prediction, node["product"], expected, step.get("context"), evidence.get("predictor_provenance")
            )
            step_view.update(hypotheses=hypotheses, forward_status=forward_status,
                             reason=reason or prediction.get("reason"))
            if reason in (
                "bounded top-three predictor provenance not recorded", "retained raw top-three beam count differs",
                "raw input or predictor provenance not recorded", "step context not recorded",
                "raw input reactants or reagents differ", "raw predictor provenance binding differs",
                "raw prompt differs from recorded reactants and reagents",
                "invalid raw hypothesis also claims a product",
                "forward status contradicts retained raw hypotheses",
            ):
                item["issues"].append(f"Step {ident} raw binding or beam inventory incomplete: {reason}")
        except (KeyError, TypeError, ValueError, AttributeError) as exc:
            item["issues"].append(f"Step {ident} raw evidence incomplete or inconsistent: {exc}")
    item["route_binding_sha256"] = binding
    item["predictor_provenance"] = evidence.get("predictor_provenance")
    item["evidence_status"] = "recorded_bound_steps" if not item["issues"] else "incomplete_or_inconsistent"
    item["forward_review_status"] = (
        "all_recorded_exact_matches"
        if not item["issues"] and all(s["forward_status"] == "recorded_exact_match" for s in step_views)
        else "contains_recorded_nonmatch"
        if not item["issues"] and any(s["forward_status"] == "recorded_complete_nonmatch" for s in step_views)
        else "incomplete_or_unknown"
    )


def _one(record):
    item = {
        "route_id": record.get("route_id"),
        "selection_status": record.get("selection_status"),
        "reason": record.get("reason") or record.get("reasons"),
        "recorded_comparison": {
            "design": record.get("design"),
            "hit_tanimoto": record.get("tanimoto", {}).get("hit_tanimoto")
            if isinstance(record.get("tanimoto"), dict) else None,
            "design_tanimoto": record.get("tanimoto", {}).get("design_tanimoto")
            if isinstance(record.get("tanimoto"), dict) else None,
            "polar_recall_upper_bound": record.get("polar_capacity", {}).get("recall_upper_bound")
            if isinstance(record.get("polar_capacity"), dict) else None,
            "polar_capacity_status": record.get("polar_capacity", {}).get("status")
            if isinstance(record.get("polar_capacity"), dict) else None,
            "qualification_status": record.get("qualification", {}).get("status")
            if isinstance(record.get("qualification"), dict) else None,
            "qualification_reasons": record.get("qualification", {}).get("reasons")
            if isinstance(record.get("qualification"), dict) else None,
        },
        "identity_status": "unassessed",
        "evidence_status": "unassessed",
        "forward_review_status": "unassessed",
        "sources": [],
        "steps": [],
        "issues": [],
    }
    try:
        _need(isinstance(item["route_id"], str) and item["route_id"], "route ID missing or malformed")
        data = _graph(record)
        graph, reactions, leaves, _, blocks, products, _ = data
        item["identity_status"] = "consistent_graph_and_verification"
        item["product"] = products[graph["root"]]
        item["sources"] = [
            {"node_id": ident, "smiles": products[ident], "catalog_record": blocks[ident]}
            for ident in leaves
        ]
        item["intermediates"] = [
            {"node_id": ident, "product": node["product"]}
            for ident, node in reactions.items()
            if ident != graph["root"]
        ]
        item["steps"] = [
            {
                "node_id": step["node_id"],
                "rule_id": step["rule_id"],
                "reactants": step["reactants"],
                "expected_product": step["product"],
                "forward_status": "unassessed",
                "reason": "stepwise forward evidence not recorded",
                "hypotheses": [],
            }
            for step in data[-1].values()
        ]
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        item["issues"].append(str(exc))
        item["identity_status"] = "inconsistent_or_incomplete"
        return item
    try:
        _evidence(record, data, item)
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        item["issues"].append(str(exc))
        item["evidence_status"] = "incomplete_or_inconsistent"
        item["forward_review_status"] = "unassessed"
        for step in item["steps"]:
            step.update(forward_status="unassessed", reason="evidence correspondence incomplete",
                        hypotheses=[])
    return item


def build(routes, selected_ids):
    """Review at most LIMIT non-primary records in their original saved order."""
    if not isinstance(routes, list):
        return {"total": 0, "displayed": 0, "omitted": 0, "routes": [], "issue": "routes list unavailable"}
    selected = {ident for ident in selected_ids if isinstance(ident, str)}
    unselected = [route for route in routes if not isinstance(route, dict) or not isinstance(route.get("route_id"), str) or route["route_id"] not in selected]
    seen = Counter(route.get("route_id") for route in unselected if isinstance(route, dict) and isinstance(route.get("route_id"), str))
    displayed = []
    for route in unselected[:LIMIT]:
        item = _one(route) if isinstance(route, dict) else {"route_id": None, "issues": ["route record malformed"]}
        if isinstance(route, dict) and isinstance(route.get("route_id"), str) and seen[route["route_id"]] > 1:
            item["issues"].append("duplicate unselected route ID")
            item["identity_status"] = "inconsistent_or_incomplete"
            item["forward_review_status"] = "unassessed"
        displayed.append(item)
    return {"total": len(unselected), "displayed": len(displayed), "omitted": len(unselected) - len(displayed), "routes": displayed}
