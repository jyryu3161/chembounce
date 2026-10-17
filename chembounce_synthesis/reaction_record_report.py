"""Standalone optional reaction-record CLI and saved-route-bound evidence reports.

No catalog, reaction execution, model or scoring imports. Rendering only displays
captured reports; newly verified matches always go through raw bundle replay.
"""

import argparse
import html
import json
from pathlib import Path
import shutil
import sys

from . import reaction_records as records


def _same(left, right):
    a, b = records.smiles_identity(left), records.smiles_identity(right)
    return a["key"] is not None and a["key"] == b["key"]


def _saved_queries(run_dir, route_id, contexts):
    """Join saved final result, graph and verification; this does not replay chemistry."""
    raw_routes = records.read_local(run_dir, "routes.json")
    raw_results = records.read_local(run_dir, "results.json")
    routes_doc, results_doc = records.json_read(raw_routes), records.json_read(raw_results)
    records.require(
        isinstance(routes_doc, dict)
        and isinstance(results_doc, dict)
        and isinstance(routes_doc.get("routes"), list)
        and isinstance(results_doc.get("results"), list)
        and all(isinstance(row, dict) for row in routes_doc["routes"] + results_doc["results"]),
        "saved_documents",
        "Expected routes/results objects containing record lists",
    )
    routes = [r for r in routes_doc.get("routes", []) if r.get("route_id") == route_id]
    results = [r for r in results_doc.get("results", []) if r.get("route_id") == route_id]
    records.require(
        len(routes) == len(results) == 1,
        "saved_route_selection",
        "Select exactly one saved final result and matching route",
    )
    saved, result = routes[0], results[0]
    graph, verification = saved.get("route", {}), saved.get("verification", {})
    records.require(
        isinstance(graph, dict) and isinstance(verification, dict),
        "saved_documents",
        "Graph and verification must be objects",
    )
    nodes = graph.get("nodes", [])
    records.require(
        isinstance(nodes, list) and 0 < len(nodes) <= 256, "saved_graph", "Bounded nodes required"
    )
    by_id = {}
    for node in nodes:
        records.require(isinstance(node, dict), "saved_graph", "Nodes must be objects")
        ident = node.get("id")
        records.require(
            isinstance(ident, str) and ident and ident not in by_id, "saved_graph", "Unique node IDs required"
        )
        records.require(
            ("building_block_id" in node) != ("rule_id" in node),
            "saved_graph",
            "Expected building block or reaction node",
        )
        by_id[ident] = node
    root = graph.get("root")
    records.require(
        isinstance(root, str) and root in by_id and "rule_id" in by_id[root],
        "saved_graph",
        "Reaction root required",
    )
    root_product = by_id[root].get("product")
    records.require(
        _same(root_product, result.get("smiles")) and _same(root_product, verification.get("product")),
        "saved_final_identity",
        "Final result, graph and verification product must match exactly",
    )
    steps = verification.get("steps", [])
    records.require(
        isinstance(steps, list)
        and all(isinstance(step, dict) for step in steps)
        and all(isinstance(step.get("node_id"), str) for step in steps),
        "saved_steps",
        "Step objects with node IDs required",
    )
    steps_by_id = {s.get("node_id"): s for s in steps}
    reactions = {i for i, n in by_id.items() if "rule_id" in n}
    records.require(
        len(steps_by_id) == len(steps) and set(steps_by_id) == reactions,
        "saved_steps",
        "Every reaction must have exactly one recorded step",
    )
    records.require(
        isinstance(contexts, dict) and not set(contexts) - reactions,
        "query_contexts",
        "Context keys must be saved reaction IDs",
    )
    # Bind saved leaves through their recorded catalog ID. No live catalog access.
    blocks = verification.get("building_blocks", [])
    records.require(
        isinstance(blocks, list) and all(isinstance(block, dict) for block in blocks),
        "saved_blocks",
        "Recorded leaf objects required",
    )
    block_lookup, explicit_nodes, optional_leaf_fields_absent = {}, set(), []
    graph_leaves = {
        ident: node["building_block_id"] for ident, node in by_id.items() if "building_block_id" in node
    }
    records.require(
        all(isinstance(bid, str) and bid for bid in graph_leaves.values()),
        "saved_blocks",
        "Graph leaf catalog IDs must be nonempty strings",
    )
    for index, block in enumerate(blocks):
        bid, smiles = block.get("id"), block.get("smiles")
        records.require(
            isinstance(bid, str) and bid in graph_leaves.values() and isinstance(smiles, str),
            "saved_blocks",
            "Recorded leaf identity required",
        )
        records.require(
            bid not in block_lookup or _same(block_lookup[bid], smiles),
            "saved_blocks",
            "Conflicting saved catalog IDs",
        )
        if "node_id" in block:
            node_id = block["node_id"]
            records.require(
                isinstance(node_id, str)
                and node_id in graph_leaves
                and node_id not in explicit_nodes
                and graph_leaves[node_id] == bid,
                "saved_blocks",
                "Verification leaf node/catalog binding differs",
            )
            explicit_nodes.add(node_id)
        else:
            optional_leaf_fields_absent.append(f"verification/building_blocks/{index}/node_id")
        block_lookup[bid] = smiles
    # Strict support also records per-node leaves; usable when legacy verification omits leaf records.
    support = saved.get("synthesis_evidence", {})
    records.require(
        isinstance(support, dict) and isinstance(support.get("route_audit", {}), dict),
        "saved_blocks",
        "Saved support/audit must be objects",
    )
    audit_blocks = support.get("route_audit", {}).get("building_blocks", [])
    records.require(
        isinstance(audit_blocks, list) and all(isinstance(block, dict) for block in audit_blocks),
        "saved_blocks",
        "Audited leaf objects required",
    )
    per_node = {}
    for block in audit_blocks:
        ident = block.get("node_id")
        records.require(
            isinstance(ident, str)
            and ident not in per_node
            and ident in graph_leaves
            and isinstance(block.get("catalog_record"), dict)
            and block["catalog_record"].get("id") == by_id[ident].get("building_block_id"),
            "saved_blocks",
            "Audited leaf/node binding differs",
        )
        if "smiles" in block["catalog_record"]:
            records.require(
                _same(block.get("smiles"), block["catalog_record"]["smiles"]),
                "saved_blocks",
                "Audited outer and nested catalog leaf identities disagree",
            )
        else:
            optional_leaf_fields_absent.append(f"route_audit/building_blocks/{ident}/catalog_record/smiles")
        per_node[ident] = block.get("smiles")
    visiting, visited, products, queries = set(), set(), {}, []

    def visit(ident):
        records.require(
            ident in by_id and ident not in visiting, "saved_graph", "Missing dependency or cycle"
        )
        if ident in visited:
            return
        visiting.add(ident)
        node = by_id[ident]
        if "building_block_id" in node:
            smiles = block_lookup.get(node["building_block_id"], per_node.get(ident))
            records.require(
                records.smiles_identity(smiles)["key"] is not None,
                "saved_blocks",
                "Exact saved leaf identity unavailable",
            )
            if ident in per_node:
                records.require(_same(smiles, per_node[ident]), "saved_blocks", "Saved leaf layers disagree")
            products[ident] = smiles
        else:
            refs = node.get("reactants")
            records.require(
                isinstance(refs, list) and 0 < len(refs) <= 256 and all(isinstance(ref, str) for ref in refs),
                "saved_graph",
                "Explicit ordered reactant slots required",
            )
            for ref in refs:
                visit(ref)
            step = steps_by_id[ident]
            expected = [products[ref] for ref in refs]
            supplied = step.get("reactants", [])
            records.require(
                isinstance(supplied, list)
                and len(supplied) == len(expected)
                and all(_same(a, b) for a, b in zip(expected, supplied))
                and _same(node.get("product"), step.get("product"))
                and step.get("rule_id") == node["rule_id"],
                "saved_step_binding",
                "Graph and saved reaction slots/product disagree",
            )
            products[ident] = node["product"]
            queries.append(
                {
                    "node_id": ident,
                    "product": node["product"],
                    "reactants": expected,
                    "context": contexts.get(ident),
                }
            )
        visiting.remove(ident)
        visited.add(ident)

    visit(root)
    records.require(visited == set(by_id), "saved_graph", "Unused graph nodes")
    return queries, {
        "route_id": route_id,
        "root_node_id": root,
        "final_product": root_product,
        "routes_json_sha256": records.sha(raw_routes),
        "results_json_sha256": records.sha(raw_results),
        "saved_route_sha256": records.binding(saved),
        "saved_result_sha256": records.binding(result),
        "query_contexts_sha256": records.binding(contexts),
        "optional_leaf_fields_absent": optional_leaf_fields_absent,
        "scope": "Saved-value graph/final-product consistency only; no rule replay, catalog lookup, "
        "score validation, eligibility validation or experimental authentication",
    }


def match_saved_route(bundle_path, run_dir, route_id, *, contexts=None, expected_snapshot_sha256=None):
    """Bind exact source-record comparisons to every step of one saved final route."""
    queries, saved = _saved_queries(run_dir, route_id, {} if contexts is None else contexts)
    loaded = records.load_bundle(bundle_path, expected_snapshot_sha256=expected_snapshot_sha256)
    steps = [
        {"node_id": q["node_id"], **records._match(loaded, q["product"], q["reactants"], q["context"])}
        for q in queries
    ]
    return {
        "schema": records.VERSION,
        "interpretation": records.INTERPRETATION,
        "integrity": {k: loaded[k] for k in ("bundle_sha256", "receipt", "runtime")},
        "source_as_declared": loaded["source"],
        "saved_run_binding": saved,
        "steps": steps,
        "route_scope": {
            "reported_product_matches": sum(bool(s["matches"]) for s in steps),
            "step_count": len(steps),
            "interpretation": "Counts of exact reported compound entries only; "
            "no experimental route validation or acceptance decision",
        },
    }


def render_captured_report(document):
    """Escape all source text; this does not freshly verify raw records or saved runs."""
    records.require(
        isinstance(document, dict) and document.get("schema") == records.VERSION,
        "report_schema",
        "Expected captured public reaction-record report",
    )
    esc = lambda value: html.escape(str(value), quote=True)
    rows = document.get("steps", [document["result"]] if "result" in document else [])
    parts = [
        "<!doctype html><html lang='en'><meta charset='utf-8'><title>Public reaction records</title>",
        "<style>body{font:16px system-ui;max-width:1100px;margin:2em auto;padding:1em}"
        "pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f3f5f7;padding:1em}"
        "td,th{padding:.5em;text-align:left;border-bottom:1px solid #ddd}</style>",
        "<h1>Captured public reaction-record report</h1>",
        "<p>Rendering displays captured assertions. It does not freshly verify raw bytes, chemistry or "
        "saved-run artifacts. Generate a new match to replay those bindings.</p>",
        "<p>" + esc(records.INTERPRETATION) + "</p>",
    ]
    if "saved_run_binding" in document:
        saved = document["saved_run_binding"]
        parts += [
            "<h2>Saved final product and route</h2><p>"
            + esc(saved.get("route_id"))
            + ": "
            + esc(saved.get("final_product"))
            + "</p>"
        ]
    for row in rows:
        parts += [
            "<h2>" + esc(row.get("node_id", "Product lookup")) + "</h2>",
            "<p>" + esc(row["query"]["product"]) + " — " + esc(row["status"]) + "</p>",
        ]
        for match in row["matches"]:
            parts += [
                "<h3>" + esc(match["reaction_id"]) + "</h3>",
                "<p>"
                + esc(match["claim_scope"])
                + "; reactants: "
                + esc(match["reactants"])
                + "; context: "
                + esc(match["context"])
                + "</p>",
                "<p>" + esc(match["scope_limit"]) + "</p>",
                "<table><tr><th>Reported measurement type</th><th>Source JSON value</th><th>Analysis</th></tr>",
            ]
            for measurement in match["measurements"]:
                value = {
                    k: measurement[k]
                    for k in ("percentage", "float_value", "string_value", "amount")
                    if k in measurement
                }
                parts += [
                    "<tr><td>"
                    + esc(measurement.get("type", "UNSPECIFIED"))
                    + "</td><td>"
                    + esc(json.dumps(value, sort_keys=True))
                    + "</td><td>"
                    + esc(measurement.get("analysis_key", "not reported"))
                    + "</td></tr>"
                ]
            parts += [
                "</table><details><summary>Context, limitations, wire values and presence</summary><pre>"
                + esc(json.dumps(match, indent=2, sort_keys=True, allow_nan=False))
                + "</pre></details>"
            ]
    parts += [
        "<details><summary>Full captured bindings and source record details</summary><pre>"
        + esc(json.dumps(document, indent=2, sort_keys=True, allow_nan=False))
        + "</pre></details></html>"
    ]
    return "\n".join(parts)


def write_report(document, output_dir):
    target = Path(output_dir)
    records.require(not target.exists(), "output_exists", "Use a new report directory")
    rendered = render_captured_report(document)
    raw = records.encoded(document) + b"\n"
    target.mkdir(parents=True)
    (target / "reaction-records.json").write_bytes(raw)
    (target / "index.html").write_text(rendered, encoding="utf-8")
    return {"report": str(target), "json_sha256": records.sha(raw)}


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Optional offline public reaction records; never acceptance evidence"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    example = sub.add_parser("example", help="Copy the attributed audited adverse fixture snapshot")
    example.add_argument("output")
    imp = sub.add_parser("import", help="Import selected retained ORD records (optional dependencies)")
    imp.add_argument("snapshot")
    imp.add_argument("output")
    imp.add_argument("--expected-snapshot-sha256")
    match = sub.add_parser("match", help="Replay raw records and report exact product or saved-route matches")
    match.add_argument("bundle")
    match.add_argument("output")
    selection = match.add_mutually_exclusive_group(required=True)
    selection.add_argument("--product")
    selection.add_argument("--run", help="Saved run directory containing results.json and routes.json")
    match.add_argument("--route-id")
    match.add_argument("--reactant", action="append")
    match.add_argument("--context", help="JSON object; for --run, keyed by reaction node ID")
    match.add_argument("--expected-snapshot-sha256")
    render = sub.add_parser(
        "render", help="Display captured JSON only; no fresh integrity/chemistry verification"
    )
    render.add_argument("report_json")
    render.add_argument("output")
    args = parser.parse_args(argv)
    try:
        if args.command == "example":
            records.require(not Path(args.output).exists(), "output_exists", "Use a new directory")
            shutil.copytree(records.RESOURCE, args.output)
            result = {
                "snapshot": args.output,
                "snapshot_sha256": records.sha(records.read_local(args.output, "snapshot.json")),
            }
        elif args.command == "import":
            result = records.import_snapshot(
                args.snapshot, args.output, expected_snapshot_sha256=args.expected_snapshot_sha256
            )
        elif args.command == "render":
            path = Path(args.report_json)
            result = write_report(records.json_read(records.read_local(path.parent, path.name)), args.output)
        else:
            context = None
            if args.context:
                path = Path(args.context)
                context = records.json_read(records.read_local(path.parent, path.name))
            if args.run:
                records.require(
                    args.route_id and not args.reactant,
                    "query_arguments",
                    "--run requires --route-id and uses saved reactants",
                )
                document = match_saved_route(
                    args.bundle,
                    args.run,
                    args.route_id,
                    contexts=context,
                    expected_snapshot_sha256=args.expected_snapshot_sha256,
                )
            else:
                records.require(not args.route_id, "query_arguments", "--route-id requires --run")
                document = records.match_product(
                    args.bundle,
                    args.product,
                    reactants=args.reactant,
                    context=context,
                    expected_snapshot_sha256=args.expected_snapshot_sha256,
                )
            result = write_report(document, args.output)
        print(json.dumps(result, indent=2))
        return 0
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"reaction-records: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
