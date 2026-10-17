"""Offline, standard-library-only review of recorded workflow artifacts.

Consistency checks are joins and saved-value comparisons, not scientific validation.
No supplied path, SMILES or artifact text is executed or used to load external data.
"""

import base64
import argparse
import copy
import hashlib
import html
import json
import math
import os
from pathlib import Path
import statistics
import stat
from urllib.parse import quote

from . import review_display as display
from . import evidence_outcomes

JSON_FILES = (
    "manifest.json",
    "results.json",
    "qualification.json",
    "summary.json",
    "routes.json",
    "refinement.json",
    "scores.json",
)
LINK_FILES = (
    "results.tsv",
    "poses_hit.sdf",
    "poses_design.sdf",
    "references_hit.sdf",
    "references_design.sdf",
)
MAX_FILE_BYTES = 64 * 1024 * 1024
MAX_TOTAL_BYTES = 256 * 1024 * 1024
UNKNOWN = "not recorded"


def _binding(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _finite_tree(value):
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("Nonfinite numeric value")
    if isinstance(value, dict):
        for child in value.values():
            _finite_tree(child)
    elif isinstance(value, list):
        for child in value:
            _finite_tree(child)


def _read(path, root):
    # Reject symlinks, including in the last component, instead of following assets.
    if path.is_symlink() or path.resolve().parent != root or not path.is_file():
        raise ValueError("Not a regular contained artifact")
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(fd, "rb") as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode) or before.st_size > MAX_FILE_BYTES:
            raise ValueError("Not a regular artifact within 64 MiB limit")
        data = stream.read(MAX_FILE_BYTES + 1)
        after = os.fstat(stream.fileno())
    if len(data) > MAX_FILE_BYTES:
        raise ValueError("Artifact exceeds 64 MiB limit")
    stamp = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    if stamp(before) != stamp(after) or stamp(after) != stamp(path.stat()) or path.is_symlink():
        raise ValueError("Artifact changed while reading")
    return data, stamp(after)


def _load(root):
    docs, inventory, issues, observed = {}, {}, [], {}
    total = 0
    for name in JSON_FILES + LINK_FILES + (display.ARTIFACT,):
        path = root / name
        observed[name] = None
        if not path.exists() and not path.is_symlink():
            if name in JSON_FILES:
                issues.append(f"{name}: missing artifact")
            continue
        try:
            data, stamp = _read(path, root)
            total += len(data)
            sha = hashlib.sha256(data).hexdigest()
            observed[name] = (sha, stamp)
            if total > MAX_TOTAL_BYTES:
                raise ValueError("Bundle exceeds 256 MiB limit")
            inventory[name] = {"sha256": sha, "bytes": len(data)}
            if name in JSON_FILES or name == display.ARTIFACT:
                value = json.loads(data, object_pairs_hook=_object)
                _finite_tree(value)
                if (
                    not isinstance(value, dict)
                    or type(value.get("schema_version")) is not int
                    or value["schema_version"] != 1
                ):
                    raise ValueError("Unsupported or missing schema_version (requires 1)")
                docs[name] = value
        except (OSError, ValueError, RecursionError, UnicodeError) as exc:
            issues.append(f"{name}: {exc}")
            if total > MAX_TOTAL_BYTES:
                break
    return docs, inventory, issues, observed


def _unchanged(root, observed):
    for name, expected in observed.items():
        path = root / name
        if expected is None:
            if path.exists() or path.is_symlink():
                raise ValueError(f"Input presence changed: {name}")
        else:
            data, stamp = _read(path, root)
            if (hashlib.sha256(data).hexdigest(), stamp) != expected:
                raise ValueError(f"Input changed during rendering: {name}")


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _number(value):
    try:
        return type(value) in (int, float) and math.isfinite(value)
    except OverflowError:
        return False


def _records(doc, key):
    value = doc.get(key)
    _require(isinstance(value, list) and all(isinstance(x, dict) for x in value), f"Invalid {key} list")
    return value


def _one(records, predicate, label):
    found = [r for r in records if predicate(r)]
    _require(len(found) == 1, f"Expected one {label}; found {len(found)}")
    return found[0]


def _equal(actual, expected, label):
    _require(actual == expected, f"Mismatch: {label}")


def _integer(value, label, minimum=0):
    _require(type(value) is int and value >= minimum, f"Invalid integer {label}")


def _unit(value, label, *, positive=False):
    _require(
        _number(value) and (0 < value if positive else 0 <= value) and value <= 1,
        f"Invalid [0,1] scalar {label}",
    )


def _settings(value):
    keys = {
        "num_conformers",
        "seed",
        "energy_window",
        "strain_penalty",
        "max_iterations",
        "pharmacophore_sigma",
        "reject_unspecified_stereo",
        "scoring_objective",
    }
    _require(isinstance(value, dict) and set(value) == keys, "Incomplete or unsupported conformer settings")
    for key in ("num_conformers", "max_iterations"):
        _integer(value[key], key, 1)
    _integer(value["seed"], "seed")
    _require(value["seed"] <= 2**31 - 1, "Seed outside supported range")
    for key in ("energy_window", "strain_penalty", "pharmacophore_sigma"):
        _require(_number(value[key]) and value[key] >= 0, f"Invalid setting {key}")
    _require(value["pharmacophore_sigma"] > 0, "Invalid pharmacophore sigma")
    _require(type(value["reject_unspecified_stereo"]) is bool, "Invalid stereo setting")
    _require(value["scoring_objective"] in ("joint-v1", "original-hit-v1"), "Unsupported scoring objective")


def _policy_scalars(policy):
    minimum = policy["minimum_scores"]
    _require(isinstance(minimum, dict), "Missing minimum scores")
    allowed = {
        f"{role}_{metric}" for role in ("hit", "design") for metric in ("shape", "pharmacophore", "tanimoto")
    }
    _require(set(minimum) <= allowed, "Unsupported score floor")
    for key, value in minimum.items():
        _unit(
            value,
            f"floor {key}",
            positive=policy["mode"] == "qualified"
            and key in {"hit_tanimoto", "hit_shape", "hit_pharmacophore"},
        )
    if policy["mode"] == "qualified":
        _require(
            {"hit_tanimoto", "hit_shape", "hit_pharmacophore"} <= set(minimum), "Missing qualified hit floors"
        )
        _require(
            policy["conformers"]["reject_unspecified_stereo"] is True,
            "Qualified policy contradicts strict stereo",
        )
        _equal(policy["refinement"]["threshold_policy"], "all_seeds", "qualified threshold policy")
    _settings(policy["conformers"])
    if "similarity_constraints" in policy:
        bound = policy["similarity_constraints"]["policy"]
        _unit(bound["minimum_polar_reference_recall"], "polar recall floor", positive=True)


def _final_settings(policy, refinement, row):
    protocol = policy["refinement"]
    base, final = protocol["num_conformers"], refinement["num_conformers"]
    _integer(base, "initial refinement size", 1)
    _integer(final, "final refinement size", 1)
    _integer(row["refinement_conformers"], "selected final size", 1)
    _equal(row["refinement_conformers"], final, "selected/report final size")
    maximum = protocol.get("adaptive_max_conformers")
    maximum = base if maximum is None else maximum
    growth = protocol["conformer_growth_factor"]
    _integer(maximum, "adaptive maximum", 1)
    _integer(growth, "conformer growth", 2)
    _require(base <= final <= maximum, "Final size outside recorded refinement bounds")
    sizes = [base]
    while sizes[-1] < final:
        _require(len(sizes) < 128, "Unsupported refinement stage count")
        sizes.append(min(maximum, sizes[-1] * growth))
    _equal(sizes[-1], final, "final size on recorded adaptive schedule")
    _equal(refinement["threshold_policy"], protocol["threshold_policy"], "refinement threshold policy")
    _require(protocol["threshold_policy"] in ("all_seeds", "representative"), "Unsupported threshold policy")
    stages = refinement.get("stages")
    if stages is not None:
        _require(
            isinstance(stages, list) and bool(stages) and all(isinstance(s, dict) for s in stages),
            "Malformed refinement stages",
        )
        for stage in stages:
            _integer(stage.get("num_conformers"), "saved stage size", 1)
        _equal([s.get("num_conformers") for s in stages], sizes, "saved adaptive stages")
        if "adaptive_stages" in refinement:
            _integer(refinement["adaptive_stages"], "adaptive stage count", 1)
            _equal(refinement["adaptive_stages"], len(stages), "adaptive stage count")
        for key in (
            "trials",
            "status",
            "pose_seed",
            "num_conformers",
            "threshold_policy",
            "seed_count",
            "seed_successes",
        ):
            _equal(stages[-1].get(key), refinement.get(key), f"final stage {key}")
    else:
        _require(final == base, "Adaptive final evidence lacks recorded stages")
    return {**policy["conformers"], "num_conformers": final}


def _route_steps(route, audit, support, row):
    """Compare saved graph records only; no reaction or molecular interpretation."""
    nodes = _records(route["route"], "nodes")
    ids = [n.get("id") for n in nodes]
    _require(
        all(isinstance(i, str) and i for i in ids) and len(set(ids)) == len(ids),
        "Invalid route node identifiers",
    )
    by_id = dict(zip(ids, nodes))
    reactions = {n["id"]: n for n in nodes if "rule_id" in n}
    _require(bool(reactions) and route["route"]["root"] in reactions, "No recorded reaction root")
    _integer(row["reaction_steps"], "reaction step count", 1)
    _equal(row["reaction_steps"], len(reactions), "selected reaction count")
    blocks = _records(audit, "building_blocks")
    block_ids = [b.get("node_id") for b in blocks]
    expected_blocks = {n["id"] for n in nodes if "building_block_id" in n}
    _require(
        len(block_ids) == len(set(block_ids)) and set(block_ids) == expected_blocks,
        "Audit building-block coverage differs from route",
    )
    _require(
        set(by_id) == set(reactions) | expected_blocks and not (set(reactions) & expected_blocks),
        "Unsupported route node kind",
    )
    reachable, pending = set(), [route["route"]["root"]]
    while pending:
        ident = pending.pop()
        _require(isinstance(ident, str) and ident in by_id, "Unknown route dependency")
        if ident not in reachable:
            reachable.add(ident)
            if ident in reactions:
                refs = reactions[ident]["reactants"]
                _require(isinstance(refs, list), "Invalid reaction dependencies")
                pending.extend(refs)
    _equal(reachable, set(by_id), "root-reachable route graph")
    if "strict_eligible" in audit:
        _require(audit["strict_eligible"] is True, "Audit eligibility contradicts route support")
    for block in blocks:
        identifier = by_id[block["node_id"]]["building_block_id"]
        _require(isinstance(identifier, str) and bool(identifier), "Invalid route building-block identifier")
        _equal(block["catalog_record"]["id"], identifier, "route leaf/catalog identifier")
    products = {b["node_id"]: b["smiles"] for b in blocks}
    products.update({ident: node["product"] for ident, node in reactions.items()})
    layers = [_records(route["verification"], "steps"), _records(audit, "steps"), _records(support, "steps")]
    for layer in layers:
        step_ids = [s.get("node_id") for s in layer]
        _require(
            len(step_ids) == len(set(step_ids)) and set(step_ids) == set(reactions),
            "Recorded reaction step coverage differs",
        )
    _equal(
        [s["node_id"] for s in layers[0]], [s["node_id"] for s in layers[1]], "verification/audit step order"
    )
    _equal([s["node_id"] for s in layers[1]], [s["node_id"] for s in layers[2]], "audit/support step order")
    seen = set(expected_blocks)
    for verified, audited, step in zip(*layers):
        ident = verified["node_id"]
        node = reactions[ident]
        refs = node["reactants"]
        _require(
            isinstance(refs, list)
            and bool(refs)
            and all(isinstance(ref, str) and ref in seen for ref in refs),
            "Reaction dependency/order mismatch",
        )
        expected = [products[ref] for ref in refs]
        for record in (verified, audited, step):
            _equal(record["product"], node["product"], "reaction node/product correspondence")
        _equal(verified["rule_id"], node["rule_id"], "verification rule")
        _equal(audited["rule_id"], node["rule_id"], "audit rule")
        _equal(verified["reactants"], expected, "verified reactants")
        _equal(step["reactants"], expected, "supported reactants")
        slots = _records(audited, "reactant_slots")
        _equal([s.get("node_id") for s in slots], refs, "audit reactant node slots")
        for index, slot in enumerate(slots):
            _integer(slot["slot"], "reactant slot")
            _equal(slot["slot"], index, "reactant slot order")
        _equal([s["smiles"] for s in slots], expected, "audit reactant identities")
        if "strict_eligible" in audited:
            _require(audited["strict_eligible"] is True, "Audited step eligibility contradicts support")
        _equal(step["rule_replay"], "matched", "recorded rule replay")
        _require(step["eligibility"]["passed"] is True, "Step eligibility contradicts route support")
        if "unmet_requirements" in step["eligibility"]:
            _equal(step["eligibility"]["unmet_requirements"], [], "eligible step unmet requirements")
        seen.add(ident)


def _summary_counts(summary, rows, policy):
    counters = (
        "proposals",
        "verified_routes",
        "scored_products",
        "accepted_products",
        "returned_products",
        "qualified_products",
        "screened_products",
        "screened_product_design_pairs",
        "screened_scaffold_families",
        "refinement_shortlisted_families",
        "refinement_budgeted_families",
        "refinement_attempted",
        "refinement_accepted",
        "eligible_products",
        "carried_eligible_products",
        "new_eligible_products",
        "carried_screened_products",
        "accepted_scaffold_families",
        "returned_scaffold_families",
        "synthesis_evidence_evaluated_routes",
    )
    for key in counters:
        if key in summary:
            _integer(summary[key], f"summary {key}")
    if "route_rejections" in summary:
        _require(isinstance(summary["route_rejections"], dict), "Invalid route rejection counts")
        for count in summary["route_rejections"].values():
            _integer(count, "route rejection count")
    returned = len(rows)
    _equal(summary["returned_products"], returned, "returned product count")
    # These final counters describe accepted rows before output truncation, not
    # all discovered/scored molecules or numbers of route records.
    accepted = summary["accepted_products"]
    if summary.get("status") == "no_candidates":
        _equal(accepted, 0, "terminal no-candidate accepted count")
    _require(accepted >= returned, "Accepted count below returned count")
    qualified = summary.get("qualified_products")
    if rows or summary.get("status") in ("success", "no_candidates"):
        _require(qualified is not None, "Missing final qualified count")
    if qualified is not None:
        _require(
            qualified >= sum(row.get("recommendation_status") == "qualified" for row in rows),
            "Qualified count below returned qualified count",
        )
        if summary.get("status") in ("success", "no_candidates"):
            _equal(
                qualified,
                accepted if policy["mode"] == "qualified" else 0,
                "qualified/accepted count relationship",
            )
    for lower, upper in (
        ("returned_scaffold_families", "accepted_scaffold_families"),
        ("accepted_scaffold_families", "accepted_products"),
        ("returned_scaffold_families", "returned_products"),
        ("refinement_budgeted_families", "refinement_shortlisted_families"),
    ):
        if lower in summary and upper in summary:
            _require(summary[lower] <= summary[upper], f"Counter relationship {lower}/{upper}")


def _pose(saved, smiles, label):
    _require(isinstance(saved, dict), f"Missing saved {label}")
    _equal(saved.get("smiles"), smiles, f"{label} exact identity")
    _require(
        type(saved.get("conformer_id")) is int and saved["conformer_id"] >= 0, f"Invalid {label} conformer ID"
    )
    xyz = saved.get("coordinates")
    _require(
        isinstance(xyz, list)
        and bool(xyz)
        and all(isinstance(p, list) and len(p) == 3 and all(_number(v) for v in p) for p in xyz),
        f"Invalid {label} saved coordinates",
    )
    _require(
        isinstance(saved.get("atom_order_sha256"), str) and len(saved["atom_order_sha256"]) == 64,
        f"Missing {label} atom-order binding",
    )
    _equal(
        saved.get("sha256"),
        _binding({k: v for k, v in saved.items() if k != "sha256"}),
        f"{label} saved binding",
    )


def _witness(trial, row, hit, policy):
    saved = trial["similarity_constraints"]
    witness, bound = saved["witness"], saved["policy"]
    _equal(bound["hit"], hit, "witness original hit")
    _equal(
        bound["minimum_scores"],
        {key: policy["minimum_scores"][key] for key in ("hit_tanimoto", "hit_shape", "hit_pharmacophore")},
        "witness hit-only floors",
    )
    selection = policy["similarity_constraints"]
    _equal(bound, selection["policy"], "complete bound similarity policy")
    _equal(_binding(bound), selection["sha256"], "effective similarity policy hash")
    _equal(witness["policy_sha256"], _binding(bound), "witness policy hash")
    _require(
        witness["passed"] is True and witness["joint_assignment"]["feasible"] is True,
        "Witness is not recorded as passing",
    )
    polar = witness["polar"]
    _integer(polar["reference_count"], "polar reference count")
    if polar["status"] == "not_applicable":
        _require(polar["reference_count"] == 0 and polar["recall"] is None, "Contradictory no-polar witness")
    else:
        _require(polar["reference_count"] > 0, "Applicable polar witness has no reference features")
        _unit(polar["recall"], "saved polar recall")
        _require(
            polar["status"] == "passed" and polar["recall"] >= bound["minimum_polar_reference_recall"],
            "Recorded polar recall contradicts floor",
        )
    evidence = trial["selective_pose_evidence"]
    for key, smiles in (
        ("hit_candidate", row["smiles"]),
        ("design_candidate", row["smiles"]),
        ("design_reference", row["design"]),
    ):
        _pose(evidence[key], smiles, key)
    _pose(witness["reference_pose"], hit, "hit reference")
    _equal(witness["candidate_pose_sha256"], evidence["hit_candidate"]["sha256"], "witness candidate pose")
    _equal(
        evidence["hit_candidate"]["conformer_id"], trial["candidate_conformer_id"], "candidate conformer ID"
    )
    _equal(
        evidence["design_candidate"]["conformer_id"],
        trial["candidate_conformer_id"],
        "design candidate conformer ID",
    )
    _equal(
        witness["reference_pose"]["conformer_id"],
        trial["hit"]["reference_conformer_id"],
        "hit reference conformer ID",
    )
    _equal(
        evidence["design_reference"]["conformer_id"],
        trial["design"]["reference_conformer_id"],
        "design reference conformer ID",
    )
    return witness


def _seed_row(trial):
    witness = _mapping(_mapping(trial.get("similarity_constraints")).get("witness"))
    polar = _mapping(witness.get("polar"))
    return {
        "seed": trial.get("seed", UNKNOWN),
        "status": trial.get("status", UNKNOWN),
        "reason": trial.get("reason", UNKNOWN),
        "hit_shape": _mapping(trial.get("hit")).get("shape", UNKNOWN),
        "hit_pharmacophore": _mapping(trial.get("hit")).get("pharmacophore", UNKNOWN),
        "score": trial.get("score", UNKNOWN),
        "polar_recall": polar.get("recall", UNKNOWN),
        "polar_status": polar.get("status", UNKNOWN),
        "witness_sha256": _binding(witness) if witness else UNKNOWN,
        "candidate_pose_sha256": witness.get("candidate_pose_sha256", UNKNOWN),
        "reference_pose_sha256": _mapping(witness.get("reference_pose")).get("sha256", UNKNOWN),
    }


def _forward_records(prediction, product):
    """Join recorded strict hypotheses; never decode tokens or normalize molecules."""
    raw = prediction.get("raw", {})
    _require(isinstance(raw, dict), "Malformed raw forward evidence")
    source = prediction if "predictions" in prediction else raw
    if "predictions" not in source:
        return  # Legacy detail absence is displayed explicitly, not reconstructed.
    hypotheses = _records(source, "predictions")
    _require(bool(hypotheses), "Matched forward support has no recorded hypotheses")
    ranks, matches = [], []
    for item in hypotheses:
        rank = item.get("rank")
        _integer(rank, "forward hypothesis rank", 1)
        ranks.append(rank)
        _require(type(item.get("complete")) is bool, "Invalid forward hypothesis completion")
        value = item.get("product")
        _require(
            "product" in item and (value is None or isinstance(value, str) and bool(value)),
            "Invalid recorded forward product",
        )
        _require(item["complete"] or value is None, "Incomplete forward hypothesis has a product")
        if item["complete"] and value == product:
            matches.append(rank)
    _require(len(set(ranks)) == len(ranks), "Duplicate forward hypothesis ranks")
    _equal(ranks, list(range(1, len(hypotheses) + 1)), "forward hypothesis rank order")
    _require(bool(matches), "No completed recorded hypothesis matches the exact step product")
    _equal(prediction["matched_rank"], min(matches), "forward matched hypothesis rank")
    if "predictions" in prediction and "predictions" in raw:
        _equal(_binding(raw["predictions"]), _binding(hypotheses), "raw/reported forward hypotheses")


def _candidate(row, docs, policy):
    result = {"row": row, "issues": [], "seed_rows": [], "derived_minima": {}, "route": {}, "refinement": {}}
    try:
        hit = docs["manifest.json"]["hit"]
        for key in ("smiles", "design", "route_id"):
            _require(isinstance(row.get(key), str) and bool(row[key]), f"Missing selected {key}")
        routes = _records(docs["routes.json"], "routes")
        route = _one(routes, lambda r: r.get("route_id") == row["route_id"], "primary route")
        result["route"] = route
        _equal(route["design"], row["design"], "primary route design")
        _equal(route["verification"]["product"], row["smiles"], "primary regenerated product")
        root_node = _one(
            _records(route["route"], "nodes"),
            lambda n: n.get("id") == route["route"]["root"],
            "route root node",
        )
        _equal(root_node["product"], row["smiles"], "route root exact product")

        decision = _one(
            _records(docs["qualification.json"], "decisions"),
            lambda r: r.get("route_id") == row["route_id"],
            "primary qualification decision",
        )
        _equal(route.get("qualification"), decision, "route qualification decision")
        if row.get("recommendation_status") == "qualified":
            _require(
                policy["mode"] == "qualified"
                and decision.get("qualified") is True
                and decision.get("status") == "qualified",
                "Recorded qualification disagrees",
            )
            _equal(row["qualification_policy"], policy["policy_version"], "qualification policy")
        if policy["mode"] == "qualified":
            _equal(row.get("recommendation_status"), "qualified", "selected qualified mode")
            for key in policy["minimum_scores"]:
                _require(_number(row.get(key)), f"Missing finite selected {key}")
                _equal(decision["scores"][key], row[key], f"decision score {key}")
        for key in (
            "hit_tanimoto",
            "hit_morgan",
            "design_tanimoto",
            "hit_shape",
            "hit_pharmacophore",
            "design_shape",
            "design_pharmacophore",
        ):
            if key in row:
                _unit(row[key], f"selected {key}")
        for key, cutoff in policy["minimum_scores"].items():
            _require(_number(row.get(key)) and row[key] >= cutoff, f"Selected {key} contradicts floor")
        profile = policy["calibration"]
        if row.get("recommendation_status") == "qualified":
            _equal(decision["policy_version"], policy["policy_version"], "decision policy version")
            _equal(
                decision["calibration"],
                {key: profile[key] for key in ("name", "status", "artifact_sha256")},
                "decision profile",
            )
            _equal(decision["reasons"], [], "qualified decision reasons")
            if "similarity_constraints" in policy:
                _equal(
                    decision["similarity_constraints"],
                    policy["similarity_constraints"],
                    "decision selective policy",
                )
        for key, expected in (
            ("recommendation_profile", profile.get("name") or "none"),
            ("calibration_status", profile["status"]),
            ("calibration_artifact_sha256", profile.get("artifact_sha256") or "none"),
            ("evidence_profile", policy["evidence_profile"]),
        ):
            _equal(row.get(key), expected, key)
        support = route.get("synthesis_evidence", {})
        if policy["evidence_profile"] == "computational-route-support-v1":
            _require(
                type(support["schema_version"]) is int and support["schema_version"] == 2,
                "Unsupported support schema",
            )
            _equal(support["profile"], policy["evidence_profile"], "support profile")
            _equal(support["evidence_level"], "conditional_computational_selection_support", "support scope")
            _equal(support["product"], row["smiles"], "support product")
            _equal(support["route_binding_sha256"], row["route_support_binding"], "selected support binding")
            audit = support["route_audit"]
            _equal(audit["route_binding_sha256"], support["route_binding_sha256"], "audit binding")
            _equal(audit["product"], row["smiles"], "audit product")
            _equal(audit["root"], route["route"]["root"], "audit root")
            _require(support["eligibility"]["passed"] is True, "Support is not recorded as eligible")
            for leaf in _records(audit, "building_blocks"):
                record = leaf.get("catalog_record")
                _require(isinstance(record, dict), "Malformed catalog leaf record")
                if "identifier_provenance" in record:
                    origin = record["identifier_provenance"]
                    _require(isinstance(origin, dict), "Malformed identifier provenance")
                    _records(origin, "source_records")
                for flag in ("supplier_sku_verified", "current_stock_verified"):
                    if flag in record:
                        _require(type(record[flag]) is bool, f"Malformed {flag}")
            _route_steps(route, audit, support, row)
            steps = support["steps"]
            _require(
                isinstance(steps, list) and len(steps) == len(audit["steps"]), "Support step count differs"
            )
            for step, audited in zip(steps, audit["steps"]):
                for key in ("node_id", "product", "step_binding_sha256"):
                    _equal(step[key], audited[key], f"support step {key}")
                prediction = step["forward_prediction"]
                _require(
                    prediction.get("status") == "matched" and prediction.get("matching_policy") == "exact",
                    "Forward support contradicts eligible route",
                )
                _require(
                    type(prediction.get("matched_rank")) is int and prediction["matched_rank"] > 0,
                    "Missing matched forward rank",
                )
                _forward_records(prediction, step["product"])
                _equal(
                    step["reactants"],
                    [slot["smiles"] for slot in audited["reactant_slots"]],
                    "support reactant slots",
                )
        refinement = _one(
            _records(docs["refinement.json"], "candidates"),
            lambda r: (
                r.get("product") == row["smiles"]
                and r.get("design") == row["design"]
                and r.get("status") == "success"
            ),
            "successful exact product/design refinement",
        )
        result["refinement"] = refinement
        final_settings = _final_settings(policy, refinement, row)
        trials = _records(refinement, "trials")
        result["seed_rows"] = [_seed_row(t) for t in trials]
        required = policy["refinement"]["seeds"]
        _require(
            isinstance(required, list)
            and bool(required)
            and all(type(s) is int for s in required)
            and len(set(required)) == len(required),
            "Invalid required seeds",
        )
        _equal([t.get("seed") for t in trials], required, "required seed identities/order")
        _require(all(t.get("status") == "success" for t in trials), "Incomplete required seeds")
        for trial in trials:
            _settings(trial["settings"])
            _equal(
                trial["settings"], {**final_settings, "seed": trial["seed"]}, "complete saved trial settings"
            )
            _equal(trial["scoring_objective"], final_settings["scoring_objective"], "trial scoring objective")
            for role in ("hit", "design"):
                for metric in ("shape", "pharmacophore"):
                    value = trial[role][metric]
                    _unit(value, f"seed {trial['seed']} {role}_{metric}")
                    if policy["refinement"]["threshold_policy"] == "all_seeds":
                        _require(
                            value >= policy["minimum_scores"].get(f"{role}_{metric}", 0),
                            f"Required seed {role}_{metric} contradicts floor",
                        )
            for key, identity in (
                ("candidate_smiles", row["smiles"]),
                ("hit_smiles", hit),
                ("design_smiles", row["design"]),
            ):
                _equal(trial[key], identity, f"trial {key}")
            _equal(trial["settings"]["seed"], trial["seed"], "trial settings seed")
            _equal(
                trial["settings"]["scoring_objective"],
                policy["conformers"]["scoring_objective"],
                "trial objective",
            )
            _equal(trial["settings"]["num_conformers"], row["refinement_conformers"], "final ensemble size")
            _require(
                all(
                    _number(trial[role][metric])
                    for role in ("hit", "design")
                    for metric in ("shape", "pharmacophore")
                )
                and _number(trial["score"])
                and _number(trial["candidate_relative_energy"]),
                "Invalid trial numeric metrics",
            )
            if "similarity_constraints" in policy:
                _witness(trial, row, hit, policy)
        representative = sorted(trials, key=lambda t: (t["score"], t["seed"]))[(len(trials) - 1) // 2]
        scores = [t["score"] for t in trials]
        expected = {
            "pose_seed": representative["seed"],
            "seed_count": len(trials),
            "seed_successes": len(trials),
            "seed_score_median": statistics.median(scores),
            "seed_score_min": min(scores),
            "seed_score_max": max(scores),
            "seed_score_range": max(scores) - min(scores),
            "seed_score_std": statistics.pstdev(scores),
        }
        for key, value in expected.items():
            _equal(row.get(key), value, f"representative/seed summary {key}")
            _equal(refinement.get(key), value, f"refinement summary {key}")
        for key, value in {
            "score": representative["score"],
            "relative_energy": representative["candidate_relative_energy"],
            "scoring_objective": representative["settings"]["scoring_objective"],
            "ranking_policy": policy["ranking_policy"],
            **{
                f"{role}_{metric}": representative[role][metric]
                for role in ("hit", "design")
                for metric in ("shape", "pharmacophore")
            },
        }.items():
            _equal(row.get(key), value, f"representative {key}")
        saved_score = _one(
            _records(docs["scores.json"], "scores"),
            lambda r: (
                r.get("scoring_stage") == "refined"
                and r.get("route_id") == row["route_id"]
                and r.get("product") == row["smiles"]
                and r.get("design_smiles") == row["design"]
            ),
            "selected refined score",
        )
        _equal(saved_score["hit_smiles"], hit, "score original hit")
        for key in (
            "score",
            "candidate_relative_energy",
            "settings",
            "hit",
            "design",
            "candidate_conformer_id",
        ):
            _equal(saved_score.get(key), representative.get(key), f"saved representative {key}")
        if "similarity_constraints" in policy:
            for key in ("similarity_constraints", "selective_pose_evidence", "pose_selection"):
                _require(key in representative, f"Missing representative {key}")
                _equal(saved_score.get(key), representative[key], f"saved representative {key}")
            witness = representative["similarity_constraints"]["witness"]
            for key, value in {
                "similarity_witness_sha256": _binding(witness),
                "similarity_constraint_sha256": witness["policy_sha256"],
                "similarity_constraint_policy": policy["similarity_constraints"]["policy"]["version"],
                "pose_selection_policy": representative["pose_selection_policy"],
                "polar_reference_recall": witness["polar"]["recall"],
            }.items():
                _equal(row.get(key), value, f"representative {key}")
        result["derived_minima"] = {
            "hit_shape": min(t["hit"]["shape"] for t in trials),
            "hit_pharmacophore": min(t["hit"]["pharmacophore"] for t in trials),
        }
        recalls = [s["polar_recall"] for s in result["seed_rows"] if _number(s["polar_recall"])]
        result["derived_minima"]["polar_recall"] = (
            min(recalls)
            if len(recalls) == len(trials)
            else "not applicable (all seeds have no polar reference features)"
            if all(s["polar_status"] == "not_applicable" for s in result["seed_rows"])
            else UNKNOWN
        )
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as exc:
        result["issues"].append(f"Candidate evidence incomplete or inconsistent: {exc}")
    return result


_IDENTITY_POLICY = "protonation-tautomer-family-v1"
_IDENTITY_RECORDED = (
    "original_hit",
    "final_product",
    "change_class",
    "atom_aware_scaffold_changed",
    "generic_ring_linker_topology_changed",
    "scaffold_presence_change",
)


def _nonempty_string(value, label):
    _require(isinstance(value, str) and bool(value), f"Invalid {label}")


def _identity_descriptor(manifest):
    """Validate only saved provenance copies, never the installed producer/runtime."""
    result = {"status": "not_exported", "warnings": [], "record": None}
    if "identity_comparison_policy" not in manifest:
        return result
    result["status"] = "unusable"
    result["record"] = copy.deepcopy(manifest["identity_comparison_policy"])
    try:
        descriptor = result["record"]
        _require(isinstance(descriptor, dict), "Identity descriptor is not an object")
        _require(
            type(descriptor["schema_version"]) is int and descriptor["schema_version"] == 1,
            "Unsupported identity descriptor schema",
        )
        _equal(descriptor["policy_version"], _IDENTITY_POLICY, "identity descriptor policy")
        _nonempty_string(descriptor["rdkit_version"], "identity descriptor runtime")
        for key in ("tautomer_limit", "transform_limit"):
            _require(type(descriptor[key]) is int and descriptor[key] > 0, f"Invalid identity {key}")
        for key in ("identity.py", "common.py"):
            value = descriptor["code_sha256"][key]
            _require(
                isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value),
                f"Invalid identity source hash {key}",
            )
            _equal(value, manifest["code_sha256"][key], f"saved identity source {key}")
        result["status"] = "recorded"
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        result["warnings"].append(f"Identity descriptor unavailable or unsupported: {exc}")
    if manifest.get("sources_changed_since_import"):
        result["warnings"].append(
            "Producer sources changed since import; saved copies do not authenticate source execution."
        )
    return result


def _identity_comparison(row, manifest, policy, descriptor):
    """Check implications of a supplied label; never reconstruct a missing identity/class."""
    view = {"status": "not_exported", "issues": [], "record": None}
    try:
        if "hit_identity_comparison" not in row:
            _require(
                descriptor["status"] == "not_exported",
                "Declared identity comparison missing from selected row",
            )
            return view
        view["status"] = "inconsistent_metadata"
        record = row["hit_identity_comparison"]
        _require(isinstance(record, dict), "Identity comparison is not an object")
        _require(
            type(record["schema_version"]) is int and record["schema_version"] == 1,
            "Unsupported identity comparison schema",
        )
        _equal(record["reference_role"], "original_hit", "identity reference role")
        for key, expected in (
            ("original_hit_smiles", manifest["hit"]),
            ("final_product_smiles", row["smiles"]),
        ):
            _nonempty_string(record[key], key)
            _equal(record[key], expected, f"identity {key}")
        if record["status"] == "unavailable":
            _require(
                not any(key in record for key in _IDENTITY_RECORDED),
                "Unavailable identity contains recorded-only fields",
            )
            for key in ("reason_code", "reason"):
                _nonempty_string(record[key], f"identity {key}")
            if "error_type" in record:
                _require(isinstance(record["error_type"], str), "Invalid identity error_type")
            view.update(status="unavailable", record=copy.deepcopy(record))
            return view
        _equal(record["status"], "recorded", "identity status")
        _require(
            not any(k in record for k in ("reason_code", "reason", "error_type")),
            "Recorded identity contains unavailable fields",
        )
        _equal(descriptor["status"], "recorded", "usable identity descriptor")
        saved_policy = descriptor["record"]
        hit, product = record["original_hit"], record["final_product"]
        for identity, exact in ((hit, manifest["hit"]), (product, row["smiles"])):
            _require(isinstance(identity, dict), "Identity is not an object")
            for key in ("exact_smiles", "family_smiles", "policy_version", "rdkit_version"):
                _nonempty_string(identity[key], f"identity {key}")
            _equal(identity["exact_smiles"], exact, "retained exact identity")
            for key in ("policy_version", "rdkit_version"):
                _equal(identity[key], saved_policy[key], f"retained identity {key}")
            _require(type(identity["has_scaffold"]) is bool, "Invalid has_scaffold boolean")
            for key in ("scaffold_family_smiles", "ring_skeleton_smiles"):
                _require(
                    isinstance(identity[key], str) and bool(identity[key]) == identity["has_scaffold"],
                    f"Identity scaffold presence contradicts {key}",
                )
        for key in ("family_smiles", "scaffold_family_smiles"):
            _equal(product[key], row[key], f"selected product identity {key}")
        exact_same = hit["exact_smiles"] == product["exact_smiles"]
        family_same = hit["family_smiles"] == product["family_smiles"]
        scaffold_same = hit["scaffold_family_smiles"] == product["scaffold_family_smiles"]
        generic_same = hit["ring_skeleton_smiles"] == product["ring_skeleton_smiles"]
        _require(not exact_same or family_same, "Equal exact identities have different families")
        _require(not family_same or scaffold_same, "Equal families have different scaffolds")
        _require(
            not scaffold_same or generic_same, "Equal atom-aware scaffolds have different generic graphs"
        )
        for key, expected in (
            ("atom_aware_scaffold_changed", not scaffold_same),
            ("generic_ring_linker_topology_changed", not generic_same),
        ):
            _require(type(record[key]) is bool, f"Invalid {key} boolean")
            _equal(record[key], expected, f"saved {key}")
        presence = {
            (True, True): "both_present",
            (False, False): "both_absent",
            (False, True): "gained",
            (True, False): "lost",
        }
        _equal(
            record["scaffold_presence_change"],
            presence[hit["has_scaffold"], product["has_scaffold"]],
            "saved scaffold presence",
        )
        supplied = record["change_class"]
        _require(
            isinstance(supplied, str) and supplied in ("unchanged", "state_only", "decoration", "scaffold"),
            "Missing or unsupported recorded change class",
        )
        implications = {
            "unchanged": exact_same,
            "state_only": not exact_same and family_same,
            "decoration": not exact_same and not family_same and scaffold_same,
            "scaffold": not exact_same and not family_same and not scaffold_same,
        }
        _require(implications[supplied], "Recorded change class contradicts retained identities")
        if (
            row.get("recommendation_status") == "qualified"
            and _mapping(_mapping(policy).get("constraints")).get("mode") == "hop"
        ):
            _equal(supplied, "scaffold", "qualified hop identity class")
        view.update(status="recorded", record=copy.deepcopy(record))
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        view["status"] = "inconsistent_metadata"
        view["issues"].append(f"Saved identity comparison incomplete or inconsistent: {exc}")
    return view


def _assemble(docs, inventory, issues):
    view = {
        "schema_version": 1,
        "artifacts": inventory,
        "issues": [i for i in issues if not i.startswith(display.ARTIFACT + ":")],
        "candidates": [],
        "manifest": docs.get("manifest.json", {}),
        "summary": docs.get("summary.json", {}),
        "policy": docs.get("results.json", {}).get("policy"),
    }
    try:
        rows = _records(docs["results.json"], "results")
        policy = view["policy"]
        _require(isinstance(policy, dict), "Effective policy missing")
        manifest, qualification, summary = docs["manifest.json"], docs["qualification.json"], view["summary"]
        _equal(manifest["recommendation_policy"], policy, "manifest/results policy")
        _equal(qualification["policy"], policy, "qualification/results policy")
        _equal(manifest["recommendation_profile"], policy["calibration"], "manifest profile provenance")
        _equal(manifest["config"]["minimum_scores"], policy["minimum_scores"], "config floors")
        for key in ("minimum_scores", "refinement"):
            _equal(manifest["effective_workflow"][key], policy[key], f"effective workflow {key}")
        if "similarity_constraints" in policy:
            _equal(
                manifest["effective_workflow"]["similarity_constraints"],
                policy["similarity_constraints"],
                "effective similarity policy",
            )
        for name, key in (("routes.json", "routes"), ("qualification.json", "decisions")):
            records = _records(docs[name], key)
            ids = [item.get("route_id") for item in records]
            _require(
                all(isinstance(i, str) and i for i in ids) and len(set(ids)) == len(ids),
                f"Duplicate/invalid {key} identifiers",
            )

        _equal(
            summary["recommendation_profile"],
            {k: policy["calibration"][k] for k in ("name", "status", "artifact_sha256")},
            "summary profile",
        )
        _equal(summary["recommendation_mode"], policy["mode"], "summary mode")
        _policy_scalars(policy)
        _settings(manifest["effective_conformers"])
        _equal(manifest["effective_conformers"], policy["conformers"], "manifest complete conformer settings")
        _summary_counts(summary, rows, policy)
        _require(isinstance(manifest["hit"], str) and bool(manifest["hit"]), "Original hit missing")
        route_records = _records(docs["routes.json"], "routes")
        evidence_outcomes.validate_records(route_records, _records(qualification, "decisions"))
        evidence_outcomes.validate_summary(summary, route_records)
        _equal(qualification["selected_route_ids"], [r["route_id"] for r in rows], "selected route IDs")
        for key in ("route_id", "smiles"):
            identifiers = [r[key] for r in rows]
            _require(
                all(isinstance(x, str) and x for x in identifiers)
                and len(set(identifiers)) == len(identifiers),
                f"Duplicate/invalid selected {key}",
            )
        _require(type(summary["returned_products"]) is int, "Invalid returned product count")
        _equal(summary["returned_products"], len(rows), "returned product count")
        if rows:
            _equal(summary["status"], "success", "selected run terminal status")
        else:
            _require(
                summary.get("status") in ("no_candidates", "error", "running"),
                "Empty run status contradicts selection",
            )
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        view["issues"].append(f"Bundle incomplete or inconsistent: {exc}")
    view["identity_descriptor"] = _identity_descriptor(_mapping(view["manifest"]))
    raw_rows = docs.get("results.json", {}).get("results", [])
    if isinstance(raw_rows, list):
        for row in raw_rows:
            if isinstance(row, dict):
                candidate = _candidate(row, docs, view["policy"] or {})
                candidate["identity_comparison"] = _identity_comparison(
                    row, _mapping(view["manifest"]), view["policy"], view["identity_descriptor"]
                )
                candidate["issues"].extend(candidate["identity_comparison"]["issues"])
                view["candidates"].append(candidate)
    # Optional visuals are validated independently, after all existing checks.
    descriptor = _mapping(view["manifest"]).get("review_visuals")
    visuals = [
        {"status": "not_exported", "reason": "Legacy bundle: display evidence not exported"}
        for _ in view["candidates"]
    ]
    view["visual_issues"] = [i for i in issues if i.startswith(display.ARTIFACT + ":")]
    if "review_visuals" in _mapping(view["manifest"]):
        try:
            visuals = display.validate(
                docs.get(display.ARTIFACT),
                descriptor,
                inventory.get(display.ARTIFACT),
                [c["row"] for c in view["candidates"]],
                docs.get("refinement.json", {}).get("candidates", []),
                _mapping(view["manifest"]).get("hit"),
            )
        except (KeyError, TypeError, ValueError, AttributeError, OverflowError) as exc:
            message = f"Saved display metadata inconsistent: {exc}"
            view["visual_issues"].append(message)
            visuals = [{"status": "inconsistent", "reason": message} for _ in view["candidates"]]
    elif display.ARTIFACT in inventory:
        view["visual_issues"].append("Unadvertised display artifact ignored")
    for candidate, visual in zip(view["candidates"], visuals):
        candidate["visuals"] = copy.deepcopy(visual)
    view["policy"] = _mapping(view["policy"])
    view["artifact_consistency"] = (
        "consistent_recorded_bundle"
        if not view["issues"] and not any(c["issues"] for c in view["candidates"])
        else "incomplete_or_inconsistent"
    )
    return view


def build_report(run_dir):
    """Read bounded saved bytes and return a review view; no scientific validation."""
    root = Path(run_dir).resolve(strict=True)
    docs, inventory, issues, observed = _load(root)
    view = _assemble(docs, inventory, issues)
    _unchanged(root, observed)
    return view


def _text(value):
    if value is None:
        return "unknown / not recorded"
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)
    return str(value)


def _escape(value):
    return html.escape(_text(value), quote=True)


def _table(rows, headings=("Field", "Recorded value")):
    return (
        "<table><thead><tr>"
        + "".join(f"<th>{_escape(h)}</th>" for h in headings)
        + "</tr></thead><tbody>"
        + "".join("<tr>" + "".join(f"<td>{_escape(v)}</td>" for v in row) + "</tr>" for row in rows)
        + "</tbody></table>"
    )


def _fields(value):
    return _table(sorted(value.items())) if isinstance(value, dict) else f"<p>{_escape(value)}</p>"


def _mapping(value):
    return value if isinstance(value, dict) else {}


def _dicts(value):
    return [v for v in value if isinstance(v, dict)] if isinstance(value, list) else []


def _details(title, content):
    return f"<details><summary>{_escape(title)}</summary>{content}</details>"


def _route_html(route):
    evidence = _mapping(route.get("synthesis_evidence"))
    audit = _mapping(evidence.get("route_audit"))
    parts = [
        "<h3>Primary route and computational evidence</h3>",
        _fields(
            {
                "route_id": route.get("route_id"),
                "selection_status": route.get("selection_status"),
                "route_binding_sha256": evidence.get("route_binding_sha256"),
                "evidence_level": evidence.get("evidence_level", UNKNOWN),
                "forward_outcome": evidence.get("forward_outcome", UNKNOWN),
                "eligibility": evidence.get("eligibility", UNKNOWN),
            }
        ),
    ]
    steps = evidence.get("steps", [])
    if isinstance(steps, list):
        for step in steps:
            if not isinstance(step, dict):
                continue
            parts += [
                f"<h4>Step {_escape(step.get('node_id'))}</h4>",
                _fields(
                    {
                        k: step.get(k, UNKNOWN)
                        for k in (
                            "reactants",
                            "product",
                            "step_binding_sha256",
                            "rule_replay",
                            "forward_outcome",
                            "atom_sources",
                            "eligibility",
                        )
                    }
                ),
            ]
            prediction = _mapping(step.get("forward_prediction"))
            raw = _mapping(prediction.get("raw"))
            hypotheses = prediction.get("predictions", raw.get("predictions", UNKNOWN))
            if "raw" in prediction and not isinstance(prediction["raw"], dict):
                parts.append(
                    "<p>Raw forward evidence is malformed; hypothesis detail availability is unknown.</p>"
                )
            elif "predictions" not in prediction and "predictions" not in raw:
                parts.append(
                    "<p>Forward hypothesis details not recorded; any matched status is a saved claim only.</p>"
                )
            elif "predictions" not in prediction:
                parts.append(
                    "<p>Forward hypothesis details retained in raw evidence; reported list absent.</p>"
                )
            parts += [
                "<h5>Forward support (recorded)</h5>",
                _fields(
                    {
                        k: prediction.get(k, UNKNOWN)
                        for k in ("status", "reason", "reason_code", "matched_rank", "matching_policy")
                    }
                    | {"predictions": hypotheses}
                ),
                "<h5>Conditions and observations</h5>",
                _fields(
                    {
                        "context": step.get("context", UNKNOWN),
                        "observations": step.get("observations", UNKNOWN),
                    }
                ),
            ]
    if not steps:
        parts.append("<p>Stepwise forward evidence not recorded.</p>")
    parts.append("<details><summary>Concrete replay slots, intermediates and rules</summary>")
    for step in audit.get("steps", []) if isinstance(audit.get("steps"), list) else []:
        if isinstance(step, dict):
            parts.append(
                _fields(
                    {
                        k: step.get(k, UNKNOWN)
                        for k in (
                            "node_id",
                            "reactant_slots",
                            "product",
                            "rule_id",
                            "smarts",
                            "rule_source",
                            "step_binding_sha256",
                            "atom_origin_status",
                            "stereo_status",
                            "enumeration",
                        )
                    }
                )
            )
    parts.append("</details><h4>Catalog leaves and every recorded origin</h4>")
    for leaf in audit.get("building_blocks", []) if isinstance(audit.get("building_blocks"), list) else []:
        if not isinstance(leaf, dict):
            continue
        record = _mapping(leaf.get("catalog_record"))
        parts += [
            _fields(
                {
                    "node_id": leaf.get("node_id"),
                    "smiles": leaf.get("smiles"),
                    **{
                        k: record.get(k, UNKNOWN)
                        for k in (
                            "id",
                            "source",
                            "snapshot",
                            "availability",
                            "in_stock",
                            "catalog_sha256",
                            "supplier_sku_verified",
                            "current_stock_verified",
                        )
                    },
                }
            )
        ]
        origins = _mapping(record.get("identifier_provenance"))
        parts.append(
            _fields(
                {
                    "identifier_provenance_status": origins.get("status", "unknown_legacy"),
                    "manifest_sha256": origins.get("manifest_sha256", UNKNOWN),
                }
            )
        )
        for origin in _dicts(origins.get("source_records")):
            parts.append(_fields(origin))
    parts += [
        "<details><summary>Recorded limitations and predictor provenance</summary>",
        _fields(
            {
                "limitations": audit.get("limitations", UNKNOWN),
                "enumeration_truncated": audit.get("enumeration_truncated", UNKNOWN),
                "predictor_provenance": evidence.get("predictor_provenance", UNKNOWN),
                "interpretation": evidence.get("interpretation", UNKNOWN),
            }
        ),
    ]
    parts.append("</details>")
    return "".join(parts)


def _visual_html(visual, index):
    if visual.get("status") != "available":
        return (
            "<p>Structure display: "
            + _escape(visual.get("status", "not_exported"))
            + " — "
            + _escape(visual.get("reason", "not recorded"))
            + "</p>"
        )
    parts = ["<div class='depictions'>"]
    for role, label in (("original_hit", "Original hit"), ("final_product", "Regenerated final product")):
        image = visual["depictions"][role]
        parts.append("<figure><figcaption>" + label + " — producer 2D depiction</figcaption>")
        if image["status"] == "available":
            parts.append(
                "<img alt='"
                + label
                + " exact structure' width='480' height='320' src='data:image/png;base64,"
                + _escape(image["png_base64"])
                + "'>"
            )
        else:
            parts.append("<p>Unavailable: " + _escape(image["reason"]) + "</p>")
        parts.append("</figure>")
    parts.append("</div>")
    if visual["seeds"]:
        payload = json.dumps(
            {k: visual[k] for k in ("seeds", "representative_seed")}, allow_nan=False, separators=(",", ":")
        )
        parts.append(
            "<section class='pose-viewer' data-payload='"
            + html.escape(payload, quote=True)
            + "'><h3>Saved original-hit pose overlap</h3><p>Blue: original hit; magenta: final product. One saved seed and shared camera. This view does not realign, rescore or validate stereochemistry; use the 2D depiction and exact SMILES for structure review.</p>"
        )
        parts.append(f"<label for='pose-seed-{index}'>Saved seed</label><select id='pose-seed-{index}'>")
        representative = visual["representative_seed"]
        if representative not in [s["seed"] for s in visual["seeds"]]:
            parts.append(
                "<option selected value='"
                + _escape(representative)
                + "'>Representative not retained</option>"
            )
        for seed in visual["seeds"]:
            selected = " selected" if seed["seed"] == representative else ""
            parts.append(
                f"<option{selected} value='{_escape(seed['seed'])}'>{_escape(seed['seed'])} — {_escape(seed['pair']['status'])}</option>"
            )
        parts.append(
            "</select><label><input type='checkbox' checked>Original hit</label><label><input type='checkbox' checked>Final product</label><button type='button' data-action='reset'>Reset camera</button><button type='button' data-zoom='1.2'>Zoom in</button><button type='button' data-zoom='.833333333333'>Zoom out</button><canvas tabindex='0' aria-label='Shared saved-pose viewer: arrow keys rotate, plus and minus zoom, zero resets'></canvas><p>Drag/touch to rotate; shift-drag to pan; wheel or buttons to zoom. Keyboard: arrows, +/−, 0. Bonds show connectivity, not stereochemical certification.</p><p class='pose-note' role='status'></p><noscript>Interactive camera unavailable without JavaScript. Saved seed evidence remains in the tables below.</noscript></section>"
        )
    else:
        parts.append("<p>Saved 3D overlap unavailable: no retained selective pose trials.</p>")
    return "".join(parts)


def render_html(view, *, artifact_links=None):
    """Render a build_report view. Caller links are restricted to approved basenames."""
    links = artifact_links or {}
    for name, href in links.items():
        _require(
            name in JSON_FILES + LINK_FILES + (display.ARTIFACT,) and href == quote(name),
            "Unsafe artifact link",
        )
    policy = _mapping(view.get("policy"))
    profile = _mapping(policy.get("calibration"))
    bound = _mapping(_mapping(policy.get("similarity_constraints")).get("policy"))
    minimum = _mapping(policy.get("minimum_scores"))
    summary = view["summary"]
    parts = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width'><title>Candidate review</title><style>body{font:16px system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#172332}table{border-collapse:collapse;width:100%;margin:1rem 0;table-layout:fixed}th,td{border:1px solid #cbd3dc;padding:.5rem;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#eef2f6}article{border-top:3px solid #637a90;margin-top:2rem;padding-top:1rem}.notice{background:#fff4d6;padding:1rem}code{overflow-wrap:anywhere}.depictions{display:flex;flex-wrap:wrap;gap:1rem}.depictions figure{margin:0;flex:1;min-width:220px}.depictions img{width:100%;max-width:480px;height:auto}.pose-viewer canvas{width:100%;height:360px;touch-action:none;display:block}.pose-viewer button,.pose-viewer select{padding:.5rem;margin:.3rem}.pose-note{overflow-wrap:anywhere;font-size:.85rem}details{border:1px solid #cbd3dc;padding:.7rem;margin:1rem 0}summary{cursor:pointer;font-weight:600}@media(max-width:600px){body{font-size:14px}th,td{padding:.3rem}.scroll{overflow-x:auto}.scroll table{min-width:520px}}</style></head><body><h1>Candidate review</h1>",
        "<p class='notice'>Workflow-recorded evidence, not scientific revalidation. Computational route support does not establish experimental synthesis, activity, reaction conditions, supplier SKU or current stock.</p>",
        _table(
            [
                ("Run status", summary.get("status", UNKNOWN)),
                ("Artifact consistency", view["artifact_consistency"]),
                ("Profile", profile.get("name", UNKNOWN)),
                ("Profile status", profile.get("status", UNKNOWN)),
                ("Original-hit Tanimoto floor", minimum.get("hit_tanimoto", UNKNOWN)),
                ("Original-hit shape floor", minimum.get("hit_shape", UNKNOWN)),
                ("Original-hit pharmacophore floor", minimum.get("hit_pharmacophore", UNKNOWN)),
                ("Polar recall floor", bound.get("minimum_polar_reference_recall", UNKNOWN)),
            ]
        ),
    ]
    for issue in view["issues"]:
        parts.append(f"<p class='notice'>{_escape(issue)}</p>")
    for issue in view.get("visual_issues", []):
        parts.append(f"<p class='notice'>{_escape(issue)}</p>")
    for warning in _mapping(view.get("identity_descriptor")).get("warnings", []):
        parts.append(f"<p class='notice'>{_escape(warning)}</p>")
    for index, candidate in enumerate(view["candidates"], 1):
        row = candidate["row"]
        comparison = _mapping(candidate.get("identity_comparison"))
        identity = _mapping(comparison.get("record"))
        original_identity = _mapping(identity.get("original_hit"))
        consistent = view["artifact_consistency"] == "consistent_recorded_bundle"
        parts += [
            f"<article><h2>Candidate {index}</h2>",
            _visual_html(candidate.get("visuals", {}), index),
            _table(
                [
                    (
                        "Workflow-recorded status" if consistent else "Unconfirmed recorded claim",
                        row.get("recommendation_status", UNKNOWN),
                    ),
                    ("Original hit — exact SMILES", view["manifest"].get("hit", UNKNOWN)),
                    ("Regenerated product — exact SMILES", row.get("smiles", UNKNOWN)),
                    ("Secondary design — exact SMILES", row.get("design", UNKNOWN)),
                    ("Product normalized scaffold", row.get("scaffold_family_smiles", UNKNOWN)),
                    (
                        "Original-hit normalized scaffold",
                        original_identity.get("scaffold_family_smiles", UNKNOWN),
                    ),
                    ("Identity comparison availability", comparison.get("status", "not_exported")),
                    ("Unavailable comparison reason", identity.get("reason", UNKNOWN)),
                    (
                        "Scaffold change type",
                        identity.get("change_class", UNKNOWN),
                    ),
                    (
                        "Atom-aware scaffold changed (recorded)",
                        identity.get("atom_aware_scaffold_changed", UNKNOWN),
                    ),
                    (
                        "Generic topology comparison — ring/linker graph changed (recorded)",
                        identity.get("generic_ring_linker_topology_changed", UNKNOWN),
                    ),
                    ("Scaffold presence (recorded)", identity.get("scaffold_presence_change", UNKNOWN)),
                ]
            ),
        ]
        for issue in candidate["issues"]:
            parts.append(f"<p class='notice'>{_escape(issue)}</p>")
        parts += [
            "<h3>Representative exported pose</h3>",
            _table(
                [
                    ("Seed", row.get("pose_seed", UNKNOWN)),
                    ("Original-hit Tanimoto", row.get("hit_tanimoto", UNKNOWN)),
                    ("Original-hit shape", row.get("hit_shape", UNKNOWN)),
                    ("Original-hit pharmacophore", row.get("hit_pharmacophore", UNKNOWN)),
                    ("Polar recall", _recall_text(row.get("polar_reference_recall"), candidate)),
                    ("Score", row.get("score", UNKNOWN)),
                ]
            ),
            "<h3>Required final seeds</h3><p>Each seed has its own same-pose witness. Different seeds need not share a pose. Observed seed spread is not calibrated uncertainty.</p>",
        ]
        seeds = candidate["seed_rows"]
        parts.append(
            f"<p>Requested seeds: {_escape(_mapping(policy.get('refinement')).get('seeds', UNKNOWN))}</p>"
        )
        if seeds:
            columns = ("seed", "status", "hit_shape", "hit_pharmacophore", "polar_recall", "score")
            cells = [
                [
                    (
                        "not applicable"
                        if s.get("polar_status") == "not_applicable" and s.get(k) is None
                        else s.get(k, UNKNOWN)
                    )
                    for k in columns
                ]
                for s in seeds
            ]
            parts.append(
                "<div class='scroll'>"
                + _table(
                    cells,
                    ("Seed", "Recorded status", "Hit shape", "Hit pharmacophore", "Polar recall", "Score"),
                )
                + "</div>"
            )
        parts += [
            "<h4>Derived arithmetic: minima across complete saved final seeds</h4>",
            _fields(
                candidate["derived_minima"] or {"minima": "unavailable; incomplete or inconsistent evidence"}
            ),
            _route_html(candidate["route"]),
        ]
        bindings = _fields(
            {
                k: row.get(k, UNKNOWN)
                for k in (
                    "route_id",
                    "alternative_route_ids",
                    "pose_seed",
                    "relative_energy",
                    "design_shape",
                    "design_pharmacophore",
                    "scoring_objective",
                    "ranking_policy",
                    "similarity_witness_sha256",
                )
            }
        )
        bindings += "<p>Alternative routes are separate records; their support is not borrowed for this primary route.</p>"
        for seed in seeds:
            bindings += _fields(
                {
                    k: seed.get(k, UNKNOWN)
                    for k in (
                        "seed",
                        "reason",
                        "polar_status",
                        "witness_sha256",
                        "candidate_pose_sha256",
                        "reference_pose_sha256",
                    )
                }
            )
        parts += [
            _details("Saved original-hit identity comparison", _fields(identity)),
            _details("Pose and route bindings (recorded)", bindings),
            "</article>",
        ]
    if not view["candidates"]:
        parts.append(
            "<h2>No selected candidates</h2><p>This is not proof that no suitable molecule exists. Only recorded run diagnostics are shown.</p>"
        )
    parts += [
        "<h2>Recorded diagnostics</h2><p>route_rejections counts route records, not unique molecules. Product, family and proposal counts retain their recorded names. Missing or unevaluated work is not a chemical failure.</p>",
        _fields(
            {
                k: summary.get(k, UNKNOWN)
                for k in (
                    "status",
                    "no_candidates_reason",
                    "error_type",
                    "error",
                    "returned_products",
                    "accepted_products",
                    "qualified_products",
                    "deadline_reached",
                    "limits_reached",
                    "route_rejections",
                    "forward_evidence_outcomes",
                )
            }
        ),
    ]
    failures = _dicts(summary.get("failures"))
    if failures:
        parts.append("<h3>Recorded failures</h3>")
        parts.extend(_fields(item) for item in failures)
    parts += [
        _details("Complete recorded run diagnostics", _fields(summary)),
        _details("Effective policy and provenance", _fields(policy)),
        _details(
            "Saved identity descriptor (not source authentication)",
            _fields(_mapping(_mapping(view.get("identity_descriptor")).get("record"))),
        ),
        "<details><summary>Rendered artifact bytes and local files</summary><p>Hashes identify rendered bytes; they do not authenticate a self-consistent bundle. The file inventory identifies SDF/TSV artifacts for review; they are not parsed or scientifically validated by this report.</p>",
    ]
    for name, item in sorted(view["artifacts"].items()):
        label = f"<a href='{links[name]}'>{_escape(name)}</a>" if name in links else _escape(name)
        parts.append(
            f"<p>{label} — {_escape(item['bytes'])} bytes — SHA256 <code>{item['sha256']}</code></p>"
        )
    parts.append("</details>")
    if any(
        c.get("visuals", {}).get("status") == "available" and c["visuals"].get("seeds")
        for c in view["candidates"]
    ):
        script = (Path(__file__).parent / "data" / "review_viewer.js").read_text(encoding="utf-8")
        parts.append("<script>" + script + "</script>")
    parts.append("</body></html>\n")
    output = "".join(parts)
    style = output.split("<style>", 1)[1].split("</style>", 1)[0]
    hashes = lambda text: base64.b64encode(hashlib.sha256(text.encode()).digest()).decode()
    scripts = "'sha256-" + hashes(script) + "'" if "script" in locals() else "'none'"
    csp = f"default-src 'none'; img-src data:; style-src 'sha256-{hashes(style)}'; script-src {scripts}; connect-src 'none'; object-src 'none'; base-uri 'none'; form-action 'none'"
    return output.replace(
        "<head>", '<head><meta http-equiv="Content-Security-Policy" content="' + csp + '">', 1
    )


def _recall_text(value, candidate):
    if value is None and any(
        s.get("seed") == candidate["row"].get("pose_seed") and s.get("polar_status") == "not_applicable"
        for s in candidate["seed_rows"]
    ):
        return "not applicable (no polar reference features)"
    return value if value is not None else UNKNOWN


def write_report(run_dir, output=None):
    """Create a new report; refuse overwrite and changes to the saved input snapshot.

    Output within the run directory receives safe local artifact links. An output
    elsewhere shows filenames/hashes without links. Returns the consistency view.
    """
    root = Path(run_dir).resolve(strict=True)
    target = Path(output) if output is not None else root / "review.html"
    target = target.absolute()
    if target.exists() or target.is_symlink() or target.name in JSON_FILES + LINK_FILES:
        raise ValueError("Report output must be a new file, not an input artifact")
    docs, inventory, issues, observed = _load(root)
    view = _assemble(docs, inventory, issues)
    links = {name: quote(name) for name in inventory} if target.parent.resolve() == root else {}
    rendered = render_html(view, artifact_links=links)
    _unchanged(root, observed)
    created = False
    try:
        with target.open("x", encoding="utf-8", newline="\n") as stream:
            created = True
            stream.write(rendered)
            stream.flush()
            os.fsync(stream.fileno())
        _unchanged(root, observed)
    except BaseException:
        if created:
            target.unlink(missing_ok=True)
        raise
    return view


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        view = write_report(args.run_dir, args.output)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"Report not written: {exc}\n")
    print(f"Artifact consistency: {view['artifact_consistency']}; no scientific revalidation")
    return 0 if view["artifact_consistency"] == "consistent_recorded_bundle" else 3


if __name__ == "__main__":
    raise SystemExit(main())
