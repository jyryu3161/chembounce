"""Validate the complete frozen grid, then independently summarize forward evidence.

No model inference is performed. Incomplete/invalid inputs produce validation
diagnostics only. The live audit program, caches and production code are read-only.
"""

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import statistics
import sys

import numpy as np


RUNTIME = Path(__file__).resolve().parent
BASE = RUNTIME.parent
ROOT = RUNTIME.parents[2]
sys.path.insert(0, str(ROOT))
PLAN_SHA256 = "385be03fee3d2df6cfc04975cb63b42717088eabe62077bdbd3c6db9ca5db529"
# Canonical sorted-JSON code_sha256 mapping from the original frozen campaign
# record, anchored independently of mutable copies in campaign manifests.
FROZEN_CODE_INVENTORY_SHA256 = "c3a926f289accdbea57c84923b82ba5863e9bbb0a5b87a73f2081e798b5cfe07"
FROZEN_CODE_FILE_COUNT = 39
FROZEN_IMPLEMENTATION_COMMIT = "d972124d587e647ef48455cf1a038e6be96e5b04"
SPLITS = ("scaffold", "time", "target")
PRIMARY = "attempt_with_common_forward_matched_route"
INFRA = "infrastructure_interruption"
COUNTS = (
    "candidate_occurrences",
    "no_route_supplied",
    "route_supplied",
    "route_assessed",
    "route_all_steps_evaluated",
    "route_all_steps_matched",
    "route_assessed_incomplete",
    "route_evaluated_mismatched",
    "route_audit_failed",
)
COVERAGE = ("route_supplied", "route_assessed", "route_all_steps_evaluated", "route_all_steps_matched")
BOUNDED = (
    (PRIMARY,)
    + tuple(f"attempt_with_{name}" for name in COVERAGE)
    + tuple(f"{name}_output_slot_fraction" for name in COVERAGE)
)
CONDITIONAL = "matched_fraction_given_evaluable_routes"
METRICS = COUNTS + BOUNDED + (CONDITIONAL,)
TERMINAL = {
    "success",
    "no_candidates",
    "timeout",
    "memory_limit",
    "budget_exceeded",
    "worker_error",
    "missing_response",
    "invalid_response",
    "launch_error",
    "interrupted",
    INFRA,
}


class ValidationError(ValueError):
    def __init__(self, message, details=None):
        self.details = details or {}
        super().__init__(message)


def require(condition, message, details=None):
    if not condition:
        raise ValidationError(message, details)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def content_key(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def read_json(path, hashes=None, label=None):
    path = Path(path)
    raw = path.read_bytes()
    if hashes is not None:
        hashes[label or str(path)] = hashlib.sha256(raw).hexdigest()

    def reject(value):
        raise ValidationError(f"Nonfinite JSON constant in {path}: {value}")

    return json.loads(raw, parse_constant=reject)


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def run_key(row):
    return row["method_id"], row["query_id"], row["seed"]


def occurrence_key(row):
    return row["split"], row["method_id"], row["query_id"], row["generation_seed"], row["rank"]


def require_frozen_inventory(inventory, label):
    require(
        isinstance(inventory, dict)
        and len(inventory) == FROZEN_CODE_FILE_COUNT
        and content_key(inventory) == FROZEN_CODE_INVENTORY_SHA256,
        f"Frozen code inventory missing or differs from original 39-file inventory: {label}",
    )


def validate_grid(base, audit_dir, plan):
    """Only JSON/hash/grid checks here: no RDKit, catalog, rules or model loading."""
    base, audit_dir = Path(base), Path(audit_dir)
    hashes, panels, protocols, rows, paths, grids = {}, {}, {}, {}, {}, {}
    provenance_path = base / "frozen-campaign-code-provenance.json"
    provenance = read_json(provenance_path, hashes)
    inventory = provenance.get("code_sha256")
    require_frozen_inventory(inventory, str(provenance_path))
    require(
        provenance.get("file_count") == FROZEN_CODE_FILE_COUNT
        and provenance.get("implementation_commit") == FROZEN_IMPLEMENTATION_COMMIT,
        "Frozen code provenance count/implementation metadata mismatch",
    )
    require(set(plan["expected_grid"]) == set(SPLITS), "Plan must cover all three splits")
    for name, sha in plan["source_sha256"].items():
        require(digest(base / name) == sha, f"Frozen plan source hash mismatch: {name}")
        hashes[str(base / name)] = sha
    for split in SPLITS:
        panel = read_json(base / f"panels/direct-human-{split}/panel.json")
        records = {r["id"]: r for r in panel["compounds"]}
        targets = {q["id"]: records[q["compound_id"]]["target_id"] for q in panel["queries"]}
        require(len(targets) == len(panel["queries"]), "Duplicate panel query")
        require(targets == plan["expected_grid"][split]["queries"], "Query/target mapping differs from plan")
        panels[split] = panel
        expected, actual_rows, split_paths = set(), [], {}
        for group, spec in plan["expected_grid"][split]["groups"].items():
            protocol_path = base / spec["protocol"]
            protocol = read_json(protocol_path)
            protocols[split, group] = protocol
            require([m["id"] for m in protocol["methods"]] == spec["methods"], "Method plan mismatch")
            require(protocol["seeds"] == spec["seeds"], "Generation seed plan mismatch")
            group_grid = {(m, q, seed) for m in spec["methods"] for q in targets for seed in spec["seeds"]}
            require(
                not expected.intersection(group_grid), "Duplicate planned method/query/seed across groups"
            )
            expected.update(group_grid)
            folder = f"comparison-{group}" if split == "scaffold" else f"comparison-{split}-{group}"
            run_path = base / folder / "runs.json"
            if not run_path.exists():
                continue
            manifest_path = run_path.parent / "manifest.json"
            manifest = read_json(manifest_path, hashes)
            require_frozen_inventory(manifest.get("code_sha256"), folder)
            require(
                manifest["protocol_sha256"] == digest(protocol_path), f"Campaign protocol mismatch: {folder}"
            )
            require(
                manifest["dataset"]["dataset_sha256"]
                == digest(base / f"panels/direct-human-{split}/panel.json"),
                f"Campaign panel mismatch: {folder}",
            )
            group_rows = read_json(run_path, hashes)
            require(isinstance(group_rows, list), "runs.json must be a list")
            for row in group_rows:
                require(type(row["seed"]) is int, "Generation seed must be integer")
                key = run_key(row)
                require(key in group_grid, f"Unexpected attempt in {folder}: {key}")
                require(row["target_id"] == targets[row["query_id"]], "Run target/query mismatch")
                require(row["status"] in TERMINAL, f"Unknown terminal status: {row['status']}")
                split_paths[key] = run_path.parent
            actual_rows.extend(group_rows)
        actual = Counter(run_key(r) for r in actual_rows)
        missing = sorted(expected - actual.keys())
        duplicates = [list(k) for k, n in actual.items() if n != 1]
        require(
            len(expected) == plan["expected_attempts"][split], "Expected attempt total disagrees with plan"
        )
        grids[split] = {
            "expected": len(expected),
            "observed": len(actual_rows),
            "missing": missing,
            "duplicates": duplicates,
            "complete": not missing and not duplicates,
        }
        rows[split], paths[split] = actual_rows, split_paths
    require(
        all(g["complete"] for g in grids.values()),
        "Incomplete/invalid global attempt grid; no comparative estimates emitted",
        {"grid": grids},
    )
    protocol_path = audit_dir / "protocol.json"
    require(digest(protocol_path) == plan["audit_protocol_sha256"], "Common forward protocol hash mismatch")
    protocol = read_json(protocol_path, hashes)
    require(protocol["expected_attempts"] == plan["expected_attempts"], "Audit attempt plan mismatch")
    progress = read_json(audit_dir / "progress.json", hashes)
    require(progress.get("complete") is True, "Common forward assessment is incomplete")
    require(progress["observed_attempts"] == plan["expected_attempts"], "Audit observed grid count mismatch")
    required_folders = {str(p.relative_to(base)) for split in paths.values() for p in split.values()}
    require(set(progress["input_run_sha256"]) == required_folders, "Missing/extraneous audit run-file hash")
    for folder, expected_sha in progress["input_run_sha256"].items():
        require(
            hashes[str(base / folder / "runs.json")] == expected_sha, f"Audit input hash mismatch: {folder}"
        )
    return {
        "rows": rows,
        "paths": paths,
        "panels": panels,
        "protocols": protocols,
        "grids": grids,
        "protocol": protocol,
        "progress": progress,
        "hashes": hashes,
        "frozen_code_sha256": inventory,
    }


def exact_smiles(smiles, *, mixture=False):
    """Preserve exact state/stereo; multi-component predictions remain whole."""
    from rdkit import Chem

    require(isinstance(smiles, str) and bool(smiles.strip()), "Missing SMILES")
    mol = Chem.MolFromSmiles(smiles)
    require(
        mol is not None and mol.GetNumAtoms() > 0 and all(a.GetAtomicNum() for a in mol.GetAtoms()),
        f"Invalid SMILES: {smiles}",
    )
    require(mixture or len(Chem.GetMolFrags(mol)) == 1, "Disconnected candidate/reactant")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def validate_step(step, replay, protocol):
    require(step["node_id"] == replay["node_id"], "Evidence/replay step node mismatch")
    require(step["product"] == replay["product"], "Evidence/replay step product mismatch")
    require(step["reactants"] == replay["reactants"], "Evidence/replay reactants mismatch")
    require(step["rule_replay"] == "matched", "Evidence does not declare rule replay")
    require(
        step["context"] == {"status": "incomplete", "supplied": {}, "validated": False},
        "Common reactant-only context changed",
    )
    require(
        step["observations"]["records"] == [] and step["observations"]["status"] == "unknown",
        "Unexpected empirical observations in common audit",
    )
    prediction = step["forward_prediction"]
    require(
        prediction["input_coverage"]
        == {
            "conditions_supplied": False,
            "conditions_used_by_predictor": False,
            "incomplete_context": True,
            "reagents_supplied": False,
        },
        "Prediction input context coverage mismatch",
    )
    values = prediction.get("predictions", [])
    require(len(values) <= protocol["model"]["top_k"], "Prediction exceeds pinned top-k")
    canonical, matched_ranks = [], []
    for rank, value in enumerate(values, 1):
        require(type(value["rank"]) is int and value["rank"] == rank, "Prediction rank mismatch")
        require(type(value["complete"]) is bool, "Prediction completion must be boolean")
        require(isinstance(value["raw_product"], str), "Raw predicted product must be a string")
        product = None
        if value["complete"]:
            try:
                product = exact_smiles(value["raw_product"], mixture=True)
            except ValidationError:
                pass
        require(value["product"] == product, "Prediction raw/canonical product mismatch")
        canonical.append(product)
        if product == step["product"]:
            matched_ranks.append(rank)
    status = prediction["status"]
    require(status in {"matched", "mismatched", "unknown", "failed"}, "Unknown prediction status")
    if any(canonical):
        require(
            status == ("matched" if matched_ranks else "mismatched"), "Prediction agreement status mismatch"
        )
        require(
            prediction["matched_rank"] == (min(matched_ranks) if matched_ranks else None),
            "Prediction matched rank mismatch",
        )
        require(prediction["matching_policy"] == "exact", "Nonexact product matching policy")
    else:
        require(status in {"unknown", "failed"}, "Unevaluable predictions promoted to matched/mismatched")
        require(bool(values) or bool(prediction.get("reason")), "Missing incomplete prediction explanation")
    return status, content_key([sorted(step["reactants"]), [], {}])


def validate_cache(cached, raw_route, product, protocol, replay):
    if cached["status"] == "audit_failed":
        require(isinstance(cached.get("error"), str) and bool(cached["error"]), "Audit failure lacks reason")
        require(
            cached.get("forward_all_steps_matched") is None and cached.get("all_steps_evaluated") is None,
            "Audit failure promoted to evaluated/matched",
        )
        return {"status": "audit_failed", "matched": False, "evaluated": False, "input_keys": []}
    require(cached["status"] == "assessed", "Unknown route cache status")
    try:
        verification = replay(raw_route)
    except (ValueError, RuntimeError) as exc:
        raise ValidationError(f"Retained assessed route cannot be replayed: {exc}") from exc
    require(
        verification == cached["verification"], "Fresh raw route replay differs from retained verification"
    )
    require(verification["product"] == product, "Replayed route product differs from candidate")
    require(
        verification["catalog_sha256"] == protocol["catalog_sha256"]
        and verification["rules_sha256"] == protocol["rules_sha256"],
        "Route replay asset hash mismatch",
    )
    evidence = cached["evidence"]
    require(evidence["product"] == product, "Evidence root product mismatch")
    require(
        evidence["predictor_provenance"] == protocol["model"], "Evidence predictor/model provenance mismatch"
    )
    require(evidence["policy"] == protocol["evidence_policy"]["policy"], "Evidence policy mismatch")
    steps = evidence["steps"]
    require(bool(steps) and len(steps) == len(verification["steps"]), "Missing/extraneous evidence step")
    states, input_keys = [], []
    for step, verified in zip(steps, verification["steps"]):
        status, input_key = validate_step(step, verified, protocol)
        states.append(status)
        input_keys.append(input_key)
    matched = all(s == "matched" for s in states)
    evaluated = all(s in {"matched", "mismatched"} for s in states)
    require(
        cached["forward_all_steps_matched"] is matched and cached["all_steps_evaluated"] is evaluated,
        "Cached route summary disagrees with prediction steps",
    )
    require(evidence["eligibility"]["passed"] is matched, "Evidence eligibility disagrees with pinned policy")
    return {"status": "assessed", "matched": matched, "evaluated": evaluated, "input_keys": input_keys}


def resolve_replay_inputs(base, snapshot, source_root=None, model_directory=None):
    """Relocate files while preserving the original protocol bytes and pins."""
    source_root = Path(source_root or ROOT).resolve()
    inventory = snapshot.get("frozen_code_sha256")
    require_frozen_inventory(inventory, "validated snapshot")
    for name, sha in inventory.items():
        path = source_root / name
        require(path.is_file(), f"Frozen campaign source file missing: {name}")
        require(digest(path) == sha, f"Frozen campaign source hash mismatch: {name}")
        snapshot["hashes"][str(path)] = sha
    protocol = snapshot["protocol"]
    for name, sha in protocol["source_sha256"].items():
        require(digest(source_root / name) == sha, f"Audit implementation source hash mismatch: {name}")
        snapshot["hashes"][str(source_root / name)] = sha
    model = protocol["model"]
    model_root = Path(model_directory or model["model_path"]).resolve()
    for name, sha in model["model_files_sha256"].items():
        path = model_root / name
        require(digest(path) == sha, f"Model artifact hash mismatch: {name}")
        snapshot["hashes"][str(path)] = sha
    assets = next(
        m["assets"] for p in snapshot["protocols"].values() for m in p["methods"] if m["id"] == "pipeline"
    )
    asset_paths = {}
    for name in ("catalog", "rules"):
        path = Path(base) / assets[name]["artifact"]
        require(
            digest(path) == assets[name]["sha256"] == protocol[f"{name}_sha256"],
            f"Replay {name} artifact hash mismatch",
        )
        snapshot["hashes"][str(path)] = assets[name]["sha256"]
        asset_paths[name] = path
    snapshot["resolved_inputs"] = {
        "source_root": str(source_root),
        "model_directory": str(model_root),
        "original_model_path_provenance": model["model_path"],
        "frozen_code_inventory_sha256": FROZEN_CODE_INVENTORY_SHA256,
        "frozen_code_file_count": len(inventory),
        "replay_assets": {k: str(v) for k, v in asset_paths.items()},
    }
    return source_root, asset_paths


def build_replay(base, snapshot, source_root=None, model_directory=None):
    """Called once, after complete JSON/hash/grid validation; never loads a model."""
    source_root, asset_paths = resolve_replay_inputs(base, snapshot, source_root, model_directory)
    sys.path.insert(0, str(source_root))
    try:
        from chembounce_synthesis import catalog as catalog_module, reactions, routes
    finally:
        sys.path.pop(0)
    for name, module in tuple(sys.modules.items()):
        if name.startswith("chembounce_synthesis") and getattr(module, "__file__", None):
            require(
                Path(module.__file__).resolve().is_relative_to(source_root),
                "Imported synthesis module differs from validated source root; use a fresh interpreter",
            )
    Catalog, Rules, verify_route = catalog_module.Catalog, reactions.Rules, routes.verify_route
    catalog, rules = Catalog(asset_paths["catalog"]), Rules(asset_paths["rules"])
    return lambda route: verify_route(route, catalog, rules)


def validate_occurrences(base, audit_dir, plan, snapshot, replay):
    hashes = snapshot["hashes"]
    occurrences = read_json(Path(audit_dir) / "occurrences.json", hashes)
    require(isinstance(occurrences, list), "Occurrences must be a list")
    require(
        all(
            type(r["rank"]) is int and r["rank"] > 0 and type(r["generation_seed"]) is int
            for r in occurrences
        ),
        "Occurrence rank/seed must be integers",
    )
    indexed = {occurrence_key(row): row for row in occurrences}
    require(len(indexed) == len(occurrences), "Duplicate audited occurrence")
    require(
        snapshot["progress"]["candidate_occurrences"] == len(occurrences),
        "Progress occurrence count mismatch",
    )
    require(
        snapshot["progress"]["status_counts"] == dict(Counter(r["status"] for r in occurrences)),
        "Progress occurrence status count mismatch",
    )
    route_files = list((Path(audit_dir) / "routes").glob("*.json"))
    require(
        snapshot["progress"]["unique_cached_routes"] == len(route_files),
        "Progress route cache count mismatch",
    )
    seen, cache, annotated = set(), {}, defaultdict(list)
    for split in SPLITS:
        panel = snapshot["panels"][split]
        query_indices = {q["id"]: i for i, q in enumerate(panel["queries"])}
        for run in snapshot["rows"][split]:
            rk = run_key(run)
            folder = snapshot["paths"][split][rk]
            attempt = (
                folder / run["method_id"] / f"query-{query_indices[run['query_id']]:04d}-seed-{run['seed']}"
            )
            require(
                read_json(attempt / "evaluation.json", hashes) == run,
                "Attempt evaluation differs from runs.json",
            )
            accepted = run["status"] in {"success", "no_candidates"}
            response_sha = run.get("response_sha256")
            if accepted:
                require(
                    isinstance(response_sha, str)
                    and len(response_sha) == 64
                    and all(c in "0123456789abcdef" for c in response_sha),
                    "Accepted response requires response_sha256",
                )
            if response_sha is not None:
                response_path = attempt / "response.json"
                require(digest(response_path) == response_sha, "Raw response hash mismatch")
                hashes[str(response_path)] = response_sha
                if accepted:
                    require(
                        read_json(response_path) == run["response"], "Accepted response differs from raw file"
                    )
            candidates = run["evaluation"]["candidates"]
            require(
                run["status"] == "success" or not candidates,
                "Failed/no-output attempt has accepted candidates",
            )
            raw_candidates = run.get("response", {}).get("candidates", [])
            raw_ranks = set()
            for rank, candidate in enumerate(candidates, 1):
                original_rank = candidate["rank"]
                require(
                    type(original_rank) is int and 1 <= original_rank <= len(raw_candidates),
                    "Evaluation candidate raw rank is invalid",
                )
                require(original_rank not in raw_ranks, "Duplicate raw candidate rank in evaluation")
                raw_ranks.add(original_rank)
                raw = raw_candidates[original_rank - 1]
                require(
                    exact_smiles(raw["smiles"]) == candidate["smiles"],
                    "Raw/evaluation candidate identity mismatch",
                )
                ok = (split, *rk, rank)
                require(ok in indexed, "Missing audited candidate occurrence")
                seen.add(ok)
                occurrence = indexed[ok]
                require(
                    occurrence["target_id"] == run["target_id"]
                    and occurrence["smiles"] == candidate["smiles"],
                    "Occurrence target/candidate identity mismatch",
                )
                route = raw.get("route")
                if not route:
                    require(
                        occurrence["status"] == "no_route_supplied"
                        and occurrence.get("forward_all_steps_matched") is None
                        and occurrence.get("all_steps_evaluated") is None
                        and "route_cache_key" not in occurrence,
                        "No-route candidate promoted to route evidence",
                    )
                    outcome = {
                        "status": "no_route_supplied",
                        "matched": False,
                        "evaluated": False,
                        "input_keys": [],
                    }
                    route_key = None
                else:
                    route_key = content_key([plan["audit_protocol_sha256"], route, candidate["smiles"]])
                    require(
                        occurrence.get("route_cache_key") == route_key,
                        "Route cache key differs from raw route/product/protocol",
                    )
                    if route_key not in cache:
                        cached = read_json(Path(audit_dir) / "routes" / f"{route_key}.json", hashes)
                        cache[route_key] = (
                            cached,
                            validate_cache(cached, route, candidate["smiles"], snapshot["protocol"], replay),
                        )
                    cached, outcome = cache[route_key]
                    require(
                        occurrence["status"] == cached["status"]
                        and occurrence.get("forward_all_steps_matched")
                        is cached.get("forward_all_steps_matched")
                        and occurrence.get("all_steps_evaluated") is cached.get("all_steps_evaluated")
                        and occurrence.get("error") == cached.get("error"),
                        "Occurrence/cache result mismatch",
                    )
                annotated[split].append({**occurrence, **outcome, "route_key": route_key})
    require(seen == indexed.keys(), "Extraneous audited occurrence")
    return dict(annotated)


def mean(values):
    values = [v for v in values if v is not None]
    return statistics.mean(values) if values else None


def bootstrap(values, config):
    values = np.array([v for v in values if v is not None], dtype=float)
    result = {
        "mean": float(values.mean()) if len(values) else None,
        "target_clusters": len(values),
        "ci95": None,
    }
    if len(values) > 1:
        draws = (
            np.random.default_rng(config["seed"])
            .choice(values, (config["draws"], len(values)), replace=True)
            .mean(axis=1)
        )
        result["ci95"] = np.quantile(draws, config["quantiles"]).tolist()
    return result


def summarize(rows, targets, names, config, *, bounds=False, planned_attempts=None):
    planned_attempts = len(rows) if planned_attempts is None else planned_attempts
    per_query = {
        q: {name: mean(r["metrics"].get(name) for r in rows if r["query_id"] == q) for name in names}
        for q in targets
    }
    per_target = {
        t: {name: mean(per_query[q][name] for q in targets if targets[q] == t) for name in names}
        for t in sorted(set(targets.values()))
    }
    eligibility = {
        name: {
            "attempts": sum(r["metrics"].get(name) is not None for r in rows),
            "queries": sum(q[name] is not None for q in per_query.values()),
            "targets": sum(t[name] is not None for t in per_target.values()),
        }
        for name in names
    }
    for values in eligibility.values():
        values.update(
            excluded_attempts=planned_attempts - values["attempts"],
            excluded_queries=len(targets) - values["queries"],
            excluded_targets=len(set(targets.values())) - values["targets"],
        )
    return {
        "attempts": len(rows),
        "planned_attempts": planned_attempts,
        "excluded_attempts": planned_attempts - len(rows),
        "metric_eligibility": eligibility,
        "per_query_attempts": {q: sum(r["query_id"] == q for r in rows) for q in targets},
        "per_query": per_query,
        "per_target": per_target,
        "target_macro": {
            name: (
                mean(t[name] for t in per_target.values())
                if bounds
                else bootstrap((t[name] for t in per_target.values()), config)
            )
            for name in names
        },
    }


def paired_rows(left, right, *, complete_case=False):
    a = {(r["query_id"], r["seed"]): r for r in left}
    b = {(r["query_id"], r["seed"]): r for r in right}
    require(
        a.keys() == b.keys() and len(a) == len(left) and len(b) == len(right),
        "Paired query/seed grids differ",
    )
    rows = []
    for key, row in a.items():
        other = b[key]
        if complete_case and INFRA in (row["status"], other["status"]):
            continue
        values = {
            name: (
                row["metrics"][name] - other["metrics"][name]
                if row["metrics"][name] is not None and other["metrics"][name] is not None
                else None
            )
            for name in METRICS
        }
        rows.append({**row, "metrics": values})
    return rows


def aggregate_split(runs, occurrences, targets, caps, config):
    """Independent seed/query/target aggregation; callable on synthetic records."""
    by_attempt = defaultdict(list)
    for row in occurrences:
        by_attempt[row["method_id"], row["query_id"], row["generation_seed"]].append(row)
    methods = defaultdict(list)
    occurrence_methods = defaultdict(list)
    for run in runs:
        output = by_attempt[run_key(run)]
        require(run["status"] == "success" or not output, "Failed attempt has forward occurrences")
        c = Counter(r["status"] for r in output)
        evaluated, matched = sum(r["evaluated"] for r in output), sum(r["matched"] for r in output)
        metrics = {
            "candidate_occurrences": len(output),
            "no_route_supplied": c["no_route_supplied"],
            "route_supplied": len(output) - c["no_route_supplied"],
            "route_assessed": c["assessed"],
            "route_all_steps_evaluated": evaluated,
            "route_all_steps_matched": matched,
            "route_assessed_incomplete": c["assessed"] - evaluated,
            "route_evaluated_mismatched": evaluated - matched,
            "route_audit_failed": c["audit_failed"],
        }
        cap = caps[run["method_id"]]
        require(
            0 <= len(output) <= cap and min(metrics.values()) >= 0, "Invalid coverage counts/candidate cap"
        )
        metrics[PRIMARY] = int(matched > 0)
        metrics.update({f"attempt_with_{name}": int(metrics[name] > 0) for name in COVERAGE})
        metrics.update({f"{name}_output_slot_fraction": metrics[name] / cap for name in COVERAGE})
        metrics[CONDITIONAL] = matched / evaluated if evaluated else None
        methods[run["method_id"]].append(
            {
                "query_id": run["query_id"],
                "target_id": run["target_id"],
                "seed": run["seed"],
                "status": run["status"],
                "metrics": metrics,
            }
        )
        occurrence_methods[run["method_id"]].extend(output)
    result, complete, bounds, paired, paired_complete, paired_bounds = {}, {}, {}, {}, {}, {}
    for method, rows in sorted(methods.items()):
        counts = {name: sum(r["metrics"][name] for r in rows) for name in COUNTS}
        output = occurrence_methods[method]
        input_keys = [key for row in output for key in row["input_keys"]]
        result[method] = {
            "status_counts": dict(Counter(r["status"] for r in rows)),
            "coverage_counts": counts,
            "unique_supplied_route_keys": len({r["route_key"] for r in output if r["route_key"]}),
            "route_step_input_occurrences": len(input_keys),
            "unique_reactant_context_input_keys": len(set(input_keys)),
            "pooled_occurrence_matched_fraction_selected_subset": (
                counts["route_all_steps_matched"] / counts["route_all_steps_evaluated"]
                if counts["route_all_steps_evaluated"]
                else None
            ),
            **summarize(rows, targets, METRICS, config),
        }
        observed = [r for r in rows if r["status"] != INFRA]
        complete[method] = {
            "excluded_infrastructure_attempts": len(rows) - len(observed),
            **summarize(observed, targets, METRICS, config, planned_attempts=len(rows)),
        }
        bound_rows = [
            {
                **r,
                "metrics": {
                    f"{name}_{side}": (value if r["status"] != INFRA else missing)
                    for name in BOUNDED
                    for side, missing in [("lower", 0), ("upper", 1)]
                    for value in [r["metrics"][name]]
                },
            }
            for r in rows
        ]
        bound_names = [f"{name}_{side}" for name in BOUNDED for side in ("lower", "upper")]
        bounds[method] = summarize(bound_rows, targets, bound_names, config, bounds=True)
    if "pipeline" in methods:
        for method, rows in sorted(methods.items()):
            if method == "pipeline":
                continue
            primary = methods["pipeline"]
            paired[method] = summarize(paired_rows(primary, rows), targets, METRICS, config)
            paired_complete[method] = summarize(
                paired_rows(primary, rows, complete_case=True),
                targets,
                METRICS,
                config,
                planned_attempts=len(primary),
            )
            lookup = {(r["query_id"], r["seed"]): r for r in rows}
            differences = []
            for row in primary:
                other = lookup[row["query_id"], row["seed"]]
                values = {}
                for name in BOUNDED:
                    low_a, high_a = (0, 1) if row["status"] == INFRA else (row["metrics"][name],) * 2
                    low_b, high_b = (0, 1) if other["status"] == INFRA else (other["metrics"][name],) * 2
                    values[f"{name}_lower"], values[f"{name}_upper"] = low_a - high_b, high_a - low_b
                differences.append({**row, "metrics": values})
            paired_bounds[method] = summarize(differences, targets, bound_names, config, bounds=True)
    input_targets, input_queries = defaultdict(set), defaultdict(set)
    for occurrence in occurrences:
        for key in occurrence["input_keys"]:
            input_targets[key].add(occurrence["target_id"])
            input_queries[key].add(occurrence["query_id"])
    dependence = {
        "unique_reactant_context_input_keys": len(input_targets),
        "input_keys_shared_across_targets": sum(len(values) > 1 for values in input_targets.values()),
        "input_keys_shared_across_queries": sum(len(values) > 1 for values in input_queries.values()),
        "uncertainty_policy": "Occurrences/routes/reactant reuse are not independent samples; target-cluster intervals are descriptive, with cross-target shared-input dependence and few-target limitations.",
    }
    return {
        "methods": result,
        "paired_pipeline_minus_comparator": paired,
        "dependence_diagnostics": dependence,
        "infrastructure_identification_bounds_not_ci": {
            "methods": bounds,
            "paired_pipeline_minus_comparator": paired_bounds,
        },
        "complete_case_secondary_selection_biased": {
            "methods": complete,
            "paired_pipeline_minus_comparator": paired_complete,
        },
    }


def analyze(base, audit_dir, plan, *, validate_only=False, source_root=None, model_directory=None):
    snapshot = validate_grid(base, audit_dir, plan)
    replay = build_replay(base, snapshot, source_root, model_directory)
    occurrences = validate_occurrences(base, audit_dir, plan, snapshot, replay)
    validated = {
        "passed": True,
        "grid": snapshot["grids"],
        "source_sha256": snapshot["hashes"],
        "total_attempts": sum(len(rows) for rows in snapshot["rows"].values()),
        "resolved_inputs": snapshot["resolved_inputs"],
        "audit_model_provenance": snapshot["protocol"]["model"],
        "analysis_runtime": {
            "python": sys.version,
            "numpy": np.__version__,
            "rdkit": importlib.metadata.version("rdkit"),
            "forward_inference_performed": False,
        },
        "analysis_helper_sha256": digest(Path(__file__)),
        "scientific_comparison_emitted": not validate_only,
    }
    if validate_only:
        return validated, None
    results = {}
    for split in SPLITS:
        settings = plan["expected_grid"][split]
        caps = {m: cap for spec in settings["groups"].values() for m, cap in spec["output_caps"].items()}
        results[split] = aggregate_split(
            snapshot["rows"][split], occurrences.get(split, []), settings["queries"], caps, plan["bootstrap"]
        )
    return validated, {
        "analysis_plan": plan,
        "validation": validated,
        "splits": results,
        "interpretation": plan["interpretation"],
        "dependence": plan["dependence"],
        "conditional_metric": plan["conditional_metric"],
        "infrastructure": plan["infrastructure"],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--audit-dir", type=Path)
    parser.add_argument("--plan", type=Path, default=RUNTIME / "forward-analysis-plan.json")
    parser.add_argument(
        "--source-root",
        type=Path,
        default=ROOT,
        help="Frozen source checkout; all pinned source hashes must still match",
    )
    parser.add_argument(
        "--model-directory",
        type=Path,
        help="Relocated model files; every pinned model file hash must still match",
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument(
        "--validation-output", type=Path, default=RUNTIME / "forward-analysis-validation.json"
    )
    parser.add_argument("--output", type=Path, default=RUNTIME / "forward-analysis.json")
    args = parser.parse_args(argv)
    try:
        require(digest(args.plan) == PLAN_SHA256, "Analysis plan differs from fixed recorded plan")
        plan = read_json(args.plan)
        validation, result = analyze(
            args.base,
            args.audit_dir or args.base / "common-forward-audit-cpu",
            plan,
            validate_only=args.validate_only,
            source_root=args.source_root,
            model_directory=args.model_directory,
        )
    except (ValidationError, OSError, KeyError, TypeError, json.JSONDecodeError) as exc:
        validation = {
            "passed": False,
            "error": f"{type(exc).__name__}: {exc}",
            "details": getattr(exc, "details", {}),
            "scientific_comparison_emitted": False,
            "checked_utc": datetime.now(timezone.utc).isoformat(),
        }
        save(args.validation_output, validation)
        print(json.dumps({"passed": False, "error": str(exc), "scientific_comparison_emitted": False}))
        return 2
    validation["analysis_plan_sha256"] = digest(args.plan)
    save(args.validation_output, validation)
    if result is not None:
        save(args.output, result)
    print(
        json.dumps(
            {
                "passed": True,
                "total_attempts": validation["total_attempts"],
                "scientific_comparison_emitted": result is not None,
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
