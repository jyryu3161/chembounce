"""CPU-only, immutable assessment of the separately declared 270-attempt study.

Invoke this file directly with the prepared runner Python. No generation,
checkpoint loading, fragment-DB API, retry or campaign scheduling is provided.
"""

import argparse
from collections import Counter, defaultdict
import copy
import hashlib
import importlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import statistics
import sys
import time
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare_comparator_followup as preparation

read_json, digest, encode = preparation.read_json, preparation.digest, preparation.encode
ACCEPTED = {"success", "no_candidates"}
STATUSES = ACCEPTED | {
    "timeout",
    "memory_limit",
    "launch_error",
    "worker_error",
    "missing_response",
    "invalid_response",
    "budget_exceeded",
    "interrupted",
    "infrastructure_interruption",
}
INFRA = "infrastructure_interruption"
RETENTION_FIELDS = (
    "query_scaffold_subgraph_retained",
    "query_generic_ring_skeleton_subgraph_retained",
    "core_retaining_scaffold_expansion",
    "query_scaffold_not_contained",
)
BOUNDED = (
    "attempt_with_hop",
    "attempt_with_verified_route",
    "attempt_with_valid_output",
    "attempt_with_common_3d",
    "common_3d_output_slot_coverage",
    "route_supplied_output_slot_coverage",
    *(f"attempt_with_{name}" for name in RETENTION_FIELDS),
    *(f"{name}_assessed_output_slot_coverage" for name in RETENTION_FIELDS),
)
UNITS = {
    "pipeline_current": {
        "synthesis product/design 3D screens (refinement separately reported)",
        "returned design scores (internal generator work also in resource artifact)",
    },
    "crem_train_partial_exo": {"CReM API-yielded products"},
    "mol2mol_pretrained": {"official Mol2Mol raw sampled sequences"},
}


def fingerprint(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def key(row):
    return row["split"], row["method_id"], row["query_id"], row["seed"]


def occurrence_key(row):
    return (*key(row), row["rank"], row["smiles"])


def same(actual, expected, message):
    """Strict JSON structure with tolerance only for finite floating arithmetic."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            raise ValueError(message)
        for name in expected:
            same(actual[name], expected[name], message + ":" + name)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise ValueError(message)
        for a, b in zip(actual, expected):
            same(a, b, message)
    elif isinstance(expected, float):
        if (
            isinstance(actual, bool)
            or not isinstance(actual, (int, float))
            or not math.isfinite(actual)
            or not math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12)
        ):
            raise ValueError(message)
    elif type(actual) is not type(expected) or actual != expected:
        raise ValueError(message)


def finite(value, label, *, nullable=False, minimum=0):
    if nullable and value is None:
        return value
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or (minimum is not None and value < minimum)
    ):
        raise ValueError(f"Invalid finite measurement: {label}")
    return value


def file_inventory(root):
    root = Path(root).resolve()
    result = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError("Result artifacts may not be symlinks")
        if path.is_file():
            result[str(path.relative_to(root))] = digest(path)
    return result


def runtime_packages():
    return {
        d.metadata["Name"].lower().replace("_", "-"): d.version for d in importlib.metadata.distributions()
    }


def expected_runtime(study):
    python = study["plan"]["commands"][0]["argv"][0]
    return python, study["plan"]["runtimes"][python]


def new_output(path, study):
    raw = Path(path)
    if raw.exists() or raw.is_symlink():
        raise ValueError("Assessment output must be new; no overwrite or resume")
    decision = read_json(study["plan"]["parent_decision"]["artifact"])
    amendment = read_json(study["plan"]["amendment"]["artifact"])
    result = preparation.validate_output_location(raw, decision, amendment)
    protected = [study["root"], *(Path(c["argv"][-1]) for c in study["plan"]["commands"])]
    if any(result.is_relative_to(Path(p).resolve()) for p in protected):
        raise ValueError("Assessment output cannot modify prepared plans or generation artifacts")
    return result


def publish(output, files, provenance, *, already_created=False):
    if not already_created:
        output.mkdir(parents=True, exist_ok=False)
    for relative, data in files.items():
        path = output / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x") as stream:
            stream.write(encode(data))
    manifest = {
        "schema_version": 1,
        **provenance,
        "files_sha256": file_inventory(output),
        "assessment_code_sha256": digest(__file__),
    }
    with (output / "manifest.json").open("x") as stream:
        stream.write(encode(manifest))
    return digest(output / "manifest.json")


def verify_archive(path, manifest_sha256):
    path = Path(path).resolve()
    if digest(path / "manifest.json") != manifest_sha256:
        raise ValueError("Assessment manifest differs from parent-pinned SHA-256")
    manifest = read_json(path / "manifest.json")
    if manifest.get("assessment_code_sha256") != digest(__file__):
        raise ValueError("Assessment archive uses another assessment implementation")
    actual = file_inventory(path)
    actual.pop("manifest.json")
    if actual != manifest["files_sha256"]:
        raise ValueError("Assessment archive files or membership changed")
    return manifest


def load_study(root, manifest_sha256):
    root = Path(root).resolve()
    preparation.verify_prepared(root, manifest_sha256)
    plan = read_json(root / "execution-plan.json")
    grid = read_json(root / "grid.json")
    expected = set()
    for split in preparation.SPLITS:
        public = read_json(root / f"public-queries/{split}.json")
        if len(public) != preparation.COUNTS[split]:
            raise ValueError("Public query count differs from declared study")
        expected.update(
            (split, method, q["id"], seed)
            for method in preparation.METHODS
            for q in public
            for seed in preparation.SEEDS
        )
    if len(grid) != 270 or len({key(r) for r in grid}) != 270 or {key(r) for r in grid} != expected:
        raise ValueError("Require the exact 270-slot supplementary grid")
    commands = plan["commands"]
    if len(commands) != 9 or {(c["split"], c["method_id"]) for c in commands} != {
        (s, m) for s in preparation.SPLITS for m in preparation.METHODS
    }:
        raise ValueError("Require exactly nine supplementary campaigns")
    analysis = read_json(root / "analysis-plan.json")
    addendum = analysis["parent_scaffold_retention_addendum"]
    if (
        addendum["sha256"] != preparation.RETENTION_ADDENDUM_SHA
        or plan["input_files_sha256"].get(addendum["artifact"]) != addendum["sha256"]
        or digest(addendum["artifact"]) != addendum["sha256"]
    ):
        raise ValueError("Scaffold-retention addendum is not bound to the parent-approved study")
    same(
        analysis["secondary_scaffold_retention"],
        read_json(addendum["artifact"]),
        "Scaffold-retention definitions differ",
    )
    return {
        "root": root,
        "anchor": manifest_sha256,
        "plan": plan,
        "grid": grid,
        "common": read_json(root / "common-3d-plan.json"),
        "analysis": analysis,
    }


def collect_generation(study):
    """No chemistry is imported: first gate every campaign and exact attempt join."""
    records, campaigns = [], []
    for command in study["plan"]["commands"]:
        split, method = command["split"], command["method_id"]
        root = Path(command["argv"][-1]).resolve()
        expected = {key(g): g for g in study["grid"] if g["split"] == split and g["method_id"] == method}
        inventory = file_inventory(root)
        for required in ("runs.json", "manifest.json", "campaign.json"):
            if required not in inventory:
                raise ValueError("Generation campaign is incomplete")
        manifest, completion = read_json(root / "manifest.json"), read_json(root / "campaign.json")
        protocol = Path(command["argv"][command["argv"].index("--protocol") + 1])
        panel = Path(command["argv"][command["argv"].index("--dataset") + 1])
        if manifest["protocol_sha256"] != digest(protocol) or manifest["dataset"]["dataset_sha256"] != digest(
            panel
        ):
            raise ValueError("Generation protocol/panel provenance differs")
        if manifest["code_sha256"] != study["plan"]["source"]["files_sha256"]:
            raise ValueError("Generation did not use the tested 41-file scientific source")
        _, runtime = expected_runtime(study)
        if (
            manifest.get("python") != runtime["version"]
            or any(runtime["packages"].get(k) != v for k, v in manifest.get("versions", {}).items())
            or set(manifest.get("versions", {})) != {"rdkit", "numpy", "scipy", "psutil"}
        ):
            raise ValueError("Generation runtime provenance differs")
        rows = read_json(root / "runs.json")
        if (
            len(rows) != len(expected)
            or completion.get("status") != "completed"
            or completion.get("attempted_runs") != len(expected)
        ):
            raise ValueError("Refuse incomplete campaign comparisons")
        seen = set()
        for row in rows:
            tagged = {**row, "split": split}
            ident = key(tagged)
            if ident not in expected or ident in seen or row["method_id"] != method:
                raise ValueError("Duplicate or unplanned generation slot")
            seen.add(ident)
            slot = expected[ident]
            if row["target_id"] != slot["target_id"] or row.get("status") not in STATUSES:
                raise ValueError("Attempt target/status does not match plan")
            directory = Path(slot["output_dir"]).resolve()
            if not directory.is_relative_to(root):
                raise ValueError("Attempt output escaped its campaign")
            if read_json(directory / "evaluation.json") != row:
                raise ValueError("runs.json differs from per-attempt evaluation.json")
            response_path = directory / "response.json"
            actual_response_hash = digest(response_path) if response_path.is_file() else None
            if row.get("response_sha256") != actual_response_hash:
                raise ValueError("Raw response hash mismatch")
            if row["status"] in ACCEPTED:
                if actual_response_hash is None or read_json(response_path) != row["response"]:
                    raise ValueError("Accepted response differs from preserved raw response")
                candidates = row["response"].get("candidates")
                if (
                    not isinstance(candidates, list)
                    or len(candidates) > 3
                    or bool(candidates) != (row["status"] == "success")
                ):
                    raise ValueError("Accepted response/status/output-cap mismatch")
            elif row["evaluation"].get("candidates"):
                raise ValueError("Failed attempt contains accepted candidates")
            records.append(
                {
                    "slot": slot,
                    "run": tagged,
                    "request": read_json(directory / "request.json"),
                    "protocol": read_json(protocol),
                    "protocol_path": str(protocol.resolve()),
                    "panel": str(panel),
                    "directory": str(directory),
                }
            )
        if seen != set(expected):
            raise ValueError("Incomplete method/query/seed grid")
        campaigns.append({"split": split, "method_id": method, "root": str(root), "files_sha256": inventory})
    if len(records) != 270 or len({key(r["run"]) for r in records}) != 270:
        raise ValueError("All 270 unique attempts are required before assessment")
    return records, campaigns


def freeze_generation(study, output):
    output = new_output(output, study)
    _, campaigns = collect_generation(study)
    snapshot = {"study_manifest_sha256": study["anchor"], "expected_attempts": 270, "campaigns": campaigns}
    # A changing result tree must not produce a completed snapshot.
    for c in campaigns:
        if file_inventory(c["root"]) != c["files_sha256"]:
            raise ValueError("Generation artifacts changed during snapshot")
    return publish(
        output,
        {"generation.json": snapshot},
        {"stage": "generation_snapshot", "complete": True, "study_manifest_sha256": study["anchor"]},
    )


def load_generation(study, archive, anchor):
    manifest = verify_archive(archive, anchor)
    if (
        manifest.get("stage") != "generation_snapshot"
        or not manifest.get("complete")
        or manifest["study_manifest_sha256"] != study["anchor"]
    ):
        raise ValueError("Generation snapshot belongs to another study or is incomplete")
    snapshot = read_json(Path(archive) / "generation.json")
    records, campaigns = collect_generation(study)
    same(campaigns, snapshot["campaigns"], "Generation files changed after parent-pinned snapshot")
    return records


def scientific_api(study):
    """Load only the pinned CPU evaluator implementation, never model adapters."""
    source = study["plan"]["source"]
    root = Path(source["root"]).resolve()
    expected_python = Path(study["plan"]["commands"][0]["argv"][0]).absolute()
    if Path(sys.executable).absolute() != expected_python:
        raise ValueError("Use the exact prepared runner Python for assessment")
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "" or any(
        os.environ.get(k) != "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
    ):
        raise ValueError("Assessment requires CUDA hidden and single-thread CPU environment")
    sys.path.insert(0, str(root))
    for name in ("benchmarks", "benchmarks.research", "chembounce_synthesis"):
        module = sys.modules.get(name)
        if (
            module is not None
            and getattr(module, "__file__", None)
            and not Path(module.__file__).resolve().is_relative_to(root)
        ):
            raise ValueError("Invoke this file directly; another project checkout was already imported")
    protocol = importlib.import_module("benchmarks.research.protocol")
    runner = importlib.import_module("benchmarks.research.runner")
    evaluation = importlib.import_module("benchmarks.research.evaluation")
    controls = importlib.import_module("benchmarks.research.controls")
    identity = importlib.import_module("chembounce_synthesis.identity")
    chem = importlib.import_module("rdkit.Chem")
    api = SimpleNamespace(
        load_dataset=protocol.load_dataset,
        public_request=runner.public_request,
        evidence_verifier=runner.evidence_verifier,
        evaluate=evaluation.evaluate_candidates,
        identity=identity.chemical_identity,
        similarity=controls.similarity_scores,
        stability=evaluation.score_stability,
        interval=evaluation.cluster_interval,
        mol_from_smiles=chem.MolFromSmiles,
    )
    verify_science(study)
    return api


def verify_science(study):
    source = study["plan"]["source"]
    actual = preparation.source_snapshot(source["root"], source["commit"], preparation.Pins())
    if actual != source:
        raise ValueError("Scientific source changed during assessment")
    if runtime_packages() != expected_runtime(study)[1]["packages"]:
        raise ValueError("Assessment distribution inventory differs from prepared CPU evaluator runtime")
    verify_import_origins(source, sys.modules)


def verify_import_origins(source, modules):
    """Check every loaded project module, including top-level helper modules."""
    root = Path(source["root"]).resolve()
    top_level = {Path(p).stem for p in source["files_sha256"] if "/" not in p}
    for name, module in list(modules.items()):
        if name in top_level | {"chembounce_synthesis", "benchmarks.research"} or name.startswith(
            ("chembounce_synthesis.", "benchmarks.research.")
        ):
            origin = getattr(module, "__file__", None)
            if origin:
                path = Path(origin).resolve()
                if not path.is_relative_to(root) or source["files_sha256"].get(
                    str(path.relative_to(root))
                ) != digest(path):
                    raise ValueError(f"Evaluator import outside pinned scientific source: {name}")


def raw_accounting(record, query_exact=None):
    """Preserve native counts; partial/unknown work is never converted to zero."""
    run, directory = record["run"], Path(record["directory"])
    accepted = run["status"] in ACCEPTED
    disk = None
    try:
        disk = read_json(directory / "response.json")
    except (ValueError, OSError):
        pass
    result = {
        "accepted_response": accepted,
        "raw_response_present": (directory / "response.json").is_file(),
        "raw_response_parseable": isinstance(disk, dict),
        "evaluations": None,
        "evaluation_unit": None,
        "raw_products_or_sequences": None,
        "raw_count_scope": "unavailable",
        "selection_counts": None,
        "internal_work_comparable_between_methods": False,
    }
    if isinstance(disk, dict):
        result.update(
            evaluations=disk.get("evaluations"),
            evaluation_unit=disk.get("evaluation_unit"),
            raw_returned_candidates=len(disk["candidates"])
            if isinstance(disk.get("candidates"), list)
            else None,
        )
    if accepted:
        count = result["evaluations"]
        cap = 12 if run["method_id"] == "pipeline_current" else 256
        if (
            type(count) is not int
            or not 0 <= count <= cap
            or result["evaluation_unit"] not in UNITS[run["method_id"]]
        ):
            raise ValueError("Accepted native-work accounting differs from declared method/cap")
    elif type(result["evaluations"]) is not int or result["evaluations"] < 0:
        result["invalid_native_count_declaration"] = result["evaluations"]
        result["evaluations"] = None
    method = run["method_id"]
    if method == "crem_train_partial_exo":
        path = directory / "response.raw.jsonl"
        raw = []
        malformed = 0
        if path.is_file():
            for line in path.read_text().splitlines():
                try:
                    row = json.loads(line)
                    if row.get("index") != len(raw) or not isinstance(row.get("raw_smiles"), str):
                        raise ValueError("Malformed CReM raw row")
                    raw.append(row)
                except (ValueError, TypeError, AttributeError):
                    malformed += 1
            result.update(
                raw_products_or_sequences=len(raw),
                raw_count_scope="complete" if accepted else "preserved_completed_rows_lower_bound",
                malformed_raw_lines=malformed,
            )
        if accepted:
            metadata = disk.get("metadata", {})
            if malformed or len(raw) != disk["evaluations"] or metadata.get("raw_products") != len(raw):
                raise ValueError("CReM raw-product accounting mismatch")
            decisions = metadata.get("decisions", [])
            if len(decisions) != len(raw) or any(
                d.get("index") != i or d.get("raw_smiles") != raw[i]["raw_smiles"]
                for i, d in enumerate(decisions)
            ):
                raise ValueError("CReM raw/selection row correspondence differs")
            counts = dict(Counter(d.get("reason", d["status"]) for d in decisions))
            if counts != metadata.get("postfilter_counts", {}):
                raise ValueError("CReM selection counts differ")
            result["selection_counts"] = counts
            selected = sorted(
                (d for d in decisions if d["status"] == "selected"),
                key=lambda d: (-d["count_tanimoto"], d["exact_smiles"], d["index"]),
            )
            if [d["exact_smiles"] for d in selected] != [c["smiles"] for c in disk["candidates"]]:
                raise ValueError("CReM selected rows differ from returned candidates")
    elif method == "mol2mol_pretrained":
        path, selected_path = directory / "mol2mol-raw.json", directory / "mol2mol-selection.json"
        raw = None
        if path.is_file():
            try:
                raw = read_json(path)
                result.update(
                    raw_products_or_sequences=len(raw["rows"]),
                    raw_count_scope="complete" if accepted else "preserved_completed_batches_lower_bound",
                )
            except (ValueError, KeyError, TypeError):
                if accepted:
                    raise ValueError("Accepted Mol2Mol raw artifact is invalid")
        if accepted:
            metadata = disk.get("metadata", {})
            if (
                raw is None
                or len(raw["rows"]) != disk["evaluations"]
                or metadata.get("raw_sha256") != digest(path)
                or metadata.get("selection_sha256") != digest(selected_path)
            ):
                raise ValueError("Mol2Mol raw/selection provenance mismatch")
            if raw.get("seed") != run["seed"] or (
                query_exact is not None and raw.get("input_smiles") != query_exact
            ):
                raise ValueError("Mol2Mol raw samples belong to another query/seed")
            decisions = read_json(selected_path)
            if len(decisions) != len(raw["rows"]) or any(
                d.get("raw_index") != i or any(d.get(k) != v for k, v in raw["rows"][i].items())
                for i, d in enumerate(decisions)
            ):
                raise ValueError("Mol2Mol raw/selection row correspondence differs")
            counts = dict(Counter(d["selection_status"] for d in decisions))
            if counts != metadata.get("selection_counts"):
                raise ValueError("Mol2Mol selection counts differ")
            result["selection_counts"] = counts
            selected = sorted(
                (d for d in decisions if d["selection_status"] == "selected"),
                key=lambda d: (d["native_nll"], d["canonical_smiles"], d["raw_index"]),
            )
            if [(d["raw_index"], d["canonical_smiles"], d["native_nll"]) for d in selected] != [
                (c["raw_index"], c["smiles"], c["native_nll"]) for c in disk["candidates"]
            ]:
                raise ValueError("Mol2Mol selected rows differ from returned candidates")
    return result


def scaffold_retention(query_smiles, candidate_smiles, api):
    """Secondary, asymmetric query-within-candidate graph annotations only."""
    result = {name: None for name in RETENTION_FIELDS}
    result.update(
        query_scaffold_heavy_atoms=None,
        candidate_scaffold_heavy_atoms=None,
        atom_bond_status="unassessed",
        generic_status="unassessed",
        reason=None,
    )
    try:
        query, candidate = api.identity(query_smiles), api.identity(candidate_smiles)
        q = api.mol_from_smiles(query.scaffold_family_smiles) if query.scaffold_family_smiles else None
        c = (
            api.mol_from_smiles(candidate.scaffold_family_smiles)
            if candidate.scaffold_family_smiles
            else None
        )
        result["query_scaffold_heavy_atoms"] = q.GetNumHeavyAtoms() if q is not None else 0
        result["candidate_scaffold_heavy_atoms"] = c.GetNumHeavyAtoms() if c is not None else 0
        if not query.has_scaffold or not candidate.has_scaffold or q is None or c is None:
            result["reason"] = "Either input lacks a valid normalized ring scaffold"
            return result
        retained = bool(c.HasSubstructMatch(q, useChirality=False))
        hop = (
            candidate.family_smiles != query.family_smiles
            and candidate.scaffold_family_smiles != query.scaffold_family_smiles
        )
        result.update(
            query_scaffold_subgraph_retained=retained,
            atom_bond_status="assessed",
            core_retaining_scaffold_expansion=hop
            and retained
            and c.GetNumHeavyAtoms() > q.GetNumHeavyAtoms(),
            query_scaffold_not_contained=hop and not retained,
        )
        qg = api.mol_from_smiles(query.ring_skeleton_smiles) if query.ring_skeleton_smiles else None
        cg = api.mol_from_smiles(candidate.ring_skeleton_smiles) if candidate.ring_skeleton_smiles else None
        if qg is not None and cg is not None and qg.GetNumHeavyAtoms() and cg.GetNumHeavyAtoms():
            result.update(
                query_generic_ring_skeleton_subgraph_retained=bool(
                    cg.HasSubstructMatch(qg, useChirality=False)
                ),
                generic_status="assessed",
            )
        else:
            result["reason"] = "Generic ring skeleton unavailable"
    except (ValueError, RuntimeError) as exc:
        result = {name: None for name in RETENTION_FIELDS} | {
            "query_scaffold_heavy_atoms": None,
            "candidate_scaffold_heavy_atoms": None,
            "atom_bond_status": "unassessed",
            "generic_status": "unassessed",
            "reason": f"{type(exc).__name__}: {exc}",
        }
    return result


def recompute(study, records, api):
    datasets, verifiers, rebuilt, occurrences = {}, {}, [], []
    for record in records:
        run, slot, protocol = record["run"], record["slot"], record["protocol"]
        split = run["split"]
        if split not in datasets:
            dataset, _ = api.load_dataset(record["panel"])
            datasets[split] = dataset, {q["id"]: q for q in dataset["queries"]}
        dataset, queries = datasets[split]
        query = queries[run["query_id"]]
        method = copy.deepcopy(protocol["methods"][0])
        # Match runner.load_protocol: resolve the protocol itself first, then
        # resolve file/directory assets against its physical parent directory.
        protocol_base = Path(record["protocol_path"]).resolve().parent
        method["resolved_assets"] = {
            k: str((protocol_base / (v["directory"] if "directory" in v else v["artifact"])).resolve())
            for k, v in method["assets"].items()
        }
        expected_request = api.public_request(
            dataset, query, method, run["seed"], protocol["budget"], slot["output_dir"]
        )
        request = copy.deepcopy(record["request"])
        finite(request.pop("deadline_monotonic"), "deadline")
        if request.pop("finalization_margin_seconds") != 2:
            raise ValueError("Generation finalization margin differs")
        same(request, expected_request, "Method request/private-field/settings provenance differs")
        verifier_key = fingerprint(
            {k: v for k, v in method["resolved_assets"].items() if k in {"catalog", "rules"}}
        )
        if verifier_key not in verifiers:
            verifiers[verifier_key] = api.evidence_verifier(method["resolved_assets"])
        evaluated = api.evaluate(
            dataset, query, run["response"].get("candidates", []), run["status"], verifiers[verifier_key]
        )
        same(run["evaluation"], evaluated, "Recomputed chemical metrics/occurrences differ")
        public_query = expected_request["query"]
        query_exact = api.identity(public_query["smiles"]).exact_smiles
        raw = raw_accounting(record, query_exact)
        resources = run["resources"]
        for field in ("wall_seconds", "sampled_peak_rss_mb", "sampled_cpu_seconds_lower_bound"):
            finite(resources.get(field), field, nullable=run["status"] == INFRA)
        if run["status"] != INFRA and resources.get("includes_import_and_model_loading") is not True:
            raise ValueError("Missing cold-start resource scope")
        item = {
            "split": split,
            "method_id": run["method_id"],
            "query_id": run["query_id"],
            "seed": run["seed"],
            "target_id": run["target_id"],
            "status": run["status"],
            "evaluation": evaluated,
            "raw_accounting": raw,
            "resources": resources,
            "response_sha256": run["response_sha256"],
        }
        rebuilt.append(item)
        accepted_raw = run["response"].get("candidates", []) if run["status"] in ACCEPTED else []
        for candidate in evaluated["candidates"]:
            source_candidate = accepted_raw[candidate["rank"] - 1]
            if api.identity(source_candidate["smiles"]).exact_smiles != candidate["smiles"]:
                raise ValueError("Raw candidate rank/identity correspondence differs")
            occurrences.append(
                {
                    **{k: item[k] for k in ("split", "method_id", "query_id", "target_id", "seed")},
                    "rank": candidate["rank"],
                    "smiles": candidate["smiles"],
                    "query_exact": query_exact,
                    "candidate": candidate,
                    "route_supplied": bool(source_candidate.get("route")),
                    "serialized_route_sha256": fingerprint(source_candidate["route"])
                    if source_candidate.get("route")
                    else None,
                    "scaffold_retention": scaffold_retention(query_exact, candidate["smiles"], api),
                }
            )
    return rebuilt, occurrences


def cache_payload(study, occurrence):
    return {
        "query_exact": occurrence["query_exact"],
        "candidate_exact": occurrence["smiles"],
        "source": study["plan"]["source"],
        "conformers": study["common"]["conformers"],
        "seeds": study["common"]["seeds"],
        "identity_policy": "protonation-tautomer-family-v1",
        "rdkit": "2025.09.6",
    }


def score_pair(payload, api):
    started = time.monotonic()
    try:
        scores, trials = api.similarity(
            payload["query_exact"], payload["candidate_exact"], payload["seeds"], payload["conformers"]
        )
        error = None
    except (ValueError, RuntimeError) as exc:
        scores, error = {}, f"{type(exc).__name__}: {exc}"
        trials = [
            {"seed": s, "status": "unassessed_due_to_pair_error", "reason": error} for s in payload["seeds"]
        ]
    return {
        "payload": payload,
        "cache_key": fingerprint(payload),
        "scores": scores,
        "trials": trials,
        "all_3d_seeds_succeeded": all(t["status"] == "success" for t in trials),
        "pair_error": error,
        "wall_seconds": time.monotonic() - started,
    }


def validate_score(scored, payload, candidate, api):
    same(scored["payload"], payload, "Common-score exact-pair/source/config provenance differs")
    if scored["cache_key"] != fingerprint(payload):
        raise ValueError("Common-score cache key differs")
    trials = scored["trials"]
    if len(trials) != 3 or [t["seed"] for t in trials] != preparation.SEEDS:
        raise ValueError("Common scoring must retain all three fixed seed outcomes")
    if any(t["status"] not in {"success", "failed", "unassessed_due_to_pair_error"} for t in trials):
        raise ValueError("Unrecognized scoring outcome")
    complete = all(t["status"] == "success" for t in trials)
    if scored["all_3d_seeds_succeeded"] is not complete:
        raise ValueError("Common scoring success flag differs from seed coverage")
    if scored.get("pair_error") is not None:
        if scored["scores"] or any(t["status"] != "unassessed_due_to_pair_error" for t in trials):
            raise ValueError("Pair-level evaluator error cannot contain assessed successful scores")
    elif any(t["status"] == "unassessed_due_to_pair_error" for t in trials):
        raise ValueError("Unassessed scoring seeds require a recorded pair-level error")
    for t in trials:
        if t["status"] == "success":
            for name in ("hit_shape", "hit_pharmacophore", "score"):
                finite(t["metrics"][name], name, minimum=None)
            if not all(0 <= t["metrics"][name] <= 1 for name in ("hit_shape", "hit_pharmacophore")):
                raise ValueError("Shape/pharmacophore scores must lie in [0,1]")
        elif not t.get("reason"):
            raise ValueError("Scoring failure is missing its reason")
    if scored.get("pair_error") is None:
        same(
            scored["scores"]["bit_tanimoto"],
            candidate["hit_bit_tanimoto"],
            "Common bit similarity differs from direct recomputation",
        )
        same(
            scored["scores"]["count_tanimoto"],
            candidate["hit_count_tanimoto"],
            "Common count similarity differs from direct recomputation",
        )
    for name in ("hit_shape", "hit_pharmacophore", "score"):
        if complete:
            values = [t["metrics"][name] for t in trials]
            same(scored["scores"].get(name), float(statistics.median(values)), "Common seed median differs")
            same(
                scored["scores"].get(name + "_stability"),
                api.stability(values),
                "Common seed stability differs",
            )
        elif name in scored["scores"] or name + "_stability" in scored["scores"]:
            raise ValueError("Incomplete 3D pair cannot publish successful conditional aggregate")


def score_generation(study, snapshot, snapshot_sha, output, api):
    output = new_output(output, study)
    records = load_generation(study, snapshot, snapshot_sha)
    rebuilt, occurrences = recompute(study, records, api)
    output.mkdir(parents=True, exist_ok=False)
    (output / "cache").mkdir()
    with (output / "context.json").open("x") as stream:
        stream.write(
            encode(
                {
                    "study_manifest_sha256": study["anchor"],
                    "generation_manifest_sha256": snapshot_sha,
                    "common_plan_sha256": digest(study["root"] / "common-3d-plan.json"),
                    "assessment_code_sha256": digest(__file__),
                    "complete": False,
                    "completion_requires_final_manifest": True,
                }
            )
        )
    files, cache, joined = {"recomputed-attempts.json": rebuilt}, {}, []
    for occurrence in occurrences:
        payload = cache_payload(study, occurrence)
        cache_key = fingerprint(payload)
        if cache_key not in cache:
            cache[cache_key] = score_pair(payload, api)
            validate_score(cache[cache_key], payload, occurrence["candidate"], api)
            # Completed pairs survive an interrupted CPU assessment. They are
            # diagnostic until a complete manifest exists; there is no auto-resume.
            with (output / "cache" / f"{cache_key}.json").open("x") as stream:
                stream.write(encode(cache[cache_key]))
        validate_score(cache[cache_key], payload, occurrence["candidate"], api)
        joined.append(
            {
                k: occurrence[k]
                for k in (
                    "split",
                    "method_id",
                    "query_id",
                    "target_id",
                    "seed",
                    "rank",
                    "smiles",
                    "query_exact",
                )
            }
            | {"cache_key": cache_key}
        )
    files["occurrences.json"] = joined
    files["coverage.json"] = {
        "attempts": len(rebuilt),
        "valid_candidate_occurrences": len(joined),
        "unique_exact_pairs": len(cache),
        "successful_3d_occurrences": sum(cache[o["cache_key"]]["all_3d_seeds_succeeded"] for o in joined),
        "interpretation": "All valid returned candidate occurrences retained; shared exact pairs reuse computation, not independent observations",
    }
    load_generation(study, snapshot, snapshot_sha)
    verify_science(study)
    return publish(
        output,
        files,
        {
            "stage": "common_3d",
            "complete": True,
            "study_manifest_sha256": study["anchor"],
            "generation_manifest_sha256": snapshot_sha,
            "common_plan_sha256": digest(study["root"] / "common-3d-plan.json"),
            "source": study["plan"]["source"],
            "execution": "CPU only, sequential, generation budgets excluded",
            "runtime_python": sys.executable,
            "runtime_packages": runtime_packages(),
        },
        already_created=True,
    )


def load_scores(study, directory, anchor, generation_anchor, occurrences, api):
    manifest = verify_archive(directory, anchor)
    if (
        manifest.get("stage") != "common_3d"
        or manifest.get("complete") is not True
        or manifest["study_manifest_sha256"] != study["anchor"]
        or manifest["generation_manifest_sha256"] != generation_anchor
    ):
        raise ValueError("Common scoring archive does not bind this completed generation snapshot")
    if (
        manifest["common_plan_sha256"] != digest(study["root"] / "common-3d-plan.json")
        or manifest["source"] != study["plan"]["source"]
    ):
        raise ValueError("Common scoring source/config provenance differs")
    python, runtime = expected_runtime(study)
    if manifest.get("runtime_python") != python or manifest.get("runtime_packages") != runtime["packages"]:
        raise ValueError("Common scoring runtime provenance differs")
    rows = read_json(Path(directory) / "occurrences.json")
    expected = {occurrence_key(o): o for o in occurrences}
    if (
        len(expected) != len(occurrences)
        or len(rows) != len(expected)
        or len({occurrence_key(o) for o in rows}) != len(rows)
        or {occurrence_key(o) for o in rows} != set(expected)
    ):
        raise ValueError("Missing, duplicate, or extraneous common-scoring candidate occurrence")
    result, cache = {}, {}
    for row in rows:
        occurrence = expected[occurrence_key(row)]
        if row["query_exact"] != occurrence["query_exact"] or row["target_id"] != occurrence["target_id"]:
            raise ValueError("Common-score query/target correspondence differs")
        payload = cache_payload(study, occurrence)
        cache_key = fingerprint(payload)
        if row["cache_key"] != cache_key:
            raise ValueError("Occurrence uses wrong exact-pair cache key")
        if cache_key not in cache:
            cache[cache_key] = read_json(Path(directory) / "cache" / f"{cache_key}.json")
        validate_score(cache[cache_key], payload, occurrence["candidate"], api)
        result[occurrence_key(row)] = cache[cache_key]
    actual_caches = {p.stem for p in (Path(directory) / "cache").glob("*.json")}
    if actual_caches != set(cache):
        raise ValueError("Scoring cache contains omitted or extraneous exact pairs")
    coverage = read_json(Path(directory) / "coverage.json")
    expected_coverage = {
        "attempts": 270,
        "valid_candidate_occurrences": len(rows),
        "unique_exact_pairs": len(cache),
        "successful_3d_occurrences": sum(result[occurrence_key(r)]["all_3d_seeds_succeeded"] for r in rows),
    }
    if any(coverage.get(k) != v for k, v in expected_coverage.items()):
        raise ValueError("Common scoring coverage counts differ from occurrence/cache joins")
    return result


def mean(values):
    observed = [v for v in values if v is not None]
    return statistics.mean(observed) if observed else None


def attempt_metrics(attempts, occurrences, scores):
    grouped = defaultdict(list)
    for o in occurrences:
        grouped[key(o)].append(o)
    result = []
    for attempt in attempts:
        rows = grouped[key(attempt)]
        successful = [o for o in rows if scores[occurrence_key(o)]["all_3d_seeds_succeeded"]]
        route_supplied = sum(o["route_supplied"] for o in rows)
        metrics = dict(attempt["evaluation"]["metrics"])
        metrics.update(
            attempt_with_hop=int(metrics["hop_families"] > 0),
            attempt_with_valid_output=int(bool(rows)),
            attempt_with_verified_route=int(metrics["verified_rule_routes"] > 0),
            attempt_with_common_3d=int(bool(successful)),
            common_3d_successful_candidates=len(successful),
            common_3d_failed_candidates=len(rows) - len(successful),
            common_3d_output_slot_coverage=len(successful) / 3,
            route_supplied_occurrences=route_supplied,
            no_route_supplied_occurrences=len(rows) - route_supplied,
            route_supplied_output_slot_coverage=route_supplied / 3,
            route_replay_verified_occurrences=sum(
                o["candidate"]["route_evidence"]["level"] == "catalog_and_rule_verified" for o in rows
            ),
            route_rejected_occurrences=sum(
                o["candidate"]["route_evidence"]["level"] == "route_rejected" for o in rows
            ),
            wall_seconds=attempt["resources"].get("wall_seconds"),
            hit_bit_tanimoto_mean_conditional=mean(o["candidate"]["hit_bit_tanimoto"] for o in rows),
            hit_count_tanimoto_mean_conditional=mean(o["candidate"]["hit_count_tanimoto"] for o in rows),
        )
        for name in ("hit_shape", "hit_pharmacophore", "score"):
            metrics[f"common_{name}_mean_success_conditional"] = mean(
                scores[occurrence_key(o)]["scores"][name] for o in successful
            )
            metrics[f"common_{name}_range_mean_success_conditional"] = mean(
                scores[occurrence_key(o)]["scores"][name + "_stability"]["range"] for o in successful
            )
        for name in RETENTION_FIELDS:
            qualifying = [o for o in rows if o["scaffold_retention"][name] is True]
            assessed = sum(o["scaffold_retention"][name] is not None for o in rows)
            metrics.update(
                {
                    name + "_occurrences": len(qualifying),
                    name + "_unique_families": len({o["candidate"]["family"] for o in qualifying}),
                    name + "_assessed_occurrences": assessed,
                    name + "_unassessed_occurrences": len(rows) - assessed,
                    name + "_assessable_fraction_conditional": assessed / len(rows) if rows else None,
                    name + "_assessed_output_slot_coverage": assessed / 3,
                    "attempt_with_" + name: int(bool(qualifying)),
                }
            )
        declarations = [o["candidate"].get("declared_external_synthesis_evidence") or {} for o in rows]
        steps = [d.get("steps", []) for d in declarations]
        metrics["declared_forward_assessed_occurrences"] = sum(
            bool(s)
            and all(t.get("forward_prediction", {}).get("status") in {"matched", "mismatched"} for t in s)
            for s in steps
        )
        metrics["declared_forward_all_steps_matched_occurrences"] = sum(
            bool(s) and all(t.get("forward_prediction", {}).get("status") == "matched" for t in s)
            for s in steps
        )
        result.append(
            {
                **{k: attempt[k] for k in ("split", "method_id", "query_id", "target_id", "seed", "status")},
                "metrics": metrics,
            }
        )
    return result


def hierarchical(values, slots, interval):
    """Missing conditional scores remain missing; targets always receive equal weight."""
    queries = sorted({s["query_id"] for s in slots})
    targets = sorted({s["target_id"] for s in slots})
    by_query = {
        q: mean(values.get((s["query_id"], s["seed"])) for s in slots if s["query_id"] == q) for q in queries
    }
    per_target = {
        t: mean(
            by_query[q] for q in queries if any(s["query_id"] == q and s["target_id"] == t for s in slots)
        )
        for t in targets
    }
    observed = {
        q: sum(values.get((s["query_id"], s["seed"])) is not None for s in slots if s["query_id"] == q)
        for q in queries
    }
    estimate = interval([v for v in per_target.values() if v is not None], seed=20261016, resamples=10000)
    return {
        "per_query": by_query,
        "per_target": per_target,
        "target_macro": estimate,
        "availability": {
            "planned_slots": len(slots),
            "observed_slots": sum(observed.values()),
            "missing_slots": len(slots) - sum(observed.values()),
            "observed_seeds_per_query": observed,
            "planned_queries": len(queries),
            "observed_queries": sum(v is not None for v in by_query.values()),
            "planned_targets": len(targets),
            "observed_targets": sum(v is not None for v in per_target.values()),
        },
    }


def comparisons(rows, grid, interval):
    output = {}
    for split in preparation.SPLITS:
        planned = [g for g in grid if g["split"] == split and g["method_id"] == preparation.METHODS[0]]
        lookup = {
            m: {(r["query_id"], r["seed"]): r for r in rows if r["split"] == split and r["method_id"] == m}
            for m in preparation.METHODS
        }
        required = {(s["query_id"], s["seed"]) for s in planned}
        if any(set(d) != required for d in lookup.values()):
            raise ValueError("Analysis is missing matched method/query/seed slots")
        names = sorted(next(iter(lookup[preparation.METHODS[0]].values()))["metrics"])
        methods, paired = {}, {}
        for method, local in lookup.items():
            if any(set(r["metrics"]) != set(names) for r in local.values()):
                raise ValueError("Metric schema differs between attempts")
            summary = {
                "attempts": len(local),
                "status_counts": dict(Counter(r["status"] for r in local.values())),
                "metrics": {},
                "complete_case_excluding_infrastructure_only": {},
                "infrastructure_counterfactual_bounds": {},
            }
            for metric in names:
                values = {k: r["metrics"][metric] for k, r in local.items()}
                summary["metrics"][metric] = hierarchical(values, planned, interval)
                complete = {k: r["metrics"][metric] for k, r in local.items() if r["status"] != INFRA}
                summary["complete_case_excluding_infrastructure_only"][metric] = hierarchical(
                    complete, planned, interval
                )
                if metric in BOUNDED:
                    low = {k: 0.0 if r["status"] == INFRA else r["metrics"][metric] for k, r in local.items()}
                    high = {
                        k: 1.0 if r["status"] == INFRA else r["metrics"][metric] for k, r in local.items()
                    }
                    summary["infrastructure_counterfactual_bounds"][metric] = [
                        hierarchical(v, planned, interval)["target_macro"]["mean"] for v in (low, high)
                    ]
            methods[method] = summary
        primary = lookup[preparation.METHODS[0]]
        for method in preparation.METHODS[1:]:
            other = lookup[method]
            paired[method] = {
                "metrics": {},
                "complete_case_excluding_either_infrastructure": {},
                "infrastructure_counterfactual_bounds": {},
            }
            for metric in names:
                paired_values, complete, low, high = {}, {}, {}, {}
                for k in required:
                    a, b = primary[k], other[k]
                    av, bv = a["metrics"][metric], b["metrics"][metric]
                    if av is not None and bv is not None:
                        paired_values[k] = av - bv
                        if a["status"] != INFRA and b["status"] != INFRA:
                            complete[k] = av - bv
                    if metric in BOUNDED:
                        alo, ahi = (0.0, 1.0) if a["status"] == INFRA else (av, av)
                        blo, bhi = (0.0, 1.0) if b["status"] == INFRA else (bv, bv)
                        low[k], high[k] = alo - bhi, ahi - blo
                paired[method]["metrics"][metric] = hierarchical(paired_values, planned, interval)
                paired[method]["complete_case_excluding_either_infrastructure"][metric] = hierarchical(
                    complete, planned, interval
                )
                if metric in BOUNDED:
                    paired[method]["infrastructure_counterfactual_bounds"][metric] = [
                        hierarchical(v, planned, interval)["target_macro"]["mean"] for v in (low, high)
                    ]
        output[split] = {
            "methods": methods,
            "paired_pipeline_current_minus_comparator": paired,
            "counterfactual_bounds_are_confidence_intervals": False,
            "complete_case_is_selection_biased": True,
        }
    return output


def analyze_generation(study, snapshot, snapshot_sha, scores_dir, scores_sha, output, api):
    output = new_output(output, study)
    records = load_generation(study, snapshot, snapshot_sha)
    rebuilt, occurrences = recompute(study, records, api)
    scores = load_scores(study, scores_dir, scores_sha, snapshot_sha, occurrences, api)
    same(
        read_json(Path(scores_dir) / "recomputed-attempts.json"),
        rebuilt,
        "Common-score chemistry assessment differs from fresh recomputation",
    )
    rows = attempt_metrics(rebuilt, occurrences, scores)
    result = comparisons(rows, study["grid"], api.interval)
    unique = {
        split: {
            method: {
                "families": len(
                    {
                        o["candidate"]["family"]
                        for o in occurrences
                        if o["split"] == split and o["method_id"] == method
                    }
                ),
                "exact_query_candidate_pairs": len(
                    {
                        (o["query_exact"], o["smiles"])
                        for o in occurrences
                        if o["split"] == split and o["method_id"] == method
                    }
                ),
                "serialized_routes": len(
                    {
                        o["serialized_route_sha256"]
                        for o in occurrences
                        if o["split"] == split and o["method_id"] == method and o["route_supplied"]
                    }
                ),
            }
            for method in preparation.METHODS
        }
        for split in preparation.SPLITS
    }
    accounting = {}
    retention_summary = {}
    for split in preparation.SPLITS:
        accounting[split] = {}
        retention_summary[split] = {}
        for method in preparation.METHODS:
            selected = [a for a in rebuilt if a["split"] == split and a["method_id"] == method]
            native = [a["raw_accounting"] for a in selected]
            accounting[split][method] = {
                "attempts": len(selected),
                "accepted_attempts": sum(a["status"] in ACCEPTED for a in selected),
                "native_unit_counts": dict(Counter(a["evaluation_unit"] or "unavailable" for a in native)),
                "native_evaluations_observed": sum(type(a["evaluations"]) is int for a in native),
                "raw_count_scope_counts": dict(Counter(a["raw_count_scope"] for a in native)),
                "raw_response_unavailable_or_unparseable": sum(
                    not a["raw_response_parseable"] for a in native
                ),
                "resource_scope": "Runner sampled process-tree RSS/CPU; GPU memory not uniformly measured. Shared host, not a controlled hardware speed benchmark",
                "counts_are_comparable_between_methods": False,
            }
            local = [o for o in occurrences if o["split"] == split and o["method_id"] == method]
            retention_summary[split][method] = {
                "attempts": len(selected),
                "valid_candidate_occurrences": len(local),
                "annotations": {
                    name: {
                        "qualifying_occurrences": sum(o["scaffold_retention"][name] is True for o in local),
                        "unique_qualifying_families": len(
                            {o["candidate"]["family"] for o in local if o["scaffold_retention"][name] is True}
                        ),
                        "assessed_occurrences": sum(o["scaffold_retention"][name] is not None for o in local),
                        "unassessed_occurrences": sum(o["scaffold_retention"][name] is None for o in local),
                        "attempts_with_qualifying_output": len(
                            {key(o) for o in local if o["scaffold_retention"][name] is True}
                        ),
                    }
                    for name in RETENTION_FIELDS
                },
            }
    load_generation(study, snapshot, snapshot_sha)
    verify_archive(scores_dir, scores_sha)
    verify_science(study)
    return publish(
        output,
        {
            "comparison.json": result,
            "attempts.json": rows,
            "recomputed-attempts.json": rebuilt,
            "unique-identities.json": unique,
            "raw-accounting-summary.json": accounting,
            "scaffold-retention-summary.json": retention_summary,
            "scaffold-retention-occurrences.json": [
                {
                    k: o[k]
                    for k in (
                        "split",
                        "method_id",
                        "query_id",
                        "target_id",
                        "seed",
                        "rank",
                        "query_exact",
                        "smiles",
                        "scaffold_retention",
                    )
                }
                | {"family": o["candidate"]["family"]}
                for o in occurrences
            ],
            "analysis-plan.json": study["analysis"],
        },
        {
            "stage": "final_comparison",
            "complete_requires_parent_review": True,
            "study_manifest_sha256": study["anchor"],
            "generation_manifest_sha256": snapshot_sha,
            "common_scores_manifest_sha256": scores_sha,
            "expected_attempts": 270,
            "interpretation": "Engineering delivery and retrospective structural proxies; no activity/synthesis probability claim. No-route is unassessed coverage. Declared forward evidence is not fresh independent inference and shares pipeline selection bias",
        },
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "score", "analyze"))
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--prepared-sha256", required=True)
    parser.add_argument("--generation", type=Path)
    parser.add_argument("--generation-sha256")
    parser.add_argument("--scores", type=Path)
    parser.add_argument("--scores-sha256")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    study = load_study(args.prepared, args.prepared_sha256)
    if args.action == "freeze":
        anchor = freeze_generation(study, args.output)
    else:
        if not args.generation or not args.generation_sha256:
            parser.error("A parent-pinned completed generation snapshot is required")
        api = scientific_api(study)
        if args.action == "score":
            anchor = score_generation(study, args.generation, args.generation_sha256, args.output, api)
        else:
            if not args.scores or not args.scores_sha256:
                parser.error("A parent-pinned completed common-scoring archive is required")
            anchor = analyze_generation(
                study,
                args.generation,
                args.generation_sha256,
                args.scores,
                args.scores_sha256,
                args.output,
                api,
            )
    print(encode({"manifest_sha256": anchor, "generation_executed": False}), end="")


if __name__ == "__main__":
    main()
