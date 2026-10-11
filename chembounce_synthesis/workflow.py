"""End-to-end synthesis-aware scaffold hopping and hit expansion."""

import csv
import importlib.metadata
import json
import logging
import math
import subprocess
import time
from dataclasses import asdict
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold

from .catalog import Catalog
from .common import FP, canonical, digest, molecule, positive, write_json
from .conformers import ConformerSettings, ensemble, reference_from_sdf, require_alignment, score_candidate
from .reactions import Rules
from .routes import SearchLimits, expand_route, resolve_proposal, verify_route
from .synllama import SavedProposals, SynLlama

LOGGER = logging.getLogger(__name__)
COLUMNS = [
    "rank",
    "smiles",
    "design",
    "score",
    "pareto_front",
    "hit_shape",
    "hit_pharmacophore",
    "design_shape",
    "design_pharmacophore",
    "hit_morgan",
    "relative_energy",
    "reaction_steps",
    "route_id",
    "alternative_route_ids",
]


def scaffold(smiles):
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(molecule(smiles)), isomericSmiles=True)


def check_constraints(product, hit, design, mode="hop", core_smiles=None, preserve_design_scaffold=False):
    mol = molecule(product)
    if canonical(product) == canonical(hit) and mode == "hop":
        raise ValueError("Reconstruction returned the original hit")
    if mode == "hop":
        if scaffold(product) == scaffold(hit):
            raise ValueError("Reconstructed product did not change the original Murcko scaffold")
    elif mode == "expand":
        if canonical(product) == canonical(design):
            raise ValueError("Hit expansion returned the unchanged design")
        if scaffold(product) != scaffold(design):
            raise ValueError("Hit expansion changed the selected design's Murcko scaffold")
    else:
        raise ValueError("mode must be hop or expand")
    if preserve_design_scaffold and scaffold(product) != scaffold(design):
        raise ValueError("Reconstruction did not preserve the selected design scaffold")
    if core_smiles and not mol.HasSubstructMatch(molecule(core_smiles), useChirality=True):
        raise ValueError("Reconstructed product lost the required core/anchor substructure")


def _config(config):
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise ValueError("Synthesis configuration requires schema_version=1")
    permitted = {
        "schema_version",
        "catalog",
        "rules",
        "provider",
        "search",
        "conformers",
        "mode",
        "core_smiles",
        "preserve_design_scaffold",
        "max_designs",
        "max_scored_products",
        "max_results",
        "expansion",
        "reference_sdf",
        "minimum_scores",
        "property_thresholds",
    }
    if set(config) - permitted:
        raise ValueError(f"Unknown synthesis configuration keys: {sorted(set(config) - permitted)}")
    for required in ("catalog", "rules", "provider"):
        if required not in config:
            raise ValueError(f"Missing synthesis configuration: {required}")
    for key in ("provider", "search", "conformers", "expansion", "minimum_scores", "property_thresholds"):
        if key in config and not isinstance(config[key], dict):
            raise ValueError(f"{key} must be a JSON object")
    if "preserve_design_scaffold" in config and not isinstance(config["preserve_design_scaffold"], bool):
        raise ValueError("preserve_design_scaffold must be true or false")
    if config.get("core_smiles") is not None:
        canonical(config["core_smiles"])
    for key, default in (("max_designs", 20), ("max_scored_products", 256), ("max_results", 50)):
        positive(config.get(key, default), key)
    if config.get("mode", "hop") not in {"hop", "expand"}:
        raise ValueError("mode must be hop or expand")
    minimum = config.get("minimum_scores", {})
    for key, value in minimum.items():
        if (
            key not in {"hit_shape", "hit_pharmacophore", "design_shape", "design_pharmacophore"}
            or not isinstance(value, (int, float))
            or isinstance(value, bool)
            or not math.isfinite(value)
            or not 0 <= value <= 1
        ):
            raise ValueError("minimum_scores must contain known score names with values in [0, 1]")
    exp = config.get("expansion", {})
    if set(exp) - {"rounds", "max_routes", "locked_building_block_ids"}:
        raise ValueError("Unknown expansion configuration key")
    if (
        isinstance(exp.get("rounds", 0), bool)
        or not isinstance(exp.get("rounds", 0), int)
        or exp.get("rounds", 0) < 0
    ):
        raise ValueError("expansion rounds must be a nonnegative integer")
    positive(exp.get("max_routes", 64), "expansion.max_routes")
    locks = exp.get("locked_building_block_ids", [])
    if not isinstance(locks, list) or not all(isinstance(value, str) and value for value in locks):
        raise ValueError("locked_building_block_ids must be a list of nonempty catalog IDs")
    # Reuse the existing public property vocabulary, validating before any expensive work.
    from chembounce import _validate_parameters

    _validate_parameters(0.5, 0.3, 42, config.get("property_thresholds", {}))
    return config


def load_config(filename):
    path = Path(filename).resolve()
    config = _config(json.loads(path.read_text()))
    for key in ("catalog", "rules", "reference_sdf"):
        if config.get(key):
            config[key] = str((path.parent / config[key]).resolve())
    provider = config["provider"] = dict(config["provider"])
    for key in ("model_path", "prompt_path", "path"):
        if provider.get(key):
            provider[key] = str((path.parent / provider[key]).resolve())
    return config


def _provider(config):
    config = dict(config)
    kind = config.pop("kind", None)
    if kind == "saved":
        if set(config) != {"path"}:
            raise ValueError("saved provider requires only path")
        return SavedProposals(config["path"])
    if kind == "synllama":
        if set(config) - {"model_path", "prompt_path", "device", "samples", "max_new_tokens", "seed"}:
            raise ValueError("Unknown SynLlama provider option")
        return SynLlama(**config)
    raise ValueError("Provider kind must be synllama or saved")


def _source_snapshot():
    return {p.name: digest(p) for p in sorted(Path(__file__).parent.glob("*.py"))}


# Snapshot at import, as with the main workflow, to detect edits in a running process.
_SOURCE_HASHES = _source_snapshot()


def _manifest(config, hit, designs):
    try:
        root = Path(__file__).parent
        revision = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, stderr=subprocess.DEVNULL, text=True
        ).strip()
        dirty = bool(
            subprocess.check_output(
                ["git", "status", "--porcelain"], cwd=root, stderr=subprocess.DEVNULL, text=True
            ).strip()
        )
    except (OSError, subprocess.CalledProcessError):
        revision, dirty = None, None
    return {
        "schema_version": 1,
        "hit": hit,
        "designs": designs,
        "config": config,
        "code_sha256": _SOURCE_HASHES,
        "sources_changed_since_import": _source_snapshot() != _SOURCE_HASHES,
        "git_revision": revision,
        "git_dirty": dirty,
        "versions": {name: importlib.metadata.version(name) for name in ("rdkit", "numpy", "scipy")},
        "chemical_state_policy": "Explicit input states; no implicit protonation, tautomer, or stereoisomer enumeration",
        "verification_scope": "Catalog snapshot membership and registered reaction replay; no experimental synthesis or current-stock guarantee",
        "ranking_scope": "Heuristic, bounded search; joint shape/pharmacophore scores are not calibrated activity probabilities",
    }


def _pareto(rows):
    remaining = list(range(len(rows)))
    level = 0
    points = [
        (
            r["hit_shape"],
            r["hit_pharmacophore"],
            r["design_shape"],
            r["design_pharmacophore"],
            -r["relative_energy"],
            -r["reaction_steps"],
        )
        for r in rows
    ]
    while remaining:
        front = [
            i
            for i in remaining
            if not any(
                all(a >= b for a, b in zip(points[j], points[i]))
                and any(a > b for a, b in zip(points[j], points[i]))
                for j in remaining
                if j != i
            )
        ]
        for i in front:
            rows[i]["pareto_front"] = level
        removed = set(front)
        remaining = [i for i in remaining if i not in removed]
        level += 1


def run(hit, designs, config, output_dir):
    """Write complete evidence for each verified route, failure, and accepted pose.

    A fresh output directory is required. Configuration paths passed directly to
    this API are relative to the current directory; load_config resolves file paths.
    """
    out = Path(output_dir)
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Synthesis output directory is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    summary = {
        "schema_version": 1,
        "status": "running",
        "proposals": 0,
        "verified_routes": 0,
        "scored_products": 0,
        "accepted_products": 0,
        "failures": [],
        "search_reports": [],
        "limits_reached": [],
    }
    write_json(out / "summary.json", summary)
    route_records, score_records, raw_records = [], [], []
    rows, poses, score_cache = {}, {}, {}
    try:
        config = _config(config)
        hit = canonical(hit)
        all_designs = list(dict.fromkeys(canonical(s) for s in designs))
        designs = all_designs[: config.get("max_designs", 20)]
        if len(designs) < len(all_designs):
            summary["limits_reached"].append("max_designs")
        summary["designs_supplied"], summary["designs_used"] = len(all_designs), len(designs)
        if config.get("core_smiles") and not molecule(hit).HasSubstructMatch(
            molecule(config["core_smiles"]), useChirality=True
        ):
            raise ValueError("The original hit does not contain core_smiles")
        limits = SearchLimits(**config.get("search", {}))
        settings = ConformerSettings(**config.get("conformers", {}))
        if designs:
            require_alignment()
        catalog, rules = Catalog(config["catalog"]), Rules(config["rules"])
        manifest = _manifest(config, hit, designs)
        manifest["catalog_sha256"], manifest["rules_sha256"] = catalog.sha256, rules.sha256
        manifest["catalog_import"] = catalog.import_provenance
        manifest["effective_search"], manifest["effective_conformers"] = asdict(limits), asdict(settings)
        manifest["effective_workflow"] = {
            "mode": config.get("mode", "hop"),
            "core_smiles": config.get("core_smiles"),
            "preserve_design_scaffold": config.get("preserve_design_scaffold", False),
            "max_designs": config.get("max_designs", 20),
            "max_scored_products": config.get("max_scored_products", 256),
            "max_results": config.get("max_results", 50),
            "expansion": {
                "rounds": 0,
                "max_routes": 64,
                "locked_building_block_ids": [],
                **config.get("expansion", {}),
            },
            "minimum_scores": config.get("minimum_scores", {}),
            "property_thresholds": config.get("property_thresholds", {}),
        }
        hit_pose = (
            reference_from_sdf(config["reference_sdf"], hit, settings)
            if config.get("reference_sdf")
            else None
        )
        if config.get("reference_sdf"):
            manifest["reference_sdf_sha256"] = digest(config["reference_sdf"])
        provider = _provider(config["provider"]) if designs else None
        manifest["provider"] = provider.provenance if provider else None
        write_json(out / "manifest.json", manifest)
        from utils import _get_molecular_prop_extd_, determ_props

        for design in designs:
            LOGGER.info("Projecting design %s onto catalog reaction space", design)
            proposals = provider.propose(design)
            for proposal in proposals:
                proposal_id = len(raw_records)
                raw_records.append({"design": design, **proposal})
                write_json(
                    out / "proposals.json",
                    {"schema_version": 1, "provider": provider.provenance, "proposals": raw_records},
                )
                summary["proposals"] += 1
                try:
                    projected, report = resolve_proposal(proposal["response"], catalog, rules, limits)
                    summary["search_reports"].append(
                        {"proposal_id": proposal_id, "stage": "projection", **asdict(report)}
                    )
                    if not projected:
                        summary["failures"].append(
                            {
                                "stage": "projection",
                                "proposal_id": proposal_id,
                                "reason": "No catalog-supported reaction route survived",
                            }
                        )
                except (ValueError, KeyError, TypeError) as exc:
                    summary["failures"].append(
                        {"stage": "proposal", "proposal_id": proposal_id, "reason": str(exc)}
                    )
                    continue
                candidates = list(projected)
                expansion = config.get("expansion", {})
                if expansion.get("rounds", 0):
                    for item in projected:
                        expanded, report = expand_route(item["route"], catalog, rules, limits, **expansion)
                        candidates.extend(expanded)
                        summary["search_reports"].append(
                            {"proposal_id": proposal_id, "stage": "expansion", **asdict(report)}
                        )
                for item in candidates:
                    route_id = f"route_{len(route_records):06d}"
                    # Final independent verification includes expansion products.
                    verification = verify_route(item["route"], catalog, rules, limits.reaction_product_limit)
                    product = verification["product"]
                    record = {
                        "route_id": route_id,
                        "proposal_id": proposal_id,
                        "design": design,
                        **item,
                        "verification": verification,
                        "selection_status": "pending",
                    }
                    route_records.append(record)
                    summary["verified_routes"] += 1
                    try:
                        check_constraints(
                            product,
                            hit,
                            design,
                            config.get("mode", "hop"),
                            config.get("core_smiles"),
                            config.get("preserve_design_scaffold", False),
                        )
                        props = _get_molecular_prop_extd_(molecule(product))
                        accepted, reasons = determ_props(
                            props, config.get("property_thresholds", {}), w_all_failed_reason=True
                        )
                        if not accepted:
                            raise ValueError("Product property filter: " + "; ".join(reasons))
                    except ValueError as exc:
                        record.update(selection_status="constraint_rejected", reason=str(exc))
                        continue
                    key = (product, design)
                    if key not in score_cache:
                        if summary["scored_products"] >= config.get("max_scored_products", 256):
                            record["selection_status"] = "not_scored_budget"
                            if "max_scored_products" not in summary["limits_reached"]:
                                summary["limits_reached"].append("max_scored_products")
                            continue
                        summary["scored_products"] += 1
                        try:
                            score_cache[key] = score_candidate(product, hit, design, settings, hit_pose)
                        except (ValueError, RuntimeError) as exc:
                            score_cache[key] = str(exc)
                    result = score_cache[key]
                    if isinstance(result, str):
                        record.update(selection_status="3d_failed", reason=result)
                        continue
                    info, hit_aligned, design_aligned = result
                    scores = {
                        f"{ref}_{metric}": info[ref][metric]
                        for ref in ("hit", "design")
                        for metric in ("shape", "pharmacophore")
                    }
                    if any(scores[k] < v for k, v in config.get("minimum_scores", {}).items()):
                        record.update(selection_status="3d_threshold_rejected", scores=scores)
                        continue
                    record["selection_status"] = "accepted"
                    score_records.append(
                        {
                            "route_id": route_id,
                            "product": product,
                            "design": design,
                            "properties": props,
                            **info,
                        }
                    )
                    row = {
                        "smiles": product,
                        "design": design,
                        "score": info["score"],
                        **scores,
                        "hit_morgan": DataStructs.TanimotoSimilarity(
                            FP.GetFingerprint(molecule(hit)), FP.GetFingerprint(molecule(product))
                        ),
                        "relative_energy": info["candidate_relative_energy"],
                        "reaction_steps": len(verification["steps"]),
                        "route_id": route_id,
                        "alternative_route_ids": [],
                    }
                    if product in rows:
                        previous = rows[product]
                        if (row["score"], -row["reaction_steps"]) > (
                            previous["score"],
                            -previous["reaction_steps"],
                        ):
                            row["alternative_route_ids"] = [
                                previous["route_id"],
                                *previous["alternative_route_ids"],
                            ]
                        else:
                            previous["alternative_route_ids"].append(route_id)
                            continue
                    rows[product] = row
                    poses[product] = (info, hit_aligned, design_aligned)
        ranked = list(rows.values())
        _pareto(ranked)
        ranked.sort(key=lambda r: (r["pareto_front"], -r["score"], r["reaction_steps"], r["smiles"]))
        if len(ranked) > config.get("max_results", 50):
            summary["limits_reached"].append("max_results")
        selected = ranked[: config.get("max_results", 50)]
        with (
            (out / "results.tsv").open("w", newline="") as stream,
            Chem.SDWriter(str(out / "poses_hit.sdf")) as hit_writer,
            Chem.SDWriter(str(out / "poses_design.sdf")) as design_writer,
            Chem.SDWriter(str(out / "references_hit.sdf")) as hit_ref_writer,
            Chem.SDWriter(str(out / "references_design.sdf")) as design_ref_writer,
        ):
            writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t")
            writer.writeheader()
            for rank, row in enumerate(selected, 1):
                row["rank"] = rank
                writer.writerow(row | {"alternative_route_ids": json.dumps(row["alternative_route_ids"])})
                info, hit_aligned, design_aligned = poses[row["smiles"]]
                for mol, sdf in ((hit_aligned, hit_writer), (design_aligned, design_writer)):
                    mol = Chem.Mol(mol)
                    mol.SetProp("_Name", row["route_id"])
                    mol.SetProp("smiles", row["smiles"])
                    mol.SetProp("score", str(row["score"]))
                    mol.SetProp("verification", "catalog_and_rule_verified; experimental_validation=false")
                    sdf.write(mol, confId=info["candidate_conformer_id"])
                for role, reference, sdf in (
                    ("hit", hit_pose if hit_pose is not None else ensemble(hit, settings), hit_ref_writer),
                    ("design", ensemble(row["design"], settings), design_ref_writer),
                ):
                    mol = Chem.Mol(reference.mol)
                    mol.SetProp("_Name", row["route_id"])
                    mol.SetProp("reference_role", role)
                    cid = info[role]["reference_conformer_id"]
                    mol.SetIntProp("source_conformer_id", cid)
                    sdf.write(mol, confId=cid)
        summary.update(
            status="success" if selected else "no_candidates",
            accepted_products=len(rows),
            returned_products=len(selected),
        )
        return selected
    except Exception as exc:
        summary.update(status="error", error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        summary["elapsed_seconds"] = time.monotonic() - start
        summary["route_rejections"] = {
            status: sum(r["selection_status"] == status for r in route_records)
            for status in sorted({r["selection_status"] for r in route_records})
        }
        write_json(out / "routes.json", {"schema_version": 1, "routes": route_records})
        write_json(out / "scores.json", {"schema_version": 1, "scores": score_records})
        write_json(out / "summary.json", summary)
