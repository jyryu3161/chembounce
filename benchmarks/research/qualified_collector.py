"""Integrity attestation of pinned trusted workflow artifacts, never new inference.

This does not turn arbitrary supplied tokens into proof of a live model call.
The source-bound child performed live cache validation; this collector checks its
saved identities, route replay, actual tensors, seed metrics and export agreement.
"""

import csv
from collections import Counter
from dataclasses import asdict, dataclass, field
import importlib.metadata
import math

import numpy as np
from pathlib import Path
import platform
import statistics
import time
from types import MappingProxyType, SimpleNamespace

from rdkit import Chem, DataStructs, rdBase

from benchmarks.research.protocol import digest
from benchmarks.research.qualified_adapter import (
    read,
    common_config,
    benchmark_profile,
    selective_benchmark_profile,
)
from chembounce_synthesis.common import FP, canonical, molecule
from chembounce_synthesis.catalog import Catalog, _import_manifest_present
from chembounce_synthesis.forward_support import (
    GENERATION_SOURCE_SHA256,
    PADDING_POLICY,
    PROFILE,
    RUNTIME,
    PinnedReactionT5Forward,
    validate_bundle,
)
from chembounce_synthesis.identity import chemical_identity
from chembounce_synthesis.profiles import (
    FIXED_CONFORMERS,
    METRICS,
    PROFILE_SHA256,
    load_profile,
)
from chembounce_synthesis.qualification import _validate_score_info
from chembounce_synthesis.conformers import ConformerSettings
from chembounce_synthesis.route_audit import audit_route


@dataclass(frozen=True, init=False)
class BoundCatalog:
    """Runner-owned full snapshot with only read-only collector lookup exposed.

    The complete Catalog remains owned here, but no molecule, fingerprint,
    nearest-neighbor or reaction-result cache is available through this API.
    This object is never global, serialized, or reused across campaigns.
    """

    path: Path
    sha256: str
    _sidecar_sha256: object
    _source: object
    _runtime: tuple
    _collector_sha256: str
    _catalog: object = field(repr=False, compare=False)
    _records: object = field(repr=False, compare=False)

    def __init__(self, filename, expected_sha256, expected_source, *, sidecar_sha256):
        object.__setattr__(self, "path", Path(filename).resolve())
        object.__setattr__(self, "sha256", expected_sha256)
        object.__setattr__(self, "_sidecar_sha256", sidecar_sha256)
        source = {name: expected_source[name] for name in ("catalog.py", "common.py")}
        object.__setattr__(self, "_source", MappingProxyType(source))
        object.__setattr__(self, "_runtime", (platform.python_version(), rdBase.rdkitVersion))
        object.__setattr__(self, "_collector_sha256", digest(__file__))
        self.validate(filename, expected_source)
        catalog = Catalog(self.path)  # Full eager validation, never selective parsing.
        if catalog.sha256 != expected_sha256:
            raise ValueError("Catalog changed during shared construction")
        object.__setattr__(self, "_catalog", catalog)
        # Frozen BuildingBlock values and a private immutable mapping prevent a
        # lookup caller from modifying the retained authoritative record set.
        object.__setattr__(self, "_records", MappingProxyType(dict(catalog.blocks)))
        self.validate(filename, expected_source)

    def validate(self, filename, expected_source):
        if Path(filename).resolve() != self.path or digest(self.path) != self.sha256:
            raise ValueError("Shared catalog path or checksum differs from the pinned snapshot")
        sidecar = self.path.with_suffix(self.path.suffix + ".manifest.json")
        actual = digest(sidecar) if _import_manifest_present(sidecar) else None
        if actual != self._sidecar_sha256:
            raise ValueError("Shared catalog sidecar presence/checksum changed")
        source_root = Path(__file__).resolve().parents[2] / "chembounce_synthesis"
        if digest(__file__) != self._collector_sha256:
            raise ValueError("Shared catalog collector source changed")
        if any(
            expected_source.get(name) != value or digest(source_root / name) != value
            for name, value in self._source.items()
        ) or self._runtime != (platform.python_version(), rdBase.rdkitVersion):
            raise ValueError("Shared catalog source/runtime binding changed")

    def get(self, ident):
        block = self._records.get(ident)
        if block is None or block.availability == "unavailable":
            raise ValueError(f"Building block is absent or unavailable: {ident}")
        return block

    def evidence(self, ident):
        expected = asdict(self.get(ident)) | {"catalog_sha256": self.sha256}
        evidence = self._catalog.evidence(ident)
        if not isinstance(evidence, dict) or any(
            evidence.get(key) != value for key, value in expected.items()
        ):
            raise ValueError("Shared catalog evidence differs from its frozen records")
        return evidence

    def _check_import_unchanged(self):
        self._catalog._check_import_unchanged()

    def binding(self):
        return {
            "catalog_path": str(self.path),
            "catalog_sha256": self.sha256,
            "catalog_manifest_path": str(self.path.with_suffix(self.path.suffix + ".manifest.json")),
            "catalog_manifest_sha256": self._sidecar_sha256,
            "source_sha256": dict(self._source),
            "collector_source_sha256": self._collector_sha256,
            "python": self._runtime[0],
            "rdkit": self._runtime[1],
            "records": len(self._records),
        }


class TokenReplay(PinnedReactionT5Forward):
    """Pinned tokenizer/config only. No weights loaded and no evaluator/cache created."""

    def __init__(self, model_path, provenance):
        ignored = validate_bundle(model_path)
        path = Path(model_path).resolve()
        generation = Path(
            importlib.metadata.distribution("transformers").locate_file("transformers/generation/utils.py")
        )
        runtime = {name: importlib.metadata.version(name) for name in RUNTIME}
        if (
            runtime != RUNTIME
            or platform.python_version() != "3.10.12"
            or digest(generation) != GENERATION_SOURCE_SHA256
        ):
            raise ValueError("Collector numerical runtime differs from pinned inference runtime")
        import chembounce_synthesis.forward_support as source
        from transformers import AutoTokenizer

        expected = {
            "model_path": str(path),
            "profile": PROFILE,
            "runtime": RUNTIME,
            "padding_policy": PADDING_POLICY,
            "adapter_sha256": digest(source.__file__),
            "base_adapter_sha256": digest(Path(source.__file__).with_name("evidence.py")),
            "ignored_non_model_metadata_sha256": ignored,
            "model_files_sha256": {p.name: digest(p) for p in sorted(path.iterdir()) if p.is_file()},
            "top_k": 3,
            "max_new_tokens": 256,
            "max_input_tokens": 512,
            "device": "cpu",
            "do_sample": False,
            "torch_threads": 1,
        }
        if any(provenance.get(key) != value for key, value in expected.items()):
            raise ValueError("Saved forward provenance differs from actual pinned assets/runtime")
        if provenance.get("generation_source", {}).get("sha256") != GENERATION_SOURCE_SHA256:
            raise ValueError("Saved generation-source binding is missing")
        self.tokenizer = AutoTokenizer.from_pretrained(
            str(path), local_files_only=True, trust_remote_code=False
        )
        self.model = SimpleNamespace(config=SimpleNamespace(**read(path / "config.json")))
        self.top_k, self.max_new_tokens, self.max_input_tokens = 3, 256, 512
        self.provenance = provenance

    def predict_raw(self, prepared):
        raise RuntimeError("Collector is tokenizer-only; inference is forbidden")


def verify_route_support(record, catalog, rules, limit, decoder):
    from chembounce_synthesis.evidence_outcomes import validate_recorded

    evidence = record["synthesis_evidence"]
    validate_recorded(evidence)
    audit = audit_route(record["route"], catalog, rules, limit)
    if (
        evidence.get("profile") != PROFILE
        or not audit["strict_eligible"]
        or evidence.get("route_audit") != audit
        or evidence.get("route_binding_sha256") != audit["route_binding_sha256"]
        or evidence.get("product") != audit["product"]
        or evidence.get("predictor_provenance") != decoder.provenance
        or len(evidence.get("steps", [])) != len(audit["steps"])
    ):
        raise ValueError("Selected route support differs from fresh concrete replay")
    if record["verification"] != audit["verification"]:
        raise ValueError("Selected verification differs from authoritative replay")
    for step, saved in zip(audit["steps"], evidence["steps"]):
        if saved["node_id"] != step["node_id"] or saved["step_binding_sha256"] != step["step_binding_sha256"]:
            raise ValueError("Selected step order/binding changed")
        prepared = decoder.prepare(
            [r["smiles"] for r in step["reactant_slots"]],
            saved.get("context", {}).get("supplied", {}).get("reagents", []),
        )
        prediction = saved["forward_prediction"]
        beams = decoder.validate_raw(prediction["raw"], prepared)
        matches = [b["rank"] for b in beams if b["complete"] and b["product"] == step["product"]]
        if (
            not matches
            or prediction.get("status") != "matched"
            or prediction.get("matched_rank") != min(matches)
            or prediction.get("predictions") != beams
        ):
            raise ValueError("Selected step lacks completed exact raw forward support")
    return audit


def _selective_outcome(summary, rows, records, refinements, qualification):
    """Cross-check saved outcome/denominators, including genuine empty runs."""
    from chembounce_synthesis.evidence_outcomes import validate_records, validate_summary

    validate_records(records, qualification["decisions"])
    validate_summary(summary, records)
    route_ids = [r["route_id"] for r in records]
    accepted = summary.get("accepted_products")
    if (
        len(route_ids) != len(set(route_ids))
        or qualification["decisions"] != [r["qualification"] for r in records]
        or any(
            type(r["qualification"].get("qualified")) is not bool
            or (r["qualification"].get("status") == "qualified") != r["qualification"]["qualified"]
            or (r["selection_status"] == "qualified") != r["qualification"]["qualified"]
            for r in records
        )
        or summary.get("route_rejections") != dict(Counter(r["selection_status"] for r in records))
        or summary.get("status") != ("success" if rows else "no_candidates")
        or summary.get("returned_products") != len(rows)
        or type(accepted) is not int
        or accepted < len(rows)
        or summary.get("qualified_products") != accepted
        or summary.get("refinement_accepted") != accepted
        or summary.get("refinement_attempted") != sum("alternative_index" in r for r in refinements)
        or any(r["route_id"] not in route_ids for r in refinements)
        or summary.get("refinement_budget", {}).get("candidate_conformer_requests")
        != sum(r["candidate_conformer_requests"] for r in refinements)
        or (not rows and (accepted != 0 or any(r["qualification"]["qualified"] for r in records)))
    ):
        raise ValueError("Selective outcome/route/refinement denominator or decision mismatch")


def _representative_fields(representative, trials, policy, design_tanimoto, hit_tanimoto):
    """Derive export metadata from retained validated trials, never new poses."""
    scores = [t["score"] for t in trials]
    return {
        "score": representative["score"],
        "relative_energy": representative["candidate_relative_energy"],
        "scoring_objective": representative["settings"]["scoring_objective"],
        "ranking_policy": policy["ranking_policy"],
        "scoring_stage": "refined",
        "qualification_policy": policy["policy_version"],
        "evidence_profile": policy["evidence_profile"],
        "pose_seed": sorted(trials, key=lambda t: (t["score"], t["seed"]))[1]["seed"],
        "refinement_conformers": representative["settings"]["num_conformers"],
        "hit_tanimoto": hit_tanimoto,
        "hit_morgan": hit_tanimoto,
        "design_tanimoto": design_tanimoto,
        **{
            f"{role}_{metric}": representative[role][metric]
            for role in ("hit", "design")
            for metric in ("shape", "pharmacophore")
        },
        "seed_score_median": statistics.median(scores),
        "seed_score_min": min(scores),
        "seed_score_max": max(scores),
        "seed_score_range": max(scores) - min(scores),
        "seed_score_std": statistics.pstdev(scores),
        "seed_count": len(trials),
        "seed_successes": len(trials),
        "ranking_score": min(scores),
        "ranking_metrics": {
            f"{role}_{metric}": min(t[role][metric] for t in trials)
            for role in ("hit", "design")
            for metric in ("shape", "pharmacophore")
        },
        "convergence_status": "not_assessed",
    }


def checked_products(
    directory,
    response,
    query,
    expected_source,
    *,
    decoder_factory=TokenReplay,
    deadline=None,
    bound_catalog=None,
    benchmark_protocol=None,
):
    """Authenticate the unchanged benchmark envelope before shared artifact checks."""
    if deadline is not None and time.monotonic() >= deadline:
        raise TimeoutError("Outer budget exhausted during collector integrity checks")
    directory = Path(directory)
    request = read(directory / "request.json")
    if request.get("benchmark_protocol") != benchmark_protocol:
        raise ValueError("Common request benchmark profile differs from frozen runner authority")
    if request["query"] != query:
        raise ValueError("Common request differs from planned public query")
    bridge = directory / "selected-primary-bridge" if request["settings"]["arm"] == "baseline" else None
    from chembounce_synthesis.workflow import _config

    return checked_workflow_artifacts(
        directory,
        response,
        query,
        expected_source,
        expected_config=_config(common_config(request, bridge)),
        profile_name=benchmark_profile(request.get("benchmark_protocol")),
        decoder_factory=decoder_factory,
        deadline=deadline,
        bound_catalog=bound_catalog,
    )


def checked_workflow_artifacts(
    directory,
    response,
    query,
    expected_source,
    *,
    expected_config,
    profile_name,
    decoder_factory=TokenReplay,
    deadline=None,
    bound_catalog=None,
):
    """Shared checks after a caller authenticates its own explicit config authority.

    This does not authenticate a benchmark or development request by itself.
    Both callers must bind their complete expected config, query and source map.
    """
    from chembounce_synthesis.reactions import Rules
    from chembounce_synthesis.workflow import check_constraints

    def check_deadline():
        if deadline is not None and time.monotonic() >= deadline:
            raise TimeoutError("Outer budget exhausted during collector integrity checks")

    check_deadline()
    directory = Path(directory)
    files = response["metadata"]["files_sha256"]
    required = {
        "manifest.json",
        "results.json",
        "qualification.json",
        "routes.json",
        "refinement.json",
        "summary.json",
    }
    if not {"synthesis/" + name for name in required} <= files.keys():
        raise ValueError("Incomplete common artifact manifest")
    if any(
        not (directory / name).resolve().is_relative_to(directory.resolve())
        or digest(directory / name) != value
        for name, value in files.items()
    ):
        raise ValueError("Common artifacts changed after execution")
    out = directory / "synthesis"
    manifest, result = read(out / "manifest.json"), read(out / "results.json")
    if manifest["sources_changed_since_import"] or manifest["code_sha256"] != expected_source:
        raise ValueError("Workflow source differs from frozen collector authority")
    if manifest["hit"] != canonical(query["smiles"]):
        raise ValueError("Common hit differs from planned query")
    config, rows, policy = manifest["config"], result["results"], result["policy"]
    from chembounce_synthesis.forward_support import binding

    if binding(expected_config) != binding(config):
        raise ValueError("Effective workflow config differs from authenticated expected config")
    if (
        len(rows) > 3
        or [{"smiles": r["smiles"], "route_id": r["route_id"]} for r in rows] != response["candidates"]
    ):
        raise ValueError("Response selected identities differ from final JSON")
    name = profile_name
    selective = selective_benchmark_profile(name)
    similarity = None
    if selective:
        from chembounce_synthesis import selective_scoring
        from chembounce_synthesis.qualification import effective_policy

        similarity = selective_scoring.configured_constraints(config, manifest["hit"])
        expected_policy = effective_policy(config)
        expected_policy["similarity_constraints"] = {"policy": similarity.policy, "sha256": similarity.sha256}
        if (
            policy != expected_policy
            or manifest.get("recommendation_policy") != expected_policy
            or manifest.get("recommendation_profile") != expected_policy["calibration"]
            or any(
                expected_source.get(Path(path).name) != value
                for path, value in expected_policy["calibration"]["current_numerical_runtime"][
                    "source_sha256"
                ].items()
            )
        ):
            raise ValueError("Selective policy/provenance differs from frozen effective protocol/source")
    if (
        policy["mode"] != "qualified"
        or policy["evidence_profile"] != PROFILE
        or policy["conformers"] != FIXED_CONFORMERS
        or policy["calibration"]["artifact_sha256"] != PROFILE_SHA256[name]
    ):
        raise ValueError("Common final policy is not the frozen qualified protocol")
    qualification = read(out / "qualification.json")
    if qualification["policy"] != policy or qualification["selected_route_ids"] != [
        r["route_id"] for r in rows
    ]:
        raise ValueError("Live qualification selection/policy mismatch")
    if policy["minimum_scores"] != load_profile(name)["minimum_scores"]:
        raise ValueError("Common minima differ from the frozen calibration artifact")
    if bound_catalog is not None:
        if type(bound_catalog) is not BoundCatalog:
            raise ValueError("Collector requires a runner-owned bound catalog")
        bound_catalog.validate(config["catalog"], expected_source)
    from chembounce_synthesis.evidence_outcomes import validate_records, validate_summary

    records = read(out / "routes.json")["routes"]
    summary = read(out / "summary.json")
    validate_records(records, qualification["decisions"])
    validate_summary(summary, records)
    if not rows and not selective:
        return []  # Legacy emptiness remains valid after checking any optional diagnostic fields.
    refinements = read(out / "refinement.json")["candidates"]
    if selective:
        _selective_outcome(summary, rows, records, refinements, qualification)
    exports = {
        "results.tsv",
        "poses_hit.sdf",
        "poses_design.sdf",
        "references_hit.sdf",
        "references_design.sdf",
    }
    if not {"synthesis/" + name for name in exports} <= files.keys():
        raise ValueError("Selected output lacks complete bound exports")
    with (out / "results.tsv").open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        table = list(reader)
        if selective:
            from chembounce_synthesis.workflow import COLUMNS, SELECTIVE_COLUMNS

            if reader.fieldnames != COLUMNS + SELECTIVE_COLUMNS:
                raise ValueError("Selective TSV header differs from the frozen export schema")
    if len(table) != len(rows):
        raise ValueError("TSV cardinality differs from JSON")
    if not rows:
        if any((out / name).read_text().strip() for name in exports - {"results.tsv"}):
            raise ValueError("Empty selective result has nonempty pose/reference exports")
        check_deadline()
        return []
    poses = {
        name: list(Chem.SDMolSupplier(str(out / name), removeHs=False))
        for name in ("poses_hit.sdf", "poses_design.sdf", "references_hit.sdf", "references_design.sdf")
    }
    if any(len(values) != len(rows) or any(m is None for m in values) for values in poses.values()):
        raise ValueError("SDF cardinality/chemistry differs from JSON")
    selected, decoder = [], None
    catalog = bound_catalog if bound_catalog is not None else Catalog(config["catalog"])
    rules = Rules(config["rules"])  # Fresh per call: no cross-attempt Rules.apply cache.
    for index, row in enumerate(rows):
        check_deadline()
        if row.get("rank") != index + 1 or row.get("recommendation_status") != "qualified":
            raise ValueError("Selected output lacks live qualification")
        found = [r for r in records if r["route_id"] == row["route_id"]]
        if len(found) != 1 or found[0]["qualification"].get("qualified") is not True:
            raise ValueError("Selected route lacks unique live qualification decision")
        record = found[0]
        if record["qualification"] not in qualification["decisions"]:
            raise ValueError("Live qualification decision differs across artifacts")
        if decoder is None:
            decoder = decoder_factory(
                config["synthesis_evidence"]["predictor"]["model_path"],
                record["synthesis_evidence"]["predictor_provenance"],
            )
        audit = verify_route_support(
            record, catalog, rules, config["search"]["reaction_product_limit"], decoder
        )
        if audit["product"] != row["smiles"] or audit["route_binding_sha256"] != row["route_support_binding"]:
            raise ValueError("Row identity differs from selected route support")
        constraint = policy["constraints"]
        check_constraints(
            row["smiles"],
            manifest["hit"],
            row["design"],
            **{key: constraint[key] for key in ("mode", "core_smiles", "preserve_design_scaffold")},
        )
        matches = [
            r
            for r in refinements
            if r["product"] == row["smiles"] and r["design"] == row["design"] and r.get("status") == "success"
        ]
        if len(matches) != 1:
            raise ValueError("Selected product lacks unique completed refinement")
        report = matches[0]
        trials = report["trials"]
        if [t.get("seed") for t in trials] != [42, 43, 44] or report["threshold_policy"] != "all_seeds":
            raise ValueError("Final seed protocol changed or incomplete")
        for trial in trials:
            _validate_score_info(trial, ConformerSettings(**(FIXED_CONFORMERS | {"seed": trial["seed"]})))
            if (
                trial.get("status") != "success"
                or trial.get("candidate_smiles") != row["smiles"]
                or trial.get("hit_smiles") != manifest["hit"]
                or trial.get("design_smiles") != row["design"]
            ):
                raise ValueError("Seed failure or scoring identity mismatch")
        if selective:
            if (
                record.get("design") != row["design"]
                or record["qualification"].get("similarity_constraints") != policy["similarity_constraints"]
                or record["qualification"].get("calibration")
                != {key: policy["calibration"][key] for key in ("name", "status", "artifact_sha256")}
                or report.get("num_conformers") != 8
                or report.get("seed_count") != 3
                or report.get("seed_successes") != 3
                or report.get("failed_thresholds") != []
            ):
                raise ValueError("Selective decision or complete seed report mismatch")
            for trial in trials:
                selective_scoring.validate_result(
                    selective_scoring.restore_result(trial),
                    row["smiles"],
                    ConformerSettings(**(FIXED_CONFORMERS | {"seed": trial["seed"]})),
                    hit=manifest["hit"],
                    design=row["design"],
                    similarity=similarity,
                )
        fp = FP.GetFingerprint(molecule(row["smiles"]))
        metrics = {
            "hit_tanimoto": DataStructs.TanimotoSimilarity(fp, FP.GetFingerprint(molecule(manifest["hit"]))),
            **{"hit_" + key: min(t["hit"][key] for t in trials) for key in ("shape", "pharmacophore")},
        }
        if metrics["hit_tanimoto"] != row["hit_tanimoto"]:
            raise ValueError("Saved original-hit fingerprint score is inconsistent")
        margins = {key: metrics[key] - policy["minimum_scores"][key] for key in METRICS}
        if any(not math.isfinite(v) or v < 0 for v in margins.values()):
            raise ValueError("A required robust original-hit minimum failed")
        representative = next((t for t in trials if t["seed"] == row["pose_seed"]), None)
        if representative is None or any(
            row["hit_" + key] != representative["hit"][key] for key in ("shape", "pharmacophore")
        ):
            raise ValueError("Exported score is not the representative seed's pose score")
        selective_fields, representative_fields = {}, {}
        if selective:
            from chembounce_synthesis.similarity_constraints import binding

            representative_fields = _representative_fields(
                representative,
                trials,
                policy,
                DataStructs.TanimotoSimilarity(fp, FP.GetFingerprint(molecule(row["design"]))),
                metrics["hit_tanimoto"],
            )
            if any(
                row.get(key) != value or isinstance(row.get(key), bool)
                for key, value in representative_fields.items()
            ):
                raise ValueError("Selective final row differs from retained representative metrics/protocol")
            selective_fields = selective_scoring.result_fields(representative)
            if (
                any(row.get(key) != value for key, value in selective_fields.items())
                or row["pose_seed"] != sorted(trials, key=lambda t: (t["score"], t["seed"]))[1]["seed"]
                or row["score"] != representative["score"]
                or row.get("recommendation_profile") != name
                or row.get("calibration_status") != policy["calibration"]["status"]
                or row.get("calibration_artifact_sha256") != PROFILE_SHA256[name]
            ):
                raise ValueError("Selective representative identity/policy/witness mismatch")
        for key in (
            "smiles",
            "route_id",
            "evidence_profile",
            "route_support_binding",
            "recommendation_status",
        ):
            if table[index][key] != str(row[key]):
                raise ValueError("TSV selected binding differs from JSON")
        if selective:
            for key in (
                "pose_seed",
                "recommendation_profile",
                "calibration_status",
                "calibration_artifact_sha256",
                "scoring_objective",
                "ranking_policy",
                *representative_fields,
            ):
                if key == "ranking_metrics":
                    import json

                    if json.loads(table[index].get(key, "null")) != row[key]:
                        raise ValueError("TSV selective ranking metrics differ from JSON")
                    continue
                if table[index].get(key) != str(row[key]):
                    raise ValueError("TSV selective representative/profile differs from JSON")
        for key, value in selective_fields.items():
            if table[index].get(key) != (str(value) if value is not None else ""):
                raise ValueError("TSV selective evidence differs from JSON")
        for pose_name, values in poses.items():
            mol = values[index]
            if pose_name.startswith("poses") and mol.GetProp("smiles") != row["smiles"]:
                raise ValueError("SDF declared product differs from JSON")
            if not mol.GetConformer().Is3D() or not all(
                math.isfinite(float(v)) for xyz in mol.GetConformer().GetPositions() for v in xyz
            ):
                raise ValueError("SDF pose has invalid coordinates")
            if pose_name.startswith("poses") and any(
                not math.isclose(float(mol.GetProp(key)), row[key], rel_tol=1e-12, abs_tol=1e-12)
                for key in ("score", "hit_tanimoto", "hit_shape", "hit_pharmacophore")
            ):
                raise ValueError("SDF pose scores differ from JSON")
            expected = (
                row["smiles"]
                if pose_name.startswith("poses")
                else manifest["hit"]
                if pose_name == "references_hit.sdf"
                else row["design"]
            )
            if (
                canonical(Chem.MolToSmiles(mol)) != expected
                or mol.GetProp("_Name") != row["route_id"]
                or mol.GetProp("route_support_binding") != row["route_support_binding"]
                or int(mol.GetProp("pose_seed")) != row["pose_seed"]
            ):
                raise ValueError("SDF selected identity/route/pose differs from JSON")
            if selective:
                shared = ("scoring_objective", "ranking_policy", "evidence_profile")
                candidate_text = ("scoring_stage", "recommendation_status", "qualification_policy")
                for key in shared + (candidate_text if pose_name.startswith("poses") else ()):
                    if not mol.HasProp(key) or mol.GetProp(key) != row[key]:
                        raise ValueError("SDF representative objective/status differs from final evidence")
                if pose_name.startswith("poses"):
                    if (
                        not mol.HasProp("verification")
                        or mol.GetProp("verification")
                        != "catalog_and_rule_verified; experimental_validation=false"
                    ):
                        raise ValueError("SDF verification scope differs from computational evidence")
                    for key in (
                        "score",
                        "hit_tanimoto",
                        "design_tanimoto",
                        "hit_shape",
                        "hit_pharmacophore",
                        "design_shape",
                        "design_pharmacophore",
                        "seed_score_median",
                        "seed_score_min",
                        "seed_score_max",
                    ):
                        if not mol.HasProp(key) or not math.isclose(
                            float(mol.GetProp(key)), representative_fields[key], rel_tol=1e-12, abs_tol=1e-12
                        ):
                            raise ValueError("SDF representative metric differs from retained evidence")
                else:
                    role = "hit" if pose_name == "references_hit.sdf" else "design"
                    if (
                        not mol.HasProp("source_conformer_id")
                        or mol.GetProp("source_conformer_id")
                        != str(representative[role]["reference_conformer_id"])
                        or not mol.HasProp("reference_role")
                        or mol.GetProp("reference_role") != role
                    ):
                        raise ValueError("SDF reference role/conformer differs from retained evidence")
                for key, value in selective_fields.items():
                    if not mol.HasProp(key) or mol.GetProp(key) != (
                        str(value) if value is not None else "not_applicable"
                    ):
                        raise ValueError("SDF selective policy/witness differs from JSON")
                saved_poses = representative["selective_pose_evidence"]
                saved = {
                    "poses_hit.sdf": saved_poses["hit_candidate"],
                    "poses_design.sdf": saved_poses["design_candidate"],
                    "references_hit.sdf": representative["similarity_constraints"]["witness"][
                        "reference_pose"
                    ],
                    "references_design.sdf": saved_poses["design_reference"],
                }[pose_name]
                from chembounce_synthesis.similarity_constraints import indexed_topology

                # SDF writes four decimals; exact hashes apply to retained JSON only.
                if (
                    indexed_topology(mol) != saved["atom_order_sha256"]
                    or not np.allclose(
                        mol.GetConformer().GetPositions(), saved["coordinates"], rtol=0, atol=0.000051
                    )
                    or any(
                        mol.GetProp(key) != row[key]
                        for key in (
                            "recommendation_profile",
                            "calibration_status",
                            "calibration_artifact_sha256",
                        )
                    )
                ):
                    raise ValueError("SDF retained representative geometry/provenance mismatch")
        identity = chemical_identity(row["smiles"])
        selective_evidence = {}
        if selective:
            polar = [t["similarity_constraints"]["witness"]["polar"]["recall"] for t in trials]
            selective_evidence = {
                "selective_evidence": {
                    "profile": name,
                    "profile_artifact_sha256": PROFILE_SHA256[name],
                    "policy_sha256": similarity.sha256,
                    "pose_seed": row["pose_seed"],
                    "all_seed_min_polar_reference_recall": min(polar)
                    if all(v is not None for v in polar)
                    else None,
                    "polar_status": "not_applicable" if all(v is None for v in polar) else "assessed",
                    "seed_witness_sha256": {
                        str(t["seed"]): binding(t["similarity_constraints"]["witness"]) for t in trials
                    },
                    "seed_pose_evidence_sha256": {
                        str(t["seed"]): binding(t["selective_pose_evidence"]) for t in trials
                    },
                }
            }
        selected.append(
            {
                "rank": index + 1,
                "smiles": identity.exact_smiles,
                "family": identity.family_smiles,
                "scaffold": identity.scaffold_family_smiles,
                "generic_scaffold": identity.ring_skeleton_smiles,
                "bit_morgan_bits": list(fp.GetOnBits()),
                "metrics": metrics,
                "margins": margins,
                **selective_evidence,
            }
        )
    if bound_catalog is not None:
        bound_catalog.validate(config["catalog"], expected_source)
    check_deadline()
    return selected
