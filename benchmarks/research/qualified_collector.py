"""Integrity attestation of pinned trusted workflow artifacts, never new inference.

This does not turn arbitrary supplied tokens into proof of a live model call.
The source-bound child performed live cache validation; this collector checks its
saved identities, route replay, actual tensors, seed metrics and export agreement.
"""

import csv
import importlib.metadata
import math
from pathlib import Path
import platform
import time
from types import SimpleNamespace

from rdkit import Chem, DataStructs

from benchmarks.research.protocol import digest
from benchmarks.research.qualified_adapter import read, common_config
from chembounce_synthesis.common import FP, canonical, molecule
from chembounce_synthesis.forward_support import (
    GENERATION_SOURCE_SHA256,
    PADDING_POLICY,
    PROFILE,
    RUNTIME,
    PinnedReactionT5Forward,
    validate_bundle,
)
from chembounce_synthesis.identity import chemical_identity
from chembounce_synthesis.profiles import FIXED_CONFORMERS, METRICS, PROFILE_SHA256, load_profile
from chembounce_synthesis.qualification import _validate_score_info
from chembounce_synthesis.conformers import ConformerSettings
from chembounce_synthesis.route_audit import audit_route


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
    evidence = record["synthesis_evidence"]
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


def checked_products(
    directory, response, query, expected_source, *, decoder_factory=TokenReplay, deadline=None
):
    """Return collector-computed chemistry/robust scores from complete consistent exports."""
    from chembounce_synthesis.catalog import Catalog
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
    request = read(directory / "request.json")
    if request["query"] != query:
        raise ValueError("Common request differs from planned public query")
    bridge = directory / "selected-primary-bridge" if request["settings"]["arm"] == "baseline" else None
    from chembounce_synthesis.workflow import _config
    from chembounce_synthesis.forward_support import binding

    if binding(_config(common_config(request, bridge))) != binding(config):
        raise ValueError("Effective workflow config differs from frozen request plus declared arm overrides")
    if (
        len(rows) > 3
        or [{"smiles": r["smiles"], "route_id": r["route_id"]} for r in rows] != response["candidates"]
    ):
        raise ValueError("Response selected identities differ from final JSON")
    if (
        policy["mode"] != "qualified"
        or policy["evidence_profile"] != PROFILE
        or policy["conformers"] != FIXED_CONFORMERS
        or policy["calibration"]["artifact_sha256"] != PROFILE_SHA256["qualified-development-v1"]
    ):
        raise ValueError("Common final policy is not the frozen qualified protocol")
    qualification = read(out / "qualification.json")
    if qualification["policy"] != policy or qualification["selected_route_ids"] != [
        r["route_id"] for r in rows
    ]:
        raise ValueError("Live qualification selection/policy mismatch")
    if policy["minimum_scores"] != load_profile("qualified-development-v1")["minimum_scores"]:
        raise ValueError("Common minima differ from the frozen calibration artifact")
    if not rows:
        return []
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
        table = list(csv.DictReader(stream, delimiter="\t"))
    if len(table) != len(rows):
        raise ValueError("TSV cardinality differs from JSON")
    poses = {
        name: list(Chem.SDMolSupplier(str(out / name), removeHs=False))
        for name in ("poses_hit.sdf", "poses_design.sdf", "references_hit.sdf", "references_design.sdf")
    }
    if any(len(values) != len(rows) or any(m is None for m in values) for values in poses.values()):
        raise ValueError("SDF cardinality/chemistry differs from JSON")
    records = read(out / "routes.json")["routes"]
    refinements = read(out / "refinement.json")["candidates"]
    selected, decoder = [], None
    catalog, rules = Catalog(config["catalog"]), Rules(config["rules"])
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
        for key in (
            "smiles",
            "route_id",
            "evidence_profile",
            "route_support_binding",
            "recommendation_status",
        ):
            if table[index][key] != str(row[key]):
                raise ValueError("TSV selected binding differs from JSON")
        for name, values in poses.items():
            mol = values[index]
            if name.startswith("poses") and mol.GetProp("smiles") != row["smiles"]:
                raise ValueError("SDF declared product differs from JSON")
            if not mol.GetConformer().Is3D() or not all(
                math.isfinite(float(v)) for xyz in mol.GetConformer().GetPositions() for v in xyz
            ):
                raise ValueError("SDF pose has invalid coordinates")
            if name.startswith("poses") and any(
                not math.isclose(float(mol.GetProp(key)), row[key], rel_tol=1e-12, abs_tol=1e-12)
                for key in ("score", "hit_tanimoto", "hit_shape", "hit_pharmacophore")
            ):
                raise ValueError("SDF pose scores differ from JSON")
            expected = (
                row["smiles"]
                if name.startswith("poses")
                else manifest["hit"]
                if name == "references_hit.sdf"
                else row["design"]
            )
            if (
                canonical(Chem.MolToSmiles(mol)) != expected
                or mol.GetProp("_Name") != row["route_id"]
                or mol.GetProp("route_support_binding") != row["route_support_binding"]
                or int(mol.GetProp("pose_seed")) != row["pose_seed"]
            ):
                raise ValueError("SDF selected identity/route/pose differs from JSON")
        identity = chemical_identity(row["smiles"])
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
            }
        )
    check_deadline()
    return selected
