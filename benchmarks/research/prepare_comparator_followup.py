"""Prepare, but never launch, the separately declared 270-attempt comparison.

Only finalized, hashed method fragments are accepted. No placeholder protocol,
model import, DB construction, benchmark subprocess, retry or resume is provided.
"""

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess


DECISION_SHA = "d18a3fc280201f933b9c2f7a0ae4151c8bdf165cf9577e947c8d16deb8606faf"
AMENDMENT_SHA = "7f43d5b06f5525005ba851bffcfbbf70f9f362ac8916ef1c81a0f99a1634654f"
ORIGINAL_PROTOCOLS = {
    "scaffold": (
        "main-protocol-gpu0.json",
        "6b8ee2ad718e23fe95f370679755e8202687f8d36059c2e58e2dee6f3b3ba7c7",
    ),
    "time": (
        "secondary-time-protocol-gpu0.json",
        "28909a8f7bf8f8c38cc2f6548eabe663b787d3b6cd28e05dd78353080ed55840",
    ),
    "target": (
        "secondary-target-protocol-gpu0.json",
        "02665eb636361a5c3e4a43b614a5a850cb405b57c45321b573c406a00918fb34",
    ),
}
SPLITS = ("scaffold", "time", "target")
COUNTS = {"scaffold": 12, "time": 12, "target": 6}
METHODS = ("pipeline_current", "crem_train_partial_exo", "mol2mol_pretrained")
SEEDS = [42, 43, 44]
BUDGET = {"wall_seconds": 360, "max_memory_mb": 24000, "max_candidates": 3, "max_evaluations": 256}
MOL2MOL_SETTINGS = {
    "protocol": "mol2mol-pretrained-v1-proposal",
    "raw_sample_cap": 256,
    "batch_size": 32,
    "temperature": 1.0,
    "strategy": "multinomial",
    "input_randomization": False,
    "device": "cuda:0",
    "filter_margin_seconds": 15.0,
    "batch_reserve_multiplier": 1.25,
}
THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "USE_TF": "0",
    "PYTHONNOUSERSITE": "1",
    "PYTHONDONTWRITEBYTECODE": "1",
    "PYTHONHASHSEED": "0",
    "TOKENIZERS_PARALLELISM": "false",
    "HF_HUB_OFFLINE": "1",
    "TRANSFORMERS_OFFLINE": "1",
    "HF_HUB_DISABLE_TELEMETRY": "1",
}
CONFORMERS = {
    "energy_window": 5.0,
    "max_iterations": 500,
    "num_conformers": 8,
    "pharmacophore_sigma": 1.0,
    "reject_unspecified_stereo": True,
    "strain_penalty": 0.05,
}
CREM_FILES = (
    "benchmarks/external/crem_baseline.py",
    "chembounce_synthesis/common.py",
    "chembounce_synthesis/identity.py",
    "fingerprint_db.py",
    "benchmarks/research/protocol.py",
    "benchmarks/research/runner.py",
    "benchmarks/research/evaluation.py",
    "chemistry.py",
    "chembounce_synthesis/exploration.py",
    "chembounce_synthesis/__init__.py",
    "benchmarks/research/__init__.py",
)
MOL2MOL_FILES = {
    **{
        "adapter_source_" + name: "benchmarks/external/mol2mol/" + name + ".py"
        for name in ("__init__", "adapter", "observation", "prepare", "selection", "support", "worker")
    },
    **{
        "common_source_" + name: "chembounce_synthesis/" + name + ".py"
        for name in ("__init__", "common", "identity")
    },
}


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def nonfinite(value):
        raise ValueError(f"Nonfinite JSON number: {value}")

    return json.loads(Path(path).read_text(), object_pairs_hook=unique, parse_constant=nonfinite)


def encode(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


class Pins:
    """Every input is checked again immediately before publishing any outputs."""

    def __init__(self):
        self.files = {}
        self.directories = {}
        self.python_trees = {}
        self.upstream_commits = {}
        self.runtimes = {}

    def file(self, path, expected=None):
        path = Path(path).absolute()
        if not path.is_file():
            raise ValueError(f"Required finalized file is absent: {path}")
        actual = self.files.get(str(path)) or digest(path)
        if expected is not None and actual != expected:
            raise ValueError(f"SHA-256 mismatch: {path}")
        if str(path) in self.files and self.files[str(path)] != actual:
            raise ValueError(f"Input changed while preparing: {path}")
        self.files[str(path)] = actual
        return {"artifact": str(path), "sha256": actual}

    def spec(self, spec, base):
        if set(spec) != {"artifact", "sha256"}:
            raise ValueError("Artifact requires exactly artifact and sha256")
        return self.file(Path(base) / spec["artifact"], spec["sha256"])

    def recheck(self):
        for path, expected in self.files.items():
            if digest(path) != expected:
                raise ValueError(f"Input changed before publication: {path}")
        for directory, expected in self.directories.items():
            actual = {str(p.relative_to(directory)) for p in Path(directory).rglob("*") if p.is_file()}
            if actual != set(expected):
                raise ValueError(f"Directory asset membership changed: {directory}")
        for directory, expected in self.python_trees.items():
            actual = {str(p.relative_to(directory)) for p in Path(directory).rglob("*.py") if p.is_file()}
            if actual != set(expected):
                raise ValueError(f"Upstream Python source membership changed: {directory}")
        for root, commit in self.upstream_commits.items():
            check_commit(root, commit)


def check_commit(root, commit):
    actual = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    if actual != commit:
        raise ValueError(f"Upstream source commit differs: {root}")
    subprocess.run(["git", "diff", "--quiet", commit, "--"], cwd=root, check=True)


def pin_inventory(root, hashes, pins):
    root = Path(root).resolve()
    if not isinstance(hashes, dict) or not hashes:
        raise ValueError("Prepared source inventory must be nonempty")
    for relative, sha in hashes.items():
        path = (root / relative).resolve()
        if Path(relative).is_absolute() or ".." in Path(relative).parts or not path.is_relative_to(root):
            raise ValueError("Source inventory path escapes its declared root")
        pins.file(path, sha)


def pin_upstream(root, commit, hashes, pins, *, generated=None, all_tracked=False):
    root = Path(root).resolve()
    check_commit(root, commit)
    tracked = set(
        subprocess.check_output(
            ["git", "ls-tree", "-r", "--name-only", commit],
            cwd=root,
            text=True,
        ).splitlines()
    )
    required = tracked if all_tracked else {p for p in tracked if p.endswith(".py")}
    if not required <= set(hashes) or (all_tracked and set(hashes) != required):
        raise ValueError("Prepared upstream inventory omits committed source files")
    generated = generated or {}
    if set(generated) - {"crem/_version.py"}:
        raise ValueError("Unexpected generated upstream source")
    pin_inventory(root, hashes, pins)
    if generated:
        pin_inventory(root, generated, pins)
    actual_python = {str(p.relative_to(root)) for p in root.rglob("*.py") if p.is_file()}
    if actual_python != {p for p in hashes.keys() | generated.keys() if p.endswith(".py")}:
        raise ValueError("Unpinned upstream Python source exists")
    pins.python_trees[str(root)] = sorted(actual_python)
    pins.upstream_commits[str(root)] = commit


def checked_assets(assets, base, pins):
    result = {}
    for name, spec in assets.items():
        if "directory" not in spec:
            result[name] = pins.spec(spec, base)
            continue
        if set(spec) != {"directory", "files_sha256"}:
            raise ValueError("Directory asset requires an exact file inventory")
        directory = (Path(base) / spec["directory"]).resolve()
        paths = {str(p.relative_to(directory)): p for p in directory.rglob("*") if p.is_file()}
        if not paths or set(paths) != set(spec["files_sha256"]):
            raise ValueError(f"Directory asset membership differs: {directory}")
        for relative, path in paths.items():
            pins.file(path, spec["files_sha256"][relative])
        pins.directories[str(directory)] = sorted(paths)
        result[name] = {"directory": str(directory), "files_sha256": spec["files_sha256"]}
    return result


def source_snapshot(root, commit, pins):
    root = Path(root).resolve()
    actual = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    if actual != commit:
        raise ValueError("Tested runner checkout is not at its declared commit")
    paths = sorted(
        {
            *root.glob("*.py"),
            *(root / "benchmarks/research").glob("*.py"),
            *(root / "chembounce_synthesis").glob("*.py"),
        }
    )
    source = {}
    for path in paths:
        relative = str(path.relative_to(root))
        committed = subprocess.check_output(["git", "show", f"{commit}:{relative}"], cwd=root)
        source[relative] = pins.file(path, hashlib.sha256(committed).hexdigest())["sha256"]
    if len(source) != 41:
        raise ValueError("Expected exactly 41 tested scientific source files")
    return {"root": str(root), "commit": commit, "files_sha256": source}


def runtime_snapshot(python, pins, rdkit_version="2025.9.6", import_crem_from=None):
    # Keep the venv spelling: resolving its interpreter symlink loses isolation.
    python = Path(python).absolute()
    binary = pins.file(python)
    config = pins.file(python.parent.parent / "pyvenv.cfg")
    if "include-system-site-packages = false" not in Path(config["artifact"]).read_text().lower():
        raise ValueError("Runtime must exclude system site packages")
    if not os.access(python, os.X_OK):
        raise ValueError("Runtime Python is not executable")
    code = (
        "import importlib.metadata as m,json,sys\n"
        "if len(sys.argv)>1:\n"
        " import crem\n"
        " from pathlib import Path\n"
        " assert Path(crem.__file__).resolve().is_relative_to(Path(sys.argv[1]).resolve())\n"
        "print(json.dumps({'prefix':sys.prefix,'base_prefix':sys.base_prefix,"
        "'version':sys.version,'packages':dict(sorted((d.metadata['Name'].lower().replace('_','-'),"
        "d.version) for d in m.distributions()))}))"
    )
    report = json.loads(
        subprocess.check_output(
            [str(python), "-I", "-c", code, *([str(import_crem_from)] if import_crem_from else [])],
            text=True,
            timeout=30,
            env={**os.environ, **THREAD_ENV, "CUDA_VISIBLE_DEVICES": ""},
        )
    )
    if Path(report["prefix"]).resolve() != python.parent.parent.resolve():
        raise ValueError("Runtime is not the declared isolated environment")
    if report["prefix"] == report["base_prefix"]:
        raise ValueError("A dedicated virtual environment is required")
    if report["packages"].get("rdkit", "").replace(".09.", ".9.") != rdkit_version:
        raise ValueError(f"Runtime requires declared RDKit {rdkit_version}")
    return {
        "python": binary,
        "venv_config": config,
        "rdkit_policy": rdkit_version,
        "metadata_probe_crem_root": str(import_crem_from) if import_crem_from else None,
        **report,
        "scope": "Interpreter/config hashes and distribution versions; not hashes of all installed binaries",
    }


def bind_runtime(
    python, pins, *, versions=None, binary_sha=None, rdkit_version="2025.9.6", import_crem_from=None
):
    python = str(Path(python).absolute())
    if python not in pins.runtimes:
        pins.runtimes[python] = runtime_snapshot(python, pins, rdkit_version, import_crem_from)
    result = pins.runtimes[python]
    if result["rdkit_policy"] != rdkit_version:
        raise ValueError("Runtime reused with a conflicting RDKit policy")
    if result["metadata_probe_crem_root"] != (str(import_crem_from) if import_crem_from else None):
        raise ValueError("Runtime reused with a conflicting import-only metadata probe")
    if versions is not None:
        normalized = {k.lower().replace("_", "-"): v for k, v in versions.items()}
        if len(normalized) != len(versions) or result["packages"] != normalized:
            raise ValueError("Runtime distribution map differs from its prepared bundle/manifest")
    if binary_sha is not None:
        pins.file(python, binary_sha)
    return result


def validate_output_location(output, decision, amendment):
    output = Path(output).resolve()
    root = Path(decision["source_root"]).resolve()
    permitted = root / "benchmark-results/research-validation"
    if output.is_relative_to(root) and not output.is_relative_to(permitted):
        raise ValueError(
            "Root campaign writes are restricted to its ignored research-validation result subtree"
        )
    protected = [
        Path(p).resolve() for p in amendment["execution_policy"]["protected_resources"] if p.startswith("/")
    ]
    protected.append(Path(amendment["current_pipeline"]["source_worktree"]).resolve())
    if any(output.is_relative_to(p) for p in protected):
        raise ValueError("Output would modify a protected worktree")
    return output


def public_queries(panel, split):
    compounds = {r["id"]: r for r in panel["compounds"]}
    if len(compounds) != len(panel["compounds"]):
        raise ValueError("Duplicate compound IDs")
    queries = []
    for query in panel["queries"]:
        if set(query) - {"id", "compound_id"}:
            raise ValueError("This ligand-only follow-up permits no query reference files or extra fields")
        record = compounds[query["compound_id"]]
        queries.append({"id": query["id"], "smiles": record["smiles"], "target_id": record["target_id"]})
    if len(queries) != COUNTS[split] or len({q["id"] for q in queries}) != len(queries):
        raise ValueError("Panel query grid differs")
    return queries


def external_method(spec, base, method_id, split, decision, panel_sha, pins):
    fragment = pins.spec(spec, base)
    path = Path(fragment["artifact"])
    method = read_json(path)
    if set(method) - {"id", "kind", "command", "settings", "budget_overrides", "assets", "source_binding"}:
        raise ValueError("Private or unknown method fragment fields")
    if method.get("id") != method_id or method.get("kind") != "external":
        raise ValueError("External method identity differs")
    if method.get("budget_overrides") != {"max_evaluations": 256}:
        raise ValueError("External native work cap must be 256")
    command = method.get("command", [])
    suffix = ["sample"] if method_id == METHODS[1] else []
    if len(command) != 6 + len(suffix) or command[2:] != suffix + [
        "--request",
        "{request}",
        "--response",
        "{response}",
    ]:
        raise ValueError("External command must be the declared direct adapter invocation")
    if not all(Path(p).is_absolute() for p in command[:2]):
        raise ValueError("External interpreter and adapter must use absolute paths")
    assets = method.get("assets", {})
    if method_id == METHODS[1]:
        required = {"crem_db", "crem_manifest"}
        sources = {"source__" + p.replace("/", "__").replace(".", "_"): p for p in CREM_FILES}
        source_root = Path(command[1]).parents[2]
    else:
        required, sources = {"mol2mol_bundle", "mol2mol_model"}, MOL2MOL_FILES
        source_root = Path(command[1]).parents[3]
    if set(assets) != required | set(sources):
        raise ValueError("External assets are restricted to declared model/DB and source files")
    if any("directory" in s for s in assets.values()):
        raise ValueError("External fragment assets must be explicit files")
    method["assets"] = checked_assets(assets, path.parent, pins)
    for name, relative in sources.items():
        if Path(method["assets"][name]["artifact"]).resolve() != (source_root / relative).resolve():
            raise ValueError("External source asset path does not match its allowlist")
    if str(Path(command[1]).resolve()) not in {
        str(Path(s["artifact"]).resolve()) for s in method["assets"].values()
    }:
        raise ValueError("Adapter command source is not hash-bound as an asset")
    if method_id == METHODS[1]:
        manifest = read_json(method["assets"]["crem_manifest"]["artifact"])
        if method["settings"] != {
            "policy": "crem-train-only-partial-exo-v1",
            "manifest_sha256": method["assets"]["crem_manifest"]["sha256"],
        }:
            raise ValueError("CReM settings or public manifest hash differ")
        if manifest.get("complete") is not True or manifest.get("labels_exposed_to_sampler") is not False:
            raise ValueError("CReM DB is not a completed public sampler asset")
        if manifest.get("input_hashes", {}).get("dataset") != panel_sha:
            raise ValueError(f"CReM database does not belong to the {split} panel")
        if manifest.get("inventory_counts", {}).get("nontrain_partition_family_overlap") != 0:
            raise ValueError("CReM DB includes held-out families")
        if manifest.get("db_sha256") != method["assets"]["crem_db"]["sha256"]:
            raise ValueError("CReM DB/public manifest mismatch")
        mutation = copy.deepcopy(decision["crem"]["mutation_proposal"])
        mutation["min_rel_size"] = mutation.pop("min_rel")
        mutation["max_rel_size"] = mutation.pop("max_rel")
        mutation.update(
            radius=3, symmetry_fixes=False, return_rxn=False, return_rxn_freq=False, return_mol=False
        )
        if manifest.get("mutation_settings") != mutation:
            raise ValueError("CReM mutation settings were changed")
        build = manifest.get("build_settings", {})
        required_build = {
            "set_name": "train",
            "radii": [3],
            "ncpu": 1,
            "max_heavy_atoms": 15,
            "keep_stereo": True,
            "frag_mode": "both_optimal",
        }
        if any(build.get(k) != v for k, v in required_build.items()):
            raise ValueError("CReM train-only DB construction settings differ")
        if manifest.get("upstream", {}).get("commit") != decision["crem"]["upstream_commit"]:
            raise ValueError("CReM upstream source differs")
        upstream = manifest["upstream"]
        pin_upstream(
            upstream["root"],
            upstream["commit"],
            upstream["tracked_files_sha256"],
            pins,
            generated=upstream["generated_files_sha256"],
            all_tracked=True,
        )
        implementation = manifest["implementation"]
        if Path(implementation["root"]).resolve() != source_root.resolve() or set(
            implementation["files_sha256"]
        ) != set(CREM_FILES):
            raise ValueError("CReM implementation inventory differs from exact source allowlist")
        pin_inventory(implementation["root"], implementation["files_sha256"], pins)
        prepared_runtime = manifest["runtime"]
        if Path(prepared_runtime["executable"]).absolute() != Path(command[0]).absolute():
            raise ValueError("CReM command runtime differs from its build environment")
        runtime = bind_runtime(
            command[0], pins, versions=prepared_runtime["distributions"], import_crem_from=upstream["root"]
        )
        if (
            runtime["version"] != prepared_runtime["python"]
            or runtime["prefix"] != prepared_runtime["prefix"]
        ):
            raise ValueError("CReM prepared interpreter version/prefix changed")
        for origin in prepared_runtime["origins"].values():
            pins.file(origin["path"], origin["sha256"])
            if Path(origin["path"]).stat().st_size != origin["bytes"]:
                raise ValueError("CReM prepared import-origin size changed")
    else:
        if method["settings"] != MOL2MOL_SETTINGS:
            raise ValueError("Mol2Mol declared sampling/selection settings were changed")
        bundle = read_json(method["assets"]["mol2mol_bundle"]["artifact"])
        if bundle.get("upstream_commit") != decision["mol2mol"]["upstream_commit"]:
            raise ValueError("Mol2Mol upstream commit differs")
        if (
            bundle.get("model_sha256") != method["assets"]["mol2mol_model"]["sha256"]
            or bundle.get("model_sha512") != decision["mol2mol"]["official_sha512"]
        ):
            raise ValueError("Mol2Mol official checkpoint differs")
        if (
            digest(method["assets"]["mol2mol_model"]["artifact"], "sha512")
            != decision["mol2mol"]["official_sha512"]
        ):
            raise ValueError("Mol2Mol checkpoint bytes do not match official SHA-512")
        if Path(bundle.get("common_python", "")).absolute() != Path(command[0]).absolute():
            raise ValueError("Mol2Mol common interpreter differs from prepared bundle")
        pin_upstream(bundle["upstream"], bundle["upstream_commit"], bundle["source_sha256"], pins)
        if (
            Path(bundle["adapter_root"]).resolve() != Path(command[1]).resolve().parent
            or Path(bundle["common_root"]).resolve() != source_root.resolve()
        ):
            raise ValueError("Mol2Mol prepared project roots differ")
        pin_inventory(bundle["adapter_root"], bundle["adapter_sha256"], pins)
        pin_inventory(bundle["common_root"], bundle["common_sha256"], pins)
        bind_runtime(
            command[0], pins, versions=bundle["common_versions"], binary_sha=bundle["common_python_sha256"]
        )
        bind_runtime(
            bundle["sampler_python"],
            pins,
            versions=bundle["sampler_versions"],
            binary_sha=bundle["sampler_python_sha256"],
            rdkit_version="2023.9.5",
        )
    return method, fragment


def analysis_plan(grid, source):
    return {
        "schema_version": 1,
        "status": "specified_not_executed",
        "expected_attempts": 270,
        "design": "Post hoc supplementary comparison; not preregistered; original 630 remain separate",
        "split_attempts": {s: COUNTS[s] * 9 for s in SPLITS},
        "expected_grid": grid,
        "gate": "All 270 unique planned slots, completed common scoring, response/source hashes and candidate occurrence joins verified before comparative estimates",
        "primary": "Engineering delivery endpoint: attempt_with_hop = int(evaluation.metrics.hop_families > 0). This measures delivery of at least one structurally changed output, not biological or synthetic utility",
        "evaluator_metrics": [
            "run_success",
            "candidate_count",
            "valid_candidate_count",
            "valid_fraction",
            "unique_families",
            "duplicate_chemical_states",
            "unchanged_query_families",
            "hop_families",
            "unique_hop_scaffolds",
            "unique_hop_ring_skeletons",
            "reference_active_scaffolds",
            "recovered_active_scaffolds",
            "active_scaffold_recall_proxy",
            "active_scaffold_precision_proxy",
            "exact_known_active_matches",
            "known_inactive_family_matches",
            "decoy_family_matches",
            "verified_rule_routes",
        ],
        "structural_definitions": {
            "identity": "Common RDKit 2025.09.6 protonation-tautomer-family-v1 policy; family ignores stereo/isotopes, exact canonical isomeric identity retains them",
            "hop": "Candidate has a ring scaffold, differs from input chemical family, and differs from normalized Murcko scaffold family; generic ring-skeleton change is separately counted",
            "deduplication": "First accepted valid occurrence per chemical family, preserving returned rank; invalid/duplicate raw-output accounting is retained separately",
            "reference_active_scaffolds": "Unique normalized active scaffold families among same-target test partition, excluding input family and input scaffold and requiring a ring scaffold",
            "recovered_active_scaffolds": "Distinct returned hop scaffold families intersecting reference_active_scaffolds",
            "active_scaffold_recall_proxy": "recovered_active_scaffolds / reference_active_scaffolds; null when reference set is empty. No Tanimoto >= 0.5 gate or other similarity threshold is applied",
            "active_scaffold_precision_proxy": "Number of unique returned hop families whose scaffold matches an active reference / hop_families; zero if no hops",
            "exact_known_active_matches": "Unique returned families represented by first-ranked valid output whose exact canonical isomeric SMILES matches an active same-target test record",
            "known_inactive_family_matches": "Unique returned families present in inactive same-target test records; this is family-level, not exact inactive identity",
            "interpretation": "Scaffold/family recovery is a retrospective structural proxy, not measured or predicted activity of new products; earlier control prefilter recall at T >= 0.5 is a different endpoint",
        },
        "similarity": {
            "candidate_fields": ["hit_bit_tanimoto", "hit_count_tanimoto"],
            "bit": "Morgan radius 2, 2048 bits, includeChirality=True",
            "count": "Count Morgan radius 2, includeChirality=False",
            "common_3d": "Fixed common CONFORMERS and seeds 42/43/44 in common-3d-plan.json; ligand-only exact canonical query/product inputs; all-seed success required, median over seeds for hit_shape/hit_pharmacophore/score and score_stability summaries",
            "aggregation": "Attempt mean of valid candidate occurrence similarities, then seed/query/target hierarchy; no-output/failed-score attempts are null for conditional means and remain in coverage denominators",
            "paired_missing": "For conditional score contrasts, pair only identical query/seed slots with observed scores for both methods; report excluded slot/query/target counts. Bounded output indicators retain every planned slot",
            "score_scope": "The common 3D score retains the pinned conformer energy/strain term; it is not a calibrated probability or pure shape score",
        },
        "secondary": [
            "distinct families and chemical/generic scaffolds",
            "exact known active and known inactive family recovery",
            "held-out active-scaffold recovery proxy",
            "input bit/count fingerprint similarity",
            "common 3D success coverage and conditional shape/pharmacophore/joint scores",
            "route-supplied/replay/forward coverage",
            "native work counts and measured resources",
        ],
        "hierarchy": "Within each split: seed mean per query, query mean per target, equal target macro; paired method differences at identical query/seed slots before aggregation",
        "comparisons": [
            "pipeline_current minus crem_train_partial_exo",
            "pipeline_current minus mol2mol_pretrained",
        ],
        "bootstrap": {
            "unit": "target",
            "draws": 10000,
            "seed": 20261016,
            "quantiles": [0.025, 0.975],
            "single_target": "no CI",
            "multiplicity": "descriptive unadjusted intervals; no confirmatory significance claim",
        },
        "failures": "Keep every planned attempted outcome, including timeout, invalid, no_candidates and infrastructure interruption; terminal failures accept zero partial candidates; no retry/backfill",
        "infrastructure": "Primary realized accepted return is zero; separately show bounded [0,1] missing-counterfactual sensitivity, not CI, and paired complete-case sensitivity excluding infrastructure only",
        "missing_scores": "Failed/unassessed 3D or route evidence is coverage loss, never a fabricated similarity zero or chemical synthesis failure; conditional means retain observed counts",
        "route_policy": "No route supplied is unassessed; a CReM graph-edit reaction is not a synthetic route. Forward support uses the same predictor as pipeline selection, with selection bias; not independent experimental validation",
        "dependence": "Repeated seeds, molecules, routes and cached pairs are not independent samples; report distinct identities/keys and occurrence denominators",
        "reuse": {
            "source": source,
            "functions": [
                "benchmarks.research.evaluation.evaluate_candidates",
                "benchmarks.research.controls.similarity_scores",
                "benchmarks.research.evaluation.cluster_interval",
            ],
            "note": "Preparation implements no scientific metrics; final paired analysis adapter still requires independent assessment",
        },
        "limits": [
            "Pretrained chemical inventory incomplete/unknown",
            "Methods use different training libraries",
            "No activity, synthesis yield/probability or journal-readiness claim",
            "Shared host is not a controlled hardware speed benchmark",
        ],
    }


def prepare(decision_path, amendment_path, fragments_path, output):
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise ValueError("Output must be new; immutable plans cannot be overwritten or resumed")
    run_root = output.with_name(output.name + "-runs")
    if run_root.exists() or run_root.is_symlink():
        raise ValueError("Planned run output must be absent; no reuse or backfill")
    pins = Pins()
    pins.file(__file__)
    decision_spec = pins.file(decision_path, DECISION_SHA)
    amendment_spec = pins.file(amendment_path, AMENDMENT_SHA)
    decision, amendment = read_json(decision_path), read_json(amendment_path)
    output = validate_output_location(output, decision, amendment)
    run_root = validate_output_location(run_root, decision, amendment)
    base = Path(decision_path).resolve().parent
    if amendment["previous_decision"]["sha256"] != DECISION_SHA:
        raise ValueError("Amendment does not bind the parent decision")
    if (
        amendment["supplementary_grid"]["methods"] != list(METHODS)
        or amendment["supplementary_grid"]["seeds"] != SEEDS
    ):
        raise ValueError("Amendment grid differs")
    fragments_spec = pins.file(fragments_path)
    fragments = read_json(fragments_path)
    if set(fragments) != {"crem_train_partial_exo", "mol2mol_pretrained"} or set(
        fragments[METHODS[1]]
    ) != set(SPLITS):
        raise ValueError("Supply three completed CReM fragments and one Mol2Mol fragment")
    current = amendment["current_pipeline"]
    source = source_snapshot(current["source_worktree"], current["source_commit"], pins)
    python = current["runtime"]
    bind_runtime(python, pins)
    documents, grid, commands, protocols = {}, [], [], {}
    for split in SPLITS:
        panel_info = decision["panels"][split]
        panel = pins.file(base / panel_info["path"], panel_info["sha256"])
        panel_data = read_json(panel["artifact"])
        queries = public_queries(panel_data, split)
        documents[f"public-queries/{split}.json"] = queries
        original_name, original_sha = ORIGINAL_PROTOCOLS[split]
        original_spec = pins.file(base / original_name, original_sha)
        original = read_json(original_spec["artifact"])
        if original["seeds"] != SEEDS:
            raise ValueError("Original pipeline seeds differ")
        matches = [m for m in original["methods"] if m["id"] == "pipeline"]
        if (
            len(matches) != 1
            or matches[0].get("kind") != "pipeline"
            or matches[0].get("budget_overrides") != {"max_evaluations": 12}
        ):
            raise ValueError("Original pipeline method/evaluation cap differs")
        pipeline = copy.deepcopy(matches[0])
        pipeline["id"] = METHODS[0]
        pipeline["assets"]["fingerprint_db"] = {
            "directory": current["fingerprint_db"]["path"],
            "files_sha256": current["fingerprint_db"]["files_sha256"],
        }
        pipeline["assets"] = checked_assets(pipeline["assets"], base, pins)
        for inventory in pipeline.get("inventories", []):
            if "artifact" in inventory:
                inventory.update(pins.spec({k: inventory[k] for k in ("artifact", "sha256")}, base))
        for method_id in METHODS:
            fragment = original_spec
            if method_id == METHODS[0]:
                method = pipeline
            else:
                spec = fragments[method_id][split] if method_id == METHODS[1] else fragments[method_id]
                method, fragment = external_method(
                    spec, Path(fragments_path).parent, method_id, split, decision, panel["sha256"], pins
                )
            name = f"protocols/{split}-{method_id}.json"
            protocol = {
                "schema_version": 1,
                "budget": BUDGET,
                "seeds": SEEDS,
                "methods": [method],
                "search_space_audit_max_records": original.get("search_space_audit_max_records", 5000),
                "campaign": {
                    "id": "comparator-followup-270",
                    "split": split,
                    "scope": "Post hoc supplementary comparison, separate from original 630",
                    "phase": "prepared_requires_parent_activation",
                    "parent_fragment": fragment,
                    "native_work_unit": "pipeline proposal evaluations"
                    if method_id == METHODS[0]
                    else (
                        "CReM API-yielded products" if method_id == METHODS[1] else "Mol2Mol raw sampled rows"
                    ),
                },
            }
            documents[name] = protocol
            protocols[f"{split}/{method_id}"] = name
            run_output = run_root / split / method_id
            gpu = method_id != METHODS[1]
            commands.append(
                {
                    "split": split,
                    "method_id": method_id,
                    "expected_attempts": len(queries) * 3,
                    "cwd": source["root"],
                    "argv": [
                        python,
                        "-m",
                        "benchmarks.research.runner",
                        "--dataset",
                        panel["artifact"],
                        "--protocol",
                        str(output / name),
                        "--output",
                        str(run_output),
                    ],
                    "environment_overrides": {
                        **THREAD_ENV,
                        "PYTHONPATH": source["root"],
                        "CUDA_VISIBLE_DEVICES": amendment["execution_policy"]["single_gpu"] if gpu else "",
                    },
                    "resource_class": "exclusive_physical_gpu0" if gpu else "cpu_single_thread",
                    "execution_order": len(commands),
                    "no_retry_or_resume": True,
                }
            )
            grid.extend(
                {
                    "split": split,
                    "method_id": method_id,
                    "query_id": q["id"],
                    "target_id": q["target_id"],
                    "query_index": index,
                    "seed": seed,
                    "output_dir": str(run_output / method_id / f"query-{index:04d}-seed-{seed}"),
                }
                for index, q in enumerate(queries)
                for seed in SEEDS
            )
    if len(grid) != 270 or len({(g["split"], g["method_id"], g["query_id"], g["seed"]) for g in grid}) != 270:
        raise ValueError("The exact 270-slot grid is required")
    documents["grid.json"] = grid
    documents["analysis-plan.json"] = analysis_plan(grid, source)
    documents["common-3d-plan.json"] = {
        "status": "specified_not_executed",
        "function": "benchmarks.research.controls.similarity_scores",
        "source": source,
        "conformers": CONFORMERS,
        "seeds": SEEDS,
        "resource_scope": "Separate post-generation pass outside each method generation budget",
        "occurrences": "One assessment record per accepted candidate rank/SMILES and attempt; retain all failure records",
        "cache": "Exact query/candidate pair plus source/config/seed key; reuse is not independent evidence",
        "coverage": "Report all accepted candidate occurrences, successful 3D, failed 3D and per-attempt output-slot coverage; means conditional on all scoring seeds succeeding",
        "basic_runner_computes_3d": False,
    }
    documents["execution-plan.json"] = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "prepared_requires_parent_activation",
        "execution_authorized": False,
        "parent_decision": decision_spec,
        "amendment": amendment_spec,
        "fragments": fragments_spec,
        "expected_attempts": 270,
        "protocols": protocols,
        "source": source,
        "runtimes": pins.runtimes,
        "commands": commands,
        "scheduling": "All nine commands sequential; never overlap GPU methods; GPU0 UUID only",
        "activation_requirements": amendment["required_before_launch"],
        "privacy": "Runner is trusted evaluator with panel access; its public_request whitelist passes id/smiles/target_id only. This is information-flow discipline, not OS filesystem isolation",
        "native_work_counts": "12 pipeline proposal evaluations vs 256 API products or model samples; units are not universally comparable",
        "input_files_sha256": pins.files,
        "input_directory_membership": pins.directories,
        "input_python_membership": pins.python_trees,
        "upstream_commits": pins.upstream_commits,
    }
    # Missing final assets and every validation error above leave no executable plan behind.
    pins.recheck()
    output.mkdir(parents=True, exist_ok=False)
    for relative, value in documents.items():
        path = output / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x") as stream:
            stream.write(encode(value))
    manifest = {
        "schema_version": 1,
        "status": "prepared_requires_parent_activation",
        "execution_authorized": False,
        "files_sha256": {p: digest(output / p) for p in documents},
    }
    with (output / "manifest.json").open("x") as stream:
        stream.write(encode(manifest))
    return manifest


def verify_prepared(output, manifest_sha256):
    output = Path(output)
    if digest(output / "manifest.json") != manifest_sha256:
        raise ValueError("Prepared manifest differs from its externally retained SHA-256")
    manifest = read_json(output / "manifest.json")
    expected = set(manifest["files_sha256"]) | {"manifest.json"}
    actual = {str(p.relative_to(output)) for p in output.rglob("*") if p.is_file()}
    if actual != expected:
        raise ValueError("Prepared directory contents changed; use a separate location for run artifacts")
    for relative, sha in manifest["files_sha256"].items():
        if digest(output / relative) != sha:
            raise ValueError(f"Prepared output was modified: {relative}")
    plan = read_json(output / "execution-plan.json")
    for path, sha in plan["input_files_sha256"].items():
        if digest(path) != sha:
            raise ValueError(f"Prepared input changed: {path}")
    for directory, expected in plan["input_directory_membership"].items():
        actual = {str(p.relative_to(directory)) for p in Path(directory).rglob("*") if p.is_file()}
        if actual != set(expected):
            raise ValueError(f"Prepared directory asset changed: {directory}")
    for directory, expected in plan["input_python_membership"].items():
        actual = {str(p.relative_to(directory)) for p in Path(directory).rglob("*.py") if p.is_file()}
        if actual != set(expected):
            raise ValueError(f"Prepared upstream Python membership changed: {directory}")
    for root, commit in plan["upstream_commits"].items():
        check_commit(root, commit)
    source = plan["source"]
    if source_snapshot(source["root"], source["commit"], Pins()) != source:
        raise ValueError("Prepared scientific source snapshot differs")
    for python, expected in plan["runtimes"].items():
        if (
            runtime_snapshot(python, Pins(), expected["rdkit_policy"], expected["metadata_probe_crem_root"])
            != expected
        ):
            raise ValueError("Prepared runtime distribution inventory differs")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("prepare")
    for name in ("decision", "amendment", "fragments", "output"):
        command.add_argument("--" + name, type=Path, required=True)
    check = sub.add_parser("verify")
    check.add_argument("--output", type=Path, required=True)
    check.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args(argv)
    if args.command == "prepare":
        result = prepare(args.decision, args.amendment, args.fragments, args.output)
    else:
        result = verify_prepared(args.output, args.manifest_sha256)
    print(
        encode(
            {
                "status": result["status"],
                "execution_authorized": False,
                "manifest_sha256": digest(args.output / "manifest.json"),
            }
        ),
        end="",
    )


if __name__ == "__main__":
    main()
