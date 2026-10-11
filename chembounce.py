#!/usr/bin/env python3
"""ChemBounce 2.0: reproducible and diagnosable scaffold hopping."""

import csv
import hashlib
import importlib.metadata
import json
import logging
import math
from pathlib import Path
import platform
import subprocess
from collections import Counter
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator

import chemistry
import cli
import cost_estimation
import scaffoldgraph_fragmenter as sg
import utils

__version__ = "2.0.0.dev0"
LOGGER = logging.getLogger(__name__)
RESULT_FEATURES = [
    "Fragment_no",
    "Original scaffold",
    "Replaced scaffold",
    "Final structure",
    "Standardized final structure",
    "Tanimoto Similarity",
    "Electron shape Similarity",
    "CATS2D dist",
    "QED",
    "SAscore",
    "logP",
    "MW",
    "H_Donors",
    "H_Acceptors",
]


def remove_substructures(mol, pattern_mol):
    return [Chem.MolToSmiles(part) for part in chemistry.cut_scaffold(mol, pattern_mol)]


def replace_molecule(
    target_mol, pattern_mol, replace_mol, max_n=10, _merge_structure_top_n_=100, *, prepared=None
):
    return chemistry.replace(
        target_mol, pattern_mol, replace_mol, max_n, _merge_structure_top_n_, prepared=prepared
    )


def _iter_library(filename, diagnostics=None):
    stats = {"path": str(filename), "valid_records": 0, "rejected_records": 0, "examples": []}
    with open(filename) as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            try:
                with rdBase.BlockLogs():
                    mol, _ = utils.init_mol(line.strip())
                if len(Chem.GetMolFrags(mol)) != 1:
                    raise ValueError("Multiple components")
            except (utils.SmilesInputError, ValueError) as exc:
                stats["rejected_records"] += 1
                if len(stats["examples"]) < 10:
                    stats["examples"].append({"line": number, "smiles": line.strip(), "reason": str(exc)})
                continue
            stats["valid_records"] += 1
            yield mol
    if diagnostics is not None:
        diagnostics["text_library"] = stats


def read_fragments(target_smiles, fragment_file, diagnostics=None):
    """Read valid connected scaffolds and report unusable library records."""
    return list(_iter_library(fragment_file, diagnostics))


def search_similar_scaffolds(
    original_scaffold,
    fragments_DB,
    low_mem=False,
    tqdm_quiet=False,
    scaffold_top_n=None,
    threshold=0.3,
    radius=2,
    fp_size=2048,
):
    """Streaming reference implementation using the same bits as the fast path."""
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size)
    query = generator.GetFingerprint(original_scaffold)
    original = Chem.MolToSmiles(original_scaffold)
    scores = {}
    for entry in fragments_DB:
        mol, smiles = utils.init_mol(entry) if isinstance(entry, str) else (entry, Chem.MolToSmiles(entry))
        if mol is None or smiles == original:
            continue
        score = DataStructs.TanimotoSimilarity(query, generator.GetFingerprint(mol))
        if score >= threshold:
            scores[smiles] = score
            if scaffold_top_n and len(scores) > 2 * scaffold_top_n:
                scores = dict(sorted(scores.items(), key=lambda item: (-item[1], item[0]))[:scaffold_top_n])
    # Canonical SMILES tie order is shared by freshly generated databases.
    ranked = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    if scaffold_top_n is not None:
        ranked = ranked[:scaffold_top_n]
    return pd.Series(dict(ranked), name="Tanimoto Similarity", dtype=float)


def _read_replacements(filename):
    """Accept one SMILES per line or exported SMILES/score TSV with a header."""
    with Path(filename).open() as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    parsed = {}
    for number, row in enumerate(rows, 1):
        if not row or not row[0].strip():
            continue
        if number == 1 and row[0].strip().lower() in ("smiles", "scaffold", "unnamed: 0"):
            continue
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(row[0])
        if mol is None:
            raise ValueError(f"Invalid replacement SMILES at {filename}:{number}: {row[0]}")
        mol, smiles = utils.init_mol(row[0])
        if len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError(f"Multicomponent replacement at {filename}:{number}")
        score = float(row[1]) if len(row) > 1 else -float(number)
        if not math.isfinite(score):
            raise ValueError(f"Non-finite replacement priority at {filename}:{number}")
        parsed[smiles] = max(score, parsed.get(smiles, -math.inf))
    return pd.Series(dict(sorted(parsed.items(), key=lambda x: (-x[1], x[0]))), name="Priority", dtype=float)


def _json_write(path, data):
    """Atomically replace a JSON report, rejecting nonstandard NaN/Infinity."""
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temp.replace(path)


def _hash_file(path):
    path = Path(path)
    if path.is_dir():
        return {p.name: _hash_file(p) for p in sorted(path.iterdir()) if p.is_file()}
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# Snapshot source files when this module is loaded, so later edits to an active
# notebook/server checkout cannot silently replace the recorded implementation.
def _source_hashes():
    root = Path(__file__).parent
    paths = sorted(root.glob("*.py")) + sorted((root / "chembounce_synthesis").glob("*.py"))
    return {str(p.relative_to(root)): _hash_file(p) for p in paths}


_SOURCE_HASHES = _source_hashes()


def _changed_sources():
    current = _source_hashes()
    return sorted(
        name
        for name in current.keys() | _SOURCE_HASHES.keys()
        if current.get(name) != _SOURCE_HASHES.get(name)
    )


def _provenance(parameters, sources):
    versions = {}
    for name in ("rdkit", "numpy", "pandas", "scipy", "scaffoldgraph", "oddt"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    try:
        root = Path(__file__).parent
        revision = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True, stderr=subprocess.DEVNULL
        ).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=root, text=True))
    except (OSError, subprocess.CalledProcessError):
        revision, dirty = None, None
    return {
        "schema_version": 1,
        "chembounce_version": __version__,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": revision,
        "git_dirty": dirty,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "dependencies": versions,
        "parameters": parameters,
        "sources": {str(Path(s).resolve()): _hash_file(s) for s in sources},
        "source_code_sha256": _SOURCE_HASHES,
        "source_files_changed_after_import": _changed_sources(),
        "methods": {
            "scaffold_search": "Morgan binary radius=2 bits=2048 unless database metadata overrides",
            "candidate_similarity": "Morgan sparse counts radius=2; chirality disabled",
            "replacement": "first substructure match, deterministic bounded beam",
            "shape": "ETKDGv3 seed; UFF maxIters=1000; heavy-atom ElectroShape with Gasteiger charges",
            "standardization": "RDKit Cleanup; input protonation retained during search",
        },
    }


def get_frags_cands(
    target_smiles,
    target_mol,
    fragments_DB=None,
    overall_max_n=None,
    frag_max_n=None,
    scaffold_top_n=None,
    cand_max_n__rplc=None,
    _merge_structure_top_n_=100,
    murcko_frag_itr_rnd=1,
    _search_scf_thr_=0.3,
    fragments=None,
    replace_scaffold_files=None,
    output_dir="./output",
    low_mem=False,
    tqdm_quiet=False,
    use_fingerprint_db=True,
    scaffold_db_file=None,
    fingerprint_db_file=None,
    trust_legacy_db=False,
    diagnostics=None,
):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    fragments, files = fragments or [], replace_scaffold_files or []
    if isinstance(files, (str, Path)):
        files = [files]
    frags = (
        [utils.init_mol(s)[0] for s in fragments]
        if fragments
        else sg.get_all_murcko_fragments(
            target_mol, break_fused_rings=False, iteration_round=murcko_frag_itr_rnd
        )
    )
    for frag in frags:
        if not target_mol.HasSubstructMatch(frag, useChirality=True):
            raise ValueError(f"Fragment is not contained in input: {Chem.MolToSmiles(frag)}")
    if files and len(files) not in (1, len(frags)):
        raise ValueError("Supply one replacement file for all fragments or one per fragment")
    pd.DataFrame(
        {"SMILES": [Chem.MolToSmiles(f) for f in frags]},
        index=pd.Index(range(1, len(frags) + 1), name="Fragment_no"),
    ).to_csv(out / "fragment_info.tsv", sep="\t")
    limits = utils._scaffold_no_reassign_(
        overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_, len(frags)
    )
    top_n = limits[2]
    fp_db = None
    # Explicit user libraries take precedence over the bundled fingerprint DB.
    if (
        frags
        and not files
        and use_fingerprint_db
        and (fingerprint_db_file or (scaffold_db_file is None and fragments_DB is None))
    ):
        try:
            fp_db = utils.load_fingerprint_db(fingerprint_db_file, trust_legacy=trust_legacy_db)
        except FileNotFoundError:
            if fingerprint_db_file:
                raise
            LOGGER.warning("Default fingerprint database missing; using scaffold text database")
    if fp_db is not None and scaffold_db_file is not None:
        recorded_hash = fp_db.get("source_sha256")
        if recorded_hash is None:
            LOGGER.warning("Legacy DB has no source hash; correspondence to --scaffold-db cannot be verified")
        elif str(recorded_hash) != _hash_file(scaffold_db_file):
            raise ValueError("Fingerprint database source hash does not match --scaffold-db")
    source = None
    if frags and not files and fp_db is None and fragments_DB is None:
        source = (
            Path(scaffold_db_file)
            if scaffold_db_file
            else Path(__file__).parent / "data/Scaffolds_processed.txt"
        )
        if not source.is_file():
            raise FileNotFoundError(
                f"Scaffold database missing: {source}. Run bash install.sh or supply --scaffold-db."
            )
        # Stream the text library anew for each query in low-memory mode.
        if not low_mem:
            fragments_DB = read_fragments(target_smiles, source, diagnostics)
    if diagnostics is not None:
        diagnostics["database"] = {
            "mode": "replacement_files" if files else ("fingerprints" if fp_db else "text"),
            "path": fp_db["path"] if fp_db else (str(source) if source else None),
            "radius": fp_db["radius"] if fp_db else 2,
            "fp_size": fp_db["fp_size"] if fp_db else 2048,
        }
    for number, frag in enumerate(frags, 1):
        if files:
            scores = _read_replacements(files[0] if len(files) == 1 else files[number - 1]).iloc[:top_n]
        elif fp_db is not None:
            scores = utils.search_similar_scaffolds_fast(frag, fp_db, top_n, _search_scf_thr_)
        elif source is not None and low_mem:
            scores = search_similar_scaffolds(
                frag, _iter_library(source, diagnostics), scaffold_top_n=top_n, threshold=_search_scf_thr_
            )
        else:
            scores = search_similar_scaffolds(
                frag, fragments_DB, scaffold_top_n=top_n, threshold=_search_scf_thr_
            )
        scores.rename_axis("SMILES").to_csv(
            out / f"replace_scaffold_list.fragment_{number:04d}.tsv", sep="\t"
        )
    if diagnostics is not None and fp_db is not None:
        diagnostics["database"]["excluded_legacy_records"] = len(fp_db.get("_legacy_rejected_rows", ()))
    return (frags, *limits)


def make_scaffold_hopping(
    target_smiles,
    target_mol,
    frags,
    core_mol=None,
    tanimoto_threshold=0.5,
    overall_max_n=None,
    frag_max_n=None,
    scaffold_top_n=None,
    cand_max_n__rplc=None,
    _merge_structure_top_n_=100,
    output_dir="./output",
    candidate_thresholds=None,
    seed=42,
    calculate_shape=True,
    diagnostics=None,
):
    out = Path(output_dir)
    counters = Counter()
    report = diagnostics if diagnostics is not None else {}
    report["fragments"] = []
    saved = {target_smiles}
    seen = set()
    cats = utils.cats_module.CATS2D()
    target_cats = np.asarray(cats.getCATs2D(target_mol))
    query_fp = chemistry.count_fingerprint(target_mol)
    output_path = out / "overall_result.txt"
    with output_path.open("w", newline="") as overall:
        writer = csv.writer(overall, delimiter="\t")
        writer.writerow(RESULT_FEATURES)
        for number, frag in enumerate(frags, 1):
            if overall_max_n is not None and counters["accepted"] >= overall_max_n:
                break
            local = Counter()
            original = Chem.MolToSmiles(frag)
            prepared = chemistry.prepare_replacement(target_mol, frag)
            replacements = _read_replacements(out / f"replace_scaffold_list.fragment_{number:04d}.tsv")
            with (out / f"fragment_result_{number:04d}.txt").open("w", newline="") as stream:
                frag_writer = csv.writer(stream, delimiter="\t")
                frag_writer.writerow(RESULT_FEATURES)
                for replacement in replacements.index[:scaffold_top_n]:
                    if (frag_max_n is not None and local["accepted"] >= frag_max_n) or (
                        overall_max_n is not None and counters["accepted"] >= overall_max_n
                    ):
                        break
                    local["scaffolds_tested"] += 1
                    replacement_mol, replacement = utils.init_mol(replacement)
                    candidates = replace_molecule(
                        target_mol,
                        frag,
                        replacement_mol,
                        cand_max_n__rplc or 10,
                        _merge_structure_top_n_,
                        prepared=prepared,
                    )
                    if not candidates:
                        local["no_valid_attachments"] += 1
                    for _, candidate in candidates:
                        if (frag_max_n is not None and local["accepted"] >= frag_max_n) or (
                            overall_max_n is not None and counters["accepted"] >= overall_max_n
                        ):
                            break
                        local["generated"] += 1
                        if candidate in saved or candidate in seen:
                            local["duplicate_or_input"] += 1
                            continue
                        seen.add(candidate)
                        mol, candidate = utils.init_mol(candidate)
                        score = DataStructs.TanimotoSimilarity(query_fp, chemistry.count_fingerprint(mol))
                        if score < tanimoto_threshold:
                            local["similarity_rejected"] += 1
                            continue
                        if core_mol is not None and not mol.HasSubstructMatch(core_mol, useChirality=True):
                            local["core_rejected"] += 1
                            continue
                        props = utils._get_molecular_prop_extd_(mol)
                        passed, reasons = utils.determ_props(props, candidate_thresholds or {}, True)
                        if not passed:
                            local["property_rejected"] += 1
                            for reason in reasons:
                                local["filter:" + reason.split(":")[0]] += 1
                            continue
                        standardized, _ = utils.get_standard_smiles(candidate)
                        distance = float(np.linalg.norm(target_cats - np.asarray(cats.getCATs2D(mol))))
                        shape_score = float("nan")
                        if calculate_shape:
                            try:
                                shape_score = utils.calc_electron_shape(target_smiles, candidate, seed)
                            except (ValueError, RuntimeError, ImportError) as exc:
                                local["shape_failed"] += 1
                                report.setdefault("shape_errors", [])
                                if len(report["shape_errors"]) < 10:
                                    report["shape_errors"].append({"smiles": candidate, "reason": str(exc)})
                        row = [
                            number,
                            original,
                            replacement,
                            candidate,
                            standardized,
                            score,
                            shape_score,
                            distance,
                            props["qed"],
                            props["sa"],
                            props["logp"],
                            props["mw"],
                            props["n_hdonor"],
                            props["n_hacceptor"],
                        ]
                        writer.writerow(row)
                        frag_writer.writerow(row)
                        saved.add(candidate)
                        local["accepted"] += 1
                        counters["accepted"] += 1
            counters.update({key: value for key, value in local.items() if key != "accepted"})
            report["fragments"].append(
                {
                    "number": number,
                    "smiles": original,
                    "available_scaffolds": len(replacements),
                    "counts": dict(local),
                }
            )
    report["source_files_changed_after_import"] = _changed_sources()
    report["fragment_count_total"] = len(frags)
    report["fragment_count_screened"] = len(report["fragments"])
    report["global_limit_reached"] = overall_max_n is not None and counters["accepted"] >= overall_max_n
    report["counts"] = dict(counters)
    report["status"] = "success" if counters["accepted"] else "no_candidates"
    if not counters["accepted"]:
        if not frags:
            reason = "No ring scaffolds found. Supply --fragments for an acyclic input."
        elif counters["property_rejected"]:
            reason = "Candidates failed property filters. For nucleotide-like inputs consider --wo_lipinski or explicit limits; filters are never relaxed automatically."
        elif counters["similarity_rejected"]:
            reason = "Candidates failed the similarity threshold. Inspect replacement lists and consider an explicit threshold change."
        elif counters["core_rejected"]:
            reason = "Candidates did not preserve the requested core."
        else:
            reason = "No distinct valid candidates within the search limits. Inspect replacement lists and attachment counts."
        report["message"] = reason
        LOGGER.warning(reason)
    _json_write(out / "run_summary.json", report)
    return str(output_path)


def _validate_parameters(tanimoto_threshold, search_threshold, seed, thresholds):
    for value in (tanimoto_threshold, search_threshold):
        if not isinstance(value, (float, int)) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("Similarity thresholds must be finite numbers between 0 and 1")
    if not isinstance(seed, int) or not 0 <= seed <= 2147483647:
        raise ValueError("seed must be an integer between 0 and 2147483647")
    for key, bounds in thresholds.items():
        if key not in ("sa", "qed", "mw", "logp", "n_hdonor", "n_hacceptor") or len(bounds) != 2:
            raise ValueError(f"Invalid property filter: {key}")
        lo, hi = bounds
        for value in bounds:
            if value is not None and (not isinstance(value, (int, float)) or math.isnan(value)):
                raise ValueError(f"Invalid {key} bound")
        if lo is not None and hi is not None and lo > hi:
            raise ValueError(f"{key} minimum exceeds maximum")


@cost_estimation.EstDecoMem
def chembounce(
    target_smiles,
    fragments_DB=None,
    core_smiles=None,
    tanimoto_threshold=0.5,
    overall_max_n=None,
    frag_max_n=None,
    scaffold_top_n=None,
    cand_max_n__rplc=None,
    _merge_structure_top_n_=100,
    output_dir="./output",
    low_mem=False,
    tqdm_quiet=False,
    fragments=None,
    replace_scaffold_files=None,
    candidate_thresholds=None,
    murcko_frag_itr_rnd=1,
    _search_scf_thr_=0.3,
    use_fingerprint_db=True,
    scaffold_db_file=None,
    fingerprint_db_file=None,
    seed=42,
    calculate_shape=True,
    trust_legacy_db=False,
    overwrite=False,
):
    """Return (results DataFrame, resource metrics); record provenance and diagnostics.

    Python API and CLI both apply Lipinski filters unless thresholds={} is supplied.
    Existing run artifacts require overwrite=True, protecting earlier results.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if not overwrite and any(
        (out / name).exists() for name in ("run_summary.json", "overall_result.txt", "run_manifest.json")
    ):
        raise FileExistsError(f"Output already contains a run: {out}. Use --overwrite or another directory.")
    if overwrite:
        for pattern in (
            "fragment_result_*.txt",
            "replace_scaffold_list.fragment_*.tsv",
            "fragment_info.tsv",
            "overall_result.txt",
            "run_summary.json",
            "run_manifest.json",
            "resource_cost.json",
        ):
            for old in out.glob(pattern):
                old.unlink()
    report = {"schema_version": 1}
    try:
        thresholds = utils._default_thrs_() if candidate_thresholds is None else candidate_thresholds
        _validate_parameters(tanimoto_threshold, _search_scf_thr_, seed, thresholds)
        utils._scaffold_no_reassign_(
            overall_max_n, frag_max_n, scaffold_top_n, cand_max_n__rplc, _merge_structure_top_n_
        )
        if not isinstance(murcko_frag_itr_rnd, int) or murcko_frag_itr_rnd < 0:
            raise ValueError("Fragment depth must be a nonnegative integer")
        target_mol, canonical = utils.init_mol(target_smiles)
        if len(Chem.GetMolFrags(target_mol)) != 1:
            raise utils.SmilesInputError(
                target_smiles, "Multi-component input: supply one compound without counterions"
            )
        core_mol = utils.init_mol(core_smiles)[0] if core_smiles is not None else None
        if core_mol is not None and not target_mol.HasSubstructMatch(core_mol, useChirality=True):
            raise utils.SmilesInputError(core_smiles, "Core SMILES is not contained in the input")
        report["input_smiles"] = canonical
        report["input_properties"] = utils._get_molecular_prop_extd_(target_mol)
        _, input_failures = utils.determ_props(report["input_properties"], thresholds, True)
        report["input_filter_failures"] = input_failures
        if input_failures:
            LOGGER.warning("Input exceeds active candidate filters: %s", "; ".join(input_failures))
        if fragments_DB is not None:
            fragments_DB = list(fragments_DB)
        args = dict(
            target_smiles=canonical,
            target_mol=target_mol,
            fragments_DB=fragments_DB,
            overall_max_n=overall_max_n,
            frag_max_n=frag_max_n,
            scaffold_top_n=scaffold_top_n,
            cand_max_n__rplc=cand_max_n__rplc,
            _merge_structure_top_n_=_merge_structure_top_n_,
            murcko_frag_itr_rnd=murcko_frag_itr_rnd,
            _search_scf_thr_=_search_scf_thr_,
            fragments=fragments,
            replace_scaffold_files=replace_scaffold_files,
            output_dir=out,
            low_mem=low_mem,
            tqdm_quiet=tqdm_quiet,
            use_fingerprint_db=use_fingerprint_db,
            scaffold_db_file=scaffold_db_file,
            fingerprint_db_file=fingerprint_db_file,
            trust_legacy_db=trust_legacy_db,
            diagnostics=report,
        )
        frags, overall, per_frag, top_n, per_replacement, beam = get_frags_cands(**args)
        parameters = {
            k: v
            for k, v in args.items()
            if k not in ("target_mol", "fragments_DB", "diagnostics", "output_dir")
        }
        parameters.update(
            input_smiles=target_smiles,
            core_smiles=core_smiles,
            seed=seed,
            calculate_shape=calculate_shape,
            tanimoto_threshold=tanimoto_threshold,
            candidate_thresholds={
                k: [v if v is None or math.isfinite(v) else str(v) for v in bounds]
                for k, bounds in thresholds.items()
            },
            effective_limits={
                "overall": overall,
                "per_fragment": per_frag,
                "scaffolds": top_n,
                "per_replacement": per_replacement,
                "beam": beam,
            },
        )
        # Normalize paths for JSON and hash exactly the actual libraries in use.
        parameters = json.loads(json.dumps(parameters, default=str))
        sources = (
            list(replace_scaffold_files or [])
            if not isinstance(replace_scaffold_files, (str, Path))
            else [replace_scaffold_files]
        )
        if report["database"]["path"]:
            sources.append(report["database"]["path"])
        manifest = _provenance(parameters, sources)
        if fragments_DB is not None:
            values = [s if isinstance(s, str) else Chem.MolToSmiles(s) for s in fragments_DB]
            manifest["in_memory_library_sha256"] = hashlib.sha256("\n".join(values).encode()).hexdigest()
        _json_write(out / "run_manifest.json", manifest)
        path = make_scaffold_hopping(
            canonical,
            target_mol,
            frags,
            core_mol,
            tanimoto_threshold,
            overall,
            per_frag,
            top_n,
            per_replacement,
            beam,
            out,
            thresholds,
            seed,
            calculate_shape,
            report,
        )
        return pd.read_csv(path, sep="\t")
    except Exception as exc:
        report.update(status="error", error_type=type(exc).__name__, message=str(exc))
        _json_write(out / "run_summary.json", report)
        raise


def main(argv=None):
    parser = cli.argument_parser()
    options = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    from chembounce_synthesis.options import apply_expansion_options, has_synthesis_overrides

    synthesis_config = None
    if options.synllama is True and not options.synthesis_config:
        LOGGER.error("--synllama requires --synthesis-config with provider.kind='synllama'")
        return 2
    if options.synllama is not False and options.synthesis_config:
        from chembounce_synthesis.workflow import _config, load_config

        try:
            synthesis_config = _config(
                apply_expansion_options(load_config(options.synthesis_config), options)
            )
            if options.synllama is True and synthesis_config["provider"].get("kind") != "synllama":
                raise ValueError(
                    "--synllama requires provider.kind='synllama'; omit it to replay saved proposals"
                )
            if options.core_smiles:
                if synthesis_config.get("core_smiles") not in (None, options.core_smiles):
                    raise ValueError("core_smiles conflicts with the synthesis configuration")
                synthesis_config["core_smiles"] = options.core_smiles
        except (ValueError, OSError, KeyError, TypeError) as exc:
            LOGGER.error("Synthesis configuration: %s", exc)
            return 2
    elif options.synllama is not False and has_synthesis_overrides(options):
        LOGGER.error("Synthesis search/refinement options require --synthesis-config")
        return 2
    thresholds = utils._default_thrs_(not options.wo_lipinski)
    for metric, prefix in (
        ("qed", "qed"),
        ("sa", "sa"),
        ("logp", "logp"),
        ("mw", "mw"),
        ("n_hdonor", "h_donor"),
        ("n_hacceptor", "h_acceptor"),
    ):
        lo, hi = getattr(options, prefix + "_min"), getattr(options, prefix + "_max")
        if lo is not None or hi is not None:
            previous = thresholds.get(metric, (None, None))
            thresholds[metric] = (
                lo if lo is not None else previous[0],
                hi if hi is not None else previous[1],
            )
    try:
        results, cost = chembounce(
            options.input_smiles,
            output_dir=options.output_dir,
            core_smiles=options.core_smiles,
            overall_max_n=options.overall_max_n,
            frag_max_n=options.frag_max_n,
            scaffold_top_n=options.scaffold_top_n,
            cand_max_n__rplc=options.cand_max_n__rplc,
            _merge_structure_top_n_=options.beam_width,
            tanimoto_threshold=options.tanimoto_threshold,
            fragments=options.fragments,
            replace_scaffold_files=options.replace_scaffold_files,
            candidate_thresholds=thresholds,
            low_mem=options.low_mem,
            use_fingerprint_db=options.use_fingerprint_db,
            scaffold_db_file=options.scaffold_db_file,
            fingerprint_db_file=options.fingerprint_db_file,
            seed=options.seed,
            calculate_shape=not options.no_shape,
            trust_legacy_db=options.trust_legacy_db,
            overwrite=options.overwrite,
            _search_scf_thr_=options.search_threshold,
            murcko_frag_itr_rnd=options.fragment_depth,
        )
    except (ValueError, OSError, utils.SmilesInputError) as exc:
        LOGGER.error("%s", exc)
        return 2
    _json_write(Path(options.output_dir) / "resource_cost.json", cost)
    if synthesis_config is not None:
        from chembounce_synthesis.workflow import run

        try:
            config = synthesis_config
            # Projection can change molecular properties. Reapply the intersection
            # of hopping and synthesis bounds to the regenerated final products.
            bounds = dict(config.get("property_thresholds", {}))
            for key, (lo, hi) in thresholds.items():
                extra_lo, extra_hi = bounds.get(key, (None, None))
                lows = [v for v in (lo, extra_lo) if v is not None and math.isfinite(v)]
                highs = [v for v in (hi, extra_hi) if v is not None and math.isfinite(v)]
                bounds[key] = [max(lows) if lows else None, min(highs) if highs else None]
            config["property_thresholds"] = bounds
            verified = run(
                options.input_smiles,
                results["Final structure"].tolist(),
                config,
                Path(options.output_dir) / "synthesis",
            )
        except (ValueError, OSError, KeyError, TypeError) as exc:
            LOGGER.error("Synthesis workflow: %s", exc)
            return 2
        print(
            f"Verified synthesis products: {len(verified)}; results: {Path(options.output_dir) / 'synthesis/results.tsv'}"
        )
        return 0 if verified else 3
    print(
        f"Found hopped structures: {len(results)}; results: {Path(options.output_dir) / 'overall_result.txt'}"
    )
    print(f"Diagnostics: {Path(options.output_dir) / 'run_summary.json'}")
    return 0 if len(results) else 3


if __name__ == "__main__":
    raise SystemExit(main())
