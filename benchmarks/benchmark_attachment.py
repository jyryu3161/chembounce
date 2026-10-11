#!/usr/bin/env python3
"""Controlled attachment comparison and beam ablation, independent of search.

Both versions receive the same RDKit molecules, full Murcko query and four
replacement choices. Each call runs in a disposable process with a wall budget.
This is a correctness/stress experiment, not an activity or synthesis benchmark.
Requires the legacy dependencies when --legacy-tag is supplied.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def worker(connection, code, case, beam, max_n):
    """Import one implementation in a fresh interpreter, then signal readiness."""
    sys.path.insert(0, str(code))
    import rdkit
    from rdkit import Chem, rdBase

    os.environ["RDBASE"] = str(Path(rdkit.__file__).parent)
    import chembounce

    molecules = [Chem.MolFromSmiles(case[key]) for key in ("input", "fragment", "replacement")]
    connection.send({"ready": True})
    started = time.perf_counter()
    try:
        with rdBase.BlockLogs():
            candidates = chembounce.replace_molecule(*molecules, max_n, beam)
        result = {"status": "success" if candidates else "empty", "candidates": candidates}
    except Exception as exc:
        result = {"status": "error", "error_type": type(exc).__name__, "message": str(exc)}
    result["seconds"] = time.perf_counter() - started
    connection.send(result)
    connection.close()


def run_case(code, case, beam, max_n, timeout):
    context = multiprocessing.get_context("spawn")
    parent, child = context.Pipe(duplex=False)
    process = context.Process(target=worker, args=(child, code, case, beam, max_n))
    process.start()
    child.close()
    try:
        if not parent.poll(60):
            return {"status": "import_timeout"}
        try:
            assert parent.recv() == {"ready": True}
            if parent.poll(timeout):
                return parent.recv()
            return {"status": "timeout", "timeout_seconds": timeout}
        except EOFError:
            return {"status": "worker_crash", "exitcode": process.exitcode}
    finally:
        if process.is_alive():
            process.terminate()
        process.join()
        parent.close()


def build_cases(panel):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    cases = []
    for record in panel:
        mol = Chem.MolFromSmiles(record["smiles"])
        if mol is None or len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError(f"Invalid panel entry: {record['name']}")
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        Chem.RemoveStereochemistry(scaffold)
        if not scaffold.GetNumAtoms():
            raise ValueError(f"Panel entry has no Murcko scaffold: {record['name']}")
        fragment = Chem.MolToSmiles(scaffold)
        # Identity is an assembly control, not a scaffold-hopping candidate.
        choices = {
            "identity_control": fragment,
            "pyridine": "c1ccncc1",
            "piperidine": "C1CCNCC1",
            "morpholine": "C1COCCN1",
        }
        for label, replacement in choices.items():
            cases.append(
                {
                    "name": record["name"],
                    "choice": label,
                    "input": Chem.MolToSmiles(mol),
                    "fragment": fragment,
                    "replacement": replacement,
                }
            )
    return cases


def assess(case, result, max_n):
    """Independent output checks, with retained fragment queries from the boundary cut.

    Wildcard endpoints permit a changed neighbour while checking specified stereo
    in retained atoms. This is a necessary substructure check, not atom-mapped proof
    for repeated identical substituents. Heavy atom/isotope accounting is also checked.
    """
    from rdkit import Chem, DataStructs
    from rdkit.Chem.Scaffolds import MurckoScaffold
    import chemistry

    target = Chem.MolFromSmiles(case["input"])
    replacement = Chem.MolFromSmiles(case["replacement"])
    outside = chemistry.cut_scaffold(target, Chem.MolFromSmiles(case["fragment"]))
    params = Chem.AdjustQueryParameters.NoAdjustments()
    params.makeDummiesQueries = True
    # Without degree constraints a terminal methyl query can match an aromatic
    # ipso carbon and the wildcard can match that methyl, falsely accepting an
    # inverted stereocenter by swapping the two neighbours.
    params.adjustDegree = True
    params.adjustDegreeFlags = Chem.ADJUST_IGNOREDUMMIES
    queries = []

    def atoms(mol):
        return Counter((a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms() if a.GetAtomicNum() > 1)

    expected_atoms = atoms(replacement)
    for part in outside:
        expected_atoms.update(atoms(part))
        for atom in part.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetIsotope(0)
        query = Chem.AdjustQueryProperties(part, params)
        assert target.HasSubstructMatch(query, useChirality=True), "Invalid retention query"
        queries.append(query)

    def scaffold(mol):
        core = MurckoScaffold.GetScaffoldForMol(mol)
        Chem.RemoveStereochemistry(core)
        return Chem.MolToSmiles(core)

    original_core = scaffold(target)
    original_fp = chemistry.count_fingerprint(target)
    failures = Counter()
    products = []
    changed = 0
    for score, smiles in result.get("candidates", []):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            failures["invalid_smiles"] += 1
            continue
        canonical = Chem.MolToSmiles(mol)
        products.append(canonical)
        if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
            failures["disconnected_or_dummy"] += 1
        if atoms(mol) != expected_atoms:
            failures["heavy_atom_or_isotope_loss"] += 1
        if not all(mol.HasSubstructMatch(query, useChirality=True) for query in queries):
            failures["retained_substructure_or_stereo"] += 1
        reference = DataStructs.TanimotoSimilarity(original_fp, chemistry.count_fingerprint(mol))
        if abs(score - reference) > 1e-12:
            failures["score_mismatch"] += 1
        changed += scaffold(mol) != original_core
    if len(products) != len(set(products)):
        failures["duplicate_products"] += 1
    if len(result.get("candidates", [])) > max_n:
        failures["output_cap"] += 1
    return {
        "checks_failed": dict(failures),
        "valid_unique": len(set(products)),
        "changed_murcko": changed,
        "input_recovered": case["input"] in products,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, default=Path(__file__).with_name("publication_molecules.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--legacy-tag", default="1.0-original")
    parser.add_argument(
        "--cases-file", type=Path, help="Exact input/fragment/replacement cases instead of a panel"
    )
    parser.add_argument("--beams", type=int, nargs="+", default=[10, 50, 100])
    parser.add_argument("--max-n", type=int, default=10)
    parser.add_argument("--timeout", type=float, default=10)
    parser.add_argument("--names", nargs="+")
    args = parser.parse_args()
    if min(args.beams) < 1 or args.max_n < 1 or args.timeout <= 0:
        parser.error("Limits must be positive")
    if args.output.exists():
        parser.error("Output exists; choose a new file to preserve previous evidence")
    sys.path.insert(0, str(ROOT))
    import rdkit

    panel = json.loads(args.panel.read_text())["molecules"]
    if args.names:
        panel = [row for row in panel if row["name"] in args.names]
        if {row["name"] for row in panel} != set(args.names):
            parser.error("Unknown panel name")
    cases = json.loads(args.cases_file.read_text()) if args.cases_file else build_cases(panel)
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "rdkit": rdkit.__version__,
        "platform": platform.platform(),
        "panel_sha256": digest(args.panel),
        "cases_file_sha256": digest(args.cases_file) if args.cases_file else None,
        "script_sha256": digest(__file__),
        "max_n": args.max_n,
        "timeout_seconds_after_import": args.timeout,
        "cases": cases,
        "current_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "git_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT)),
        "source_sha256": {p.name: digest(p) for p in sorted(ROOT.glob("*.py"))},
        "results": [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as directory:
        configurations = []
        if args.legacy_tag:
            commit = subprocess.check_output(
                ["git", "rev-parse", args.legacy_tag + "^{commit}"], cwd=ROOT, text=True
            ).strip()
            names = subprocess.check_output(
                ["git", "ls-tree", "--name-only", commit], cwd=ROOT, text=True
            ).splitlines()
            for name in names:
                if name.endswith(".py"):
                    (Path(directory) / name).write_bytes(
                        subprocess.check_output(["git", "show", commit + ":" + name], cwd=ROOT)
                    )
            report["legacy_commit"] = commit
            report["legacy_compatibility_shim"] = "RDBASE set to installed pip RDKit; tagged source unchanged"
            configurations.append((args.legacy_tag, directory, 50))
        configurations.extend(("2.0", str(ROOT), beam) for beam in args.beams)
        for version, code, beam in configurations:
            for index, case in enumerate(cases):
                result = run_case(code, case, beam, args.max_n, args.timeout)
                result.update(assess(case, result, args.max_n))
                result.update(version=version, beam_width=beam, case=index)
                report["results"].append(result)
                args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
                print(
                    f"{version} beam={beam} {case['name']}/{case['choice']}: {result['status']} "
                    f"n={result['valid_unique']} failures={result['checks_failed']}",
                    flush=True,
                )


if __name__ == "__main__":
    main()
