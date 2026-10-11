"""ADP regressions for the narrowly backported 1.x attachment fix."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import chembounce as cb
import utils

ADP = "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N"


@pytest.fixture(scope="module")
def database(tmp_path_factory):
    smiles = ["C1CCNC1", "C1CCOC1", "c1ccncc1", "c1ncc2cc[nH]c2n1", "c1ncc2[nH]cnc2n1"]
    bits = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for index, value in enumerate(smiles):
        fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(value), 2, nBits=2048)
        DataStructs.ConvertToNumpyArray(fp, bits[index])
    path = tmp_path_factory.mktemp("library") / "scaffolds.npz"
    np.savez(path, fingerprints=bits, smiles=np.asarray(smiles), radius=2, fp_size=2048)
    return str(path)


def test_all_adp_scaffolds_assemble():
    mol, _ = utils.init_mol(ADP)
    fragments = cb.sg.get_all_murcko_fragments(mol, False, 1)
    assert len(fragments) == 3
    for fragment in fragments:
        candidates = cb.replace_molecule(mol, fragment, Chem.MolFromSmiles("c1ccncc1"), 3, 30)
        assert candidates
        for _, smiles in candidates:
            product = Chem.MolFromSmiles(smiles)
            assert product is not None and len(Chem.GetMolFrags(product)) == 1
            assert all(atom.GetAtomicNum() for atom in product.GetAtoms())


def test_ribose_stereo_and_phosphates_survive():
    mol, canonical = utils.init_mol(ADP)
    assert canonical == Chem.MolToSmiles(Chem.MolFromSmiles(ADP))
    core = Chem.MolFromSmiles("O=P(O)(O)OP(=O)(O)OC[C@H]1OC[C@H](O)[C@@H]1O")
    assert mol.HasSubstructMatch(core, useChirality=True)
    candidates = cb.replace_molecule(
        mol, Chem.MolFromSmiles("c1ncc2nc[nH]c2n1"), Chem.MolFromSmiles("c1ccncc1"), 3, 30
    )
    assert candidates
    assert all(Chem.MolFromSmiles(s).HasSubstructMatch(core, useChirality=True) for _, s in candidates)


@pytest.mark.parametrize("lipinski", [False, True])
def test_full_workflow_reports_results_or_filter_rejection(database, tmp_path, lipinski):
    result, _ = cb.chembounce(
        ADP,
        fingerprint_db_file=database,
        output_dir=str(tmp_path),
        frag_max_n=2,
        scaffold_top_n=2,
        cand_max_n__rplc=2,
        _merge_structure_top_n_=30,
        candidate_thresholds=utils._default_thrs_(lipinski=lipinski),
    )
    summary = json.loads((tmp_path / "run_summary.json").read_text())
    assert summary["counts"].get("assembly_failed", 0) == 0
    assert summary["counts"]["accepted"] == len(result)
    if lipinski:
        assert result.empty
        assert summary["status"] == "no_candidates"
        assert summary["counts"]["property_rejected"] > 0
        assert "--wo_lipinski" in summary["message"]
    else:
        assert not result.empty and summary["status"] == "success"
        assert result["Final structure"].is_unique
        assert result["Tanimoto Similarity"].ge(0.5).all()
        assert result["Electron shape Similarity"].notna().any()


def test_assembly_errors_are_reported(database, tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise ValueError("deliberate attachment error")

    monkeypatch.setattr(cb, "replace_molecule", fail)
    result, _ = cb.chembounce(
        ADP,
        fingerprint_db_file=database,
        output_dir=str(tmp_path),
        frag_max_n=1,
        scaffold_top_n=1,
        cand_max_n__rplc=1,
    )
    assert result.empty
    summary = json.loads((tmp_path / "run_summary.json").read_text())
    assert summary["counts"]["assembly_failed"] > 0
    assert "deliberate attachment error" in summary["assembly_errors"][0]["reason"]


def test_cli_adp_can_produce_final_output(database, tmp_path):
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "chembounce.py"),
            "-i",
            ADP,
            "-o",
            str(tmp_path),
            "--fingerprint-db",
            database,
            "--wo_lipinski",
            "--scaffold_top_n",
            "2",
            "-n",
            "2",
            "--cand_max_n__rplc",
            "2",
        ],
        check=False,
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert not pd.read_csv(tmp_path / "overall_result.txt", sep="\t").empty
