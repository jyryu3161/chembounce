"""Actual source strings whose first canonical serialization is not a fixed point."""

import json
from unittest.mock import Mock

import numpy as np
import pytest
from rdkit import Chem, rdBase

import fingerprint_db as fp
import generate_fingerprints as builder


CONVERGENT_SOURCES = ["[C@]12C[C@H](CNC1=O)C2", "[C@]12CC[C@](c3ccccc3)(C1)C2"]
CYCLE_SOURCES = [
    "[C@]1(c2ccccc2)C[C@@H](NC(=O)c2nc[nH]c2)C1",
    "[C@]1(c2ccccc2)C[C@H](N2CCCC2=O)C1",
    "[C@]1(c2ccccc2)C[C@@H](N2CCCC2=O)C1",
]
REAL_CASES = pytest.mark.skipif(
    rdBase.rdkitVersion != "2025.09.6", reason="Empirical serialization cycles specific to RDKit 2025.09.6"
)


@REAL_CASES
@pytest.mark.parametrize("suffix", ["", ".npz"])
def test_actual_builder_sources_stabilize_or_record_original_failed_inputs(tmp_path, suffix):
    source = tmp_path / "source.smi"
    source.write_text("\n".join(CONVERGENT_SOURCES + CYCLE_SOURCES) + "\n")
    output = tmp_path / ("db" + suffix)
    assert builder.generate_fingerprints(source, output) == 2
    data = fp.load_fingerprint_db(output)
    assert data["canonicalization_policy"] == "rdkit-roundtrip-fixed-point-v1"
    assert data["canonicalization_max_roundtrips"] == 8
    assert data["canonicalization_adjusted_lines"] == 2
    assert data["rejected_lines"] == 3
    assert json.loads(str(data["rejection_counts_json"])) == {"canonicalization_unstable": 3}
    assert json.loads(str(data["canonicalization_failure_counts_json"])) == {"cycle": 3}
    failures = json.loads(str(data["canonicalization_failure_examples_json"]))
    assert [row["smiles"] for row in failures] == CYCLE_SOURCES
    assert [row["line"] for row in failures] == [3, 4, 5]
    assert all(row["reason"] == "cycle" and row["roundtrips"] == 2 for row in failures)
    for index, smiles in enumerate(data["smiles"]):
        parsed = Chem.MolFromSmiles(smiles)
        assert Chem.MolToSmiles(parsed) == smiles
        np.testing.assert_array_equal(
            data["fingerprints"][index], np.packbits(fp.generate_query_fingerprint(parsed))
        )
        assert smiles not in fp.search_similar_scaffolds_fast(parsed, data, 1, 0).index


def test_serialized_molecule_is_used_for_bits_and_exact_states_remain_distinct(tmp_path, monkeypatch):
    states = [
        "N[C@H](C)[13C](=O)[O-]",
        "N[C@@H](C)[13C](=O)[O-]",
        "N[C@H](C)C(=O)[O-]",
        "N[C@H](C)[13C](=O)O",
        "F/C=C/F",
        "F/C=C\\F",
        "O=c1[nH]cccc1",
        "Oc1ncccc1",
    ]
    source = tmp_path / "source.smi"
    source.write_text("\n".join(states) + "\n")
    fingerprint = builder.generate_query_fingerprint
    atom_orders = []

    def inspect_molecule(mol, radius, fp_size):
        serialized = Chem.MolToSmiles(mol)
        expected = Chem.MolFromSmiles(serialized)
        assert Chem.MolToSmiles(mol, canonical=False) == Chem.MolToSmiles(expected, canonical=False)
        atom_orders.append(serialized)
        return fingerprint(mol, radius, fp_size)

    monkeypatch.setattr(builder, "generate_query_fingerprint", inspect_molecule)
    output = tmp_path / "db"
    assert builder.generate_fingerprints(source, output) == len(states)
    data = fp.load_fingerprint_db(output)
    assert set(data["smiles"]) == {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in states}
    assert len(atom_orders) == len(states)
    assert set(atom_orders) == set(data["smiles"])


def test_noncyclic_nonconvergence_has_an_explicit_small_work_bound(monkeypatch):
    # A writer that never repeats cannot bypass the iteration budget.
    mol = Chem.MolFromSmiles("C")
    writer = Mock(side_effect=["C", "CC", "CCC", "CCCC", "CCCCC"])
    parser = Mock(wraps=Chem.MolFromSmiles)
    monkeypatch.setattr(fp.Chem, "MolToSmiles", writer)
    monkeypatch.setattr(fp.Chem, "MolFromSmiles", parser)
    with pytest.raises(ValueError, match="limit") as caught:
        fp.stable_canonical_molecule(mol, max_roundtrips=3)
    assert caught.value.reason == "limit" and caught.value.roundtrips == 3
    assert parser.call_count == 3


@REAL_CASES
@pytest.mark.parametrize("source", CYCLE_SOURCES)
def test_unresolved_query_is_rejected_before_scoring_without_mutating_caller(monkeypatch, source):
    query = Chem.MolFromSmiles(source)
    for index, atom in enumerate(query.GetAtoms(), 1):
        atom.SetAtomMapNum(index)
    before = Chem.MolToSmiles(query)
    scorer = Mock(side_effect=AssertionError("Unresolved query cannot be ranked"))
    monkeypatch.setattr(fp, "calculate_bulk_tanimoto", scorer)
    data = {
        "rdkit_version": rdBase.rdkitVersion,
        "radius": 2,
        "fp_size": 1,
        "canonical_sorted_unique": True,
        "canonicalization_policy": "rdkit-roundtrip-fixed-point-v1",
        "smiles": np.array(["CC"]),
        "fingerprints": np.ones((1, 1), dtype=np.uint8),
    }
    with pytest.raises(ValueError, match="cycle"):
        fp.search_similar_scaffolds_fast(query, data, 1, 0)
    assert Chem.MolToSmiles(query) == before
    scorer.assert_not_called()


@REAL_CASES
@pytest.mark.parametrize("version", [None, rdBase.rdkitVersion])
def test_old_claims_do_not_bypass_cycle_rejection_and_refill(version):
    strings = sorted(
        [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in CYCLE_SOURCES] + ["c1ccccc1", "c1ccncc1"]
    )
    data = {
        "radius": 2,
        "fp_size": 1,
        "canonical_sorted_unique": True,
        "smiles": np.array(strings),
        "fingerprints": np.ones((len(strings), 1), dtype=np.uint8),
    }
    if version is not None:
        data["rdkit_version"] = version
    result = fp.search_similar_scaffolds_fast(Chem.MolFromSmiles("N"), data, 2, 0)
    assert list(result.index) == ["c1ccccc1", "c1ccncc1"]
    assert len(data["_legacy_rejected_rows"]) == 3
    assert all(row["reason"] == "cycle" for row in data["_canonicalization_rejected_rows"].values())


@REAL_CASES
def test_all_cyclic_inputs_fail_without_publishing_partial_database(tmp_path):
    source, output = tmp_path / "source.smi", tmp_path / "db"
    source.write_text("\n".join(CYCLE_SOURCES) + "\n")
    with pytest.raises(ValueError, match="No valid connected scaffolds") as caught:
        builder.generate_fingerprints(source, output)
    assert '"cycle": 3' in str(caught.value)
    assert all(source_smiles in str(caught.value) for source_smiles in CYCLE_SOURCES)
    assert not output.exists()
