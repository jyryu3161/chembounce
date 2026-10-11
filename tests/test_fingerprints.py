import numpy as np
import pytest
from rdkit import Chem, DataStructs
import fingerprint_db as fp
from generate_fingerprints import generate_fingerprints, verify_fingerprints
from chembounce import search_similar_scaffolds

SMILES = ["CCO", "CCCO", "CCN", "c1ccncc1", "c1ccccc1", "CCCl"]


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("top_n", [None, 1, 3])
def test_bulk_matches_rdkit_with_ties_and_empty_vectors(packed, top_n):
    rng = np.random.default_rng(2)
    bits = rng.integers(0, 2, (31, 64), dtype=np.uint8)
    bits[0] = 0
    bits[2] = bits[1]
    for query in (bits[0], bits[1]):
        expected = [
            DataStructs.TanimotoSimilarity(
                DataStructs.CreateFromBitString("".join(map(str, query))),
                DataStructs.CreateFromBitString("".join(map(str, row))),
            )
            for row in bits
        ]
        order = np.lexsort((np.arange(len(bits)), -np.array(expected)))[:top_n]
        actual_i, actual_s = fp.calculate_bulk_tanimoto(
            np.packbits(query) if packed else query,
            np.packbits(bits, axis=1) if packed else bits,
            packed=packed,
            top_n=top_n,
            chunk_size=4,
        )
        np.testing.assert_array_equal(actual_i, order)
        np.testing.assert_allclose(actual_s, np.array(expected)[order])


@pytest.mark.parametrize("suffix", [".npz", ""])
def test_generate_verify_and_search_equivalence(tmp_path, suffix):
    source, output = tmp_path / "scaffolds.smi", tmp_path / ("db" + suffix)
    source.write_text("\n".join(SMILES + ["OCC", "bad", "CC.O"]))
    assert generate_fingerprints(source, output, batch_size=2) == len(SMILES)
    assert verify_fingerprints(output) == len(SMILES)
    database = fp.load_fingerprint_db(output)
    assert database["fingerprints"].shape == (len(SMILES), 256)
    if not suffix:
        assert isinstance(database["fingerprints"], np.memmap)
    query = Chem.MolFromSmiles("CCO")
    actual = fp.search_similar_scaffolds_fast(query, database, 3, 0)
    expected = search_similar_scaffolds(query, SMILES, scaffold_top_n=3, threshold=0)
    assert list(actual.index) == list(expected.index)
    np.testing.assert_allclose(actual, expected)
    assert "CCO" not in actual.index


def test_different_structures_with_same_fingerprint_are_retained():
    # Small fp_size deliberately exposes bit collisions.
    data = {
        "radius": 2,
        "fp_size": 1,
        "smiles": np.array(["CC", "CCC"]),
        "fingerprints": np.ones((2, 1), dtype=np.uint8),
        "packed": False,
    }
    result = fp.search_similar_scaffolds_fast(Chem.MolFromSmiles("CC"), data, 1, 1)
    assert list(result.index) == ["CCC"]
    assert result.iloc[0] == 1


def test_top_k_does_not_lose_slots_to_duplicates_or_identity():
    data = {
        "radius": 2,
        "fp_size": 1,
        "smiles": np.array(["CC", "CCC", "CCC", "CCCC"]),
        "fingerprints": np.ones((4, 1), dtype=np.uint8),
        "packed": False,
    }
    result = fp.search_similar_scaffolds_fast(Chem.MolFromSmiles("CC"), data, 2, 0)
    assert list(result.index) == ["CCC", "CCCC"]


def test_legacy_pickle_requires_explicit_trust(tmp_path):
    path = tmp_path / "legacy.npz"
    np.savez(
        path,
        smiles=np.array(["CC"], dtype=object),
        fingerprints=np.ones((1, 2), dtype=np.uint8),
        radius=2,
        fp_size=2,
    )
    with pytest.raises(ValueError, match="explicit trust"):
        fp.load_fingerprint_db(path)
    assert fp.load_fingerprint_db(path, trust_legacy=True)["smiles"][0] == "CC"


def test_corrupted_database_rejected(tmp_path):
    path = tmp_path / "bad.npz"
    np.savez(path, smiles=np.array(["CC"]), fingerprints=np.ones((2, 3), dtype=np.uint8), radius=2, fp_size=3)
    with pytest.raises(ValueError, match="dimensions"):
        fp.load_fingerprint_db(path)


def test_verifier_actually_checks_small_databases(tmp_path):
    path = tmp_path / "bad.npz"
    np.savez(
        path,
        smiles=np.array(["CC"]),
        fingerprints=np.zeros((1, 2048), dtype=np.uint8),
        radius=2,
        fp_size=2048,
    )
    with pytest.raises(ValueError, match="does not match"):
        verify_fingerprints(path)


@pytest.mark.parametrize(
    "metadata",
    [
        {"schema_version": 999},
        {"bitorder": "little"},
        {"include_chirality": True},
        {"fp_size": 3.5},
        {"radius": -1},
    ],
)
def test_incompatible_metadata_rejected(tmp_path, metadata):
    path = tmp_path / "bad.npz"
    data = dict(smiles=np.array(["CC"]), fingerprints=np.ones((1, 3), dtype=np.uint8), radius=2, fp_size=3)
    data.update(metadata)
    np.savez(path, **data)
    with pytest.raises(ValueError):
        fp.load_fingerprint_db(path)


def test_generation_records_excluded_records_without_pickle(tmp_path):
    import json

    source, output = tmp_path / "source.smi", tmp_path / "db.npz"
    source.write_text("CC\nCC\n\nCC.O\nCC*\nbad\n")
    generate_fingerprints(source, output)
    data = fp.load_fingerprint_db(output)
    assert data["rows"] == 1
    assert data["duplicate_lines"] == 1
    assert data["blank_lines"] == 1
    assert data["rejected_lines"] == 3
    assert json.loads(str(data["rejection_counts_json"])) == {
        "invalid_smiles": 1,
        "multiple_components": 1,
        "dummy_atoms": 1,
    }
    with np.load(output, allow_pickle=False) as archive:
        for key in archive.files:
            archive[key]


def test_legacy_invalid_hits_and_equivalent_smiles_do_not_consume_top_k():
    strings = ["CC*", "CC.CC", "C(C)", "CCC", "C(CC)", "CCCC"]
    data = {
        "radius": 2,
        "fp_size": 1,
        "smiles": np.array(strings),
        "fingerprints": np.ones((6, 1), dtype=np.uint8),
        "packed": False,
    }
    result = fp.search_similar_scaffolds_fast(Chem.MolFromSmiles("CC"), data, 2, 0)
    assert list(result.index) == ["CCC", "CCCC"]
    assert len(data["_legacy_rejected_rows"]) == 2
