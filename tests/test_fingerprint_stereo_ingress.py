"""Synthetic exact-state boundaries; no 3D, models or production database scan."""

import json
from unittest.mock import Mock

import numpy as np
import pytest
from rdkit import Chem, rdBase

import fingerprint_db as fp
import generate_fingerprints as builder
from chembounce_synthesis.common import canonical, parse_absolute_smiles

PLAIN = "C[C@H](F)[C@@H](O)Cl"
ATROP = "CC1=CC=CC(I)=C1N1C(C)=CC=C1Br"
UNSUPPORTED = [PLAIN + " |o1:1,3|", PLAIN + " |&1:1,3|", ATROP + " |wU:7.7|", ATROP + " |wD:7.7|"]
SUPPORTED = [PLAIN, "[H][C@](F)(Cl)Br |a:1|", "[2H][C@](F)(Cl)Br |a:1|", "F/C=C/F", "F/C=C\\F"]


@pytest.fixture(autouse=True)
def forbid_3d(monkeypatch):
    from rdkit.Chem import AllChem, rdMolAlign

    def forbidden(*args, **kwargs):
        raise AssertionError("Fingerprint fixtures must not embed or align")

    monkeypatch.setattr(AllChem, "EmbedMultipleConfs", forbidden)
    monkeypatch.setattr(AllChem, "EmbedMolecule", forbidden)
    monkeypatch.setattr(rdMolAlign, "AlignMol", forbidden)


@pytest.mark.parametrize("suffix", ["", ".npz"])
def test_mixed_builder_preserves_state_bits_ordinals_and_rejections(tmp_path, suffix):
    source, output = tmp_path / "tiny.smi", tmp_path / ("db" + suffix)
    # A plain/ABS duplicate must not be invented from any unsupported row.
    texts = UNSUPPORTED + SUPPORTED + [PLAIN + " |a:1,3|", ""]
    source.write_text("\n".join(texts) + "\n")
    before = source.read_bytes()
    assert builder.generate_fingerprints(source, output, fp_size=33, batch_size=2) == len(SUPPORTED)
    data = fp.load_fingerprint_db(output)
    assert data["rejected_lines"] == 4 and data["duplicate_lines"] == 1 and data["blank_lines"] == 1
    assert json.loads(str(data["rejection_counts_json"])) == {"unsupported_molecular_representation": 4}
    rejected = json.loads(str(data["rejection_examples_json"]))
    assert [r["line"] for r in rejected] == [1, 2, 3, 4]
    assert [r["smiles"] for r in rejected] == UNSUPPORTED
    assert all(r["detail"].startswith("Unsupported") for r in rejected)
    assert set(data["smiles"]) == {canonical(s) for s in SUPPORTED}
    for i, smiles in enumerate(data["smiles"]):
        # Independent RDKit generator oracle: settings remain chirality=False.
        generator = Chem.rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=33)
        bits = generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(str(smiles)))
        np.testing.assert_array_equal(data["fingerprints"][i], np.packbits(bits))
    assert builder.verify_fingerprints(output, sample_size=100) == len(SUPPORTED)
    assert source.read_bytes() == before


@pytest.mark.parametrize("text", UNSUPPORTED)
def test_raw_text_rejection_precedes_runtime_parser_and_does_not_publish(tmp_path, monkeypatch, text):
    source, output = tmp_path / "tiny.smi", tmp_path / "db"
    source.write_text(text + "\n")
    parser = Mock(wraps=Chem.MolFromSmiles)
    monkeypatch.setattr(Chem, "MolFromSmiles", parser)
    with pytest.raises(ValueError, match="unsupported_molecular_representation") as caught:
        builder.generate_fingerprints(source, output)
    assert text in str(caught.value)
    assert not output.exists()
    parser.assert_not_called()  # especially old RDKit must not discard wU/wD first


@pytest.mark.parametrize("text", UNSUPPORTED[:2])
def test_raw_mol_query_and_fixed_point_reject_before_work_without_mutation(monkeypatch, text):
    raw = Chem.MolFromSmiles(text)
    before = Chem.MolToCXSmiles(raw)
    writer, generator = (
        Mock(side_effect=AssertionError("no serialization")),
        Mock(side_effect=AssertionError("no bits")),
    )
    monkeypatch.setattr(Chem, "MolToSmiles", writer)
    monkeypatch.setattr(fp, "_generator", generator)
    for operation in [fp.stable_canonical_molecule, fp.generate_query_fingerprint]:
        with pytest.raises(ValueError, match="nonabsolute enhanced"):
            operation(raw)
    assert Chem.MolToCXSmiles(raw) == before
    writer.assert_not_called()
    generator.assert_not_called()


def test_raw_atrop_mol_rejects_where_runtime_retains_state():
    mol = Chem.MolFromSmiles(UNSUPPORTED[2])
    if not any(str(b.GetStereo()).startswith("STEREOATROP") for b in mol.GetBonds()):
        pytest.skip("Only raw atrop enum fixture unavailable; text rejection runs on all runtimes")
    for operation in [fp.stable_canonical_molecule, fp.generate_query_fingerprint]:
        with pytest.raises(ValueError, match="atropisomeric"):
            operation(mol)


def collision_db(strings, packed, claimed):
    # Deliberate one-bit collisions isolate selected-row qualification/refill.
    strings = sorted(strings)
    return {
        "smiles": np.array(strings),
        "fingerprints": np.full((len(strings), 1), 128 if packed else 1, dtype=np.uint8),
        "radius": 2,
        "fp_size": 1,
        "packed": packed,
        **(
            {
                "rdkit_version": rdBase.rdkitVersion,
                "canonical_sorted_unique": True,
                "canonicalization_policy": fp.CANONICALIZATION_POLICY,
            }
            if claimed
            else {}
        ),
    }


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("claimed", [False, True])
@pytest.mark.parametrize("text", UNSUPPORTED + ["C ordinary-name"])
def test_legacy_and_claimed_fast_path_cannot_return_unsupported_annotations(packed, claimed, text):
    data = collision_db([text, "N", "NN"], packed, claimed)
    result = fp.search_similar_scaffolds_fast(Chem.MolFromSmiles("O"), data, 2, 0)
    assert list(result.index) == ["N", "NN"]
    rejected = data["_representation_rejected_rows"]
    assert set(rejected) == {list(data["smiles"]).index(text)}
    assert data["_legacy_rejected_rows"] == set(rejected)


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("claimed", [False, True])
@pytest.mark.parametrize("query", [PLAIN, "[H][C@](F)(Cl)Br"])
def test_supported_abs_query_exclusion_duplicate_normalization_and_refill(packed, claimed, query):
    annotation = "a:1,3" if query == PLAIN else "a:1"
    alternatives = ["N", "NN"] if query == PLAIN else ["[Na+]", "[Ne]"]
    strings = [query + " |" + annotation + "|", canonical(query), *alternatives]
    data = collision_db(strings, packed, claimed)
    result = fp.search_similar_scaffolds_fast(parse_absolute_smiles(query), data, 2, 0)
    assert list(result.index) == alternatives
    assert not data.get("_representation_rejected_rows")


@pytest.mark.parametrize("packed", [False, True])
def test_abs_and_plain_candidate_duplicate_cannot_consume_refilled_slot(packed):
    data = collision_db([PLAIN, PLAIN + " |a:1,3|", "N"], packed, True)
    result = fp.search_similar_scaffolds_fast(Chem.MolFromSmiles("O"), data, 2, 0)
    assert list(result.index) == [PLAIN, "N"]


@pytest.mark.parametrize("text", UNSUPPORTED)
@pytest.mark.parametrize("packed", [False, True])
def test_sampled_verifier_rejects_unsupported_original_text(tmp_path, text, packed):
    path = tmp_path / "tiny.npz"
    np.savez(
        path,
        smiles=np.array([text]),
        fingerprints=np.array([[128 if packed else 1]], dtype=np.uint8),
        radius=2,
        fp_size=1,
        packed=packed,
    )
    with pytest.raises(ValueError, match="unsupported database SMILES at row 0"):
        builder.verify_fingerprints(path, sample_size=1)


@pytest.mark.parametrize("version", ["v2000", "v3000"])
def test_fixed_sdf_unconsumed_wedging_rejected_before_bits_or_serialization(version):
    from pathlib import Path

    path = Path(__file__).parent / "synthesis/fixtures" / ("atrop-up-" + version + ".sdf")
    mol = list(Chem.SDMolSupplier(str(path), removeHs=False))[0]
    before = Chem.MolToCXSmiles(mol)
    for operation in [fp.stable_canonical_molecule, fp.generate_query_fingerprint]:
        with pytest.raises(ValueError, match="atropisomeric|wedging"):
            operation(mol)
    assert Chem.MolToCXSmiles(mol) == before


@pytest.mark.parametrize("text", SUPPORTED + ["[13CH3:8][C@H:2](F)Cl |a:1|"])
def test_supported_raw_abs_and_mapped_states_keep_bit_values_and_caller(text):
    mol = parse_absolute_smiles(text)
    before = Chem.MolToCXSmiles(mol)
    expected_bits = Chem.rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=33).GetFingerprintAsNumPy(
        mol
    )
    np.testing.assert_array_equal(fp.generate_query_fingerprint(mol, 2, 33), expected_bits)
    exact, normalized, rounds = fp.stable_canonical_molecule(mol)
    assert exact == canonical(text)
    assert Chem.MolToSmiles(normalized) == exact and rounds >= 1
    assert Chem.MolToCXSmiles(mol) == before
