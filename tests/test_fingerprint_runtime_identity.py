"""Exact identities are canonical under the running RDKit, not an asset's release."""

import json
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem, rdBase

import chembounce
import fingerprint_db as fp
from generate_fingerprints import verify_fingerprints


# Literal entry from the RDKit 2026.09.1 mw250 asset, row 3167. RDKit 2025.09.6
# uses the opposite slash spelling for the same E/Z state, not its stereoisomer.
STORED_SLASH = r"C(=C\c1ccccc1)/CCCCc1ccccc1"
ALTERNATE_SLASH = r"C(=C/c1ccccc1)\CCCCc1ccccc1"


def database(strings, version, packed, *, collisions=False):
    strings = sorted(strings)
    size = 1 if collisions else 2048
    bits = np.stack([fp.generate_query_fingerprint(Chem.MolFromSmiles(s), 2, size) for s in strings])
    data = {
        "schema_version": 2,
        "radius": 2,
        "fp_size": size,
        "smiles": np.asarray(strings),
        "fingerprints": np.packbits(bits, axis=1) if packed else bits,
        "packed": packed,
        "canonical_sorted_unique": True,
    }
    if version is not None:
        data["rdkit_version"] = version
    return data


def write_database(path, data):
    if path.suffix == ".npz":
        np.savez(path, **data)
    else:
        path.mkdir(parents=True)
        np.save(path / "smiles.npy", data["smiles"])
        np.save(path / "fingerprints.npy", data["fingerprints"])
        (path / "metadata.json").write_text(
            json.dumps({k: v for k, v in data.items() if not isinstance(v, np.ndarray)})
        )


@pytest.mark.parametrize("packed", [False, True])
def test_missing_version_real_alternate_slash_identity_is_excluded(packed):
    query = Chem.MolFromSmiles(ALTERNATE_SLASH)
    canonical = Chem.MolToSmiles(query)
    assert Chem.MolToSmiles(Chem.MolFromSmiles(STORED_SLASH)) == canonical
    opposite = Chem.MolToSmiles(Chem.MolFromSmiles(r"C(=C/c1ccccc1)/CCCCc1ccccc1"))
    assert opposite != canonical
    data = database([STORED_SLASH, opposite, "c1ccncc1"], None, packed)
    actual = fp.search_similar_scaffolds_fast(query, data, 2, 0)
    assert canonical not in [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in actual.index]
    assert list(actual.index) == [opposite, "c1ccncc1"]
    # Achiral fingerprints collide for these different specified stereoisomers.
    assert actual.iloc[0] == 1.0


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("top_n", [3, None])
def test_missing_version_runtime_duplicates_and_query_do_not_consume_slots(packed, top_n):
    data = database(["C(C)C", "C(C)CC", "C(CO)", "CCC", "CCCC", "CCCCC", "CO"], None, packed, collisions=True)
    query = Chem.MolFromSmiles("CCO")
    actual = fp.search_similar_scaffolds_fast(query, data, top_n, 0)
    assert list(actual.index) == ["CCC", "CCCC", "CCCCC", "CO"][:top_n]
    assert actual.tolist() == [1.0] * len(actual)
    pd.testing.assert_series_equal(actual, fp.search_similar_scaffolds_fast(query, data, top_n, 0))


@pytest.mark.parametrize("version", [None, rdBase.rdkitVersion])
def test_mapped_query_exact_identity_uses_copy_and_preserves_stereo_isotope_charge(version):
    query = Chem.MolFromSmiles("[13CH3][C@H](O)[NH3+]")
    for index, atom in enumerate(query.GetAtoms(), 1):
        atom.SetAtomMapNum(index)
    before = Chem.MolToSmiles(query)
    original_maps = [a.GetAtomMapNum() for a in query.GetAtoms()]
    variants = ["[13CH3][C@H](O)[NH3+]", "[13CH3][C@@H](O)[NH3+]", "C[C@H](O)[NH3+]", "[13CH3][C@H](O)N"]
    variants = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in variants]
    data = database(variants, version, True)
    actual = fp.search_similar_scaffolds_fast(query, data, 4, 0)
    assert set(actual.index) == set(variants[1:])
    assert Chem.MolToSmiles(query) == before
    assert [a.GetAtomMapNum() for a in query.GetAtoms()] == original_maps


@pytest.mark.parametrize("packed", [False, True])
def test_matching_runtime_keeps_string_fast_path_without_parsing_hits(monkeypatch, packed):
    data = database(["CC", "CCC", "CCCC", "CCCCC"], rdBase.rdkitVersion, packed, collisions=True)
    data["canonicalization_policy"] = fp.CANONICALIZATION_POLICY
    query = Chem.MolFromSmiles("CC")
    parse = Mock(wraps=Chem.MolFromSmiles)
    monkeypatch.setattr(fp.Chem, "MolFromSmiles", parse)
    actual = fp.search_similar_scaffolds_fast(query, data, 2, 0)
    assert list(actual.index) == ["CCC", "CCCC"]
    # Only the query is stabilized. Certified stored hits need no reparsing.
    assert [call.args[0] for call in parse.call_args_list] == ["CC"]


def test_missing_version_database_only_parses_selected_hits(monkeypatch):
    data = database(["C" * n for n in range(1, 41)], None, False, collisions=True)
    query = Chem.MolFromSmiles("N")
    parse = Mock(wraps=Chem.MolFromSmiles)
    monkeypatch.setattr(fp.Chem, "MolFromSmiles", parse)
    actual = fp.search_similar_scaffolds_fast(query, data, 3, 0)
    assert list(actual.index) == ["C", "CC", "CCC"]
    # One query parse, then input parsing + fixed-point verification per hit.
    assert [call.args[0] for call in parse.call_args_list] == ["N", "C", "C", "CC", "CC", "CCC", "CCC"]


@pytest.mark.parametrize("packed", [False, True])
@pytest.mark.parametrize("canonical_claim", [False, True])
def test_known_version_mismatch_rejected_before_scoring(monkeypatch, packed, canonical_claim):
    data = database([STORED_SLASH, "c1ccncc1"], "other-rdkit-release", packed)
    data["canonical_sorted_unique"] = canonical_claim
    score = Mock(side_effect=AssertionError("Do not rank known cross-version fingerprints"))
    monkeypatch.setattr(fp, "calculate_bulk_tanimoto", score)
    with pytest.raises(ValueError, match="RDKit.*Regenerate.*new output path"):
        fp.search_similar_scaffolds_fast(Chem.MolFromSmiles(ALTERNATE_SLASH), data, 1, 0)
    score.assert_not_called()


@pytest.mark.parametrize("directory", [False, True])
def test_passing_sample_does_not_unlock_foreign_version_ranking(tmp_path, directory):
    data = database([STORED_SLASH, "c1ccncc1"], "other-rdkit-release", True)
    path = tmp_path / ("db" if directory else "db.npz")
    write_database(path, data)
    assert verify_fingerprints(path, sample_size=2) == 2
    loaded = fp.load_fingerprint_db(path)
    assert str(loaded["rdkit_version"]) == "other-rdkit-release"
    # The loader remains available for read-only cross-version verification.
    with pytest.raises(ValueError, match="RDKit.*Regenerate"):
        fp.search_similar_scaffolds_fast(Chem.MolFromSmiles(ALTERNATE_SLASH), loaded, 1, 0)


@pytest.mark.parametrize("available", ["runtime", "generic", "legacy"])
def test_default_prefers_runtime_specific_database_and_explicit_path_wins(tmp_path, monkeypatch, available):
    monkeypatch.setattr(fp, "__file__", str(tmp_path / "fingerprint_db.py"))
    root = tmp_path / "data"
    root.mkdir()
    paths = {
        "runtime": root / f"fingerprints_mw250_rdkit_{rdBase.rdkitVersion}",
        "generic": root / "fingerprints_mw250_v2",
        "legacy": root / "scaffold_fingerprints_mw250.npz",
    }
    data = database(["CC", "CCC"], rdBase.rdkitVersion, True)
    # Install the selected preference and all lower-precedence alternatives.
    names = list(paths)
    for name in names[names.index(available) :]:
        write_database(paths[name], data)
    assert fp.load_fingerprint_db()["path"] == str(paths[available].resolve())
    assert fp.load_fingerprint_db(paths["legacy"])["path"] == str(paths["legacy"].resolve())


@pytest.mark.parametrize("policy", ["legacy", "diverse"])
def test_generator_dispatch_receives_runtime_distinct_priority_hits(tmp_path, monkeypatch, policy):
    data = database([STORED_SLASH, "c1ccncc1", "c1ccccc1"], None, True)
    path = tmp_path / "db.npz"
    np.savez(path, **data)
    query = Chem.MolFromSmiles(ALTERNATE_SLASH)
    original = Chem.MolToSmiles(query)
    if policy == "diverse":
        actual_diverse = chembounce._search_diverse_scaffolds

        def inspect_priority(*args, **kwargs):
            priority = kwargs["priority"]
            assert len(priority) == 2
            assert original not in [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in priority.index]
            return actual_diverse(*args, **kwargs)

        monkeypatch.setattr(chembounce, "_search_diverse_scaffolds", inspect_priority)
    output = tmp_path / policy
    chembounce.get_frags_cands(
        original,
        query,
        fragments=[original],
        scaffold_top_n=2,
        _search_scf_thr_=0,
        fingerprint_db_file=str(path),
        output_dir=output,
        search_policy=policy,
    )
    hits = pd.read_csv(output / "replace_scaffold_list.fragment_0001.tsv", sep="\t")
    assert len(hits) == 2
    assert original not in [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in hits["SMILES"]]


@pytest.mark.parametrize("policy", ["legacy", "diverse"])
def test_generator_policies_reject_known_version_mismatch(tmp_path, monkeypatch, policy):
    path = tmp_path / "db.npz"
    write_database(path, database(["c1ccccc1", "c1ccncc1"], "other-rdkit-release", True))
    score = Mock(side_effect=AssertionError("No ranking across known different RDKit versions"))
    monkeypatch.setattr(fp, "calculate_bulk_tanimoto", score)
    with pytest.raises(ValueError, match="RDKit.*Regenerate"):
        chembounce.get_frags_cands(
            "c1ccccc1",
            Chem.MolFromSmiles("c1ccccc1"),
            fragments=["c1ccccc1"],
            scaffold_top_n=2,
            fingerprint_db_file=str(path),
            output_dir=tmp_path / policy,
            search_policy=policy,
        )
    score.assert_not_called()
