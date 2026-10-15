"""Search-policy decisions and bounded exploration independent of 3D claims."""

import json

import pytest
from rdkit import Chem

import chembounce as cb
import chemistry
from chembounce_synthesis.exploration import hash_reservoir, select_diverse
from chembounce_synthesis.identity import scaffold_family_key


def test_hash_sample_order_independent_and_bounded_with_duplicates():
    values = ["CC", "CCC", "CCCC", "CO", "CCN", "c1ccccc1"]
    first = hash_reservoir(values + values, 3, smiles=lambda s: s, seed=42)
    assert first == hash_reservoir(list(reversed(values)), 3, smiles=lambda s: s, seed=42)
    assert len(set(first)) == 3


def test_reserved_lane_retains_lower_2d_scaffold_despite_higher_scoring_decoys():
    values = ["Cc1ccncc1", "CCc1ccncc1", "CCCc1ccncc1", "C1CCOC1"]
    scores = dict(zip(values, [0.9, 0.8, 0.7, 0.1]))
    fp = lambda s: chemistry.count_fingerprint(Chem.MolFromSmiles(s))
    chosen, report = select_diverse(
        values,
        2,
        score=scores.get,
        smiles=lambda s: s,
        fingerprint=fp,
        threshold=0.5,
        fraction=0.5,
        group=scaffold_family_key,
    )
    assert set(chosen) == {values[0], values[-1]}
    assert report["below_threshold_selected"] == 1
    assert chosen[0] == values[-1]  # First-N downstream caps do not erase the reserved lane.


def test_zero_exploration_is_an_explicit_ablation():
    chosen, report = select_diverse(
        ["CC", "CCC"],
        1,
        score=lambda s: len(s) / 3,
        smiles=lambda s: s,
        fingerprint=lambda s: chemistry.count_fingerprint(Chem.MolFromSmiles(s)),
        fraction=0,
    )
    assert chosen == ["CCC"]
    assert report["exploration_selected"] == 0


def test_diverse_pipeline_preserves_constraints_but_does_not_hard_drop_low2d(tmp_path):
    library = ["c1ccncc1", "C1CCCCC1", "c1ccoc1"]
    args = dict(
        fragments_DB=library,
        use_fingerprint_db=False,
        fragments=["c1ccccc1"],
        candidate_thresholds={},
        calculate_shape=False,
        overall_max_n=2,
        frag_max_n=2,
        scaffold_top_n=3,
        cand_max_n__rplc=2,
        _merge_structure_top_n_=8,
        tanimoto_threshold=0.9,
        _search_scf_thr_=0.9,
    )
    old, _ = cb.chembounce("CCc1ccccc1", output_dir=tmp_path / "legacy", **args)
    assert old.empty
    rows, _ = cb.chembounce(
        "CCc1ccccc1",
        output_dir=tmp_path / "diverse",
        search_policy="diverse",
        exploration_fraction=0.5,
        search_candidate_budget=6,
        core_smiles="CC",
        **args,
    )
    assert len(rows) == 2
    assert rows["Tanimoto Similarity"].lt(0.9).all()
    assert all(
        Chem.MolFromSmiles(s).HasSubstructMatch(Chem.MolFromSmiles("CC")) for s in rows["Final structure"]
    )
    report = json.loads((tmp_path / "diverse/run_summary.json").read_text())
    assert report["search_budget"]["candidate_requests"] <= 6
    assert report["counts"]["generated"] <= 6
    assert report["candidate_selection"]["below_threshold_selected"] == 2
    decisions = json.loads((tmp_path / "diverse/search_decisions.json").read_text())["decisions"]
    assert len([r for r in decisions if r["status"] == "selected"]) == 2
    assert all(r["lane"] == "exploration" for r in decisions if r["status"] == "selected")
    manifest = json.loads((tmp_path / "diverse/run_manifest.json").read_text())
    assert manifest["parameters"]["search_policy"] == "diverse"
    assert manifest["parameters"]["search_candidate_budget"] == 6


def test_failed_attachments_still_charge_global_search_budget(tmp_path, monkeypatch):
    calls = []

    def impossible(*args, **kwargs):
        calls.append(args[3])
        return []

    monkeypatch.setattr(cb, "replace_molecule", impossible)
    cb.chembounce(
        "CCc1ccccc1",
        output_dir=tmp_path,
        fragments_DB=["c1ccncc1", "c1ccoc1", "C1CCCCC1"],
        use_fingerprint_db=False,
        candidate_thresholds={},
        calculate_shape=False,
        search_policy="diverse",
        search_candidate_budget=3,
        cand_max_n__rplc=2,
        scaffold_top_n=3,
        overall_max_n=2,
    )
    assert calls == [2, 1]
    report = json.loads((tmp_path / "run_summary.json").read_text())
    assert report["search_budget"]["candidate_requests"] == 3
    assert report["search_budget"]["unvisited_replacements"] == 1


def test_diverse_beam_keeps_alternative_attachment_instead_of_only_top2d():
    target = Chem.MolFromSmiles("COc1cc(Cl)ccc1N")
    pattern, replacement = Chem.MolFromSmiles("c1ccccc1"), Chem.MolFromSmiles("c1ccncc1")
    legacy = chemistry.replace(target, pattern, replacement, 2, 2)
    report = {}
    diverse = chemistry.replace(
        target,
        pattern,
        replacement,
        2,
        2,
        search_policy="diverse",
        exploration_fraction=0.5,
        seed=42,
        diagnostics=report,
    )
    assert len(diverse) <= 2
    assert diverse == chemistry.replace(
        target, pattern, replacement, 2, 2, search_policy="diverse", exploration_fraction=0.5, seed=42
    )
    assert {s for _, s in diverse} != {s for _, s in legacy}
    assert any(stage["exploration_selected"] for stage in report["beam_stages"])
    assert all(stage["selected"] <= 2 for stage in report["beam_stages"])
    assert max(score for score, _ in diverse) <= max(score for score, _ in legacy)


@pytest.mark.parametrize(
    "option",
    [
        {"search_policy": "invented"},
        {"exploration_fraction": -1},
        {"exploration_fraction": float("nan")},
        {"search_candidate_budget": 0},
        {"search_candidate_budget": True},
    ],
)
def test_invalid_research_settings_fail_before_search(tmp_path, option):
    with pytest.raises(ValueError):
        cb.chembounce("CCc1ccccc1", output_dir=tmp_path, **option)


@pytest.mark.parametrize("database_mode", ["memory", "text", "stream", "fingerprints"])
def test_diverse_scaffold_selection_agrees_across_library_backends(tmp_path, database_mode):
    from generate_fingerprints import generate_fingerprints

    library = ["c1ccncc1", "C1CCCCC1", "c1ccoc1", "c1cc[nH]c1", "C1CCNCC1"]
    source = tmp_path / "scaffolds.smi"
    source.write_text("\n".join(reversed(library)) + "\n")
    options = {"use_fingerprint_db": False}
    if database_mode == "memory":
        options["fragments_DB"] = library
    elif database_mode in {"text", "stream"}:
        options.update(scaffold_db_file=source, low_mem=database_mode == "stream")
    else:
        database = tmp_path / "fingerprints.npz"
        generate_fingerprints(source, database)
        options.update(use_fingerprint_db=True, fingerprint_db_file=database)
    out = tmp_path / "out"
    cb.chembounce(
        "CCc1ccccc1",
        output_dir=out,
        fragments=["c1ccccc1"],
        search_policy="diverse",
        exploration_fraction=0.5,
        search_candidate_budget=2,
        scaffold_top_n=3,
        overall_max_n=2,
        candidate_thresholds={},
        calculate_shape=False,
        **options,
    )
    selected = cb._read_replacements(out / "replace_scaffold_list.fragment_0001.tsv", preserve_order=True)
    expected, _ = cb._search_diverse_scaffolds(Chem.MolFromSmiles("c1ccccc1"), library, 3, 0.3, 42, 0.5)
    assert selected.index.tolist() == expected.index.tolist()
    assert selected.tolist() == pytest.approx(expected.tolist())


def test_zero_exploration_never_admits_below_threshold_on_unused_priority_slots():
    chosen, report = select_diverse(
        ["CC", "CCC"],
        2,
        score=lambda s: 0.9 if s == "CC" else 0.1,
        smiles=lambda s: s,
        fingerprint=lambda s: chemistry.count_fingerprint(Chem.MolFromSmiles(s)),
        threshold=0.5,
        fraction=0,
    )
    assert chosen == ["CC"]
    assert report["below_threshold_selected"] == 0


def test_library_identity_failure_excludes_only_that_candidate(monkeypatch):
    import chembounce_synthesis.identity as identity

    monkeypatch.setattr(identity, "TAUTOMER_LIMIT", 2)
    unstable = "O=c1[nH]c(=O)[nH]c(=O)[nH]1"
    with pytest.raises(ValueError, match="incomplete"):
        identity.chemical_identity(unstable)
    chosen, report = cb._search_diverse_scaffolds(
        Chem.MolFromSmiles("c1ccccc1"), [unstable, "c1ccncc1"], 2, 0.3, 42, 0.5
    )
    assert chosen.index.tolist() == ["c1ccncc1"]
    assert report["identity_or_structure_rejected"] == 1
    assert "incomplete" in report["rejection_examples"][0]["reason"]
    with pytest.raises(ValueError, match="Query scaffold identity failure"):
        cb._search_diverse_scaffolds(Chem.MolFromSmiles(unstable), ["c1ccncc1"], 2, 0.3, 42, 0.5)


def test_generated_identity_failure_retains_valid_other_candidate(tmp_path, monkeypatch):
    import chembounce_synthesis.identity as identity

    monkeypatch.setattr(identity, "TAUTOMER_LIMIT", 2)
    unstable = "O=c1[nH]c(=O)[nH]c(=O)[nH]1"
    monkeypatch.setattr(
        cb, "replace_molecule", lambda *args, **kwargs: [(0.1, unstable), (0.7, "CCc1ccncc1")]
    )
    rows, _ = cb.chembounce(
        "CCc1ccccc1",
        output_dir=tmp_path,
        fragments_DB=["c1ccncc1"],
        use_fingerprint_db=False,
        fragments=["c1ccccc1"],
        search_policy="diverse",
        search_candidate_budget=2,
        cand_max_n__rplc=2,
        scaffold_top_n=1,
        overall_max_n=2,
        calculate_shape=False,
        candidate_thresholds={},
    )
    assert rows["Final structure"].tolist() == ["CCc1ccncc1"]
    report = json.loads((tmp_path / "run_summary.json").read_text())
    assert report["counts"]["identity_rejected"] == 1
    decisions = json.loads((tmp_path / "search_decisions.json").read_text())["decisions"]
    assert decisions[0]["status"] == "identity_rejected"
    assert "incomplete" in decisions[0]["reason"]
