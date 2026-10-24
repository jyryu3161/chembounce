"""Synthetic helper controls only; no saved source pool or reaction run."""

import time

import pytest
from rdkit import Chem, DataStructs
from rdkit.Chem import rdChemReactions

from chembounce_synthesis.route_first_marker import Meter, UnknownMarker, departure_spec, marked_view, max_marked_similarity, source_views

RULE036 = rdChemReactions.ReactionFromSmarts(
    "[#6;$(C=C-[#6]),$(c:c):1][Br,I].[Cl,Br,I][c:2]>>[c:2][#6:1]"
)
BROAD = rdChemReactions.ReactionFromSmarts("[#6:1][Cl].[Br][c:2]>>[c:2][#6:1]")


def meter(**kwargs):
    return Meter(time.monotonic() + 5, **kwargs)


def test_fixed_rule_template_and_cl_br_marker_invariance_with_spectator():
    assert departure_spec(RULE036, 0) == (1, 0)
    assert departure_spec(RULE036, 1) == (0, 1)
    m = meter()
    cl = Chem.MolFromSmiles("Clc1cc(F)ccc1")
    br = Chem.MolFromSmiles("Brc1cc(F)ccc1")
    cl_before, br_before = Chem.MolToMolBlock(cl), Chem.MolToMolBlock(br)
    cl_sites, cl_views = source_views(cl, RULE036, 1, m)
    br_sites, br_views = source_views(br, RULE036, 1, m)
    assert len(cl_sites) == len(br_sites) == len(cl_views) == len(br_views) == 1
    assert set(cl_views) == set(br_views)
    assert next(iter(cl_views.values())) == next(iter(br_views.values()))
    assert "F" in next(iter(cl_views))
    assert (Chem.MolToMolBlock(cl), Chem.MolToMolBlock(br)) == (cl_before, br_before)

    spectator = Chem.MolFromSmiles("Brc1cc(Cl)ccc1")
    _, views = source_views(spectator, RULE036, 0, m)
    assert len(views) == 1 and "Cl" in next(iter(views))


def test_multisite_max_records_all_sites_without_claiming_reaction_site():
    source = Chem.MolFromSmiles("Clc1cc(OC)ccc1Cl")
    source_before = Chem.MolToMolBlock(source)
    m = meter()
    sites, views = source_views(source, RULE036, 1, m)
    assert len(sites) == len(views) == 2
    guide = m.fingerprint(marked_view(source, sites[0]["departing_atom_index"], sites[0]["attachment_atom_index"], m))
    assert max_marked_similarity(guide, views, m) == 1.0
    assert m.tanimoto_attempts == 2
    assert min(DataStructs.TanimotoSimilarity(guide, fp) for fp in views.values()) < 1.0
    assert len({item["departing_atom_index"] for item in sites}) == 2
    assert Chem.MolToMolBlock(source) == source_before

    # Symmetry can yield two possible sites yet only one projected view. The
    # score retains multiplicity and cannot prove which site a rule would use.
    symmetric = Chem.MolFromSmiles("Clc1ccc(Cl)cc1")
    symmetric_sites, symmetric_views = source_views(symmetric, RULE036, 1, m)
    assert len(symmetric_sites) == 2 and len(symmetric_views) == 1


def test_sentinel_and_unsafe_second_site_make_whole_source_unknown():
    many = Chem.MolFromSmiles("ClC(Cl)" + "C(Cl)" * 65 + "C")
    assert many is not None and many.GetNumAtoms() <= 256
    with pytest.raises(UnknownMarker, match="65_match_sentinel"):
        source_views(many, BROAD, 0, meter())

    mixed = Chem.MolFromSmiles("Clc1ccc([C@H](Cl)F)cc1")
    before = Chem.MolToMolBlock(mixed)
    charged = meter()
    with pytest.raises(UnknownMarker, match="unsafe_concrete_departure"):
        source_views(mixed, BROAD, 0, charged)
    assert charged.fingerprints >= 1  # First safe view was charged, never returned as a partial score.
    assert Chem.MolToMolBlock(mixed) == before


def test_nonterminal_stereo_charge_and_isotope_reject():
    nonterminal = rdChemReactions.ReactionFromSmarts("[C:1]O[C:3].[Cl][c:2]>>[c:2][C:1]")
    with pytest.raises(UnknownMarker, match="unproved_single_halogen_template"):
        departure_spec(nonterminal, 0)
    element_change = rdChemReactions.ReactionFromSmarts(
        "[C:1][Cl].[Br][c:2]>>[c:2][N:1]"
    )
    with pytest.raises(UnknownMarker, match="endpoint_element_not_proven_same"):
        departure_spec(element_change, 0)
    stereo = Chem.MolFromSmiles("Cl[C@H](F)C")
    with pytest.raises(UnknownMarker, match="unsafe_concrete_departure"):
        marked_view(stereo, 0, 1, meter())
    isotope = Chem.MolFromSmiles("[37Cl]c1ccccc1")
    with pytest.raises(UnknownMarker, match="unsafe_concrete_departure"):
        marked_view(isotope, 0, 1, meter())
    charged = Chem.MolFromSmiles("Clc1ccccc1")
    charged.GetAtomWithIdx(0).SetFormalCharge(1)
    with pytest.raises(UnknownMarker, match="unsafe_concrete_departure"):
        marked_view(charged, 0, 1, meter())


def test_budget_failure_invalidates_complete_score_and_keeps_original():
    source = Chem.MolFromSmiles("Clc1cc(OC)ccc1Cl")
    original = Chem.MolToMolBlock(source)
    m = meter(fp_limit=1)
    with pytest.raises(UnknownMarker, match="fingerprint_cap"):
        source_views(source, RULE036, 1, m)
    assert m.fingerprints == 1
    with pytest.raises(UnknownMarker, match="match_request_cap"):
        source_views(source, RULE036, 1, meter(match_limit=0))
    with pytest.raises(UnknownMarker, match="deadline"):
        source_views(source, RULE036, 1, Meter(time.monotonic() - 1))
    scores = meter(score_limit=0)
    with pytest.raises(UnknownMarker, match="tanimoto_cap"):
        max_marked_similarity(meter().fingerprint(source), {"source": meter().fingerprint(source)}, scores)
    assert scores.tanimoto_attempts == 0
    assert Chem.MolToMolBlock(source) == original
