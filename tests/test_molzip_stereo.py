"""Synthetic molzip compatibility; empty internal ABS is not external authority."""

from types import SimpleNamespace

import pytest
from rdkit import Chem

import chemistry
from chembounce_synthesis.common import require_absolute_stereo


def groups(mol):
    return [
        (
            str(g.GetGroupType()),
            tuple(a.GetIdx() for a in g.GetAtoms()),
            tuple(b.GetIdx() for b in g.GetBonds()) if hasattr(g, "GetBonds") else (),
        )
        for g in mol.GetStereoGroups()
    ]


@pytest.mark.parametrize(
    "text",
    [
        "COc1ccc(N)cc1",
        "C[C@H](O)c1ccccc1 |a:1|",
        "[13CH3][C@H](O)c1ccc(N)cc1 |a:1|",
        "C/C=C/c1ccc(N)cc1",
    ],
)
def test_actual_multistage_assembly_preserves_state_and_removes_only_internal_empty_abs(text, monkeypatch):
    target = Chem.MolFromSmiles(text)
    pattern = Chem.MolFromSmiles("c1ccccc1")
    before = Chem.MolToCXSmiles(target)
    original = chemistry._attach
    observed = []

    def inspect(state, site, label):
        state_before = Chem.MolToCXSmiles(state)
        product = original(state, site, label)
        assert Chem.MolToCXSmiles(state) == state_before
        assert not any(chemistry._empty_absolute_group(g) for g in product.GetStereoGroups())
        require_absolute_stereo(product)
        observed.append((label, groups(product)))
        return product

    monkeypatch.setattr(chemistry, "_attach", inspect)
    results = chemistry.replace(target, pattern, pattern, 100, 100)
    assert Chem.MolToSmiles(target) in {s for _, s in results}
    assert observed
    assert Chem.MolToCXSmiles(target) == before
    if " |a:" in text:
        assert any(any(g[0] == "STEREO_ABSOLUTE" and g[1] for g in state) for _, state in observed)
    if "ccc(N)" in text:
        assert len({label for label, _ in observed}) >= 2


def empty_group_mol(kind):
    mol = Chem.RWMol(Chem.MolFromSmiles("CC"))
    try:
        empty = Chem.CreateStereoGroup(kind, mol, [])
    except ValueError as exc:
        # Newer runtimes prohibit this artifact at construction; that is a
        # positive capability outcome, not a skipped ordinary assembly test.
        assert str(exc) == "New StereoGroup must contain at least one atom or bond."
        return None
    mol.SetStereoGroups([empty])
    return mol.GetMol()


@pytest.mark.parametrize(
    "kind",
    [Chem.StereoGroupType.STEREO_ABSOLUTE, Chem.StereoGroupType.STEREO_OR, Chem.StereoGroupType.STEREO_AND],
)
def test_external_empty_group_still_rejected(kind):
    raw = empty_group_mol(kind)
    if raw is None:
        # Still execute the unchanged guard's empty-group decision on a tiny
        # group-protocol sentinel; no claim this is a real new-runtime Mol.
        group = SimpleNamespace(GetGroupType=lambda: kind, GetAtoms=lambda: [])
        with pytest.raises(ValueError, match="Unsupported"):
            require_absolute_stereo(SimpleNamespace(GetStereoGroups=lambda: [group]))
        return
    before = groups(raw)
    with pytest.raises(ValueError, match="Unsupported"):
        require_absolute_stereo(raw)
    with pytest.raises(ValueError, match="Unsupported"):
        chemistry.replace(raw, Chem.MolFromSmiles("C"), Chem.MolFromSmiles("N"))
    assert groups(raw) == before


@pytest.mark.parametrize(
    "kind,atoms,bonds,expected",
    [
        (Chem.StereoGroupType.STEREO_ABSOLUTE, [], [], True),
        (Chem.StereoGroupType.STEREO_ABSOLUTE, [object()], [], False),
        (Chem.StereoGroupType.STEREO_ABSOLUTE, [], [object()], False),
        (Chem.StereoGroupType.STEREO_OR, [], [], False),
        (Chem.StereoGroupType.STEREO_AND, [], [], False),
    ],
)
def test_empty_group_predicate_requires_abs_and_zero_atom_and_bond_members(kind, atoms, bonds, expected):
    group = SimpleNamespace(GetGroupType=lambda: kind, GetAtoms=lambda: atoms, GetBonds=lambda: bonds)
    assert chemistry._empty_absolute_group(group) is expected


@pytest.mark.parametrize("kind", [Chem.StereoGroupType.STEREO_OR, Chem.StereoGroupType.STEREO_AND])
def test_internal_cleanup_cannot_erase_nonabsolute_group(monkeypatch, kind):
    # A controlled molzip output tests cleanup itself, with the real shared guard.
    tag = "o1:1" if kind == Chem.StereoGroupType.STEREO_OR else "&1:1"
    result = Chem.MolFromSmiles("F[C@H](Cl)Br |" + tag + "|")
    monkeypatch.setattr(Chem, "molzip", lambda *args: result)
    state = Chem.MolFromSmiles("C.[1*]O")
    product = chemistry._attach(state, 0, 1)
    assert groups(product) == groups(result)
    with pytest.raises(ValueError, match="nonabsolute enhanced"):
        require_absolute_stereo(product)


def test_internal_empty_abs_cleanup_does_not_mutate_original_product(monkeypatch):
    raw = empty_group_mol(Chem.StereoGroupType.STEREO_ABSOLUTE)
    if raw is None:
        # Actual cleanup branch is exercised on the minimum runtime that can
        # produce this artifact; newer constructors reject it before assembly.
        return
    before = groups(raw)
    monkeypatch.setattr(Chem, "molzip", lambda *args: raw)
    product = chemistry._attach(Chem.MolFromSmiles("C.[1*]O"), 0, 1)
    assert groups(raw) == before and groups(product) == []
    assert Chem.MolToSmiles(product) == Chem.MolToSmiles(raw) == "CC"
    require_absolute_stereo(product)


@pytest.mark.parametrize(
    "kind",
    [Chem.StereoGroupType.STEREO_ABSOLUTE, Chem.StereoGroupType.STEREO_OR, Chem.StereoGroupType.STEREO_AND],
)
def test_cleanup_copies_nonempty_group_members_without_changing_ownership_or_state(monkeypatch, kind):
    raw = Chem.RWMol(Chem.MolFromSmiles("F[C@H](Cl)Br"))
    members = [Chem.CreateStereoGroup(kind, raw, [1])]
    try:
        members.insert(0, Chem.CreateStereoGroup(Chem.StereoGroupType.STEREO_ABSOLUTE, raw, []))
    except ValueError as exc:
        assert str(exc) == "New StereoGroup must contain at least one atom or bond."
    raw.SetStereoGroups(members)
    raw = raw.GetMol()
    original = groups(raw)
    monkeypatch.setattr(Chem, "molzip", lambda *args: raw)
    product = chemistry._attach(Chem.MolFromSmiles("C.[1*]O"), 0, 1)
    assert groups(raw) == original
    assert groups(product) == [(str(kind), (1,), ())]
    atom = product.GetStereoGroups()[0].GetAtoms()[0]
    assert atom.GetOwningMol().GetNumAtoms() == product.GetNumAtoms() == 4
    assert Chem.MolToSmiles(atom.GetOwningMol()) == Chem.MolToSmiles(raw)
    if kind == Chem.StereoGroupType.STEREO_ABSOLUTE:
        require_absolute_stereo(product)
    else:
        with pytest.raises(ValueError, match="nonabsolute enhanced"):
            require_absolute_stereo(product)
