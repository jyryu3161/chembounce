"""Synthetic-only selection and fallback controls; no saved source rescoring."""

import json
import time
from types import SimpleNamespace

import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from chembounce_synthesis import route_first_marker_second as selector
from chembounce_synthesis.catalog import BuildingBlock
from chembounce_synthesis.common import molecule
from chembounce_synthesis.reactions import Rules
from chembounce_synthesis.route_first import RouteFirstLimits, RouteFirstSearch


RULE = rdChemReactions.ReactionFromSmarts(
    "[#6;$(C=C-[#6]),$(c:c):1][Br,I].[Cl,Br,I][c:2]>>[c:2][#6:1]"
)
FIRST = SimpleNamespace(id="first-id", smiles="reserved-first")
INCUMBENT = SimpleNamespace(id="incumbent-id", smiles="Clc1ccccc1")
FAVORED = SimpleNamespace(id="favored-id", smiles="Brc1cc(F)ccc1")
ORDINARY = ((-0.5, 0, INCUMBENT.smiles), (-0.4, 0, FAVORED.smiles))
MOLS = {block.smiles: Chem.MolFromSmiles(block.smiles) for block in (INCUMBENT, FAVORED)}
BLOCKS = {block.smiles: block for block in (INCUMBENT, FAVORED)}


def setup(monkeypatch, *, reserve=None, ranked=ORDINARY, meter=None):
    counters = {"reserved_fp": 0, "eligibility": 0}

    def charge():
        counters["reserved_fp"] += 1
        return True if reserve is None else reserve(counters["reserved_fp"])

    meter = meter or selector.ChargedMeter(time.monotonic() + 5, charge, lambda: True)
    chooser = selector.Selector(meter)
    assert chooser.static_departure_supported(RULE, "rule", 1)
    guide = Chem.MolFromSmiles("Clc1cc(F)ccc1")
    guide_view = selector.marker.marked_view(guide, 0, 1, meter)
    guide_fp = meter.fingerprint(guide_view)
    monkeypatch.setattr(selector, "proven_guide_fingerprint", lambda *_: guide_fp)

    def eligibility(_smiles):
        counters["eligibility"] += 1
        return True

    args = dict(context_key=(("target", "role", False), "rule", "target", ("guide",), 1, "outer", 0),
                legacy_pair=(FIRST, INCUMBENT), ordinary_rows=ranked,
                pool_complete=True, registry={"by_guide": {}}, guide=("guide",),
                reaction=RULE, rule_id="rule", source_slot=1, tags={}, molecules=MOLS, blocks=BLOCKS,
                eligible=eligibility)
    return chooser, args, counters


def test_complementarity_changes_only_second_and_preserves_source_mols(monkeypatch):
    before = {smiles: Chem.MolToMolBlock(mol) for smiles, mol in MOLS.items()}
    chooser, args, counters = setup(monkeypatch)
    pair, receipt = chooser.choose(**args)
    assert pair == (FIRST, FAVORED)
    assert receipt["status"] == "complete_replacement"
    assert [row[0] for row in receipt["frozen_candidates"]] == [INCUMBENT.smiles, FAVORED.smiles]
    assert receipt["scores"][1]["marked_tanimoto"] > receipt["scores"][0]["marked_tanimoto"]
    assert receipt["scores"][1]["site_count"] == len(receipt["scores"][1]["site_to_view"])
    assert receipt["scores"][1]["distinct_projected_views"] >= 1
    assert receipt["scores"][1]["chosen_reaction_site"] is None
    assert counters["reserved_fp"] == chooser.meter.fingerprints
    assert all(Chem.MolToMolBlock(MOLS[k]) == value for k, value in before.items())


def test_incumbent_wins_equal_score_and_cache_is_separate(monkeypatch):
    chooser, args, counters = setup(monkeypatch)
    twin = SimpleNamespace(id="twin-id", smiles="Brc1ccccc1")
    args["ordinary_rows"] = ((-0.5, 0, INCUMBENT.smiles), (-0.4, 0, twin.smiles))
    args["molecules"] = {**MOLS, twin.smiles: Chem.MolFromSmiles(twin.smiles)}
    args["blocks"] = {**BLOCKS, twin.smiles: twin}
    pair, receipt = chooser.choose(**args)
    assert pair == args["legacy_pair"]
    assert receipt["status"] == "complete_noop"
    prior = (counters["reserved_fp"], counters["eligibility"])
    assert chooser.choose(**args)[0] == pair
    assert prior == (counters["reserved_fp"], counters["eligibility"])
    assert len(chooser.selected) == 1


def test_same_smiles_different_catalog_id_cannot_reuse_completed_pair(monkeypatch):
    chooser, args, _ = setup(monkeypatch)
    _, first_receipt = chooser.choose(**args)
    assert first_receipt["status"] == "complete_replacement"
    alternate = SimpleNamespace(id="different-incumbent-id", smiles=INCUMBENT.smiles)
    args["legacy_pair"] = (FIRST, alternate)
    _, second_receipt = chooser.choose(**args)
    assert second_receipt["status"] == "complete_replacement"
    assert len(chooser.selected) == 2


@pytest.mark.parametrize("unknown", ["eligibility", "budget", "registry", "truncation", "expiry"])
def test_unknown_atomically_restores_incoming_pair(monkeypatch, unknown):
    reserve = (lambda count: count <= 1) if unknown == "budget" else None
    chooser, args, _ = setup(monkeypatch, reserve=reserve)
    if unknown == "eligibility":
        args["eligible"] = lambda smiles: None if smiles == FAVORED.smiles else True
    elif unknown == "registry":
        monkeypatch.undo()
        args["registry"] = {"incomplete": "snapshot_budget", "by_guide": {}}
    elif unknown == "truncation":
        args["ordinary_rows"] = args["ordinary_rows"] * 4097
    elif unknown == "expiry":
        chooser.meter.deadline = time.monotonic() - 1
    pair, receipt = chooser.choose(**args)
    assert pair == args["legacy_pair"]
    assert receipt["status"] == "legacy_unknown"
    assert not chooser.selected


def test_definite_negative_excluded_and_first_fifteen_only(monkeypatch):
    chooser, args, counters = setup(monkeypatch)
    args["eligible"] = lambda smiles: smiles != FAVORED.smiles
    pair, receipt = chooser.choose(**args)
    assert pair == args["legacy_pair"]
    assert len(receipt["scores"]) == 1
    assert chooser.candidate_attempts == 2
    assert counters["reserved_fp"] == chooser.meter.fingerprints

    many = tuple((-0.5 + index / 1000, 0, f"source-{index}") for index in range(30))
    frozen = selector.frozen_candidates("reserved-first", "source-0", many)
    assert len(frozen) == 16 and frozen[0][0] == "source-0"
    assert frozen[-1][0] == "source-15"


def test_existing_tag_metadata_must_agree_without_compile():
    tags = {(1, 0): 3}
    assert selector.existing_tags({("rule", False): (object(), tags),
                                   ("rule", True): (object(), dict(tags))}, "rule") == tags
    with pytest.raises(selector.UnknownMarker, match="missing_or_conflicting"):
        selector.existing_tags({("rule", False): (object(), tags),
                                ("rule", True): (object(), {(1, 0): 9})}, "rule")


def test_exact_guide_alternatives_all_required_and_consistent(monkeypatch):
    clean = Chem.MolFromSmiles("Clc1cc(F)ccc1")
    other = Chem.MolFromSmiles("Brc1cc(F)ccc1")
    guide = (Chem.MolToSmiles(clean, isomericSmiles=True),)
    registry = {"incomplete": None, "raw_cap_saturated": False,
                "by_guide": {guide: [((clean,), (clean,)), ((clean,), (clean,))]}}
    meter = selector.ChargedMeter(time.monotonic() + 5, lambda: True, lambda: True)
    monkeypatch.setattr(selector.marker, "guide_view", lambda raw, *_args: raw)
    selector.proven_guide_fingerprint(registry, guide, RULE, 0, {}, meter)
    assert meter.fingerprints == 1  # Duplicate exact alternatives share one projected view.
    registry["by_guide"][guide][1] = ((other,), (clean,))
    with pytest.raises(selector.UnknownMarker, match="projected_view_conflict"):
        selector.proven_guide_fingerprint(registry, guide, RULE, 0, {}, meter)
    registry["raw_cap_saturated"] = True
    with pytest.raises(selector.UnknownMarker, match="raw_outcome_cap_saturated"):
        selector.proven_guide_fingerprint(registry, guide, RULE, 0, {}, meter)


@pytest.mark.parametrize("missing", ["incomplete", "raw_cap_saturated", "by_guide"])
def test_missing_raw_registry_completeness_fields_fail_closed(monkeypatch, missing):
    clean = Chem.MolFromSmiles("Clc1cc(F)ccc1")
    guide = (Chem.MolToSmiles(clean, isomericSmiles=True),)
    registry = {"incomplete": None, "raw_cap_saturated": False,
                "by_guide": {guide: [((clean,), (clean,))]}}
    del registry[missing]
    meter = selector.ChargedMeter(time.monotonic() + 5, lambda: True, lambda: True)
    monkeypatch.setattr(selector.marker, "guide_view", lambda raw, *_args: raw)
    with pytest.raises(selector.UnknownMarker):
        selector.proven_guide_fingerprint(registry, guide, RULE, 0, {}, meter)
    assert meter.fingerprints == 0


def test_any_unsafe_source_view_restores_pair_without_partial_score(monkeypatch):
    chooser, args, _ = setup(monkeypatch)
    original = selector.marker.source_views

    def source_views(mol, *rest):
        if Chem.MolToSmiles(mol) == Chem.MolToSmiles(MOLS[FAVORED.smiles]):
            raise selector.UnknownMarker("unsafe_second_view")
        return original(mol, *rest)

    monkeypatch.setattr(selector.marker, "source_views", source_views)
    pair, receipt = chooser.choose(**args)
    assert pair == args["legacy_pair"]
    assert receipt["status"] == "legacy_unknown"
    assert chooser.meter.fingerprints > 1  # Incumbent work was charged but cannot be used partially.


def test_failed_contexts_are_bounded_and_metadata_visits_charged(monkeypatch):
    chooser, args, _ = setup(monkeypatch)
    args["eligible"] = lambda _smiles: None
    for index in range(selector.MAX_CONTEXTS):
        args["context_key"] = ("distinct", index)
        pair, receipt = chooser.choose(**args)
        assert pair == args["legacy_pair"] and receipt["status"] == "legacy_unknown"
        assert receipt["failed_stage"] == "source_eligibility"
        assert receipt["failed_candidate"] == INCUMBENT.smiles
    assert chooser.attempted_contexts == selector.MAX_CONTEXTS
    assert chooser.ranked_metadata_rows == 2 * selector.MAX_CONTEXTS
    args["context_key"] = ("distinct", selector.MAX_CONTEXTS)
    _, receipt = chooser.choose(**args)
    assert receipt["reason"] == "context_cap"
    assert chooser.ranked_metadata_rows == 2 * selector.MAX_CONTEXTS


def test_shared_fp_reservation_survives_post_reservation_deadline():
    clock = [time.monotonic() + 5]
    meter = None

    def reserve():
        meter.deadline = time.monotonic() - 1
        return True

    meter = selector.ChargedMeter(clock[0], reserve, lambda: True)
    with pytest.raises(selector.UnknownMarker, match="deadline"):
        meter.fingerprint(Chem.MolFromSmiles("C"))
    assert meter.global_fingerprint_reservations == 1
    assert meter.fingerprints == 0


def test_shared_and_local_role_caps_prevent_native_match_request():
    source = Chem.MolFromSmiles("Clc1cc(F)ccc1")
    attempts = []

    def deny():
        attempts.append("shared")
        return False

    meter = selector.ChargedMeter(time.monotonic() + 5, lambda: True, deny)
    with pytest.raises(selector.UnknownMarker, match="shared_role_match_budget"):
        selector.marker.source_views(source, RULE, 1, meter)
    assert attempts == ["shared"]
    assert meter.match_requests == 0 and meter.shared_role_match_reservations == 0
    meter = selector.ChargedMeter(time.monotonic() + 5, lambda: True, deny)
    meter.match_limit = 0
    attempts.clear()
    with pytest.raises(selector.UnknownMarker, match="match_request_cap"):
        selector.marker.source_views(source, RULE, 1, meter)
    assert attempts == [] and meter.match_requests == 0


def test_unsupported_template_is_cached_before_context_attempt(monkeypatch):
    chooser, args, _ = setup(monkeypatch)
    unsupported = rdChemReactions.ReactionFromSmarts(
        "[C:1]O[C:3].[Cl][c:2]>>[c:2][C:1]"
    )
    assert chooser.static_departure_supported(unsupported, "unsupported", 0) is False
    before = chooser.static_template_checks
    assert chooser.static_departure_supported(unsupported, "unsupported", 0) is False
    assert chooser.static_template_checks == before
    args.update(rule_id="unsupported", reaction=unsupported, source_slot=0)
    pair, receipt = chooser.choose(**args)
    assert pair == args["legacy_pair"] and receipt["reason"] == "static_departure_preflight_not_proven"
    assert chooser.attempted_contexts == 0


def production_fixture(tmp_path, *, marked, alternative="Brc1ccc(Cl)cc1"):
    """One small synthetic reverse context, with no loaded catalog or model."""
    path = tmp_path / "rules.json"
    path.write_text(json.dumps({"schema_version": 1, "rules": [{
        "id": "synthetic-coupling", "smarts": rdChemReactions.ReactionToSmarts(RULE),
        "source": "synthetic",
    }]}))
    rules = Rules(path)
    smiles = ("Brc1ccccc1", "Clc1ccccc1", alternative)
    blocks = {value: BuildingBlock(value, value, "synthetic", "2026") for value in smiles}
    catalog = SimpleNamespace(
        sha256="synthetic", mols={value: molecule(value) for value in smiles},
        blocks=blocks, by_smiles={value: {value} for value in smiles},
    )
    limits = RouteFirstLimits(prefer_continuation_second=True,
                              prefer_marked_continuation_second=marked)
    search = RouteFirstSearch("Clc1ccc(C=CC)cc1", catalog, rules, limits)
    search._phase = 2
    rule = rules.get("synthetic-coupling")
    guide = search._reverse(rule, search.hit, False, record_centers=True)[0]
    key = (guide[1], rule.roles[1], False)
    incoming = (blocks[smiles[0]], blocks[smiles[1]])
    search._nearest_cache[key] = incoming
    search._continuation_ranked[key] = (
        ((-0.5, 0, smiles[1]), (-0.4, 0, smiles[2])), True,
    )
    kwargs = {"rule_id": rule.id, "role_slot": 1,
              "continuation_context": (rule, 1, 1, rule.id, search.hit, guide)}
    return search, rule, guide, incoming, kwargs


def test_production_hook_opt_in_changes_only_second_with_exact_raw_and_shared_budgets(tmp_path):
    off, _, old_guide, old_pair, old_kwargs = production_fixture(tmp_path, marked=False)
    on, rule, guide, incoming, kwargs = production_fixture(tmp_path, marked=True)
    assert guide == old_guide
    assert off._nearest(old_guide[1], rule.roles[1], **old_kwargs) == old_pair
    assert "marked_continuation_second" not in off.report
    before_reverse = on.report["reverse_attempts"]
    selected = on._nearest(guide[1], rule.roles[1], **kwargs)
    assert selected[0] == incoming[0]
    assert selected[1].smiles == "Brc1ccc(Cl)cc1"
    assert on._nearest_cache[(guide[1], rule.roles[1], False)] == incoming
    assert on.report["reverse_attempts"] == before_reverse
    work = on.report["marked_continuation_second"]
    assert work["context_attempts"] == work["completed_contexts"] == 1
    assert work["shared_generation_fp_reservations"] == work["helper_fp_attempts"]
    assert work["site_match_requests"] == work["shared_site_role_reservations"]
    assert work["eligibility_role_checks"] + work["site_match_requests"] <= 1536
    assert work["events"][0]["receipt"]["scores"][1]["chosen_reaction_site"] is None
    spent = (on.report["fingerprint_calls"], on.report["context_role_checks_total"])
    assert on._nearest(guide[1], rule.roles[1], **kwargs) == selected
    assert spent == (on.report["fingerprint_calls"], on.report["context_role_checks_total"])


def test_production_hook_unsafe_view_missing_proof_and_role_cap_restore_pair(tmp_path):
    unsafe, rule, guide, incoming, kwargs = production_fixture(
        tmp_path, marked=True, alternative="[37Cl]c1ccc(Br)cc1"
    )
    assert unsafe._nearest(guide[1], rule.roles[1], **kwargs) == incoming
    assert unsafe.report["marked_continuation_second"]["completed_contexts"] == 0
    missing, rule, guide, incoming, kwargs = production_fixture(tmp_path, marked=True)
    del missing._continuation_raw[(rule.id, missing.hit, False)]["incomplete"]
    assert missing._nearest(guide[1], rule.roles[1], **kwargs) == incoming
    assert missing.report["marked_continuation_second"]["context_attempts"] == 0
    capped, rule, guide, incoming, kwargs = production_fixture(tmp_path, marked=True)
    capped._context_role_checks = 65536
    assert capped._nearest(guide[1], rule.roles[1], **kwargs) == incoming
    assert capped.report["context_role_checks_total"] == 0
    assert capped.report["marked_continuation_second"]["completed_contexts"] == 0


def test_marked_policy_requires_explicit_boolean_and_existing_continuation():
    assert RouteFirstLimits().prefer_marked_continuation_second is False
    assert RouteFirstLimits(prefer_continuation_second=True,
                            prefer_marked_continuation_second=True).prefer_marked_continuation_second
    with pytest.raises(ValueError, match="requires prefer_continuation_second"):
        RouteFirstLimits(prefer_marked_continuation_second=True)
    with pytest.raises(ValueError, match="must be a boolean"):
        RouteFirstLimits(prefer_continuation_second=True, prefer_marked_continuation_second=1)
