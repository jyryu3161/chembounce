"""Bounded marked-site ordinary-second selection for route-first search.

The caller supplies already retained raw outcomes, completed ordinary ranks and
catalog molecules. This module never retrieves a new pool or runs a reaction.
"""

from rdkit import Chem

from . import route_first_marker as marker

UnknownMarker = marker.UnknownMarker

MAX_CONTEXTS = 16
MAX_CANDIDATE_ATTEMPTS = 256
MAX_CANDIDATES_PER_CONTEXT = 16
MAX_NEW_FINGERPRINTS = 1024
MAX_SITE_MATCH_REQUESTS = 1024
MAX_STRUCTURAL_UNITS = 262144
MAX_RANKED_METADATA_ROWS = MAX_CONTEXTS * 8192
MAX_STATIC_TEMPLATE_CHECKS = 128


class ChargedMeter(marker.Meter):
    """Reserve every projected FP in the caller's existing generation budget."""

    def __init__(self, deadline, reserve_generation_fingerprint, reserve_role_match):
        self.reserve_role_match = reserve_role_match
        self._match_requests = 0
        self.shared_role_match_reservations = 0
        super().__init__(deadline, fp_limit=MAX_NEW_FINGERPRINTS,
                         match_limit=MAX_SITE_MATCH_REQUESTS,
                         parse_limit=0, score_limit=MAX_NEW_FINGERPRINTS,
                         units_limit=MAX_STRUCTURAL_UNITS)
        self.reserve_generation_fingerprint = reserve_generation_fingerprint
        self.global_fingerprint_reservations = 0
        self.guide_clean_serialization_attempts = 0

    @property
    def match_requests(self):
        return self._match_requests

    @match_requests.setter
    def match_requests(self, value):
        difference = value - self._match_requests
        if difference not in (0, 1):
            raise UnknownMarker("role_match_meter_increment_invalid")
        if difference:
            # Frozen marker calls this setter before each native full-role
            # GetSubstructMatches request. Both the local and original search
            # budgets must reserve the request before chemistry work starts.
            if self.reserve_role_match() is not True:
                raise UnknownMarker("shared_role_match_budget")
            self.shared_role_match_reservations += 1
        self._match_requests = value

    def fingerprint(self, mol):
        self.check()
        if self.fingerprints >= self.fp_limit:
            raise UnknownMarker("fingerprint_cap")
        # A failed or late FP attempt remains charged by the caller. This hook
        # must check report FP + accept FP against the generation ceiling.
        if self.reserve_generation_fingerprint() is not True:
            raise UnknownMarker("shared_generation_fingerprint_ceiling")
        self.global_fingerprint_reservations += 1
        return super().fingerprint(mol)


def existing_tags(compiled_reverse_cache, rule_id):
    """Read only pre-existing tagged templates, never compile or reverse again."""
    found = [tags for (ident, _nonring), (_reverse, tags) in compiled_reverse_cache.items()
             if ident == rule_id]
    if not found or any(tags != found[0] for tags in found[1:]):
        raise UnknownMarker("tagged_reverse_tags_missing_or_conflicting")
    return found[0]


def proven_guide_fingerprint(registry, guide, reaction, slot, tags, meter):
    """Require every exact raw alternative to produce the same safe marked view."""
    if (not isinstance(registry, dict) or "incomplete" not in registry
            or registry["incomplete"] is not None):
        raise UnknownMarker("raw_registry_missing_or_incomplete")
    if registry.get("raw_cap_saturated") is not False:
        raise UnknownMarker("raw_outcome_cap_saturated")
    if not isinstance(registry.get("by_guide"), dict):
        raise UnknownMarker("exact_guide_registry_missing")
    alternatives = registry["by_guide"].get(guide, ())
    if not alternatives:
        raise UnknownMarker("exact_guide_raw_alternatives_missing")
    views = {}
    for raw_parts, clean_parts in alternatives:
        meter.check()
        if len(raw_parts) != len(guide) or len(clean_parts) != len(guide) or slot not in range(len(guide)):
            raise UnknownMarker("guide_precursor_count_or_slot_mismatch")
        for index, clean in enumerate(clean_parts):
            meter.charge_structure(clean)
            meter.guide_clean_serialization_attempts += 1
            if Chem.MolToSmiles(clean, isomericSmiles=True) != guide[index]:
                raise UnknownMarker("exact_emitted_clean_guide_mismatch")
        meter.charge_structure(*reaction.GetReactants(), *reaction.GetProducts())
        view = marker.guide_view(raw_parts[slot], clean_parts[slot], reaction, slot, tags, meter)
        meter.projected_serialization_attempts += 1
        smiles = Chem.MolToSmiles(view, isomericSmiles=True)
        views.setdefault(smiles, view)
    if len(views) != 1:
        raise UnknownMarker("raw_alternative_projected_view_conflict")
    meter.check()
    return meter.fingerprint(next(iter(views.values())))


def frozen_candidates(first_smiles, incumbent_smiles, ordinary_rows, *, charge_row=None):
    """Incumbent plus at most 15 other original ordinary ranks, before scores."""
    if len(ordinary_rows) > 8192:
        raise UnknownMarker("ordinary_ranked_pool_limit")
    ranks = {}
    for rank, row in enumerate(ordinary_rows, 1):
        if charge_row is not None:
            charge_row()
        if not isinstance(row, tuple) or len(row) != 3 or row[1] != 0 or not isinstance(row[2], str):
            raise UnknownMarker("ordinary_ranked_pool_malformed")
        ranks.setdefault(row[2], rank)
    if first_smiles == incumbent_smiles or incumbent_smiles not in ranks:
        raise UnknownMarker("ordinary_second_or_first_reservation_unproved")
    chosen = [(incumbent_smiles, ranks[incumbent_smiles], True)]
    for smiles, rank in ranks.items():
        if smiles not in (first_smiles, incumbent_smiles):
            chosen.append((smiles, rank, False))
            if len(chosen) == MAX_CANDIDATES_PER_CONTEXT:
                break
    return tuple(chosen)


class Selector:
    """A separate completed-context cache; unknown never caches a replacement."""

    def __init__(self, meter):
        self.meter = meter
        self.selected = {}
        self.candidate_attempts = 0
        self.attempted_contexts = 0
        self.completed_contexts = 0
        self.ranked_metadata_rows = 0
        self.static_template_checks = 0
        self.static_template_cache = {}

    def static_departure_supported(self, reaction, rule_id, slot):
        """Cache generic type applicability before spending a score context."""
        key = (rule_id, slot)
        if key in self.static_template_cache:
            return self.static_template_cache[key]
        self.meter.check()
        if self.static_template_checks >= MAX_STATIC_TEMPLATE_CHECKS:
            raise UnknownMarker("static_template_check_cap")
        self.static_template_checks += 1
        self.meter.charge_structure(*reaction.GetReactants(), *reaction.GetProducts())
        try:
            marker.departure_spec(reaction, slot)
        except UnknownMarker:
            self.static_template_cache[key] = False
        else:
            self.static_template_cache[key] = True
        return self.static_template_cache[key]

    def choose(self, *, context_key, legacy_pair, ordinary_rows, pool_complete,
               registry, guide, reaction, rule_id, source_slot, tags, molecules, blocks, eligible):
        """Return (pair, receipt). `eligible(smiles)` is True/False/None.

        It must use the unchanged full outer role and strict source stereo checks;
        None means unknown and atomically restores the incoming pair.
        """
        stage = "preflight"
        active_candidate = None

        def old(reason, candidates=()):
            return legacy_pair, {"status": "legacy_unknown", "reason": reason,
                                 "frozen_candidates": candidates, "scores": [],
                                 "failed_stage": stage, "failed_candidate": active_candidate,
                                 "candidate_attempts_total": self.candidate_attempts,
                                 "context_attempts_total": self.attempted_contexts,
                                 "ranked_metadata_rows_total": self.ranked_metadata_rows,
                                 "shared_fp_reservations_total": self.meter.global_fingerprint_reservations,
                                 "helper_fp_attempts_total": self.meter.fingerprints,
                                 "site_match_requests_total": self.meter.match_requests,
                                 "shared_role_reservations_total": self.meter.shared_role_match_reservations,
                                 "structural_units_total": self.meter.structural_units}

        try:
            self.meter.check()
            if len(legacy_pair) != 2 or not pool_complete:
                raise UnknownMarker("complete_ordinary_pool_or_top2_unavailable")
            first, incumbent = legacy_pair
            first_smiles, incumbent_smiles = first.smiles, incumbent.smiles
            if not all(isinstance(block.id, str) and block.id for block in legacy_pair):
                raise UnknownMarker("concrete_pair_catalog_ids_missing")
            key = (context_key, (first.id, first_smiles), (incumbent.id, incumbent_smiles))
            if key in self.selected:
                return self.selected[key]
            if self.attempted_contexts >= MAX_CONTEXTS:
                raise UnknownMarker("context_cap")
            if self.static_template_cache.get((rule_id, source_slot)) is not True:
                raise UnknownMarker("static_departure_preflight_not_proven")
            self.attempted_contexts += 1

            def charge_ranked_row():
                if self.ranked_metadata_rows >= MAX_RANKED_METADATA_ROWS:
                    raise UnknownMarker("ranked_metadata_row_cap")
                self.ranked_metadata_rows += 1

            stage = "freeze_original_ranks"
            candidates = frozen_candidates(first_smiles, incumbent_smiles, ordinary_rows,
                                           charge_row=charge_ranked_row)
            stage = "exact_raw_guide_proof"
            guide_fp = proven_guide_fingerprint(registry, guide, reaction, source_slot, tags, self.meter)
            scores = []
            for smiles, rank, is_incumbent in candidates:
                active_candidate = smiles
                stage = "source_eligibility"
                self.meter.check()
                if self.candidate_attempts >= MAX_CANDIDATE_ATTEMPTS:
                    raise UnknownMarker("candidate_attempt_cap")
                self.candidate_attempts += 1
                allowed = eligible(smiles)
                self.meter.check()
                if allowed is None or type(allowed) is not bool:
                    raise UnknownMarker("outer_role_or_stereo_unknown")
                if not allowed:
                    continue
                mol = molecules.get(smiles)
                if mol is None:
                    raise UnknownMarker("existing_catalog_molecule_missing")
                stage = "complete_source_site_projection"
                self.meter.charge_structure(*reaction.GetReactants(), *reaction.GetProducts())
                sites, views = marker.source_views(mol, reaction, source_slot, self.meter)
                stage = "complete_view_scoring"
                score = marker.max_marked_similarity(guide_fp, views, self.meter)
                scores.append((score, is_incumbent, rank, smiles, sites, len(views)))
            if not scores:
                raise UnknownMarker("no_definitely_eligible_source")
            self.meter.check()
            stage = "publish_completed_choice"
            best = min(scores, key=lambda row: (-row[0], not row[1], row[2], row[3]))
            if best[3] not in blocks:
                raise UnknownMarker("existing_catalog_block_missing")
            selected = legacy_pair if best[3] == incumbent_smiles else (first, blocks[best[3]])
            receipt = {"status": "complete_noop" if selected == legacy_pair else "complete_replacement",
                       "reason": None, "frozen_candidates": candidates,
                       "scores": [{"smiles": smiles, "ordinary_rank": rank,
                                   "incumbent": incumbent_flag, "marked_tanimoto": score,
                                   "site_count": len(sites), "distinct_projected_views": views_count,
                                   "site_to_view": sites, "chosen_reaction_site": None}
                                  for score, incumbent_flag, rank, smiles, sites, views_count in scores],
                       "chosen_smiles": best[3]}
            receipt.update(candidate_attempts_total=self.candidate_attempts,
                           context_attempts_total=self.attempted_contexts,
                           ranked_metadata_rows_total=self.ranked_metadata_rows,
                           shared_fp_reservations_total=self.meter.global_fingerprint_reservations,
                           helper_fp_attempts_total=self.meter.fingerprints,
                           site_match_requests_total=self.meter.match_requests,
                           shared_role_reservations_total=self.meter.shared_role_match_reservations,
                           structural_units_total=self.meter.structural_units)
            self.completed_contexts += 1
            self.selected[key] = (selected, receipt)
            return selected, receipt
        except (UnknownMarker, ValueError, RuntimeError, AttributeError, KeyError, IndexError) as exc:
            return old(str(exc), locals().get("candidates", ()))
