"""Opt-in, bounded route-first retrieval. Reverse templates are hints, never evidence."""

import hashlib
import heapq
import itertools
import math
import time
from dataclasses import asdict, dataclass

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdChemReactions

from .common import FP, absolute_smiles, canonical, molecule, positive
from .route_audit import audit_route, template_maps
from .routes import verify_route

VERSION = "catalog-rule-route-first-v2"


class HitStereoAnchors:
    """Bounded local hit anchors; compile only representable, unique source mappings."""

    def __init__(self, hit):
        self.source = molecule(hit)
        if self.source.GetNumAtoms() > 128:
            raise ValueError("Hit stereo anchor molecule exceeds bounded comparison size")
        tetra = [
            (atom.GetIdx(),)
            for atom in self.source.GetAtoms()
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        ]
        bonds = [
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            for bond in self.source.GetBonds()
            if bond.GetStereo() != Chem.BondStereo.STEREONONE
        ]
        if len(tetra) + len(bonds) > 16:
            raise ValueError("Too many specified hit stereo anchors")
        if any(
            self.source.GetBondBetweenAtoms(*pair).GetStereo()
            not in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ)
            for pair in bonds
        ):
            raise ValueError("Unsupported specified hit bond stereo")
        self.queries = [self._compile(roots) for roots in tetra + bonds]

    def _compile(self, roots):
        included, frontier = set(roots), set(roots)
        rings = self.source.GetRingInfo().AtomRings()
        for radius in range(1, 9):
            previous_included = set(included)
            frontier = {
                neighbor.GetIdx()
                for index in frontier
                for neighbor in self.source.GetAtomWithIdx(index).GetNeighbors()
            } - included
            included.update(frontier)
            for _ in range(len(rings) + 1):
                previous = len(included)
                for ring in rings:
                    if included.intersection(ring):
                        included.update(ring)
                if len(included) == previous:
                    break
            frontier = included - previous_included
            if len(included) > 64:
                break
            marked = Chem.Mol(self.source)
            for position, index in enumerate(roots):
                marked.GetAtomWithIdx(index).SetAtomMapNum(1000000 + position)
            fragment = Chem.MolFragmentToSmiles(
                marked, atomsToUse=sorted(included), canonical=True, isomericSmiles=True
            )
            query = Chem.MolFromSmiles(fragment)
            if query is None:
                continue
            positions = [
                next(
                    (atom.GetIdx() for atom in query.GetAtoms() if atom.GetAtomMapNum() == 1000000 + i), None
                )
                for i in range(len(roots))
            ]
            if None in positions:
                continue
            if len(roots) == 1:
                defined = query.GetAtomWithIdx(positions[0]).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
            else:
                bond = query.GetBondBetweenAtoms(*positions)
                defined = bond is not None and bond.GetStereo() in (
                    Chem.BondStereo.STEREOE,
                    Chem.BondStereo.STEREOZ,
                )
            if not defined:
                continue
            for position in positions:
                query.GetAtomWithIdx(position).SetAtomMapNum(0)
            matches = self.source.GetSubstructMatches(
                query, useChirality=False, uniquify=False, maxMatches=65
            )
            locations = {tuple(match[position] for position in positions) for match in matches}
            chiral_matches = self.source.GetSubstructMatches(
                query, useChirality=True, uniquify=False, maxMatches=65
            )
            chiral_locations = {tuple(match[position] for position in positions) for match in chiral_matches}
            if (
                len(matches) < 65
                and len(chiral_matches) < 65
                and locations == {tuple(roots)}
                and chiral_locations == {tuple(roots)}
            ):
                return query, tuple(positions), radius
            if not frontier:
                break
        raise ValueError("Hit stereo anchor is not uniquely representable within bounded local context")

    def check(self, product):
        candidate = molecule(product)
        if candidate.GetNumAtoms() > 128:
            raise ValueError("Hit stereo anchor molecule exceeds bounded comparison size")
        used_tetra, used_bonds = set(), set()
        for query, positions, _radius in self.queries:
            achiral_matches = candidate.GetSubstructMatches(
                query, useChirality=False, uniquify=False, maxMatches=65
            )
            if len(achiral_matches) >= 65:
                raise ValueError("Hit stereo anchor match bound reached")
            locations = {tuple(match[position] for position in positions) for match in achiral_matches}
            if len(locations) != 1:
                raise ValueError("Hit stereo anchor is missing or has ambiguous correspondence")
            chiral_matches = candidate.GetSubstructMatches(
                query, useChirality=True, uniquify=False, maxMatches=65
            )
            if len(chiral_matches) >= 65:
                raise ValueError("Hit stereo anchor match bound reached")
            if locations != {tuple(match[position] for position in positions) for match in chiral_matches}:
                raise ValueError("Hit stereo anchor is inverted")
            mapped = next(iter(locations))
            used = used_tetra if len(positions) == 1 else used_bonds
            mapped = mapped if len(positions) == 1 else tuple(sorted(mapped))
            if mapped in used:
                raise ValueError("Two hit stereo anchors map to one product anchor")
            used.add(mapped)


def check_hit_stereo_anchors(hit, product):
    """Conservative standalone control; workflows compile anchors only once."""
    HitStereoAnchors(hit).check(product)


@dataclass(frozen=True)
class RouteFirstLimits:
    version: str = VERSION
    max_rule_inspections: int = 512
    max_rule_scans_per_target: int = 128
    max_rules: int = 16
    max_reverse_outcomes: int = 4
    max_outer_guides: int = 8
    max_catalog_index_rows: int = 250000
    max_catalog_rows_per_query: int = 8192
    max_catalog_row_visits: int = 262144
    max_fingerprint_calls: int = 65536
    building_block_top_k: int = 2
    max_combinations: int = 512
    max_products: int = 1024
    max_routes: int = 64
    beam_width: int = 12
    max_steps: int = 2
    max_inner_rules: int = 4
    fallback_rules: int = 2
    reaction_product_limit: int = 16

    def __post_init__(self):
        if self.version != VERSION:
            raise ValueError(f"route_first.version must be {VERSION}")
        for key, value in asdict(self).items():
            if key != "version":
                positive(value, f"route_first.{key}")
        if self.max_steps not in (1, 2):
            raise ValueError("route_first.max_steps must be 1 or 2")
        if self.max_rules > self.max_rule_inspections or self.fallback_rules > self.max_rules:
            raise ValueError("route_first rule budgets are inconsistent")
        ceilings = {
            "max_rule_inspections": 1024,
            "max_reverse_outcomes": 16,
            "max_catalog_index_rows": 500000,
            "max_catalog_rows_per_query": 32768,
            "max_catalog_row_visits": 524288,
            "max_fingerprint_calls": 262144,
            "max_combinations": 4096,
            "max_products": 8192,
            "max_routes": 256,
            "beam_width": 64,
        }
        for key, ceiling in ceilings.items():
            if getattr(self, key) > ceiling:
                raise ValueError(f"route_first.{key} exceeds the implementation ceiling")


class RouteFirstSearch:
    def __init__(self, hit, catalog, rules, limits, deadline_monotonic=None, accept=None, hit_anchors=None):
        self.hit = canonical(hit)
        self.hit_mol = molecule(self.hit)
        self.hit_fp = FP.GetFingerprint(self.hit_mol)
        self.catalog, self.rules, self.limits = catalog, rules, limits
        self.deadline = deadline_monotonic
        self.accept = accept or (lambda product: None)
        self.report = {
            "version": VERSION,
            "limits": asdict(limits),
            "catalog_sha256": catalog.sha256,
            "rules_sha256": rules.sha256,
            "rule_inspections": 0,
            "rule_template_checks": 0,
            "rule_template_exclusions": [],
            "rules_selected": [],
            "reverse_attempts": 0,
            "reverse_failures": 0,
            "reverse_outcomes": 0,
            "reverse_fallbacks": 0,
            "retrieval_guides_normalized": 0,
            "retrieval_guides_unusable": 0,
            "catalog_index_rows_scanned": 0,
            "catalog_index_hashes": 0,
            "anchor_index_queries": min(2, len(hit_anchors.queries)) if hit_anchors else 0,
            "anchor_index_queries_unindexed": max(0, len(hit_anchors.queries) - 2) if hit_anchors else 0,
            "anchor_index_match_checks": 0,
            "anchor_index_matches": [],
            "anchor_index_postings_retained": [],
            "anchor_index_postings_truncated": [],
            "anchor_index_incomplete": False,
            "guide_anchor_checks": 0,
            "anchor_posting_metadata_visits": 0,
            "anchor_role_checks": 0,
            "anchor_role_compatible": 0,
            "catalog_row_visits": 0,
            "catalog_bucket_rankings": 0,
            "fingerprint_calls": 1,
            "admission_fingerprint_calls": 0,
            "role_queries": 0,
            "combinations": 0,
            "products_generated": 0,
            "route_audit_attempts": 0,
            "routes_verified": 0,
            "product_rejections": {},
            "limits_reached": [],
            "selection_policy": "Product-template hit match first; bounded reverse precursor guides; explicit hit-role fallback",
            "retrieval_scope": "One bounded catalog feature index; each guide ranks charged bucket metadata and scans at most the nearest 8192 catalog rows; exact target ID is additionally considered",
            "evidence_scope": "Reverse results only guide retrieval. Forward rule replay and later route audit/evidence are authoritative.",
        }
        self._nearest_cache = {}
        self._rule_selection_cache = {}
        self._eligible_rules = None
        self._reverse_cache = {}
        self._retrieval_target_cache = {}
        self._fp_cache = {self.hit: self.hit_fp}
        self._catalog_index = None
        self._catalog_pools = {}
        self._anchor_queries = tuple(item[0] for item in hit_anchors.queries[:2]) if hit_anchors else ()
        self._anchor_postings = ()
        self._anchor_scores = {}
        self._guide_anchor_cache = {}
        self._phase = 1
        self.routes = []
        self._route_keys = set()
        self._pending = []

    def _stop(self, name=None):
        if self.deadline is not None and time.monotonic() >= self.deadline:
            name = "deadline_monotonic"
        if name is not None and name not in self.report["limits_reached"]:
            self.report["limits_reached"].append(name)
        return name is not None

    def _rank(self, smiles, target=None):
        if self._stop():
            return 0
        fp = self._fingerprint(smiles)
        reference = self.hit_fp if target is None else self._fingerprint(target)
        return (
            -DataStructs.TanimotoSimilarity(reference, fp) if fp is not None and reference is not None else 0
        )

    def _fingerprint(self, smiles):
        if smiles in self._fp_cache:
            return self._fp_cache[smiles]
        phase_limit = (
            self.limits.max_fingerprint_calls // 2
            if self._phase == 1 and self.limits.max_steps == 2
            else self.limits.max_fingerprint_calls
        )
        if self.report["fingerprint_calls"] + getattr(self.accept, "fingerprint_calls", 0) >= phase_limit:
            self._stop(
                "depth_one_fingerprint_budget"
                if self._phase == 1 and self.limits.max_steps == 2
                else "max_fingerprint_calls"
            )
            return None
        self.report["fingerprint_calls"] += 1
        fp = FP.GetFingerprint(molecule(smiles))
        self._fp_cache[smiles] = fp
        return fp

    @staticmethod
    def _features(mol):
        heavy = mol.GetNumHeavyAtoms()
        hetero = sum(atom.GetAtomicNum() not in (1, 6) for atom in mol.GetAtoms())
        aromatic = sum(atom.GetIsAromatic() for atom in mol.GetAtoms())
        rings = mol.GetRingInfo().NumRings()
        return (heavy // 4, hetero // 2, aromatic // 4, min(rings, 8))

    def _index(self):
        if self._catalog_index is not None:
            return self._catalog_index
        buckets = {}
        reserve = min(256, max(1, self.limits.max_catalog_rows_per_query // 4))
        postings = [[] for _ in self._anchor_queries]
        match_totals = [0 for _ in self._anchor_queries]
        for smiles, mol in itertools.islice(self.catalog.mols.items(), self.limits.max_catalog_index_rows):
            if self._stop():
                break
            self.report["catalog_index_rows_scanned"] += 1
            self.report["catalog_index_hashes"] += 1
            key = self._features(mol)
            score = int.from_bytes(hashlib.sha256(smiles.encode()).digest()[:8], "big")
            buckets.setdefault(key, []).append((score, smiles))
            for anchor_index, anchor in enumerate(self._anchor_queries):
                if self.report["anchor_index_match_checks"] >= 500000:
                    self.report["anchor_index_incomplete"] = True
                    break
                self.report["anchor_index_match_checks"] += 1
                if mol.HasSubstructMatch(anchor, useChirality=True):
                    match_totals[anchor_index] += 1
                    heap = postings[anchor_index]
                    item = (-score, smiles)
                    if len(heap) < reserve:
                        heapq.heappush(heap, item)
                    elif item > heap[0]:
                        heapq.heapreplace(heap, item)
        if len(self.catalog.mols) > self.limits.max_catalog_index_rows:
            self._stop("max_catalog_index_rows")
        for rows in buckets.values():
            rows.sort()
        self._anchor_postings = tuple(
            tuple(smiles for _score, smiles in sorted(heap, reverse=True)) for heap in postings
        )
        self._anchor_scores = {
            smiles: -negative_score for heap in postings for negative_score, smiles in heap
        }
        self.report["anchor_index_matches"] = match_totals
        self.report["anchor_index_postings_retained"] = [len(posting) for posting in self._anchor_postings]
        self.report["anchor_index_postings_truncated"] = [
            total > len(posting) for total, posting in zip(match_totals, self._anchor_postings)
        ]
        if self.report["catalog_index_rows_scanned"] < min(
            len(self.catalog.mols), self.limits.max_catalog_index_rows
        ):
            self.report["anchor_index_incomplete"] = True
        self._catalog_index = buckets
        self.report["catalog_index_size"] = sum(map(len, buckets.values()))
        self.report["catalog_index_buckets"] = len(buckets)
        return buckets

    def _pool(self, target):
        if target in self._catalog_pools:
            return self._catalog_pools[target]
        buckets = self._index()
        if self._stop():
            return ()
        phase_limit = (
            self.limits.max_catalog_row_visits // 2
            if self._phase == 1 and self.limits.max_steps == 2
            else self.limits.max_catalog_row_visits
        )
        if self.report["catalog_row_visits"] + len(buckets) > phase_limit:
            self._stop("max_catalog_row_visits")
            return ()
        self.report["catalog_row_visits"] += len(buckets)
        self.report["catalog_bucket_rankings"] += len(buckets)
        features = self._features(molecule(target))

        def distance(key):
            return (
                4 * abs(key[0] - features[0])
                + 3 * abs(key[1] - features[1])
                + 2 * abs(key[2] - features[2])
                + 2 * abs(key[3] - features[3]),
                key,
            )

        selected = []
        for key in sorted(buckets, key=distance):
            remaining = self.limits.max_catalog_rows_per_query - len(selected)
            if remaining <= 0:
                break
            selected.extend(smiles for _, smiles in buckets[key][:remaining])
        pool = tuple(selected)
        self._catalog_pools[target] = pool
        return pool

    def _eligible_rule_ids(self):
        if self._eligible_rules is not None:
            return self._eligible_rules
        eligible = []
        bounded_ids = sorted(self.rules.rules)[: self.limits.max_rule_scans_per_target]
        for ident in bounded_ids:
            if self._stop():
                break
            if self.report["rule_inspections"] >= self.limits.max_rule_inspections:
                self._stop("max_rule_inspections")
                break
            self.report["rule_inspections"] += 1
            self.report["rule_template_checks"] += 1
            rule = self.rules.get(ident)
            _, issues = template_maps(rule.reaction)
            if any(issues.values()):
                self.report["rule_template_exclusions"].append(
                    {"rule_id": ident, "template_audit": issues}
                )
            else:
                eligible.append(ident)
        # A partial precheck is usable, but never admits uninspected rules.
        self._eligible_rules = tuple(eligible)
        return self._eligible_rules

    def _select_rules(self, target, maximum, fallback=0):
        key = (target, maximum, fallback)
        if key in self._rule_selection_cache:
            return self._rule_selection_cache[key]
        mol = molecule(target)
        primary, secondary = [], []
        eligible_ids = self._eligible_rule_ids()
        # For inner searches, this is the same ordering used by the original
        # full scan's primary sort, so the first maximum matches are final.
        if not fallback:
            eligible_ids = tuple(sorted(
                eligible_ids,
                key=lambda ident: (
                    -self.rules.get(ident).reaction.GetProducts()[0].GetNumHeavyAtoms(), ident
                ),
            ))
        complete = True
        for ident in eligible_ids:
            if self._stop():
                complete = False
                break
            if self.report["rule_inspections"] >= self.limits.max_rule_inspections:
                self._stop("max_rule_inspections")
                complete = False
                break
            self.report["rule_inspections"] += 1
            rule = self.rules.get(ident)
            product_template = rule.reaction.GetProducts()[0]
            if mol.HasSubstructMatch(product_template, useChirality=True):
                primary.append(rule)
                if not fallback and len(primary) >= maximum:
                    break
            elif fallback and any(
                mol.HasSubstructMatch(q, useChirality=True) for q in rule.reaction.GetReactants()
            ):
                secondary.append(rule)
        primary.sort(key=lambda rule: (-rule.reaction.GetProducts()[0].GetNumHeavyAtoms(), rule.id))
        secondary.sort(key=lambda rule: rule.id)
        selected = [(rule, False) for rule in primary[:maximum]] + [
            (rule, True) for rule in secondary[: max(0, min(fallback, maximum - len(primary)))]
        ]
        if complete and not self._stop():
            self._rule_selection_cache[key] = selected
        return selected

    def _reverse(self, rule, target, fallback):
        key = (rule.id, target, fallback)
        if key in self._reverse_cache:
            return self._reverse_cache[key]
        guides = []
        if not fallback and not self._stop():
            self.report["reverse_attempts"] += 1
            try:
                lhs, rhs = rule.smarts.split(">>")
                reverse = rdChemReactions.ReactionFromSmarts(rhs + ">>" + lhs)
                with rdBase.BlockLogs():
                    reverse.Initialize()
                    outcomes = reverse.RunReactants(
                        (molecule(target),), maxProducts=self.limits.max_reverse_outcomes
                    )
                for outcome in outcomes:
                    if len(outcome) != len(rule.roles):
                        continue
                    try:
                        parts = []
                        for mol in outcome:
                            Chem.SanitizeMol(mol)
                            parts.append(canonical(absolute_smiles(mol, isomericSmiles=True)))
                        guides.append(tuple(parts))
                    except (ValueError, RuntimeError):
                        self.report["reverse_failures"] += 1
            except (ValueError, RuntimeError, AttributeError):
                self.report["reverse_failures"] += 1
        if (
            not guides
            and not fallback
            and any(
                molecule(target).HasSubstructMatch(q, useChirality=True) for q in rule.reaction.GetReactants()
            )
        ):
            fallback = True
            self.report["reverse_fallbacks"] += 1
        if fallback:
            # Coarse, explicitly capped fallback. Unmatched roles use a small
            # retrieval query so a pair of nearly whole-hit blocks is not forced.
            mol = molecule(target)
            guides = [
                tuple(
                    target if mol.HasSubstructMatch(q, useChirality=True) else "C"
                    for q in rule.reaction.GetReactants()
                )
            ]
        guides = sorted(set(guides))[: self.limits.max_reverse_outcomes]
        self.report["reverse_outcomes"] += len(guides)
        self._reverse_cache[key] = guides
        return guides

    def _retrieval_target(self, target):
        """Erase reverse-SMARTS dummy/radical artifacts for retrieval only.

        The original rule, catalog block and forward replay are never changed.
        An unconstructive guide falls back to its original literal query.
        """
        if target in self._retrieval_target_cache:
            return self._retrieval_target_cache[target]
        try:
            source = molecule(target)
            edited = Chem.RWMol(source)
            changed = False
            for atom in edited.GetAtoms():
                if atom.GetNumRadicalElectrons():
                    atom.SetNumRadicalElectrons(0)
                    atom.SetNumExplicitHs(0)
                    atom.SetNoImplicit(False)
                    changed = True
            for index in reversed(range(edited.GetNumAtoms())):
                if edited.GetAtomWithIdx(index).GetAtomicNum() == 0:
                    edited.RemoveAtom(index)
                    changed = True
            if edited.GetNumAtoms() == 0:
                raise ValueError("Guide has no concrete atoms")
            clean = edited.GetMol()
            Chem.SanitizeMol(clean)
            if len(Chem.GetMolFrags(clean)) != 1:
                raise ValueError("Guide is disconnected after placeholder removal")
            result = canonical(Chem.MolToSmiles(clean, isomericSmiles=True))
            if changed:
                self.report["retrieval_guides_normalized"] += 1
        except (ValueError, RuntimeError):
            self.report["retrieval_guides_unusable"] += 1
            result = target
        self._retrieval_target_cache[target] = result
        return result

    def _guide_anchors(self, target):
        if target in self._guide_anchor_cache:
            return self._guide_anchor_cache[target]
        matched = []
        if self._stop():
            return ()
        mol = molecule(target)
        for index, query in enumerate(self._anchor_queries):
            if self._stop():
                break
            if self.report["guide_anchor_checks"] >= 2048:
                self._stop("guide_anchor_check_budget")
                break
            self.report["guide_anchor_checks"] += 1
            if mol.HasSubstructMatch(query, useChirality=True):
                matched.append(index)
        result = tuple(matched)
        self._guide_anchor_cache[target] = result
        return result

    def _nearest(self, target, role):
        key = (target, role)
        if key in self._nearest_cache:
            return self._nearest_cache[key]
        self.report["role_queries"] += 1
        query = Chem.MolFromSmarts(role)
        retrieval_target = self._retrieval_target(target)
        # Radical hydrogens in reversed templates can retain useful local
        # positional features in the fingerprint. Dummy atoms cannot, so only
        # those use the normalized fingerprint; both remain retrieval hints.
        fingerprint_target = (
            retrieval_target
            if any(atom.GetAtomicNum() == 0 for atom in molecule(target).GetAtoms())
            else target
        )
        fp = self._fingerprint(fingerprint_target)
        if fp is None:
            return ()
        ranked = []
        pool = self._pool(retrieval_target)
        active_anchors = self._guide_anchors(target)
        anchor_strength = {}
        for index in active_anchors:
            for smiles in self._anchor_postings[index]:
                anchor_strength[smiles] = anchor_strength.get(smiles, 0) + 1
        complete = True
        phase_limit = (
            self.limits.max_catalog_row_visits // 2
            if self._phase == 1 and self.limits.max_steps == 2
            else self.limits.max_catalog_row_visits
        )
        if self.report["catalog_row_visits"] + len(anchor_strength) > phase_limit:
            self._stop("max_catalog_row_visits")
            return ()
        self.report["catalog_row_visits"] += len(anchor_strength)
        self.report["anchor_posting_metadata_visits"] += len(anchor_strength)
        ordered_anchors = sorted(
            anchor_strength,
            key=lambda smiles: (-anchor_strength[smiles], self._anchor_scores[smiles], smiles),
        )
        chosen_anchors = []
        posting_attempts = 0
        posting_limit = min(256, self.limits.max_catalog_rows_per_query)
        for smiles in ordered_anchors[:posting_limit]:
            if self._stop():
                complete = False
                break
            if self.report["catalog_row_visits"] >= phase_limit:
                self._stop("max_catalog_row_visits")
                complete = False
                break
            self.report["catalog_row_visits"] += 1
            self.report["anchor_role_checks"] += 1
            posting_attempts += 1
            mol = self.catalog.mols[smiles]
            if mol.HasSubstructMatch(query, useChirality=True):
                self.report["anchor_role_compatible"] += 1
                chosen_anchors.append(smiles)
        remaining = self.limits.max_catalog_rows_per_query - posting_attempts
        base_pool = tuple(smiles for smiles in pool if smiles not in chosen_anchors)[:remaining]
        if target in self.catalog.mols and target not in chosen_anchors and target not in base_pool:
            base_pool = (target,) + base_pool[: max(0, remaining - 1)] if remaining else ()
        for smiles in itertools.chain(chosen_anchors, base_pool):
            if self._stop():
                complete = False
                break
            if smiles not in chosen_anchors:
                if self.report["catalog_row_visits"] >= phase_limit:
                    self._stop("max_catalog_row_visits")
                    complete = False
                    break
                self.report["catalog_row_visits"] += 1
                mol = self.catalog.mols[smiles]
                if not mol.HasSubstructMatch(query, useChirality=True):
                    continue
            candidate_fp = self._fingerprint(smiles)
            if candidate_fp is None:
                complete = False
                break
            score = DataStructs.TanimotoSimilarity(fp, candidate_fp)
            ranked.append((-score, -anchor_strength.get(smiles, 0), smiles))
        anchored = sorted((row for row in ranked if row[1] < 0), key=lambda row: (row[1], row[0], row[2]))
        ordinary = sorted((row for row in ranked if row[1] == 0), key=lambda row: (row[0], row[2]))
        selected = []
        if anchored:
            selected.append(anchored.pop(0))
        if ordinary and len(selected) < self.limits.building_block_top_k:
            selected.append(ordinary.pop(0))
        selected.extend((anchored + ordinary)[: self.limits.building_block_top_k - len(selected)])
        blocks = tuple(
            self.catalog.blocks[sorted(self.catalog.by_smiles[smiles])[0]]
            for _score, _anchors, smiles in selected
        )
        if complete:
            self._nearest_cache[key] = blocks
        return blocks

    def _stage(self, rule, guide, prefix, replace=None, budget=None, product_budget=None):
        options = []
        for slot, (target, role) in enumerate(zip(guide, rule.roles)):
            if replace is not None and replace[0] == slot:
                options.append(replace[1])
            else:
                options.append(
                    tuple(
                        {
                            "root": f"{prefix}b{slot}",
                            "product": block.smiles,
                            "nodes": [
                                {
                                    "id": f"{prefix}b{slot}",
                                    "building_block_id": block.id,
                                }
                            ],
                        }
                        for block in self._nearest(target, role)
                    )
                )
        states = []
        for combination in itertools.product(*options):
            if self._stop():
                break
            if self.report["combinations"] >= self.limits.max_combinations or (
                budget is not None and self.report["combinations"] >= budget
            ):
                self._stop("max_combinations" if budget is None else "depth_one_combination_budget")
                break
            self.report["combinations"] += 1
            reactants = tuple(item["product"] for item in combination)
            products, capped = self.rules.apply(rule.id, reactants, self.limits.reaction_product_limit)
            if capped:
                self._stop("reaction_product_limit")
            for product in products:
                if self.report["products_generated"] >= (
                    product_budget if product_budget is not None else self.limits.max_products
                ):
                    self._stop("depth_one_product_budget" if product_budget is not None else "max_products")
                    break
                self.report["products_generated"] += 1
                nodes = [node for item in combination for node in item["nodes"]]
                nodes.append(
                    {
                        "id": f"{prefix}r",
                        "rule_id": rule.id,
                        "reactants": [item["root"] for item in combination],
                        "product": product,
                    }
                )
                states.append({"root": f"{prefix}r", "product": product, "nodes": nodes})
            if self.report["products_generated"] >= (
                product_budget if product_budget is not None else self.limits.max_products
            ):
                break
        return states

    def _retain(self, states):
        # Keep the bounded forward outcomes until both search phases finish.
        # Product admission and route audit then spend their budget on the
        # closest eligible exact products, rather than the first rule's output.
        for state in states:
            if self._stop():
                return
            self._pending.append(state)

    def _admit(self):
        if self._stop():
            return
        unique_products = sorted({state["product"] for state in self._pending})
        product_ranks = {product: self._rank(product) for product in unique_products}
        if self._stop() or "max_fingerprint_calls" in self.report["limits_reached"]:
            return
        eligible = set()
        for product in sorted(unique_products, key=lambda product: (product_ranks[product], product)):
            if self._stop():
                return
            if (
                self.report["fingerprint_calls"] + getattr(self.accept, "fingerprint_calls", 0)
                >= self.limits.max_fingerprint_calls
            ):
                self._stop("max_fingerprint_calls")
                break
            try:
                self.accept(product)
            except (ValueError, KeyError, TypeError) as exc:
                reason = f"{type(exc).__name__}: {exc}"
                self.report["product_rejections"][reason] = (
                    self.report["product_rejections"].get(reason, 0) + 1
                )
                continue
            finally:
                self.report["admission_fingerprint_calls"] = getattr(self.accept, "fingerprint_calls", 0)
            eligible.add(product)
        if self._stop():
            return
        states = sorted(
            (state for state in self._pending if state["product"] in eligible),
            key=lambda state: (
                product_ranks[state["product"]],
                state["product"],
                sum("rule_id" in node for node in state["nodes"]),
                str(state["nodes"]),
            ),
        )
        first, alternatives, seen_products = [], [], set()
        for state in states:
            if state["product"] in seen_products:
                alternatives.append(state)
            else:
                first.append(state)
                seen_products.add(state["product"])
        for state in itertools.chain(first, alternatives):
            if self._stop():
                return
            if self.report["route_audit_attempts"] >= self.limits.max_routes:
                self._stop("max_routes")
                return
            route = {"schema_version": 1, "root": state["root"], "nodes": state["nodes"]}
            key = str(route)
            if key in self._route_keys:
                continue
            self.report["route_audit_attempts"] += 1
            try:
                verification = verify_route(
                    route, self.catalog, self.rules, self.limits.reaction_product_limit
                )
                audit = audit_route(route, self.catalog, self.rules, self.limits.reaction_product_limit)
                if not audit["strict_eligible"]:
                    raise ValueError("Full route atom/stereo audit is not strict eligible")
            except (ValueError, KeyError, TypeError) as exc:
                reason = f"{type(exc).__name__}: {exc}"
                self.report["product_rejections"][reason] = (
                    self.report["product_rejections"].get(reason, 0) + 1
                )
                continue
            self._route_keys.add(key)
            self.routes.append(
                {
                    "route": route,
                    "verification": verification,
                    "generation": {
                        "policy": VERSION,
                        "reference": "self_product_design",
                        "route_audit_binding_sha256": audit["route_binding_sha256"],
                    },
                }
            )
            self.report["routes_verified"] += 1

    def run(self):
        outer = self._select_rules(self.hit, self.limits.max_rules, fallback=self.limits.fallback_rules)
        self.report["rules_selected"] = [
            {"rule_id": rule.id, "fallback": fallback} for rule, fallback in outer
        ]
        depth_one_budget = (
            self.limits.max_combinations
            if self.limits.max_steps == 1
            else max(1, self.limits.max_combinations // 2)
        )
        plan_candidates = []
        for rule, fallback in outer:
            if self._stop():
                break
            for guide in self._reverse(rule, self.hit, fallback):
                if self._stop():
                    break
                normalized = tuple(self._retrieval_target(part) for part in guide)
                artifacts = sum(a != b for a, b in zip(guide, normalized))
                size_delta = abs(
                    sum(molecule(part).GetNumHeavyAtoms() for part in normalized)
                    - self.hit_mol.GetNumHeavyAtoms()
                )
                plan_candidates.append((rule, guide, fallback, artifacts, size_delta))
        plan_candidates.sort(
            key=lambda item: (
                item[2],
                -item[0].reaction.GetProducts()[0].GetNumHeavyAtoms(),
                item[3],
                item[4],
                item[0].id,
                item[1],
            )
        )
        plans = []
        seen_rules = set()
        for only_new_rules in (True, False):
            for rule, guide, _fallback, _artifacts, _size_delta in plan_candidates:
                if len(plans) >= self.limits.max_outer_guides:
                    break
                if only_new_rules and rule.id in seen_rules:
                    continue
                if (rule, guide) in plans:
                    continue
                plans.append((rule, guide))
                seen_rules.add(rule.id)
        if len(plan_candidates) > len(plans):
            self._stop("max_outer_guides")
        self.report["outer_plans_selected"] = [
            {"rule_id": rule.id, "guide": list(guide)} for rule, guide in plans
        ]
        depth_one_product_budget = (
            self.limits.max_products if self.limits.max_steps == 1 else max(1, self.limits.max_products // 2)
        )
        for rule, guide in plans:
            if self._stop() or self.report["combinations"] >= depth_one_budget:
                break
            states = self._stage(
                rule,
                guide,
                "o",
                budget=depth_one_budget,
                product_budget=depth_one_product_budget,
            )
            self._retain(states)
        if self.limits.max_steps == 2 and not self._stop():
            self._phase = 2
            for outer_rule, guide in plans:
                for slot, target in enumerate(guide):
                    if self._stop() or self.report["combinations"] >= self.limits.max_combinations:
                        break
                    inner_rules = self._select_rules(target, self.limits.max_inner_rules)
                    for inner_rule, _ in inner_rules:
                        for inner_guide in self._reverse(inner_rule, target, False):
                            inner = self._stage(inner_rule, inner_guide, "i")
                            inner = sorted(
                                inner,
                                key=lambda state: (self._rank(state["product"], target), state["product"]),
                            )[: self.limits.beam_width]
                            if not inner:
                                continue
                            states = self._stage(outer_rule, guide, "o", replace=(slot, inner))
                            self._retain(states)
                            if self._stop() or self.report["combinations"] >= self.limits.max_combinations:
                                break

        self._admit()
        self.report["total_fingerprint_calls"] = (
            self.report["fingerprint_calls"] + self.report["admission_fingerprint_calls"]
        )
        self.report["distinct_exact_products"] = len(
            {item["verification"]["product"] for item in self.routes}
        )
        self.report["depth_counts"] = {
            str(depth): sum(len(item["verification"]["steps"]) == depth for item in self.routes)
            for depth in (1, 2)
        }
        return self.routes, self.report


def generate(hit, catalog, rules, config, *, deadline_monotonic=None, accept=None, hit_anchors=None):
    if not isinstance(config, dict):
        raise ValueError("route_first must be an object")
    if set(config) - set(RouteFirstLimits.__dataclass_fields__):
        raise ValueError("Unknown route_first option")
    limits = RouteFirstLimits(**config)
    if deadline_monotonic is not None and (
        isinstance(deadline_monotonic, bool)
        or not isinstance(deadline_monotonic, (int, float))
        or not math.isfinite(deadline_monotonic)
    ):
        raise ValueError("Invalid route-first deadline")
    return RouteFirstSearch(hit, catalog, rules, limits, deadline_monotonic, accept, hit_anchors).run()
