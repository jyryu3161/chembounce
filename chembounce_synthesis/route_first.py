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
from .route_audit import audit_route
from .routes import verify_route

VERSION = "catalog-rule-route-first-v1"


def check_hit_stereo_anchors(hit, product):
    """Require an unambiguous, stereo-matched radius-two motif at each hit center.

    This deliberately rejects changed/ambiguous local neighborhoods. The atom
    correspondence comes from the rooted local substructure, not a global MCS
    or a potentially relabelled CIP R/S string.
    """
    source, candidate = molecule(hit), molecule(product)
    if source.GetNumAtoms() > 128 or candidate.GetNumAtoms() > 128:
        raise ValueError("Hit stereo anchor molecule exceeds bounded comparison size")
    centers = [
        atom.GetIdx() for atom in source.GetAtoms() if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    ]
    stereo_bonds = [bond for bond in source.GetBonds() if bond.GetStereo() != Chem.BondStereo.STEREONONE]
    if len(centers) + len(stereo_bonds) > 16:
        raise ValueError("Too many specified hit stereo anchors")

    def verify_anchor(roots, radius):
        included, frontier = set(roots), set(roots)
        for _ in range(radius):
            frontier = {
                neighbor.GetIdx()
                for index in frontier
                for neighbor in source.GetAtomWithIdx(index).GetNeighbors()
            } - included
            included.update(frontier)
        rings = source.GetRingInfo().AtomRings()
        for _ in range(len(rings) + 1):
            previous = len(included)
            for ring in rings:
                if included.intersection(ring):
                    included.update(ring)
            if len(included) == previous:
                break
        if len(included) > 64:
            raise ValueError("Hit stereo anchor exceeds bounded local ring context")
        marked = Chem.Mol(source)
        for position, index in enumerate(roots):
            marked.GetAtomWithIdx(index).SetAtomMapNum(1000000 + position)
        fragment = Chem.MolFragmentToSmiles(
            marked, atomsToUse=sorted(included), canonical=True, isomericSmiles=True
        )
        query = Chem.MolFromSmiles(fragment)
        if query is None:
            raise ValueError("Hit stereo anchor could not be represented")
        positions = [
            next((a.GetIdx() for a in query.GetAtoms() if a.GetAtomMapNum() == 1000000 + i), None)
            for i in range(len(roots))
        ]
        if None in positions:
            raise ValueError("Hit stereo anchor correspondence is ambiguous")
        if (
            len(roots) == 1
            and query.GetAtomWithIdx(positions[0]).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
        ):
            raise ValueError("Hit stereo anchor lost defined tetrahedral geometry in local context")
        if len(roots) == 2:
            query_bond = query.GetBondBetweenAtoms(*positions)
            if query_bond is None or query_bond.GetStereo() not in (
                Chem.BondStereo.STEREOE,
                Chem.BondStereo.STEREOZ,
            ):
                raise ValueError("Hit stereo anchor lost defined double-bond geometry in local context")
        for position in positions:
            query.GetAtomWithIdx(position).SetAtomMapNum(0)
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

    for center in centers:
        verify_anchor((center,), 2)
    for bond in stereo_bonds:
        if bond.GetStereo() not in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
            raise ValueError("Unsupported specified hit bond stereo")
        verify_anchor((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()), 1)


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
    def __init__(self, hit, catalog, rules, limits, deadline_monotonic=None, accept=None):
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
            "rules_selected": [],
            "reverse_attempts": 0,
            "reverse_failures": 0,
            "reverse_outcomes": 0,
            "reverse_fallbacks": 0,
            "catalog_index_rows_scanned": 0,
            "catalog_index_hashes": 0,
            "catalog_row_visits": 0,
            "fingerprint_calls": 1,
            "role_queries": 0,
            "combinations": 0,
            "products_generated": 0,
            "routes_verified": 0,
            "product_rejections": {},
            "limits_reached": [],
            "selection_policy": "Product-template hit match first; bounded reverse precursor guides; explicit hit-role fallback",
            "retrieval_scope": "One bounded deterministic stratified hash reservoir over catalog rows, plus exact target ID lookup; role scans/fingerprint work separately charged",
            "evidence_scope": "Reverse results only guide retrieval. Forward rule replay and later route audit/evidence are authoritative.",
        }
        self._nearest_cache = {}
        self._reverse_cache = {}
        self._fp_cache = {self.hit: self.hit_fp}
        self._catalog_pool = None
        self._phase = 1
        self.routes = []
        self._route_keys = set()

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
        if self.report["fingerprint_calls"] >= phase_limit:
            self._stop("max_fingerprint_calls")
            return None
        self.report["fingerprint_calls"] += 1
        fp = FP.GetFingerprint(molecule(smiles))
        self._fp_cache[smiles] = fp
        return fp

    def _pool(self):
        if self._catalog_pool is not None:
            return self._catalog_pool
        size = self.limits.max_catalog_rows_per_query
        quota = max(1, (size + 3) // 4)
        buckets = [[] for _ in range(4)]
        for smiles, mol in itertools.islice(self.catalog.mols.items(), self.limits.max_catalog_index_rows):
            if self._stop():
                break
            self.report["catalog_index_rows_scanned"] += 1
            self.report["catalog_index_hashes"] += 1
            heavy = mol.GetNumHeavyAtoms()
            bucket = 0 if heavy <= 8 else 1 if heavy <= 16 else 2 if heavy <= 28 else 3
            score = int.from_bytes(hashlib.sha256(smiles.encode()).digest()[:8], "big")
            heap = buckets[bucket]
            item = (-score, smiles)
            if len(heap) < quota:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)
        if len(self.catalog.mols) > self.limits.max_catalog_index_rows:
            self._stop("max_catalog_index_rows")
        self._catalog_pool = tuple(sorted(smiles for heap in buckets for _, smiles in heap)[:size])
        self.report["catalog_index_size"] = len(self._catalog_pool)
        return self._catalog_pool

    def _select_rules(self, target, maximum, fallback=0):
        mol = molecule(target)
        primary, secondary = [], []
        for ident in sorted(self.rules.rules)[: self.limits.max_rule_scans_per_target]:
            if self._stop():
                break
            if self.report["rule_inspections"] >= self.limits.max_rule_inspections:
                self._stop("max_rule_inspections")
                break
            self.report["rule_inspections"] += 1
            rule = self.rules.get(ident)
            product_template = rule.reaction.GetProducts()[0]
            if mol.HasSubstructMatch(product_template, useChirality=True):
                primary.append(rule)
            elif fallback and any(
                mol.HasSubstructMatch(q, useChirality=True) for q in rule.reaction.GetReactants()
            ):
                secondary.append(rule)
        primary.sort(key=lambda rule: (-rule.reaction.GetProducts()[0].GetNumHeavyAtoms(), rule.id))
        secondary.sort(key=lambda rule: rule.id)
        return [(rule, False) for rule in primary[:maximum]] + [
            (rule, True) for rule in secondary[: max(0, min(fallback, maximum - len(primary)))]
        ]

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

    def _nearest(self, target, role):
        key = (target, role)
        if key in self._nearest_cache:
            return self._nearest_cache[key]
        self.report["role_queries"] += 1
        query = Chem.MolFromSmarts(role)
        fp = self._fingerprint(target)
        if fp is None:
            return ()
        ranked = []
        pool = self._pool()
        exact = (target,) if target in self.catalog.mols and target not in pool else ()
        complete = True
        for smiles in itertools.chain(exact, pool):
            if self._stop():
                complete = False
                break
            phase_limit = (
                self.limits.max_catalog_row_visits // 2
                if self._phase == 1 and self.limits.max_steps == 2
                else self.limits.max_catalog_row_visits
            )
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
            ranked.append((-score, smiles))
        blocks = tuple(
            self.catalog.blocks[sorted(self.catalog.by_smiles[smiles])[0]]
            for _, smiles in sorted(ranked)[: self.limits.building_block_top_k]
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
        for state in states:
            if self._stop():
                return
            depth = sum("rule_id" in node for node in state["nodes"])
            if len(self.routes) >= self.limits.max_routes:
                self._stop("max_routes")
                return
            if (
                depth == 1
                and self.limits.max_steps == 2
                and sum(len(item["verification"]["steps"]) == 1 for item in self.routes)
                >= max(1, self.limits.max_routes // 2)
            ):
                self._stop("depth_one_route_budget")
                return
            route = {"schema_version": 1, "root": state["root"], "nodes": state["nodes"]}
            key = str(route)
            if key in self._route_keys:
                continue
            try:
                self.accept(state["product"])
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
        plans = []
        depth_one_product_budget = (
            self.limits.max_products if self.limits.max_steps == 1 else max(1, self.limits.max_products // 2)
        )
        for rule, fallback in outer:
            for guide in self._reverse(rule, self.hit, fallback):
                if self._stop():
                    break
                if len(plans) >= self.limits.max_outer_guides:
                    self._stop("max_outer_guides")
                    break
                plans.append((rule, guide))
                states = self._stage(
                    rule,
                    guide,
                    "o",
                    budget=depth_one_budget,
                    product_budget=depth_one_product_budget,
                )
                self._retain(states)
                if (
                    self.report["combinations"] >= depth_one_budget
                    or len(plans) >= self.limits.max_outer_guides
                ):
                    break
            if (
                self.report["combinations"] >= depth_one_budget
                or len(plans) >= self.limits.max_outer_guides
                or self._stop()
            ):
                break
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

        def order(item):
            return (
                self._rank(item["verification"]["product"]),
                item["verification"]["product"],
                str(item["route"]),
            )

        by_depth = {
            depth: sorted(
                (item for item in self.routes if len(item["verification"]["steps"]) == depth),
                key=order,
            )
            for depth in (1, 2)
        }
        first = []
        seen = set()
        for depth in (1, 2):
            for item in by_depth[depth]:
                product = item["verification"]["product"]
                if product not in seen:
                    first.append(item)
                    seen.add(product)
                    break
        remainder = sorted((item for item in self.routes if item not in first), key=order)
        self.routes = (first + remainder)[: self.limits.max_routes]
        self.report["distinct_exact_products"] = len(
            {item["verification"]["product"] for item in self.routes}
        )
        self.report["depth_counts"] = {
            str(depth): sum(len(item["verification"]["steps"]) == depth for item in self.routes)
            for depth in (1, 2)
        }
        return self.routes, self.report


def generate(hit, catalog, rules, config, *, deadline_monotonic=None, accept=None):
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
    return RouteFirstSearch(hit, catalog, rules, limits, deadline_monotonic, accept).run()
