"""Deterministic, bounded similarity exploitation and structural exploration.

This is a search allocation heuristic, not a claim of activity or 3D similarity.
The reserved lane is selected without a mandatory 2D similarity floor.
"""

import hashlib
import heapq
import math

from rdkit import DataStructs


def validate_policy(policy, fraction, seed, candidate_budget=4096):
    if not isinstance(policy, str) or policy not in {"legacy", "diverse"}:
        raise ValueError("search_policy must be legacy or diverse")
    if (
        isinstance(fraction, bool)
        or not isinstance(fraction, (int, float))
        or not math.isfinite(fraction)
        or not 0 <= fraction <= 1
    ):
        raise ValueError("exploration_fraction must be a finite number in [0, 1]")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed <= 2**31 - 1:
        raise ValueError("seed must be an integer in [0, 2**31-1]")
    if isinstance(candidate_budget, bool) or not isinstance(candidate_budget, int) or candidate_budget < 1:
        raise ValueError("search_candidate_budget must be a positive integer")


def stable_key(smiles, seed):
    return int.from_bytes(hashlib.sha256(f"{seed}\0{smiles}".encode()).digest()[:8], "big")


def hash_reservoir(items, limit, *, smiles, seed):
    """Order-independent bounded sample; duplicate identities consume no slots."""
    heap, retained = [], set()
    for item in items:
        key = smiles(item)
        if key in retained:
            continue
        rank = stable_key(key, seed)
        if len(heap) < limit:
            heapq.heappush(heap, (-rank, key, item))
            retained.add(key)
        elif rank < -heap[0][0]:
            _, removed, _ = heapq.heapreplace(heap, (-rank, key, item))
            retained.remove(removed)
            retained.add(key)
    return [item for _, _, item in sorted(heap, key=lambda row: (-row[0], row[1]))]


def select_diverse(
    items,
    limit,
    *,
    score,
    smiles,
    fingerprint,
    fraction=0.25,
    seed=42,
    threshold=None,
    group=None,
    can_select=None,
):
    """Reserve slots for farthest-point exploration from a bounded input pool.

    Low-similarity entries have first access to reserved slots when a threshold
    is supplied. Farthest-point distance is to already selected fingerprints;
    ties are reproducible hashed SMILES. Group novelty precedes distance when
    supplied (e.g. standardized scaffold families). External caps are checked
    for every selection and updated through the returned selection callback.
    """
    if limit <= 0 or not items:
        return [], {"selected": 0, "exploration_selected": 0, "below_threshold_selected": 0}
    unique = {smiles(item): item for item in items}
    ranked = sorted(unique.values(), key=lambda item: (-score(item), smiles(item)))
    reserved = min(limit, math.ceil(limit * fraction))
    priority = [item for item in ranked if threshold is None or score(item) >= threshold]
    selected, selected_ids, fps, groups = [], set(), [], set()
    lanes = {}
    fp_cache = {smiles(item): fingerprint(item) for item in ranked}
    ordered_keys = list(fp_cache)
    all_fps = list(fp_cache.values())
    max_similarities = {key: 0.0 for key in ordered_keys}
    tie_keys = {key: stable_key(key, seed) for key in ordered_keys}
    group_keys = {smiles(item): group(item) for item in ranked} if group is not None else {}

    def add(item, lane):
        key = smiles(item)
        if key in selected_ids or (can_select is not None and not can_select(item, selected)):
            return False
        selected.append(item)
        selected_ids.add(key)
        fps.append(fp_cache[key])
        for other_key, similarity in zip(
            ordered_keys, DataStructs.BulkTanimotoSimilarity(fp_cache[key], all_fps)
        ):
            max_similarities[other_key] = max(max_similarities[other_key], similarity)
        if group is not None:
            groups.add(group_keys[key])
        lanes[key] = lane
        return True

    for item in priority:
        if len(selected) >= limit - reserved:
            break
        add(item, "similarity")
    # Fill reserved slots plus any priority slots with no eligible high-2D item.
    exploration_selected = 0
    available = [item for item in ranked if smiles(item) not in selected_ids]
    while fraction > 0 and available and len(selected) < limit:
        eligible = [item for item in available if can_select is None or can_select(item, selected)]
        if not eligible:
            break
        below = [item for item in eligible if threshold is not None and score(item) < threshold]
        pool = below or eligible

        def novelty(item):
            key = smiles(item)
            distance = 1.0 - max_similarities[key]
            return (group is not None and group_keys[key] not in groups, distance, -tie_keys[key], key)

        chosen = max(pool, key=novelty)
        add(chosen, "exploration")
        exploration_selected += 1
        chosen_id = smiles(chosen)
        available = [item for item in available if smiles(item) != chosen_id]
    # Interleave the reserved lane throughout the output rather than leave all
    # exploration at the end where downstream first-N budgets would erase it.
    high = [item for item in selected if lanes[smiles(item)] == "similarity"]
    explore = [item for item in selected if lanes[smiles(item)] == "exploration"]
    ordered = []
    for index in range(len(selected)):
        expected_exploration = math.ceil((index + 1) * len(explore) / len(selected))
        emitted_exploration = sum(lanes[smiles(item)] == "exploration" for item in ordered)
        lane = explore if explore and emitted_exploration < expected_exploration else high
        if not lane:
            lane = explore or high
        ordered.append(lane.pop(0))
    return ordered, {
        "pool_size": len(ranked),
        "selected": len(ordered),
        "reserved_exploration_slots": reserved,
        "exploration_selected": exploration_selected,
        "below_threshold_selected": sum(
            threshold is not None and score(item) < threshold for item in ordered
        ),
        "lanes": {smiles(item): lanes[smiles(item)] for item in ordered},
    }
