"""Bounded multi-seed reevaluation with an auditable representative pose."""

import math
import statistics
from dataclasses import asdict, dataclass, replace

from .common import positive
from .conformers import score_candidate


@dataclass(frozen=True)
class RefinementSettings:
    # Opt in explicitly so historical configurations retain their scoring protocol.
    top_n: int = 0
    num_conformers: int = 24
    seeds: tuple = (42, 43, 44)

    def __post_init__(self):
        if isinstance(self.top_n, bool) or not isinstance(self.top_n, int) or self.top_n < 0:
            raise ValueError("refinement.top_n must be a nonnegative integer")
        positive(self.num_conformers, "refinement.num_conformers")
        if (
            not isinstance(self.seeds, (list, tuple))
            or not self.seeds
            or any(
                isinstance(s, bool) or not isinstance(s, int) or not 0 <= s <= 2**31 - 1 for s in self.seeds
            )
            or len(set(self.seeds)) != len(self.seeds)
        ):
            raise ValueError("refinement.seeds must contain distinct integers in [0, 2**31-1]")


def refine_candidate(smiles, hit, design, settings, refinement, hit_pose=None):
    """Return median-score representative and every seed's diagnostics.

    Require success on every requested seed. A failed seed cannot silently make
    a molecule look better. For even counts the lower of tied median-distance
    scores wins. Exported components always describe one real pose, not a vector
    assembled from separate seed maxima or medians.
    """
    trials, results = [], []
    for seed in refinement.seeds:
        seeded = replace(settings, num_conformers=refinement.num_conformers, seed=seed)
        try:
            result = score_candidate(smiles, hit, design, seeded, hit_pose)
            if not math.isfinite(result[0]["score"]):
                raise ValueError("Nonfinite refinement score")
            trials.append({"seed": seed, "status": "success", **result[0]})
            results.append((seed, result))
        except (ValueError, RuntimeError) as exc:
            trials.append({"seed": seed, "status": "failed", "reason": str(exc), "settings": asdict(seeded)})
    report = {"trials": trials, "seed_count": len(trials), "seed_successes": len(results)}
    if len(results) != len(trials):
        return None, {**report, "status": "failed"}
    values = [result[0]["score"] for _, result in results]
    median = statistics.median(values)
    seed, representative = sorted(results, key=lambda r: (r[1][0]["score"], r[0]))[(len(results) - 1) // 2]
    return representative, {
        **report,
        "status": "success",
        "pose_seed": seed,
        "seed_score_median": median,
        "seed_score_min": min(values),
        "seed_score_max": max(values),
    }
