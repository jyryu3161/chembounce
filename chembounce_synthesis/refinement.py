"""Bounded multi-seed reevaluation with an auditable representative pose."""

import math
import statistics
import time
from dataclasses import asdict, dataclass, replace

from .common import positive
from .conformers import score_candidate

POSE_METRICS = tuple(f"{ref}_{metric}" for ref in ("hit", "design") for metric in ("shape", "pharmacophore"))


@dataclass(frozen=True)
class RefinementSettings:
    # Opt in explicitly so historical configurations retain their scoring protocol.
    top_n: int = 0
    num_conformers: int = 24
    seeds: tuple = (42, 43, 44)
    threshold_policy: str = "all_seeds"
    max_score_range: float | None = None
    max_metric_range: float | None = None
    adaptive_max_conformers: int | None = None
    conformer_growth_factor: int = 2
    borderline_margin: float = 0.0
    max_total_conformers: int | None = None
    ranking_statistic: str = "median"

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
        if not isinstance(self.threshold_policy, str) or self.threshold_policy not in {
            "all_seeds",
            "representative",
        }:
            raise ValueError("refinement.threshold_policy must be all_seeds or representative")
        if not isinstance(self.ranking_statistic, str) or self.ranking_statistic not in {"median", "minimum"}:
            raise ValueError("refinement.ranking_statistic must be median or minimum")
        for key in ("max_score_range", "max_metric_range", "borderline_margin"):
            value = getattr(self, key)
            if value is None and key in {"max_score_range", "max_metric_range"}:
                continue
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0
            ):
                raise ValueError(f"refinement.{key} must be a finite nonnegative number")
        if (self.max_score_range is not None or self.max_metric_range is not None) and len(self.seeds) < 2:
            raise ValueError("Refinement range tolerance requires at least two seeds")
        if self.adaptive_max_conformers is not None:
            positive(self.adaptive_max_conformers, "refinement.adaptive_max_conformers")
            if self.adaptive_max_conformers < self.num_conformers:
                raise ValueError("refinement.adaptive_max_conformers must be >= num_conformers")
        positive(self.conformer_growth_factor, "refinement.conformer_growth_factor")
        if self.conformer_growth_factor < 2:
            raise ValueError("refinement.conformer_growth_factor must be >= 2")
        if self.max_total_conformers is not None:
            positive(self.max_total_conformers, "refinement.max_total_conformers")
        if (
            self.adaptive_max_conformers or self.num_conformers
        ) > self.num_conformers and self.max_total_conformers is None:
            raise ValueError("Adaptive refinement requires refinement.max_total_conformers")


@dataclass
class RefinementBudget:
    """Logical candidate-conformer requests; reference generation is not included.

    Shared across products and design alternatives. Reserve an entire seed block
    before starting it, avoiding partially sampled or selectively omitted seeds.
    """

    limit: int | None = None
    requested: int = 0
    score_calls: int = 0
    reserved_for_others: int = 0

    def reserve(self, num_conformers, seed_count):
        cost = num_conformers * seed_count
        if self.limit is not None and self.requested + cost + self.reserved_for_others > self.limit:
            return False
        self.requested += cost
        return True

    def snapshot(self):
        return {
            "limit": self.limit,
            "candidate_conformer_requests": self.requested,
            "score_calls": self.score_calls,
            "reserved_for_others": self.reserved_for_others,
        }


def _trial_metrics(info):
    return {
        f"{ref}_{metric}": info[ref][metric]
        for ref in ("hit", "design")
        for metric in ("shape", "pharmacophore")
        if ref in info and metric in info[ref]
    }


def _summarize(trials, results, refinement, minimum_scores):
    report = {"trials": trials, "seed_count": len(trials), "seed_successes": len(results)}
    if len(results) != len(trials):
        return None, {**report, "status": "failed", "convergence_status": "failed_seeds"}
    values = [result[0]["score"] for _, result in results]
    median = statistics.median(values)
    seed, representative = sorted(results, key=lambda r: (r[1][0]["score"], r[0]))[(len(results) - 1) // 2]
    metrics = [_trial_metrics(result[0]) for _, result in results]
    metric_ranges = {
        key: {"min": min(m[key] for m in metrics), "max": max(m[key] for m in metrics)}
        for key in POSE_METRICS
        if all(key in m for m in metrics)
    }
    if all("candidate_relative_energy" in result[0] for _, result in results):
        energies = [result[0]["candidate_relative_energy"] for _, result in results]
        metric_ranges["relative_energy"] = {"min": min(energies), "max": max(energies)}
    score_range = max(values) - min(values)
    metric_unstable = refinement.max_metric_range is not None and any(
        values["max"] - values["min"] > refinement.max_metric_range
        for key, values in metric_ranges.items()
        if key in POSE_METRICS
    )
    if refinement.max_score_range is None and refinement.max_metric_range is None:
        convergence = "not_assessed" if len(values) > 1 else "insufficient_seeds"
    else:
        score_unstable = refinement.max_score_range is not None and score_range > refinement.max_score_range
        convergence = "unstable" if score_unstable or metric_unstable else "converged"
    assessed = {key: value for key, value in minimum_scores.items() if key in POSE_METRICS}
    # All requested criteria must be present. Missing metrics cannot pass silently.
    rep_metrics = _trial_metrics(representative[0])
    failed_thresholds = (
        [
            key
            for key, value in assessed.items()
            if (key not in metric_ranges or metric_ranges[key]["min"] < value)
        ]
        if refinement.threshold_policy == "all_seeds"
        else [key for key, value in assessed.items() if key not in rep_metrics or rep_metrics[key] < value]
    )
    borderline = any(
        key in metric_ranges
        and metric_ranges[key]["min"] - refinement.borderline_margin
        <= value
        <= metric_ranges[key]["max"] + refinement.borderline_margin
        for key, value in assessed.items()
    )
    return representative, {
        **report,
        "status": "success",
        "pose_seed": seed,
        "seed_score_median": median,
        "seed_score_min": min(values),
        "seed_score_max": max(values),
        "seed_score_range": score_range,
        "seed_score_std": statistics.pstdev(values),
        "ranking_score": min(values) if refinement.ranking_statistic == "minimum" else median,
        "metric_ranges": metric_ranges,
        "convergence_status": convergence,
        "convergence_scope": "Observed dispersion across requested seeds at the final ensemble size; not asymptotic convergence or a confidence interval",
        "threshold_policy": refinement.threshold_policy,
        "failed_thresholds": failed_thresholds,
        "borderline": borderline,
    }


def refine_candidate(
    smiles,
    hit,
    design,
    settings,
    refinement,
    hit_pose=None,
    minimum_scores=None,
    budget=None,
    pose_constraints=None,
    deadline_monotonic=None,
    result_validator=None,
):
    """Refine reproducibly; all exported pose components belong to one real pose.

    Each adaptive stage reruns the same seeds with a larger ensemble. Only the
    final complete stage supplies acceptance/ranking; earlier trials remain in
    ``stages``. Observed seed dispersion is not a calibrated confidence interval.
    """
    if deadline_monotonic is not None and (
        isinstance(deadline_monotonic, bool)
        or not isinstance(deadline_monotonic, (int, float))
        or not math.isfinite(deadline_monotonic)
    ):
        raise ValueError("deadline_monotonic must be a finite monotonic-clock timestamp")
    minimum_scores = minimum_scores or {}
    budget = budget if budget is not None else RefinementBudget(refinement.max_total_conformers)
    before = budget.requested
    stages = []
    num_conformers = refinement.num_conformers
    maximum = refinement.adaptive_max_conformers or num_conformers
    result = None
    while True:
        if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
            report = (
                {
                    "status": "deadline_exhausted",
                    "trials": [],
                    "seed_count": len(refinement.seeds),
                    "seed_successes": 0,
                    "convergence_status": "not_evaluated",
                }
                if not stages
                else {**stages[-1], "status": "deadline_exhausted"}
            )
            report["reason"] = "Deadline reached before the next complete refinement stage"
            result = None
            break
        if not budget.reserve(num_conformers, len(refinement.seeds)):
            report = (
                {
                    "status": "budget_exhausted",
                    "trials": [],
                    "seed_count": 0,
                    "seed_successes": 0,
                    "convergence_status": "not_evaluated",
                }
                if not stages
                else {**stages[-1], "status": "budget_exhausted"}
            )
            report = {
                **report,
                "reason": "Cannot fund the next complete seed block without exceeding the global budget or reserved initial evaluations",
                "unfunded_candidate_conformer_requests": num_conformers * len(refinement.seeds),
            }
            result = None
            break
        trials, results = [], []
        deadline_interrupted = False
        for seed_index, seed in enumerate(refinement.seeds):
            if deadline_monotonic is not None and time.monotonic() >= deadline_monotonic:
                remaining = refinement.seeds[seed_index:]
                budget.requested -= num_conformers * len(remaining)
                trials.extend(
                    {
                        "seed": pending_seed,
                        "status": "not_evaluated",
                        "reason": "Deadline reached before this seed",
                        "settings": asdict(
                            replace(settings, num_conformers=num_conformers, seed=pending_seed)
                        ),
                    }
                    for pending_seed in remaining
                )
                deadline_interrupted = True
                break
            seeded = replace(settings, num_conformers=num_conformers, seed=seed)
            budget.score_calls += 1
            try:
                scoring_options = (
                    {"pose_constraints": pose_constraints} if pose_constraints is not None else {}
                )
                candidate = score_candidate(smiles, hit, design, seeded, hit_pose, **scoring_options)
                if result_validator is not None:
                    result_validator(candidate, smiles, seeded)
                if not math.isfinite(candidate[0]["score"]) or not all(
                    math.isfinite(v) for v in _trial_metrics(candidate[0]).values()
                ):
                    raise ValueError("Nonfinite refinement score")
                trials.append({"seed": seed, "status": "success", **candidate[0]})
                results.append((seed, candidate))
            except (ValueError, RuntimeError) as exc:
                trials.append(
                    {"seed": seed, "status": "failed", "reason": str(exc), "settings": asdict(seeded)}
                )
        result, report = _summarize(trials, results, refinement, minimum_scores)
        report["num_conformers"] = num_conformers
        if deadline_interrupted:
            report.update(
                status="deadline_exhausted",
                convergence_status="incomplete_seed_block",
                reason="Deadline reached with unevaluated requested seeds; partial block cannot be accepted",
            )
            stages.append(report)
            result = None
            break
        if stages and "seed_score_median" in stages[-1] and "seed_score_median" in report:
            report["stage_score_median_change"] = (
                report["seed_score_median"] - stages[-1]["seed_score_median"]
            )
            report["stage_metric_min_changes"] = {
                key: values["min"] - stages[-1]["metric_ranges"][key]["min"]
                for key, values in report["metric_ranges"].items()
                if key in stages[-1]["metric_ranges"]
            }
        stages.append(report)
        needs_more = (
            result is None
            or report.get("convergence_status") == "unstable"
            or report.get("borderline", False)
        )
        if not needs_more or num_conformers >= maximum:
            if result is not None and report["convergence_status"] == "unstable":
                result, report = None, {**report, "status": "unstable"}
            elif result is not None and report["failed_thresholds"]:
                result, report = None, {**report, "status": "threshold_rejected"}
            break
        num_conformers = min(maximum, num_conformers * refinement.conformer_growth_factor)
    return result, {
        **report,
        "stages": stages,
        "adaptive_stages": len(stages),
        "candidate_conformer_requests": budget.requested - before,
        "global_budget": budget.snapshot(),
    }
