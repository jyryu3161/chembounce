import numpy as np
import pytest

from chembounce_synthesis.workflow import _pareto


def row(shape, tanimoto):
    return dict(
        hit_shape=shape,
        hit_pharmacophore=0.6,
        design_shape=0.8,
        design_pharmacophore=0.6,
        relative_energy=0,
        reaction_steps=1,
        hit_tanimoto=tanimoto,
        design_tanimoto=0.5,
    )


def test_pareto_keeps_tradeoffs_and_ties_on_the_same_front():
    rows = [row(0.9, 0.8), row(0.8, 0.9), row(0.8, 0.8), row(0.7, 0.8), row(0.8, 0.8)]
    _pareto(rows)
    assert [r["pareto_front"] for r in rows] == [0, 0, 1, 2, 1]


def test_pareto_fully_ordered_pool_is_independent_of_input_order():
    order = np.random.default_rng(42).permutation(256)
    rows = [row(i / 256, i / 256) for i in order]
    _pareto(rows)
    assert [r["pareto_front"] for r in rows] == list(255 - order)
    _pareto([])


@pytest.mark.parametrize("invalid", [float("nan"), float("inf"), -float("inf")])
def test_pareto_rejects_nonfinite_metrics(invalid):
    with pytest.raises(ValueError, match="finite"):
        _pareto([row(invalid, 0.5)])
