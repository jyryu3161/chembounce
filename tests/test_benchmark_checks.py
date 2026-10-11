"""Negative controls for the retained-fragment checks used in the benchmark."""

from rdkit import Chem, DataStructs
import chemistry
from benchmarks.benchmark_attachment import assess


def test_retained_fragment_check_rejects_inverted_stereo():
    case = {"input": "C[C@H](O)c1ccccc1", "fragment": "c1ccccc1", "replacement": "c1ccncc1"}
    target = Chem.MolFromSmiles(case["input"])
    for product, should_fail in [("C[C@H](O)c1ccncc1", False), ("C[C@@H](O)c1ccncc1", True)]:
        score = DataStructs.TanimotoSimilarity(
            chemistry.count_fingerprint(target), chemistry.count_fingerprint(Chem.MolFromSmiles(product))
        )
        result = assess(case, {"candidates": [(score, product)]}, 10)
        assert bool(result["checks_failed"].get("retained_substructure_or_stereo")) == should_fail


def test_retained_fragment_check_rejects_lost_isotope():
    case = {"input": "[13CH3]c1ccccc1", "fragment": "c1ccccc1", "replacement": "c1ccncc1"}
    result = assess(case, {"candidates": [(0.0, "Cc1ccncc1")]}, 10)
    assert result["checks_failed"]["heavy_atom_or_isotope_loss"] == 1
    assert result["checks_failed"]["retained_substructure_or_stereo"] == 1
