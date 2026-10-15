"""Identical reference objects save alignment work without changing two-pose output."""

from dataclasses import replace
from unittest.mock import patch

import numpy as np
import pytest
from rdkit import Chem

pytest.importorskip("rdkit.Chem.rdShapeAlign", reason="Optional synthesis alignment requires a newer RDKit")

from chembounce_synthesis import conformers as module
from chembounce_synthesis.common import canonical

SETTINGS = module.ConformerSettings(num_conformers=6, max_iterations=1000)


def assert_same_results(actual, expected):
    assert actual[0] == expected[0]
    for observed, baseline in zip(actual[1:], expected[1:]):
        assert canonical(Chem.MolToSmiles(observed)) == canonical(Chem.MolToSmiles(baseline))
        assert observed.GetPropsAsDict() == baseline.GetPropsAsDict()
        assert observed.GetNumConformers() == baseline.GetNumConformers()
        for left, right in zip(observed.GetConformers(), baseline.GetConformers()):
            assert left.GetId() == right.GetId()
            np.testing.assert_array_equal(left.GetPositions(), right.GetPositions())


def evaluate_with_trace(product, reference, *, shared, constraints=None):
    original = module.ensemble
    # Separate wrappers retain exactly the same molecules, energies, diagnostics
    # and conformer order but force the original two-overlay branch. No scoring
    # formula or expected score is copied into this regression.
    factory = original if shared else lambda *args, **kwargs: replace(original(*args, **kwargs))
    with (
        patch.object(module, "ensemble", factory),
        patch.object(module, "_overlay", wraps=module._overlay) as calls,
    ):
        result = module.score_candidate(product, reference, reference, SETTINGS, pose_constraints=constraints)
    return result, calls.call_args_list


@pytest.mark.parametrize(
    "product,reference",
    [
        ("CC(=O)Nc1ccncc1", "CC(=O)Nc1ccccc1"),
        ("CCCCCCO", "CCCCCO"),
        ("C[C@H](O)c1ccncc1", "C[C@H](O)c1ccccc1"),
        ("[2H]C([2H])([2H])c1ccncc1", "Cc1ccccc1"),
    ],
)
def test_reuse_matches_two_real_overlays_for_every_score_diagnostic_and_pose(product, reference):
    baseline, double = evaluate_with_trace(product, reference, shared=False)
    actual, single = evaluate_with_trace(product, reference, shared=True)
    conformer_count = len(module.ensemble(canonical(product), SETTINGS).relative_energies)
    assert len(double) == 2 * conformer_count
    assert len(single) == conformer_count
    assert all(left.args[0] is not right.args[0] for left, right in zip(double[::2], double[1::2]))
    assert_same_results(actual, baseline)
    if product == "CCCCCCO":
        assert conformer_count > 1  # Exercise selection across multiple candidate conformations.


def test_returned_reference_comparisons_and_poses_do_not_share_mutable_output():
    product, reference = "CCCCCCO", "CCCCCO"
    result, _ = evaluate_with_trace(product, reference, shared=True)
    info, hit_pose, design_pose = result
    cid = info["candidate_conformer_id"]
    before = hit_pose.GetConformer(cid).GetPositions().copy()
    shape = info["hit"]["shape"]
    info["design"]["shape"] = -1
    design_pose.SetProp("annotation", "design only")
    design_pose.GetConformer(cid).SetAtomPosition(0, (100, 100, 100))
    assert info["hit"]["shape"] == shape
    assert not hit_pose.HasProp("annotation")
    np.testing.assert_array_equal(hit_pose.GetConformer(cid).GetPositions(), before)
    assert hit_pose is not design_pose and info["hit"] is not info["design"]


def test_identical_smiles_bound_and_generated_references_still_receive_separate_overlays(tmp_path):
    reference, product = "CCCCCO", "CCCCCCO"
    generated = module.ensemble(canonical(reference), SETTINGS)
    molecule = Chem.Mol(generated.mol)
    cid = next(iter(generated.relative_energies))
    for atom in range(molecule.GetNumAtoms()):
        point = molecule.GetConformer(cid).GetAtomPosition(atom)
        molecule.GetConformer(cid).SetAtomPosition(atom, point + Chem.rdGeometry.Point3D(20, 5, -3))
    path = tmp_path / "bound.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(molecule, confId=cid)
    bound = module.reference_from_sdf(path, reference, SETTINGS)
    assert bound is not generated and bound.bound_pose
    with patch.object(module, "_overlay", wraps=module._overlay) as calls:
        info, hit_pose, design_pose = module.score_candidate(
            product, reference, reference, SETTINGS, hit_pose=bound
        )
    count = len(module.ensemble(canonical(product), SETTINGS).relative_energies)
    assert calls.call_count == 2 * count
    assert all(call.args[0] is bound for call in calls.call_args_list[::2])
    assert all(call.args[0] is generated for call in calls.call_args_list[1::2])
    assert info["ensemble_diagnostics"]["hit"] == {}
    assert info["ensemble_diagnostics"]["design"] == generated.failures
    selected = info["candidate_conformer_id"]
    assert not np.allclose(
        hit_pose.GetConformer(selected).GetPositions(), design_pose.GetConformer(selected).GetPositions()
    )


class PoseAssessment:
    """Test-only annotating assessor, including rejection before design alignment."""

    provenance = {"provider": "test fixture"}

    def __init__(self, reject):
        self.reject = reject
        self.evaluated = []

    def validate_reference(self, pose):
        pass

    def evaluate(self, pose, cid):
        self.evaluated.append(cid)
        # The scorer accepts an assessment object; annotations must not leak into
        # the independent design pose when its alignment is reused.
        pose.SetProp("pocket_annotation", str(cid))
        return {"passed": cid not in self.reject, "candidate_conformer_id": cid}


def test_shared_reference_keeps_pocket_rejection_order_diagnostics_and_pose_isolation():
    product, reference = "CCCCCCO", "CCCCCO"
    ids = list(module.ensemble(canonical(product), SETTINGS).relative_energies)
    assert len(ids) > 1
    old_assessor, assessor = PoseAssessment({ids[0]}), PoseAssessment({ids[0]})
    baseline, double = evaluate_with_trace(product, reference, shared=False, constraints=old_assessor)
    actual, single = evaluate_with_trace(product, reference, shared=True, constraints=assessor)
    assert assessor.evaluated == old_assessor.evaluated == ids
    assert len(single) == len(ids) and len(double) == 2 * len(ids) - 1
    assert actual[0]["pocket_selection"]["rejected_conformations"] == [
        {"candidate_conformer_id": ids[0], "passed": False}
    ]
    assert actual[1].HasProp("pocket_annotation")
    assert not actual[2].HasProp("pocket_annotation")
    assert_same_results(actual, baseline)


def test_shared_reference_does_not_turn_pocket_rejection_into_success():
    product, reference = "CCCCCCO", "CCCCCO"
    ids = set(module.ensemble(canonical(product), SETTINGS).relative_energies)
    messages = []
    for shared in (False, True):
        with pytest.raises(ValueError, match="No candidate conformation satisfies") as error:
            evaluate_with_trace(product, reference, shared=shared, constraints=PoseAssessment(ids))
        messages.append(str(error.value))
    assert messages[0] == messages[1]
