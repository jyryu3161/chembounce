"""Deterministic, depth-limited Murcko scaffold enumeration.

Based on ScaffoldGraph (Scott and Chan, Bioinformatics, 2020,
https://doi.org/10.1093/bioinformatics/btaa219).
"""

from rdkit import Chem, rdBase
from scaffoldgraph import get_murcko_scaffold
from scaffoldgraph.core.fragment import MurckoRingFragmenter, MurckoRingSystemFragmenter
from scaffoldgraph.core.scaffold import Scaffold


def get_all_murcko_fragments(mol, break_fused_rings=True, iteration_round=1):
    if iteration_round < 0:
        raise ValueError("Fragment depth must be nonnegative")
    root = get_murcko_scaffold(mol)
    if root.GetNumAtoms() == 0:
        return []
    Chem.RemoveStereochemistry(root)
    fragmenter = MurckoRingFragmenter() if break_fused_rings else MurckoRingSystemFragmenter()
    current = [Scaffold(root)]
    found = {Chem.MolToSmiles(root): root}
    depth = 0
    while current and depth < iteration_round:
        next_level = []
        for child in current:
            Chem.GetSymmSSSR(child.mol)
            # ScaffoldGraph probes an uninitialized RingInfo and repairs it via
            # GetSymmSSSR internally. Silence that known RDKit diagnostic only
            # while it enumerates parents; Python exceptions still propagate.
            with rdBase.BlockLogs():
                parents = fragmenter.fragment(child)
            for parent in parents:
                smiles = Chem.MolToSmiles(parent.mol)
                if smiles and smiles not in found:
                    found[smiles] = parent.mol
                    next_level.append(parent)
        current = sorted(next_level, key=lambda x: Chem.MolToSmiles(x.mol))
        depth += 1
    return [found[key] for key in sorted(found)]
