from scaffoldgraph import *
from scaffoldgraph.core.fragment import *
# from scaffoldgraph.core.scaffold import Scaffold
# from scaffoldgraph.utils import suppress_rdlogger

@suppress_rdlogger()
def get_all_murcko_fragments(mol, break_fused_rings=True,iteration_round:int=float('inf')):
    if break_fused_rings:
        fragmenter = MurckoRingFragmenter()
    else:
        fragmenter = MurckoRingSystemFragmenter()
    mol = get_murcko_scaffold(mol)
    rdmolops.RemoveStereochemistry(mol)
    scaffold = Scaffold(mol)
    parents = {scaffold}

    def recursive_generation(child):
        it_cnt = 0
        for parent in fragmenter.fragment(child):
            it_cnt += 1
            if parent in parents:
                continue
            parents.add(parent)
            if it_cnt > iteration_round:
                break
            recursive_generation(parent)

    recursive_generation(scaffold)
    return [f.mol for f in parents]
