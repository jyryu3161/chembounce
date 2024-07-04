#! usr/bin/env python3
"""
Modified module of scaffoldgraph
github)
https://github.com/UCLCheminformatics/ScaffoldGraph
citation)
@article{10.1093/bioinformatics/btaa219,
    author = {Scott, Oliver B and Chan, A W Edith},
    title = "{ScaffoldGraph: an open-source library for the generation and analysis of molecular scaffold networks and scaffold trees}",
    journal = {Bioinformatics},
    year = {2020},
    month = {03},
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btaa219},
    url = {https://doi.org/10.1093/bioinformatics/btaa219},
    note = {btaa219}
    eprint = {https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btaa219/32984904/btaa219.pdf},
}
"""
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
