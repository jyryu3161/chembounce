#!/usr/bin/env python3
"""CLI options"""
import argparse

#### CLI options ####

# Argument of parser
def argument_parser():
    parser = argparse.ArgumentParser()    
    parser.add_argument('-o', '--output_dir', required=True, help="Output location")
    parser.add_argument('-i', '--input_smiles', required=True, help="Input SMILES, the target molecular structure")
    # Max numbers to test
    parser_cnt = parser.add_argument_group(title='iteration number options',description='Limit numbers of iterations')
    parser_cnt.add_argument('--overall_max_n', required=False, default=None, type=int,
                        help="Maximal number of scaffold-hopped candidates for overall fragments")
    parser_cnt.add_argument('-n','--frag_max_n', required=False, default=1000, type=int,
                        help="Maximal number of scaffold-hopped candidates for a fragment")
    parser_cnt.add_argument('--scaffold_top_n', required=False, default=None, type=int,
                        help="Number of scaffolds to test for a fragment")
    parser_cnt.add_argument('--cand_max_n__rplc', required=False, default=10, type=int,
                        help="Maximal number of candidates for a replaced scaffold")
#     parser_cnt.add_argument('--merge_structure_top_n', required=False, default=100, help="Number of top fragments to test", type=int)
    parser_thr = parser.add_argument_group(title='threshold options',description='Thresholds for candidates')
    parser_thr.add_argument('-t', '--tanimoto_threshold', required=False, default=0.5, type=float,
                        help="Similarity threshold, between 0 and 1: used to exclude irrelated molecular structure, based on the similarity between the original structure and scaffold-hopped one. Default is 0.5")
    # Thresholds on molecular properties
    parser_thr.add_argument('--qed_min', required=False, default=None, help="Minimal QED score", type=float)
    parser_thr.add_argument('--qed_max', required=False, default=None, help="Maximal QED score", type=float)
    parser_thr.add_argument('--sa_min', required=False, default=None, help="Minimal SAscore (0.0 ~ 10.0)", type=float)
    parser_thr.add_argument('--sa_max', required=False, default=None, help="Maximal SAscore (0.0 ~ 10.0)", type=float)
    parser_thr.add_argument('--logp_min', required=False, default=None, help="Minimal logP limit", type=float)
    parser_thr.add_argument('--logp_max', required=False, default=None, help="Maximal logP limit", type=float)
    parser_thr.add_argument('--mw_min', required=False, default=None, help="Minimal molecular weight (0.0 ~ )", type=float)
    parser_thr.add_argument('--mw_max', required=False, default=None, help="Maximal molecular weight (0.0 ~ )", type=float)
    parser_thr.add_argument('--h_donor_min', required=False, default=None, help="Minimal H donor limit", type=float)
    parser_thr.add_argument('--h_donor_max', required=False, default=None, help="Maximal H donor limit", type=float)
    parser_thr.add_argument('--h_acceptor_min', required=False, default=None, help="Minimal H acceptor limit", type=float)
    parser_thr.add_argument('--h_acceptor_max', required=False, default=None, help="Maximal H acceptor limit", type=float)
    parser_thr.add_argument('--wo_lipinski', required=False, action='store_true', default=False,
                        help="Turn off  Lipinski\'s rule of five: logp_max=5, qed_max=500, h_donor_max=5, h_acceptor_max=10")
    # Minor options for process
    parser.add_argument('-l', '--low_mem', required=False, action='store_true', default=False,
                        help="Low memory mode")
    parser.add_argument('-f', '--fragments', required=False, action='append', default=[],
                        help="Specific fragment SMILES of the input structure to replace. For multiple fragment, repeatedly impose this option: e.g. --fragments SMILES_A --fragments SMILES_B")
    parser.add_argument('--replace_scaffold_files', required=False, action='append', default=[],
                        help="Replace scaffold file, for a specific fragment")
    return parser
