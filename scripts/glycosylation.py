#!/usr/bin/python

import argparse, json
import numpy as np
from Bio import Phylo, AlignIO, SeqIO
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="find where sequences are glycosylated",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="FASTA file of aligned sequences")
    parser.add_argument("--tree", required=True, help="NWK file, phylogenetic tree")
    parser.add_argument("--output", required=True, help="JSON file, places where glycosylation has occured")

    args = parser.parse_args()
    

def glycosylation_count(total_aa_seq, glyc_mask=None):
    if glyc_mask is None:
        glyc_mask = np.ones(len(total_aa_seq), dtype=bool)

    total_aa_seq_masked = "".join([aa if mask else 'X'
                                   for (mask, aa) in zip(glyc_mask, total_aa_seq)])

    return len(re.findall('N[^P][ST][^P]', total_aa_seq_masked))


with open (args.alignment) as f:
    file = SeqIO.parse(f, 'fasta')  

tree = args.tree
alignment = args.alignment
T = Phylo.read(tree, 'newick')

glyc_json = {}
aln = {s.name:str(s.seq) for s in AlignIO.read(alignment, 'fasta')}
root_seq = aln[T.root.name]
root_glyc = glycosylation_count(root_seq)
for n in T.find_clades(order='preorder'):
    glyc_json[n.name] = {'glyc':glycosylation_count(aln[n.name]) - root_glyc}

with open(args.output, 'wt') as fh:
    json.dump({'nodes':glyc_json}, fh)
