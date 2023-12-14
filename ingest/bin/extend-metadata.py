from Bio import SeqIO
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

NEXTCLADE_JOIN_COLUMN_NAME = 'seqName'
VALUE_MISSING_DATA = '?'

column_map = {
    "clade": "clade",
    "lineage": "lineage",
    "coverage": "genome_coverage",
    "totalMissing": "missing_data",
    "totalSubstitutions": "divergence",
    "totalNonACGTNs": "nonACGTN"
}

coordinates = {'a':{'G':[4652, 5617], 'F':[5697,7421]},
                'b':{'G':[4646, 5578], 'F':[5676,7400]}}

def coverage(target, total):
    if total[0]>target[1] or total[1]<target[0]:
        # to overlap
        return 0
    elif total[0]<=target[0] and total[1]>=target[1]:
        # total overlap
        return 1
    elif total[0]>target[0] and total[1]<target[1]:
        # total contained in target
        return (total[1]-total[0])/(target[1]-target[0])
    elif total[0]>target[0] and total[1]>target[1]:
        # overlap with total to the right of target
        return (target[1]-total[0])/(target[1]-target[0])
    else:
        # overlap with total to the left of target
        return (total[1]-target[0])/(target[1]-target[0])


if __name__=="__main__":
    import argparse, sys
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata")
    parser.add_argument("--nextclade")
    parser.add_argument("--id-field")
    parser.add_argument("--virus-type")
    parser.add_argument("--output", default=sys.stdout)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, index_col=args.id_field,
                           sep='\t', low_memory=False, na_filter = False)

    # Read and rename clade column to be more descriptive
    clades = pd.read_csv(args.nextclade, index_col=NEXTCLADE_JOIN_COLUMN_NAME,
                         sep='\t', low_memory=False, na_filter = False) \
            .rename(columns=column_map)

    # Concatenate on columns
    result = pd.merge(
        metadata, clades,
        left_index=True,
        right_index=True,
        how='left'
    )

    for gene in coordinates[args.virus_type]:
        def get_coverage(d):
            try:
                return coverage(coordinates[args.virus_type][gene], [int(d.alignmentStart), int(d.alignmentEnd)])
            except:
                print('missing alignment for ',d.name)
                return np.nan

        result[f"{gene}_coverage"] = result.apply(get_coverage, axis=1)

    result.to_csv(args.output, index_label=args.id_field, sep='\t')
