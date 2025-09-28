import argparse
import sys
import pandas as pd

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

if __name__=="__main__":
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

    for gene in ["F", "G"]:
        def get_coverage(d):
            cov_dict = {k:float(v) for k,v in (x.split(":") for x in d['cdsCoverage'].split(",") if ":" in x)}
            return cov_dict.get(gene, 0.0)

        result[f"{gene}_coverage"] = result.apply(get_coverage, axis=1)

    result.to_csv(args.output, index_label=args.id_field, sep='\t')
