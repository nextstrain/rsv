import pandas as pd
import argparse


def metadata_deduplication(old, new):
    df  = pd.read_csv(old, sep='\t')
    df1 = df.drop_duplicates(subset=['strain'], keep='first')
    df1.to_csv(new,sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make new reference depending on whether the entire genome or only part is to be used for the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata-original", required=True, help="original metadata file, tsv format")
    parser.add_argument("--metadata-output", required=True, help="deduplicated metadata file, tsv format")
    args = parser.parse_args()
    metadata_deduplication(args.metadata_original, args.metadata_output)
