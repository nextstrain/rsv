import pandas as pd
import argparse

def metadatadeduplication(old, new):
    df  = pd.read_csv(old, sep='\t')
    df1 = df.drop_duplicates(subset=['strain'], keep='first')
    df1.to_csv(new,sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make new reference depending on whether the entire genome or only part is to be used for the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadataoriginal", required=True, help="original metadata file, tsv format")
    parser.add_argument("--metadataoutput", required=True, help="deduplicated metadata file, tsv format")
    args = parser.parse_args()
    metadatadeduplication(args.metadataoriginal, args.metadataoutput)
