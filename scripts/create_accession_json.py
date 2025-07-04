import json
import argparse
import pandas as pd
from Bio import Phylo



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create json file with node attributes for PPX and genbank accession numbers including links to sequences and data use terms.")
    parser.add_argument('--metadata', type=str, help='Path to the metadata file in TSV format.')
    parser.add_argument('--tree', type=str, help='path to the tree file in Newick format.')
    parser.add_argument('--output-node-data', type=str, help='Path to the output JSON file.')


    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, index_col='accession', sep='\t')
    tips = [n.name for n in Phylo.read(args.tree, "newick").get_terminals()]

    node_data = {}
    for tip in tips:
        print(tip)
        if tip in metadata.index:
            row = metadata.loc[tip]
            node_data[tip] = {
                'PPX_accession': {'value': row.get('accessionVersion', ''), 'url': f"https://pathoplexus.org/seq/{tip}"},
                'dataUseTerms': {'value': row.get('dataUseTerms', ''), 'url': row.get('dataUseTermsUrl', '')},
                'dataUseTerms_URL': {"value": f"{row.get('dataUseTermsUrl', '')}"}
            }
            if row.get('dataUseTerms') == "RESTRICTED":
                node_data[tip]['restrictedUntil'] = {"value": row.get('dataUseTermsRestrictedUntil', '')}
            if not pd.isna(row['genbank_accession']):
                node_data[tip]['insdcAccession'] = {"value": row['genbank_accession'], 'url': f"https://www.ncbi.nlm.nih.gov/nuccore/{row['genbank_accession']}"}


    with open(args.output_node_data, 'w') as f:
        json.dump({'nodes':node_data}, f)

