import argparse
import pandas as pd



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate INSDC accession numbers to PPX format.")
    parser.add_argument("--metadata", type=str, required=True, help="Path to the metadata file.")
    parser.add_argument("--accessions", type=str, required=True, help="Path to the file containing INSDC accession numbers.")
    parser.add_argument("--output", type=str, required=True, help="Path to the output file.")

    args = parser.parse_args()

    # Read the metadata file
    metadata = pd.read_csv(args.metadata, sep="\t")
    # LOOKUP dict
    lookup = {k.split('.')[0]:v for k, v in zip(metadata["INSDC_accession"], metadata["accession"]) if pd.notna(k)}

    # Read the input file
    with open(args.accessions, "r") as f:
        accessions = f.read().splitlines()

    # Translate accessions
    translated_accessions = []
    for acc in accessions:
        if acc in lookup:
            translated_accessions.append(lookup[acc])
        else:
            translated_accessions.append(acc)  # Keep original if not found

    # Write the output file
    with open(args.output, "w") as f:
        for acc in translated_accessions:
            f.write(f"{acc}\n")
    print(f"Translated accessions written to {args.output}")