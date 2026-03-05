#!/usr/bin/env python3
"""
Merge F protein scores into metadata TSV file.
This script takes the original metadata and adds F protein scores as a new column.
"""

import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Merge F protein scores into metadata")
    parser.add_argument("--metadata", required=True, help="Input metadata TSV file")
    parser.add_argument(
        "--scores", required=True, help="Input CSV file with F protein scores"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file with merged metadata and scores",
    )
    parser.add_argument(
        "--strain-id-field",
        required=True,
        help='Name of the strain ID field in metadata (e.g., "accession")',
    )

    args = parser.parse_args()

    strain_id = args.strain_id_field
    metadata = pd.read_csv(args.metadata, sep="\t")
    scores = pd.read_csv(args.scores, sep="\t")

    assert strain_id in metadata.columns, f"{strain_id=}, {metadata.columns=}"
    assert strain_id in scores.columns, f"{strain_id=}, {scores.columns=}"
    assert set(scores[strain_id]).issubset(metadata[strain_id])

    assert set(metadata.columns).intersection(scores.columns) == {strain_id}

    metadata.merge(scores, on=strain_id, validate="one_to_one").to_csv(
        args.output, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
