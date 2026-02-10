"""Enrich for sequences with high antibody escape scores.

This script selects sequences with high escape scores while maintaining diversity
by limiting the number of sequences with identical F protein mutations or identical
top escape mutations within each group.
"""

import argparse

import augur.dates

import Bio.SeqIO

import pandas as pd


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--input-sequences", required=True, help="Input FASTA file with all sequences"
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Metadata TSV file with escape scores (enhanced metadata from pre-subsample scoring)",
    )
    parser.add_argument(
        "--output-sequences",
        required=True,
        help="Output FASTA file with selected sequences",
    )
    parser.add_argument(
        "--strain-id",
        required=True,
        help="Column name for strain identifier in metadata",
    )
    parser.add_argument(
        "--escape-col",
        required=True,
        help="Column name in metadata with escape scores (e.g., 'Nirsevimab-Fab_total_escape')",
    )
    parser.add_argument(
        "--nseqs", type=int, required=True, help="Number of sequences to select"
    )
    parser.add_argument(
        "--group-by",
        nargs="+",
        required=True,
        help="Metadata columns to group by (e.g., country year)",
    )
    parser.add_argument(
        "--max-identical-f-prot-muts",
        type=int,
        required=True,
        help="Maximum number of sequences with identical F protein mutations per group",
    )
    parser.add_argument(
        "--max-identical-max-escape-mut",
        type=int,
        required=True,
        help="Maximum number of sequences with identical top escape mutation per group",
    )

    args = parser.parse_args()

    print(f"Would select {args.nseqs=} sequences from {args.input_sequences=}")
    print(f"Using escape column: {args.escape_col=}")
    print(f"Grouping by: {args.group_by=}")
    print(
        f"Max identical F protein mutations per group: {args.max_identical_f_prot_muts=}"
    )
    print(
        f"Max identical max escape mutation per group: {args.max_identical_max_escape_mut=}"
    )

    metadata = pd.read_csv(args.metadata, sep="\t").assign(
        year=lambda x: x["date"].map(lambda d: augur.dates.get_year_month_day(d)[0])
    )
    assert args.strain_id in metadata.columns, f"{args.strain_id}, {metadata.columns=}"
    assert (
        args.escape_col in metadata.columns
    ), f"{args.escape_col}, {metadata.columns=}"
    assert set(args.group_by).issubset(
        metadata.columns
    ), f"{args.group_by=}, {metadata.columns=}"
    assert metadata[args.escape_col].notnull().all()
    print(f"Read metadata for {len(metadata)=} strains")

    # get column with max escape for this antibody
    if args.escape_col.endswith("_max_escape"):
        antibody = args.escape_col[: -len("_max_escape")]
    elif args.escape_col.endswith("_total_escape"):
        antibody = args.escape_col[: -len("_total_escape")]
    else:
        raise ValueError(f"Cannot parse antibody from {args.escape_col=}")
    max_escape_mut_col = f"{antibody}_max_escape_mutation"
    assert (
        max_escape_mut_col in metadata.columns
    ), f"{max_escape_mut_col=}, {metadata.columns=}"

    # we will drop grouping variables of NA, so need to make sure empty strings
    # are not that for these columns
    for col in ["mutations_from_DMS_strain", max_escape_mut_col]:
        metadata[col] = metadata[col].fillna("")

    metadata = metadata.groupby(
        args.group_by + ["mutations_from_DMS_strain"], as_index=False, dropna=True
    ).head(n=args.max_identical_f_prot_muts)
    print(f"Retained {len(metadata)=} after {args.max_identical_f_prot_muts=}")

    metadata = (
        metadata.sort_values(args.escape_col, ascending=False)
        .groupby(args.group_by + [max_escape_mut_col], as_index=False, dropna=True)
        .head(n=args.max_identical_max_escape_mut)
    )
    print(f"Retained {len(metadata)=} after {args.max_identical_max_escape_mut=}")

    metadata = metadata.sort_values(args.escape_col, ascending=False).head(n=args.nseqs)
    print(f"Retained {len(metadata)=} after {args.nseqs=}")

    strains_to_keep = {strain: 0 for strain in metadata[args.strain_id]}
    seqs = []
    for seq in Bio.SeqIO.parse(args.input_sequences, "fasta"):
        if seq.id in strains_to_keep:
            seqs.append(seq)
            strains_to_keep[seq.id] += 1

    assert args.strain_id != "n"
    strains_to_keep = (
        pd.Series(strains_to_keep).rename_axis(args.strain_id).rename("n").reset_index()
    )
    if not (strains_to_keep["n"] == 1).all():
        raise ValueError(
            "Did not find exactly one sequence for all strains to keep:\n"
            + str(strains_to_keep.query("n != 1"))
        )

    assert len(metadata) == len(seqs) <= args.nseqs
    Bio.SeqIO.write(seqs, args.output_sequences, "fasta")
    print(f"Wrote output to {args.output_sequences}")


if __name__ == "__main__":
    main()
