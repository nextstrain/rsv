#!/usr/bin/env python3
"""Score F protein sequences based on DMS antibody escape data.

This script can score sequences in two modes:
1. fasta: Score sequences from a FASTA file and output to TSV
2. tree: Score tree nodes using cumulative mutations and output to Auspice JSON
"""

import argparse
import json

import Bio.SeqIO
import pandas as pd
from Bio import Phylo


class SequenceScorer:
    """Assign scores to lists of mutations, mutation can be between any two amino acids.

    Parameters
    ----------
    mutation_effects_df : pandas.DataFrame
        Should have columns `site`, `wildtype`, `mutant`, `mutation_effect`.
    only_positive_escape : bool
        If True, clip mutation effects to be >= 0.
    """

    def __init__(self, mutation_effects_df, only_positive_escape):
        mutation_effects_df = mutation_effects_df[
            ["site", "wildtype", "mutant", "mutation_effect"]
        ].dropna()
        assert len(mutation_effects_df) == len(
            mutation_effects_df[["site", "mutant"]].drop_duplicates()
        )
        self.sites = sorted(set(mutation_effects_df["site"]))
        assert len(self.sites) == len(
            mutation_effects_df[["site", "wildtype"]].drop_duplicates()
        )
        self.wts = mutation_effects_df.set_index("site")["wildtype"].to_dict()
        if only_positive_escape:
            mutation_effects_df = mutation_effects_df.assign(
                mutation_effect=lambda x: x["mutation_effect"].clip(lower=0)
            )
        self.effects = {
            site: site_df.set_index("mutant")["mutation_effect"].to_dict()
            for site, site_df in mutation_effects_df.groupby("site")
        }
        for site, wt in self.wts.items():
            assert wt not in self.effects[site] or self.effects[site][wt] == 0
            self.effects[site][wt] = 0.0

    def total_escape(self, muts):
        """Returns total escape for list of `muts` as `(wildtype, site, mutant)`."""
        pheno = 0.0
        for wt, site, m in muts:
            if (
                (site in self.effects)
                and (wt in self.effects[site])
                and (m in self.effects[site])
            ):
                pheno += self.effects[site][m] - self.effects[site][wt]
        return pheno

    def max_mut_effect(self, muts):
        """Returns `(effect, mutation)` for max effect mutation > 0 for list of muts."""
        effect = 0
        mutation = ""
        for wt, site, m in muts:
            if (
                (site in self.effects)
                and (wt in self.effects[site])
                and (m in self.effects[site])
            ):
                mut_effect = self.effects[site][m] - self.effects[site][wt]
                if mut_effect > effect:
                    effect = mut_effect
                    mutation = f"{wt}{site}{m}"
        return (effect, mutation)


def load_dms_scores(dms_scores_file, antibodies):
    """Load and validate DMS scores CSV file.

    Parameters
    ----------
    dms_scores_file : str
        Path to CSV file with DMS scores.
    antibodies : list of str
        List of antibody column names to extract.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: wildtype, site, mutant, and antibody columns.
    """
    print(f"Reading DMS scores from {dms_scores_file}")
    dms_scores = pd.read_csv(dms_scores_file)
    req_cols = ["wildtype", "site", "mutant"]
    if set(req_cols).intersection(antibodies):
        raise ValueError(f"{antibodies=} cannot contain {req_cols}")
    if not set(antibodies + req_cols).issubset(dms_scores.columns):
        raise ValueError(f"DMS scores file lacks {req_cols} or {antibodies=}")
    dms_scores = dms_scores[req_cols + antibodies]
    print(f"Read scores for {antibodies} for {len(dms_scores)=} mutations")
    return dms_scores


def score_sequences_from_fasta(
    sequences_file, dms_scores, antibodies, max_muts, only_positive_escape
):
    """Score sequences from a FASTA file.

    Parameters
    ----------
    sequences_file : str
        Path to FASTA file with aligned sequences.
    dms_scores : pandas.DataFrame
        DMS scores with columns: site, wildtype, mutant, and antibody columns.
    antibodies : list of str
        List of antibody names to score.
    max_muts : int
        Maximum number of mutations allowed per sequence.
    only_positive_escape : bool
        If True, only consider positive escape values.

    Returns
    -------
    pandas.DataFrame
        DataFrame with scores for each sequence.
    """
    print(f"Reading protein sequences from {sequences_file}")
    seqs = {}
    for s in Bio.SeqIO.parse(sequences_file, "fasta"):
        if s.id in seqs:
            raise ValueError(f"Duplicate {s.id=} in {s}")
        seqs[s.id] = str(s.seq)
    assert seqs, f"no sequences in {sequences_file}"
    seqlens = list(set(len(s) for s in seqs.values()))
    if len(seqlens) != 1:
        raise ValueError(f"Sequences not all the same length:\n{seqlens=}")
    seqlen = seqlens[0]
    print(f"Read {len(seqs)=} sequences with {seqlen=}")

    # Get DMS parent sequence with X for missing amino acids
    assert dms_scores["site"].nunique() == len(
        dms_scores[["site", "wildtype"]].drop_duplicates()
    )
    dms_wts = dms_scores.set_index("site")["wildtype"].to_dict()
    dms_seq = "".join(
        dms_wts[r] if r in dms_wts else "X" for r in range(1, seqlen + 1)
    ).upper()
    assert len(dms_seq) == seqlen

    # Create dataframe including tuples with mutations
    seqs_df = (
        pd.Series(seqs)
        .rename("sequence")
        .rename_axis("accession")
        .reset_index()
        .assign(
            sequence=lambda x: x["sequence"].str.upper(),
            mutations=lambda x: x["sequence"].map(
                lambda s: tuple(
                    (wt, r, mut)
                    for r, (wt, mut) in enumerate(zip(dms_seq, s, strict=True), start=1)
                    if (wt != mut) and not {wt, mut}.intersection({"X", "-"})
                )
            ),
            n_mutations=lambda x: x["mutations"].map(len),
        )
    )

    if any(seqs_df["n_mutations"] > max_muts):
        raise ValueError(
            f"Some sequences have > {max_muts=} mutations:\n"
            + str(
                seqs_df[["accession", "n_mutations", "mutations"]].query(
                    "n_mutations > @max_muts"
                )
            )
        )

    # Calculate effects on data frame of mutations then merge back into sequence data frame;
    # this is more efficient as it avoids duplicate calculations on identical proteins
    print("Scoring sequences...")
    muts_df = seqs_df[["mutations"]].drop_duplicates()
    dms_cols = []
    for ab in antibodies:
        scorer = SequenceScorer(
            dms_scores[["site", "wildtype", "mutant", ab]].rename(
                columns={ab: "mutation_effect"}
            ),
            only_positive_escape,
        )
        muts_df[f"{ab}_total_escape"] = muts_df["mutations"].map(scorer.total_escape)
        muts_df[f"{ab}_max_escape"], muts_df[f"{ab}_max_escape_mutation"] = zip(
            *muts_df["mutations"].map(scorer.max_mut_effect)
        )
        dms_cols += [
            f"{ab}_total_escape",
            f"{ab}_max_escape",
            f"{ab}_max_escape_mutation",
        ]

    df = seqs_df.merge(muts_df, on="mutations", validate="many_to_one").assign(
        mutations_from_DMS_strain=lambda x: x["mutations"].map(
            lambda ms: ",".join(f"{wt}{r}{m}" for (wt, r, m) in ms)
        )
    )[["accession"] + dms_cols + ["mutations_from_DMS_strain"]]
    assert len(df) == len(seqs_df)

    return df


def parse_aa_mutations(aa_muts_json, gene):
    """Parse amino acid mutations from augur ancestral JSON for a specific gene.

    Parameters
    ----------
    aa_muts_json : str
        Path to aa_muts JSON file from augur translate.
    gene : str
        Gene name to extract mutations for.

    Returns
    -------
    dict
        Mapping of node name to list of (wildtype, site, mutant) tuples
        representing branch mutations.
    """
    with open(aa_muts_json) as f:
        aa_data = json.load(f)

    node_mutations = {}

    for node, node_data in aa_data["nodes"].items():
        if gene in node_data.get("aa_muts", {}):
            # Parse mutations like "A123B" into (wildtype, site, mutant)
            mutations = []
            for mut_str in node_data["aa_muts"][gene]:
                wt = mut_str[0]
                mutant = mut_str[-1]
                site = int(mut_str[1:-1])
                mutations.append((wt, site, mutant))
            node_mutations[node] = mutations
        else:
            # No mutations for this node
            node_mutations[node] = []

    return node_mutations


def compute_cumulative_mutations(aa_muts_json, tree_newick, gene):
    """Compute cumulative mutations from root for each node.

    Parameters
    ----------
    aa_muts_json : str
        Path to aa_muts JSON file from augur translate.
    tree_newick : str
        Path to Newick tree file from augur refine.
    gene : str
        Gene name to compute mutations for.

    Returns
    -------
    dict
        Mapping of node name to list of (wildtype, site, mutant) tuples
        representing cumulative mutations from root.
    """
    # Get mutations along branches
    branch_mutations = parse_aa_mutations(aa_muts_json, gene)

    # Load tree structure from Newick file
    tree = Phylo.read(tree_newick, "newick")

    # Build parent map from tree
    parent_map = {}
    for clade in tree.find_clades():
        for child in clade.clades:
            assert child.name not in parent_map, child.name
            parent_map[child.name] = clade.name

    # Compute cumulative mutations for each node
    cumulative_mutations = {}

    def get_cumulative_muts(node_name):
        if node_name in cumulative_mutations:
            return cumulative_mutations[node_name]
        if node_name not in parent_map:
            assert node_name == tree.root.name
            cumulative_mutations[node_name] = []
            return []
        parent_muts = get_cumulative_muts(parent_map[node_name])
        branch_muts = branch_mutations.get(node_name, [])
        # Combine mutations: parent cumulative + branch mutations
        # Need to handle if a mutation reverts a previous one
        muts_dict = {}
        for wt, site, mut in parent_muts:
            muts_dict[site] = (wt, site, mut)
        for wt, site, mut in branch_muts:
            if site in muts_dict:
                orig_wt = muts_dict[site][0]
                if mut == orig_wt:
                    del muts_dict[site]
                else:
                    muts_dict[site] = (orig_wt, site, mut)
            else:
                muts_dict[site] = (wt, site, mut)
        cumulative_mutations[node_name] = list(muts_dict.values())
        return cumulative_mutations[node_name]

    for node_name in branch_mutations.keys():
        get_cumulative_muts(node_name)

    return cumulative_mutations


def score_tree_nodes(
    tree_newick, aa_muts_json, gene, dms_scores, antibodies, only_positive_escape
):
    """Score tree nodes based on cumulative mutations.

    Parameters
    ----------
    tree_newick : str
        Path to Newick tree file from augur refine.
    aa_muts_json : str
        Path to aa_muts JSON file from augur translate.
    gene : str
        Gene name to score.
    dms_scores : pandas.DataFrame
        DMS scores with columns: site, wildtype, mutant, and antibody columns.
    antibodies : list of str
        List of antibody names to score.
    only_positive_escape : bool
        If True, only consider positive escape values.

    Returns
    -------
    dict
        Auspice node data JSON with scores for each node.
    """
    print(f"Reading amino acid mutations from {aa_muts_json}")
    print(f"Reading tree from {tree_newick}")

    # Compute cumulative mutations using Newick tree and aa_muts JSON
    cumulative_mutations = compute_cumulative_mutations(aa_muts_json, tree_newick, gene)

    print(f"Computed cumulative mutations for {len(cumulative_mutations)} nodes")

    # Create node data structure for Auspice
    node_data = {"nodes": {}}

    # Compute scores for each node
    print("Computing F scores for each node...")
    for ab in antibodies:
        print(f"  Processing {ab}...")
        scorer = SequenceScorer(
            dms_scores[["site", "wildtype", "mutant", ab]].rename(
                columns={ab: "mutation_effect"}
            ),
            only_positive_escape,
        )

        for node_name, mutations in cumulative_mutations.items():
            if node_name not in node_data["nodes"]:
                node_data["nodes"][node_name] = {}

            # Compute total escape score and top escape mutation
            total_escape = scorer.total_escape(mutations)
            max_escape, max_escape_mutation = scorer.max_mut_effect(mutations)

            # Store as node attribute
            node_data["nodes"][node_name][f"{ab}_total_escape"] = total_escape
            node_data["nodes"][node_name][f"{ab}_max_escape"] = max_escape
            node_data["nodes"][node_name][
                f"{ab}_max_escape_mutation"
            ] = max_escape_mutation

    return node_data


def main():
    """Main entry point with subcommand parsing."""
    parser = argparse.ArgumentParser(
        description="Score F protein sequences based on DMS antibody escape data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
        help="Scoring mode",
        metavar="{fasta,tree}",
    )

    # Subcommand for scoring sequences from FASTA
    fasta_parser = subparsers.add_parser(
        "fasta",
        help="Score sequences from FASTA file and output to TSV",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    fasta_parser.add_argument(
        "--sequences",
        required=True,
        help="Input FASTA file with aligned F protein sequences",
    )
    fasta_parser.add_argument(
        "--output", required=True, help="Output TSV file with scores for each sequence"
    )
    fasta_parser.add_argument(
        "--dms-scores",
        required=True,
        help="Input CSV with per-mutation DMS antibody escape effects (columns: site, wildtype, mutant, <antibody names>)",
    )
    fasta_parser.add_argument(
        "--dms-antibodies",
        required=True,
        nargs="+",
        help="List of antibody names (column names in --dms-scores file) to score",
    )
    fasta_parser.add_argument(
        "--only-positive-escape",
        required=True,
        choices=["true", "false", "True", "False"],
        help="Only consider positive escape values (clip negative values to 0)",
    )

    # Subcommand for scoring tree nodes
    tree_parser = subparsers.add_parser(
        "tree",
        help="Score tree nodes using cumulative mutations and output Auspice JSON",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    tree_parser.add_argument(
        "--tree-newick",
        required=True,
        help="Input Newick tree file from augur refine",
    )
    tree_parser.add_argument(
        "--aa-muts",
        required=True,
        help="Input amino acid mutations JSON file from augur translate",
    )
    tree_parser.add_argument(
        "--gene",
        required=True,
        help="Gene name to compute scores for (e.g., 'F')",
    )
    tree_parser.add_argument(
        "--output",
        required=True,
        help="Output JSON file with Auspice node data containing scores",
    )
    tree_parser.add_argument(
        "--dms-scores",
        required=True,
        help="Input CSV with per-mutation DMS antibody escape effects (columns: site, wildtype, mutant, <antibody names>)",
    )
    tree_parser.add_argument(
        "--dms-antibodies",
        required=True,
        nargs="+",
        help="List of antibody names (column names in --dms-scores file) to score",
    )
    tree_parser.add_argument(
        "--only-positive-escape",
        required=True,
        choices=["true", "false", "True", "False"],
        help="Only consider positive DMS escape values (clip negative values to 0)",
    )

    args = parser.parse_args()

    # Convert string to boolean
    only_positive_escape = args.only_positive_escape.lower() == "true"

    # Load DMS scores (common to both modes)
    dms_scores = load_dms_scores(args.dms_scores, args.dms_antibodies)

    if args.command == "fasta":
        # Score sequences from FASTA
        max_muts = 100  # raise error if any sequence has more than this many mutations
        scores = score_sequences_from_fasta(
            args.sequences,
            dms_scores,
            args.dms_antibodies,
            max_muts,
            only_positive_escape,
        )
        print(f"Writing scores to {args.output}")
        scores.to_csv(args.output, index=False, sep="\t", float_format="%.5g")
        print(f"Successfully scored {len(scores)} sequences")

    elif args.command == "tree":
        # Score tree nodes
        node_data = score_tree_nodes(
            args.tree_newick,
            args.aa_muts,
            args.gene,
            dms_scores,
            args.dms_antibodies,
            only_positive_escape,
        )
        print(f"Writing node data to {args.output}")
        with open(args.output, "w") as f:
            json.dump(node_data, f, indent=2)
        print(f"Successfully computed F scores for {len(node_data['nodes'])} nodes")

    else:
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
