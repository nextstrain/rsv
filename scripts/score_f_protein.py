#!/usr/bin/env python3
"""Score protein sequences from antibody-escape DMS data and output scores to TSV."""


import argparse

import Bio.SeqIO

import pandas as pd


class SequenceScorer:
    """Assign scores to lists of mutations, mutation can be between any two amino acids.

    Parameters
    ----------
    mutation_effects_df : pandas.DataFrame
        Should have columns `site`, `wildtype`, `mutant`, `mutation_effect`.

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

    def phenotype(self, muts):
        """Returns phenotype for list of `muts` as `(wildtype, site, mutant)`."""
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


def score_sequences(seqs, dms_scores, dms_antibodies, max_muts, only_positive_escape):
    """Score the total and top escape from each antibody for all sequences."""

    # get DMS parentl sequence with X for missing amino acids
    seqlen = list(set(len(s) for s in seqs.values()))
    assert len(seqlen) == 1
    seqlen = seqlen[0]
    assert dms_scores["site"].nunique() == len(
        dms_scores[["site", "wildtype"]].drop_duplicates()
    )
    dms_wts = dms_scores.set_index("site")["wildtype"].to_dict()
    dms_seq = "".join(
        dms_wts[r] if r in dms_wts else "X" for r in range(1, seqlen + 1)
    ).upper()
    assert len(dms_seq) == seqlen

    # create dataframe including tuples with mutations
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

    # calculate effects on data frame of mutations then merge back into sequence data frame;
    # this is more efficient as it avoids duplicate calculations on identical proteins
    muts_df = seqs_df[["mutations"]].drop_duplicates()
    dms_cols = []
    for ab in dms_antibodies:
        scorer = SequenceScorer(
            dms_scores[["site", "wildtype", "mutant", ab]].rename(
                columns={ab: "mutation_effect"}
            ),
            only_positive_escape,
        )
        muts_df[f"{ab}_total_escape"] = muts_df["mutations"].map(scorer.phenotype)
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


def main():
    parser = argparse.ArgumentParser(description="Score F protein sequences")
    parser.add_argument(
        "--sequences",
        required=True,
        help="Input FASTA file with aligned F protein sequences",
    )
    parser.add_argument(
        "--dms-scores", required=True, help="Input CSV with per-mutation effects"
    )
    parser.add_argument("--output", required=True, help="Output TSV file with scores")
    parser.add_argument(
        "--dms-antibodies",
        required=True,
        help="List antibodies with escape scores as columns in DMS scores",
        nargs="+",
    )
    parser.add_argument(
        "--only-positive-escape",
        required=True,
        choices=["true", "false", "True", "False"],
        help="Only consider positive escape values (true or false)",
    )

    args = parser.parse_args()

    # Convert string to boolean
    only_positive_escape = args.only_positive_escape.lower() == "true"

    max_muts = 100  # raise error if any sequence has more than this many mutations

    print(f"Reading protein sequences from {args.sequences}")
    seqs = {}
    for s in Bio.SeqIO.parse(args.sequences, "fasta"):
        if s.id in seqs:
            raise ValueError(f"Duplicate {s.id=} in {s}")
        seqs[s.id] = str(s.seq)
    assert seqs, "no sequences in {args.sequences}"
    seqlens = list(set(len(s) for s in seqs.values()))
    if len(seqlens) != 1:
        raise ValueError(f"Sequences not all the same length:\n{seqlens=}")
    seqlen = seqlens[0]
    print(f"Read {len(seqs)=} sequences with {seqlen=}")

    print(f"Reading DMS scores from {args.dms_scores}")
    dms_scores = pd.read_csv(args.dms_scores)
    dms_antibodies = args.dms_antibodies
    req_cols = ["wildtype", "site", "mutant"]
    if set(req_cols).intersection(dms_antibodies):
        raise ValueError(f"{dms_antibodies=} cannot contain {req_cols}")
    if not set(dms_antibodies + req_cols).issubset(dms_scores.columns):
        raise ValueError(
            f"{dms_antibodies.columns=} lacks {req_cols} or {dms_antibodies=}"
        )
    dms_scores = dms_scores[req_cols + dms_antibodies]
    print(f"Read scores for {dms_antibodies} for {len(dms_scores)=} mutations")

    print("Scoring sequences...")
    scores = score_sequences(
        seqs, dms_scores, dms_antibodies, max_muts, only_positive_escape
    )
    print(f"Writing scores to {args.output}")
    scores.to_csv(args.output, index=False, sep="\t", float_format="%.5g")


if __name__ == "__main__":
    main()
