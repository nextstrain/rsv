from Bio import SeqIO
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from sort import sequence_to_int_array


#this function finds the coverage of the F and G gene for each sequence in each alignment. The genome annotations used are those of the respective alignment reference.

def G_and_F_coverage(alignment, indices_G, indices_F):
    gap_symbol = 45
    G_coverage, F_coverage = (dict() for i in range(2))
    for seq in SeqIO.parse(alignment, "fasta"):

        seq_array_G = sequence_to_int_array(seq.seq[indices_G[0]:indices_G[1]], fill_gaps=False)

        seq_array_F = sequence_to_int_array(seq.seq[indices_F[0]:indices_F[1]], fill_gaps=False)

        F_coverage[seq.id] = np.mean(seq_array_F!=gap_symbol)
        G_coverage[seq.id] = np.mean(seq_array_G!=gap_symbol)

    return(G_coverage, F_coverage)


if __name__=="__main__":
    alignments = ["data/a/1_sequences.aligned.fasta", "data/a/2_sequences.aligned.fasta", "data/a/3_sequences.aligned.fasta",
                  "data/b/1_sequences.aligned.fasta", "data/b/2_sequences.aligned.fasta", "data/b/3_sequences.aligned.fasta"]

    start_and_end_G = [ [4682, 5653], [4637,5602], [4633, 5229], [4713, 5645],[4688, 5641], [4615, 5578]]

    start_and_end_F = [[5733, 7457],[5682,7401], [5606, 7349], [5747, 7467], [5718, 7442], [5628, 7529]]

    column = ['a_1', 'a_2', 'a_3', 'b_1', 'b_2', 'b_3']

    accessions = [pd.read_csv("data/a/1_metadata.tsv", sep='\t').accession,
                    pd.read_csv("data/a/2_metadata.tsv", sep='\t').accession,
                    pd.read_csv("data/a/3_metadata.tsv", sep='\t').accession,
                    pd.read_csv("data/b/1_metadata.tsv", sep='\t').accession,
                    pd.read_csv("data/b/2_metadata.tsv", sep='\t').accession,
                    pd.read_csv("data/b/3_metadata.tsv", sep='\t').accession]

    everything = set().union(*accessions)

    # add the G and F coverage into the dataframe
    G_coverage, F_coverage = (defaultdict(list) for i in range(2))
    for filename, G_s, F_s in zip(alignments, start_and_end_G, start_and_end_F):
        coverages_G = G_and_F_coverage(filename, G_s, F_s)[0]
        coverages_F = G_and_F_coverage(filename, G_s, F_s)[1]

        for accession in everything:
            G_coverage[accession].append(coverages_G.get(accession, 0.00))
            F_coverage[accession].append(coverages_F.get(accession, 0.00))

    metadata_files = ["data/b/metadata_no_covg.tsv", "data/a/metadata_no_covg.tsv"]
    outputs = ["data/b/metadata.tsv", "data/a/metadata.tsv"]

    # this part of the script reads the already separated a and b metadata and adds the F and G covg values to the correct row based on accession
    for metadata_fname, output in zip(metadata_files, outputs):
        metadata = pd.read_csv(metadata_fname, sep='\t')
        F_covg, G_covg = [], []

        for acc in metadata['accession']:
            F_covg.append(max(F_coverage[acc]))
            G_covg.append(max(G_coverage[acc]))

        coverages = pd.DataFrame({'F coverage':F_covg, 'G coverage': G_covg})
        m = pd.DataFrame(metadata.join(coverages))
        m.to_csv(output, sep='\t')
