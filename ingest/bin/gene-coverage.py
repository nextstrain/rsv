from Bio import SeqIO
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import defaultdict



#same as in the sort function

def sequence_to_int_array(s, fill_value=110, fill_gaps=True):
    seq = np.frombuffer(str(s).lower().encode('utf-8'), dtype=np.int8).copy()
    if fill_gaps:
        seq[(seq!=97) & (seq!=99) & (seq!=103) & (seq!=116)] = fill_value
    else:
        seq[(seq!=97) & (seq!=99) & (seq!=103) & (seq!=116) & (seq!=45)] = fill_value
    return seq


#this function finds the coverage of the F and G gene for each sequence in each alignment. The genome annotations used are those of the respective alignment reference.

def G_and_F_coverage(alignment, indices_G, indices_F):
    G_coverage, F_coverage = (dict() for i in range(2))
    for seq in SeqIO.parse(alignment, "fasta"):

        seq_array_G = sequence_to_int_array(seq.seq[indices_G[0]:indices_G[1]], fill_gaps=False)
        other_array_G = sequence_to_int_array(('-'*(indices_G[1]-indices_G[0])), fill_gaps=False)

        seq_array_F = sequence_to_int_array(seq.seq[indices_F[0]:indices_F[1]], fill_gaps=False)
        other_array_F = sequence_to_int_array(('-'*(indices_F[1]-indices_F[0])), fill_gaps=False)

        F_coverage[seq.id] = np.mean(seq_array_F!=other_array_F)
        G_coverage[seq.id] = np.mean(seq_array_G!=other_array_G)
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

    G_coverage, F_coverage = (defaultdict(list) for i in range(2))



    # add the G and F coverage into the dataframe 

    for filename, G_s, F_s in zip(alignments, start_and_end_G, start_and_end_F):
        coverages_G = G_and_F_coverage(filename, G_s, F_s)[0]
        coverages_F = G_and_F_coverage(filename, G_s, F_s)[1]

    for accession in everything:
                G_coverage[accession].append(coverages_G.get(accession, 0.001))
                F_coverage[accession].append(coverages_F.get(accession, 0.001))

    G_coverage_df = pd.DataFrame(G_coverage, index = column).T
    F_coverage_df = pd.DataFrame(F_coverage, index = column).T

    max_F_coverage = F_coverage_df.max(axis=1)
    max_G_coverage = G_coverage_df.max(axis=1)
    F_df = pd.DataFrame({'Accession': max_F_coverage.index, 'Coverage': max_F_coverage.values})
    G_df = pd.DataFrame({'Accession': max_G_coverage.index, 'Coverage': max_G_coverage.values})


    metadata_files = ["data/b/metadata_no_covg.tsv", "data/a/metadata_no_covg.tsv"]
    outputs = ["data/b/metadata.tsv", "data/a/metadata.tsv"]

    # this part of the script reads the already separated a and b metadata and adds the F and G covg values to the correct row based on accession

    for metadata_, output in zip(metadata_files, outputs):
        metadata = pd.read_csv(metadata_, sep='\t')
        F_covg, G_covg = ([] for i in range(2))

        for acc in metadata['accession']:
                Fcovg = (F_df.loc[F_df['Accession'] == acc, 'Coverage'].iloc[0])
                Gcovg = (G_df.loc[G_df['Accession'] == acc, 'Coverage'].iloc[0])
                F_covg.append(Fcovg)
                G_covg.append(Gcovg)

        coverages = pd.DataFrame({'F coverage':F_covg, 'G coverage': G_covg})
        m = pd.DataFrame(metadata.join(coverages))
        m.to_csv(output, sep='\t')