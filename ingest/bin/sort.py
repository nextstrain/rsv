import numpy as np
import pandas as pd
from Bio import SeqIO
from functools import reduce

dfs, a_list, b_list, a_sequences, b_sequences, listofdictionaries, a_metadata, b_metadata = ([] for i in range(8))
to_a = ['a_new', 'a_middle', 'a_old']
to_b = ['b_new', 'b_middle', 'b_old']
seq = ["data/a/1_sequences.aligned.fasta", "data/a/2_sequences.aligned.fasta", "data/a/3_sequences.aligned.fasta", "data/b/1_sequences.aligned.fasta", "data/b/2_sequences.aligned.fasta", "data/b/3_sequences.aligned.fasta"]
references = ["config/a_1_reference.fasta", "config/a_2_reference.fasta", "config/a_3_reference.fasta", "config/b_1_reference.fasta", "config/b_2_reference.fasta", "config/b_3_reference.fasta"]
accessions = [pd.read_csv("data/a/1_metadata.tsv", sep='\t').accession, pd.read_csv("data/a/2_metadata.tsv", sep='\t').accession, pd.read_csv("data/a/3_metadata.tsv", sep='\t').accession, pd.read_csv("data/b/1_metadata.tsv", sep='\t').accession, pd.read_csv("data/b/2_metadata.tsv", sep='\t').accession, pd.read_csv("data/b/3_metadata.tsv", sep='\t').accession] 
everything = (set().union(*accessions))

for file_, reference in zip(seq, references):
    dictionary = dict()
    accessions_, values = ([] for i in range(2))
    file = SeqIO.parse(file_,"fasta")
    ref_array = np.array(SeqIO.read(reference,"fasta").seq)

    for record in file:
        record_array = np.array(record.seq)
        good_indices = (ref_array!='-')&(record_array!='-')
        mean = np.mean(record_array[good_indices]==ref_array[good_indices])
        accessions_.append(record.id)
        values.append(mean)

    for accession in everything:
        if accession not in accessions_:
            accessions_.append(accession)
            values.append(0.001)

    dictionary['accession'] = accessions_
    dictionary['value'] = values
    listofdictionaries.append(dictionary)

for i, column in zip(listofdictionaries, (to_a+to_b)):
    df = pd.DataFrame(i).groupby('accession', group_keys=False).apply(lambda x: x.loc[x.value.idxmax()]).iloc[:,1:].reset_index().sort_values('accession').rename(columns={"value": column})
    dfs.append(df)
merged_df = reduce(lambda l, r: pd.merge(l, r, on='accession'), dfs)

df_drop = merged_df.iloc[:, 1:]
first_column = (merged_df.iloc[:, 0])
df_values = (df_drop.idxmax(axis=1)).rename("values")
result = pd.concat([first_column,df_values], axis=1)

for acc, a_or_b in zip(result['accession'], result['values']):
    if a_or_b in to_a: a_list.append(acc)
    elif a_or_b in to_b: b_list.append(acc)

for record in SeqIO.parse("data/sequences.fasta", "fasta"):
    if record.id in a_list: a_sequences.append(record)
    if record.id in b_list: b_sequences.append(record)

SeqIO.write(a_sequences, "data/a/sequences_notdedup.fasta","fasta")
SeqIO.write(b_sequences, "data/b/sequences_notdedup.fasta", "fasta")

tsv_file = pd.read_csv("data/metadata.tsv", sep="\t")
a_metadata = pd.DataFrame(data =tsv_file.loc[tsv_file['accession'].isin(a_list)], columns=tsv_file.columns)
b_metadata = pd.DataFrame(data=tsv_file.loc[tsv_file['accession'].isin(b_list)], columns=tsv_file.columns)
a_metadata.to_csv('data/a/metadata_notdedup.tsv', sep="\t")
b_metadata.to_csv('data/b/metadata_notdedup.tsv', sep="\t")
