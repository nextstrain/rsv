import numpy as np
import pandas as pd
from Bio import SeqIO
import csv

listofa, lista, ids_a, listofb, ids_b, listb, alltherest, list1, list2 = ([] for i in range(9))

seq_a = ["data/a/newsequences.aligned.fasta", "data/a/middlesequences.aligned.fasta", "data/a/oldsequences.aligned.fasta"]
unaligned_a = ["data/a/newsequences.fasta", "data/a/middlesequences.fasta", "data/a/oldsequences.fasta"]
unaligned_b = ["data/b/newsequences.fasta", "data/b/middlesequences.fasta", "data/b/oldsequences.fasta"]
seq_b = ["data/b/newsequences.aligned.fasta", "data/b/middlesequences.aligned.fasta", "data/b/oldsequences.aligned.fasta"]
references_a = ["config/anewreference.fasta", "config/amiddlereference.fasta", "config/aoldreference.fasta"]
references_b = ["config/bnewreference.fasta", "config/bmiddlereference.fasta", "config/boldreference.fasta"]
accessions_a = [pd.read_csv("data/a/newmetadata.tsv", sep='\t').accession, pd.read_csv("data/a/middlemetadata.tsv", sep='\t').accession, pd.read_csv("data/a/oldmetadata.tsv", sep='\t').accession] 
accessions_b = [pd.read_csv("data/b/newmetadata.tsv", sep='\t').accession, pd.read_csv("data/b/middlemetadata.tsv", sep='\t').accession, pd.read_csv("data/b/oldmetadata.tsv", sep='\t').accession]

everything_a = (set().union(*accessions_a))
everything_b = (set().union(*accessions_b))
everything = everything_a.union(everything_b)

for file_a, file_b,  a_fasta, b_fasta in zip(seq_a, seq_b, references_a, references_b):

    sequences = [file_a, file_b]
    both = list(everything_a.intersection(everything_b))

    for file in sequences:
        for record in SeqIO.parse(file, "fasta"):
            if record.id not in both:
                if record.id in everything_a: listofa.append(record.id)
                if record.id in everything_b: listofb.append(record.id)
            else: alltherest.append(record.id)

    refseq_a = SeqIO.read(a_fasta, "fasta").seq
    seq2 = np.array(refseq_a)
    refseq_b = SeqIO.read(b_fasta, "fasta").seq
    seq3 = np.array(refseq_b)

    for record in SeqIO.parse(file_a,"fasta"):
        if record.id in alltherest:
            seq1 = np.array(record.seq)
            good_indices = (refseq_a!='-')&(seq1!='-')
            a = np.mean(seq1[good_indices]==seq2[good_indices])
            lista.append(a)
            ids_a.append(record.id)

    for record in SeqIO.parse(file_b,"fasta"):
        if record.id in alltherest:
            seq1 = np.array(record.seq)
            good_indices = (refseq_b!='-')&(seq1!='-')
            b = np.mean(seq1[good_indices]==seq3[good_indices])
            listb.append(b)
            ids_b.append(record.id)

dict_a, dict_b = (dict() for i in range(2))

dict_a["a"] = lista
dict_a["id a"] = ids_a
dict_b["b"] = listb
dict_b["id b"] = ids_b

df_b = pd.DataFrame(dict_b)
df_a = pd.DataFrame(dict_a)
df_a = (df_a.drop_duplicates(subset=["id a"])).sort_values("id a")
df_b = (df_b.drop_duplicates(subset=["id b"])).sort_values("id b")

for i, j, a, b in zip(df_a["id a"], df_b["id b"], df_a["a"], df_b["b"]):
    if i == j:
        if a > b: listofa.append(i)
        if b > a: listofb.append(i)

#writing new sequences files for a and b 

for file_a, file_b in zip(unaligned_a, unaligned_b):
    for record in SeqIO.parse(file_a,"fasta"):
        if record.id in listofa: list1.append(record)

    for record in SeqIO.parse(file_b, "fasta"):
        if record.id in listofb: list2.append(record)
        
SeqIO.write(list1, "data/a/sequences_notdedup.fasta","fasta")
SeqIO.write(list2, "data/b/sequences_notdedup.fasta", "fasta")

for type in ["a", "b"]:
    listofid = []
    metadata_ = []
    with open(f"data/{type}/sequences_notdedup.fasta") as file:
        parse_file = SeqIO.parse(file, "fasta")
        for i in parse_file:
            listofid.append(i.id)

#writing new metadata files for a and b

    with open("data/metadata.tsv") as file:
        tsv_file = csv.reader(file, delimiter="\t")
        row1= next(tsv_file)
        metadata_.append(row1)
        for line in tsv_file:
            if line[0] in listofid: metadata_.append(line)

    with open(f"data/{type}/metadata_notdedup.tsv", 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for i in metadata_: writer.writerow(i)