from Bio import SeqIO
import pandas as pd
import argparse

def sequencesandmetadata(sortedalignment,allmetadata, allsequences,newmetadata, newsequences):
    seq_new, listofid =([] for i in range(2))
    alignment = SeqIO.parse(sortedalignment, "fasta")

    for record in alignment: listofid.append(record.id)

    tsv_file = pd.read_csv(allmetadata, sep="\t")
    metadata = pd.DataFrame(data =tsv_file.loc[tsv_file['accession'].isin(listofid)], columns=tsv_file.columns)
    metadata.to_csv(newmetadata, sep="\t")
    sequences = SeqIO.parse(allsequences, "fasta")
    for record in sequences:
        if record.id in listofid: seq_new.append(record)
    SeqIO.write(seq_new, newsequences,"fasta")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="write separate files for sequences and metadata for a and b",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("--sortedalignment", required=True, help="FASTA file of aligned and sorted sequences")
    parser.add_argument("--allmetadata", required=True, help="all metadata input file")
    parser.add_argument("--allsequences", required=True, help="all sequences input file")
    parser.add_argument("--metadata", required=True, help="output metadata file for a or b, tsv")
    parser.add_argument("--sequences", required=True, help="output sequences file for a or b, FASTA")
    args = parser.parse_args()
        
    sequencesandmetadata(args.sortedalignment, args.allmetadata, args.allsequences,  args.metadata, args.sequences)
