from Bio import SeqIO
import csv
import argparse

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


def sequencesandmetadata(sortedalignment,allmetadata, newmetadata, allsequences, newsequences):

    lista = []
    seqa =[]
    listofid = []

    with open(sortedalignment) as file:
        a = SeqIO.parse(file, "fasta")
        for i in a:
            listofid.append(i.id)

    with open(allmetadata) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        row1= next(tsv_file)
        lista.append(row1)
        for line in tsv_file:
            if line[0] in listofid:
                lista.append(line)

    with open(newmetadata, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for a in lista:
            writer.writerow(a)

    with open(allsequences) as file:
        b = SeqIO.parse(file, "fasta")
        for i in b:
            if(i.id) in listofid:
                seqa.append(i)
    SeqIO.write(seqa, newsequences,"fasta")

print(sequencesandmetadata(args.sortedalignment, args.allmetadata, args.metadata, args.allsequences, args.sequences))
