import numpy as np
import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="sort into only sequences more than 90 percent similar to reference file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="FASTA file of aligned sequences")
    parser.add_argument("--referencefile", required=True, help="genbank reference file")
    parser.add_argument("--output", required=True, help="output sorted alignment fasta, more than 0.9 similar to reference")
    args = parser.parse_args()
    

def sorted_file(referencefile, alignedfile):
    
    refseq = SeqIO.parse(args.referencefile, "genbank")
    seq2 = np.array(refseq.__next__())
    list1 =[]
    for record in SeqIO.parse(args.alignment,"fasta"):

        seq1 = np.array(record.seq)
        good_indices = (refseq!='N')&(seq1!='N')
        a = np.mean(seq1[good_indices]==seq2[good_indices])

        if a >0.9: list1.append(record)
        else: continue
    SeqIO.write(list1, args.output,"fasta")
print(sorted_file(args.referencefile, args.alignment))


   
