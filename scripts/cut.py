from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import Seq
import argparse

def cut(oldalignment, newalignment, referencefile, gene, min_length=0):
        list1 = []
        list2 = []
        if gene == 'genome': gene = 'G'

        ref = SeqIO.read(referencefile, "genbank")
        for feature in ref.features:
            if feature.type =='gene':
                a =str((list(feature.qualifiers.items())[0])[-1])[2:-2]
                if a == gene:
                    print(a)
                    startofgene = int(list(feature.location)[0])
                    endofgene =  int(list(feature.location)[-1])+1
                    print(startofgene, endofgene, (endofgene-startofgene))

        alignment =SeqIO.parse(oldalignment, 'fasta')
        for entry in alignment:
            newrecord = SeqRecord(Seq(entry.seq[startofgene:endofgene]), id=entry.id, description=entry.description)
            if len(newrecord.seq.ungap('-')) >= min_length:
                list1.append(newrecord)

        for record in list1:
            sequence =str(record.seq)
            sequence = sequence.replace("-", "")
            if sequence != "":
                list2.append(record)
        SeqIO.write(list2, newalignment, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make new reference depending on whether the entire genome or only part is to be used for the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--oldalignment", required=True, help="fasta input with alignment of entire genome")
    parser.add_argument("--slicedalignment", required=True, help="fasta output with alignment of relevant part of geome")
    parser.add_argument("--reference", required=True, help="reference genbank file of entire genome")
    parser.add_argument("--min-length", default=0, type=int, help="minimal length of output sequence")
    parser.add_argument("--gene", required=True, help="build name, G or F or genome")
    args = parser.parse_args()

    cut(args.oldalignment, args.slicedalignment, args.reference, args.gene, args.min_length)
