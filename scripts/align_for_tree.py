from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import shutil
import argparse

def alignfortree(realign, align, reference, newoutput, build):
    records = []
    if build != "genome":
        shutil.copy(realign, newoutput)
    else:
        realigned = {s.id:s for s in SeqIO.parse(realign, "fasta")}
        original = SeqIO.parse(align, "fasta")
        ref = SeqIO.read(reference, "genbank")

        for feature in ref.features:
                if feature.type =='gene':
                    a =str((list(feature.qualifiers.items())[0])[-1])[2:-2]
                    if a == "G":
                        startofgene = int(list(feature.location)[0])
                        endofgene =  int(list(feature.location)[-1])+1
                        break

        for record_original in original:
            sequence_to_insert = realigned.get(record_original.id, None)
            if sequence_to_insert is None:
                sequence_to_insert = '-' * (endofgene - startofgene)
            else:
                sequence_to_insert = sequence_to_insert.seq

            record_for_tree = record_original.seq.replace(record_original.seq[startofgene:endofgene], sequence_to_insert)
            newrecord = SeqRecord(record_for_tree, id=record_original.id, description=record_original.description)
            records.append(newrecord)

        SeqIO.write(records, newoutput, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="insert the G gene realignment into the alignment of the entire genome",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--realign", required=True, help="multiple sequence aligned G gene input FASTA")
    parser.add_argument("--original", required=True, help="entire genome pairwise aligned FASTA")
    parser.add_argument("--reference", required=True, help="original reference GENBANK file of entire genome")
    parser.add_argument("--output", help="output FASTA file")
    parser.add_argument("--gene", help="specify either G or genome")
    args = parser.parse_args()

    alignfortree(args.realign, args.original, args.reference, args.output, args.gene)
