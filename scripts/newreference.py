from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, Seq
import shutil
import argparse

def new_reference(greference, referencefile, newfile, gene='genome'):
    ref = SeqIO.read(f"{referencefile}.gbk", "genbank")

    for feature in ref.features:
        if feature.type =='gene':
            a =str((list(feature.qualifiers.items())[0])[-1])[2:-2]
            if a == 'G':
                startofgene = int(list(feature.location)[0])
                endofgene =  int(list(feature.location)[-1])+1
                finish = endofgene-startofgene
    sequence_object = Seq(ref.seq[startofgene:endofgene])

    record = SeqRecord(sequence_object, id=ref.id, name=ref.name, description=ref.description, annotations=ref.annotations)

    for feature in ref.features:
        if feature.type == 'source':
            feature1 =  SeqFeature(FeatureLocation(start=0, end=finish), type='source', qualifiers=feature.qualifiers)
    record.features.append(feature1)

    for feature in ref.features:
        if feature.type == 'CDS':
            a =str((list(feature.qualifiers.items())[0])[-1])[2:-2]
            if a == gene:
                startofgene = int(list(feature.location)[0])
                endofgene =  int(list(feature.location)[-1])+1
                finish = endofgene-startofgene
                feature2 = SeqFeature(FeatureLocation(start=0, end=finish), type='CDS', qualifiers=feature.qualifiers)
                record.features.append(feature2)

    reffasta = SeqIO.read(f"{referencefile}.fasta", 'fasta')
    newrecord = SeqRecord(Seq(reffasta.seq[startofgene:endofgene]), id=reffasta.id, description=reffasta.description)
    SeqIO.write(newrecord, f"{greference}.fasta", "fasta")

    if gene=='genome':
        shutil.copy(f'{referencefile}.gbk', f'{newfile}.gbk')
        shutil.copy(f'{referencefile}.fasta', f'{newfile}.fasta')
    else:
        SeqIO.write(record, f"{newfile}.gbk", 'genbank')
        reffasta = SeqIO.read(f"{referencefile}.fasta", 'fasta')
        newrecord = SeqRecord(Seq(reffasta.seq[startofgene:endofgene]), id=reffasta.id, description=reffasta.description)
        SeqIO.write(newrecord, f"{newfile}.fasta", "fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make new reference depending on whether the entire genome or only part is to be used for the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--greference", required=True, help="fasta output file with reference sequences")
    parser.add_argument("--reference", required=True, help="GenBank file with reference sequences")
    parser.add_argument("--output", required=True, help="GenBank new reference file")
    parser.add_argument("--gene", help="add gene name, otherwise entire genome will be used")
    args = parser.parse_args()

    new_reference(args.greference, args.reference, args.output, args.gene)

