from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse
import sys

def new_reference(referencefile, outgenbank, outfasta, gene):
    ref = SeqIO.read(referencefile, "genbank")
    startofgene = None
    endofgene = None
    for feature in ref.features:
        if feature.type == 'source':
            ref_source_feature = feature
        if feature.type =='gene' or feature.type == 'CDS':
            a = list(feature.qualifiers.items())[0][-1][0]
            if a == gene:
                startofgene = int(list(feature.location)[0])
                endofgene =  int(list(feature.location)[-1])+1

    # If user provides a --gene 'some name' that is not found, error out as this may indicate that
    # the gene name is misspelled or the user may be using the wrong GenBank file.
    if(startofgene is None and endofgene is None):
        print(f"ERROR: No '{gene}' was found under 'gene' or 'CDS' features in the GenBank file.", file=sys.stderr)
        sys.exit(1)

    record = ref[startofgene:endofgene]
    source_feature =  SeqFeature(FeatureLocation(start=0, end=len(record)), type='source',
                                 qualifiers=ref_source_feature.qualifiers)
    record.features.append(source_feature)

    SeqIO.write(record, outgenbank, 'genbank')
    SeqIO.write(record, outfasta, "fasta")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make new reference based on a gene",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--reference", required=True, help="GenBank file with reference sequences")
    parser.add_argument("--output-fasta", required=True, help="FASTA new reference file")
    parser.add_argument("--output-genbank", required=True, help="GenBank new reference file")
    parser.add_argument("--gene", required=True, help="gene name")
    args = parser.parse_args()

    new_reference(args.reference, args.output_genbank, args.output_fasta, args.gene)

