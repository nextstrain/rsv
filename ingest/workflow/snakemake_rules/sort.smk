"""
This part of the workflow handles sorting downloaded sequences and metadata
into a and b by aligning them to reference sequences.

It produces output files as

    metadata = "data/{type}/metadata.tsv"

    sequences = "data/{type}/sequences.fasta"

"""

TIME = ['1', '2','3']


rule sort:
    input:
        sequences = rules.transform.output.sequences
    output:
        "data/a/sequences.fasta",
        "data/b/sequences.fasta"
    shell:
        '''
        nextclade3 sort {input.sequences} --output-dir tmp
        seqkit rmdup tmp/nextstrain/rsv/b/sequences.fasta > data/b/sequences.fasta
        seqkit rmdup tmp/nextstrain/rsv/a/sequences.fasta > data/a/sequences.fasta
        rm -r tmp
        '''

rule metadata:
    input:
        metadata = rules.transform.output.metadata,
        sequences = "data/{type}/sequences.fasta"
    output:
        metadata = "data/{type}/metadata.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO

        strains = [s.id for s in SeqIO.parse(input.sequences, 'fasta')]
        d = pd.read_csv(input.metadata, sep='\t', index_col='accession').loc[strains].drop_duplicates()
        d.to_csv(output.metadata)
