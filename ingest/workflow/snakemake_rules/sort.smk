"""
This part of the workflow handles sorting downloaded sequences and metadata
into a and b by aligning them to reference sequences.

It produces output files as

    metadata = "data/{type}/metadata.tsv"

    sequences = "data/{type}/sequences.fasta"

"""

rule sort:
    input:
        sequences = rules.curate.output.sequences
    output:
        "data/a/sequences.fasta",
        "data/b/sequences.fasta"
    shell:
        '''
        seqkit rmdup {input.sequences} | \
        nextclade3 sort - --output-dir tmp
        mv tmp/nextstrain/rsv/b/sequences.fasta  data/b/sequences.fasta
        mv tmp/nextstrain/rsv/a/sequences.fasta  data/a/sequences.fasta
        rm -r tmp
        '''

rule metadata:
    input:
        metadata = "data/metadata.tsv",
        sequences = "data/{type}/sequences.fasta"
    output:
        metadata = "data/{type}/metadata_raw.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO

        strains = [s.id for s in SeqIO.parse(input.sequences, 'fasta')]
        d = pd.read_csv(input.metadata, sep='\t', index_col='accession').loc[strains].drop_duplicates()
        d.to_csv(output.metadata, sep='\t')

rule nextclade_dataset:
    output:
        ref_a = "rsv-a_nextclade/reference.fasta",
        ref_b = "rsv-b_nextclade/reference.fasta"
    params:
        dataset_a = "nextstrain/rsv/a/EPI_ISL_412866",
        dataset_b = "nextstrain/rsv/b/EPI_ISL_1653999"
    shell:
        """
        nextclade3 dataset get -n {params.dataset_a} --output-dir rsv-a_nextclade
        nextclade3 dataset get -n {params.dataset_b} --output-dir rsv-b_nextclade
        """

rule nextclade:
    input:
        sequences = "data/{type}/sequences.fasta",
        ref = "rsv-a_nextclade/reference.fasta"
    output:
        nextclade = "data/{type}/nextclade.tsv"
    params:
        dataset = "rsv-{type}_nextclade",
        output_columns = "seqName clade qc.overallScore qc.overallStatus totalMissing alignmentScore  alignmentStart  alignmentEnd  coverage cdsCoverage dynamic"
    threads: 8
    shell:
        """
        nextclade3 run -D {params.dataset}  -j {threads} \
                          --output-columns-selection {params.output_columns} \
                          --output-tsv {output.nextclade} \
                          {input.sequences}
        """

rule extend_metadata:
    input:
        nextclade = "data/{type}/nextclade.tsv",
        metadata = "data/{type}/metadata_raw.tsv"
    output:
        metadata = "data/{type}/metadata.tsv"
    shell:
        """
        python3 bin/extend-metadata.py --metadata {input.metadata} \
                                       --id-field accession \
                                       --virus-type {wildcards.type} \
                                       --nextclade {input.nextclade} \
                                       --output {output.metadata}
        """


rule extract_open_data:
    input:
        metadata = "data/{type}/metadata.tsv",
        sequences = "data/{type}/sequences.fasta"
    output:
        metadata = "data/{type}/metadata_open.tsv",
        sequences = "data/{type}/sequences_open.fasta"
    shell:
        """
        augur filter --metadata {input.metadata} \
                     --sequences {input.sequences} \
                     --metadata-id-columns accession \
                     --include-where "dataUseTerms=OPEN" \
                     --output-metadata {output.metadata} \
                     --output-sequences {output.sequences}
        """

rule extract_restricted_data:
    input:
        metadata = "data/{type}/metadata.tsv",
        sequences = "data/{type}/sequences.fasta"
    output:
        metadata = "data/{type}/metadata_restricted.tsv",
        sequences = "data/{type}/sequences_restricted.fasta"
    shell:
        """
        augur filter --metadata {input.metadata} \
                     --sequences {input.sequences} \
                     --metadata-id-columns accession \
                     --include-where "dataUseTerms=RESTRICTED" \
                     --output-metadata {output.metadata} \
                     --output-sequences {output.sequences}
        """
