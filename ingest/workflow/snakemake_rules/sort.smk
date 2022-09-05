"""
This part of the workflow handles sorting downloaded sequences and metadata
into a and b by aligning them to reference sequences. 

It produces output files as

    metadata = "data/{type}/metadata.tsv"

    sequences = "data/{type}/sequences.fasta"

"""

TIME = ['new', 'middle','old']

rule align:
    input: 
        sequences = rules.transform.output.sequences,
        reference = "config/{type}{time}reference.fasta"
    output:
        alignment = "data/{type}/{time}sequences.aligned.fasta"
    shell:
        """
        nextalign run \
        --reference {input.reference} \
        --output-fasta {output.alignment} \
        {input.sequences}
        """

rule metadataandsequences:
    input:
        alignment = rules.align.output.alignment,
        metadata = rules.transform.output.metadata,
        sequences = rules.transform.output.sequences
    output:
        metadata = "data/{type}/{time}metadata.tsv",
        sequences = "data/{type}/{time}sequences.fasta"

    shell:
        """
        python bin/metadataandsequences.py \
        --sortedalignment {input.alignment} \
        --allmetadata {input.metadata} \
        --allsequences {input.sequences} \
        --metadata {output.metadata} \
        --sequences {output.sequences}
        """

rule sort:
    input:
        alignment_a = expand("data/a/{time}sequences.aligned.fasta", time=TIME),
        alignment_b = expand("data/b/{time}sequences.aligned.fasta", time=TIME),
        reference_a = expand("config/a{time}reference.fasta", time=TIME),
        reference_b = expand("config/b{time}reference.fasta", time=TIME),
        metadata_b = expand("data/b/{time}metadata.tsv", time=TIME),
        metadata_a = expand("data/a/{time}metadata.tsv", time=TIME)
    output:
        sequences_a = "data/a/sequences_notdedup.fasta",
        metadata_a = "data/a/metadata_notdedup.tsv",
        sequences_b = "data/b/sequences_notdedup.fasta",
        metadata_b = "data/b/metadata_notdedup.tsv"
    shell:
        """
        python bin/sort.py 
        """

rule deduplication:
    input:
        sequences_a = rules.sort.output.sequences_a,
        metadata_a = rules.sort.output.metadata_a,
        sequences_b = rules.sort.output.sequences_b,
        metadata_b = rules.sort.output.metadata_b     
    output:
        dedup_seq_a = "data/a/sequences.fasta",
        dedup_metadata_a = "data/a/metadata.tsv",
        dedup_seq_b = "data/b/sequences.fasta",
        dedup_metadata_b = "data/b/metadata.tsv"
    shell:
        """
        seqkit rmdup < {input.sequences_a} > {output.dedup_seq_a}
        seqkit rmdup < {input.sequences_b} > {output.dedup_seq_b}

        python bin/metadatadedup.py \
            --metadataoriginal {input.metadata_a} \
            --metadataoutput {output.dedup_metadata_a}

        python bin/metadatadedup.py \
            --metadataoriginal {input.metadata_b} \
            --metadataoutput {output.dedup_metadata_b}
        """