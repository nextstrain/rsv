"""
This part of the workflow handles sorting downloaded sequences and metadata
into a and b by aligning them to reference sequences.

It produces output files as

    metadata = "data/{type}/metadata.tsv"

    sequences = "data/{type}/sequences.fasta"

"""

TIME = ['1', '2','3']

rule align:
    input:
        sequences = rules.transform.output.sequences,
        reference = "config/{type}_{time}_reference.fasta"
    output:
        alignment = "data/{type}/{time}_sequences.aligned.fasta"
    threads: 2
    shell:
        """
        nextalign run -j {threads} --silent \
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
        metadata = "data/{type}/{time}_metadata.tsv",
        sequences = "data/{type}/{time}_sequences.fasta"

    shell:
        """
        python bin/sequencesandmetadata.py \
        --sortedalignment {input.alignment} \
        --allmetadata {input.metadata} \
        --allsequences {input.sequences} \
        --metadata {output.metadata} \
        --sequences {output.sequences}
        """

rule sort:
    input:
        allsequences = "data/sequences.fasta",
        metadata = "data/metadata.tsv",
        alignment_a = expand("data/a/{time}_sequences.aligned.fasta", time=TIME),
        alignment_b = expand("data/b/{time}_sequences.aligned.fasta", time=TIME),
        reference_a = expand("config/a_{time}_reference.fasta", time=TIME),
        reference_b = expand("config/b_{time}_reference.fasta", time=TIME),
        metadata_b = expand("data/b/{time}_metadata.tsv", time=TIME),
        metadata_a = expand("data/a/{time}_metadata.tsv", time=TIME)
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
        dedup_metadata_a = "data/a/metadata_no_covg.tsv",
        dedup_seq_b = "data/b/sequences.fasta",
        dedup_metadata_b = "data/b/metadata_no_covg.tsv"
    shell:
        """
        seqkit rmdup < {input.sequences_a} > {output.dedup_seq_a}
        seqkit rmdup < {input.sequences_b} > {output.dedup_seq_b}

        python bin/metadata_dedup.py \
            --metadata-original {input.metadata_a} \
            --metadata-output {output.dedup_metadata_a}

        python bin/metadata_dedup.py \
            --metadata-original {input.metadata_b} \
            --metadata-output {output.dedup_metadata_b}
        """

rule coverage:
    input:
        alignment_a = expand("data/a/{time}_sequences.aligned.fasta", time=TIME),
        alignment_b = expand("data/b/{time}_sequences.aligned.fasta", time=TIME),
        metadata_b = expand("data/b/{time}_metadata.tsv", time=TIME),
        metadata_a = expand("data/a/{time}_metadata.tsv", time=TIME),
        dedup_metadata_a = rules.deduplication.output.dedup_metadata_a,
        dedup_metadata_b = rules.deduplication.output.dedup_metadata_b
    output:
        metadata_a = "data/a/metadata.tsv",
        metadata_b = "data/b/metadata.tsv"
    shell:
        """
        python bin/gene-coverage.py
        """