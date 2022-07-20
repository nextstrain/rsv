
rule download:
    message: "downloading sequences and metadata from data.nextstrain.org"
    output:
        metadata =  "data/{a_or_b}/metadata.tsv.gz",
        sequences = "data/{a_or_b}/sequences.fasta.xz"
    params:
        metadata_url = "http://data.nextstrain.org/files/workflows/rsv/{a_or_b}/metadata.tsv.gz",
        sequence_url = "http://data.nextstrain.org/files/workflows/rsv/{a_or_b}/sequences.fasta.xz"
    shell:
        """
        curl -fsSL --compressed {params.metadata_url:q} --output {output.metadata}
        curl -fsSL --compressed {params.sequence_url:q} --output {output.sequences}
        """

rule decompress:
    message: "decompressing sequences and metadata"
    input:
        sequences = "data/{a_or_b}/sequences.fasta.xz",
        metadata = "data/{a_or_b}/metadata.tsv.gz"
    output:
        sequences = "data/{a_or_b}/sequences.fasta",
        metadata = "data/{a_or_b}/metadata.tsv"
    shell:
        """
        gzip --decompress --keep {input.metadata}
        xz --decompress --keep {input.sequences}
        """

