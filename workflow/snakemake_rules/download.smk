rule download:
    message: "downloading sequences and metadata from data.nextstrain.org"
    output:
        metadata = "data/{type}_metadata.tsv.gz",
        sequences = "data/{type}_sequences.fasta.xz"
    params:
        metadata_url = "http://data.nextstrain.org/files/workflows/rsv/metadata.tsv.gz",
        sequence_url = "http://data.nextstrain.org/files/workflows/rsv/sequences.fasta.xz"
    shell:
        """
        curl -fsSL --compressed {params.metadata_url:q} --output {output.metadata}
        curl -fsSL --compressed {params.sequence_url:q} --output {output.sequences}
        """

rule decompress:
    message: "decompressing sequences and metadata"
    input:
        sequences = "data/{type}_sequences.fasta.xz",
        metadata = "data/{type}_metadata.tsv.gz"
    output:
        sequences = "data/sequences.fasta",
        metadata = "data/metadata.tsv"
    shell:
        """
        gzip --decompress --keep {input.metadata}
        xz --decompress --keep {input.sequences}
        """
         
