typ = config['rsv']

rule download:
    message: "downloading sequences and metadata from data.nextstrain.org"
    output:
        metadata = "data/" +typ + "/metadata.tsv.gz",
        sequences = "data/"  + typ + "/sequences.fasta.xz"
    params:
        metadata_url = "http://data.nextstrain.org/files/workflows/rsv/" + typ + "/metadata.tsv.gz",
        sequence_url = "http://data.nextstrain.org/files/workflows/rsv/" + typ + "/sequences.fasta.xz"
    shell:
        """
        curl -fsSL --compressed {params.metadata_url:q} --output {output.metadata}
        curl -fsSL --compressed {params.sequence_url:q} --output {output.sequences}
        """

rule decompress:
    message: "decompressing sequences and metadata"
    input:
        sequences = "data/" + typ + "/sequences.fasta.xz",
        metadata = "data/" + typ + "/metadata.tsv.gz"
    output:
        sequences = "data/"+ typ +"/sequences.fasta",
        metadata = "data/"+ typ +"/metadata.tsv"
    shell:
        """
        gzip --decompress --keep {input.metadata}
        xz --decompress --keep {input.sequences}
        """
         
