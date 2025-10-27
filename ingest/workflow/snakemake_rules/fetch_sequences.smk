"""
This part of the workflow handles fetching sequences from various sources.
Uses `config.sources` to determine which sequences to include in final output.

Currently only fetches sequences from GenBank, but other sources can be
defined in the config. If adding other sources, add a new rule upstream
of rule `fetch_all_sequences` to create the file `data/{source}.ndjson` or the
file must exist as a static file in the repo.

Produces final output as

    sequences_ndjson = "data/sequences.ndjson"

"""


rule download_ppx_seqs:
    output:
        sequences= "data/{subtype}/ppx_sequences.fasta",
    params:
        sequences_url=lambda w: config["ppx_fetch"][w.subtype]["seqs"],
    # Allow retries in case of network errors
    retries: 5
    shell:
        """
        curl -fsSL {params.sequences_url:q} -o {output.sequences}
        """

rule download_ppx_meta:
    output:
        metadata= "data/{subtype}/ppx_metadata.csv"
    params:
        metadata_url=lambda w: config["ppx_fetch"][w.subtype]["meta"],
        fields = ",".join(config["ppx_metadata_fields"])
    # Allow retries in case of network errors
    retries: 5
    shell:
        """
        curl -fsSL '{params.metadata_url}&fields={params.fields}' -o {output.metadata}
        """

rule format_ppx_ndjson:
    input:
        sequences="data/{subtype}/ppx_sequences.fasta",
        metadata="data/{subtype}/ppx_metadata.csv"
    output:
        ndjson="data/{subtype}/ppx.ndjson"
    log:
        "logs/{subtype}/format_ppx_ndjson.txt"
    shell:
        """
        augur curate passthru \
            --metadata {input.metadata} \
            --fasta {input.sequences} \
            --seq-id-column accessionVersion \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            2> logs/{wildcards.subtype}/format_ppx_ndjson.txt > {output.ndjson}
        """


rule fetch_all_ppx_sequences:
    input:
        "data/a/ppx.ndjson",
        "data/b/ppx.ndjson",
    output:
        sequences_ndjson="data/sequences.ndjson",
    shell:
        """
        cat {input} > {output.sequences_ndjson}
        """
