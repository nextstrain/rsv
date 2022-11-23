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

rule fetch_from_genbank:
    output:
        csv = "data/genbank.csv"
    params:
        URL_a = config['fetch']['genbank_url']['a'],
        URL_b = config['fetch']['genbank_url']['b'],
        URL_general = config['fetch']['genbank_url']['general']
    conda: config["conda_environment"]
    shell:
        """
        curl "{params.URL_a}" --fail --silent --show-error --http1.1 \
             --header 'User-Agent: https://github.com/nextstrain/rsv (hello@nextstrain.org)' >> {output}
        curl "{params.URL_b}" --fail --silent --show-error --http1.1 \
             --header 'User-Agent: https://github.com/nextstrain/rsv (hello@nextstrain.org)' >> {output}
        curl "{params.URL_general}" --fail --silent --show-error --http1.1 \
             --header 'User-Agent: https://github.com/nextstrain/rsv (hello@nextstrain.org)' >> {output}
        """

rule csv_to_ndjson:
    input:
        csv = rules.fetch_from_genbank.output.csv
    output:
        ndjson = "data/genbank.ndjson"
    shell:
        """
        python bin/csv-to-ndjson.py \
            --input {input.csv} \
            --output {output.ndjson}
        """


rule fetch_all_sequences:
    input:
        all_sources = "data/genbank.ndjson"
    output:
        sequences_ndjson = "data/sequences.ndjson"
    shell:
        """
        cat {input.all_sources} > {output.sequences_ndjson}
        """
