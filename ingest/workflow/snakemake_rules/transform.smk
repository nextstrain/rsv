"""
This part of the workflow handles transforming the data into standardized
formats and expects input file

    sequences_ndjson = "data/sequences.ndjson"

This will produce output files as

    metadata = "data/metadata.tsv"
    sequences = "data/sequences.fasta"

Parameters are expected to be defined in `config.transform`.
"""

rule fetch_general_geolocation_rules:
    output:
        general_geolocation_rules = "data/general-geolocation-rules.tsv"
    params:
        geolocation_rules_url = config['transform']['geolocation_rules_url']
    shell:
        """
        curl {params.geolocation_rules_url} > {output.general_geolocation_rules}
        """

rule concat_geolocation_rules:
    input:
        general_geolocation_rules = "data/general-geolocation-rules.tsv",
        local_geolocation_rules = config['transform']['local_geolocation_rules']
    output:
        all_geolocation_rules = "data/all-geolocation-rules.tsv"
    shell:
        """
        cat {input.general_geolocation_rules} {input.local_geolocation_rules} >> {output.all_geolocation_rules}
        """




rule transform:
    input:
        sequences_ndjson = "data/sequences.ndjson",
        all_geolocation_rules = "data/all-geolocation-rules.tsv",
        annotations = config['transform']['annotations'],
    output:
        metadata = "data/curated_metadata.tsv",
        sequences = "data/sequences.fasta"
    log:
        "logs/transform.txt"
    params:
        field_map = config['transform']['field_map'],
        strain_regex = config['transform']['strain_regex'],
        strain_backup_fields = config['transform']['strain_backup_fields'],
        date_fields = config['transform']['date_fields'],
        expected_date_formats = config['transform']['expected_date_formats'],
        genbank_location_field=config["transform"]["genbank_location_field"],
        articles = config['transform']['titlecase']['articles'],
        abbreviations = config['transform']['titlecase']['abbreviations'],
        titlecase_fields = config['transform']['titlecase']['fields'],
        authors_field = config['transform']['authors_field'],
        authors_default_value = config['transform']['authors_default_value'],
        abbr_authors_field = config['transform']['abbr_authors_field'],
        annotations_id = config['transform']['annotations_id'],
        id_field = config['transform']['id_field'],
        sequence_field = config['transform']['sequence_field']
    shell:
        """
        (cat {input.sequences_ndjson} \
            | augur curate rename \
                --field-map {params.field_map} \
            | augur curate normalize-strings \
            | augur curate transform-strain-name \
                --strain-regex {params.strain_regex} \
                --backup-fields {params.strain_backup_fields} \
            | augur curate format-dates \
                --date-fields {params.date_fields} \
                --expected-date-formats {params.expected_date_formats} \
            | augur curate parse-genbank-location \
                --location-field {params.genbank_location_field} \
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields} \
                --articles {params.articles} \
                --abbreviations {params.abbreviations} \
            | augur curate abbreviate-authors \
                --authors-field {params.authors_field} \
                --default-value {params.authors_default_value} \
                --abbr-authors-field {params.abbr_authors_field} \
            | augur curate apply-geolocation-rules \
                --geolocation-rules {input.all_geolocation_rules} \
            | augur curate apply-record-annotations \
                --annotations {input.annotations} \
                --id-field {params.annotations_id} \
                --output-fasta {output.sequences} \
                --output-metadata {output.metadata} \
                --output-id-field {params.id_field} \
                --output-seq-field {params.sequence_field} ) 2>> {log}
        """

rule subset_metadata:
    input:
        metadata = "data/curated_metadata.tsv",
    output:
        subset_metadata="data/metadata.tsv",
    params:
        metadata_fields=",".join(config["transform"]["metadata_columns"]),
    shell:
        """
        tsv-select -H -f {params.metadata_fields} \
            {input.metadata} > {output.subset_metadata}
        """