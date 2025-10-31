"""
This part of the workflow handles transforming the data into standardized
formats and expects input file

    sequences_ndjson = "data/sequences.ndjson"

This will produce output files as

    metadata = "data/metadata.tsv"
    sequences = "data/sequences.fasta"

Parameters are expected to be defined in `config.transform`.
"""


rule curate:
    input:
        sequences_ndjson = "data/sequences.ndjson",
        geolocation_rules=config["curate"]["local_geolocation_rules"],
        annotations = config['curate']['annotations'],
    output:
        metadata = "data/curated_metadata.tsv",
        sequences = "data/sequences.fasta"
    log:
        "logs/curate.txt"
    params:
        field_map = config['curate']['ppx_field_map'],
        strain_regex = config['curate']['strain_regex'],
        strain_backup_fields = config['curate']['strain_backup_fields'],
        date_fields = config['curate']['date_fields'],
        expected_date_formats = config['curate']['expected_date_formats'],
        articles = config['curate']['titlecase']['articles'],
        abbreviations = config['curate']['titlecase']['abbreviations'],
        titlecase_fields = config['curate']['titlecase']['fields'],
        authors_field = config['curate']['authors_field'],
        authors_default_value = config['curate']['authors_default_value'],
        abbr_authors_field = config['curate']['abbr_authors_field'],
        ppx_division_field = config['curate']['ppx_division_field'],
        location_field = config['curate']['location_field'],
        annotations_id = config['curate']['annotations_id'],
        id_field = config['curate']['id_field'],
        sequence_field = config['curate']['sequence_field']
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
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields} \
                --articles {params.articles} \
                --abbreviations {params.abbreviations} \
            | augur curate abbreviate-authors \
                --authors-field {params.authors_field} \
                --default-value {params.authors_default_value} \
                --abbr-authors-field {params.abbr_authors_field} \
            | ./bin/parse-ppx-division \
                --division-field {params.ppx_division_field:q} \
                --location-field {params.location_field:q} \
            | augur curate apply-geolocation-rules \
                --geolocation-rules {input.geolocation_rules} \
            | python ./bin/curate-urls.py \
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
        metadata_fields=",".join(config["curate"]["ppx_metadata_columns"]),
    shell:
        """
        tsv-select -H -f {params.metadata_fields} \
            {input.metadata} > {output.subset_metadata}
        """
