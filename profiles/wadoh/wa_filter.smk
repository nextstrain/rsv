rule filter_wa:
    message:
        """
        filtering Washington sequences
        """
    input:
        sequences="data/{a_or_b}/sequences.fasta",
        reference="config/{a_or_b}reference.gbk",
        metadata="data/{a_or_b}/metadata.tsv",
        sequence_index=rules.index_sequences.output,
        exclude=config["exclude"],
    output:
        sequences=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/filtered_wa.fasta",
    params:
        #group_by=config["filter"]["group_by"],
        min_coverage=lambda w: f'{w.build_name}_coverage>{config["filter"]["min_coverage"].get(w.build_name, 10000)}',
        min_length=lambda w: config["filter"]["min_length"].get(w.build_name, 10000),
        subsample_max_sequences=lambda w: config["filter"][
            "subsample_max_sequences"
        ].get(w.build_name, 1000),
        strain_id=config["strain_id_field"],
        min_date=lambda w: config["filter"]["resolutions"][w.resolution]["min_date"],
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --exclude {input.exclude} \
            --exclude-where 'qc.overallStatus=bad' \
            --min-date {params.min_date} \
            --min-length {params.min_length} \
            --output {output.sequences} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --query "({params.min_coverage}) & missing_data<1000 & division == 'Washington' "
        """


rule filter_contextual:
    message:
        """
        filtering sequences
        """
    input:
        sequences="data/{a_or_b}/sequences.fasta",
        reference="config/{a_or_b}reference.gbk",
        metadata="data/{a_or_b}/metadata.tsv",
        sequence_index=rules.index_sequences.output,
        exclude=config["exclude"],
    output:
        sequences=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/filtered_background.fasta",
    params:
        group_by=config["filter"]["group_by"],
        min_coverage=lambda w: f'{w.build_name}_coverage>{config["filter"]["min_coverage"].get(w.build_name, 10000)}',
        min_length=lambda w: config["filter"]["min_length"].get(w.build_name, 10000),
        subsample_max_sequences=lambda w: config["filter"]["subsample_max_sequences"].get(w.build_name, 3000),
        strain_id=config["strain_id_field"],
        max_date=lambda w: config["filter"]["resolutions"][w.resolution]["min_date"],
        min_date=lambda w: config["filter"]["resolutions"][w.resolution][
            "background_min_date"
        ],
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --exclude {input.exclude} \
            --exclude-where 'qc.overallStatus=bad' 'qc.overallStatus=mediocre'\
            --min-date {params.min_date} \
            --min-length {params.min_length} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --query "({params.min_coverage}) & missing_data<1000 & division != 'Washington' "
        """


rule combine_all_samples:
    input:
        subsamples=lambda w: (
            [
                rules.filter_wa.output.sequences,
                rules.filter_contextual.output.sequences,
            ]
            if "background_min_date" in config["filter"]["resolutions"][w.resolution]
            else [rules.filter_wa.output.sequences]
        ),
    output:
        sequences=build_dir + "/{a_or_b}/{build_name}/{resolution}/filtered.fasta",
    shell:
        """
        cat {input.subsamples} | seqkit rmdup > {output}
        """

ruleorder: filter_wa > filter_recent
ruleorder: filter_contextual > filter_background
ruleorder: combine_all_samples > combine_samples
