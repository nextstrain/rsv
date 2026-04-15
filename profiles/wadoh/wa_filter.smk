rule filter_recent_wa:
    message:
        """
        filtering Washington sequences
        """
    input:
        sequences="results/{a_or_b}/sequences.fasta",
        reference="config/{a_or_b}reference.gbk",
        metadata="results/{a_or_b}/metadata.tsv",
        sequence_index=rules.index_sequences.output,
        exclude=config["exclude"],
    output:
        sequences=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/filtered_recent.fasta",
    params:
#        group_by=config["filter"]["group_by"],
        min_coverage=lambda w: f'{w.build_name}_coverage>{config["filter"]["min_coverage"].get(w.build_name, 10000)}',
        min_length=lambda w: config["filter"]["min_length"].get(w.build_name, 10000),
#        subsample_max_sequences=lambda w: config["filter"][
#            "subsample_max_sequences"
#        ].get(w.build_name, 1000),
        strain_id=config["strain_id_field"],
#        min_date=lambda w: config["filter"]["resolutions"][w.resolution]["min_date"],
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --exclude {input.exclude} \
            --exclude-where 'qc.overallStatus=bad' \
            --min-date 12Y \
            --min-length {params.min_length} \
            --output {output.sequences} \
            --subsample-max-sequences 10000 \
            --query '({params.min_coverage}) & missing_data<1000 & division=="Washington"'
        """


rule filter_background_wa:
    message:
        """
        filtering background sequences
        """
    input:
        sequences="results/{a_or_b}/sequences.fasta",
        reference="config/{a_or_b}reference.gbk",
        metadata="results/{a_or_b}/metadata.tsv",
        sequence_index=rules.index_sequences.output,
        include="config/include_{a_or_b}.txt",
        exclude=config["exclude"],
    output:
        sequences=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/filtered_background.fasta",
        metadata=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/filtered_background_metadata.tsv",
    params:
        group_by=config["filter"]["group_by"],
        min_coverage=lambda w: f'{w.build_name}_coverage>{config["filter"]["min_coverage"].get(w.build_name, 10000)}',
        min_length=lambda w: config["filter"]["min_length"].get(w.build_name, 10000),
#        subsample_max_sequences=lambda w: int(
#            config["filter"]["subsample_max_sequences"].get(w.build_name, 1000)
#        ),
        strain_id=config["strain_id_field"],
#        max_date=lambda w: config["filter"]["resolutions"][w.resolution]["min_date"],
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
            --include {input.include} \
            --exclude {input.exclude} \
            --exclude-where 'qc.overallStatus=bad' 'qc.overallStatus=mediocre' \
            --min-date {params.min_date} \
            --min-length {params.min_length} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --group-by {params.group_by} \
            --subsample-max-sequences 3000 \
            --query '({params.min_coverage}) & missing_data<1000 & division!="Washington"'
        """

rule combine_samples_wa:
    input:
        subsamples=lambda w: (
            (
                [
                    rules.filter_recent_wa.output.sequences,
                    rules.filter_background_wa.output.sequences,
                ]
                if "background_min_date" in config["filter"]["resolutions"][w.resolution]
                else [rules.filter_recent_wa.output.sequences]
            )
            # potentially add sequences sampled to include maximum escape sequences
            + (
                [
                    f"{build_dir}/{w.a_or_b}/{w.build_name}/{w.resolution}/filtered_{antibody}_{scoretype}.fasta"
                    for antibody in config["f_dms_antibodies"]
                    for scoretype in ["total_escape", "max_escape"]
                ]
                if w.build_name in config["enrich_antibody_escape"]
                else []
            )
        ),
    output:
        sequences=build_dir + "/{a_or_b}/{build_name}/{resolution}/filtered.fasta",
    log:
        "logs/combine_samples_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/combine_samples_{a_or_b}_{build_name}_{resolution}.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.subsamples} | seqkit rmdup > {output}
        """

ruleorder: filter_recent_wa > filter_recent
ruleorder: filter_background_wa > filter_background
ruleorder: combine_samples_wa > combine_samples
