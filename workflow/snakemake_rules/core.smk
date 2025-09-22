"""
This part of the workflow expects input files
            sequences = "data/sequences.fasta"
            metadata = "data/metadata.tsv"
"""


rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences="data/{a_or_b}/sequences.fasta",
    output:
        sequence_index=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/sequence_index.tsv",
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """


rule newreference:
    message:
        """
        Making new reference
        """
    input:
        oldreference="config/{a_or_b}reference.gbk",
    output:
        newreferencegbk=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/{gene}_reference.gbk",
        newreferencefasta=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/{gene}_reference.fasta",
    params:
        gene=lambda w: w.gene,
    shell:
        """
        python scripts/newreference.py \
            --reference {input.oldreference} \
            --output-genbank {output.newreferencegbk} \
            --output-fasta {output.newreferencefasta} \
            --gene {params.gene}
        """


rule filter_recent:
    message:
        """
        filtering sequences
        """
    input:
        sequences="data/{a_or_b}/sequences.fasta",
        metadata="data/{a_or_b}/metadata.tsv",
        sequence_index=rules.index_sequences.output,
        exclude=config["exclude"],
    output:
        sequences=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/filtered_recent.fasta",
    params:
        group_by=config["filter"]["group_by"],
        min_coverage=lambda w: f'{w.build_name}_coverage>{config["filter"]["min_coverage"].get(w.build_name, 10000)}',
        min_length=lambda w: config["filter"]["min_length"].get(w.build_name, 10000),
        subsample_max_sequences=lambda w: config["filter"][
            "subsample_max_sequences"
        ].get(w.build_name, 1000),
        strain_id=config["strain_id_field"],
        min_date=lambda w: config["filter"]["resolutions"][w.resolution]["min_date"],
        exclude_where=lambda w: " ".join([f"'{item}'" for item in config["filter"]["exclude_where"]["recent"]]),
        missing_data_threshold=config["filter"]["missing_data_threshold"],
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where} \
            --min-date {params.min_date} \
            --min-length {params.min_length} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --query '({params.min_coverage}) & missing_data<{params.missing_data_threshold}'
        """


rule filter_background:
    message:
        """
        filtering sequences
        """
    input:
        sequences="data/{a_or_b}/sequences.fasta",
        metadata="data/{a_or_b}/metadata.tsv",
        sequence_index=rules.index_sequences.output,
        include="config/include_{a_or_b}.txt",
        exclude=config["exclude"],
    output:
        sequences=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/filtered_background.fasta",
    params:
        group_by=config["filter"]["group_by"],
        min_coverage=lambda w: f'{w.build_name}_coverage>{config["filter"]["min_coverage"].get(w.build_name, 10000)}',
        min_length=lambda w: config["filter"]["min_length"].get(w.build_name, 10000),
        subsample_max_sequences=lambda w: int(
            config["filter"]["subsample_max_sequences"].get(w.build_name, 1000)
        )
        // 10,
        strain_id=config["strain_id_field"],
        max_date=lambda w: config["filter"]["resolutions"][w.resolution]["min_date"],
        min_date=lambda w: config["filter"]["resolutions"][w.resolution][
            "background_min_date"
        ],
        exclude_where=lambda w: " ".join([f"'{item}'" for item in config["filter"]["exclude_where"]["background"]]),
        missing_data_threshold=config["filter"]["missing_data_threshold"],
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --include {input.include} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --min-length {params.min_length} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --query '({params.min_coverage}) & missing_data<{params.missing_data_threshold}'
        """


rule combine_samples:
    input:
        subsamples=lambda w: (
            [
                rules.filter_recent.output.sequences,
                rules.filter_background.output.sequences,
            ]
            if "background_min_date" in config["filter"]["resolutions"][w.resolution]
            else [rules.filter_recent.output.sequences]
        ),
    output:
        sequences=build_dir + "/{a_or_b}/{build_name}/{resolution}/filtered.fasta",
    shell:
        """
        cat {input.subsamples} | seqkit rmdup > {output}
        """


rule get_nextclade_dataset:
    message:
        """
        fetching nextclade dataset
        """
    output:
        dataset="results/nextclade_rsv-{a_or_b}.zip",
    params:
        ds_name=lambda w: (
            "nextstrain/rsv/a/EPI_ISL_412866"
            if w.a_or_b == "a"
            else "nextstrain/rsv/b/EPI_ISL_1653999"
        ),
    shell:
        """
        nextclade3 dataset get -n {params.ds_name} --output-zip {output.dataset}
        """


rule genome_align:
    message:
        """
        Aligning sequences to the reference
        """
    input:
        sequences=rules.combine_samples.output.sequences,
        dataset=rules.get_nextclade_dataset.output.dataset,
    output:
        alignment=build_dir + "/{a_or_b}/{build_name}/{resolution}/sequences.aligned.fasta",
        translations=directory(build_dir + "/{a_or_b}/{build_name}/{resolution}/translations"),
    params:
        genes=lambda w: config["cds"][w.build_name],
    threads: 4
    log:
        "logs/align_{a_or_b}_{build_name}_{resolution}.txt",
    shell:
        """
        nextclade3 run -j {threads}\
            {input.sequences} \
            -D {input.dataset} \
            --output-fasta {output.alignment} \
            --cds-selection {params.genes} \
            --output-translations "{output.translations}/{{cds}}.fasta" 2>&1 | tee {log} \
        """


# cut out the G-Gene for alignment refinement
rule cut:
    input:
        oldalignment=rules.genome_align.output.alignment,
        reference="config/{a_or_b}reference.gbk",
    output:
        slicedalignment=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/{gene}_slicedalignment.fasta",
    params:
        gene=lambda w: w.gene,
    shell:
        """
        python scripts/cut.py \
            --oldalignment {input.oldalignment} \
            --slicedalignment {output.slicedalignment} \
            --reference {input.reference} \
            --gene {params.gene}
        """


# align the G gene with mafft
rule realign:
    input:
        slicedalignment=rules.cut.output.slicedalignment,
        reference=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/{gene}_reference.fasta",
    output:
        realigned=build_dir + "/{a_or_b}/{build_name}/{resolution}/{gene}_aligned.fasta",
    threads: 4
    shell:
        """
        augur align --nthreads {threads} \
            --sequences {input.slicedalignment} \
            --reference-sequence {input.reference} \
            --output {output.realigned}
        """


rule hybrid_align:
    input:
        original=rules.genome_align.output.alignment,
        G_alignment=build_dir + "/{a_or_b}/{build_name}/{resolution}/G_aligned.fasta",
        reference="config/{a_or_b}reference.gbk",
    output:
        hybrid_alignment=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/hybrid_alignment.fasta",
    params:
        gene=lambda w: w.build_name,
    shell:
        """
        python scripts/align_for_tree.py \
            --realign {input.G_alignment} \
            --original {input.original} \
            --reference {input.reference} \
            --output {output.hybrid_alignment} \
            --gene {params.gene}
        """


def get_alignment(w):
    if w.build_name == "genome":
        return rules.hybrid_align.output.hybrid_alignment
    else:
        return (
            build_dir
            + f"/{w.a_or_b}/{w.build_name}/{w.resolution}/{w.build_name}_aligned.fasta"
        )


rule tree:
    message:
        "Building tree"
    input:
        alignment=get_alignment,
    output:
        tree=build_dir + "/{a_or_b}/{build_name}/{resolution}/tree_raw.nwk",
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --tree-builder-args '-ninit 10 -n 4 -czb' \
            --nthreads {threads}
        """


rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree=rules.tree.output.tree,
        alignment=get_alignment,
        metadata="data/{a_or_b}/metadata.tsv",
    output:
        tree=build_dir + "/{a_or_b}/{build_name}/{resolution}/tree.nwk",
        node_data=build_dir + "/{a_or_b}/{build_name}/{resolution}/branch_lengths.json",
    params:
        coalescent=config["refine"]["coalescent"],
        clock_filter_iqd=config["refine"]["clock_filter_iqd"],
        date_inference=config["refine"]["date_inference"],
        strain_id=config["strain_id_field"],
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --timetree \
            --stochastic-resolve \
            --use-fft \
            --clock-filter-iqd {params.clock_filter_iqd}
        """


def _get_build_distance_map_config(wildcards):
    distance_config = distance_map_config[
        (distance_map_config["a_or_b"] == wildcards.a_or_b)
        & (distance_map_config["build_name"] == wildcards.build_name)
        & (distance_map_config["resolution"] == wildcards.resolution)
    ]
    if distance_config.shape[0] > 0:
        return distance_config
    else:
        return None

def _get_distance_comparisons_by_lineage_and_segment(wildcards):
    config = _get_build_distance_map_config(wildcards)
    return " ".join(config.loc[:, "compare_to"].values)


def _get_distance_attributes_by_lineage_and_segment(wildcards):
    config = _get_build_distance_map_config(wildcards)
    return " ".join(config.loc[:, "attribute"].values)


def _get_distance_maps_by_lineage_and_segment(wildcards):
    distance_config = _get_build_distance_map_config(wildcards)
    if wildcards.build_name != "G":
        return [
            "config/distance_maps/{wildcards.build_name}/{distance_map}.json".format(wildcards=wildcards, distance_map=distance_map)
            for distance_map in distance_config.loc[:, "distance_map"].values
        ]
    else:
        return ""


rule distances:
    input:
        tree=rules.refine.output.tree,
        distance_maps = _get_distance_maps_by_lineage_and_segment,
        translations_done=build_dir + "/{a_or_b}/{build_name}/{resolution}/translations.done",
    output:
        distances= build_dir + "/{a_or_b}/{build_name}/{resolution}/distances.json"
    params:
        genes=lambda w: config["cds"][w.build_name],
        alignments=lambda w: [f"{build_dir}/{w.a_or_b}/{w.build_name}/{w.resolution}/translations/{gene}_withInternalNodes.fasta" for gene in config["cds"][w.build_name]],
        comparisons=_get_distance_comparisons_by_lineage_and_segment,
        attribute_names=_get_distance_attributes_by_lineage_and_segment,
    resources:
        mem_mb=8000,
        time="00:30:00",
    shell:
        """
        augur distance \
            --alignment {params.alignments} \
            --tree {input.tree} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --output {output.distances} 2>&1 | tee {log}
        """

rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - inferring ambiguous mutations
        """
    input:
        tree=rules.refine.output.tree,
        alignment=get_alignment,
        translations=rules.genome_align.output.translations,
        root_sequence=build_dir + "/{a_or_b}/{build_name}/{resolution}/{build_name}_reference.gbk",
    output:
        node_data=build_dir + "/{a_or_b}/{build_name}/{resolution}/nt_muts.json",
        translations_done=build_dir + "/{a_or_b}/{build_name}/{resolution}/translations.done",
    params:
        inference=config["ancestral"]["inference"],
        genes=lambda w: config["cds"][w.build_name],
        output_translations=lambda w: build_dir + f"/{w.a_or_b}/{w.build_name}/{w.resolution}/translations/%GENE_withInternalNodes.fasta",
        input_translations=lambda w: build_dir + f"/{w.a_or_b}/{w.build_name}/{w.resolution}/translations/%GENE.fasta",
    log:
        "logs/ancestral_{a_or_b}_{build_name}_{resolution}.txt",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --annotation {input.root_sequence} \
            --root-sequence {input.root_sequence} \
            --genes {params.genes} \
            --translations "{params.input_translations}" \
            --output-translations "{params.output_translations}" \
            --inference {params.inference} 2>&1 | tee {log} && touch {output.translations_done}
        """


rule translate:
    message:
        "Translating amino acid sequences"
    input:
        tree=rules.refine.output.tree,
        node_data=rules.ancestral.output.node_data,
        reference=build_dir + "/{a_or_b}/{build_name}/{resolution}/{build_name}_reference.gbk",
    output:
        node_data=build_dir + "/{a_or_b}/{build_name}/{resolution}/aa_muts.json",
    params:
        alignment_file_mask=build_dir + "/{a_or_b}/{build_name}/{resolution}/aligned_%GENE.fasta",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
            --alignment-output {params.alignment_file_mask}
        """


rule traits:
    input:
        tree=rules.refine.output.tree,
        metadata="data/{a_or_b}/metadata.tsv",
    output:
        node_data=build_dir + "/{a_or_b}/{build_name}/{resolution}/traits.json",
    log:
        "logs/{a_or_b}/traits_{build_name}_{resolution}_rsv.txt",
    params:
        columns=config["traits"]["columns"],
        strain_id=config["strain_id_field"],
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule frequencies:
    input:
        tree = rules.refine.output.tree,
        metadata = "data/{a_or_b}/metadata.tsv"
    output:
        frequencies = build_dir + "/{a_or_b}/{build_name}/{resolution}/frequencies.json"
    params:
        min_date_arg = lambda w: f"--min-date {config['frequencies']['resolutions'][w.resolution]['min_date']}" if w.resolution in config["frequencies"].get('resolutions', {}) else "",
    shell:
        """
        augur frequencies \
            --tree {input.tree} \
            --method kde \
            --metadata-id-columns accession \
            {params.min_date_arg} \
            --metadata {input.metadata} \
            --output {output.frequencies}
        """