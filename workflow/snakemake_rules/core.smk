"""
This part of the workflow expects input files
            sequences = "data/sequences.fasta"
            metadata = "data/metadata.tsv"
"""



rule index_sequences:
    """
    Creating an index of sequence composition for filtering.
    """
    input:
        sequences="results/{a_or_b}/sequences.fasta",
    output:
        sequence_index=results_dir
        + "/{a_or_b}/sequence_index.tsv",
    log:
        "logs/index_sequences_{a_or_b}.txt"
    benchmark:
        "benchmarks/index_sequences_{a_or_b}.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """


rule newreference:
    """
    Making new reference
    """
    input:
        oldreference="config/{a_or_b}reference.gbk",
    output:
        newreferencegbk=results_dir
        + "/{a_or_b}/{gene}_reference.gbk",
        newreferencefasta=results_dir
        + "/{a_or_b}/{gene}_reference.fasta",
    log:
        "logs/newreference_{a_or_b}_{gene}.txt"
    benchmark:
        "benchmarks/newreference_{a_or_b}_{gene}.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/newreference.py \
            --reference {input.oldreference} \
            --output-genbank {output.newreferencegbk} \
            --output-fasta {output.newreferencefasta} \
            --gene {wildcards.gene}
        """


rule filter_recent:
    """
    filtering sequences
    """
    input:
        sequences=lambda w: f"results/{get_subtype(w.build)}/sequences.fasta",
        metadata=lambda w: f"results/{get_subtype(w.build)}/metadata.tsv",
        sequence_index=lambda w: f"results/{get_subtype(w.build)}/sequence_index.tsv",
        exclude=config["exclude"],
    output:
        sequences=results_dir + "/{build}/filtered_recent.fasta",
    log:
        "logs/{build}/filter_recent.txt"
    benchmark:
        "benchmarks/{build}/filter_recent.txt"
    params:
        group_by=config["filter"]["group_by"],
        min_coverage=lambda w: f'{get_gene(w.build)}_coverage>{config["filter"]["min_coverage"][get_gene_build(w.build)]}',
        min_length=lambda w: config["filter"]["min_length"][get_gene_build(w.build)],
        subsample_max_sequences=lambda w: config["filter"][
            "subsample_max_sequences"
        ][get_gene_build(w.build)],
        strain_id=config["strain_id_field"],
        min_date=lambda w: config["filter"]["resolutions"][get_resolution(w.build)]["min_date"],
        exclude_where=config["filter"]["exclude_where"]["recent"],
        missing_data_threshold=config["filter"]["missing_data_threshold"],
    shell:
        r"""
        exec &> >(tee {log:q})

        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where:q} \
            --min-date {params.min_date} \
            --min-length {params.min_length} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --query '({params.min_coverage}) & missing_data<{params.missing_data_threshold}'
        """


rule filter_background:
    """
    filtering sequences
    """
    input:
        sequences=lambda w: f"results/{get_subtype(w.build)}/sequences.fasta",
        metadata=lambda w: f"results/{get_subtype(w.build)}/metadata.tsv",
        sequence_index=lambda w: f"results/{get_subtype(w.build)}/sequence_index.tsv",
        include=lambda w: f"config/include_{get_subtype(w.build)}.txt",
        exclude=config["exclude"],
    output:
        sequences=results_dir
        + "/{build}/filtered_background.fasta",
        metadata=results_dir
        + "/{build}/filtered_background_metadata.tsv",
    log:
        "logs/{build}/filter_background.txt"
    benchmark:
        "benchmarks/{build}/filter_background.txt"
    params:
        group_by=config["filter"]["group_by"],
        min_coverage=lambda w: f'{get_gene(w.build)}_coverage>{config["filter"]["min_coverage"][get_gene_build(w.build)]}',
        min_length=lambda w: config["filter"]["min_length"][get_gene_build(w.build)],
        subsample_max_sequences=lambda w: int(
            config["filter"]["subsample_max_sequences"][get_gene_build(w.build)],
        )
        // 10,
        strain_id=config["strain_id_field"],
        max_date=lambda w: config["filter"]["resolutions"][get_resolution(w.build)]["min_date"],
        min_date=lambda w: config["filter"]["resolutions"][get_resolution(w.build)][
            "background_min_date"
        ],
        exclude_where=config["filter"]["exclude_where"]["background"],
        missing_data_threshold=config["filter"]["missing_data_threshold"],
        clade_prefix=lambda w: f"{get_subtype(w.build).upper()}.D",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --include {input.include} \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where:q}  \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --min-length {params.min_length} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --query '({params.min_coverage}) & missing_data<{params.missing_data_threshold} & clade.str.startswith("{params.clade_prefix}", na=False)'
        """

rule combine_samples:
    input:
        subsamples=lambda w: (
            (
                [
                    rules.filter_recent.output.sequences,
                    rules.filter_background.output.sequences,
                ]
                if "background_min_date" in config["filter"]["resolutions"][get_resolution(w.build)]
                else [rules.filter_recent.output.sequences]
            )
            # potentially add sequences sampled to include maximum escape sequences
            + (
                [
                    f"{results_dir}/{w.build}/filtered_{antibody}_{scoretype}.fasta"
                    for antibody in config["f_dms_antibodies"]
                    for scoretype in ["total_escape", "max_escape"]
                ]
                if get_gene_build(w.build) in config["enrich_antibody_escape"]
                else []
            )
        ),
    output:
        sequences=results_dir + "/{build}/filtered.fasta",
    log:
        "logs/{build}/combine_samples.txt"
    benchmark:
        "benchmarks/{build}/combine_samples.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.subsamples} | seqkit rmdup > {output}
        """


rule get_nextclade_dataset:
    """
    fetching nextclade dataset
    """
    output:
        dataset="results/nextclade_rsv-{a_or_b}.zip",
    log:
        "logs/get_nextclade_dataset_{a_or_b}.txt"
    benchmark:
        "benchmarks/get_nextclade_dataset_{a_or_b}.txt"
    params:
        ds_name=lambda w: (
            "nextstrain/rsv/a/EPI_ISL_412866"
            if w.a_or_b == "a"
            else "nextstrain/rsv/b/EPI_ISL_1653999"
        ),
    shell:
        r"""
        exec &> >(tee {log:q})

        nextclade3 dataset get -n {params.ds_name} --output-zip {output.dataset}
        """


rule filter_for_pre_subsample_alignment:
    """
    Do the quality filtering applied to each sequence set before subsampling
    """
    input:
        sequences=lambda w: f"results/{get_subtype(w.build)}/sequences.fasta",
        metadata=lambda w: f"results/{get_subtype(w.build)}/metadata.tsv",
        exclude=config["exclude"],
    output:
        sequences=results_dir + "/{build}/pre_subsample/filtered_for_alignment.fasta",
    log:
        "logs/{build}/filter_for_pre_subsample_alignment.txt"
    benchmark:
        "benchmarks/{build}/filter_for_pre_subsample_alignment.txt"
    params:
        min_coverage=lambda w: f'{get_gene(w.build).split("-")[0]}_coverage>{config["filter"]["min_coverage"][get_gene(w.build)]}',
        min_length=lambda w: config["filter"]["min_length"][get_gene(w.build)],
        strain_id=config["strain_id_field"],
        min_date=lambda w: config["filter"]["resolutions"][get_resolution(w.build)]["min_date"],
    shell:
        r"""
        exec &> >(tee {log:q})

        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --exclude {input.exclude} \
            --exclude-where 'qc.overallStatus=bad' \
            --min-length {params.min_length} \
            --min-date {params.min_date} \
            --output {output.sequences} \
            --query '({params.min_coverage}) & missing_data<1000'
        """


rule align_pre_subsample_sequences:
    """
    Aligning all pre-subsampled quality-filtered sequences
    """
    input:
        sequences=rules.filter_for_pre_subsample_alignment.output.sequences,
        dataset=lambda w: f"results/nextclade_rsv-{get_subtype(w.build)}.zip",
    output:
        alignment=results_dir + "/{build}/pre_subsample/sequences.aligned.fasta",
        translations=directory(results_dir + "/{build}/pre_subsample/translations"),
        translations_done=results_dir + "/{build}/pre_subsample/translations.done",
    params:
        genes=lambda w: get_gene(w.build),
    threads: 8
    log:
        "logs/{build}/align_pre_subsample_sequences.txt"
    benchmark:
        "benchmarks/{build}/align_pre_subsample_sequences.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        nextclade3 run -j {threads}\
            {input.sequences} \
            -D {input.dataset} \
            --output-fasta {output.alignment} \
            --cds-selection {params.genes} \
            --output-translations "{output.translations}/{{cds}}.fasta" && touch {output.translations_done}
        """


rule score_pre_subsample_f_proteins:
    """
    Computing F protein DMS scores for pre-subsampled sequences
    """
    input:
        translations_done=rules.align_pre_subsample_sequences.output.translations_done,
        dms_scores=config["f_dms_data"],
    output:
        scores=results_dir + "/{build}/pre_subsample/f_protein_scores.tsv",
    log:
        "logs/{build}/score_pre_subsample_f_proteins.txt"
    benchmark:
        "benchmarks/{build}/score_pre_subsample_f_proteins.txt"
    params:
        f_sequences=lambda w: results_dir + f"/{w.build}/pre_subsample/translations/F.fasta",
        dms_antibodies=lambda w: " ".join(shlex.quote(ab) for ab in config["f_dms_antibodies"]),
        only_positive_escape=config["dms_only_positive_escape"],
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/score_f_sequences.py fasta \
            --sequences {params.f_sequences} \
            --dms-scores {input.dms_scores} \
            --output {output.scores} \
            --dms-antibodies {params.dms_antibodies} \
            --only-positive-escape {params.only_positive_escape}
        """


rule add_f_scores_to_pre_subsample_metadata:
    """
    Adding F protein scores to pre-subsampled metadata
    """
    input:
        original_metadata=lambda w: f"results/{get_subtype(w.build)}/metadata.tsv",
        f_scores=rules.score_pre_subsample_f_proteins.output.scores,
    output:
        enhanced_metadata=results_dir + "/{build}/pre_subsample/metadata_with_scores.tsv",
    log:
        "logs/{build}/add_f_scores_to_pre_subsample_metadata.txt"
    benchmark:
        "benchmarks/{build}/add_f_scores_to_pre_subsample_metadata.txt"
    params:
        strain_id=config["strain_id_field"],
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/merge_f_scores.py \
            --metadata {input.original_metadata} \
            --scores {input.f_scores} \
            --output {output.enhanced_metadata} \
            --strain-id-field {params.strain_id}
        """


rule enrich_antibody_escape:
    """
    Get sequences with high antibody escape to add to tree via custom filtering rule.
    """
    wildcard_constraints:
        antibody="|".join(re.escape(antibody) for antibody in config["f_dms_antibodies"]),
        scoretype="total_escape|max_escape",
    input:
        sequences=lambda w: f"results/{get_subtype(w.build)}/sequences.fasta",
        metadata=rules.add_f_scores_to_pre_subsample_metadata.output.enhanced_metadata,
    output:
        sequences=results_dir + "/{build}/filtered_{antibody}_{scoretype}.fasta"
    log:
        "logs/{build}/enrich_antibody_escape_{antibody}_{scoretype}.txt"
    benchmark:
        "benchmarks/{build}/enrich_antibody_escape_{antibody}_{scoretype}.txt"
    params:
        strain_id=config["strain_id_field"],
        escape_col=lambda w: f"{w.antibody}_{w.scoretype}",
        nseqs=lambda w: config["enrich_antibody_escape"][get_gene_build(w.build)]["nseqs_per_antibody_scoretype"],
        group_by=lambda w: " ".join(shlex.quote(g) for g in config["enrich_antibody_escape"][get_gene_build(w.build)]["group_by"]),
        max_identical_f_prot_muts=lambda w: config["enrich_antibody_escape"][get_gene_build(w.build)]["max_identical_f_prot_muts"],
        max_identical_max_escape_mut=lambda w: config["enrich_antibody_escape"][get_gene_build(w.build)]["max_identical_max_escape_mut"],
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/enrich_antibody_escape.py \
            --input-sequences {input.sequences} \
            --metadata {input.metadata} \
            --output-sequences {output.sequences} \
            --strain-id {params.strain_id} \
            --escape-col {params.escape_col} \
            --nseqs {params.nseqs} \
            --group-by {params.group_by} \
            --max-identical-f-prot-muts {params.max_identical_f_prot_muts} \
            --max-identical-max-escape-mut {params.max_identical_max_escape_mut}
        """


rule genome_align:
    """
    Aligning sequences to the reference
    """
    input:
        sequences=rules.combine_samples.output.sequences,
        dataset=lambda w: f"results/nextclade_rsv-{get_subtype(w.build)}.zip",
    output:
        alignment=results_dir + "/{build}/sequences.aligned.fasta",
        translations=directory(results_dir + "/{build}/translations"),
    params:
        genes=lambda w: get_gene(w.build),
    threads: 4
    log:
        "logs/{build}/genome_align.txt"
    benchmark:
        "benchmarks/{build}/genome_align.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        nextclade3 run -j {threads}\
            {input.sequences} \
            -D {input.dataset} \
            --output-fasta {output.alignment} \
            --cds-selection {params.genes} \
            --output-translations "{output.translations}/{{cds}}.fasta"
        """


# cut out the G-Gene for alignment refinement
rule cut:
    input:
        oldalignment=rules.genome_align.output.alignment,
        reference=lambda w: f"config/{get_subtype(w.build)}reference.gbk",
    output:
        slicedalignment=results_dir
        + "/{build}/{gene}_slicedalignment.fasta",
    log:
        "logs/{build}/cut_{gene}.txt"
    benchmark:
        "benchmarks/{build}/cut_{gene}.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/cut.py \
            --oldalignment {input.oldalignment} \
            --slicedalignment {output.slicedalignment} \
            --reference {input.reference} \
            --gene {wildcards.gene}
        """


# align the G gene with mafft
rule realign:
    input:
        slicedalignment=rules.cut.output.slicedalignment,
        reference=lambda w: results_dir
        + f"/{get_subtype(w.build)}/{w.gene}_reference.fasta",
    output:
        realigned=results_dir + "/{build}/{gene}_aligned.fasta",
    log:
        "logs/{build}/realign_{gene}.txt"
    benchmark:
        "benchmarks/{build}/realign_{gene}.txt"
    threads: 2
    shell:
        r"""
        exec &> >(tee {log:q})

        augur align --nthreads {threads} \
            --sequences {input.slicedalignment} \
            --reference-sequence {input.reference} \
            --output {output.realigned}
        """


rule hybrid_align:
    input:
        original=rules.genome_align.output.alignment,
        G_alignment=lambda w: results_dir + f"/{w.build}/G_aligned.fasta",
        reference=lambda w: f"config/{get_subtype(w.build)}reference.gbk",
    output:
        hybrid_alignment=results_dir
        + "/{build}/hybrid_alignment.fasta",
    log:
        "logs/{build}/hybrid_align.txt"
    benchmark:
        "benchmarks/{build}/hybrid_align.txt"
    params:
        gene_name=lambda w: get_gene_build(w.build),
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/align_for_tree.py \
            --realign {input.G_alignment} \
            --original {input.original} \
            --reference {input.reference} \
            --output {output.hybrid_alignment} \
            --build {params.gene_name}
        """


def get_alignment(w):
    gene_build = get_gene_build(w.build)
    if gene_build == "genome":
        return rules.hybrid_align.output.hybrid_alignment
    else:
        gene = get_gene(w.build)
        return (
            results_dir
            + f"/{w.build}/{gene}_aligned.fasta"
        )


def get_reference(w):
    gene_build = get_gene_build(w.build)
    subtype = get_subtype(w.build)
    if gene_build == "genome":
        return f"config/{subtype}reference.gbk"
    else:
        gene = get_gene(w.build)
        return results_dir + f"/{subtype}/{gene}_reference.gbk"


rule tree:
    """
    Building tree
    """
    input:
        alignment=get_alignment,
    output:
        tree=results_dir + "/{build}/tree_raw.nwk",
    log:
        "logs/{build}/tree.txt"
    benchmark:
        "benchmarks/{build}/tree.txt"
    threads: 2
    shell:
        r"""
        exec &> >(tee {log:q})

        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --tree-builder-args '-ninit 10 -n 4 -czb' \
            --nthreads {threads}
        """


rule refine:
    """
    Refining tree
      - estimate timetree
      - use {params.coalescent} coalescent timescale
      - estimate {params.date_inference} node dates
    """
    input:
        tree=rules.tree.output.tree,
        alignment=get_alignment,
        metadata=lambda w: f"results/{get_subtype(w.build)}/metadata.tsv",
    output:
        tree=results_dir + "/{build}/tree.nwk",
        node_data=results_dir + "/{build}/branch_lengths.json",
    log:
        "logs/{build}/refine.txt"
    benchmark:
        "benchmarks/{build}/refine.txt"
    params:
        coalescent=config["refine"]["coalescent"],
        clock_filter_iqd=config["refine"]["clock_filter_iqd"],
        date_inference=config["refine"]["date_inference"],
        strain_id=config["strain_id_field"],
    shell:
        r"""
        exec &> >(tee {log:q})

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
        (distance_map_config["a_or_b"] == get_subtype(wildcards.build))
        & (distance_map_config["build_name"] == get_gene(wildcards.build))
        & (distance_map_config["resolution"] == get_resolution(wildcards.build))
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
    if (gene := get_gene(wildcards.build)) != "G":
        return [
            f"config/distance_maps/{gene}/{distance_map}.json"
            for distance_map in distance_config.loc[:, "distance_map"].values
        ]
    else:
        return ""


rule distances:
    input:
        tree=rules.refine.output.tree,
        distance_maps = _get_distance_maps_by_lineage_and_segment,
        translations_done=results_dir + "/{build}/translations.done",
    output:
        distances= results_dir + "/{build}/distances.json"
    params:
        genes=lambda w: get_gene(w.build),
        alignments=lambda w: [f"{results_dir}/{w.build}/translations/{get_gene(w.build)}_withInternalNodes.fasta"],
        comparisons=_get_distance_comparisons_by_lineage_and_segment,
        attribute_names=_get_distance_attributes_by_lineage_and_segment,
    log:
        "logs/{build}/distances.txt"
    benchmark:
        "benchmarks/{build}/distances.txt"
    resources:
        mem_mb=8000,
        time="00:30:00",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur distance \
            --alignment {params.alignments} \
            --tree {input.tree} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --output {output.distances}
        """


rule ancestral:
    """
    Reconstructing ancestral sequences and mutations
      - inferring ambiguous mutations
    """
    input:
        tree=rules.refine.output.tree,
        alignment=get_alignment,
        translations=rules.genome_align.output.translations,
        root_sequence=get_reference,
    output:
        node_data=results_dir + "/{build}/nt_muts.json",
        translations_done=results_dir + "/{build}/translations.done",
    params:
        inference=config["ancestral"]["inference"],
        genes=lambda w: get_gene(w.build),
        output_translations=lambda w: results_dir + f"/{w.build}/translations/%GENE_withInternalNodes.fasta",
        input_translations=lambda w: results_dir + f"/{w.build}/translations/%GENE.fasta",
    log:
        "logs/{build}/ancestral.txt"
    benchmark:
        "benchmarks/{build}/ancestral.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --annotation {input.root_sequence} \
            --root-sequence {input.root_sequence} \
            --genes {params.genes} \
            --translations "{params.input_translations}" \
            --output-translations "{params.output_translations}" \
            --inference {params.inference} && touch {output.translations_done}
        """


rule translate:
    """
    Translating amino acid sequences
    """
    input:
        tree=rules.refine.output.tree,
        node_data=rules.ancestral.output.node_data,
        reference=get_reference,
    output:
        node_data=results_dir + "/{build}/aa_muts.json",
    log:
        "logs/{build}/translate.txt"
    benchmark:
        "benchmarks/{build}/translate.txt"
    params:
        alignment_file_mask=results_dir + "/{build}/aligned_%GENE.fasta",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
            --alignment-output {params.alignment_file_mask}
        """


rule compute_f_scores_node_data:
    """
    Computing F protein antibody escape scores for all tree nodes
    """
    input:
        tree_newick=rules.refine.output.tree,
        aa_muts=rules.translate.output.node_data,
        f_scores=config["f_dms_data"],
    output:
        f_scores_node_data=results_dir + "/{build}/f_scores.json"
    params:
        gene="F",
        f_antibodies=lambda w: " ".join(shlex.quote(ab) for ab in config["f_dms_antibodies"]),
        only_positive_escape=config["dms_only_positive_escape"],
    log:
        "logs/{build}/compute_f_scores_node_data.txt"
    benchmark:
        "benchmarks/{build}/compute_f_scores_node_data.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/score_f_sequences.py tree \
            --tree-newick {input.tree_newick} \
            --aa-muts {input.aa_muts} \
            --gene {params.gene} \
            --dms-scores {input.f_scores} \
            --dms-antibodies {params.f_antibodies} \
            --only-positive-escape {params.only_positive_escape} \
            --output {output.f_scores_node_data}
        """


rule traits:
    input:
        tree=rules.refine.output.tree,
        metadata=lambda w: f"results/{get_subtype(w.build)}/metadata.tsv",
    output:
        node_data=results_dir + "/{build}/traits.json",
    log:
        "logs/{build}/traits.txt"
    benchmark:
        "benchmarks/{build}/traits.txt"
    params:
        columns=config["traits"]["columns"],
        strain_id=config["strain_id_field"],
    shell:
        r"""
        exec &> >(tee {log:q})

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
        metadata = lambda w: f"results/{get_subtype(w.build)}/metadata.tsv"
    output:
        frequencies = results_dir + "/{build}/frequencies.json"
    log:
        "logs/{build}/frequencies.txt"
    benchmark:
        "benchmarks/{build}/frequencies.txt"
    params:
        min_date_arg = lambda w: f"--min-date {config['frequencies']['resolutions'][get_resolution(w.build)]['min_date']}",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur frequencies \
            --tree {input.tree} \
            --method kde \
            --metadata-id-columns accession \
            {params.min_date_arg} \
            --metadata {input.metadata} \
            --output {output.frequencies}
        """
