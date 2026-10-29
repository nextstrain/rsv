"""
This part of the workflow expects input files
            sequences = "data/sequences.fasta"
            metadata = "data/metadata.tsv"
"""
from augur.subsample import get_referenced_files

rule filter:
    input:
        sequences="results/{a_or_b}/sequences.fasta",
        metadata="results/{a_or_b}/metadata.tsv",
        config="results/{a_or_b}/config/filter.yaml",
        referenced_files=lambda w: get_referenced_files(f"results/{w.a_or_b}/config/filter.yaml"),
    output:
        sequences="results/{a_or_b}/sequences_filtered.fasta",
        metadata="results/{a_or_b}/metadata_filtered.tsv",
    log:
        "logs/filter_{a_or_b}.txt"
    benchmark:
        "benchmarks/filter_{a_or_b}.txt"
    params:
        strain_id=config["strain_id_field"],
    shell:
        r"""
        exec &> >(tee {log:q})

        augur subsample \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --config {input.config} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
        """


rule index_sequences:
    """
    Creating an index of sequence composition for filtering.
    """
    input:
        sequences="results/{a_or_b}/sequences_filtered.fasta",
    output:
        sequence_index=build_dir
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
        newreferencegbk=build_dir
        + "/{a_or_b}/{gene}_reference.gbk",
        newreferencefasta=build_dir
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


rule cut:
    wildcard_constraints:
        build_name="|".join(b for b in config["cds"] if b != "genome"),
    input:
        oldalignment=build_dir + "/{a_or_b}/sequences_aligned.fasta",
        reference="config/{a_or_b}reference.gbk",
    output:
        slicedalignment=build_dir + "/{a_or_b}/{build_name}/sequences_aligned.fasta",
    log:
        "logs/cut_{a_or_b}_{build_name}.txt"
    benchmark:
        "benchmarks/cut_{a_or_b}_{build_name}.txt"
    params:
        gene=lambda w: config["cds"][w.build_name],
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/cut.py \
            --oldalignment {input.oldalignment} \
            --slicedalignment {output.slicedalignment} \
            --reference {input.reference} \
            --gene {params.gene}
        """


def get_pre_subsample_alignment(w):
    if w.build_name == "genome":
        return build_dir + f"/{w.a_or_b}/sequences_aligned.fasta"
    else:
        return rules.cut.output.slicedalignment


rule subsample:
    input:
        sequences=get_pre_subsample_alignment,
        metadata="results/{a_or_b}/metadata_filtered.tsv",
        sequence_index=rules.index_sequences.output,
        config="results/{a_or_b}/{build_name}/{resolution}/subsample_config.yaml",
        referenced_files=lambda w: get_referenced_files(f"results/{w.a_or_b}/{w.build_name}/{w.resolution}/subsample_config.yaml"),
    output:
        sequences=build_dir + "/{a_or_b}/{build_name}/{resolution}/subsampled.fasta",
    log:
        "logs/subsample_{a_or_b}_{build_name}_{resolution}.txt",
    benchmark:
        "benchmarks/subsample_{a_or_b}_{build_name}_{resolution}.txt",
    params:
        strain_id=config["strain_id_field"],
    shell:
        r"""
        exec &> >(tee {log:q})

        augur subsample \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --config {input.config} \
            --output-sequences {output.sequences}
        """

rule combine_samples:
    input:
        subsamples=lambda w: (
            [rules.subsample.output.sequences]
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


rule filter_for_f_antibody_escape:
    """
    Do the quality filtering applied to each sequence set before subsampling
    """
    input:
        sequences="results/{a_or_b}/sequences_filtered.fasta",
        metadata="results/{a_or_b}/metadata_filtered.tsv",
    output:
        sequences=build_dir + "/{a_or_b}/{build_name}/{resolution}/pre_subsample/filtered.fasta",
    log:
        "logs/filter_for_f_antibody_escape_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/filter_for_f_antibody_escape_{a_or_b}_{build_name}_{resolution}.txt"
    params:
        min_coverage=f'F_coverage>{config["filter_for_f_antibody_escape"]["min_coverage"]}',
        min_length=config["filter_for_f_antibody_escape"]["min_length"],
        strain_id=config["strain_id_field"],
        min_date=lambda w: config["filter_for_f_antibody_escape"]["resolutions"][w.resolution]["min_date"],
    shell:
        r"""
        exec &> >(tee {log:q})

        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --min-length {params.min_length} \
            --min-date {params.min_date} \
            --output {output.sequences} \
            --query '{params.min_coverage}'
        """


rule filter_f_translations:
    """
    Filtering F translations to pre-subsampled quality-filtered sequences
    """
    input:
        translations_done=build_dir + "/{a_or_b}/translations.done",
        filtered_sequences=rules.filter_for_f_antibody_escape.output.sequences,
    output:
        f_sequences=build_dir
        + "/{a_or_b}/{build_name}/{resolution}/pre_subsample/translations/F.fasta",
    log:
        "logs/filter_f_translations_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/filter_f_translations_{a_or_b}_{build_name}_{resolution}.txt"
    params:
        f_translations=lambda w: f"results/{w.a_or_b}/translations/F.fasta",
    shell:
        r"""
        exec &> >(tee {log:q})

        seqkit grep -f <(seqkit seq -n -i {input.filtered_sequences}) {params.f_translations} > {output.f_sequences}
        """


rule score_pre_subsample_f_proteins:
    """
    Computing F protein DMS scores for pre-subsampled sequences
    """
    input:
        f_sequences=rules.filter_f_translations.output.f_sequences,
        dms_scores=config["f_dms_data"],
    output:
        scores=build_dir + "/{a_or_b}/{build_name}/{resolution}/pre_subsample/f_protein_scores.tsv",
    log:
        "logs/score_pre_subsample_f_proteins_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/score_pre_subsample_f_proteins_{a_or_b}_{build_name}_{resolution}.txt"
    params:
        dms_antibodies=lambda w: " ".join(shlex.quote(ab) for ab in config["f_dms_antibodies"]),
        only_positive_escape=config["dms_only_positive_escape"],
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/score_f_sequences.py fasta \
            --sequences {input.f_sequences} \
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
        original_metadata="results/{a_or_b}/metadata_filtered.tsv",
        f_scores=rules.score_pre_subsample_f_proteins.output.scores,
    output:
        enhanced_metadata=build_dir + "/{a_or_b}/{build_name}/{resolution}/pre_subsample/metadata_with_scores.tsv",
    log:
        "logs/add_f_scores_to_pre_subsample_metadata_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/add_f_scores_to_pre_subsample_metadata_{a_or_b}_{build_name}_{resolution}.txt"
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
        sequences=rules.cut.output.slicedalignment,
        metadata=rules.add_f_scores_to_pre_subsample_metadata.output.enhanced_metadata,
    output:
        sequences=build_dir + "/{a_or_b}/{build_name}/{resolution}/filtered_{antibody}_{scoretype}.fasta"
    log:
        "logs/enrich_antibody_escape_{a_or_b}_{build_name}_{resolution}_{antibody}_{scoretype}.txt"
    benchmark:
        "benchmarks/enrich_antibody_escape_{a_or_b}_{build_name}_{resolution}_{antibody}_{scoretype}.txt"
    params:
        strain_id=config["strain_id_field"],
        escape_col=lambda w: f"{w.antibody}_{w.scoretype}",
        nseqs=lambda w: config["enrich_antibody_escape"][w.build_name]["nseqs_per_antibody_scoretype"],
        group_by=lambda w: " ".join(shlex.quote(g) for g in config["enrich_antibody_escape"][w.build_name]["group_by"]),
        max_identical_f_prot_muts=lambda w: config["enrich_antibody_escape"][w.build_name]["max_identical_f_prot_muts"],
        max_identical_max_escape_mut=lambda w: config["enrich_antibody_escape"][w.build_name]["max_identical_max_escape_mut"],
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


rule nextclade_align:
    """
    Aligning all filtered sequences to the reference
    """
    input:
        sequences="results/{a_or_b}/sequences_filtered.fasta",
        dataset=rules.get_nextclade_dataset.output.dataset,
    output:
        alignment=build_dir + "/{a_or_b}/sequences_aligned.fasta",
        translations=directory(build_dir + "/{a_or_b}/translations"),
        translations_done=build_dir + "/{a_or_b}/translations.done",
    params:
        genes="G,F",
    # This is a bottleneck step with at most two (one per subtype) running at a time.
    threads: max(1, workflow.cores // 2)
    log:
        "logs/nextclade_align_{a_or_b}.txt"
    benchmark:
        "benchmarks/nextclade_align_{a_or_b}.txt"
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


# cut out the G-Gene for alignment refinement
rule cut_G:
    input:
        oldalignment=build_dir + "/{a_or_b}/genome/{resolution}/filtered.fasta",
        reference="config/{a_or_b}reference.gbk",
    output:
        slicedalignment=build_dir
        + "/{a_or_b}/genome/{resolution}/G_slicedalignment.fasta",
    log:
        "logs/cut_{a_or_b}_genome_{resolution}_G.txt"
    benchmark:
        "benchmarks/cut_{a_or_b}_genome_{resolution}_G.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/cut.py \
            --oldalignment {input.oldalignment} \
            --slicedalignment {output.slicedalignment} \
            --reference {input.reference} \
            --gene G
        """


# align the G gene with mafft
rule realign_G:
    wildcard_constraints:
        build_name="genome|G",
    input:
        slicedalignment=lambda w: (
            rules.cut_G.output.slicedalignment
            if w.build_name == "genome"
            else rules.combine_samples.output.sequences
        ),
        reference=build_dir
        + "/{a_or_b}/G_reference.fasta",
    output:
        realigned=build_dir + "/{a_or_b}/{build_name}/{resolution}/G_aligned.fasta",
    log:
        "logs/realign_{a_or_b}_{build_name}_{resolution}_G.txt"
    benchmark:
        "benchmarks/realign_{a_or_b}_{build_name}_{resolution}_G.txt"
    threads: 2
    shell:
        r"""
        exec &> >(tee {log:q})

        augur align --nthreads {threads} \
            --sequences {input.slicedalignment} \
            --reference-sequence {input.reference} \
            --output {output.realigned}
        """


rule hybrid_align_genome:
    input:
        original=build_dir + "/{a_or_b}/genome/{resolution}/filtered.fasta",
        G_alignment=build_dir + "/{a_or_b}/genome/{resolution}/G_aligned.fasta",
        reference="config/{a_or_b}reference.gbk",
    output:
        hybrid_alignment=build_dir
        + "/{a_or_b}/genome/{resolution}/hybrid_alignment.fasta",
    log:
        "logs/hybrid_align_genome_{a_or_b}_{resolution}.txt"
    benchmark:
        "benchmarks/hybrid_align_genome_{a_or_b}_{resolution}.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/align_for_tree.py \
            --realign {input.G_alignment} \
            --original {input.original} \
            --reference {input.reference} \
            --output {output.hybrid_alignment}
        """


def get_alignment(w):
    if w.build_name == "genome":
        return rules.hybrid_align_genome.output.hybrid_alignment
    elif w.build_name == "G":
        return rules.realign_G.output.realigned
    else:
        return rules.combine_samples.output.sequences


def get_reference(w):
    if w.build_name == "genome":
        return f"config/{w.a_or_b}reference.gbk"
    else:
        gene = config["cds"][w.build_name]
        return build_dir + f"/{w.a_or_b}/{gene}_reference.gbk"


rule tree:
    """
    Building tree
    """
    input:
        alignment=get_alignment,
    output:
        tree=build_dir + "/{a_or_b}/{build_name}/{resolution}/tree_raw.nwk",
    log:
        "logs/tree_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/tree_{a_or_b}_{build_name}_{resolution}.txt"
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
        metadata="results/{a_or_b}/metadata_filtered.tsv",
    output:
        tree=build_dir + "/{a_or_b}/{build_name}/{resolution}/tree.nwk",
        node_data=build_dir + "/{a_or_b}/{build_name}/{resolution}/branch_lengths.json",
    log:
        "logs/refine_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/refine_{a_or_b}_{build_name}_{resolution}.txt"
    params:
        coalescent=config["refine"]["coalescent"],
        clock_filter_iqd=config["refine"]["clock_filter_iqd"],
        date_inference=config["refine"]["date_inference"],
        strain_id=config["strain_id_field"],
        divergence_units=config["refine"]["divergence_units"],
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
            --clock-filter-iqd {params.clock_filter_iqd} \
            --divergence-units {params.divergence_units}
        """


def _get_build_distance_map_config(wildcards):
    build_gene = config["cds"][wildcards.build_name]
    distance_config = distance_map_config[
        (distance_map_config["a_or_b"] == wildcards.a_or_b)
        & (distance_map_config["build_name"] == build_gene)
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
        build_gene = config["cds"][wildcards.build_name]
        return [
            f"config/distance_maps/{build_gene}/{distance_map}.json"
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
    log:
        "logs/distances_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/distances_{a_or_b}_{build_name}_{resolution}.txt"
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
        translations=rules.nextclade_align.output.translations,
        root_sequence=get_reference,
    output:
        node_data=build_dir + "/{a_or_b}/{build_name}/{resolution}/nt_muts.json",
        translations_done=build_dir + "/{a_or_b}/{build_name}/{resolution}/translations.done",
    params:
        inference=config["ancestral"]["inference"],
        genes=lambda w: config["cds"][w.build_name],
        output_translations=lambda w: build_dir + f"/{w.a_or_b}/{w.build_name}/{w.resolution}/translations/%GENE_withInternalNodes.fasta",
        input_translations=lambda w: build_dir + f"/{w.a_or_b}/translations/%GENE.fasta",
    log:
        "logs/ancestral_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/ancestral_{a_or_b}_{build_name}_{resolution}.txt"
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
        node_data=build_dir + "/{a_or_b}/{build_name}/{resolution}/aa_muts.json",
    log:
        "logs/translate_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/translate_{a_or_b}_{build_name}_{resolution}.txt"
    params:
        alignment_file_mask=build_dir + "/{a_or_b}/{build_name}/{resolution}/aligned_%GENE.fasta",
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
        f_scores_node_data=build_dir + "/{a_or_b}/{build_name}/{resolution}/f_scores.json"
    params:
        gene="F",
        f_antibodies=lambda w: " ".join(shlex.quote(ab) for ab in config["f_dms_antibodies"]),
        only_positive_escape=config["dms_only_positive_escape"],
    log:
        "logs/compute_f_scores_node_data_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/compute_f_scores_node_data_{a_or_b}_{build_name}_{resolution}.txt"
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
        metadata="results/{a_or_b}/metadata_filtered.tsv",
    output:
        node_data=build_dir + "/{a_or_b}/{build_name}/{resolution}/traits.json",
    log:
        "logs/traits_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/traits_{a_or_b}_{build_name}_{resolution}.txt"
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
        metadata = "results/{a_or_b}/metadata_filtered.tsv"
    output:
        frequencies = build_dir + "/{a_or_b}/{build_name}/{resolution}/frequencies.json"
    log:
        "logs/frequencies_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/frequencies_{a_or_b}_{build_name}_{resolution}.txt"
    params:
        min_date_arg = lambda w: f"--min-date {config['frequencies']['resolutions'][w.resolution]['min_date']}",
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
