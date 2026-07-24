def get_node_data(w):
    node_data = [rules.refine.output.node_data,
                    rules.traits.output.node_data,
                    rules.ancestral.output.node_data,
                    rules.translate.output.node_data]
    gene_build = get_gene_build(w.build)
    gene = get_gene(w.build)
    if gene_build in config["genesforglycosylation"]:
        node_data.append(rules.glycosylation.output.glycosylations)
    if gene_build == "genome":
        node_data.append(rules.clades_consortium.output.node_data)
    if gene == "F":
        node_data.append(rules.distances.output.distances)
        node_data.append(rules.compute_f_scores_node_data.output.f_scores_node_data)

    return node_data

rule generate_f_dms_antibody_auspice_config:
    """
    Generating auspice config for F DMS antibody colorings
    """
    input:
        f_scores_node_data = rules.compute_f_scores_node_data.output.f_scores_node_data
    output:
        auspice_config = "results/{build}/auspice_config_f_dms_antibodies.json"
    log:
        "logs/{build}/generate_f_dms_antibody_auspice_config.txt"
    benchmark:
        "benchmarks/{build}/generate_f_dms_antibody_auspice_config.txt"
    params:
        antibodies = config["f_dms_antibodies"],
        continuous_scale = " ".join(
            shlex.quote(hexcode)
            for hexcode in
            ["#440154", "#472d7b", "#3b528b", "#2c728e", "#21918c", "#28ae80", "#5ec962", "#addc30", "#fde725"]
        )
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/generate_f_dms_antibody_auspice_config.py \
            --antibodies {params.antibodies} \
            --continuous-scale {params.continuous_scale} \
            --node-data {input.f_scores_node_data} \
            --output {output.auspice_config}
        """

rule colors:
    input:
        color_schemes = "config/color_schemes.tsv",
        color_orderings = "config/color_orderings.tsv",
        metadata = "results/{a_or_b}/metadata.tsv",
    output:
        colors = "results/{a_or_b}/colors.tsv"
    log:
        "logs/colors_{a_or_b}.txt"
    benchmark:
        "benchmarks/colors_{a_or_b}.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        python scripts/assign-colors.py \
            --color-schemes {input.color_schemes} \
            --ordering {input.color_orderings} \
            --metadata {input.metadata} \
            --output {output.colors}
        """

def auspice_configs(wildcards):
    """
    Returns a list of auspice configs which change depending on the wildcards.
    This is preferable to using command-line args to define per-build changes as
    these tend to clobber the data in the provided config JSON.
    (Note: multiple config file support requires Augur v29.1.0)
    """
    configs = [config["files"]["auspice_config"]] # base config
    gene_build = get_gene_build(wildcards.build)
    gene = get_gene(wildcards.build)
    if gene_build != 'genome':
        configs.append(config['files']['auspice_config_additional_colorings'])
    if gene == "F":
        configs.append(
            f"results/{wildcards.build}/auspice_config_f_dms_antibodies.json"
        )
    if gene_build == "F-antibody-escape":
        configs.append(config['files']['auspice_config_f_antibody_escape'])
    elif gene_build not in ["genome", "G"]:
        configs.append(config['files']['auspice_config_non-genome_builds'])
    return configs

rule export:
    """
    Exporting data files for auspice
    """
    input:
        tree = rules.refine.output.tree,
        metadata = lambda w: f"results/{get_subtype(w.build)}/metadata.tsv",
        node_data = get_node_data,
        colors = lambda w: f"results/{get_subtype(w.build)}/colors.tsv",
        auspice_config = auspice_configs,
        description = config["description"]
    output:
        auspice_json = results_dir + "/{build}/tree.json"
    log:
        "logs/{build}/export.txt"
    benchmark:
        "benchmarks/{build}/export.txt"
    params:
        title = lambda w: f"RSV-{get_subtype(w.build).upper()} phylogeny",
        strain_id=config["strain_id_field"],
    shell:
        r"""
        exec &> >(tee {log:q})

        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --node-data {input.node_data} \
            --title {params.title:q} \
            --description {input.description} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence-inline \
            --validation warn \
            --output {output.auspice_json}
        """


rule final_strain_name:
    input:
        auspice_json= rules.export.output.auspice_json,
        metadata = lambda w: f"results/{get_subtype(w.build)}/metadata.tsv",
        frequencies = results_dir + "/{build}/frequencies.json"
    output:
        auspice_json=results_dir + "/{build}/tree_renamed.json",
        freq_json= results_dir + "/{build}/tip-frequencies.json"
    log:
        "logs/{build}/final_strain_name.txt"
    benchmark:
        "benchmarks/{build}/final_strain_name.txt"
    params:
        strain_id=config["strain_id_field"],
        display_strain_field=config["display_strain_field"],
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/set_final_strain_name.py --metadata {input.metadata} \
                --metadata-id-columns {params.strain_id} \
                --input-auspice-json {input.auspice_json} \
                --input-frequency-json {input.frequencies} \
                --display-strain-name {params.display_strain_field} \
                --output-auspice-json {output.auspice_json}\
                --output-frequencies-json {output.freq_json}
        """


rule rename_and_ready_for_nextclade:
    input:
        auspice_json= rules.final_strain_name.output.auspice_json,
        pathogen_json= "nextclade/config/pathogen.json",
    output:
        auspice_json= results_dir + "/{build}/rsv.json",
    log:
        "logs/{build}/rename_and_ready_for_nextclade.txt"
    benchmark:
        "benchmarks/{build}/rename_and_ready_for_nextclade.txt"
    params:
        accession= lambda w: config["nextclade_attributes"][get_subtype(w.build)]["accession"],
        name= lambda w: config["nextclade_attributes"][get_subtype(w.build)]["name"],
        ref_name= lambda w: config["nextclade_attributes"][get_subtype(w.build)]["reference_name"]
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/rename_and_nextclade.py  \
                --input-auspice-json {input.auspice_json} \
                --pathogen-json {input.pathogen_json} \
                --reference {params.ref_name:q} \
                --build-name {params.name:q} \
                --reference-accession {params.accession:q} \
                --output {output.auspice_json}
        """


rule copy_export:
    input:
        auspice_json = lambda w: f"results/{w.build_with_underscores.replace('_', '/')}/rsv.json",
        tip_freq = lambda w: f"results/{w.build_with_underscores.replace('_', '/')}/tip-frequencies.json",
    output:
        auspice_json = "auspice/rsv_{build_with_underscores}.json",
        tip_freq = "auspice/rsv_{build_with_underscores}_tip-frequencies.json",
    shell:
        """
        cp {input.auspice_json} {output.auspice_json}
        cp {input.tip_freq} {output.tip_freq}
        """
