def get_node_data(w):
    node_data = [rules.refine.output.node_data,
                    rules.traits.output.node_data,
                    rules.ancestral.output.node_data,
                    rules.translate.output.node_data]
    if w.build_name in config["genesforglycosylation"]:
        node_data.append(rules.glycosylation.output.glycosylations)
    if w.build_name == "genome":
        node_data.append(rules.clades_consortium.output.node_data)
    if w.build_name not in ["G"]:
        node_data.append(rules.distances.output.distances)
        node_data.append(rules.compute_f_scores_node_data.output.f_scores_node_data)

    return node_data

rule generate_f_dms_antibody_auspice_config:
    message: "Generating auspice config for F DMS antibody colorings"
    input:
        f_scores_node_data = rules.compute_f_scores_node_data.output.f_scores_node_data
    output:
        auspice_config = "results/{a_or_b}/{build_name}/{resolution}/auspice_config_f_dms_antibodies.json"
    params:
        antibodies = config["f_dms_antibodies"],
        continuous_scale = " ".join(
            shlex.quote(hexcode)
            for hexcode in
            ["#440154", "#472d7b", "#3b528b", "#2c728e", "#21918c", "#28ae80", "#5ec962", "#addc30", "#fde725"]
        )
    shell:
        """
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
        metadata = "data/{a_or_b}/metadata.tsv",
    output:
        colors = "results/{a_or_b}/{build_name}/{resolution}/colors.tsv"
    shell:
        """
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
    if wildcards.build_name!='genome':
        configs.append(config['files']['auspice_config_additional_colorings'])
    if wildcards.build_name not in ["G"]:
        configs.append(
            f"results/{wildcards.a_or_b}/{wildcards.build_name}/{wildcards.resolution}/auspice_config_f_dms_antibodies.json"
        )
    return configs

rule export:
    message: "Exporting data files for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = "data/{a_or_b}/metadata.tsv",
        node_data = get_node_data,
        colors = rules.colors.output.colors,
        auspice_config = auspice_configs,
        description = config["description"]
    output:
        auspice_json =  build_dir + "/{a_or_b}/{build_name}/{resolution}/tree.json"
    params:
        title = lambda w: f"RSV-{w.a_or_b.upper()} phylogeny",
        strain_id=config["strain_id_field"],
    shell:
        """
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
        metadata = "data/{a_or_b}/metadata.tsv",
        frequencies = build_dir + "/{a_or_b}/{build_name}/{resolution}/frequencies.json"
    output:
        auspice_json=build_dir + "/{a_or_b}/{build_name}/{resolution}/tree_renamed.json",
        freq_json= "auspice/rsv_{a_or_b}_{build_name}_{resolution}_tip-frequencies.json"
    params:
        strain_id=config["strain_id_field"],
        display_strain_field=config["display_strain_field"],
    shell:
        """
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
        auspice_json= "auspice/rsv_{a_or_b}_{build_name}_{resolution}.json",
    params:
        accession= lambda w: config["nextclade_attributes"][w.a_or_b]["accession"],
        name= lambda w: config["nextclade_attributes"][w.a_or_b]["name"],
        ref_name= lambda w: config["nextclade_attributes"][w.a_or_b]["reference_name"]
    shell:
        """
        python3 scripts/rename_and_nextclade.py  \
                --input-auspice-json {input.auspice_json} \
                --pathogen-json {input.pathogen_json} \
                --reference {params.ref_name:q} \
                --build-name {params.name:q} \
                --reference-accession {params.accession:q} \
                --output {output.auspice_json}
        """
