def get_node_data(w):
    node_data = [rules.refine.output.node_data,
                    rules.traits.output.node_data,
                    rules.ancestral.output.node_data,
                    rules.translate.output.node_data]
    if w.build_name in config["genesforglycosylation"]:
        node_data.append(rules.glycosylation.output.glycosylations)
    if w.build_name == "genome":
        node_data.append(rules.clades_genome.output.node_data)
    if w.build_name in ["genome", "G"]:
        node_data.append(rules.clades.output.node_data)
    return node_data

rule colors:
    input:
        color_schemes = "config/color_schemes.tsv",
        color_orderings = "config/color_orderings.tsv",
        metadata = "data/{a_or_b}/metadata.tsv",
    output:
        colors = "results/{a_or_b}/{build_name}/colors.tsv"
    shell:
        """
        python scripts/assign-colors.py \
            --color-schemes {input.color_schemes} \
            --ordering {input.color_orderings} \
            --metadata {input.metadata} \
            --output {output.colors}
        """

rule export:
    message: "Exporting data files for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.filter.input.metadata,
        node_data = get_node_data,
        colors = rules.colors.output.colors,
        auspice_config = config["files"]["auspice_config"],
        description = config["description"]
    output:
        auspice_json =  build_dir + "/{a_or_b}/{build_name}/tree.json",
        root_sequence = build_dir + "/{a_or_b}/{build_name}/tree_root-sequence.json"
    params:
    	title = lambda w: f"RSV-{w.a_or_b.upper()} phylogeny"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --title {params.title:q} \
            --description {input.description} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence \
            --output {output.auspice_json}
        """

rule final_strain_name:
    input:
        auspice_json= rules.export.output.auspice_json,
        metadata = rules.filter.input.metadata,
    output:
        auspice_json=build_dir + "/{a_or_b}/{build_name}/tree_renamed.json"
    params:
        display_strain_field=lambda w: config.get("display_strain_field", "strain"),
    shell:
        """
        python3 scripts/set_final_strain_name.py --metadata {input.metadata} \
                --input-auspice-json {input.auspice_json} \
                --display-strain-name {params.display_strain_field} \
                --output {output.auspice_json}
        """


rule rename_clade_labels:
    input:
        auspice_json= rules.final_strain_name.output.auspice_json,
        root_sequence= rules.export.output.root_sequence
    output:
        auspice_json= "auspice/rsv_{a_or_b}_{build_name}.json",
        root_sequence= "auspice/rsv_{a_or_b}_{build_name}_root-sequence.json"
    shell:
        """
        python3 scripts/clade_names.py  \
                --input-auspice-json {input.auspice_json} \
                --output {output.auspice_json}

        cp {input.root_sequence} {output.root_sequence}
        """
