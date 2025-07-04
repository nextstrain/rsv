def get_node_data(w):
    node_data = [rules.refine.output.node_data,
                    rules.traits.output.node_data,
                    rules.ancestral.output.node_data,
                    rules.translate.output.node_data]
    if w.build_name in config["genesforglycosylation"]:
        node_data.append(rules.glycosylation.output.glycosylations)
    if w.build_name == "genome":
        node_data.append(rules.clades_consortium.output.node_data)
    if w.build_name in ["F", "genome"]:
        node_data.append(rules.distances.output.distances)

    return node_data

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

rule create_accession_json:
    input:
        metadata = "data/{a_or_b}/metadata.tsv",
        tree = rules.refine.output.tree,
    output:
        accession_json = build_dir + "/{a_or_b}/{build_name}/{resolution}/accession.json"
    shell:
        """
        python3 scripts/create_accession_json.py \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --output {output.accession_json}
        """


rule export:
    message: "Exporting data files for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = "data/{a_or_b}/metadata.tsv",
        node_data = get_node_data,
        colors = rules.colors.output.colors,
        auspice_config = config["files"]["auspice_config"],
        description = config["description"]
    output:
        auspice_json =  build_dir + "/{a_or_b}/{build_name}/{resolution}/tree_prelim.json"
    params:
        title = lambda w: f"RSV-{w.a_or_b.upper()} phylogeny",
        strain_id=config["strain_id_field"],
        metadata_colors = lambda w: '' if w.build_name=='genome' else f"--color-by-metadata clade"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} {params.metadata_colors} \
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

# add PPX specific fields to the exported tree that augur export v2 does not support
# for genbank only, rename output file of rule export to tree.json
rule patch_exported_tree:
    input:
        auspice_json = build_dir + "/{a_or_b}/{build_name}/{resolution}/tree_prelim.json",
        accessions = rules.create_accession_json.output.accession_json,
    output:
        auspice_json = build_dir + "/{a_or_b}/{build_name}/{resolution}/tree.json",
    run:
        import json

        with open(input.auspice_json) as f:
            tree = json.load(f)

        with open(input.accessions) as f:
            accessions = json.load(f)['nodes']


        def update_node_data(node):
            name = node["name"]
            if name in accessions:
                node["node_attrs"] |= accessions[name]
            if "children" in node:
                for child in node["children"]:
                    update_node_data(child)

        update_node_data(tree["tree"])

        with open(output.auspice_json, "w") as f:
            json.dump(tree, f, indent=2)



rule final_strain_name:
    input:
        auspice_json= build_dir + "/{a_or_b}/{build_name}/{resolution}/tree.json",
        metadata = "data/{a_or_b}/metadata.tsv",
        frequencies = build_dir + "/{a_or_b}/{build_name}/{resolution}/frequencies.json"
    output:
        auspice_json=build_dir + "/{a_or_b}/{build_name}/{resolution}/tree_renamed.json",
        freq_json= "auspice/rsv_{a_or_b}_{build_name}_{resolution}_tip-frequencies.json"
    params:
        strain_id=config["strain_id_field"],
        display_strain_field=config.get("display_strain_field", "strain"),
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
        auspice_json= build_dir + "/{a_or_b}/{build_name}/{resolution}/tree_renamed.json",
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
