def get_node_data(w):
    node_data = [rules.refine.output.node_data,
                    rules.traits.output.node_data,
                    rules.ancestral.output.node_data,
                    rules.translate.output.node_data]
    if config["gene"]=='G':
        node_data.append(rules.glycosylation.output.node_data)
    return node_data


rule export:
    message: "Exporting data files for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.traits.input.metadata,
        node_data = get_node_data,
        auspice_config = config["files"]["auspice_config"]
    output:
        auspice_json = "auspice/{build_name}.json"
    params:
    	title = "RSV-A phylogeny"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --title {params.title:q} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """
