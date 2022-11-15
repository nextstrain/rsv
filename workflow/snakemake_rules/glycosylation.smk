rule glycosylation:
    input:
        tree = rules.refine.output.tree,
        translations = rules.translate.output.node_data
    output:
        glycosylations = build_dir + "/{a_or_b}/{build_name}/glyc.json"
    params:
        aa_data = build_dir + "/{a_or_b}/{build_name}/aligned_{build_name}.fasta"
    shell:
     """
     python scripts/glycosylation.py \
     --alignment {params.aa_data} \
     --tree {input.tree} \
     --output {output.glycosylations}
     """
