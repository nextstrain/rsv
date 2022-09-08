rule glycosylation:
    input:
        tree = rules.refine.output.tree, 
        alignment = rules.translate.output.aa_data
    output:
        glycosylations = build_dir + "/{a_or_b}/{build_name}/glyc.json"

    shell:
     """
     python scripts/glycosylation.py \
     --alignment {input.alignment} \
     --tree {input.tree} \
     --output {output.glycosylations}
     """
