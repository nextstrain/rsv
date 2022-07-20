rule glycosylation:
    input:
        tree = rules.refine.output.tree,
    output:
        glycosylations = build_dir + "/{build_name}/{typ}/glyc.json"
    params:
        alignment = rules.translate.params.aa_data

    shell:
     """
     python scripts/glycosylation.py \
     --alignment {params.alignment} \
     --tree {input.tree} \
     --output {output.glycosylations}
     """
     
