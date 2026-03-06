rule glycosylation:
    input:
        tree = rules.refine.output.tree,
        translations = rules.translate.output.node_data
    output:
        glycosylations = build_dir + "/{a_or_b}/{build_name}/{resolution}/glyc.json"
    log:
        "logs/glycosylation_{a_or_b}_{build_name}_{resolution}.txt"
    benchmark:
        "benchmarks/glycosylation_{a_or_b}_{build_name}_{resolution}.txt"
    params:
        aa_data = build_dir + "/{a_or_b}/{build_name}/{resolution}/aligned_{build_name}.fasta"
    shell:
     r"""
     exec &> >(tee {log:q})

     python scripts/glycosylation.py \
     --alignment {params.aa_data} \
     --tree {input.tree} \
     --output {output.glycosylations}
     """
