rule glycosylation:
    input:
        tree = rules.refine.output.tree,
        translations = rules.translate.output.node_data
    output:
        glycosylations = results_dir + "/{build}/glyc.json"
    log:
        "logs/{build}/glycosylation.txt"
    benchmark:
        "benchmarks/{build}/glycosylation.txt"
    params:
        aa_data = lambda w: results_dir + f"/{w.build}/aligned_{get_gene(w.build)}.fasta"
    shell:
     r"""
     exec &> >(tee {log:q})

     python scripts/glycosylation.py \
     --alignment {params.aa_data} \
     --tree {input.tree} \
     --output {output.glycosylations}
     """
