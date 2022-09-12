rule clades_genome:
    message:
        "adding clades based on the entire genome"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = "config/clades_genome_{a_or_b}.tsv"
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/clades_genome.json"
    log:
        "logs/{a_or_b}/clades_genome_{build_name}.txt"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data tmp.json 2>&1 | tee {log}
        sed -e 's/clade_membership/genome_clade/g' tmp.json |
        sed -e 's/clade_annotation/genome_clade_annotation/g' > {output.node_data}
        """
