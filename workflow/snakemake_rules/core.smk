'''
This part of the workflow expects input files
            sequences = "data/sequences.fasta"
            metadata = "data/metadata.tsv"
'''


rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences = "data/{a_or_b}/sequences.fasta"
    output:
        sequence_index = build_dir + "/{a_or_b}/{build_name}/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule newreference:
    message:
        """
        Making new reference
        """
    input:
        oldreference = "config/{a_or_b}reference.gbk"
    output:
        newreferencegbk = build_dir + "/{a_or_b}/{build_name}/newreference.gbk",
        newreferencefasta = build_dir + "/{a_or_b}/{build_name}/newreference.fasta",
        greference = build_dir + "/{a_or_b}/{build_name}/greference.fasta"
    params:
        gene = lambda w: w.build_name,
        newreference = build_dir + "/{a_or_b}/{build_name}/newreference",
        oldreference = 'config/{a_or_b}reference',
        greference = build_dir + "/{a_or_b}/{build_name}/greference"
    shell:
        """
        python scripts/newreference.py \
            --greference {params.greference} \
            --reference {params.oldreference} \
            --output {params.newreference} \
            --gene {params.gene}
        """


rule filter:
    message:
        """
        Aligning sequences to {input.reference}
            - gaps relative to reference are considered real
        """
    input:
        sequences = "data/{a_or_b}/sequences.fasta",
        reference = "config/{a_or_b}reference.gbk",
        metadata = "data/{a_or_b}/metadata.tsv",
        sequence_index = rules.index_sequences.output
    output:
    	sequences = build_dir + "/{a_or_b}/{build_name}/filtered.fasta"
    params:
    	group_by = config["filter"]["group_by"],
    	min_length = lambda w: config["filter"]["min_length"].get(w.build_name, 10000),
    	subsample_max_sequences = config["filter"]["subsample_max_sequences"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --min-length {params.min_length}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = "config/{a_or_b}reference.fasta"
    output:
        alignment = build_dir + "/{a_or_b}/{build_name}/sequences.aligned.fasta"
    shell:
        """
        nextalign run \
            --reference {input.reference} \
            --output-fasta {output.alignment} \
            {input.sequences}
        """

rule cut:
    input:
        oldalignment = rules.align.output.alignment,
        reference = "config/{a_or_b}reference.gbk"
    output:
        newalignment = build_dir + "/{a_or_b}/{build_name}/newalignment.fasta"
    shell:
        """
        python scripts/cut.py \
            --oldalignment {input.oldalignment} \
            --newalignment {output.newalignment} \
            --reference {input.reference} \
        """

rule realign:
    input:
        newalignment = rules.cut.output.newalignment,
        reference = rules.newreference.output.greference
    output:
        realigned = build_dir + "/{a_or_b}/{build_name}/realigned.fasta"
    shell:
        """
        augur align \
            --sequences {input.newalignment} \
            --reference-sequence {input.reference} \
            --output {output.realigned}
        """

rule alignment_for_tree:
    input:
        realigned = rules.realign.output.realigned,
        original = rules.align.output.alignment,
        reference = "config/{a_or_b}reference.gbk"
    output:
        aligned_for_tree = build_dir + "/{a_or_b}/{build_name}/alignment_for_tree.fasta"
    params:
        gene = lambda w: w.build_name
    shell:
        """
        python scripts/align_for_tree.py \
            --realign {input.realigned} \
            --original {input.original} \
            --reference {input.reference} \
            --output {output.aligned_for_tree} \
            --gene {params.gene}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.alignment_for_tree.output.aligned_for_tree
    output:
        tree = build_dir + "/{a_or_b}/{build_name}/tree_raw.nwk"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads 4
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.alignment_for_tree.output.aligned_for_tree,
        metadata = rules.filter.input.metadata
    output:
        tree = build_dir + "/{a_or_b}/{build_name}/tree.nwk",
        node_data = build_dir + "/{a_or_b}/{build_name}/branch_lengths.json"
    params:
    	coalescent = config["refine"]["coalescent"],
    	clock_filter_iqd = config["refine"]["clock_filter_iqd"],
    	date_inference = config["refine"]["date_inference"]
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --date-confidence \
            --timetree \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - inferring ambiguous mutations
        """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.alignment_for_tree.output.aligned_for_tree
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/nt_muts.json"
    params:
    	inference = config["ancestral"]["inference"]
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = rules.newreference.output.newreferencegbk
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/aa_muts.json",
        aa_data = build_dir + "/{a_or_b}/{build_name}/alignedG.fasta"
    params:
    	alignment_file_mask = build_dir + "/{a_or_b}/{build_name}/aligned%GENE.fasta"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
            --alignment-output {params.alignment_file_mask}
        """

rule traits:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.filter.input.metadata
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/traits.json"
    log:
        "logs/{a_or_b}/traits_{build_name}_rsv.txt"
    params:
    	columns = config["traits"]["columns"]
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule clades:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = "config/clades_G_{a_or_b}.tsv"
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/clades_G.json"
    log:
        "logs/{a_or_b}/clades_{build_name}.txt"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """