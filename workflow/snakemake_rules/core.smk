'''
This part of the workflow expects input files
            sequences = "data/sequences.fasta"
            metadata = "data/metadata.tsv"
'''

build_name = config['build_name']
typ = config['rsv']


rule sequencesdeduplicated:
    input:
        allsequences = "data/rsv" + typ + "/sequences.fasta"
    output:
        sequences = build_dir + "/{build_name}/sequencesdedup.fasta"
    shell:
     """
     seqkit rmdup < {input.allsequences} > {output.sequences}
     """

rule metadatadeduplicated:
    input:
        metadata = "data/rsv" + typ + "/metadata.tsv"
    output:
        metadata = "data/rsv" + typ + "/metadatadedup.tsv"
    shell:
        """
        python scripts/metadatadedup.py \
            --metadataoriginal {input.metadata} \
            --metadataoutput {output.metadata}
        """

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences = rules.sequencesdeduplicated.output.sequences
    output:
        sequence_index = build_dir + "/{build_name}/sequence_index.tsv"
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
        oldreference = "config/" + reference +".gbk"
    output:
        newreferencegbk = build_dir + "/{build_name}/newreference.gbk",
        newreferencefasta = build_dir + "/{build_name}/newreference.fasta"
    params:
        gene = config['gene'],
        newreference = build_dir + "/{build_name}/newreference",
        oldreference = 'config/' + reference
    shell:
        """
        python scripts/newreference.py \
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
        sequences = build_dir + "/{build_name}/sequencesdedup.fasta",
        reference = "config/" + reference + ".gbk",
        metadata = rules.metadatadeduplicated.output.metadata,
        sequence_index = rules.index_sequences.output
    output:
    	sequences = build_dir + "/{build_name}/filtered.fasta"
    params:
    	group_by = config["filter"]["group_by"],
    	min_length =config["filter"]["min_length"],
    	sequences_per_group = config["filter"]["sequences_per_group"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-length {params.min_length}
        """


rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = rules.newreference.output.newreferencefasta
    output:
        alignment = build_dir + "/{build_name}/sequences.aligned.fasta"
    shell:
        """
        nextalign run \
            --reference {input.reference} \
            --output-fasta {output.alignment} \
            {input.sequences}
        """

rule sorted:
    input:
        sequences = rules.align.output.alignment,
        reference = rules.newreference.output.newreferencegbk
    output:
        alignedandsorted = build_dir + "/{build_name}/alignedandsorted.fasta"
    shell:
        """
        python scripts/sort.py \
            --alignment {input.sequences} \
            --reference {input.reference} \
            --output {output.alignedandsorted}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.sorted.output.alignedandsorted
    output:
        tree = build_dir + "/{build_name}/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
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
        alignment = rules.sorted.output.alignedandsorted,
        metadata = rules.metadatadeduplicated.output
    output:
        tree = build_dir + "/{build_name}/tree.nwk",
        node_data = build_dir + "/{build_name}/branch_lengths.json"
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
        alignment = rules.align.output.alignment
    output:
        node_data = build_dir + "/{build_name}/nt_muts.json"
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
        node_data = build_dir + "/{build_name}/aa_muts.json"
    params:
    	alignment_file_mask = build_dir + "/{build_name}/alignedandsorted%GENE.fasta",
        aa_data = expand(build_dir + "/{{build_name}}/alignedandsorted{gene}.fasta", gene=config["genesforglycosylation"])
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
        metadata = rules.metadatadeduplicated.output
    output:
        node_data = build_dir + "/{build_name}/traits.json"
    log:
        "logs/traits_{build_name}_rsv.txt"
    conda:
        config["conda_environment"]
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

