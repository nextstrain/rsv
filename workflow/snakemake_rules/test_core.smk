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
        newreferencegbk = build_dir + "/{a_or_b}/{build_name}/{gene}_reference.gbk",
        newreferencefasta = build_dir + "/{a_or_b}/{build_name}/{gene}_reference.fasta",
    params:
        gene = lambda w: w.gene,
    shell:
        """
        python scripts/newreference.py \
            --reference {input.oldreference} \
            --output-genbank {output.newreferencegbk} \
            --output-fasta {output.newreferencefasta} \
            --gene {params.gene}
        """


rule filter:
    message:
        """
        filtering sequences
        """
    input:
        sequences = "data/{a_or_b}/sequences.fasta",
        reference = "config/{a_or_b}reference.gbk",
        metadata = "data/{a_or_b}/metadata.tsv",
        sequence_index = rules.index_sequences.output,
        exclude = config['exclude']
    output:
    	sequences = build_dir + "/{a_or_b}/{build_name}/filtered.fasta"
    params:
        group_by = config["filter"]["group_by"],
        min_coverage = lambda w: f'{w.build_name}_coverage>{config["filter"]["min_coverage"].get(w.build_name, 10000)}',
        subsample_max_sequences = lambda w: config["filter"]["subsample_max_sequences"].get(w.build_name, 1000),
        strain_id=config["strain_id_field"],
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --exclude {input.exclude} \
            --exclude-where 'qc.overallStatus=bad' \
            --output sample_strains_state.txt \
            --group-by {params.group_by} \
            --subsample-max-sequences 500 \
            --query "{params.min_coverage} & division == 'Washington'"
            
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --exclude {input.exclude} \
            --exclude-where 'qc.overallStatus=bad' \
            --output sample_strains_country.txt \
            --group-by {params.group_by} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --query "{params.min_coverage} & division != 'Washington'"
            
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --exclude-all \
            --include sample_strains_state.txt \
                      sample_strains_country.txt \
            --output {output.sequences} \
        """

rule genome_align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = build_dir + "/{a_or_b}/{build_name}/genome_reference.fasta"
    output:
        alignment = build_dir + "/{a_or_b}/{build_name}/sequences.aligned.fasta"
    threads: 4
    shell:
        """
        nextclade3 run -j {threads}\
            --input-ref {input.reference} \
            --output-fasta {output.alignment} \
            {input.sequences}
        """

# cut out the G-Gene for alignment refinement
rule cut:
    input:
        oldalignment = rules.genome_align.output.alignment,
        reference = "config/{a_or_b}reference.gbk"
    output:
        slicedalignment = build_dir + "/{a_or_b}/{build_name}/{gene}_slicedalignment.fasta"
    params:
        gene = lambda w: w.gene
    shell:
        """
        python scripts/cut.py \
            --oldalignment {input.oldalignment} \
            --slicedalignment {output.slicedalignment} \
            --reference {input.reference} \
            --gene {params.gene}
        """

# align the G gene with mafft
rule realign:
    input:
        slicedalignment = rules.cut.output.slicedalignment,
        reference = build_dir + "/{a_or_b}/{build_name}/{gene}_reference.fasta"
    output:
        realigned = build_dir + "/{a_or_b}/{build_name}/{gene}_aligned.fasta"
    threads: 4
    shell:
        """
        augur align --nthreads {threads} \
            --sequences {input.slicedalignment} \
            --reference-sequence {input.reference} \
            --output {output.realigned}
        """


rule hybrid_align:
    input:
        original = rules.genome_align.output.alignment,
        G_alignment = build_dir + "/{a_or_b}/{build_name}/G_aligned.fasta",
        reference = "config/{a_or_b}reference.gbk"
    output:
        hybrid_alignment = build_dir + "/{a_or_b}/{build_name}/hybrid_alignment.fasta"
    params:
        gene = lambda w: w.build_name
    shell:
        """
        python scripts/align_for_tree.py \
            --realign {input.G_alignment} \
            --original {input.original} \
            --reference {input.reference} \
            --output {output.hybrid_alignment} \
            --gene {params.gene}
        """

def get_alignment(w):
    if w.build_name == "genome":
        return rules.hybrid_align.output.hybrid_alignment
    else:
        return build_dir + f"/{w.a_or_b}/{w.build_name}/{w.build_name}_aligned.fasta"

rule tree:
    message: "Building tree"
    input:
        alignment = get_alignment
    output:
        tree = build_dir + "/{a_or_b}/{build_name}/tree_raw.nwk"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads}
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
        alignment =get_alignment,
        metadata = rules.filter.input.metadata
    output:
        tree = build_dir + "/{a_or_b}/{build_name}/tree.nwk",
        node_data = build_dir + "/{a_or_b}/{build_name}/branch_lengths.json"
    params:
        coalescent = config["refine"]["coalescent"],
        clock_filter_iqd = config["refine"]["clock_filter_iqd"],
        date_inference = config["refine"]["date_inference"],
        strain_id=config["strain_id_field"],
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --timetree \
            --use-fft \
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
        alignment = get_alignment
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
        reference = build_dir + "/{a_or_b}/{build_name}/{build_name}_reference.gbk",
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/aa_muts.json"
    params:
    	alignment_file_mask = build_dir + "/{a_or_b}/{build_name}/aligned_%GENE.fasta"
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
        columns = config["traits"]["columns"],
        strain_id=config["strain_id_field"],
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """