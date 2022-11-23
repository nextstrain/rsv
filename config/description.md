We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work. Please note that although data generators have generously shared data in an open fashion, that does not mean there should be free license to publish on this data. Data generators should be cited where possible and collaborations should be sought in some circumstances. Please try to avoid scooping someone else's work. Reach out if uncertain.


We maintain three views of human respiratory syncytial virus evolution for each RSV subtype:

The first is ['rsv/a/genome'](https://nextstrain.org/staging/rsv/a/genome) and ['rsv/b/genome'](https://nextstrain.org/staging/rsv/b/genome), which show evolution of full genome sequences.

The second is ['rsv/a/G'](https://nextstrain.org/staging/rsv/a/G) and ['rsv/b/G'](https://nextstrain.org/staging/rsv/b/G), which show evolution of only the G gene, which is highly variable and for which more data points are available than for the full-genome.

The second is ['rsv/a/F'](https://nextstrain.org/staging/rsv/a/F) and ['rsv/b/F'](https://nextstrain.org/staging/rsv/b/G), which show evolution of only the F gene. The F gene builds currently don't contain any clade annotations.


#### Analysis
Our bioinformatic processing workflow can be found at [github.com/nextstrain/rsv](https://github.com/nextstrain/rsv) and includes:

- sequence alignment by a combination of [nextalign](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextalign-cli.html) and [MAFFT](https://mafft.cbrc.jp/alignment/software/).
- phylogenetic reconstruction using [IQTREE](http://www.iqtree.org/)
- ancestral state reconstruction and temporal inference using [TreeTime](https://github.com/neherlab/treetime)
- clade assignment via clade definitions defined here for
    [RSV-A/genome](https://github.com/nextstrain/rsv/blob/master/config/clades_genome_a.tsv),
    [RSV-B/genome](https://github.com/nextstrain/rsv/blob/master/config/clades_genome_b.tsv),
    [RSV-A/G gene](https://github.com/nextstrain/rsv/blob/master/config/clades_G_a.tsv), and
    [RSV-B/G gene](https://github.com/nextstrain/rsv/blob/master/config/clades_G_b.tsv) to label RSV clades based on the entire genome and for just the G gene.
    These clade definitions are based on the proposed nomenclatures by [Goya et al](https://onlinelibrary.wiley.com/doi/abs/10.1111/irv.12715) and [Ramaekers et al](https://doi.org/10.1093/ve/veaa052).

#### Underlying data
We curate sequence data and metadata from [NCBI Virus (RSV)](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Human%20orthopneumovirus,%20taxid:11250), [NCBI Virus (RSV-A)](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Human%20respiratory%20syncytial%20virus%20A,%20taxid:208893), [NCBI Virus (RSV-B)](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Human%20respiratory%20syncytial%20virus%20B,%20taxid:208895) as starting point for these analyses. Curated sequences and metadata are available as flat files at:

- [RSV-A sequences](https://data.nextstrain.org/files/workflows/rsv/a/sequences.fasta.xz)
- [RSV-A metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)

- [RSV-B sequences](https://data.nextstrain.org/files/workflows/rsv/b/sequences.fasta.xz)
- [RSV-B metadata](https://data.nextstrain.org/files/workflows/rsv/b/metadata.tsv.gz)



