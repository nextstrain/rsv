We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work. Please note that although data generators have generously shared data in an open fashion, that does not mean there should be free license to publish on this data. Data generators should be cited where possible and collaborations should be sought in some circumstances. Please try to avoid scooping someone else's work. Reach out if uncertain.


We maintain two views of human respiratory syncytial virus evolution for each RSV subtype:

The first is ['rsv/a/genome'](https://nextstrain.org/staging/rsv/a/genome) and ['rsv/b/genome'](https://nextstrain.org/staging/rsv/b/genome), which show evolution of full genome sequences.

The second is ['rsv/a/G'](https://nextstrain.org/staging/rsv/a/G) and ['rsv/b/G'](https://nextstrain.org/staging/rsv/b/G), which show evolution of only the G gene, which is highly variable.



#### Analysis
Our bioinformatic processing workflow can be found at [github.com/nextstrain/rsv](https://github.com/nextstrain/rsv) and includes:

- sequence alignment by [nextalign](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextalign-cli.html)
- phylogenetic reconstruction using [IQTREE](http://www.iqtree.org/)
- ancestral state reconstruction and temporal inference using [TreeTime](https://github.com/neherlab/treetime)
- clade assignment via clade definitions defined here for [genome, RSV-A](https://github.com/nextstrain/rsv/blob/master/config/clades_genome_a.tsv), [genome, RSV-B](https://github.com/nextstrain/rsv/blob/master/config/clades_genome_b.tsv), [G gene, RSV-A](https://github.com/nextstrain/rsv/blob/master/config/clades_G_a.tsv), [G gene, RSV-B](https://github.com/nextstrain/rsv/blob/master/config/clades_G_b.tsv),to label RSV clades based on the entire genome and for just the G gene.

#### Underlying data
We curate sequence data and metadata from [NCBI Virus (RSV)](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Human%20orthopneumovirus,%20taxid:11250), [NCBI Virus (RSV-A)](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Human%20respiratory%20syncytial%20virus%20A,%20taxid:208893), [NCBI Virus (RSV-B)](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Human%20respiratory%20syncytial%20virus%20B,%20taxid:208895) as starting point for these analyses. Curated sequences and metadata are available as flat files at:

- [RSV-A sequences](https://data.nextstrain.org/files/workflows/rsv/a/sequences.fasta.xz)
- [RSV-A metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)

- [RSV-B sequences](https://data.nextstrain.org/files/workflows/rsv/b/sequences.fasta.xz)
- [RSV-B metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)



