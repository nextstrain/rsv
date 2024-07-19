We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work. Please note that although data generators have generously shared data in an open fashion, that does not mean there should be free license to publish on this data. Data generators should be cited where possible and collaborations should be sought in some circumstances. Please try to avoid scooping someone else's work. Reach out if uncertain.

We maintain three views of human respiratory syncytial virus evolution for each RSV subtype:

The first is ['rsv/a/genome'](https://nextstrain.org/staging/rsv/a/genome) and ['rsv/b/genome'](https://nextstrain.org/staging/rsv/b/genome), which show evolution of full genome sequences.

The second is ['rsv/a/G'](https://nextstrain.org/staging/rsv/a/G) and ['rsv/b/G'](https://nextstrain.org/staging/rsv/b/G), which show evolution of only the G gene, which is highly variable and for which more data points are available than for the full-genome.

The second is ['rsv/a/F'](https://nextstrain.org/staging/rsv/a/F) and ['rsv/b/F'](https://nextstrain.org/staging/rsv/b/G), which show evolution of only the F gene. The F gene builds currently don't contain any clade annotations.

#### Analysis

Our bioinformatic processing workflow can be found at [github.com/nextstrain/rsv](https://github.com/nextstrain/rsv) and includes:

- sequence alignment by a combination of [Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html) and [MAFFT](https://mafft.cbrc.jp/alignment/software/).
- phylogenetic reconstruction using [IQTREE](http://www.iqtree.org/)
- ancestral state reconstruction and temporal inference using [TreeTime](https://github.com/neherlab/treetime)
- clade assignment via clade definitions defined here:
  [RSV-A](https://raw.githubusercontent.com/rsv-lineages/lineage-designation-A/main/.auto-generated/lineage.tsv)
  [RSV-B](https://raw.githubusercontent.com/rsv-lineages/lineage-designation-A/main/.auto-generated/lineage.tsv)
  These clade definitions are based on the not-yet-published nomenclature of the RSV Genotyping Consensus Consortium.

#### Underlying data

We curate sequence data and metadata from the [NCBI Datasets command line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)
as starting point for these analyses. See our [ingest configuration file](https://github.com/nextstrain/rsv/blob/master/ingest/config/config.yaml)
for the NCBI Taxonomy IDs used to fetch the virus genomes.
Curated sequences and metadata are available as flat files at:

- [RSV-A sequences](https://data.nextstrain.org/files/workflows/rsv/a/sequences.fasta.xz)
- [RSV-A metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)

- [RSV-B sequences](https://data.nextstrain.org/files/workflows/rsv/b/sequences.fasta.xz)
- [RSV-B metadata](https://data.nextstrain.org/files/workflows/rsv/b/metadata.tsv.gz)
