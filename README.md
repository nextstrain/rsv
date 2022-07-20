# nextstrain.org/rsv

##

This is the Nextstrain build for respiratory syncytial virus (RSV). Output from this build is visible at nextstrain.org/rsv.

### Input Data

Input metadata and sequences for RSV-A and RSV-B are available on data.nextstrain.org

* [RSV-A sequences](https://data.nextstrain.org/files/workflows/rsv/a/sequences.fasta.xz)
* [RSV-A metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)

* [RSV-B sequences](https://data.nextstrain.org/files/workflows/rsv/b/sequences.fasta.xz)
* [RSV-B metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)


This data is generously shared by labs around the world. Please contact these labs first if you plan to publish using this data.
RSV sequences and metadata can be downloaded in the ```/ingest``` folder using
```nextstrain build --cpus 1 ingest``` or ```nextstrain build --cpus 1 .``` if running directly from the ```/ingest``` directory.

The ingest pipeline is based on the Nextstrain Monkeypox ingest (nextstrain.org/monkeypox/ingest). 
This will produce ingest/data/{a and b}/metadata.tsv and ingest/data/{a and b}/sequences.fasta.


###Run analysis pipeline

##

The workflow produces whole genome and G gene trees for RSV-A and RSV-B. 
To run the workflow, use ``snakemake -j4 -p --configfile config/configfile.yaml ` and ```nextstrain view auspice``` to visualise results.

##


