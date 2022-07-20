# nextstrain.org/rsv


This is the Nextstrain build for respiratory syncytial virus (RSV). Output from this build is available at nextstrain.org/rsv.

### Input Data

Input metadata and sequences for RSV-A and RSV-B are available on data.nextstrain.org

* [RSV-A sequences](https://data.nextstrain.org/files/workflows/rsv/a/sequences.fasta.xz)
* [RSV-A metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)

* [RSV-B sequences](https://data.nextstrain.org/files/workflows/rsv/b/sequences.fasta.xz)
* [RSV-B metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)


These data are generously shared by labs around the world and deposited in NCBI genbank by the authors.
Please contact these labs first if you plan to publish using these data.
RSV sequences and metadata can be downloaded in the ```/ingest``` folder using
```nextstrain build --cpus 1 ingest``` or ```nextstrain build --cpus 1 .``` if running directly from the ```/ingest``` directory.

The ingest pipeline is based on the Nextstrain Monkeypox ingest (nextstrain.org/monkeypox/ingest).
Running the ingest pipeline produces ```ingest/data/{a and b}/metadata.tsv``` and ```ingest/data/{a and b}/sequences.fasta```.


## Run Analysis Pipeline

The workflow produces whole genome and G gene trees for RSV-A and RSV-B.
To run the workflow, use ``snakemake -j4 -p --configfile config/configfile.yaml ` and ```nextstrain view auspice``` to visualise results.

## Installation

Follow the standard [installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Data Use

We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work. Please note that although data generators have generously shared data in an open fashion, that does not mean there should be free license to publish on this data. Data generators should be cited where possible and collaborations should be sought in some circumstances.
