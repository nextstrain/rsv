# nextstrain.org/rsv

This is the Nextstrain build for respiratory syncytial virus (RSV). Output from this build is available at nextstrain.org/rsv.

## Input Data

Input metadata and sequences for RSV-A and RSV-B are available via <https://data.nextstrain.org>

- [RSV-A sequences](https://data.nextstrain.org/files/workflows/rsv/a/sequences.fasta.xz)
- [RSV-A metadata](https://data.nextstrain.org/files/workflows/rsv/a/metadata.tsv.gz)

- [RSV-B sequences](https://data.nextstrain.org/files/workflows/rsv/b/sequences.fasta.xz)
- [RSV-B metadata](https://data.nextstrain.org/files/workflows/rsv/b/metadata.tsv.gz)

These data are generously shared by labs around the world and deposited in NCBI genbank by the authors.
Please contact these labs first if you plan to publish using these data.
RSV sequences and metadata can be downloaded in the `/ingest` folder using
`nextstrain build --cpus 1 ingest` or `nextstrain build --cpus 1 .` if running directly from the `/ingest` directory.

The ingest pipeline is based on the Nextstrain mpox ingest workflow (<https://github.com/nextstrain/mpox/tree/master/ingest>).
Running the ingest pipeline produces `ingest/data/{a,b}/metadata.tsv` and `ingest/data/{a,b}/sequences.fasta`.

### `ingest/vendored`

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage copies of ingest scripts in [`ingest/vendored`](./ingest/vendored), from [nextstrain/ingest](https://github.com/nextstrain/ingest). To pull new changes from the central ingest repository, first install `git subrepo`, then run:

See [ingest/vendored/README.md](./ingest/vendored/README.md#vendoring) for instructions on how to update the vendored scripts.

## Run Analysis Pipeline

The workflow produces whole genome and G gene trees for RSV-A and RSV-B.
To run the workflow, use `snakemake -j4 -p --configfile config/configfile.yaml` and `nextstrain view auspice` to visualise results.

## Installation

Follow the standard [installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Data use

We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work. Please note that although data generators have generously shared data in an open fashion, that does not mean there should be free license to publish on this data. Data generators should be cited where possible and collaborations should be sought in some circumstances.

## Update example data

[Example data](./example_data/) is used by [CI](https://github.com/nextstrain/rsv/actions/workflows/ci.yaml). It can also be used as a small subset of real-world data.

Example data should be updated every time metadata schema is changed. To update, run:

```sh
nextstrain build --docker . update_example_data -F
```


## Sending data to the `nextclade_data` repo

From within the destination directory, run
```
rsync -a <path-to>/rsv/nextclade/datasets/ .
```
