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

### Use locally ingested data

Once you have run the ingest pipeline locally you can copy the files into the top-level `data` directory so that the main phylo workflow uses these files rather than downloading from s3:

```sh
mkdir -p data/{a,b}
for i in ingest/data/*/{metadata.tsv,sequences.fasta}; do cp $i ${i#ingest/}; done
```


### `ingest/vendored`

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage copies of ingest scripts in [`ingest/vendored`](./ingest/vendored), from [nextstrain/ingest](https://github.com/nextstrain/ingest). To pull new changes from the central ingest repository, first install `git subrepo`, then run:

See [ingest/vendored/README.md](./ingest/vendored/README.md#vendoring) for instructions on how to update the vendored scripts.

## Run Analysis Pipeline

The workflow produces whole genome and G gene trees for RSV-A and RSV-B.
To run the workflow, use `snakemake -j4 -p --configfile config/configfile.yaml` and `nextstrain view auspice` to visualise results.

## Installation

Follow the standard [installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Configuration

The pipeline configuration is managed using [CUE](https://cuelang.org), which provides better maintainability than YAML anchors and aliases. The source file is `config/configfile.cue`, and it generates `config/configfile.yaml` that Snakemake uses.

### For Contributors (modifying config)

If you need to modify the configuration:

1. **Install CUE** (one time setup):
   ```sh
   micromamba create -n cue -c conda-forge cue pyyaml
   ```

2. **Edit the CUE file** (not the YAML):
   ```sh
   # Edit config/configfile.cue with your changes
   ```

3. **Regenerate YAML**:
   ```sh
   make config/configfile.yaml
   # or manually: micromamba run -n cue cue export config/configfile.cue --out yaml --outfile config/configfile.yaml
   ```

4. **Validate** (optional):
   ```sh
   make validate-config
   ```

See `make help` for other config management targets.

### For Users (running the pipeline)

No CUE installation required! Just use the committed `config/configfile.yaml` as normal. The workflow uses the YAML file directly.

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
