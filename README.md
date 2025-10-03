# nextstrain.org/rsv

This is the Nextstrain build for respiratory syncytial virus (RSV).
Output from this build is available at <https://nextstrain.org/rsv>.

## Data use

The default builds start from the public Nextstrain data that have been preprocessed
and cleaned from [Pathoplexus][] that includes [RESTRICTED data][].
Please see [Pathoplexus data use terms][].

We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work. Please note that although data generators have generously shared data in an open fashion, that does not mean there should be free license to publish on this data. Data generators should be cited where possible and collaborations should be sought in some circumstances.

## Software requirements

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html)
for Nextstrain's suite of software tools.

## Usage

If you're unfamiliar with Nextstrain builds, you may want to follow our
[Running a Pathogen Workflow guide][] first and then come back here.

### With `nextstrain build`

If you don't have a local copy of the rsv repository, use Git to download it

    git clone https://github.com/nextstrain/rsv.git

Otherwise, update your local copy of the workflow with:

    cd rsv
    git pull --ff-only origin master

Run the phylogenetic workflow workflow with

    nextstrain build .

The workflow's intermediate files will be output to `results/` and the final
outputs will be in `auspice/`.

Once you've run the build, you can view the results with:

    nextstrain view .


## Configuration

The default configuration is in [`config/configfile.yaml`](./config/configfile.yaml).
The workflow is contained in the [Snakefile](Snakefile) with included
[rules](workflow/snakemake_rules/). Each rule specifies its file inputs and outputs
and pulls its parameters from the config. There is little redirection and each
rule should be able to be reasoned with on its own.

### Default input data

The default builds start from the public Nextstrain data that have been preprocessed
and cleaned from [Pathoplexus][] that includes [RESTRICTED data][].
The [default auspice_config.json](./config/auspice_config.json) includes the
`metadata_columns` "PPX_accession", "INSDC_accession", and "restrictedUntil"
to ensure the builds adhere to the [Pathoplexus data use terms][].

```yaml
subtypes: ['a', 'b']
inputs:
  - name: ppx_with_restricted
    metadata: "https://data.nextstrain.org/files/workflows/rsv/{a_or_b}/metadata_with_restricted.tsv.gz"
    sequences: "https://data.nextstrain.org/files/workflows/rsv/{a_or_b}/sequences_with_restricted.fasta.xz"
```

Note the inputs require the `{a_or_b}` expandable field, to be replaced by the
config parameter `subtypes` values.

### Adding your own data

If you want to add your own data to the default input, specify your inputs with
the `additional_inputs` config parameter. For example, this repo has a small set
of example data that could be added to the default inputs via:

```yaml
additional_inputs:
  - name: example-data
    metadata: example_data/{a_or_b}/metadata.tsv
    sequences: example_data/{a_or_b}/sequences.fasta
```

Note that the additional inputs also require the `{a_or_b}` expandable field.
If you only have data for a single subtype, then you can do so with

```yaml
serotypes: ["a"]
additional_inputs:
  - name: private
    metadata: private/a/metadata.tsv
    sequences: private/a/sequences.fasta
```

If you want to run the builds _without_ the default data and only use your own
data, you can do so by specifying the `inputs` parameter.

```yaml
inputs:
  - name: example-data
    metadata: example_data/{a_or_b}/metadata.tsv
    sequences: example_data/{a_or_b}/sequences.fasta
```

### Using locally ingested data

Run the ingest pipeline locally with

```sh
nextstrain build ingest
```

Then you can point the phylogenetic workflow to run from the produced results with

```yaml
inputs:
    - name: local_ingest
      metadata: "ingest/data/{a_or_b}/metadata.tsv"
      sequences: "ingest/data/{a_or_b}/sequences.fasta"
```

## `shared/vendored`

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage copies of shared scripts in [`shared/vendored`](./shared/vendored), from [nextstrain/shared](https://github.com/nextstrain/shared).
See [shared/vendored/README.md](./shared/vendored/README.md#vendoring) for instructions on how to update the vendored scripts.


## Deep mutational scanning data
The trees are annotated with escape from key monoclonal antibodies as assessed by the effects of mutations measured in deep mutational scanning by [Simonich et al]() **ADD CITATION**.

## Update example data

[Example data](./example_data/) is used by [CI](https://github.com/nextstrain/rsv/actions/workflows/ci.yaml). It can also be used as a small subset of real-world data.

Example data should be updated every time metadata schema is changed. To update, run:

```sh
nextstrain build --docker . update_example_data --configfile config/chores.yaml -F
```


## Sending data to the `nextclade_data` repo

From within the destination directory, run
```
rsync -a <path-to>/rsv/nextclade/datasets/ .
```

[Installing Nextstrain guide]: https://docs.nextstrain.org/en/latest/install.html
[Pathoplexus]: https://pathoplexus.org
[Pathoplexus data use terms]: https://pathoplexus.org/about/terms-of-use/data-use-terms
[RESTRICTED data]: https://pathoplexus.org/about/terms-of-use/restricted-data
[Running a Pathogen Workflow guide]: https://docs.nextstrain.org/en/latest/tutorials/running-a-workflow.html
