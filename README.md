# Nextstrain RSV

### Input Data

RSV sequences and metadata can be downloaded in the ```/ingest``` folder using the command
```nextstrain build --cpus 1 ingest``` or ```nextstrain build --cpus 1 .``` if running it from the ```/ingest``` directory.
This will produce ingest/data/metadata.tsv and ingest/data/sequences.fasta.

To run the pipeline, use ```snakemake  --cores 1``` and ```nextstrain view auspice``` to visualise results.

