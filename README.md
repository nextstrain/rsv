# nextstrain rsv

### Input Data

RSV sequences and metadata can be downloaded in the ```/ingest``` folder using the command
```nextstrain build --cpus 1 ingest``` or ```nextstrain build --cpus 1 .``` if running directly from the ```/ingest``` directory.

This will produce ingest/data/metadata.tsv and ingest/data/sequences.fasta.

##

By default, the workflow uses whole genome sequences. To run for individual genes, change 'everything' to gene of interest (eg. 'G') in ```/config/configfile.yaml```.

To run the workflow, use ```snakemake  --cores 1``` and ```nextstrain view auspice``` to visualise results.

##


The workflow can be used for either for either rsv-A or rsv-B. To specify which, add 'a' or 'b' in the config file for rsv.
