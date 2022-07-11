configfile: "config/configfile.yaml"

build_dir = 'results'
auspice_dir = 'auspice'

rule all:
    input:
        results = "auspice/rsv.json"

include: "workflow/snakemake_rules/core.smk"

if config['gene'] == 'G':
    include: "workflow/snakemake_rules/glycosylation.smk"
else:
    include: "workflow/snakemake_rules/export.smk"
