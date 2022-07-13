configfile: "config/configfile.yaml"

build_dir = 'results'
auspice_dir = 'auspice'

if config['rsv'] == 'a':
    reference = 'areference'
if config['rsv'] == 'b':
    reference = 'breference'

rule all:
    input:
        results = "auspice/rsv.json"

include: "workflow/snakemake_rules/core.smk"

include: "workflow/snakemake_rules/export.smk"

include: "workflow/snakemake_rules/download.smk"

if config['gene'] == 'G':
    include: "workflow/snakemake_rules/glycosylation.smk"

