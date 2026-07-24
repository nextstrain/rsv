import re
import shlex

import pandas as pd
from snakemake.utils import min_version

# Minimum Snakemake version needed for the storage plugins used in remote_files.smk
min_version("8.0.0")

configfile: "config/configfile.yaml"

include: "shared/vendored/snakemake/config.smk"
include: "workflow/snakemake_rules/config.smk"

builds = config["builds"]

wildcard_constraints:
    build="|".join(builds),
    build_with_underscores="|".join(b.replace("/", "_") for b in builds),
    a_or_b="|".join(set(get_subtype(b) for b in builds)),
    gene="|".join(set(get_gene(b) for b in builds)),
    resolution="|".join(set(get_resolution(b) for b in builds)),


results_dir = "results"
auspice_dir = "auspice"

distance_map_config = pd.read_table("config/distance_maps.tsv")


rule all:
    input:
        [f"auspice/rsv_{build.replace('/', '_')}.json" for build in builds],
        [f"auspice/rsv_{build.replace('/', '_')}_tip-frequencies.json" for build in builds],


# remote_files.smk must be before merge_inputs.smk
include: "shared/vendored/snakemake/remote_files.smk"
include: "workflow/snakemake_rules/merge_inputs.smk"

include: "workflow/snakemake_rules/core.smk"
include: "workflow/snakemake_rules/export.smk"
include: "workflow/snakemake_rules/glycosylation.smk"
include: "workflow/snakemake_rules/clades.smk"


if "deploy_url" in config:

    include: "workflow/snakemake_rules/nextstrain_automation.smk"


rule clean:
    params:
        targets=["auspice", "results"],
    shell:
        """
        rm -rf {params.targets}
        """


rule clobber:
    params:
        targets=["data", "auspice", "results"],
    shell:
        """
        rm -rf {params.targets}
        rm config/clades*tsv
        """

if "custom_rules" in config:
    for rule_file in config["custom_rules"]:
        include: rule_file
