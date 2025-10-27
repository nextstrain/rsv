import re
import shlex

import pandas as pd

configfile: "config/configfile.yaml"


wildcard_constraints:
    a_or_b=r"a|b",
    build_name="|".join(config.get("builds_to_run", ["genome"])),
    resolution="|".join(config.get("resolutions_to_run", ["all-time"])),


build_dir = "results"
auspice_dir = "auspice"

distance_map_config = pd.read_table("config/distance_maps.tsv")


rule all:
    input:
        expand("auspice/rsv_{subtype}_{build}_{resolution}.json",
               subtype = config.get("subtypes",['a']),
               build = config.get("builds_to_run", ['genome']),
               resolution = config.get("resolutions_to_run", ["all-time"])),
        expand("auspice/rsv_{subtype}_{build}_{resolution}_tip-frequencies.json",
               subtype = config.get("subtypes",['a']),
               build = config.get("builds_to_run", ['genome']),
               resolution = config.get("resolutions_to_run", ["all-time"])),
        expand(build_dir + "/{subtype}/{build}/{resolution}/pre_subsample/metadata_with_scores.tsv",
               subtype=config.get("subtypes",['a']),
               build=["F-antibody-escape"],
               resolution = config.get("resolutions_to_run", ["all-time"])),


include: "workflow/snakemake_rules/chores.smk"


include: "workflow/snakemake_rules/core.smk"
include: "workflow/snakemake_rules/export.smk"
include: "workflow/snakemake_rules/download.smk"
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
