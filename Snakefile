import pandas as pd

configfile: "config/configfile.yaml"

include: "workflow/snakemake_rules/config.smk"

wildcard_constraints:
    a_or_b=r"a|b",
    build_name=r"genome|G|F",
    resolution=r"all-time|6y|3y",

build_dir = "results"
auspice_dir = "auspice"

distance_map_config = pd.read_table("config/distance_maps.tsv")


rule all:
    input:
        expand("auspice/rsv_{subtype}_{build}_{resolution}.json",
               subtype = config["subtypes"],
               build = config["builds_to_run"],
               resolution = config["resolutions_to_run"]),
        expand("auspice/rsv_{subtype}_{build}_{resolution}_tip-frequencies.json",
               subtype = config["subtypes"],
               build = config["builds_to_run"],
               resolution = config["resolutions_to_run"]),

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
