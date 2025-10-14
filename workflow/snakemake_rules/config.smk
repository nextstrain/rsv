"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/config_raw.yaml
    results/config_processed.yaml
"""
from copy import deepcopy


def main():
    global config

    # Dump the unmodified config to a YAML file.
    write_config("results/config_raw.yaml")

    process_subsample_config()

    # Move user config to its own section.
    user_config = deepcopy(config)
    config = {
        "WORKFLOW_CONFIG": {},
        "USER_CONFIG": user_config,
    }

    restructure_and_resolve_wildcards()
    # TODO: resolve filepaths

    # 3. Write the modified config to a file.
    #    This file is used by augur subsample.
    write_config("results/config_processed.yaml")


def restructure_and_resolve_wildcards():
    """Populate WORKFLOW_CONFIG into a standardized format for internal usage.

    Format is:

    WORKFLOW_CONFIG:
      <build>/<resolution>:
        <rule>:
          <config>

    This workflow doesn't use wildcards in config values (i.e. string
    placeholders), but if it ever does in the future, they should be resolved in
    this function.
    """

    global config

    user_config = config["USER_CONFIG"]

    # FIXME: add 'a_or_b' as a dimension?
    for build_name in user_config['builds_to_run']:
        for resolution in user_config['resolutions_to_run']:
            key = f"{build_name}/{resolution}"

            config["WORKFLOW_CONFIG"][key] = {}
            analysis_config = config["WORKFLOW_CONFIG"][key]

            analysis_config["subsample"] = {
                "strain_id_field": user_config["strain_id_field"],
                "config": user_config["subsample"][build_name][resolution]
            }

            analysis_config["refine"] = {
                **user_config["refine"],
                "strain_id_field": user_config["strain_id_field"],
            }

            analysis_config["ancestral"] = {
                **user_config["ancestral"],
                "genes": user_config["cds"][build_name],
            }

            analysis_config["frequencies"] = {
                **user_config["frequencies"]["resolutions"][resolution],
            }

            analysis_config["traits"] = {
                **user_config["traits"],
                "strain_id_field": user_config["strain_id_field"],
            }

            analysis_config["export"] = {
                "strain_id_field": user_config["strain_id_field"],
                "display_strain_field": user_config["display_strain_field"],
                "auspice_config": user_config["files"]["auspice_config"],
                "auspice_config_additional_colorings": user_config["files"]["auspice_config_additional_colorings"],
                "description": user_config["description"],
            }

            analysis_config["final_strain_name"] = {
                "strain_id_field": user_config["strain_id_field"],
                "display_strain_field": user_config["display_strain_field"],
            }

            analysis_config["genome_align"] = {
                "genes": user_config["cds"][build_name],
            }

            analysis_config["distances"] = {
                "genes": user_config["cds"][build_name],
            }


def get_config_value(*keys):
    """
    Return a callable to retrieve a config value from WORKFLOW_CONFIG given wildcards.
    """
    def _get(wildcards):
        build_name = wildcards.build_name if hasattr(wildcards, 'build_name') else wildcards["build_name"]
        resolution = wildcards.resolution if hasattr(wildcards, 'resolution') else wildcards["resolution"]
        key = f"{build_name}/{resolution}"
        result = config["WORKFLOW_CONFIG"][key]
        for k in keys:
            result = result[k]
        return result
    return _get


def get_subsample_config_section():
    def _get(wildcards):
        build_name = wildcards.build_name if hasattr(wildcards, 'build_name') else wildcards["build_name"]
        resolution = wildcards.resolution if hasattr(wildcards, 'resolution') else wildcards["resolution"]
        key = f"{build_name}/{resolution}"
        return ["WORKFLOW_CONFIG", key, "subsample", "config"]
    return _get


main()
