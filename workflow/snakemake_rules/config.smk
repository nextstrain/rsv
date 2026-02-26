"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/run_config.yaml
"""
from textwrap import dedent


def main():
    validate_config()
    write_subsample_config()
    write_config("results/run_config.yaml")


def validate_config():
    """
    Validate the config.

    This could be improved with a schema definition file, but for now it serves
    to provide useful error messages for common user errors and effects of
    breaking changes.
    """

    # Check for deprecated 'filter' key
    if "filter" in config:
        print(dedent(f"""\
            ERROR: The 'filter' configuration key is no longer supported.

            See 'subsample_matrix' in the default config (config/configfile.yaml) for
            an example of how to specify filtering and subsampling parameters."""))
        exit(1)


def write_subsample_config():
    matrix = config["subsample_matrix"]
    default_options = matrix.get("defaults", {})

    for a_or_b in config["subtypes"]:
        for build_name in config["builds_to_run"]:
            for resolution in config["resolutions_to_run"]:
                build = matrix["builds"][build_name]
                build_options = {k: v for k, v in build.items() if k != "sample_overrides"}
                build_sample_overrides = build.get("sample_overrides", {})
                subsample_config = {
                    "samples": {
                        sample_name: {
                            **default_options,
                            **build_options,
                            **sample_options,
                            **build_sample_overrides.get(sample_name, {}),
                        } for sample_name, sample_options in matrix["resolutions"][resolution]["samples"].items()
                    }
                }
                path = build_dir + f"/{a_or_b}/{build_name}/{resolution}/subsample_config.yaml"
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w") as f:
                    yaml.dump(subsample_config, f, sort_keys=False, Dumper=NoAliasDumper)
                print(f"Saved subsampling config to {path!r}.", file=sys.stderr)


main()
