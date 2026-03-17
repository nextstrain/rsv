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
    if "custom_subsample_matrix" in config:
        matrix = config["custom_subsample_matrix"]
    else:
        matrix = config["subsample_matrix"]
    default_options = matrix.get("defaults", {})

    # Handle `keep_all` options from inputs
    # input_sources is from merge_inputs.smk
    keep_all_samples = {}
    for input_name, input_info in input_sources.items():
        if input_info.get("keep_all", False):
            keep_all_samples[f"_workflow_keep_{input_name}"] = {
                "exclude_all": True,
                "include_where": f"input_{input_name}=1",
            }

    for a_or_b in config["subtypes"]:
        subtype = matrix["subtypes"][a_or_b]
        subtype_options = {k: v for k, v in subtype.items() if k != "sample_overrides"}
        subtype_overrides = subtype.get("sample_overrides", {})

        for build_name in config["builds_to_run"]:
            build = matrix["builds"][build_name]
            build_options = {k: v for k, v in build.items() if k != "sample_overrides"}
            build_overrides = build.get("sample_overrides", {})

            for resolution in config["resolutions_to_run"]:
                subsample_config = {"samples": {}}

                # Samples are defined at the resolution level
                samples = matrix["resolutions"][resolution]["samples"]

                if conflicting_sample_names := samples.keys() & keep_all_samples.keys():
                    print(dedent(f"""\
                        ERROR: The following sample names conflict with generated sample names:
                            {conflicting_sample_names!r}"""))
                    exit(1)

                # Merge options
                for sample_name, sample_options in samples.items():
                    # Check overrides
                    sample_subtype_overrides = subtype_overrides.get(sample_name, {})
                    sample_build_overrides = build_overrides.get(sample_name, {})
                    if conflicting_override_keys := sample_subtype_overrides.keys() & sample_build_overrides.keys():
                        print(dedent(f"""\
                            ERROR: sample_overrides for sample {sample_name!r} in subtype {a_or_b!r}
                            and build {build_name!r} both set: {conflicting_override_keys}"""))
                        exit(1)

                    # All option sources ordered by precedence (last one wins)
                    options_source = [
                        default_options,
                        subtype_options,
                        build_options,
                        sample_options,
                        sample_subtype_overrides,
                        sample_build_overrides,
                    ]

                    # Merge sources in order, overriding any previously existing options.
                    # Custom merging is implemented below for specific options.
                    merged_sample_options = {}
                    for options in options_source:
                        merged_sample_options.update(options)

                    # Custom merge for 'query'
                    queries = []
                    for options in options_source:
                        if value := options.get("query"):
                            queries.append(value)
                    if queries:
                        merged_sample_options["query"] = " & ".join(f"({q})" for q in queries)

                    # Custom merge for 'exclude_where'
                    exclude_wheres = []
                    for options in options_source:
                        if value := options.get("exclude_where"):
                            if isinstance(value, str):
                                exclude_wheres.append(value)
                            elif isinstance(value, list):
                                exclude_wheres.extend(value)
                    if exclude_wheres:
                        merged_sample_options["exclude_where"] = exclude_wheres

                    subsample_config["samples"][sample_name] = merged_sample_options

                subsample_config["samples"].update(keep_all_samples)

                path = build_dir + f"/{a_or_b}/{build_name}/{resolution}/subsample_config.yaml"
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w") as f:
                    yaml.dump(subsample_config, f, sort_keys=False, Dumper=NoAliasDumper)
                print(f"Saved subsampling config to {path!r}.", file=sys.stderr)


main()
