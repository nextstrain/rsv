
"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/run_config.yaml
"""
import os
import sys
import yaml
from textwrap import dedent


RUN_CONFIG = f"results/run_config.yaml"


def main():
    validate_config()
    write_config(RUN_CONFIG)


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

            See 'subsample' in the default config (config/configfile.yaml) for
            an example of how to specify filtering and subsampling parameters."""))
        exit(1)


def write_config(path):
    """
    Write Snakemake's 'config' variable to a file.

    This is used for the subsample rule and is generally useful for debugging.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, 'w') as f:
        yaml.dump(config, f, sort_keys=False)

    print(f"Saved current run config to {path!r}.", file=sys.stderr)


def indented_list(xs, prefix):
    return f"\n{prefix}".join(xs)


main()
