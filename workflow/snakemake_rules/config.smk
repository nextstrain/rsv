
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
    write_config(RUN_CONFIG)


def write_config(path):
    """
    Write Snakemake's 'config' variable to a file.

    This is used for the subsample rule and is generally useful for debugging.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)

    with open(path, 'w') as f:
        yaml.dump(config, f, sort_keys=False)

    print(f"Saved current run config to {path!r}.", file=sys.stderr)


main()
