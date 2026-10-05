"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/run_config.yaml
"""
import sys
from pathlib import Path

def main():
    dump_and_validate(
        "results/run_config.yaml",
        Path(workflow.basedir) / "config.schema.yaml"
    )
    write_subsample_config()


def write_subsample_config():
    # TODO: Support custom build names in the workflow and infer from
    # config["builds"].
    for a_or_b in ["a", "b"]:
        for build_name in ["genome", "G", "F", "F-antibody-escape"]:
            for resolution in ["all-time", "6y", "3y"]:
                build = f"{a_or_b}/{build_name}/{resolution}"
                if "custom_subsample" in config:
                    section = ["custom_subsample", build]
                else:
                    section = ["subsample", build]
                write_config(f"results/{build}/subsample_config.yaml", section=section)


try:
    main()
except InvalidConfigError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    exit(1)
