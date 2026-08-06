"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/run_config.yaml
"""
import sys
from augur.validate import load_json_schema_locally, validate_json, ValidateError
from pathlib import Path

def main():
    dump_and_validate(
        "results/run_config.yaml",
        Path(workflow.basedir) / "config.schema.yaml"
    )
    write_subsample_config()


# TODO: move this to nextstrain/shared
# Copied from measles with some cleanups
# <https://github.com/nextstrain/measles/blob/60634e2f/phylogenetic/rules/config.smk#L36-L57>
def dump_and_validate(dump_path, schema_path):
    """
    Write Snakemake's 'config' variable to a file, then validate it against the
    schema. Do both in the same function so that the validation output can
    easily reference the path of the dumped config for inspection.
    """
    global config

    write_config(dump_path)

    if "custom_rules" in config:
        print("WARNING: Skipping config schema validation because custom rules are defined.", file=sys.stderr)
        return

    try:
        validator = load_json_schema_locally(schema_path)
        validate_json(config, validator, dump_path)
    except ValidateError as e:
        raise InvalidConfigError(str(e)) from e


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
