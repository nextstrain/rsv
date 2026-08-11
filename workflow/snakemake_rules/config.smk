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

try:
    main()
except InvalidConfigError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    exit(1)
