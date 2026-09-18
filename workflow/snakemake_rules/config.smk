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

try:
    main()
except InvalidConfigError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    exit(1)
