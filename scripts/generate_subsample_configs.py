#!/usr/bin/env python3
"""
Generate individual subsample config files from subsample_source.yaml.

This script reads a source YAML file and generates separate config files for
each build/resolution combination.

Usage:
    python scripts/generate_subsample_configs.py <source_file> <output_dir>

Examples:
    python scripts/generate_subsample_configs.py config/subsample_source.yaml config/subsample/
    python scripts/generate_subsample_configs.py profiles/wadoh/subsample_source.yaml profiles/wadoh/subsample/
"""

import sys
import yaml
from copy import deepcopy
from pathlib import Path


def main():
    # Parse command-line arguments
    if len(sys.argv) != 3:
        print("Error: Both source_file and output_dir are required.", file=sys.stderr)
        print("Usage: python scripts/generate_subsample_configs.py <source_file> <output_dir>", file=sys.stderr)
        print("\nExamples:", file=sys.stderr)
        print("  python scripts/generate_subsample_configs.py config/subsample_source.yaml config/subsample/", file=sys.stderr)
        print("  python scripts/generate_subsample_configs.py profiles/wadoh/subsample_source.yaml profiles/wadoh/subsample/", file=sys.stderr)
        sys.exit(1)

    source_path = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])

    # Load source config
    with open(source_path) as f:
        source = yaml.safe_load(f)

    base = source["base"]
    builds = source["builds"]
    resolutions = source["resolutions"]
    builds_to_run = source["builds_to_run"]
    resolutions_to_run = source["resolutions_to_run"]

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate config for each build/resolution combination
    for build in builds_to_run:
        for resolution in resolutions_to_run:
            # Get resolution config
            resolution_config = resolutions[resolution]

            # Build the samples section
            samples = {}
            for sample_name, sample_params in resolution_config["samples"].items():
                # Merge: base -> builds -> sample-specific
                # Later values override earlier ones
                merged = merge_dicts(base, builds[build], sample_params)
                samples[sample_name] = merged

            # Create output config
            output_config = {"samples": samples}

            # Write to file
            output_path = output_dir / f"{build}_{resolution}.yaml"
            with open(output_path, "w") as f:
                yaml.dump(output_config, f, sort_keys=False, default_flow_style=False)

            print(f"Generated {output_path}")

    print(f"\nGenerated {len(builds_to_run) * len(resolutions_to_run)} subsample config files to {output_dir}/")


def merge_dicts(*dicts):
    """Merge multiple dictionaries, with later values overriding earlier ones."""
    result = {}
    for d in dicts:
        # Deep copy to ensure nested structures aren't shared
        result.update(deepcopy(d))
    return result


if __name__ == "__main__":
    main()
