#!/usr/bin/env python3
"""
Generate individual subsample config files from subsample_source.yaml.

This script uses a generic N-dimensional format that works for any number of
dimensions.

Format:
    base: {...}                    # Base parameters applied to all
    params:                        # Dimension-specific parameters
      <dimension>:                 # Dimensions are derived from these keys
        <value>: {...}             # Values are derived from nested keys
    output_name_template: "..."    # Template for output filenames
    samples:                       # Sample definitions (can also be in params)
      <sample>: {...}

Usage:
    python scripts/generate_subsample_configs.py <source_file> <output_dir>

Examples:
    python scripts/generate_subsample_configs.py config/subsample_source.yaml config/subsample/
    python scripts/generate_subsample_configs.py profiles/wadoh/subsample_source.yaml profiles/wadoh/subsample/
"""

import sys
import yaml
from copy import deepcopy
from itertools import product
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

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate configs using N-dimensional format
    if "params" in source:
        count = generate_multidimensional(source, output_dir)
    else:
        print("Error: Unrecognized config format. Expected 'params' key.", file=sys.stderr)
        print("\nThe config should have the following structure:", file=sys.stderr)
        print("  base: {...}", file=sys.stderr)
        print("  params:", file=sys.stderr)
        print("    <dimension>:", file=sys.stderr)
        print("      <value>: {...}", file=sys.stderr)
        print("  output_name_template: '...'", file=sys.stderr)
        print("  samples:", file=sys.stderr)
        print("    <sample>: {...}", file=sys.stderr)
        sys.exit(1)

    print(f"\nGenerated {count} subsample config file(s) to {output_dir}/")


def generate_multidimensional(source, output_dir):
    """
    Generate configs from N-dimensional format.

    Supports arbitrary number of dimensions via Cartesian product.
    Dimensions are automatically derived from keys in params.
    Merge order: base → dim1[v1] → dim2[v2] → ... → dimN[vN] → samples[s]
    """
    base = source.get("base", {})
    params = source["params"]

    # Derive dimensions from params keys
    # Dimension names = top-level keys under params
    # Dimension values = nested keys under each dimension
    dimensions = {dim_name: list(dim_values.keys())
                  for dim_name, dim_values in params.items()}

    output_name_template = source.get("output_name_template", "_".join(f"{{{k}}}" for k in dimensions.keys()))
    top_level_samples = source.get("samples", {"global": {}})

    count = 0

    # Get dimension names and values in order (dict maintains insertion order in Python 3.7+)
    dim_names = list(dimensions.keys())
    dim_values = [dimensions[name] for name in dim_names]

    # Generate all combinations via Cartesian product
    for combination in product(*dim_values):
        # Create context dict mapping dimension names to values
        # e.g., {build: "genome", resolution: "6y"}
        context = dict(zip(dim_names, combination))

        # Start with base parameters
        merged_params = deepcopy(base) if base else {}

        # Apply dimension-specific parameters in order
        for dim_name in dim_names:
            dim_value = context[dim_name]
            dim_params_key = f"{dim_name}"

            if dim_params_key in params and dim_value in params[dim_params_key]:
                dim_specific = params[dim_params_key][dim_value]

                # Check if this dimension defines samples (like RSV's resolution-based samples)
                if dim_specific and "samples" in dim_specific:
                    # This dimension provides sample definitions
                    # We'll handle this separately below
                    pass
                else:
                    # Regular parameters to merge
                    merged_params = merge_dicts(merged_params, dim_specific)

        # Determine which samples to use
        # Priority: dimension-specific samples > top-level samples
        samples_to_use = top_level_samples
        for dim_name in reversed(dim_names):  # Check in reverse order for highest priority
            dim_value = context[dim_name]
            dim_params_key = f"{dim_name}"
            if dim_params_key in params and dim_value in params[dim_params_key]:
                dim_specific = params[dim_params_key][dim_value]
                if dim_specific and "samples" in dim_specific:
                    samples_to_use = dim_specific["samples"]
                    break

        # Build the samples section
        samples = {}
        for sample_name, sample_params in samples_to_use.items():
            # Merge: merged_params (base + all dimensions) → sample-specific
            merged = merge_dicts(merged_params, sample_params)
            samples[sample_name] = merged

        # Create output config
        output_config = {"samples": samples}

        # Generate output name from template
        output_name = output_name_template.format(**context)

        # Write to file
        output_path = output_dir / f"{output_name}.yaml"
        with open(output_path, "w") as f:
            yaml.dump(output_config, f, sort_keys=False, default_flow_style=False)

        print(f"Generated {output_path}")
        count += 1

    return count


def merge_dicts(*dicts):
    """Merge multiple dictionaries, with later values overriding earlier ones."""
    result = {}
    for d in dicts:
        if d is not None:
            # Deep copy to ensure nested structures aren't shared
            result.update(deepcopy(d))
    return result


if __name__ == "__main__":
    main()
