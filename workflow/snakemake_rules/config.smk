
"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/config_processed.yaml
"""
import os
import sys
import yaml
from copy import deepcopy
from itertools import product
from textwrap import dedent


def main():
    write_config("results/config_raw.yaml")
    validate_config()
    process_subsample_config()
    write_config("results/config_processed.yaml")


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


def process_subsample_config():
    """
    Process the subsample config to expand matrix format into nested dicts.

    If config['subsample'] contains 'defaults' and 'matrix' keys, expands
    the N-dimensional matrix into config['subsample'][dim1][dim2]...[dimN].

    Merge order: defaults → dim1[v1] → dim2[v2] → ... → dimN[vN] → samples[s]
    """
    if not isinstance(config.get("subsample"), dict):
        # Old format (string path) or not present, skip processing
        return

    subsample_config = config["subsample"]

    if "defaults" not in subsample_config or "matrix" not in subsample_config:
        # Already expanded or different format, skip processing
        return

    defaults = subsample_config.get("defaults", {})
    matrix = subsample_config["matrix"]

    # Derive dimensions from matrix keys
    # Dimension names = top-level keys under matrix
    # Dimension values = nested keys under each dimension
    dimensions = {dim_name: list(dim_values.keys())
                  for dim_name, dim_values in matrix.items()}

    # Get dimension names and values in order
    dim_names = list(dimensions.keys())
    dim_values = [dimensions[name] for name in dim_names]

    # Build expanded config
    expanded = {}

    # Generate all combinations via Cartesian product
    for combination in product(*dim_values):
        # Create context dict mapping dimension names to values
        # e.g., {build: "genome", resolution: "6y"}
        context = dict(zip(dim_names, combination))

        # Start with defaults
        merged_params = deepcopy(defaults) if defaults else {}

        # Apply dimension-specific parameters in order
        for dim_name in dim_names:
            dim_value = context[dim_name]

            if dim_name in matrix and dim_value in matrix[dim_name]:
                dim_specific = matrix[dim_name][dim_value]

                # Check if this dimension defines samples
                if dim_specific and "samples" in dim_specific:
                    # This dimension provides sample definitions
                    # We'll handle this separately below
                    pass
                else:
                    # Regular parameters to merge
                    merged_params = merge_dicts(merged_params, dim_specific)

        # Determine which samples to use
        # Priority: last dimension with samples > earlier dimensions with samples
        samples_to_use = None
        for dim_name in reversed(dim_names):  # Check in reverse order for highest priority
            dim_value = context[dim_name]
            if dim_name in matrix and dim_value in matrix[dim_name]:
                dim_specific = matrix[dim_name][dim_value]
                if dim_specific and "samples" in dim_specific:
                    samples_to_use = dim_specific["samples"]
                    break

        # If no samples found, use empty dict with single "global" sample
        if samples_to_use is None:
            samples_to_use = {"global": {}}

        # Build the samples section
        samples = {}
        for sample_name, sample_params in samples_to_use.items():
            # Merge: merged_params (defaults + all dimensions) → sample-specific
            merged = merge_dicts(merged_params, sample_params)
            samples[sample_name] = merged

        # Create output config
        output_config = {"samples": samples}

        # Store in nested dict structure using dimension values as keys
        keys = [context[dim_name] for dim_name in dim_names]
        set_nested_dict(expanded, keys, output_config)

    # Replace subsample config with expanded version
    config["subsample"] = expanded

    print(f"Expanded subsample config matrix into {len(list(product(*dim_values)))} configurations.", file=sys.stderr)


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


def merge_dicts(*dicts):
    """Merge multiple dictionaries, with later values overriding earlier ones."""
    result = {}
    for d in dicts:
        if d is not None:
            # Deep copy to ensure nested structures aren't shared
            result.update(deepcopy(d))
    return result


def set_nested_dict(d, keys, value):
    """Set a value in a nested dict at arbitrary depth."""
    for key in keys[:-1]:
        if key not in d:
            d[key] = {}
        d = d[key]
    d[keys[-1]] = value


main()
