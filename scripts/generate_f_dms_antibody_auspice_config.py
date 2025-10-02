"""Generate auspice config for F DMS antibody colorings.

Reads node data JSON to compute min/max escape values for each antibody,
then generates an auspice config with custom viridis color scale.
"""

import argparse
import json


def get_min_max_from_node_data(node_data_file, key):
    """Extract min and max values for a given key from node data JSON."""
    with open(node_data_file) as f:
        node_data = json.load(f)

    values = []
    for node, attrs in node_data.get('nodes', {}).items():
        if key in attrs:
            val = attrs[key]
            # Handle both direct values and dict with 'value' key
            if isinstance(val, dict) and 'value' in val:
                val = val['value']
            if val is not None:
                values.append(val)

    if not values:
        return None, None
    return min(values), max(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--antibodies", nargs="+", required=True, help="List of antibody names")
    parser.add_argument("--continuous-scale", nargs="+", required=True, help="Color hex codes for viridis scale")
    parser.add_argument("--node-data", required=True, help="Node data JSON file with escape scores")
    parser.add_argument("--output", required=True, help="Output auspice config JSON file")
    args = parser.parse_args()

    colorings = []
    for antibody in args.antibodies:
        # Get min/max for total escape
        total_key = f"{antibody}_total_escape"
        total_min, total_max = get_min_max_from_node_data(args.node_data, total_key)

        # Get min/max for max escape
        max_key = f"{antibody}_max_escape"
        max_min, max_max = get_min_max_from_node_data(args.node_data, max_key)

        # Create scale with [value, color] pairs for total escape
        if total_min is not None and total_max is not None:
            num_colors = len(args.continuous_scale)
            total_scale = [
                [total_min + (total_max - total_min) * i / (num_colors - 1), color]
                for i, color in enumerate(args.continuous_scale)
            ]
            colorings.append({
                "key": total_key,
                "title": f"{antibody} Total Escape",
                "type": "continuous",
                "scale": total_scale
            })

        # Create scale for max escape
        if max_min is not None and max_max is not None:
            num_colors = len(args.continuous_scale)
            max_scale = [
                [max_min + (max_max - max_min) * i / (num_colors - 1), color]
                for i, color in enumerate(args.continuous_scale)
            ]
            colorings.append({
                "key": max_key,
                "title": f"{antibody} Max Escape",
                "type": "continuous",
                "scale": max_scale
            })

        # Max escape mutation - categorical (no scale needed)
        colorings.append({
            "key": f"{antibody}_max_escape_mutation",
            "title": f"{antibody} Max Escape Mutation",
            "type": "categorical"
        })

    config_data = {"colorings": colorings}

    with open(args.output, 'w') as f:
        json.dump(config_data, f, indent=2)


if __name__ == "__main__":
    main()
