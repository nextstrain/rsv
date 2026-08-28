#!/usr/bin/env python3
"""Automatically generate Markdown configuration documentation from config.schema.yaml and config/configfile.yaml."""

import argparse
import json
import sys
from pathlib import Path
import yaml


def format_default(val):
    if val is None or val == "-":
        return "-"
    if isinstance(val, dict):
        if len(val) > 4:
            return f"({len(val)} entries)"
        return f"`{json.dumps(val)}`"
    if isinstance(val, list):
        if len(val) > 3 and any(isinstance(x, dict) for x in val):
            return f"({len(val)} items)"
        return f"`{json.dumps(val)}`"
    if isinstance(val, bool):
        return "`true`" if val else "`false`"
    if isinstance(val, (int, float, str)):
        return f"`{val}`"
    return f"`{val}`"


def resolve_ref(spec, defs):
    """Resolve $ref pointer against $defs."""
    if isinstance(spec, dict) and "$ref" in spec:
        ref_name = spec["$ref"].split("/")[-1]
        resolved = defs.get(ref_name, spec)
        return resolved, ref_name
    return spec, None


def format_type(spec, defs):
    spec, ref_name = resolve_ref(spec, defs)
    t = spec.get("type")
    if t == "array":
        items = spec.get("items", {})
        items, item_ref = resolve_ref(items, defs)
        if item_ref:
            return f"list[[{item_ref}](#{item_ref.lower()})]"
        return f"list[{items.get('type', 'string')}]"
    if isinstance(t, list):
        return " \\| ".join(t)
    if "enum" in spec:
        return " \\| ".join([f'"{e}"' for e in spec["enum"]])
    if "anyOf" in spec:
        types = []
        for a in spec["anyOf"]:
            at = a.get("type")
            if at == "array":
                it = a.get("items", {}).get("type", "string")
                types.append(f"list[{it}]")
            elif at:
                types.append(at)
        return " \\| ".join(types) if types else "any"
    if ref_name:
        return f"[{ref_name}](#{ref_name.lower()})"
    return str(t or "object")


def get_placeholder(child_dict):
    """Extract dynamic key placeholder directly from the schema description."""
    desc = child_dict.get("description", "").strip()
    if desc.startswith("<") and desc.endswith(">"):
        return desc
    return f"<{desc}>" if desc else "<key>"


def get_container_spec(spec, defs):
    """Identify if a specification represents a nested container object/array."""
    spec, _ = resolve_ref(spec, defs)
    if "properties" in spec:
        return spec, ""
    if "additionalProperties" in spec and isinstance(spec["additionalProperties"], dict):
        child, _ = resolve_ref(spec["additionalProperties"], defs)
        if "properties" in child or "additionalProperties" in child:
            placeholder = get_placeholder(spec["additionalProperties"])
            return child, f".{placeholder}"
    if "patternProperties" in spec:
        for pat, p_spec in spec["patternProperties"].items():
            child, _ = resolve_ref(p_spec, defs)
            if "properties" in child or "additionalProperties" in child:
                placeholder = get_placeholder(p_spec)
                return child, f".{placeholder}"
    if "items" in spec and isinstance(spec["items"], dict):
        child, _ = resolve_ref(spec["items"], defs)
        if "properties" in child:
            return child, "[]"
    return None, None


def generate_docs(schema, defaults):
    defs = schema.get("$defs", {})
    sections = []

    def traverse(path_prefix, spec, depth=2, current_defaults=None, parent_desc=""):
        spec, _ = resolve_ref(spec, defs)
        props = spec.get("properties", {})
        if not props:
            return

        scalars = []
        nested = []

        for name, child in sorted(props.items()):
            full_path = f"{path_prefix}.{name}" if path_prefix else name
            child_default = current_defaults.get(name, "-") if isinstance(current_defaults, dict) else "-"
            container_spec, sub_path = get_container_spec(child, defs)

            if container_spec is not None:
                nested.append((full_path + sub_path, child, container_spec, child_default))
            else:
                p_type = format_type(child, defs)
                p_desc = child.get("description", "-").replace("\n", " ")
                p_def = format_default(child_default)
                scalars.append((full_path, p_type, p_def, p_desc))

        if scalars:
            heading = f"{'#' * depth} `{path_prefix}`" if path_prefix else f"{'#' * depth} Top-Level Settings"
            desc = parent_desc or spec.get("description", "")
            sections.append((heading, desc, scalars))

        for nested_path, orig_child, container_spec, nested_default in nested:
            child_desc = orig_child.get("description", container_spec.get("description", ""))
            traverse(nested_path, container_spec, depth=min(depth + 1, 4), current_defaults=nested_default, parent_desc=child_desc)

    traverse("", schema, depth=2, current_defaults=defaults)

    lines = [
        "<!-- [DO NOT EDIT] This file was generated automatically from cue/ definitions. -->",
        "# Configuration Reference",
        "",
        "This reference is automatically generated from [`cue/`](../cue/).",
        "",
    ]

    for heading, desc, rows in sections:
        lines.append(heading)
        lines.append("")
        if desc:
            lines.append(f"{desc.strip()}\n")
        lines.append("| Parameter Path | Type | Default | Description |")
        lines.append("| :--- | :--- | :--- | :--- |")
        for p, t, d, descr in rows:
            lines.append(f"| `{p}` | `{t}` | {d} | {descr} |")
        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Generate Markdown documentation from JSON/YAML schema and defaults."
    )
    parser.add_argument("schema", type=Path, help="Path to config schema (e.g. config.schema.yaml)")
    parser.add_argument("config", type=Path, help="Path to config file (e.g. config/configfile.yaml)")
    args = parser.parse_args()

    with open(args.schema, "r") as f:
        schema = yaml.safe_load(f)

    with open(args.config, "r") as f:
        defaults = yaml.safe_load(f)

    output = generate_docs(schema, defaults)
    print(output)


if __name__ == "__main__":
    main()
