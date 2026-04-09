"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/run_config.yaml
"""
from itertools import product
from textwrap import dedent
from typing import Any, Literal, TypedDict


DATASET_LEVELS = [
    {"name": "subtype", "values": config["valid_subtypes"]},
    {"name": "build", "values": config["valid_builds"]},
    {"name": "resolution", "values": config["valid_resolutions"]},
]


def main():
    validate_config()
    write_subsample_config()
    write_config("results/run_config.yaml")


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

    validate_subsample_config(config["subsample"], DATASET_LEVELS)


def write_subsample_config() -> None:
    """
    Write per-dataset augur subsample configs.
    """
    subsample_entries = config["subsample"]

    for dataset in get_datasets(DATASET_LEVELS):
        subsample_config = build_subsample_config(subsample_entries, dataset)

        a_or_b, build_name, resolution = dataset

        path = build_dir + f"/{a_or_b}/{build_name}/{resolution}/subsample_config.yaml"
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            yaml.dump(subsample_config, f, sort_keys=False, Dumper=NoAliasDumper)
        print(f"Saved subsampling config to {path!r}.", file=sys.stderr)


# FIXME: everything below is independent of DATASET_LEVELS and can be moved to shared/vendored/config.smk or Augur

import jsonschema
from augur.validate import ValidateError, validate_json


ExactDataset = tuple[str, ...]
"""Exact dataset values, ordered to match DATASET_LEVELS."""

SampleOptions = dict[str, Any]
"""Merged or unmerged options for a single named sample."""

class DatasetPatternPart(TypedDict):
    type: Literal["wildcard", "literal", "multivalue"]
    matches: tuple[str, ...] | None


class DatasetLevel(TypedDict):
    name: str
    values: list[str]


SUBSAMPLE_SCHEMA = {
    "type": "array",
    "items": {
        "type": "object",
        "required": ["dataset"],
        "properties": {
            "dataset": {"type": "string"},
            "samples": {
                "type": "object",
                "additionalProperties": {
                    "type": "object",
                    # FIXME: reference augur subsample config schema for exact properties
                },
            },
        },
        "additionalProperties": True,
        "minProperties": 2,
    },
}
SubsampleSchemaValidator = jsonschema.validators.validator_for(SUBSAMPLE_SCHEMA)
SubsampleSchemaValidator.check_schema(SUBSAMPLE_SCHEMA)
SUBSAMPLE_SCHEMA_VALIDATOR = SubsampleSchemaValidator(SUBSAMPLE_SCHEMA)


def validate_subsample_config(
    subsample_entries: list[dict[str, Any]],
    dataset_levels: list[DatasetLevel],
) -> None:
    try:
        validate_json(subsample_entries, SUBSAMPLE_SCHEMA_VALIDATOR, "subsample")
    except ValidateError:
        exit(1)

    # Check things that aren't covered by schema
    for index, subsample_entry in enumerate(subsample_entries, start=1):
        dataset = subsample_entry["dataset"]
        pattern_parts = parse_dataset_pattern(dataset)
        if len(pattern_parts) != len(dataset_levels):
            print(dedent(f"""\
                ERROR: Invalid subsample dataset pattern {dataset!r}.
                Expected {len(dataset_levels)} slash-separated parts matching:
                    {'/'.join(level['name'] for level in dataset_levels)}"""))
            exit(1)

        for pattern_part, level in zip(pattern_parts, dataset_levels):
            if pattern_part["type"] == "wildcard":
                continue

            invalid_values = sorted(set(pattern_part["matches"]) - set(level["values"]))
            if invalid_values:
                print(dedent(f"""\
                    ERROR: Invalid subsample dataset value(s) {invalid_values!r} in pattern {dataset!r}.
                    Expected {level['name']} values from: {level['values']}"""))
                exit(1)


def get_datasets(levels: list[DatasetLevel]) -> list[ExactDataset]:
    """
    Return all datasets requested by config, in the given levels order.
    """
    return product(*(level["values"] for level in levels))


def build_subsample_config(
    subsample_entries: list[dict[str, Any]],
    dataset: ExactDataset,
) -> dict[str, dict[str, SampleOptions]]:
    """
    Build the augur subsample config for a dataset.

    Dataset-matching entries provide option layers that are merged for each
    sample in original config order.
    """
    matching_entries = get_matching_subsample_entries(subsample_entries, dataset)
    if not matching_entries:
        dataset_label = "/".join(dataset)
        print(dedent(f"""\
            ERROR: No subsample entries matched dataset {dataset_label!r}.
            Add a matching subsample entry."""))
        exit(1)

    # Collect sample names
    layers_by_sample = {}
    for entry in matching_entries:
        for sample_name in entry.get("samples", {}):
            layers_by_sample.setdefault(sample_name, [])

    if not layers_by_sample:
        # There is implicitly one sample for this dataset
        layers_by_sample["sample"] = []

    # Add layers
    for entry in matching_entries:
        if "samples" in entry:
            # Sample-specific options
            for sample_name, sample_layer in entry.get("samples", {}).items():
                layers_by_sample[sample_name].append(sample_layer)
        else:
            # Options applied to all samples
            # FIXME: apply these even in the presence of "samples", or enforce mutually exclusiveness in the schema? currently they're silently dropped.
            for sample_name in layers_by_sample:
                layers_by_sample[sample_name].append({
                    key: value
                    for key, value in entry.items()
                    if key != "dataset"
                })

    subsample_config = {"samples": {}}
    for sample_name, layers in layers_by_sample.items():
        subsample_config["samples"][sample_name] = merge_sample_options(layers)

    return subsample_config


def get_matching_subsample_entries(
    entries: list[dict[str, Any]],
    dataset: ExactDataset,
) -> list[dict[str, Any]]:
    """
    Return matching subsample entries for an exact dataset, preserving config order.
    """
    matching_entries = []
    for entry in entries:
        if pattern_matches_dataset(entry["dataset"], dataset):
            matching_entries.append(entry)

    return matching_entries


def pattern_matches_dataset(
    pattern: str,
    dataset: ExactDataset,
) -> bool:
    """
    Return whether a dataset pattern matches an exact dataset.
    """
    pattern_parts = parse_dataset_pattern(pattern)
    if len(pattern_parts) != len(dataset):
        return False

    for pattern_part, dataset_value in zip(pattern_parts, dataset):
        if pattern_part["type"] == "wildcard":
            continue

        if dataset_value not in pattern_part["matches"]:
            return False

    return True


def merge_sample_options(layers: list[SampleOptions]) -> SampleOptions:
    """
    Merge layered sample options.

    `query` values are concatenated with `&`, and `exclude`, `include`,
    `exclude_where`, and `include_where` values are concatenated into one list.
    All other keys follow last-one-wins semantics. A null value removes any
    earlier values for the same key.
    """
    merged = {}
    queries = []
    extendables = {
        "exclude": [],
        "include": [],
        "exclude_where": [],
        "include_where": [],
    }

    for options in layers:
        for key, value in options.items():
            # null = remove
            if value is None:
                if key == "query":
                    queries.clear()
                elif key in extendables:
                    extendables[key].clear()
                else:
                    merged.pop(key, None)
                continue

            if key == "query":
                if value:
                    queries.append(value)
            elif key in extendables:
                if isinstance(value, str):
                    extendables[key].append(value)
                elif isinstance(value, list):
                    extendables[key].extend(value)
            else:
                merged[key] = value

    if queries:
        merged["query"] = " & ".join(f"({q})" for q in queries)
    for key, values in extendables.items():
        if values:
            merged[key] = values

    return merged


def parse_dataset_pattern(pattern: str) -> tuple[DatasetPatternPart, ...]:
    """
    Parse a slash-delimited dataset pattern.
    """
    return tuple(parse_dataset_pattern_part(part) for part in pattern.split("/"))


def parse_dataset_pattern_part(part: str) -> DatasetPatternPart:
    """
    Parse one part of a dataset pattern.

    Supported syntax:
    1. A literal value : 6y
    2. Multiple values : (6y|3y)
    3. All values      : *
    """
    if part == "*":
        return {"type": "wildcard", "matches": None}

    if part.startswith("(") and part.endswith(")"):
        values = tuple(part[1:-1].split("|"))
        if not values or any(not value for value in values):
            print(f"ERROR: Invalid multivalue dataset part {part!r}.")
            exit(1)
        return {"type": "multivalue", "matches": values}

    if any(char in part for char in "()|"):
        print(dedent(f"""\
            ERROR: Invalid subsample dataset part {part!r}.
            Use '*', a literal value, or a whole-part multivalue like '(6y|3y)'."""))
        exit(1)

    return {"type": "literal", "matches": (part,)}


main()
