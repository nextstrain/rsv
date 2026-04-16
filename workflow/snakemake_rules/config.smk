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
    for dataset in get_datasets(DATASET_LEVELS):
        out = build_subsample_config(config["subsample"], dataset)

        a_or_b, build_name, resolution = dataset

        path = build_dir + f"/{a_or_b}/{build_name}/{resolution}/subsample_config.yaml"
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            yaml.dump(out, f, sort_keys=False, Dumper=NoAliasDumper)
        print(f"Saved subsampling config to {path!r}.", file=sys.stderr)


# FIXME: everything below is independent of DATASET_LEVELS and can be moved to shared/vendored/config.smk or Augur

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


EXTENDABLE_OPTIONS = {
    "exclude",
    "exclude_where",
    "include",
    "include_where",
}


def validate_subsample_config(
    subsample_config: dict[str, Any],
    dataset_levels: list[DatasetLevel],
) -> None:
    if not isinstance(subsample_config, dict):
        print("ERROR: 'subsample' must be a mapping of option names to dataset patterns.")
        exit(1)

    for option, dataset_map in subsample_config.items():
        if not isinstance(dataset_map, dict):
            print(f"ERROR: subsample option {option!r} must be a mapping of dataset patterns to values.")
            exit(1)

        for pattern in dataset_map:
            pattern_parts = parse_dataset_pattern(pattern)
            if len(pattern_parts) != len(dataset_levels):
                print(dedent(f"""\
                    ERROR: Invalid subsample dataset pattern {pattern!r}.
                    Expected {len(dataset_levels)} slash-separated parts matching:
                        {'/'.join(level['name'] for level in dataset_levels)}"""))
                exit(1)

            for pattern_part, level in zip(pattern_parts, dataset_levels):
                if pattern_part["type"] == "wildcard":
                    continue

                invalid_values = sorted(set(pattern_part["matches"]) - set(level["values"]))
                if invalid_values:
                    print(dedent(f"""\
                        ERROR: Invalid subsample dataset value(s) {invalid_values!r} in pattern {pattern!r}.
                        Expected {level['name']} values from: {level['values']}"""))
                    exit(1)


def get_datasets(levels: list[DatasetLevel]) -> list[ExactDataset]:
    """
    Return all datasets requested by config, in the given levels order.
    """
    return product(*(level["values"] for level in levels))


def build_subsample_config(
    subsample_config: dict[str, Any],
    dataset: ExactDataset,
) -> dict[str, dict[str, SampleOptions]]:
    """
    Build the augur subsample config for a dataset.

    For each option, dataset-matching patterns are applied top-to-bottom:
    - A non-dict value applies to all samples.
    - A dict value is a per-sample mapping of sample name to value.

    `query` values are concatenated with `&`, and `exclude`, `include`,
    `exclude_where`, and `include_where` values are concatenated into one list.
    All other keys follow last-one-wins semantics. A null value removes any
    earlier values for the same key.

    Sample names are combined from all dict-value keys that match this dataset.
    If none appear, there is one implicit sample named 'sample'.
    """
    # Collect sample names
    sample_names = set()
    for option, dataset_map in subsample_config.items():
        for pattern, value in dataset_map.items():
            if pattern_matches_dataset(pattern, dataset) and isinstance(value, dict):
                sample_names.update(value.keys())

    if not sample_names:
        # There is implicitly one sample for this dataset
        sample_names = {"sample"}

    # Per-sample accumulators
    scalars = {sample: dict() for sample in sample_names}
    queries = {sample: list() for sample in sample_names}
    extendables = {sample: {option: list() for option in EXTENDABLE_OPTIONS} for sample in sample_names}

    for option, dataset_map in subsample_config.items():
        for pattern, pattern_value in dataset_map.items():
            if not pattern_matches_dataset(pattern, dataset):
                continue

            if isinstance(pattern_value, dict):
                # Sample-specific values
                sample_values = pattern_value.items()
            else:
                # Value applies to all samples
                sample_values = [(sample, pattern_value) for sample in sample_names]

            # Merge values in order of appearance in config
            for sample, value in sample_values:

                # null = remove
                if value is None:
                    if option == "query":
                        queries[sample].clear()
                    elif option in EXTENDABLE_OPTIONS:
                        extendables[sample][option].clear()
                    else:
                        scalars[sample].pop(option, None)

                elif option == "query":
                    queries[sample].append(value)

                elif option in EXTENDABLE_OPTIONS:
                    if isinstance(value, list):
                        extendables[sample][option].extend(value)
                    else:
                        extendables[sample][option].append(value)

                else:
                    # Last one wins
                    scalars[sample][option] = value

    result = {}
    for sample in sample_names:
        options = dict(scalars[sample])
        if queries[sample]:
            options["query"] = " & ".join(f"({q})" for q in queries[sample])
        for option, values in extendables[sample].items():
            if values:
                options[option] = values
        result[sample] = options

    return {"samples": result}


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
