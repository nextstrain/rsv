"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/run_config.yaml
"""
from itertools import product
from textwrap import dedent
from typing import Any, Literal, TypedDict


VALID_DATASET_LEVELS = [
    {"name": "subtype", "values": config["valid_subtypes"]},
    {"name": "build", "values": config["valid_builds"]},
    {"name": "resolution", "values": config["valid_resolutions"]},
]

DATASET_LEVELS_TO_RUN = [
    {"name": "subtype", "values": config["subtypes"]},
    {"name": "build", "values": config["builds_to_run"]},
    {"name": "resolution", "values": config["resolutions_to_run"]},
]

AUGUR_SUBSAMPLE_OPTIONS = {
    # Filter sample options
    "context_sample",
    "drop_sample",
    "exclude",
    "exclude_all",
    "exclude_ambiguous_dates_by",
    "exclude_where",
    "include",
    "include_where",
    "min_date",
    "max_date",
    "min_length",
    "max_length",
    "exclude_invalid",
    "non_nucleotide",
    "query",
    "query_columns",
    "group_by",
    "group_by_weights",
    "probabilistic_sampling",
    "sequences_per_group",
    "max_sequences",
    # Proximal sample options
    "focal_sample",
    "method",
    "k",
    "max_distance",
    "ignore_missing_data",
}


def main():
    validate_config()
    write_subsample_config()
    write_refine_config()
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

    if "custom_subsample" in config and "subsample_datasets" in config:
        print(dedent("""\
            ERROR: 'custom_subsample' and 'subsample_datasets' cannot both be set.

            Use 'custom_subsample' to replace the default dataset-patterned
            subsampling config, or use 'subsample_datasets' to provide raw
            Augur subsample configs for exact datasets."""))
        exit(1)

    validate_subsample_config(config["subsample"], VALID_DATASET_LEVELS, "subsample")

    if "custom_subsample" in config:
        validate_subsample_config(config["custom_subsample"], VALID_DATASET_LEVELS, "custom_subsample")

    if "subsample_datasets" in config:
        validate_subsample_datasets_config(config["subsample_datasets"], VALID_DATASET_LEVELS)


def write_subsample_config() -> None:
    """
    Write per-dataset augur subsample configs.
    """
    for dataset in get_datasets(DATASET_LEVELS_TO_RUN):
        out = get_subsample_config_for_dataset(dataset)

        a_or_b, build_name, resolution = dataset

        path = build_dir + f"/{a_or_b}/{build_name}/{resolution}/subsample_config.yaml"
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            yaml.dump(out, f, sort_keys=False, Dumper=NoAliasDumper)
        print(f"Saved subsampling config to {path!r}.", file=sys.stderr)


def write_refine_config() -> None:
    """
    Write per-dataset augur refine configs.
    """
    for dataset in get_datasets(DATASET_LEVELS_TO_RUN):
        out = build_refine_config(config["refine"], dataset)

        a_or_b, build_name, resolution = dataset

        path = build_dir + f"/{a_or_b}/{build_name}/{resolution}/refine_config.yaml"
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            yaml.dump(out, f, sort_keys=False, Dumper=NoAliasDumper)
        print(f"Saved refine config to {path!r}.", file=sys.stderr)


# FIXME: everything below is independent of VALID_DATASET_LEVELS/DATASET_LEVELS_TO_RUN and can be moved to shared/vendored/config.smk or Augur

ExactDataset = tuple[str, ...]
"""Exact dataset values, ordered to match VALID_DATASET_LEVELS."""

SampleOptions = dict[str, Any]
"""Merged options for a single named sample."""

SubsampleEvent = tuple[Literal["all", "sample"], str | None, str, Any]
"""An ordered option assignment from a matching dataset pattern."""

class DatasetPatternPart(TypedDict):
    type: Literal["wildcard", "literal", "multivalue"]
    matches: tuple[str, ...] | None


class DatasetLevel(TypedDict):
    name: str
    values: list[str]


def validate_subsample_config(
    subsample_config: dict[str, Any],
    dataset_levels: list[DatasetLevel],
    config_key: str,
) -> None:
    if not isinstance(subsample_config, dict):
        print(f"ERROR: {config_key!r} must be a mapping.")
        exit(1)

    for pattern, dataset_config in subsample_config.items():
        if not isinstance(pattern, str):
            print(f"ERROR: {config_key} dataset pattern {pattern!r} must be a string.")
            exit(1)
        _validate_dataset_pattern(pattern, dataset_levels, config_key)

        if not isinstance(dataset_config, dict):
            print(f"ERROR: {config_key} value for pattern {pattern!r} must be a mapping of sample names and/or Augur subsample options.")
            exit(1)

        for key, value in dataset_config.items():
            if not isinstance(key, str):
                print(f"ERROR: {config_key} key {key!r} in pattern {pattern!r} must be a string.")
                exit(1)

            if key in AUGUR_SUBSAMPLE_OPTIONS:
                continue

            if not isinstance(value, dict):
                print(dedent(f"""\
                    ERROR: Invalid key {key!r} in {config_key} pattern {pattern!r}.
                    Keys must be Augur subsample option names or sample names mapped to option mappings."""))
                exit(1)

            invalid_options = sorted(set(value) - AUGUR_SUBSAMPLE_OPTIONS)
            if invalid_options:
                print(dedent(f"""\
                    ERROR: Invalid Augur subsample option(s) {invalid_options!r} in sample {key!r} of {config_key} pattern {pattern!r}.
                    Valid options are: {sorted(AUGUR_SUBSAMPLE_OPTIONS)}"""))
                exit(1)


def _validate_dataset_pattern(
    pattern: str,
    dataset_levels: list[DatasetLevel],
    context: str,
) -> None:
    pattern_parts = parse_dataset_pattern(pattern)
    if len(pattern_parts) != len(dataset_levels):
        print(dedent(f"""\
            ERROR: Invalid {context} dataset pattern {pattern!r}.
            Expected {len(dataset_levels)} slash-separated parts matching:
                {'/'.join(level['name'] for level in dataset_levels)}"""))
        exit(1)

    for pattern_part, level in zip(pattern_parts, dataset_levels):
        if pattern_part["type"] == "wildcard":
            continue

        invalid_values = sorted(set(pattern_part["matches"]) - set(level["values"]))
        if invalid_values:
            print(dedent(f"""\
                ERROR: Invalid {context} dataset value(s) {invalid_values!r} in pattern {pattern!r}.
                Expected {level['name']} values from: {level['values']}"""))
            exit(1)


def validate_subsample_datasets_config(
    subsample_datasets_config: dict[str, Any],
    dataset_levels: list[DatasetLevel],
) -> None:
    if not isinstance(subsample_datasets_config, dict):
        print("ERROR: 'subsample_datasets' must be a mapping of exact dataset names to Augur subsample configs.")
        exit(1)

    for dataset_name, augur_subsample_config in subsample_datasets_config.items():
        if not isinstance(dataset_name, str):
            print(f"ERROR: subsample_datasets key {dataset_name!r} must be a string.")
            exit(1)

        dataset_parts = parse_dataset_pattern(dataset_name)
        if len(dataset_parts) != len(dataset_levels):
            print(dedent(f"""\
                ERROR: Invalid subsample_datasets key {dataset_name!r}.
                Expected an exact dataset name with {len(dataset_levels)} slash-separated parts:
                    {'/'.join(level['name'] for level in dataset_levels)}"""))
            exit(1)

        for dataset_part, level in zip(dataset_parts, dataset_levels):
            if dataset_part["type"] != "literal":
                print(dedent(f"""\
                    ERROR: Invalid subsample_datasets key {dataset_name!r}.
                    Keys must be exact dataset names and cannot use '*' or multivalue patterns."""))
                exit(1)

            invalid_values = sorted(set(dataset_part["matches"]) - set(level["values"]))
            if invalid_values:
                print(dedent(f"""\
                    ERROR: Invalid subsample_datasets value(s) {invalid_values!r} in key {dataset_name!r}.
                    Expected {level['name']} values from: {level['values']}"""))
                exit(1)

        if not isinstance(augur_subsample_config, dict):
            print(f"ERROR: subsample_datasets value for {dataset_name!r} must be an Augur subsample config mapping.")
            exit(1)


def get_datasets(levels: list[DatasetLevel]) -> list[ExactDataset]:
    """
    Return all datasets requested by config, in the given levels order.
    """
    return product(*(level["values"] for level in levels))


def get_subsample_config_for_dataset(dataset: ExactDataset) -> dict[str, Any]:
    """
    Return the Augur subsample config for a dataset.

    Raw exact-dataset configs from 'subsample_datasets' take precedence over
    the generated config. Otherwise, generate from 'custom_subsample' when set,
    or from the default/merged 'subsample' config.
    """
    dataset_name = "/".join(dataset)
    if dataset_name in config.get("subsample_datasets", {}):
        return config["subsample_datasets"][dataset_name]

    subsample_config = config.get("custom_subsample", config["subsample"])
    return build_subsample_config(subsample_config, dataset)


def build_subsample_config(
    subsample_config: dict[str, Any],
    dataset: ExactDataset,
) -> dict[str, dict[str, SampleOptions]]:
    """
    Build the augur subsample config for a dataset.

    Top-level keys are dataset patterns. Matching patterns are applied
    top-to-bottom. Inside a matching pattern:

    - Augur subsample option keys apply to all samples for this dataset.
    - Other keys are sample names mapped to Augur subsample option mappings.

    All options use last-one-wins semantics. A null value removes any earlier
    value for the same key. If no named samples are produced, there is one
    implicit sample named 'sample'.
    """
    sample_names: list[str] = []
    events: list[SubsampleEvent] = []

    for pattern, dataset_config in subsample_config.items():
        if not pattern_matches_dataset(pattern, dataset):
            continue

        for key, value in dataset_config.items():
            if key in AUGUR_SUBSAMPLE_OPTIONS:
                events.append(("all", None, key, value))
                continue

            if key not in sample_names:
                sample_names.append(key)

            for option, option_value in value.items():
                events.append(("sample", key, option, option_value))

    if not sample_names:
        sample_names = ["sample"]

    sample_options = {sample: dict() for sample in sample_names}

    for scope, sample, option, value in events:
        if scope == "all":
            target_samples = sample_names
        else:
            target_samples = [sample]

        for target_sample in target_samples:
            if value is None:
                sample_options[target_sample].pop(option, None)
            else:
                sample_options[target_sample][option] = value

    return {"samples": {sample: dict(sample_options[sample]) for sample in sample_names}}


def build_refine_config(
    refine_config: dict[str, Any],
    dataset: ExactDataset,
) -> dict[str, Any]:
    """
    Build the augur refine config for a dataset.

    Each option's value is either:
    - A scalar/list: applies directly to all datasets.
    - A dict of dataset patterns: matched top-to-bottom.

    Dataset-matching patterns are applied with last-wins semantics. A null
    value removes any earlier value for that option.
    """
    merged = {}
    for option, option_value in refine_config.items():
        if isinstance(option_value, dict):
            matched_values = [
                value
                for pattern, value in option_value.items()
                if pattern_matches_dataset(pattern, dataset)
            ]
        else:
            # Scalar/list value applies to all datasets
            matched_values = [option_value]

        for value in matched_values:
            if value is None:
                merged.pop(option, None)
            else:
                merged[option] = value
    return merged


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
