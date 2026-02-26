"""
This part of the workflow deals with configuration.

OUTPUTS:

    results/run_config.yaml
"""
from itertools import product
from textwrap import dedent
from typing import Any, Literal, TypedDict


DATASET_LEVELS = [
    {"name": "subtype", "values": config["subtypes"]},
    {"name": "build", "values": config["builds_to_run"]},
    {"name": "resolution", "values": config["resolutions_to_run"]},
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

    validate_subsample(config["subsample"], DATASET_LEVELS)


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


# NOTE: everything below is independent of DATASET_LEVELS and can be moved to shared/vendored/config.smk or Augur


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


def validate_subsample(
    subsample_entries: list[dict[str, Any]],
    dataset_levels: list[DatasetLevel],
) -> None:
    """
    Validate subsample entries.

    Each entry must be a mapping with a dataset key and one or more
    option keys.
    """
    if not isinstance(subsample_entries, list):
        print(dedent("""\
            ERROR: 'subsample' must be a list of subsample entries.

            Example:
              subsample:
                - dataset: '*/*/*'
                  query: missing_data<1000
                - dataset: '*/*/6y'
                  samples:
                    recent:
                      max_sequences: 3000
            """))
        exit(1)

    for index, subsample_entry in enumerate(subsample_entries, start=1):
        if not isinstance(subsample_entry, dict):
            print(f"ERROR: subsample entry #{index} must be a mapping.")
            exit(1)

        if "dataset" not in subsample_entry:
            print(f"ERROR: subsample entry #{index} is missing required key 'dataset'.")
            exit(1)

        dataset = subsample_entry["dataset"]
        options = {
            key: value
            for key, value in subsample_entry.items()
            if key != "dataset"
        }
        pattern_parts = parse_dataset_pattern(dataset)
        if len(pattern_parts) != len(dataset_levels):
            print(dedent(f"""\
                ERROR: Invalid subsample dataset pattern {dataset!r}.
                Expected {len(dataset_levels)} slash-separated parts matching:
                    {'/'.join(level['name'] for level in dataset_levels)}"""))
            exit(1)

        if not options:
            print(f"ERROR: subsample entry #{index} with dataset {dataset!r} has no options.")
            exit(1)


def get_datasets(levels: list[DatasetLevel]) -> product:
    """
    Return the exact datasets requested by config, in the given levels order.
    """
    return product(*(level["values"] for level in levels))


def build_subsample_config(
    subsample_entries: list[dict[str, Any]],
    dataset: ExactDataset,
) -> dict[str, dict[str, SampleOptions]]:
    """
    Build the augur subsample config for one exact dataset.

    Dataset-matching entries without `samples` provide default option layers
    that apply to every sample. Matching `samples` blocks provide or override
    per-sample options in config order.
    """
    default_options = get_default_options_for_dataset(subsample_entries, dataset)
    samples = get_samples_for_dataset(subsample_entries, dataset)

    subsample_config = {"samples": {}}
    for sample_name, sample_options in samples.items():
        layers = [
            *default_options,
            sample_options,
        ]
        subsample_config["samples"][sample_name] = merge_sample_options(layers)

    return subsample_config


def get_default_options_for_dataset(
    subsample_entries: list[dict[str, Any]],
    dataset: ExactDataset,
) -> list[SampleOptions]:
    """
    Return default option layers for one exact dataset.

    These come from matching subsample entries that do not define `samples`,
    and are applied before each sample's own options.
    """
    return [
        {
            key: value
            for key, value in entry.items()
            if key not in {"dataset", "samples"}
        }
        for entry in get_matching_subsample_entries(subsample_entries, dataset)
        if "samples" not in entry
    ]


def get_samples_for_dataset(
    subsample_entries: list[dict[str, Any]],
    dataset: ExactDataset,
) -> dict[str, SampleOptions]:
    """
    Return merged sample definitions for one exact dataset.

    Matching `samples` blocks are merged sample-by-sample in config order using
    the same semantics as other option layers.
    """
    matching_entries = get_matching_subsample_entries(subsample_entries, dataset)
    sample_layers = {}
    for entry in matching_entries:
        if "samples" not in entry:
            continue

        for sample_name, sample_options in entry["samples"].items():
            sample_layers.setdefault(sample_name, []).append(sample_options)

    if not sample_layers:
        dataset_label = "/".join(dataset)
        print(dedent(f"""\
            ERROR: No samples were defined for subsample dataset {dataset_label!r}.
            Add a matching subsample entry with a 'samples' section."""))
        exit(1)

    return {
        sample_name: merge_sample_options(layers)
        for sample_name, layers in sample_layers.items()
    }


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

    `query` values are concatenated with `&`, and `exclude_where` values are
    concatenated into one list. All other keys follow last-one-wins semantics. A
    null value removes any earlier values for the same key.
    """
    merged = {}
    queries = []
    exclude_wheres = []

    for options in layers:
        # FIXME: verify that 'options' is a dict
        for key, value in options.items():
            # null = remove
            if value is None:
                if key == "query":
                    queries.clear()
                elif key == "exclude_where":
                    exclude_wheres.clear()
                else:
                    merged.pop(key, None)
                continue

            if key == "query":
                if value:
                    queries.append(value)
            elif key == "exclude_where":
                # FIXME: do the same for exclude, include, include_where
                if isinstance(value, str):
                    exclude_wheres.append(value)
                elif isinstance(value, list):
                    exclude_wheres.extend(value)
            else:
                merged[key] = value

    if queries:
        merged["query"] = " & ".join(f"({q})" for q in queries)
    if exclude_wheres:
        merged["exclude_where"] = exclude_wheres

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
