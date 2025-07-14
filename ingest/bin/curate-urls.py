"""custom curate script to add URLs"""
import sys
import argparse
from typing import Iterable, List, Tuple, Literal, Union

from augur.curate import validate_records
from augur.io.json import dump_ndjson, load_ndjson


def urlify(value: str, url: Union[str, None]) -> str:
    # Encode via "Value: URL"
    # Other options are "Value <URL>" or make a separate column (ndjson field)
    if url:
        return f"{value}: {url}"
    return value

def run(args: argparse.Namespace, records: Iterable[dict]) -> Iterable[dict]:

    for index, record in enumerate(records):
        record = record.copy()

        # modify (ppx) accessionVersion in-place to add a ppx URL (URL refers to the un-versioned accession)
        record['PPX_accession'] = urlify(record['PPX_accession'], f"https://pathoplexus.org/seq/{record['accession']}")

        # NDJSON includes a dataUseTermsUrl which we combine with the dataUseTerms field
        record['dataUseTerms'] = urlify(record['dataUseTerms'], record['dataUseTermsUrl'])

        if record['dataUseTerms'].startswith("RESTRICTED"):
            pass                            # keep the 'restrictedUntil' column as-is
        else:
            record['restrictedUntil'] = ''  # ensure it's empty

        if record['insdcAccession']:
            record['insdcAccession'] = urlify(record['insdcAccession'], f"https://www.ncbi.nlm.nih.gov/nuccore/{record['insdcAccession']}")

        yield record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()

    records = load_ndjson(sys.stdin)

    # Validate records have the same input fields
    validated_input_records = validate_records(records, __doc__, True)

    # Run this custom curate command to get modified records
    modified_records = run(args, validated_input_records)

    # Validate modified records have the same output fields
    validated_output_records = validate_records(modified_records, __doc__, False)

    dump_ndjson(validated_output_records)
