"""custom curate script to add URLs"""
import sys
import argparse
from typing import Iterable

from augur.curate import validate_records
from augur.io.json import dump_ndjson, load_ndjson

def run(args: argparse.Namespace, records: Iterable[dict]) -> Iterable[dict]:

    for index, record in enumerate(records):
        record = record.copy()

        ppx_accession = record.get('PPX_accession', None) # versioned
        insdc_accession = record.get('INSDC_accession', None) # versioned

        # Add INSDC_accession__url and PPX_accession__url fields to NDJSON records
        record['PPX_accession__url'] = f"https://pathoplexus.org/seq/{ppx_accession}" \
            if ppx_accession \
            else ""
        record['INSDC_accession__url'] = f"https://www.ncbi.nlm.nih.gov/nuccore/{insdc_accession}" \
            if insdc_accession \
            else ""

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
