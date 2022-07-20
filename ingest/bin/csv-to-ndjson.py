#!/usr/bin/env python3
"""
Copied from "bin/csv-to-ndjson" in nextstrain/ncov-ingest:
https://github.com/nextstrain/ncov-ingest/blob/2a5f255329ee5bdf0cabc8b8827a700c92becbe4/bin/csv-to-ndjson

Convert CSV on stdin to NDJSON on stdout.
"""
import csv
import json
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="find where sequences are glycosylated",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="csv file")
    parser.add_argument("--output", required=True, help="ndjson file")

    args = parser.parse_args()

# 200 MiB; default is 128 KiB
csv.field_size_limit(200 * 1024 * 1024)

with open(args.input) as file:
    with open(args.output, 'w') as output_file:
        for row in csv.DictReader(file):
            json.dump(row, output_file, allow_nan = False, indent = None, separators = ',:')
            output_file.write("\n")
