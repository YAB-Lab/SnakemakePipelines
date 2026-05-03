#!/usr/bin/env python3
"""Natural-sort a TSV by its first column (individual ID).

Reads a tab-separated file with a header row, sorts the data rows by the
first-column value (natural sort, case-sensitive), and writes the result.

Replaces the in-script sort_file() in msg/pull_thin_tsv.py.
"""
import argparse
import csv
import sys
from pathlib import Path

from _natsort import natsort_key


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",  required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--delim", default="\t",
                   help="Field delimiter (default: tab)")
    args = p.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with open(args.input, newline="") as fin:
        reader = csv.reader(fin, delimiter=args.delim)
        try:
            header = next(reader)
        except StopIteration:
            sys.exit(f"error: empty input {args.input}")
        rows = list(reader)

    rows.sort(key=lambda row: natsort_key(row[0]))

    with open(args.output, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter=args.delim)
        writer.writerow(header)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
