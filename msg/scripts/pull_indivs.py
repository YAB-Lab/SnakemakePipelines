#!/usr/bin/env python3
"""Subset rows of a TSV by individual ID (first column).

When --indivs is "all" (case-insensitive) or empty, this is a pass-through:
the full input is copied to the output unchanged. Otherwise, --indivs is a
comma-separated list of patterns; rows whose first-column value contains
any pattern (substring match, matching legacy regex behavior) are kept.

The header row is always passed through.

Replaces the in-script pull_idds() / file_to_pulled() in msg/pull_thin_tsv.py.
"""
import argparse
import re
import shutil
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",  required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--indivs", default="all",
                   help='"all" (default) or comma-separated patterns')
    args = p.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    spec = args.indivs.strip()
    if not spec or spec.lower() == "all":
        shutil.copyfile(args.input, args.output)
        return

    patterns = [re.compile(re.escape(p.strip()))
                for p in spec.split(",") if p.strip()]

    with open(args.input) as fin, open(args.output, "w") as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            name = line.split("\t", 1)[0]
            if any(pat.search(name) for pat in patterns):
                fout.write(line)


if __name__ == "__main__":
    main()
