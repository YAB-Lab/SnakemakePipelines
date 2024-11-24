#!/usr/bin/env python

import json
import re
import os
import sys

def parse_filenames(folder_path):
    groups = {}
    pattern = re.compile(r"(.+)_(\d+)\.R\d+\.fastq\.gz")

    for filename in os.listdir(folder_path):
        match = pattern.match(filename)
        if match:
            sample_name, replicate = match.groups()
            full_sample_name = f"{sample_name}_{replicate}"

            if sample_name not in groups:
                groups[sample_name] = []

            if full_sample_name not in groups[sample_name]:
                groups[sample_name].append(full_sample_name)

    for sample in groups:
        groups[sample].sort()

    return groups

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <folder_path>")
        sys.exit(1)

    folder_path = sys.argv[1]
    groups = parse_filenames(folder_path)
    
    # Convert to JSON string with indentation
    groups_json = json.dumps(groups, indent=4)
    
    # Print in the desired format
    print("GROUPS: ", groups_json)

if __name__ == "__main__":
    main()
