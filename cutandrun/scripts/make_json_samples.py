#!/usr/bin/env python3
"""Build a samples.json from a folder of paired-end FASTQ files.

Expects file names of the form:
    <sample_name>.R1.fastq.gz
    <sample_name>.R2.fastq.gz
"""

import json
import glob
import sys
import argparse
import os


def main():
    parser = argparse.ArgumentParser(description='Make a samples.json for the CUT&RUN pipeline.')
    parser.add_argument('folder', help='Folder containing FASTQ files (*.fastq.gz)')
    parser.add_argument('-o', '--out', default='samples_pe.json',
                        help='Output JSON path (default: samples_pe.json)')
    args = parser.parse_args()

    fastqs = sorted(glob.glob(os.path.join(args.folder, '*.fastq.gz')))
    if not fastqs:
        sys.exit(f"No *.fastq.gz files found in {args.folder}")

    samples = sorted({os.path.basename(f).split('.')[0] for f in fastqs})

    pe = {}
    for s in samples:
        r1 = sorted(f for f in fastqs if f'{s}.' in os.path.basename(f) and 'R1' in os.path.basename(f))
        r2 = sorted(f for f in fastqs if f'{s}.' in os.path.basename(f) and 'R2' in os.path.basename(f))
        if r1 and r2:
            pe[s] = {'R1': r1, 'R2': r2}

    with open(args.out, 'w') as fh:
        json.dump(pe, fh, indent=4, sort_keys=True)
    print(f"Wrote {len(pe)} samples to {args.out}")


if __name__ == '__main__':
    main()
