#!/usr/bin/env python3

import json
import glob
import sys
import argparse
import os

def msg(name=None):
    return '''Usage: make_json_samples.py <folder>

    fastq file names should have the following format:
        paired-end: <sample_name>.R1.fastq.gz
                    <sample_name>.R2.fastq.gz

        single-end: <sample_name>.R1.fastq.gz
    '''

parser = argparse.ArgumentParser(description='Make a samples.json file with sample names and file names.', usage=msg())
parser.add_argument('folder', help='Folder containing FASTQ files')

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()
folder = args.folder

fastqs = glob.glob(os.path.join(folder, '*.fastq.gz'))

peFILES = {}
seFILES = {}

# Extract sample names from filenames
SAMPLES = [os.path.basename(fastq).split('.')[0] for fastq in fastqs]

for sample in SAMPLES:
    mate1 = lambda fastq: sample in fastq and 'R1' in fastq
    mate2 = lambda fastq: sample in fastq and 'R2' in fastq
    if any('R2' in s for s in sorted(filter(mate2, fastqs))):
        peFILES[sample] = {}
        peFILES[sample]['R1'] = sorted(filter(mate1, fastqs))
        peFILES[sample]['R2'] = sorted(filter(mate2, fastqs))
    else:
        seFILES[sample] = {}
        seFILES[sample]['R1'] = sorted(filter(mate1, fastqs))

js_pe = json.dumps(peFILES, indent=4, sort_keys=True)
js_se = json.dumps(seFILES, indent=4, sort_keys=True)

with open('samples_pe.json', 'w') as pe_file:
    pe_file.write(js_pe)

with open('samples_se.json', 'w') as se_file:
    se_file.write(js_se)
