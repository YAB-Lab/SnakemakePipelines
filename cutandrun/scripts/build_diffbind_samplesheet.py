"""Build a DiffBind sample sheet CSV from the YabLab CUT&RUN config.

Expected columns (DiffBind convention):
    SampleID, Tissue, Factor, Condition, Treatment, Replicate,
    bamReads, ControlID, bamControl, Peaks, PeakCaller

Invoked from the Snakemake `diffbind_samplesheet` rule.
"""

import argparse
import csv
import json
import os
from os.path import join

import yaml


def control_type_for(group, sample_controls):
    return sample_controls.get(group, 'IgG')


def control_samples_for(ctrl_type, igg, panh3):
    return igg if ctrl_type == 'IgG' else panh3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config',     required=True, help='Pipeline config.yml')
    ap.add_argument('--out',        required=True, help='Output CSV path')
    ap.add_argument('--peakcaller', default=None,
                    help="Override caller: 'macs3' or 'seacr' (default: from config)")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    out_dir = cfg['OUT_DIR']
    groups          = cfg['GROUPS']
    igg_samples     = list(cfg.get('IGG_SAMPLES', []) or [])
    panh3_samples   = list(cfg.get('PANH3_SAMPLES', []) or [])
    sample_controls = dict(cfg.get('SAMPLE_CONTROLS', {}) or {})
    seacr_mode      = cfg.get('SEACR_MODE', 'stringent')

    caller = (args.peakcaller or cfg.get('DIFFBIND_PEAKCALLER', 'macs3')).lower()
    if caller not in ('macs3', 'seacr'):
        raise SystemExit(f"Unknown peakcaller '{caller}'; expected macs3 or seacr")

    def peak_path(sample):
        if caller == 'macs3':
            return join(out_dir, 'MACS3', sample, f'{sample}_peaks.narrowPeak')
        return join(out_dir, 'SEACR', sample, f'{sample}.{seacr_mode}.bed')

    peak_format = 'narrow' if caller == 'macs3' else 'bed'

    rows = []
    for group, replicates in groups.items():
        ctrl_type    = control_type_for(group, sample_controls)
        ctrl_pool    = control_samples_for(ctrl_type, igg_samples, panh3_samples)
        has_ctrl     = bool(ctrl_pool)
        control_id   = f'{ctrl_type}_pool' if has_ctrl else ''
        control_bam  = join(out_dir, 'pooled_controls', f'{ctrl_type}.pooled.bam') if has_ctrl else ''

        for replicate_idx, sample in enumerate(replicates, start=1):
            rows.append({
                'SampleID':   sample,
                'Tissue':     '',
                'Factor':     group,
                'Condition':  group,
                'Treatment':  '',
                'Replicate':  replicate_idx,
                'bamReads':   join(out_dir, 'Bowtie2', 'target', f'{sample}.filt.bam'),
                'ControlID':  control_id,
                'bamControl': control_bam,
                'Peaks':      peak_path(sample),
                'PeakCaller': peak_format,
            })

    if not rows:
        raise SystemExit("No target samples found in GROUPS; nothing to write.")

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(args.out, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} samples to {args.out} (peakcaller={caller})")


if __name__ == '__main__':
    main()
