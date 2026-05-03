#!/usr/bin/env python3
"""F2-only X-chromosome probability slot conversion.

For each FEMALE individual on X-chromosome columns:
  new par1 prob (homozygous) = old par1 prob + old par2 prob
  new par2 prob (heterozygous) = 1 - new par1 prob

For each MALE individual on X-chromosome columns:
  new par1 prob = 1 - par2 prob   (kept hemizygous)
  new par2 prob unchanged

Autosomal columns are passed through unchanged for both sexes.

Outputs are preserved as side artifacts because downstream user tooling
(outside this Snakemake pipeline) may consume them; the rqtl conversion
itself reads the thinned files directly, not these.

Replaces f2_reduce_prob_slots() in msg/pull_thin_tsv.py.
"""
import argparse
import csv
import sys
from pathlib import Path


def normalize_sex(token):
    s = str(token).strip()
    return "1" if s in ("1", "m", "M") else "0"


def resolve_sex(num_inds, sex_all, phenofile):
    if sex_all:
        s = sex_all.strip().lower()
        if s in ("male", "m", "1"):
            return ["1"] * num_inds
        if s in ("female", "f", "0"):
            return ["0"] * num_inds
        sys.exit(f"error: invalid --sex-all '{sex_all}'")

    if phenofile:
        with open(phenofile) as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            if "sex" not in header:
                sys.exit(f"error: no 'sex' column in {phenofile}")
            idx = header.index("sex")
            sex = [normalize_sex(row[idx]) for row in reader]
        if len(sex) != num_inds:
            sys.exit(f"error: phenofile has {len(sex)} rows, "
                     f"par2 has {num_inds} individuals")
        return sex

    return ["0"] * num_inds


def count_data_rows(path):
    with open(path) as fh:
        next(fh)  # header
        return sum(1 for _ in fh)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--par1",     required=True, type=Path)
    p.add_argument("--par2",     required=True, type=Path)
    p.add_argument("--par1-out", required=True, type=Path)
    p.add_argument("--par2-out", required=True, type=Path)
    p.add_argument("--xchroms",  default="X",
                   help="Comma-separated X-chromosome names")
    p.add_argument("--sex-all",   default="")
    p.add_argument("--phenofile", type=Path)
    args = p.parse_args()

    n_inds = count_data_rows(args.par2)
    sex = resolve_sex(n_inds, args.sex_all, args.phenofile)
    xchroms = set(c.strip() for c in args.xchroms.split(",") if c.strip())

    args.par1_out.parent.mkdir(parents=True, exist_ok=True)
    args.par2_out.parent.mkdir(parents=True, exist_ok=True)

    with open(args.par2) as f2, open(args.par1) as f1, \
         open(args.par2_out, "w", newline="") as o2, \
         open(args.par1_out, "w", newline="") as o1:
        r2 = csv.reader(f2, delimiter="\t")
        r1 = csv.reader(f1, delimiter="\t")
        w2 = csv.writer(o2, delimiter="\t")
        w1 = csv.writer(o1, delimiter="\t")

        header = next(r2)
        next(r1)  # discard par1 header (must match par2 — same marker order)
        w2.writerow(header)
        w1.writerow(header)
        chroms = [c.split(":", 1)[0] for c in header]

        for ind, row2 in enumerate(r2):
            row1 = next(r1)
            out2, out1 = [], []
            for col, p2_val in enumerate(row2):
                p1_val = row1[col]
                if col == 0 or chroms[col] not in xchroms:
                    out2.append(p2_val)
                    out1.append(p1_val)
                    continue
                # X-chrom column
                if sex[ind] == "1":
                    # male: par1 = 1 - par2 (par2 unchanged)
                    out2.append(p2_val)
                    out1.append(1 - float(p2_val))
                else:
                    # female: par1 = par2 + par1, par2 = 1 - par1_new
                    new_par1 = float(p2_val) + float(p1_val)
                    out2.append(1 - new_par1)
                    out1.append(new_par1)
            w2.writerow(out2)
            w1.writerow(out1)


if __name__ == "__main__":
    main()
