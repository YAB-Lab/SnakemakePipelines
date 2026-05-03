#!/usr/bin/env python3
"""Convert thinned ancestry-prob TSVs to r/qtl CSV with hard genotype calls.

BC cross (par2 only):
  female: par2_prob > 0.5 -> BB else BA   (BA == AB)
  male:   autosomal: same as female
          X-chrom:    par2_prob > 0.5 -> BB else AA  (hemizygous)

F2 cross (par2 + par1):
  per site, genotype probs = [par1, 1 - par1 - par2, par2]  (AA, AB, BB)
  female: argmax -> AA / AB / BB
  male:   autosomal: same
          X-chrom:    AA if par1 > par2, BB if par1 < par2 (no AB)
  Equal-prob homozygotes on X (males) -> "-" (missing)

Both modes: an autosomal marker exactly equal to the autosomal prior
(or [auto_prior, 2*auto_prior, auto_prior] for f2) is reported as "-"
since it's a NA-was-imputed sentinel.

Replaces tsv2csv_bc / tsv2csv_f2 in msg/pull_thin_tsv.py.
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
        next(fh)
        return sum(1 for _ in fh)


# ---------------------------------------------------------------------------
# BC cross
# ---------------------------------------------------------------------------

def write_bc(par2_path, out_path, sex, autosome_prior, x_prior, xchroms):
    auto_sentinel = str(autosome_prior)
    x_sentinel    = str(x_prior)

    with open(par2_path) as fin, open(out_path, "w", newline="") as fout:
        r = csv.reader(fin, delimiter="\t")
        w = csv.writer(fout, delimiter=",")
        header = next(r)
        chroms = [c.split(":", 1)[0] for c in header]
        header[0] = "id"
        w.writerow(header)
        w.writerow(chroms)
        w.writerow("")

        for ind_idx, row in enumerate(r):
            ind = row[0]
            markers = row[1:]
            ind_sex = sex[ind_idx]
            out = [ind]
            for col_idx, marker in enumerate(markers):
                chrom = chroms[col_idx + 1]
                is_x  = chrom in xchroms
                sentinel = x_sentinel if is_x else auto_sentinel

                if marker == sentinel:
                    out.append("-")
                elif float(marker) > 0.5:
                    out.append("BB")
                else:
                    if ind_sex == "1" and is_x:
                        out.append("AA")  # male hemizygous
                    else:
                        out.append("BA")
            w.writerow(out)


# ---------------------------------------------------------------------------
# F2 cross
# ---------------------------------------------------------------------------

def write_f2(par2_path, par1_path, out_path, sex, autosome_prior, xchroms):
    sentinel_genos = [autosome_prior, autosome_prior * 2, autosome_prior]

    with open(par2_path) as f2, open(par1_path) as f1, \
         open(out_path, "w", newline="") as fout:
        r2 = csv.reader(f2, delimiter="\t")
        r1 = csv.reader(f1, delimiter="\t")
        w  = csv.writer(fout, delimiter=",")

        header = next(r2)
        next(r1)
        chroms_row = [c.split(":", 1)[0] for c in header]
        chroms = chroms_row[1:]   # drop the leading empty/id-column entry
        header[0] = "id"
        w.writerow(header)
        w.writerow(chroms_row)
        w.writerow("")

        for ind_idx, row_p2 in enumerate(r2):
            row_p1 = next(r1)
            ind = row_p2[0]
            m_p2 = row_p2[1:]
            m_p1 = row_p1[1:]
            ind_sex = sex[ind_idx]
            out = [ind]

            for z, p2_str in enumerate(m_p2):
                p1_str = m_p1[z]
                p1, p2 = float(p1_str), float(p2_str)
                genos = [p1, 1 - p1 - p2, p2]
                geno_idx = genos.index(max(genos))
                is_x = chroms[z] == "X" or chroms[z] in xchroms

                if ind_sex == "0" or not is_x:
                    if genos == sentinel_genos:
                        out.append("-")
                    else:
                        out.append(["AA", "AB", "BB"][geno_idx])
                else:
                    # male, X-chrom
                    if genos[0] == genos[2]:
                        out.append("-")
                    elif genos[0] > genos[2]:
                        out.append("AA")
                    else:
                        out.append("BB")
            w.writerow(out)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cross", required=True, choices=["bc", "f2"])
    p.add_argument("--par2",  required=True, type=Path)
    p.add_argument("--par1",  type=Path,
                   help="Required when --cross f2")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--xchroms", default="X",
                   help="Comma-separated X-chromosome names")
    p.add_argument("--autosome-prior", type=float, default=0.5)
    p.add_argument("--x-prior",        type=float, default=0.5)
    p.add_argument("--sex-all",   default="")
    p.add_argument("--phenofile", type=Path)
    args = p.parse_args()

    if args.cross == "f2" and not args.par1:
        sys.exit("error: --par1 is required when --cross f2")

    n_inds = count_data_rows(args.par2)
    sex = resolve_sex(n_inds, args.sex_all, args.phenofile)
    xchroms = set(c.strip() for c in args.xchroms.split(",") if c.strip())

    args.output.parent.mkdir(parents=True, exist_ok=True)

    if args.cross == "bc":
        write_bc(args.par2, args.output, sex,
                 args.autosome_prior, args.x_prior, xchroms)
    else:
        write_f2(args.par2, args.par1, args.output, sex,
                 args.autosome_prior, xchroms)


if __name__ == "__main__":
    main()
