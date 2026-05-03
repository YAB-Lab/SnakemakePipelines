#!/usr/bin/env python3
"""Thin per-chromosome markers in an ancestry-probability TSV and replace
NaNs with sex-aware priors.

Modes:
  --diffac F     keep marker i if max(|p[i]-p[i+1]|, |p[i]-p[i+2]|/2) >= F
  --num-markers N  keep the top-N markers per chromosome by that score

For an F2 cross, also align par1 markers to par2's selected positions and
apply the same NA-prior replacement (port of convert_and_thin_Par1).

Sex resolution (used for X-chrom NA replacement):
  --sex-all male|female      apply to all individuals (overrides phenofile)
  --phenofile FILE.tsv       read 'sex' column (0/F=female, 1/M=male)
  default: all female

Replaces convert_and_thin / convert_and_thin_Par1 in msg/pull_thin_tsv.py.
Bug fixes vs legacy:
  - numpy.where -> np.where typo (line 853 of legacy)
  - 'rU' file mode (deprecated) dropped
  - num_inds reference-before-assignment in sex defaulting
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Sex resolution
# ---------------------------------------------------------------------------

_SEX_TRUE = {"0": "0", "f": "0", "F": "0",
             "1": "1", "m": "1", "M": "1"}


def normalize_sex(token):
    return _SEX_TRUE.get(str(token).strip(), "0")


def resolve_sex(num_inds, sex_all, phenofile):
    """Return list of '0'/'1' (female/male), one per individual.

    Order: --sex-all overrides --phenofile; both fall through to all-female.
    """
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


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------

def read_marker_header(path):
    """Return (marker_names, chromosome_per_marker) for the first row."""
    markers = np.genfromtxt(path, max_rows=1, delimiter="\t", dtype=str)
    markers = np.delete(markers, 0, axis=0)  # drop leading empty cell
    chroms = np.array([m.split(":", 1)[0] for m in markers])
    return markers, chroms


def read_indiv_names(path):
    return np.atleast_1d(
        np.genfromtxt(path, comments="#", delimiter="\t",
                      usecols=0, skip_header=1, dtype=str))


def count_data_rows(path):
    with open(path) as fh:
        ncols = len(fh.readline().split("\t"))
        nrows = sum(1 for _ in fh)
    return nrows, ncols


# ---------------------------------------------------------------------------
# Thinning core
# ---------------------------------------------------------------------------

def select_markers(pp, diffac, num_markers, ignore_nan, rng):
    """Return a boolean mask of length pp.shape[1] selecting markers to keep.

    Always keeps the first and last column. Middle columns are scored by:
      score[i] = max(|p[i]-p[i+1]|, |p[i]-p[i+2]|/2)
    """
    n_cols = pp.shape[1]
    if n_cols <= 1:
        return np.ones(n_cols, dtype=bool)

    # NaN-as-change handling matches the legacy convention: NaN->non-NaN
    # contributes diff=1 for the adjacent comparison, diff=2/2=1 for the
    # skip-one comparison.
    if not ignore_nan:
        masked = pp.copy()
        masked[np.isnan(masked)] = 2
        d1 = np.max(np.abs(masked[:, 0:-2] - masked[:, 1:-1]), axis=0)
        d2 = np.max(np.abs(masked[:, 0:-2] - masked[:, 2:]),   axis=0)
        d1[d1 > 1] = 1
        d2[d2 > 1] = 2
    else:
        d1 = np.nanmax(np.abs(pp[:, 0:-2] - pp[:, 1:-1]), axis=0)
        d2 = np.nanmax(np.abs(pp[:, 0:-2] - pp[:, 2:]),   axis=0)
        d1[np.isnan(d1)] = 0
        d2[np.isnan(d2)] = 0

    score = np.maximum(d1, d2 / 2)

    if diffac is not None:
        keep_mid = score >= diffac
    else:
        # Top-N mode: pick threshold so we keep approximately num_markers
        # markers per chromosome (including the two pinned endpoints).
        sorted_desc = np.sort(score)[::-1]
        target_idx = max(0, min(num_markers - 3, len(sorted_desc) - 1))
        thresh = sorted_desc[target_idx]
        if thresh == 1:
            print("Warning: number of markers requested is fewer than the "
                  "number of NA<->non-NA transitions; remaining markers "
                  "will be picked at random.")
        if np.sum(score == thresh) > 1:
            keep_mid = score > thresh
            tied_idx = np.where(score == thresh)[0]
            n_needed = num_markers - int(np.sum(keep_mid)) - 2
            if n_needed > 0 and len(tied_idx) > 0:
                picked = rng.choice(tied_idx,
                                    size=min(n_needed, len(tied_idx)),
                                    replace=False)
                keep_mid[picked] = True
        else:
            keep_mid = score >= thresh

    return np.concatenate(([True], keep_mid, [True]))


def replace_nans_with_priors(pp, is_x_chrom, sex_array, autosome_prior, x_prior):
    """In-place: replace NaN entries with priors. Sex-aware on X chrom."""
    if is_x_chrom:
        male = (sex_array == "1")[:, None]
        female = (sex_array == "0")[:, None]
        male = np.broadcast_to(male, pp.shape)
        female = np.broadcast_to(female, pp.shape)
        nan_mask = np.isnan(pp)
        pp[nan_mask & male]   = x_prior * 2
        pp[nan_mask & female] = x_prior
    else:
        pp[np.isnan(pp)] = autosome_prior


# ---------------------------------------------------------------------------
# Top-level orchestration
# ---------------------------------------------------------------------------

def thin_par2(par2_in, par2_out, *, chroms_spec, xchroms, diffac, num_markers,
              ignore_nan, autosome_prior, x_prior, sex_array, seed):
    markers, marker_chroms = read_marker_header(par2_in)
    inds = read_indiv_names(par2_in)
    rng = np.random.default_rng(seed)

    if chroms_spec.lower() == "all":
        chroms = np.unique(marker_chroms).tolist()
    else:
        chroms = [c.strip() for c in chroms_spec.split(",") if c.strip()]
    xchroms_set = set(c.strip() for c in xchroms.split(",") if c.strip())

    out_pp = None
    out_markers = None

    for chrom in chroms:
        print(f"Thinning chromosome {chrom}", file=sys.stderr)
        col_idx = np.where(marker_chroms == chrom)[0]
        if col_idx.size == 0:
            continue
        marker_subset = markers[col_idx]
        # +1 to skip the leading individual-ID column on data rows
        pp = np.genfromtxt(par2_in, skip_header=1, delimiter="\t",
                           missing_values="NA", filling_values=np.nan,
                           usecols=col_idx + 1)
        if pp.ndim == 1:
            pp = pp.reshape(-1, 1)

        keep = select_markers(pp, diffac, num_markers, ignore_nan, rng)
        pp_thin = pp[:, keep]
        markers_thin = marker_subset[keep]

        replace_nans_with_priors(pp_thin, chrom in xchroms_set, sex_array,
                                 autosome_prior, x_prior)

        if out_pp is None:
            out_pp = pp_thin
            out_markers = markers_thin
        else:
            out_pp = np.concatenate((out_pp, pp_thin), axis=1)
            out_markers = np.concatenate((out_markers, markers_thin), axis=0)

    if out_pp is None:
        sys.exit("error: no chromosomes matched; nothing to write")

    print("Joining marker names to genotype matrix", file=sys.stderr)
    header_row = np.concatenate(([""], out_markers))[None, :]
    body = np.concatenate((inds[:, None], out_pp), axis=1)
    full = np.concatenate((header_row, body), axis=0)

    par2_out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(par2_out, full, delimiter="\t", newline="\n", fmt="%s")


def align_par1_to_par2(par1_in, par2_thinned, par1_out, *, xchroms,
                       autosome_prior, x_prior, sex_array):
    """Pick par1 columns that match par2's thinned marker set, apply priors."""
    markers_p2, _ = read_marker_header(par2_thinned)
    markers_p1, _ = read_marker_header(par1_in)
    p2_set = set(markers_p2.tolist())

    # Preserve par1 column order; +1 because col 0 of the data rows is the ID
    keep_idx = [0] + [i + 1 for i, m in enumerate(markers_p1) if m in p2_set]

    pp = np.genfromtxt(par1_in, delimiter="\t", usecols=keep_idx, dtype=str)

    # Sex-aware X-chrom NA replacement on the data rows only (skip header row)
    chroms_per_kept = np.array([
        markers_p1[i].split(":", 1)[0] for i in [k - 1 for k in keep_idx[1:]]
    ])
    xchroms_set = set(c.strip() for c in xchroms.split(",") if c.strip())
    is_x = np.array([c in xchroms_set for c in chroms_per_kept])

    sex_arr = np.array(sex_array)
    n_data_rows = pp.shape[0] - 1  # subtract header row

    male_rows   = (sex_arr == "1")
    female_rows = (sex_arr == "0")

    # Build a (n_data_rows, n_kept_cols) NaN mask for the data block
    data_block = pp[1:, 1:]   # strip header row + ID col
    is_na = (data_block == "NA")

    # For X cols, sex-aware priors; for autosomal, single prior
    for col_local, is_x_col in enumerate(is_x):
        if is_x_col:
            data_block[is_na[:, col_local] & male_rows,   col_local] = str(x_prior * 2)
            data_block[is_na[:, col_local] & female_rows, col_local] = str(x_prior)
        else:
            data_block[is_na[:, col_local], col_local] = str(autosome_prior)

    pp[1:, 1:] = data_block
    par1_out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(par1_out, pp, delimiter="\t", newline="\n", fmt="%s")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--par2",     required=True, type=Path)
    p.add_argument("--par2-out", required=True, type=Path)
    p.add_argument("--par1",     type=Path,
                   help="Required when --par1-out is given (F2 alignment).")
    p.add_argument("--par1-out", type=Path)
    p.add_argument("--chroms",   default="all",
                   help='"all" or comma-separated contig names')
    p.add_argument("--xchroms",  default="X",
                   help="Comma-separated X-chromosome names")
    p.add_argument("--diffac",       type=float, default=0.01)
    p.add_argument("--num-markers",  type=int, default=None,
                   help="If set, overrides --diffac")
    p.add_argument("--ignore-nan",   type=int, default=0)
    p.add_argument("--autosome-prior", type=float, default=0.5)
    p.add_argument("--x-prior",        type=float, default=0.5)
    p.add_argument("--sex-all",   default="",
                   help="'male' / 'female' to apply to all individuals")
    p.add_argument("--phenofile", type=Path,
                   help="TSV with 'sex' column (used when --sex-all is empty)")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    if args.diffac is not None and args.num_markers is not None:
        # Match legacy: diffac and num_markers are mutually exclusive.
        # We treat num_markers as the override when both are present.
        diffac = None
    else:
        diffac = args.diffac

    n_inds, _ = count_data_rows(args.par2)
    sex_array = np.array(resolve_sex(n_inds, args.sex_all, args.phenofile))

    thin_par2(
        args.par2, args.par2_out,
        chroms_spec=args.chroms,
        xchroms=args.xchroms,
        diffac=diffac,
        num_markers=args.num_markers,
        ignore_nan=bool(args.ignore_nan),
        autosome_prior=args.autosome_prior,
        x_prior=args.x_prior,
        sex_array=sex_array,
        seed=args.seed,
    )

    if args.par1_out:
        if not args.par1:
            sys.exit("error: --par1 is required when --par1-out is given")
        align_par1_to_par2(
            args.par1, args.par2_out, args.par1_out,
            xchroms=args.xchroms,
            autosome_prior=args.autosome_prior,
            x_prior=args.x_prior,
            sex_array=sex_array,
        )


if __name__ == "__main__":
    main()
