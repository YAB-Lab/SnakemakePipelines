#!/usr/bin/env python3
"""
Extract reference-allele information from paired-end BWA-MEM alignments to
two parental genomes.

Replaces the legacy Python2 ``msg/extract-ref-alleles.py``. Reads two
qname-sorted BAMs (one per parent), joins reads by qname, applies PE-aware
filters, writes filtered BAMs and per-contig ref/orth allele tables.

Output layout (preserves legacy format expected by downstream rules):

    {outdir}/aln_{indiv}_par1-filtered.bam
    {outdir}/aln_{indiv}_par2-filtered.bam
    {outdir}/refs/par1/{contig}-ref.alleles
    {outdir}/refs/par1/{contig}-orths.alleles
    {outdir}/refs/par2/{contig}-ref.alleles
    {outdir}/refs/par2/{contig}-orths.alleles
"""
import argparse
import gzip
import os
import re
import shutil
import sys
from collections import defaultdict

import pysam
from Bio import SeqIO
from Bio.Seq import Seq


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('-i', '--individual', required=True,
                   help='Sample / individual ID (used in output filenames)')
    p.add_argument('--bam-par1', required=True,
                   help='qname-sorted BAM aligned to parent 1')
    p.add_argument('--bam-par2', required=True,
                   help='qname-sorted BAM aligned to parent 2')
    p.add_argument('--parent1', required=True,
                   help='Parent 1 reference FASTA (uncompressed or .gz)')
    p.add_argument('--parent2', required=True,
                   help='Parent 2 reference FASTA (uncompressed or .gz)')
    p.add_argument('-o', '--outdir', required=True, help='Output directory')
    p.add_argument('--chroms', default='all',
                   help='Comma-separated chromosomes to keep, or "all" (default: all)')
    p.add_argument('--repeat-threshold', type=int, default=6,
                   help='Min AS-XS difference to consider a mapping non-repetitive '
                        '(default: 6)')
    p.add_argument('--verbose', action='store_true')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Reference loading
# ---------------------------------------------------------------------------

def open_fasta(path):
    if path.endswith(('.gz', '.gzip')):
        return gzip.open(path, 'rt')
    return open(path)


def load_ref(path):
    with open_fasta(path) as fh:
        return SeqIO.to_dict(SeqIO.parse(fh, 'fasta'))


# ---------------------------------------------------------------------------
# qname-merge-join over two qname-sorted BAMs
# ---------------------------------------------------------------------------

def qname_groups(bam):
    """Yield (qname, [reads]) groups from a qname-sorted BAM."""
    current_name = None
    current_reads = []
    for read in bam.fetch(until_eof=True):
        if read.query_name != current_name:
            if current_reads:
                yield current_name, current_reads
            current_name = read.query_name
            current_reads = [read]
        else:
            current_reads.append(read)
    if current_reads:
        yield current_name, current_reads


def merge_join(bam1, bam2):
    """Yield (qname, reads_par1, reads_par2) for qnames present in both BAMs."""
    it1 = qname_groups(bam1)
    it2 = qname_groups(bam2)
    g1 = next(it1, None)
    g2 = next(it2, None)
    while g1 is not None and g2 is not None:
        if g1[0] < g2[0]:
            g1 = next(it1, None)
        elif g1[0] > g2[0]:
            g2 = next(it2, None)
        else:
            yield g1[0], g1[1], g2[1]
            g1 = next(it1, None)
            g2 = next(it2, None)


# ---------------------------------------------------------------------------
# PE-aware filtering
# ---------------------------------------------------------------------------

CIGAR_INDEL_OPS = {1, 2}  # I, D


def primary_mate(reads, is_first):
    """Primary R1 (is_first=True) or R2 (is_first=False) alignment, or None."""
    flag_match = 0x40 if is_first else 0x80
    for r in reads:
        if r.is_secondary or r.is_supplementary:
            continue
        if r.flag & flag_match:
            return r
    return None


def passes_mate_filter(r1, r2, repeat_threshold):
    """Apply PE-aware quality filter to one mate aligned to par1 (r1) and par2 (r2)."""
    if r1 is None or r2 is None:
        return False
    if r1.is_unmapped or r2.is_unmapped:
        return False
    for r in (r1, r2):
        if any(op in CIGAR_INDEL_OPS for op, _ in (r.cigartuples or [])):
            return False
    # Non-repetitive: AS - XS > threshold for both alignments.
    # Missing XS in BWA-MEM means truly unique mapping — keep it.
    for r in (r1, r2):
        try:
            if (r.get_tag('AS') - r.get_tag('XS')) <= repeat_threshold:
                return False
        except KeyError:
            pass
    return True


# ---------------------------------------------------------------------------
# Reference-allele recording (ported from legacy extract-ref-alleles.py)
# ---------------------------------------------------------------------------

_DELETION_RE = re.compile(r'-')


def update_ref_read(read_id, cigar_ops, seq_forward, rev_comp_flag, padding_str):
    """Pad a read against its CIGAR; matches updateRefRead() in the legacy script.

    Returns (refseq_padded, read_padded). '_' marks deletions vs. ref, '-' marks
    insertions, 'X' marks soft/hard clip positions.
    """
    refseq = []
    updated = []
    i = 0
    j = 0
    for op, op_len in cigar_ops:
        if op == 0:  # M
            refseq.append(seq_forward[i:i + op_len])
            updated.append(seq_forward[j:j + op_len])
        elif op == 1:  # I
            refseq.append('-' * op_len)
            updated.append(seq_forward[j:j + op_len])
        elif op == 2:  # D
            refseq.append('_' * op_len)
            updated.append(padding_str * op_len)
        elif op in (4, 5):  # S, H
            refseq.append('X' * op_len)
            updated.append(seq_forward[j:j + op_len])
        elif op == 6:  # P
            pass
        else:
            raise ValueError(f'{read_id}: unsupported CIGAR op {op} (len {op_len})')
        if op != 2:
            i += op_len
            j += op_len
    return ''.join(refseq), ''.join(updated)


def record_reference_alleles(refs, orths, sim_read, read,
                             seq_forward, par1_ref_seq, par2_ref_seq,
                             contig, scaffold):
    """Walk both alignments together; populate refs/orths dicts.

    Mutates refs[par1][contig], refs[par2][contig],
    orths[par1][contig], orths[par2][contig].

    Algorithm preserved verbatim from the legacy script (lines 546-727 of
    msg/extract-ref-alleles.py); only Python2 syntax was modernised.
    """
    pos = sim_read.reference_start + 1
    ref_pos = read.reference_start + 1

    rev_comp = 0
    if sim_read.flag != read.flag and (read.is_reverse or sim_read.is_reverse):
        rev_comp = 1
        ref_pos = read.reference_start + len(read.query_sequence)

    cigar_ops = read.cigartuples or []
    cigar_ops_sim = sim_read.cigartuples or []

    if rev_comp:
        cigar_ops = tuple(reversed(cigar_ops))
        for op, op_len in cigar_ops:
            if op == 1:
                ref_pos -= op_len
            elif op == 2:
                ref_pos += op_len

    sim_ref_seq, _ = update_ref_read(sim_read.query_name, cigar_ops_sim,
                                     str(seq_forward), 0, '+')
    ref_seq, _ = update_ref_read(read.query_name, cigar_ops,
                                 str(seq_forward), rev_comp, '*')

    alleles_par1 = refs['par1'][contig]
    alleles_par2 = refs['par2'][contig]
    orths_par1 = orths['par1'][contig]
    orths_par2 = orths['par2'][contig]

    str_sim = sim_ref_seq
    str_ref = ref_seq
    len_sim = len(str_sim)
    len_ref = len(str_ref)
    idx_sim = 0
    idx_ref = 0

    while idx_sim < len_sim and idx_ref < len_ref:
        # Walk past one-sided deletions
        while (idx_sim < len_sim and idx_ref < len_ref
               and str_sim[idx_sim] == '_' and str_ref[idx_ref] != '_'):
            idx_sim += 1
            pos += 1
        while (idx_sim < len_sim and idx_ref < len_ref
               and str_ref[idx_ref] == '_' and str_sim[idx_sim] != '_'):
            idx_ref += 1
            ref_pos += -1 if rev_comp else 1

        if idx_sim >= len_sim or idx_ref >= len_ref:
            break

        s = str_sim[idx_sim]
        r = str_ref[idx_ref]
        if s == '_' and r == '-':
            idx_ref += 1
            pos += 1
        elif r == '_' and s == '-':
            idx_sim += 1
            if r != '-':
                ref_pos += -1 if rev_comp else 1

        if s not in ('_', '-') and r not in ('_', '-'):
            base1 = par1_ref_seq[pos - 1:pos]
            alleles_par1[pos] = base1
            orths_par1[pos] = f'{contig}\t{pos}\t0\t{base1}'

            base2 = par2_ref_seq[ref_pos - 1:ref_pos]
            if rev_comp:
                base2 = base2.complement()
            alleles_par2[pos] = base2
            orths_par2[pos] = f'{scaffold}\t{ref_pos}\t{rev_comp}\t{base2}'

        if s != '-':
            pos += 1
        if r != '-':
            ref_pos += -1 if rev_comp else 1

        idx_sim += 1
        idx_ref += 1


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_allele_tables(refsdir, refs, orths):
    for sp in ('par1', 'par2'):
        sp_dir = os.path.join(refsdir, sp)
        os.makedirs(sp_dir, exist_ok=True)
        for contig, d in refs[sp].items():
            if not d:
                continue
            with open(os.path.join(sp_dir, f'{contig}-ref.alleles'), 'w') as fh:
                for k, v in sorted(d.items()):
                    fh.write(f'{k}\t{v}\n')
        for contig, d in orths[sp].items():
            if not d:
                continue
            with open(os.path.join(sp_dir, f'{contig}-orths.alleles'), 'w') as fh:
                for k, v in sorted(d.items()):
                    fh.write(f'{k}\t{v}\n')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    refsdir = os.path.join(args.outdir, 'refs')
    if os.path.exists(refsdir):
        shutil.rmtree(refsdir)
    os.makedirs(refsdir)

    print(f'Loading parent 1 FASTA: {args.parent1}', flush=True)
    ref_par1 = load_ref(args.parent1)
    print(f'Loading parent 2 FASTA: {args.parent2}', flush=True)
    ref_par2 = load_ref(args.parent2)

    chroms = None if args.chroms == 'all' else set(args.chroms.split(','))

    bam1 = pysam.AlignmentFile(args.bam_par1, 'rb')
    bam2 = pysam.AlignmentFile(args.bam_par2, 'rb')

    out1_path = os.path.join(args.outdir, f'aln_{args.individual}_par1-filtered.bam')
    out2_path = os.path.join(args.outdir, f'aln_{args.individual}_par2-filtered.bam')
    out1 = pysam.AlignmentFile(out1_path, 'wb', template=bam1)
    out2 = pysam.AlignmentFile(out2_path, 'wb', template=bam2)

    refs = {'par1': defaultdict(dict), 'par2': defaultdict(dict)}
    orths = {'par1': defaultdict(dict), 'par2': defaultdict(dict)}

    n_qnames = 0
    n_kept_mates = 0
    n_seq_mismatch = 0
    for qname, reads_p1, reads_p2 in merge_join(bam1, bam2):
        n_qnames += 1
        for is_first in (True, False):
            r1 = primary_mate(reads_p1, is_first)
            r2 = primary_mate(reads_p2, is_first)
            if not passes_mate_filter(r1, r2, args.repeat_threshold):
                continue

            contig_p1 = bam1.get_reference_name(r1.reference_id)
            contig_p2 = bam2.get_reference_name(r2.reference_id)
            if chroms is not None and contig_p1 not in chroms:
                continue

            seq_p1 = Seq(r1.query_sequence)
            seq_p2 = Seq(r2.query_sequence)
            if r1.is_reverse != r2.is_reverse:
                seq_p2 = seq_p2.reverse_complement()
            if str(seq_p1) != str(seq_p2):
                n_seq_mismatch += 1
                if args.verbose:
                    print(f'Sequence mismatch for {qname}; skipping', file=sys.stderr)
                continue

            try:
                par1_ref = ref_par1[contig_p1].seq
                par2_ref = ref_par2[contig_p2].seq
            except KeyError:
                continue

            out1.write(r1)
            out2.write(r2)
            record_reference_alleles(refs, orths, r1, r2,
                                     seq_p1, par1_ref, par2_ref,
                                     contig_p1, contig_p2)
            n_kept_mates += 1

        if n_qnames % 100000 == 0:
            print(f'Processed {n_qnames:,} qnames; kept {n_kept_mates:,} mates',
                  flush=True)

    out1.close()
    out2.close()
    bam1.close()
    bam2.close()

    print()
    print(f'Total qnames in both BAMs:       {n_qnames:,}')
    print(f'Mate alignments kept:            {n_kept_mates:,}')
    print(f'Mate alignments dropped (seq):   {n_seq_mismatch:,}')

    print('\nWriting per-contig allele tables...', flush=True)
    write_allele_tables(refsdir, refs, orths)
    print('Done.')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'ERROR in extract_ref_alleles: {e}', file=sys.stderr)
        sys.exit(2)
