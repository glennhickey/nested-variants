#!/usr/bin/env python3
"""Annotate indels with tandem repeat motifs.

For each indel allele, checks whether the inserted/deleted sequence is a
tandem repeat (s^n for some minimal motif s) whose motif appears in the
flanking reference, following the criteria of Salehi Nowbandegani et al. 2025
(pantree, Fig 2a/b).

Adds INFO/TR_MOTIF (Number=A): motif string per ALT allele, or "." if not a
tandem repeat.

Dependencies: bcftools, samtools (no Python libraries beyond stdlib).
"""

import argparse
import subprocess
import sys


def minimal_motif(seq):
    """Return the shortest repeating unit s such that seq == s * n."""
    for k in range(1, len(seq) + 1):
        if len(seq) % k == 0 and seq[:k] * (len(seq) // k) == seq:
            return seq[:k]
    return seq


class RefCache:
    """Lazy-loading reference sequence cache using samtools faidx."""

    def __init__(self, ref_path):
        self.ref_path = ref_path
        self.chrom_len = {}
        self._seqs = {}  # chrom -> uppercase sequence string
        self._load_fai()

    def _load_fai(self):
        fai = self.ref_path + ".fai"
        with open(fai) as f:
            for line in f:
                parts = line.split("\t")
                self.chrom_len[parts[0]] = int(parts[1])

    def _load_chrom(self, chrom):
        if chrom in self._seqs:
            return
        result = subprocess.run(
            ["samtools", "faidx", self.ref_path, chrom],
            capture_output=True, text=True
        )
        lines = result.stdout.strip().split("\n")
        self._seqs[chrom] = "".join(lines[1:]).upper() if len(lines) > 1 else ""

    def fetch(self, chrom, start, end):
        """Fetch reference sequence (0-based half-open)."""
        if start >= end:
            return ""
        self._load_chrom(chrom)
        return self._seqs[chrom][start:end]

    def evict(self, chrom):
        """Free memory for a contig no longer needed."""
        self._seqs.pop(chrom, None)


def annotate_allele(ref_str, alt_str, chrom, pos, cache):
    """Return TR_MOTIF value for a single REF/ALT pair. pos is 1-based."""
    ref_str = ref_str.upper()
    alt_str = alt_str.upper()

    ref_len = len(ref_str)
    alt_len = len(alt_str)

    if ref_len == alt_len:
        return "."

    # Common prefix length
    prefix_len = 0
    for i in range(min(ref_len, alt_len)):
        if ref_str[i] == alt_str[i]:
            prefix_len += 1
        else:
            break

    # Extract indel sequence
    if alt_len > ref_len:
        indel_seq = alt_str[prefix_len:]
        is_deletion = False
    else:
        indel_seq = ref_str[prefix_len:]
        is_deletion = True

    if not indel_seq:
        return "."

    motif = minimal_motif(indel_seq)
    motif_len = len(motif)

    # Insert point in 0-based coordinates
    insert_point = (pos - 1) + prefix_len
    clen = cache.chrom_len.get(chrom, 0)

    # Check upstream flank
    up_start = max(0, insert_point - motif_len)
    up_end = insert_point
    if up_end - up_start == motif_len:
        upstream = cache.fetch(chrom, up_start, up_end)
        if upstream == motif:
            return motif

    # Check downstream flank
    if is_deletion:
        down_start = insert_point + len(indel_seq)
    else:
        down_start = insert_point
    down_end = min(down_start + motif_len, clen)
    if down_end - down_start == motif_len:
        downstream = cache.fetch(chrom, down_start, down_end)
        if downstream == motif:
            return motif

    return "."


def is_symbolic(alt):
    return alt.startswith("<") or alt == "*"


def main():
    parser = argparse.ArgumentParser(
        description="Annotate indels with tandem repeat motifs")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    cache = RefCache(args.ref)

    hdr_line = '##INFO=<ID=TR_MOTIF,Number=A,Type=String,Description="Minimal tandem repeat motif for indel alleles (. if not a repeat)">'

    bcf_in = subprocess.Popen(
        ["bcftools", "view", args.vcf],
        stdout=subprocess.PIPE, text=True
    )

    if args.output.endswith(".gz"):
        out_file = open(args.output, "wb")
        out_proc = subprocess.Popen(
            ["bgzip", "-c"], stdin=subprocess.PIPE,
            stdout=out_file, text=True
        )
        out_f = out_proc.stdin
    else:
        out_proc = None
        out_file = None
        out_f = open(args.output, "w")

    n_records = 0
    n_tr = 0
    header_done = False
    prev_chrom = None

    for line in bcf_in.stdout:
        if line.startswith("#"):
            if line.startswith("#CHROM") and not header_done:
                out_f.write(hdr_line + "\n")
                header_done = True
            out_f.write(line)
            continue

        n_records += 1
        fields = line.rstrip("\n").split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alts_str = fields[4]
        info = fields[7]

        # Evict previous contig when we move to a new one (save memory)
        if chrom != prev_chrom:
            if prev_chrom is not None:
                cache.evict(prev_chrom)
            prev_chrom = chrom

        alts = alts_str.split(",")

        if all(is_symbolic(a) for a in alts):
            out_f.write(line)
            continue

        motifs = []
        has_indel = False

        for alt in alts:
            if is_symbolic(alt):
                motifs.append(".")
            elif len(alt) == len(ref):
                motifs.append(".")
            else:
                has_indel = True
                m = annotate_allele(ref, alt, chrom, pos, cache)
                motifs.append(m)

        if has_indel:
            tr_val = ",".join(motifs)
            if info == ".":
                fields[7] = f"TR_MOTIF={tr_val}"
            else:
                fields[7] = info + f";TR_MOTIF={tr_val}"
            if any(m != "." for m in motifs):
                n_tr += 1

        out_f.write("\t".join(fields) + "\n")

    out_f.close()
    if out_proc:
        out_proc.wait()
    if out_file:
        out_file.close()
    bcf_in.wait()

    print(f"Processed {n_records} records, {n_tr} with tandem repeat motif(s)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
