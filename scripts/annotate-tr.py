#!/usr/bin/env python3
"""Annotate indels with tandem repeat motifs.

For each indel allele, checks whether the inserted/deleted sequence is a
tandem repeat (s^n for some minimal motif s) whose motif appears in the
flanking reference, following the criteria of Salehi Nowbandegani et al. 2025
(pantree, Fig 2a/b).

Adds INFO/TR_MOTIF (Number=A): motif string per ALT allele, or "." if not a
tandem repeat.
"""

import argparse
import sys

import pysam


def minimal_motif(seq: str) -> str:
    """Return the shortest repeating unit s such that seq == s * n."""
    for k in range(1, len(seq) + 1):
        if len(seq) % k == 0 and seq[:k] * (len(seq) // k) == seq:
            return seq[:k]
    return seq


def annotate_allele(ref: str, alt: str, chrom: str, pos: int,
                    fasta: pysam.FastaFile) -> str:
    """Return TR_MOTIF value for a single REF/ALT pair.

    pos is 1-based VCF POS.
    """
    ref = ref.upper()
    alt = alt.upper()

    ref_len = len(ref)
    alt_len = len(alt)

    # Not an indel
    if ref_len == alt_len:
        return "."

    # Find common prefix length
    prefix_len = 0
    for i in range(min(ref_len, alt_len)):
        if ref[i] == alt[i]:
            prefix_len += 1
        else:
            break

    # Extract indel sequence from the longer allele
    if alt_len > ref_len:
        # Insertion
        indel_seq = alt[prefix_len:]
        is_deletion = False
    else:
        # Deletion
        indel_seq = ref[prefix_len:]
        is_deletion = True

    if not indel_seq:
        return "."

    motif = minimal_motif(indel_seq)
    motif_len = len(motif)

    # Insert point in 0-based coordinates
    insert_point = (pos - 1) + prefix_len
    chrom_len = fasta.get_reference_length(chrom)

    # Check upstream flank
    up_start = max(0, insert_point - motif_len)
    up_end = insert_point
    if up_end - up_start == motif_len:
        upstream = fasta.fetch(chrom, up_start, up_end).upper()
        if upstream == motif:
            return motif

    # Check downstream flank
    if is_deletion:
        down_start = insert_point + len(indel_seq)
    else:
        down_start = insert_point
    down_end = min(down_start + motif_len, chrom_len)
    if down_end - down_start == motif_len:
        downstream = fasta.fetch(chrom, down_start, down_end).upper()
        if downstream == motif:
            return motif

    return "."


def is_symbolic(alt: str) -> bool:
    """True for symbolic ALTs like <DEL>, *, etc."""
    return alt.startswith("<") or alt == "*"


def main():
    parser = argparse.ArgumentParser(
        description="Annotate indels with tandem repeat motifs")
    parser.add_argument("--vcf", required=True, help="Input VCF (.vcf or .vcf.gz)")
    parser.add_argument("--ref", required=True, help="Indexed reference FASTA (.fa.gz or .fa)")
    parser.add_argument("-o", "--output", required=True, help="Output VCF (.vcf or .vcf.gz)")
    args = parser.parse_args()

    fasta = pysam.FastaFile(args.ref)
    vcf_in = pysam.VariantFile(args.vcf)

    # Add TR_MOTIF header
    vcf_in.header.info.add(
        "TR_MOTIF", "A", "String",
        "Minimal tandem repeat motif for indel alleles (. if not a repeat)")

    out_mode = "wz" if args.output.endswith(".gz") else "w"
    vcf_out = pysam.VariantFile(args.output, out_mode, header=vcf_in.header)

    n_records = 0
    n_tr = 0

    for rec in vcf_in:
        n_records += 1

        # Skip records with only symbolic ALTs
        if all(is_symbolic(str(a)) for a in rec.alts or []):
            vcf_out.write(rec)
            continue

        motifs = []
        has_indel = False

        for alt in (rec.alts or []):
            alt_str = str(alt)
            if is_symbolic(alt_str):
                motifs.append(".")
                continue

            ref_str = rec.ref
            if len(alt_str) == len(ref_str):
                # SNP or MNP
                motifs.append(".")
            else:
                has_indel = True
                m = annotate_allele(ref_str, alt_str, rec.chrom, rec.pos, fasta)
                motifs.append(m)

        if has_indel:
            rec.info["TR_MOTIF"] = tuple(motifs)
            if any(m != "." for m in motifs):
                n_tr += 1

        vcf_out.write(rec)

    vcf_out.close()
    vcf_in.close()
    fasta.close()

    print(f"Processed {n_records} records, {n_tr} with tandem repeat motif(s)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
