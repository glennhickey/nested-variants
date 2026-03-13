#!/usr/bin/env python3
"""Extract per-variant records from a pantree VCF into a standardized TSV.

Pantree VCFs use a non-standard format: off-reference variants have REF=".",
variant types are in INFO/VT, and the non-reference allele is in INFO/NR.
This script normalizes these into a TSV compatible with vcf-stats.R --tsv.

Output columns:
  CHROM POS REF ALT ref_context variant_type size size_signed nonref_af is_repeat
"""

import argparse
import subprocess
import sys


def compute_size(vt, ref, alt, nr):
    """Compute allele size and signed size from pantree fields."""
    is_offref = (ref == ".")

    if vt == "SNP":
        return 0, 0
    if vt == "MNP":
        return 0, 0

    if vt == "INS":
        if is_offref:
            # Off-ref insertion: ALT has the inserted sequence, NR="."
            sz = len(alt) if alt != "." else 0
        else:
            sz = len(alt) - len(ref)
        return abs(sz), abs(sz)

    if vt == "DEL":
        if is_offref:
            # Off-ref deletion: NR has the deleted sequence, ALT="."
            sz = len(nr) if nr != "." else 0
        else:
            sz = len(ref) - len(alt) if alt != "." else len(ref)
        return abs(sz), -abs(sz)

    if vt == "REP":
        # Replacement: compare ALT and NR lengths
        alt_len = len(alt) if alt != "." else 0
        nr_len = len(nr) if nr != "." else 0
        if is_offref:
            sz = alt_len - nr_len
        else:
            sz = len(alt) - len(ref)
        return abs(sz), sz

    # DUP, INV, or unknown
    return 0, 0


def classify_variant(vt, size, size_signed):
    """Map pantree VT + size to our variant_type vocabulary."""
    if vt == "SNP":
        return "SNP"
    if vt == "MNP":
        return "MNP"
    if vt == "REP" and size == 0:
        return "MNP"
    if vt == "INV":
        return "Other"
    if vt == "DUP":
        return "SV Insertion" if size >= 50 else "Insertion"
    if size >= 50:
        return "SV Insertion" if size_signed > 0 else "SV Deletion"
    if size_signed > 0:
        return "Insertion"
    if size_signed < 0:
        return "Deletion"
    return "Other"


def main():
    parser = argparse.ArgumentParser(
        description="Extract per-variant records from pantree VCF")
    parser.add_argument("--vcf", required=True, help="Pantree VCF (plain or bgzipped)")
    parser.add_argument("--output", "-o", required=True, help="Output TSV path")
    args = parser.parse_args()

    query_fmt = r"%CHROM\t%POS\t%REF\t%ALT\t%INFO/VT\t%INFO/NR\t%INFO/AC\t%INFO/AN\t%INFO/TR_MOTIF\n"
    proc = subprocess.Popen(
        ["bcftools", "query", "-f", query_fmt, args.vcf],
        stdout=subprocess.PIPE, text=True
    )

    n = 0
    n_offref = 0
    with open(args.output, "w") as out:
        out.write("CHROM\tPOS\tREF\tALT\tref_context\tvariant_type\tsize\tsize_signed\tnonref_af\tis_repeat\n")

        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            chrom, pos, ref, alt, vt, nr, ac_str, an_str, tr_motif = fields[:9]

            # ref_context
            is_offref = (ref == ".")
            ref_context = "Off-reference" if is_offref else "On-reference"
            if is_offref:
                n_offref += 1

            # Size
            size, size_signed = compute_size(vt, ref, alt, nr)

            # Variant type
            variant_type = classify_variant(vt, size, size_signed)

            # AF
            try:
                ac = int(ac_str)
                an = int(an_str)
                nonref_af = ac / an if an > 0 else 0.0
            except (ValueError, ZeroDivisionError):
                nonref_af = 0.0

            # TR
            is_repeat = "TRUE" if tr_motif not in (".", "") else "FALSE"

            # Use actual allele bases for REF/ALT when available (for Ts/Tv in vcf-stats.R)
            out_ref = ref
            out_alt = alt

            out.write(f"{chrom}\t{pos}\t{out_ref}\t{out_alt}\t{ref_context}\t{variant_type}\t{size}\t{size_signed}\t{nonref_af:.6f}\t{is_repeat}\n")
            n += 1

    proc.wait()

    print(f"Extracted {n} records ({n_offref} off-reference, {n - n_offref} on-reference)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
