#!/usr/bin/env python3
"""Extract compact variant stats from a VCF, streaming line by line.

Computes variant type and size from allele lengths without storing full
sequences in memory. Outputs a TSV compatible with vcf-stats.R --tsv.

Usage:
    python3 vcf-extract-stats.py input.vcf.gz > output.tsv
"""

import gzip
import sys


def classify(ref_len, alt_lens):
    """Classify variant from REF length and list of ALT lengths."""
    best_diff = 0
    for al in alt_lens:
        d = al - ref_len
        if abs(d) > abs(best_diff):
            best_diff = d
    size = abs(best_diff)
    if size == 0 and ref_len == 1:
        return "SNP", 0, best_diff
    elif size == 0:
        return "MNP", 0, best_diff
    elif size < 50 and best_diff > 0:
        return "Insertion", size, best_diff
    elif size < 50:
        return "Deletion", size, best_diff
    elif best_diff > 0:
        return "SV Insertion", size, best_diff
    else:
        return "SV Deletion", size, best_diff


def parse_info(info, key):
    """Extract a value from the INFO field."""
    for field in info.split(";"):
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return None


def main():
    vcf = sys.argv[1]
    opener = gzip.open if vcf.endswith(".gz") else open

    print("CHROM\tPOS\tref_context\tvariant_type\tsize\tsize_signed\tnonref_af\tis_repeat\ttstv")

    with opener(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            # Parse only what we need — avoid split() on the full line
            # VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO ...
            # We need CHROM(0), POS(1), REF(3), ALT(4), INFO(7)
            # Split only up to column 8 to avoid touching genotype columns
            fields = line.split("\t", 9)
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt_str = fields[4]
            info = fields[7]

            ref_len = len(ref)
            alts = [a for a in alt_str.split(",") if a and a not in ("*", ".")]
            if not alts:
                continue
            alt_lens = [len(a) for a in alts]

            vtype, size, size_signed = classify(ref_len, alt_lens)

            # On-ref vs off-ref
            ref_context = "Off-reference" if "_alt" in chrom else "On-reference"

            # AF from INFO
            af_raw = parse_info(info, "AF")
            nonref_af = ""
            if af_raw:
                try:
                    afs = [float(x) for x in af_raw.split(",") if x != "."]
                    nonref_af = min(sum(afs), 1.0) if afs else ""
                except ValueError:
                    pass

            # TR_MOTIF from INFO
            tr_raw = parse_info(info, "TR_MOTIF")
            is_repeat = "TRUE" if (tr_raw and any(c not in ".,\t" for c in tr_raw)) else "FALSE"

            # Ts/Tv for biallelic SNPs
            tstv = ""
            if vtype == "SNP" and ref_len == 1 and len(alts) == 1 and len(alts[0]) == 1:
                pair = ref + alts[0]
                tstv = "Ts" if pair in ("AG", "GA", "CT", "TC") else "Tv"

            print(f"{chrom}\t{pos}\t{ref_context}\t{vtype}\t{size}\t{size_signed}\t{nonref_af}\t{is_repeat}\t{tstv}")


if __name__ == "__main__":
    main()
