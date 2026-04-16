#!/usr/bin/env python3
"""Extract per-haplotype variant counts from a multi-sample VCF.

Splits each diploid GT (e.g. 0|1, 0/1) into two haplotype observations
and counts non-ref alleles per (haplotype, variant_type, ref_context).

Output format matches per-sample-types.tsv:
    sample  variant_type  ref_context  count

where sample carries a haplotype suffix (e.g. HG00097.1, HG00097.2).

Usage:
    python3 vcf-per-haplotype-stats.py <input.vcf.gz> > output.tsv
"""

import gzip
import re
import sys
from collections import defaultdict


def classify(ref_len, alt_lens):
    best_diff = 0
    for al in alt_lens:
        d = al - ref_len
        if abs(d) > abs(best_diff):
            best_diff = d
    size = abs(best_diff)
    if size == 0 and ref_len == 1:
        return "SNP"
    elif size == 0:
        return "MNP"
    elif size < 50 and best_diff > 0:
        return "Insertion"
    elif size < 50:
        return "Deletion"
    elif best_diff > 0:
        return "SV Insertion"
    else:
        return "SV Deletion"


def main():
    vcf = sys.argv[1]
    opener = gzip.open if vcf.endswith(".gz") else open

    samples = []
    # counts[(hap_name, variant_type, ref_context)] = int
    counts = defaultdict(int)

    with opener(vcf, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                continue

            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            ref = fields[3]
            alt_str = fields[4]

            ref_len = len(ref)
            alts = [a for a in alt_str.split(",") if a and a not in ("*", ".")]
            if not alts:
                continue
            alt_lens = [len(a) for a in alts]
            vtype = classify(ref_len, alt_lens)
            ref_context = "Off-reference" if "_alt" in chrom else "On-reference"

            sample_fields = fields[9:]

            for i, sf in enumerate(sample_fields):
                if i >= len(samples):
                    break
                gt = sf.split(":")[0]
                alleles = re.split(r"[/|]", gt)
                if len(alleles) < 2:
                    continue

                sample = samples[i]
                for hap_idx, allele in enumerate(alleles[:2]):
                    if allele == "." or allele == "0":
                        continue
                    hap_name = f"{sample}.{hap_idx + 1}"
                    counts[(hap_name, vtype, ref_context)] += 1

    # Ensure every haplotype × variant_type × ref_context has an entry
    type_order = ["SNP", "MNP", "Insertion", "Deletion", "SV Insertion", "SV Deletion"]
    ref_contexts = ["Off-reference", "On-reference"]
    all_haps = sorted(set(h for h, _, _ in counts))
    # Also include haplotypes with zero counts (samples present in VCF header)
    for s in samples:
        for h in [f"{s}.1", f"{s}.2"]:
            if h not in all_haps:
                all_haps.append(h)
    all_haps = sorted(set(all_haps))

    print("sample\tvariant_type\tref_context\tcount")
    for hap in all_haps:
        for vtype in type_order:
            for rc in ref_contexts:
                c = counts.get((hap, vtype, rc), 0)
                print(f"{hap}\t{vtype}\t{rc}\t{c}")


if __name__ == "__main__":
    main()
