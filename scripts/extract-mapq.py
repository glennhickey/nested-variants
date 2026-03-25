#!/usr/bin/env python3
"""Extract MAPQ distribution from annotated GAM (JSON) or BAM.

For GAM: reads stdin of vg view -aj output (after vg annotate -p -m).
         Classifies reads as on-ref or off-ref based on refpos path names.

For BAM: reads stdin of samtools view output (RNAME + MAPQ columns).
         Classifies reads based on RNAME matching _NNN_alt pattern.

Output: TSV with columns: mapq, ref_context, count
"""

import argparse
import json
import re
import sys
from collections import Counter


def extract_gam(infile):
    """Extract MAPQ from annotated GAM JSON stream."""
    counts = Counter()  # (mapq, ref_context) -> count
    for line in infile:
        try:
            d = json.loads(line)
        except json.JSONDecodeError:
            continue
        mapq = d.get("mapping_quality", 0)
        refpos = d.get("refpos", [])
        names = set(r.get("name", "") for r in refpos)
        has_alt = any(re.search(r"_\d+_alt$", n) for n in names)
        ctx = "Off-reference" if has_alt else "On-reference"
        counts[(mapq, ctx)] += 1
    return counts


def extract_bam(infile):
    """Extract MAPQ from samtools view output (col 1=RNAME, col 2=MAPQ)."""
    counts = Counter()
    for line in infile:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        rname = fields[0]
        try:
            mapq = int(fields[1])
        except ValueError:
            continue
        has_alt = bool(re.search(r"_\d+_alt$", rname))
        ctx = "Off-reference" if has_alt else "On-reference"
        counts[(mapq, ctx)] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description="Extract MAPQ distribution")
    parser.add_argument("--mode", choices=["gam", "bam"], required=True)
    parser.add_argument("--output", "-o", required=True)
    args = parser.parse_args()

    if args.mode == "gam":
        counts = extract_gam(sys.stdin)
    else:
        counts = extract_bam(sys.stdin)

    with open(args.output, "w") as out:
        out.write("mapq\tref_context\tcount\n")
        for (mapq, ctx), n in sorted(counts.items()):
            out.write(f"{mapq}\t{ctx}\t{n}\n")

    total = sum(counts.values())
    on_ref = sum(n for (_, ctx), n in counts.items() if ctx == "On-reference")
    off_ref = total - on_ref
    print(f"Extracted {total} reads: {on_ref} on-ref, {off_ref} off-ref",
          file=sys.stderr)


if __name__ == "__main__":
    main()
