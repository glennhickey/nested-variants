#!/usr/bin/env python3
"""Convert GFA P-line path names to match chr12.vg PanSN format.

Renames paths from plain "SAMPLE.HAP" to PanSN "SAMPLE#HAP#SAMPLE.HAP"
and sets the RS (reference sample) header tag so the reference is stored
as a REFERENCE path and haplotypes as HAPLOTYPE paths.

Usage:
    python3 fix-path-names.py [--ref CHM13] < input.gfa > output.gfa
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ref", default="CHM13",
                        help="Reference sample name (default: CHM13)")
    args = parser.parse_args()

    for line in sys.stdin:
        if line.startswith("H\t"):
            # Replace header: set GFA 1.1 and RS tag for reference sample
            sys.stdout.write(f"H\tVN:Z:1.1\tRS:Z:{args.ref}\n")
        elif line.startswith("P\t"):
            fields = line.split("\t")
            name = fields[1]
            dot = name.rfind(".")
            if dot != -1:
                sample = name[:dot]
                hap = name[dot + 1:]
                fields[1] = f"{sample}#{hap}#{name}"
            sys.stdout.write("\t".join(fields))
        else:
            sys.stdout.write(line)


if __name__ == "__main__":
    main()
