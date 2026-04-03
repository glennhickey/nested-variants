#!/usr/bin/env python3
"""Convert GFA path names and header to match chr12.vg PanSN format.

By default, renames paths from plain "SAMPLE.HAP" to PanSN
"SAMPLE#HAP#SAMPLE.HAP" and sets the RS (reference sample) header tag
so the reference is stored as a REFERENCE path.

With --header-only, only sets the RS header tag without renaming paths
(for GFAs that already have PanSN names, e.g. from vg convert -f).

Usage:
    python3 fix-path-names.py [--ref CHM13] < input.gfa > output.gfa
    python3 fix-path-names.py --header-only [--ref CHM13] < input.gfa > output.gfa
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ref", default="CHM13",
                        help="Reference sample name (default: CHM13)")
    parser.add_argument("--header-only", action="store_true",
                        help="Only set RS header tag, don't rename paths")
    args = parser.parse_args()

    for line in sys.stdin:
        if line.startswith("H\t"):
            sys.stdout.write(f"H\tVN:Z:1.1\tRS:Z:{args.ref}\n")
        elif not args.header_only and line.startswith("P\t"):
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
