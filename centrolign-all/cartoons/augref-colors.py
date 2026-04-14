#!/usr/bin/env python3
"""Generate a Bandage CSV to color augref paths in an aug.gfa.

Scheme:
  - Main augref path nodes -> light gray (backbone)
  - Top-N longest _alt paths -> distinct rainbow colors
  - Remaining _alt paths -> uniform muted pink
  - Unassigned nodes (sample-only, shouldn't exist with full augref) -> white

Usage:
    augref-colors.py <aug.gfa> <out.csv> [--top N]

Load in BandageNG via File -> Load CSV data, then switch the Colour
dropdown to "Custom colours".
"""

import argparse
import colorsys
import sys


def parse_gfa(path):
    """Return (seg_lens, paths) — seg_lens[node] = seq length, paths[name] = [node_ids]."""
    seg_lens = {}
    paths = {}
    with open(path) as f:
        for line in f:
            if line.startswith("S\t"):
                fields = line.rstrip("\n").split("\t")
                seg_lens[fields[1]] = len(fields[2])
            elif line.startswith("P\t"):
                fields = line.rstrip("\n").split("\t")
                name = fields[1]
                nodes = [tok[:-1] for tok in fields[2].split(",")]
                paths[name] = nodes
    return seg_lens, paths


def rainbow(n):
    """Return n visually-distinct hex colors via HSV."""
    out = []
    for i in range(n):
        h = (i / n) % 1.0
        # alternate saturation/value to push adjacent colors apart
        s = 0.85 if i % 2 == 0 else 0.65
        v = 0.95 if i % 2 == 0 else 0.75
        r, g, b = colorsys.hsv_to_rgb(h, s, v)
        out.append(f"#{int(r*255):02X}{int(g*255):02X}{int(b*255):02X}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("out_csv")
    ap.add_argument("--top", type=int, default=20,
                    help="Number of longest alts to give distinct colors (default 20)")
    ap.add_argument("--ref-color", default="#D3D3D3", help="Main ref color")
    ap.add_argument("--rest-color", default="#F5B7B1",
                    help="Color for alt paths outside top-N")
    args = ap.parse_args()

    seg_lens, paths = parse_gfa(args.gfa)

    main_ref = [n for n in paths if n.startswith("augref_") and "_alt" not in n]
    if len(main_ref) != 1:
        print(f"Expected exactly one main augref path, got: {main_ref}",
              file=sys.stderr)
        sys.exit(1)
    main_ref = main_ref[0]

    alts = [n for n in paths if "_alt" in n]
    alts.sort(key=lambda a: sum(seg_lens[n] for n in paths[a]), reverse=True)
    top_alts = alts[:args.top]
    rest_alts = alts[args.top:]

    print(f"Main ref path: {main_ref}", file=sys.stderr)
    print(f"Alt paths: {len(alts)} total, coloring top {len(top_alts)}",
          file=sys.stderr)

    palette = rainbow(len(top_alts))

    node_color = {}
    # Lowest priority: main ref backbone
    for n in paths[main_ref]:
        node_color[n] = args.ref_color
    # Middle priority: rest alts
    for alt in rest_alts:
        for n in paths[alt]:
            node_color[n] = args.rest_color
    # Highest priority: top alts (overwrite rest/main if overlap)
    for alt, color in zip(top_alts, palette):
        for n in paths[alt]:
            node_color[n] = color

    with open(args.out_csv, "w") as f:
        f.write("Name,Colour\n")
        for n, c in sorted(node_color.items(), key=lambda kv: int(kv[0])):
            f.write(f"{n},{c}\n")

    print(f"Wrote {len(node_color)} nodes to {args.out_csv}", file=sys.stderr)
    print("\nTop alts (longest first):", file=sys.stderr)
    for alt, color in zip(top_alts, palette):
        length = sum(seg_lens[n] for n in paths[alt])
        print(f"  {color}  {length:>10,} bp  {alt}", file=sys.stderr)


if __name__ == "__main__":
    main()
