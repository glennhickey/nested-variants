#!/usr/bin/env python3
"""Compute pairwise Jaccard similarity of GFA paths using node sets.

Usage:
    path-similarity.py <gfa> [<reference_path>]

With a reference path, prints each path's similarity to the reference,
then pairs of non-reference paths ranked by (mutual_sim - avg_ref_sim):
i.e. pairs that are similar to each other but different from the reference.
"""

import re
import sys
from itertools import combinations


def parse_gfa_paths(gfa_path):
    """Return {path_name: frozenset(node_ids)}."""
    paths = {}
    node_re = re.compile(r"(\d+)[+-]")
    with open(gfa_path) as f:
        for line in f:
            if not line.startswith("P\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            name = fields[1]
            nodes = frozenset(node_re.findall(fields[2]))
            paths[name] = nodes
    return paths


def jaccard(a, b):
    if not a and not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    gfa = sys.argv[1]
    ref = sys.argv[2] if len(sys.argv) > 2 else None

    paths = parse_gfa_paths(gfa)
    print(f"Loaded {len(paths)} paths from {gfa}", file=sys.stderr)

    if ref is None:
        ref = next(iter(paths))
        print(f"(no reference given, using first path: {ref})", file=sys.stderr)
    if ref not in paths:
        matches = [n for n in paths if ref in n]
        if len(matches) == 1:
            ref = matches[0]
        else:
            print(f"Reference '{ref}' ambiguous/missing. Candidates: {matches}",
                  file=sys.stderr)
            sys.exit(1)

    ref_nodes = paths[ref]
    others = [n for n in paths if n != ref]

    print(f"\n=== Similarity to {ref} ===")
    print("path\tjaccard_to_ref\tpath_nodes\tref_nodes\tshared")
    ref_sims = {}
    for name in sorted(others, key=lambda n: -jaccard(paths[n], ref_nodes)):
        j = jaccard(paths[name], ref_nodes)
        ref_sims[name] = j
        shared = len(paths[name] & ref_nodes)
        print(f"{name}\t{j:.4f}\t{len(paths[name])}\t{len(ref_nodes)}\t{shared}")

    print(f"\n=== Pairs: high mutual similarity, low ref similarity ===")
    print("rank_score = pair_jaccard - mean(ref_jaccard_of_pair)")
    print("rank_score\tpair_jaccard\tmean_ref_sim\tpath_a\tpath_b")
    pair_rows = []
    for a, b in combinations(others, 2):
        pair_j = jaccard(paths[a], paths[b])
        mean_ref = (ref_sims[a] + ref_sims[b]) / 2
        rank = pair_j - mean_ref
        pair_rows.append((rank, pair_j, mean_ref, a, b))
    pair_rows.sort(reverse=True)
    for rank, pj, mr, a, b in pair_rows[:15]:
        print(f"{rank:.4f}\t{pj:.4f}\t{mr:.4f}\t{a}\t{b}")


if __name__ == "__main__":
    main()
