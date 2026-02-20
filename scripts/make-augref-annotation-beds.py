#!/usr/bin/env python3
"""
make-augref-annotation-beds.py — Per-segment TSV + annotation BEDs → augref-space BEDs.

For each annotation, produces one BED file in augref coordinate space that covers
both on-reference and off-reference contigs.  This lets vcf-stats.R use a simple
interval overlap (foverlaps) to annotate every variant, including those on reference
contigs which were previously dropped.

On-ref entries:  Filter annotation BED for reference assembly rows ({ref}#0#*),
                 prepend "augref_" to contig names (same coordinate system).
Off-ref entries: From per-segment TSV, emit full-segment intervals for segments
                 where annotation overlap > threshold.

Usage:
    python scripts/make-augref-annotation-beds.py \
        --per-segment-tsv {annot-per-segment.tsv} \
        --annotation-beds bed1,bed2,... \
        --annotation-names genes,repeats,... \
        --ref S288C --min-overlap 0.5 \
        --output-dir {OUT_DIR} --output-prefix {OUT_NAME}

Outputs: {OUT_DIR}/{OUT_NAME}.augref-annot-{name}.bed
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--per-segment-tsv', required=True,
                   help='Per-segment annotation TSV from intersect-annotations.py')
    p.add_argument('--annotation-beds', required=True,
                   help='Comma-separated list of original annotation BED files')
    p.add_argument('--annotation-names', required=True,
                   help='Comma-separated list of annotation display names (same order as --annotation-beds)')
    p.add_argument('--ref', required=True,
                   help='Reference sample name (e.g. S288C, CHM13)')
    p.add_argument('--min-overlap', type=float, default=0.5,
                   help='Minimum source_overlap_frac for off-ref segments (default: 0.5)')
    p.add_argument('--output-dir', required=True,
                   help='Output directory for augref BED files')
    p.add_argument('--output-prefix', required=True,
                   help='Output file prefix (e.g. chrI.nested)')
    return p.parse_args()


def read_offref_segments(tsv_path, annot_names, min_overlap):
    """Read per-segment TSV and identify off-ref segments passing overlap threshold.

    For grouped annotations (those with _total rows), use the _total row.
    For ungrouped annotations, use rows where annotation_class == annotation.

    Returns:
        dict: annot_name -> set of (augref_path, source_len) tuples
    """
    offref = defaultdict(set)

    with open(tsv_path) as f:
        header = f.readline().strip().split('\t')
        col = {name: i for i, name in enumerate(header)}

        for line in f:
            fields = line.strip().split('\t')
            augref_path = fields[col['augref_path']]
            source_len = int(fields[col['source_len']])
            annotation = fields[col['annotation']]
            annotation_class = fields[col['annotation_class']]
            source_overlap_frac = float(fields[col['source_overlap_frac']])

            if annotation not in annot_names:
                continue

            # For grouped annotations, use _total row; for ungrouped, use self-named row
            use_row = (annotation_class == '_total') or (annotation_class == annotation)
            if not use_row:
                continue

            if source_overlap_frac >= min_overlap:
                offref[annotation].add((augref_path, source_len))

    return offref


def read_onref_entries(bed_path, ref_name):
    """Read annotation BED and extract entries for the reference assembly.

    Pre-filters with grep for efficiency on large (100+ GB) BEDs.
    Reference entries have contig names matching {ref}#0#*.
    We prepend "augref_" to produce augref-space coordinates.

    Returns:
        list of (chrom, start, end) tuples in augref space
    """
    prefix = f"{ref_name}#0#"
    entries = []

    # Pre-filter with grep to avoid scanning the entire file in Python.
    # grep returns exit code 1 when no matches — that's OK.
    env = dict(os.environ, LC_ALL="C")
    proc = subprocess.Popen(
        ["grep", f"^{prefix}", bed_path],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, env=env)

    for line in proc.stdout:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 3:
            continue
        start = int(fields[1])
        end = int(fields[2])
        entries.append((f"augref_{fields[0]}", start, end))

    proc.wait()  # exit code 1 (no matches) is fine
    return entries


def write_augref_bed(output_path, onref_entries, offref_segments):
    """Write sorted augref-space BED file combining on-ref and off-ref entries.

    On-ref entries: exact intervals from the annotation BED.
    Off-ref entries: full-segment intervals (0 to source_len).
    Uses external sort for bedtools -sorted compatibility and to avoid
    in-memory sort on large entry lists.
    """
    tmp_path = output_path + ".unsorted"
    n = 0

    with open(tmp_path, 'w') as f:
        for chrom, start, end in onref_entries:
            f.write(f"{chrom}\t{start}\t{end}\n")
            n += 1
        for augref_path, source_len in offref_segments:
            f.write(f"{augref_path}\t0\t{source_len}\n")
            n += 1

    # External sort for consistent ordering with bedtools -sorted
    subprocess.run(
        f"LC_ALL=C sort -k1,1 -k2,2n '{tmp_path}' > '{output_path}'",
        shell=True, check=True)
    os.remove(tmp_path)

    return n


def main():
    args = parse_args()

    bed_files = args.annotation_beds.split(',')
    annot_names = args.annotation_names.split(',')

    if len(bed_files) != len(annot_names):
        sys.exit(f"Error: {len(bed_files)} BED files but {len(annot_names)} names")

    os.makedirs(args.output_dir, exist_ok=True)

    # Read off-ref segments from per-segment TSV
    sys.stderr.write(f"Reading per-segment TSV: {args.per_segment_tsv}\n")
    offref = read_offref_segments(args.per_segment_tsv, set(annot_names), args.min_overlap)

    for bed_file, name in zip(bed_files, annot_names):
        sys.stderr.write(f"Processing {name} from {bed_file}...\n")

        # On-ref entries from annotation BED
        onref = read_onref_entries(bed_file, args.ref)
        sys.stderr.write(f"  On-ref entries: {len(onref)}\n")

        # Off-ref segments from per-segment TSV
        offref_segs = offref.get(name, set())
        sys.stderr.write(f"  Off-ref segments: {len(offref_segs)}\n")

        # Write combined augref BED
        output_path = os.path.join(args.output_dir,
                                   f"{args.output_prefix}.augref-annot-{name}.bed")
        n = write_augref_bed(output_path, onref, offref_segs)
        sys.stderr.write(f"  Wrote {n} entries to {output_path}\n")

    sys.stderr.write("Done.\n")


if __name__ == '__main__':
    main()
