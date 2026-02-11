#!/usr/bin/env python3
"""
Intersect off-reference coordinates with annotation BED files using bedtools.

Takes a nesting TSV file with off-reference coordinate mappings and intersects
them with annotation BED files (genes, repeats, segmental duplications, etc.)
to calculate overlap statistics for downstream visualization.
"""

import os
import sys
import subprocess
import argparse
import tempfile
from collections import defaultdict


def check_dependencies():
    """Check that required external tools are available."""
    missing = []
    for tool in ['bedtools']:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing.append(tool)
    if missing:
        sys.exit(f"Error: missing required tools: {', '.join(missing)}")


def extract_coordinates(nesting_file, offref_bed, onref_bed):
    """
    Extract both off-reference and on-reference coordinates from nesting TSV file.

    Args:
        nesting_file: Path to nesting.tsv file with columns:
                      offref_contig, offref_start, offref_end, node, onref_contig, onref_start, onref_end
        offref_bed: Output BED file for off-reference coordinates (columns 1-3)
        onref_bed: Output BED file for on-reference coordinates (columns 5-7)

    Returns:
        Tuple of (total_offref_bp, total_onref_bp, num_skipped)
    """
    total_offref_bp = 0
    total_onref_bp = 0
    num_skipped = 0

    with open(nesting_file) as infile, \
         open(offref_bed, 'w') as offref_out, \
         open(onref_bed, 'w') as onref_out:

        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) < 7:
                continue

            # Extract off-reference coordinates (columns 1-3, 0-indexed: 0-2)
            offref_contig = fields[0]
            offref_start = int(fields[1])
            offref_end = int(fields[2])

            # Extract on-reference coordinates (columns 5-7, 0-indexed: 4-6)
            onref_contig = fields[4]
            onref_start = int(fields[5])
            onref_end = int(fields[6])

            # Skip records with invalid coordinates (negative or zero-length)
            if offref_start < 0 or offref_end < 0 or onref_start < 0 or onref_end < 0:
                num_skipped += 1
                continue
            if offref_start >= offref_end or onref_start >= onref_end:
                num_skipped += 1
                continue

            offref_out.write(f'{offref_contig}\t{offref_start}\t{offref_end}\n')
            onref_out.write(f'{onref_contig}\t{onref_start}\t{onref_end}\n')

            total_offref_bp += offref_end - offref_start
            total_onref_bp += onref_end - onref_start

    return total_offref_bp, total_onref_bp, num_skipped


def run_bedtools_intersect(coords_bed, annot_bed, output_file=None):
    """
    Run bedtools intersect and return results.

    Args:
        coords_bed: BED file with off-reference coordinates
        annot_bed: BED file with annotations
        output_file: Optional output file (if None, returns to stdout)

    Returns:
        Subprocess result
    """
    cmd = ['bedtools', 'intersect', '-a', coords_bed, '-b', annot_bed, '-wa', '-wb']

    if output_file:
        with open(output_file, 'w') as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    else:
        result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        sys.stderr.write(f"bedtools intersect failed: {result.stderr}\n")

    return result


def calculate_coverage(coords_bed, annot_bed):
    """
    Calculate how many bp of coords_bed overlap with annot_bed.

    Args:
        coords_bed: BED file with off-reference coordinates
        annot_bed: BED file with annotations

    Returns:
        Total base pairs of overlap
    """
    cmd = ['bedtools', 'intersect', '-a', coords_bed, '-b', annot_bed, '-u']
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        return 0

    total_bp = 0
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if len(fields) >= 3:
            total_bp += int(fields[2]) - int(fields[1])

    return total_bp


def count_intersections(coords_bed, annot_bed):
    """
    Count number of intersection events.

    Args:
        coords_bed: BED file with off-reference coordinates
        annot_bed: BED file with annotations

    Returns:
        Number of intersections
    """
    cmd = ['bedtools', 'intersect', '-a', coords_bed, '-b', annot_bed, '-c']
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        return 0

    total_count = 0
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if len(fields) >= 4:
            total_count += int(fields[3])

    return total_count


def calculate_grouped_stats(intersect_file, group_column, total_bp):
    """
    Calculate statistics grouped by a specific column value.

    Args:
        intersect_file: Path to intersection output file
        group_column: 1-indexed column number to group by (in the annotation, after -wb)
        total_bp: Total base pairs in coordinate regions

    Returns:
        Dictionary mapping group values to their statistics
    """
    grouped_stats = defaultdict(lambda: {'overlap_bp': 0, 'num_intersections': 0})

    if not os.path.exists(intersect_file) or os.path.getsize(intersect_file) == 0:
        return {}

    with open(intersect_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue

            # First 3 fields are from coords_bed (-wa), rest are from annot_bed (-wb)
            # So annotation columns start at index 3
            annot_start_idx = 3
            group_col_idx = annot_start_idx + (group_column - 1)

            if group_col_idx >= len(fields):
                continue

            group_value = fields[group_col_idx]

            # Get the coords region size
            coord_start = int(fields[1])
            coord_end = int(fields[2])
            overlap = coord_end - coord_start

            grouped_stats[group_value]['overlap_bp'] += overlap
            grouped_stats[group_value]['num_intersections'] += 1

    # Calculate percentages
    for group_value in grouped_stats:
        overlap_bp = grouped_stats[group_value]['overlap_bp']
        grouped_stats[group_value]['percent_coverage'] = (overlap_bp / total_bp * 100) if total_bp > 0 else 0

    return dict(grouped_stats)


def process_annotations(coords_bed, annot_files, output_dir, total_bp, coord_type, group_column=None):
    """
    Process all annotation files and generate statistics.

    Args:
        coords_bed: BED file with coordinates
        annot_files: List of annotation BED file paths
        output_dir: Directory for output files
        total_bp: Total base pairs in coordinate regions
        coord_type: String describing coordinate type ('offref' or 'onref')
        group_column: Optional 1-indexed column number to group statistics by

    Returns:
        Dictionary of statistics per annotation file
    """
    stats = {}

    for annot_file in annot_files:
        annot_name = os.path.splitext(os.path.basename(annot_file))[0]
        sys.stderr.write(f"  Processing {annot_name} ({coord_type})...\n")

        # Run intersection and save detailed output
        intersect_out = os.path.join(output_dir, f'{annot_name}.{coord_type}.intersect.bed')
        run_bedtools_intersect(coords_bed, annot_file, intersect_out)

        # Calculate coverage
        overlap_bp = calculate_coverage(coords_bed, annot_file)

        # Count intersections
        num_intersections = count_intersections(coords_bed, annot_file)

        # Store overall stats
        stats[annot_name] = {
            'overlap_bp': overlap_bp,
            'percent_coverage': (overlap_bp / total_bp * 100) if total_bp > 0 else 0,
            'num_intersections': num_intersections,
            'intersect_file': intersect_out
        }

        # If grouping requested, calculate grouped statistics
        if group_column:
            grouped_stats = calculate_grouped_stats(intersect_out, group_column, total_bp)
            stats[annot_name]['grouped'] = grouped_stats

    return stats


def write_summary(offref_stats, onref_stats, total_offref_bp, total_onref_bp, output_file, grouped_output_file=None):
    """
    Write summary statistics to TSV file.

    Args:
        offref_stats: Dictionary of statistics for off-reference coordinates
        onref_stats: Dictionary of statistics for on-reference coordinates
        total_offref_bp: Total base pairs in off-reference regions
        total_onref_bp: Total base pairs in on-reference regions
        output_file: Output TSV file path
        grouped_output_file: Optional path for grouped statistics output
    """
    # Write overall summary
    with open(output_file, 'w') as out:
        out.write('annotation\tcoord_type\toverlap_bp\tpercent_coverage\tnum_intersections\n')

        # Write off-reference stats
        for annot_name in sorted(offref_stats.keys()):
            s = offref_stats[annot_name]
            out.write(f"{annot_name}\toff-reference\t{s['overlap_bp']}\t{s['percent_coverage']:.2f}\t{s['num_intersections']}\n")

        # Write on-reference stats
        for annot_name in sorted(onref_stats.keys()):
            s = onref_stats[annot_name]
            out.write(f"{annot_name}\ton-reference\t{s['overlap_bp']}\t{s['percent_coverage']:.2f}\t{s['num_intersections']}\n")

    sys.stderr.write(f"\nSummary written to {output_file}\n")

    # Write grouped summary if available
    if grouped_output_file:
        with open(grouped_output_file, 'w') as out:
            out.write('annotation\tcoord_type\tgroup\toverlap_bp\tpercent_coverage\tnum_intersections\n')

            # Write off-reference grouped stats
            for annot_name in sorted(offref_stats.keys()):
                if 'grouped' in offref_stats[annot_name]:
                    for group_value, stats in sorted(offref_stats[annot_name]['grouped'].items()):
                        out.write(f"{annot_name}\toff-reference\t{group_value}\t"
                                f"{stats['overlap_bp']}\t{stats['percent_coverage']:.2f}\t"
                                f"{stats['num_intersections']}\n")

            # Write on-reference grouped stats
            for annot_name in sorted(onref_stats.keys()):
                if 'grouped' in onref_stats[annot_name]:
                    for group_value, stats in sorted(onref_stats[annot_name]['grouped'].items()):
                        out.write(f"{annot_name}\ton-reference\t{group_value}\t"
                                f"{stats['overlap_bp']}\t{stats['percent_coverage']:.2f}\t"
                                f"{stats['num_intersections']}\n")

        sys.stderr.write(f"Grouped summary written to {grouped_output_file}\n")


def extract_coordinates_per_segment(nesting_file, offref_bed, onref_bed):
    """
    Extract coordinates from nesting TSV as BED4 (with augref_path as name field).

    Args:
        nesting_file: Path to nesting.tsv file with columns:
                      offref_contig, offref_start, offref_end, node, onref_contig, onref_start, onref_end
        offref_bed: Output BED4 file for off-reference coordinates
        onref_bed: Output BED4 file for on-reference coordinates

    Returns:
        Tuple of (total_offref_bp, total_onref_bp, num_skipped, segments)
        where segments is a list of (augref_path, source_len, ref_len) tuples
    """
    total_offref_bp = 0
    total_onref_bp = 0
    num_skipped = 0
    segments = []

    with open(nesting_file) as infile, \
         open(offref_bed, 'w') as offref_out, \
         open(onref_bed, 'w') as onref_out:

        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) < 7:
                continue

            offref_contig = fields[0]
            offref_start = int(fields[1])
            offref_end = int(fields[2])
            augref_path = fields[3]
            onref_contig = fields[4]
            onref_start = int(fields[5])
            onref_end = int(fields[6])

            if offref_start < 0 or offref_end < 0 or onref_start < 0 or onref_end < 0:
                num_skipped += 1
                continue
            if offref_start >= offref_end or onref_start >= onref_end:
                num_skipped += 1
                continue

            offref_out.write(f'{offref_contig}\t{offref_start}\t{offref_end}\t{augref_path}\n')
            onref_out.write(f'{onref_contig}\t{onref_start}\t{onref_end}\t{augref_path}\n')

            source_len = offref_end - offref_start
            ref_len = onref_end - onref_start
            total_offref_bp += source_len
            total_onref_bp += ref_len
            segments.append((augref_path, source_len, ref_len))

    return total_offref_bp, total_onref_bp, num_skipped, segments


def calculate_per_segment_coverage(coords_bed, annot_bed, group_column=None):
    """
    Calculate per-segment overlap with an annotation BED using bedtools intersect -wao.

    Args:
        coords_bed: BED4 file with coordinates (col 4 = augref_path)
        annot_bed: Annotation BED file
        group_column: 1-indexed column in annotation BED to extract class from (e.g., 6 for RepeatMasker)

    Returns:
        dict: augref_path -> {class -> overlap_bp}
    """
    cmd = ['bedtools', 'intersect', '-a', coords_bed, '-b', annot_bed, '-wao']
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        sys.stderr.write(f"bedtools intersect -wao failed: {result.stderr}\n")
        return {}

    coverage = defaultdict(lambda: defaultdict(int))

    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        # BED4 input: chrom, start, end, name, then annotation fields, then overlap_bp at end
        augref_path = fields[3]
        overlap_bp = int(fields[-1])

        if overlap_bp > 0 and group_column is not None:
            # Annotation fields start at index 4 (after the BED4 input fields)
            annot_col_idx = 4 + (group_column - 1)
            if annot_col_idx < len(fields) - 1:
                annot_class = fields[annot_col_idx]
                coverage[augref_path][annot_class] += overlap_bp
            else:
                coverage[augref_path]['unknown'] += overlap_bp
        elif overlap_bp > 0:
            coverage[augref_path]['_total'] += overlap_bp

    return dict(coverage)


def write_per_segment_summary(segments, source_results, ref_results, annot_names, annot_group_columns, output_file):
    """
    Write per-segment annotation overlap summary TSV.

    Args:
        segments: list of (augref_path, source_len, ref_len) tuples
        source_results: list of dicts (one per annotation), each: augref_path -> {class -> overlap_bp}
        ref_results: list of dicts (one per annotation), each: augref_path -> {class -> overlap_bp}
        annot_names: list of annotation display names
        annot_group_columns: list of group_column values (None or int) per annotation
        output_file: output TSV path
    """
    with open(output_file, 'w') as out:
        out.write('augref_path\tsource_len\tref_len\tannotation\tannotation_class\t'
                  'source_overlap_bp\tsource_overlap_frac\tref_overlap_bp\tref_overlap_frac\n')

        for augref_path, source_len, ref_len in segments:
            for i, annot_name in enumerate(annot_names):
                src_cov = source_results[i].get(augref_path, {})
                ref_cov = ref_results[i].get(augref_path, {})
                group_col = annot_group_columns[i]

                if group_col is not None:
                    # Grouped annotation (e.g., repeats): one row per class
                    all_classes = set(src_cov.keys()) | set(ref_cov.keys())
                    if not all_classes:
                        # No overlap at all — single zero row
                        src_frac = 0.0
                        ref_frac = 0.0
                        out.write(f'{augref_path}\t{source_len}\t{ref_len}\t{annot_name}\t{annot_name}\t'
                                  f'0\t{src_frac:.4f}\t0\t{ref_frac:.4f}\n')
                    else:
                        for cls in sorted(all_classes):
                            src_bp = src_cov.get(cls, 0)
                            ref_bp = ref_cov.get(cls, 0)
                            src_frac = min(src_bp / source_len, 1.0) if source_len > 0 else 0.0
                            ref_frac = min(ref_bp / ref_len, 1.0) if ref_len > 0 else 0.0
                            out.write(f'{augref_path}\t{source_len}\t{ref_len}\t{annot_name}\t{cls}\t'
                                      f'{src_bp}\t{src_frac:.4f}\t{ref_bp}\t{ref_frac:.4f}\n')
                else:
                    # Ungrouped annotation: single row, annotation_class = annotation name
                    src_bp = sum(src_cov.values())
                    ref_bp = sum(ref_cov.values())
                    src_frac = min(src_bp / source_len, 1.0) if source_len > 0 else 0.0
                    ref_frac = min(ref_bp / ref_len, 1.0) if ref_len > 0 else 0.0
                    out.write(f'{augref_path}\t{source_len}\t{ref_len}\t{annot_name}\t{annot_name}\t'
                              f'{src_bp}\t{src_frac:.4f}\t{ref_bp}\t{ref_frac:.4f}\n')

    sys.stderr.write(f"Per-segment summary written to {output_file}\n")


def main(command_line=None):
    parser = argparse.ArgumentParser(
        description='Intersect both off-reference and on-reference coordinates with annotation BED files. Requires: bedtools'
    )
    parser.add_argument('nesting_file',
                        help='Input nesting TSV file with coordinate mappings (cols 1-3: off-ref, cols 5-7: on-ref)')
    parser.add_argument('annotation_files', nargs='+',
                        help='One or more annotation BED files to intersect')
    parser.add_argument('--output-dir', '-o', default='intersect_results',
                        help='Output directory for results (default: intersect_results)')
    parser.add_argument('--summary', '-s', default='intersection_summary.tsv',
                        help='Summary statistics output file (default: intersection_summary.tsv)')
    parser.add_argument('--keep-coords-bed', action='store_true',
                        help='Keep the intermediate coordinate BED files')
    parser.add_argument('--offref-only', action='store_true',
                        help='Only process off-reference coordinates')
    parser.add_argument('--onref-only', action='store_true',
                        help='Only process on-reference coordinates')
    parser.add_argument('--group-by-column', '-g', type=int,
                        help='Column number (1-indexed) in annotation files to group statistics by (e.g., 7 for RepeatMasker class)')
    parser.add_argument('--grouped-summary', default='intersection_grouped.tsv',
                        help='Grouped statistics output file (default: intersection_grouped.tsv)')
    parser.add_argument('--per-segment', action='store_true',
                        help='Enable per-segment output mode')
    parser.add_argument('--per-segment-output', default=None,
                        help='Output path for per-segment TSV (default: per_segment_annotations.tsv in output dir)')
    parser.add_argument('--annotation-names', nargs='+', default=None,
                        help='Clean display names for each annotation file (same order as positional args)')

    options = parser.parse_args(command_line)

    check_dependencies()

    # Create output directory
    os.makedirs(options.output_dir, exist_ok=True)

    # Resolve annotation display names
    annot_names = options.annotation_names
    if annot_names is None:
        annot_names = [os.path.splitext(os.path.basename(f))[0] for f in options.annotation_files]
    if len(annot_names) != len(options.annotation_files):
        sys.exit(f"Error: {len(annot_names)} annotation names provided for {len(options.annotation_files)} annotation files")

    # Build per-annotation group_column list: only the annotation whose name is "repeats" gets grouping
    annot_group_columns = []
    for name in annot_names:
        if name == 'repeats' and options.group_by_column:
            annot_group_columns.append(options.group_by_column)
        else:
            annot_group_columns.append(None)

    # Determine which coordinate types to process
    process_offref = not options.onref_only
    process_onref = not options.offref_only

    # Extract coordinates to BED format
    sys.stderr.write(f"Extracting coordinates from {options.nesting_file}...\n")

    if options.keep_coords_bed:
        offref_bed = os.path.join(options.output_dir, 'offref_coords.bed')
        onref_bed = os.path.join(options.output_dir, 'onref_coords.bed')
    else:
        # Use temporary files
        offref_tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        onref_tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        offref_bed = offref_tmp.name
        onref_bed = onref_tmp.name
        offref_tmp.close()
        onref_tmp.close()

    if options.per_segment:
        total_offref_bp, total_onref_bp, num_skipped, segments = \
            extract_coordinates_per_segment(options.nesting_file, offref_bed, onref_bed)
    else:
        total_offref_bp, total_onref_bp, num_skipped = \
            extract_coordinates(options.nesting_file, offref_bed, onref_bed)

    sys.stderr.write(f"Total off-reference bp: {total_offref_bp:,}\n")
    sys.stderr.write(f"Total on-reference bp: {total_onref_bp:,}\n")
    if num_skipped > 0:
        sys.stderr.write(f"Warning: Skipped {num_skipped:,} records with invalid coordinates (negative or zero-length)\n")
    sys.stderr.write("\n")

    # Per-segment mode
    if options.per_segment:
        per_seg_output = options.per_segment_output
        if per_seg_output is None:
            per_seg_output = os.path.join(options.output_dir, 'per_segment_annotations.tsv')

        source_results = []
        ref_results = []

        for i, annot_file in enumerate(options.annotation_files):
            sys.stderr.write(f"  Per-segment intersect: {annot_names[i]}...\n")
            gc = annot_group_columns[i]
            source_results.append(calculate_per_segment_coverage(offref_bed, annot_file, group_column=gc))
            ref_results.append(calculate_per_segment_coverage(onref_bed, annot_file, group_column=gc))

        write_per_segment_summary(segments, source_results, ref_results, annot_names, annot_group_columns, per_seg_output)

    # Aggregate mode (only when not in per-segment mode, which uses BED4)
    offref_stats = {}
    onref_stats = {}

    if not options.per_segment:
        if process_offref:
            sys.stderr.write("Processing off-reference coordinates:\n")
            offref_stats = process_annotations(
                offref_bed, options.annotation_files, options.output_dir, total_offref_bp, 'offref',
                group_column=options.group_by_column
            )

        if process_onref:
            sys.stderr.write("\nProcessing on-reference coordinates:\n")
            onref_stats = process_annotations(
                onref_bed, options.annotation_files, options.output_dir, total_onref_bp, 'onref',
                group_column=options.group_by_column
            )

        # Write summary
        summary_path = os.path.join(options.output_dir, options.summary)
        grouped_summary_path = os.path.join(options.output_dir, options.grouped_summary) if options.group_by_column else None
        write_summary(offref_stats, onref_stats, total_offref_bp, total_onref_bp, summary_path, grouped_summary_path)

    # Cleanup
    if not options.keep_coords_bed:
        if os.path.exists(offref_bed):
            os.remove(offref_bed)
        if os.path.exists(onref_bed):
            os.remove(onref_bed)

    sys.stderr.write("\nDone!\n")

    if not options.per_segment:
        # Print summary to stdout as well
        print(f"\n{'Type':<15} {'Annotation':<35} {'Overlap (bp)':<15} {'Coverage %':<12} {'Intersections':<12}")
        print('-' * 95)

        if process_offref:
            for annot_name in sorted(offref_stats.keys()):
                s = offref_stats[annot_name]
                print(f"{'Off-reference':<15} {annot_name:<35} {s['overlap_bp']:<15,} {s['percent_coverage']:<12.2f} {s['num_intersections']:<12,}")

        if process_onref:
            for annot_name in sorted(onref_stats.keys()):
                s = onref_stats[annot_name]
                print(f"{'On-reference':<15} {annot_name:<35} {s['overlap_bp']:<15,} {s['percent_coverage']:<12.2f} {s['num_intersections']:<12,}")

        # Print grouped stats if available
        if options.group_by_column:
            print(f"\n\nGrouped by column {options.group_by_column}:")
            print(f"{'Type':<15} {'Annotation':<30} {'Group':<20} {'Overlap (bp)':<15} {'Coverage %':<12} {'Intersections':<12}")
            print('-' * 110)

            if process_offref:
                for annot_name in sorted(offref_stats.keys()):
                    if 'grouped' in offref_stats[annot_name]:
                        for group_value, gstats in sorted(offref_stats[annot_name]['grouped'].items()):
                            print(f"{'Off-reference':<15} {annot_name:<30} {group_value:<20} {gstats['overlap_bp']:<15,} "
                                  f"{gstats['percent_coverage']:<12.2f} {gstats['num_intersections']:<12,}")

            if process_onref:
                for annot_name in sorted(onref_stats.keys()):
                    if 'grouped' in onref_stats[annot_name]:
                        for group_value, gstats in sorted(onref_stats[annot_name]['grouped'].items()):
                            print(f"{'On-reference':<15} {annot_name:<30} {group_value:<20} {gstats['overlap_bp']:<15,} "
                                  f"{gstats['percent_coverage']:<12.2f} {gstats['num_intersections']:<12,}")


if __name__ == '__main__':
    main()
