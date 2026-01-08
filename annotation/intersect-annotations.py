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
        Tuple of (total_offref_bp, total_onref_bp)
    """
    total_offref_bp = 0
    total_onref_bp = 0

    with open(nesting_file) as infile, \
         open(offref_bed, 'w') as offref_out, \
         open(onref_bed, 'w') as onref_out:

        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) < 7:
                continue

            # Extract off-reference coordinates (columns 1-3, 0-indexed: 0-2)
            offref_contig = fields[0]
            offref_start = fields[1]
            offref_end = fields[2]

            # Extract on-reference coordinates (columns 5-7, 0-indexed: 4-6)
            onref_contig = fields[4]
            onref_start = fields[5]
            onref_end = fields[6]

            offref_out.write(f'{offref_contig}\t{offref_start}\t{offref_end}\n')
            onref_out.write(f'{onref_contig}\t{onref_start}\t{onref_end}\n')

            total_offref_bp += int(offref_end) - int(offref_start)
            total_onref_bp += int(onref_end) - int(onref_start)

    return total_offref_bp, total_onref_bp


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

    options = parser.parse_args(command_line)

    check_dependencies()

    # Create output directory
    os.makedirs(options.output_dir, exist_ok=True)

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

    total_offref_bp, total_onref_bp = extract_coordinates(options.nesting_file, offref_bed, onref_bed)

    sys.stderr.write(f"Total off-reference bp: {total_offref_bp:,}\n")
    sys.stderr.write(f"Total on-reference bp: {total_onref_bp:,}\n\n")

    # Process annotations for both coordinate sets
    offref_stats = {}
    onref_stats = {}

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
