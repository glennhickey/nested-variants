#!/usr/bin/env python3
"""
Download HPRC annotation BED files (RepeatMasker, Segmental Duplications, CAT genes)
from AWS S3 buckets using index CSVs from the HPRC GitHub repository.

Optionally adds GRCh38 reference annotations with appropriate naming prefix.
"""

import os
import sys
import subprocess
import argparse

# HPRC index URLs
RM_IDX_URL = 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/repeat_masker/repeat_masker_bed_pre_release_v0.2.index.csv'
SD_IDX_URL = 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/segdups/segdups_v1.0.csv'
CAT_IDX_URL = 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/cat/cat_genes_v1.index.csv'


def check_dependencies():
    """Check that required external tools are available."""
    missing = []
    for tool in ['wget', 'aws', 'parallel', 'bedtools', 'dos2unix']:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing.append(tool)
    if missing:
        sys.exit(f"Error: missing required tools: {', '.join(missing)}")


def run(cmd, shell=False):
    """Run a command, raising an error on failure."""
    subprocess.run(cmd, shell=shell, check=True)


def download_from_table(table_url, column, threads, out_annot_path, out_dir, gff=False):
    """
    Download annotation files listed in a CSV index table.

    Args:
        table_url: URL to the CSV index file
        column: 1-indexed column number containing S3 paths
        threads: Number of parallel downloads
        out_annot_path: Output BED file path
        out_dir: Output directory
        gff: If True, extract BED coordinates from GFF format
    """
    full_out_path = os.path.join(out_dir, out_annot_path)

    if os.path.isfile(full_out_path):
        sys.stderr.write(f'{full_out_path} exists, skipping re-gen\n')
        return full_out_path

    index_file = os.path.join(out_dir, os.path.basename(table_url))

    # Download the index CSV
    run(['wget', '-q', table_url, '-O', index_file])

    # Download all annotation files in parallel
    run(f'grep s3 {index_file} | dos2unix | awk -F "," \'{{print ${column}}}\' | '
        f'parallel -j {threads} "aws s3 cp --no-sign-request --no-progress {{}} {out_dir}/"',
        shell=True)

    # Concatenate all downloaded files
    run(['rm', '-f', full_out_path])
    with open(index_file) as table_file:
        for line in table_file:
            if 's3' not in line:
                continue
            annot_basename = os.path.basename(line.split(',')[column-1].strip())
            annot_path = os.path.join(out_dir, annot_basename)

            # Use zcat for compressed files, cat otherwise
            cat_cmd = 'zcat' if annot_path.endswith('.gz') else 'cat'
            gff_cmd = ' | cut -f1,4,5 | bedtools sort | bedtools merge' if gff else ''
            run(f'{cat_cmd} {annot_path} {gff_cmd} >> {full_out_path}', shell=True)
            os.remove(annot_path)

    # Sort the output
    run(f'sort -k1,1 -k2,2n {full_out_path} > {full_out_path}.tmp', shell=True)
    os.rename(full_out_path + '.tmp', full_out_path)

    # Clean up index file
    os.remove(index_file)

    return full_out_path


def add_local(local_path, annot_path, out_path, prefix='GRCh38#0#', merge=False, max_col=None):
    """
    Combine a local BED file with downloaded annotations.

    Args:
        local_path: Path to local GRCh38 BED file
        annot_path: Path to HPRC annotation BED file
        out_path: Output path for combined file
        prefix: Prefix to add to local file contigs (default: 'GRCh38#0#')
        merge: If True, merge overlapping intervals
        max_col: Maximum number of columns to keep (None for all)
    """
    with open(out_path, 'w') as out_file:
        # Add local annotations with prefix
        with open(local_path, 'r') as local_file:
            for line in local_file:
                fields = line.strip().split()
                if max_col:
                    fields = fields[:max_col]
                # Add prefix to contig name
                fields[0] = prefix + fields[0]
                out_file.write('\t'.join(fields) + '\n')

        # Add HPRC annotations as-is
        with open(annot_path, 'r') as annot_file:
            for line in annot_file:
                fields = line.strip().split()
                if max_col:
                    fields = fields[:max_col]
                out_file.write('\t'.join(fields) + '\n')

    # Sort and optionally merge
    merge_cmd = ' | bedtools merge' if merge else ''
    run(f'sort -k1,1 -k2,2n {out_path} {merge_cmd} > {out_path}.tmp', shell=True)
    os.rename(out_path + '.tmp', out_path)


def main(command_line=None):
    parser = argparse.ArgumentParser(
        description='Download HPRC annotation BED files. Requires: aws, parallel, bedtools, wget, dos2unix'
    )
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of parallel downloads (default: 8)')
    parser.add_argument('--output-dir', '-o', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--hg38-genes',
                        help='Path to GRCh38 genes BED file (optional)')
    parser.add_argument('--hg38-rm',
                        help='Path to GRCh38 RepeatMasker BED file (optional)')
    parser.add_argument('--hg38-sd',
                        help='Path to GRCh38 segmental duplications BED file (optional)')
    parser.add_argument('--skip-genes', action='store_true',
                        help='Skip downloading CAT genes')
    parser.add_argument('--skip-rm', action='store_true',
                        help='Skip downloading RepeatMasker')
    parser.add_argument('--skip-sd', action='store_true',
                        help='Skip downloading segmental duplications')

    options = parser.parse_args(command_line)

    check_dependencies()

    # Create output directory if needed
    os.makedirs(options.output_dir, exist_ok=True)

    # CAT genes
    if not options.skip_genes:
        genes_path = download_from_table(
            CAT_IDX_URL, 4, options.threads, 'hprc-v2-genes.bed', options.output_dir, gff=True
        )
        if options.hg38_genes:
            add_local(
                options.hg38_genes, genes_path,
                os.path.join(options.output_dir, 'hprc-v2-genes-grch38.bed'),
                merge=True, max_col=3
            )

    # RepeatMasker
    if not options.skip_rm:
        rm_path = download_from_table(
            RM_IDX_URL, 4, options.threads, 'hprc-v2-rm.bed', options.output_dir
        )
        if options.hg38_rm:
            add_local(
                options.hg38_rm, rm_path,
                os.path.join(options.output_dir, 'hprc-v2-rm-grch38.bed'),
                merge=False, max_col=6
            )

    # Segmental duplications
    if not options.skip_sd:
        sd_path = download_from_table(
            SD_IDX_URL, 4, options.threads, 'hprc-v2-sd.bed', options.output_dir
        )
        if options.hg38_sd:
            add_local(
                options.hg38_sd, sd_path,
                os.path.join(options.output_dir, 'hprc-v2-sd-grch38.bed'),
                merge=True, max_col=3
            )

    sys.stderr.write('Done!\n')


if __name__ == '__main__':
    main()
