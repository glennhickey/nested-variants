#!/usr/bin/env python3
"""
Download HPRC annotation BED files (RepeatMasker, Segmental Duplications, CAT genes,
CenSat) from AWS S3 buckets using index CSVs from the HPRC GitHub repository.

GRCh38 and CHM13 reference annotations are downloaded from UCSC by default;
use --skip-grch38 or --skip-chm13 to disable.
"""

import os
import sys
import subprocess
import argparse

# HPRC index URLs
RM_IDX_URL = 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/repeat_masker/repeat_masker_bed_hprc_r2_v1.0.index.csv'
SD_IDX_URL = 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/segdups/segdups_hprc_r2_v1.1.index.csv'
CAT_IDX_URL = 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/cat/cat_genes_hprc_r2_v1.2.index.csv'
CENSAT_IDX_URL = 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/censat/censat_hprc_r2_v1.0.index.csv'

# UCSC download URLs for reference annotations
UCSC_HG38_RMSK = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz'
UCSC_HG38_SEGDUPS = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz'
UCSC_HG38_GENES = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz'
UCSC_HS1_RMSK = 'https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.repeatMasker.out.gz'
UCSC_HS1_SEGDUPS = 'https://hgdownload.soe.ucsc.edu/gbdb/hs1/sedefSegDups/sedefSegDups.bb'
UCSC_HS1_GENES = 'https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz'

# CenSat: CHM13 has full satellite classification; GRCh38 has only centromere boundaries
# (centromeric regions are mostly gaps in GRCh38, so detailed CenSat is unavailable)
UCSC_HG38_CENTROMERES = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz'
UCSC_HS1_CENSAT = 'https://hgdownload.soe.ucsc.edu/gbdb/hs1/censat/censat.bb'


def check_dependencies(need_bigbed=False):
    """Check that required external tools are available."""
    missing = []
    for tool in ['wget', 'aws', 'parallel', 'bedtools', 'dos2unix']:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing.append(tool)
    if need_bigbed:
        if subprocess.run(['which', 'bigBedToBed'], capture_output=True).returncode != 0:
            missing.append('bigBedToBed')
    if missing:
        sys.exit(f"Error: missing required tools: {', '.join(missing)}")


def run(cmd, shell=False):
    """Run a command, raising an error on failure."""
    subprocess.run(cmd, shell=shell, check=True)


def download_ucsc_file(url, out_dir):
    """Download a file from UCSC, return local path. Skip if exists."""
    local_path = os.path.join(out_dir, os.path.basename(url))
    if os.path.isfile(local_path):
        sys.stderr.write(f'  {local_path} exists, skipping download\n')
        return local_path
    sys.stderr.write(f'  Downloading {url}...\n')
    run(['wget', '-q', url, '-O', local_path])
    return local_path


def convert_ucsc_rmsk_to_bed(rmsk_gz, out_bed):
    """
    Convert UCSC rmsk.txt.gz table to 6-column repeat BED.

    UCSC rmsk columns (tab-separated):
      bin(0) swScore(1) milliDiv(2) milliDel(3) milliIns(4) genoName(5)
      genoStart(6) genoEnd(7) genoLeft(8) strand(9) repName(10)
      repClass(11) repFamily(12) ...

    Output: chrom, start, end, repName, swScore, repClass/repFamily
    """
    if os.path.isfile(out_bed):
        sys.stderr.write(f'  {out_bed} exists, skipping conversion\n')
        return out_bed
    sys.stderr.write(f'  Converting {rmsk_gz} to BED...\n')
    run(f'zcat {rmsk_gz} | awk -F"\\t" \'{{OFS="\\t"; print $6,$7,$8,$11,$2,$12"/"$13}}\' '
        f'| sort -k1,1 -k2,2n > {out_bed}', shell=True)
    return out_bed


def convert_rm_out_to_bed(rm_out_gz, out_bed):
    """
    Convert RepeatMasker .out.gz to 6-column repeat BED.

    RM .out format (whitespace-separated, 3 header lines):
      score(0) div(1) del(2) ins(3) chrom(4) begin(5, 1-based) end(6)
      left(7) strand(8) repName(9) class/family(10) ...

    Output: chrom, start(0-based), end, repName, score, class/family
    """
    if os.path.isfile(out_bed):
        sys.stderr.write(f'  {out_bed} exists, skipping conversion\n')
        return out_bed
    sys.stderr.write(f'  Converting {rm_out_gz} to BED...\n')
    run(f'zcat {rm_out_gz} | tail -n+4 | '
        f'awk \'NF>=11 {{OFS="\\t"; print $5,$6-1,$7,$10,$1,$11}}\' '
        f'| sort -k1,1 -k2,2n > {out_bed}', shell=True)
    return out_bed


def convert_ucsc_segdups_to_bed(segdups_gz, out_bed):
    """
    Convert UCSC genomicSuperDups.txt.gz table to 3-column BED.

    UCSC genomicSuperDups columns: bin(0) chrom(1) chromStart(2) chromEnd(3) ...
    """
    if os.path.isfile(out_bed):
        sys.stderr.write(f'  {out_bed} exists, skipping conversion\n')
        return out_bed
    sys.stderr.write(f'  Converting {segdups_gz} to BED...\n')
    run(f'zcat {segdups_gz} | awk -F"\\t" \'{{OFS="\\t"; print $2,$3,$4}}\' '
        f'| sort -k1,1 -k2,2n | bedtools merge > {out_bed}', shell=True)
    return out_bed


def convert_bigbed_to_bed(bb_file, out_bed):
    """Convert BigBed to 3-column merged BED using bigBedToBed."""
    if os.path.isfile(out_bed):
        sys.stderr.write(f'  {out_bed} exists, skipping conversion\n')
        return out_bed
    sys.stderr.write(f'  Converting {bb_file} to BED...\n')
    tmp_bed = bb_file + '.tmp.bed'
    run(['bigBedToBed', bb_file, tmp_bed])
    run(f'cut -f1-3 {tmp_bed} | sort -k1,1 -k2,2n | bedtools merge > {out_bed}', shell=True)
    os.remove(tmp_bed)
    return out_bed


def convert_ucsc_centromeres_to_bed(centromeres_gz, out_bed):
    """
    Convert UCSC centromeres.txt.gz table to 3-column BED.

    UCSC centromeres columns (tab-separated):
      bin(0) chrom(1) chromStart(2) chromEnd(3) name(4)
    """
    if os.path.isfile(out_bed):
        sys.stderr.write(f'  {out_bed} exists, skipping conversion\n')
        return out_bed
    sys.stderr.write(f'  Converting {centromeres_gz} to BED...\n')
    run(f'zcat {centromeres_gz} | awk -F"\\t" \'{{OFS="\\t"; print $2,$3,$4}}\' '
        f'| sort -k1,1 -k2,2n | bedtools merge > {out_bed}', shell=True)
    return out_bed


def convert_gtf_to_bed(gtf_gz, out_bed):
    """
    Convert GTF to 3-column merged gene BED.

    GTF: chrom(0) source(1) feature(2) start(3, 1-based) end(4) ...
    Filters for 'transcript' features (UCSC ncbiRefSeq GTFs lack 'gene' rows),
    converts to 0-based BED, merges overlapping.
    """
    if os.path.isfile(out_bed):
        sys.stderr.write(f'  {out_bed} exists, skipping conversion\n')
        return out_bed
    sys.stderr.write(f'  Converting {gtf_gz} to BED...\n')
    run(f'zcat {gtf_gz} | awk -F"\\t" \'$3=="transcript" {{OFS="\\t"; print $1,$4-1,$5}}\' '
        f'| sort -k1,1 -k2,2n | bedtools merge > {out_bed}', shell=True)
    return out_bed


def download_from_table(table_url, column, threads, out_annot_path, out_dir, gff=False,
                        test_sample=None, max_col=None):
    """
    Download annotation files listed in a CSV index table.

    Each per-sample file is sorted individually, then all files are concatenated
    in contig-name order (sample prefix determines global sort). This avoids an
    expensive global sort on the combined multi-GB output.

    Args:
        table_url: URL to the CSV index file
        column: 1-indexed column number containing S3 paths
        threads: Number of parallel downloads
        out_annot_path: Output BED file path
        out_dir: Output directory
        gff: If True, extract BED coordinates from GFF format
        test_sample: If set, only download this sample's rows (both haplotypes)
        max_col: Maximum number of columns to keep (None for all)
    """
    full_out_path = os.path.join(out_dir, out_annot_path)

    if os.path.isfile(full_out_path):
        sys.stderr.write(f'{full_out_path} exists, skipping re-gen\n')
        return full_out_path

    index_file = os.path.join(out_dir, os.path.basename(table_url))

    # Download the index CSV
    run(['wget', '-q', table_url, '-O', index_file])

    # Filter to test sample if requested
    if test_sample is True:
        # Auto-pick: read first sample_id from the CSV
        with open(index_file) as f:
            for line in f:
                line = line.strip()
                if not line or 's3' not in line:
                    continue
                test_sample = line.split(',')[0]
                break
        sys.stderr.write(f'  Test mode: using sample {test_sample}\n')

    if test_sample:
        # Filter index to only matching sample rows (+ header)
        filtered = index_file + '.filtered'
        with open(index_file) as f_in, open(filtered, 'w') as f_out:
            for line in f_in:
                if 's3' not in line:
                    f_out.write(line)
                elif line.startswith(test_sample + ','):
                    f_out.write(line)
        os.rename(filtered, index_file)

    # Download all annotation files in parallel
    run(f'grep s3 {index_file} | dos2unix | awk -F "," \'{{print ${column}}}\' | '
        f'parallel -j {threads} "aws s3 cp --no-sign-request --no-progress {{}} {out_dir}/"',
        shell=True)

    # Sort each file individually, then concatenate in contig-name order.
    # Each file's contigs share a sample prefix (e.g. HG00408#1#), so
    # ordering files by their first contig gives global sort order without
    # an expensive sort on the full concatenation.
    sorted_files = []
    with open(index_file) as table_file:
        for line in table_file:
            if 's3' not in line:
                continue
            annot_basename = os.path.basename(line.split(',')[column-1].strip())
            annot_path = os.path.join(out_dir, annot_basename)
            sorted_path = annot_path + '.sorted'

            cat_cmd = 'zcat' if annot_path.endswith('.gz') else 'cat'
            # Strip UCSC track/browser headers and comments
            strip_cmd = " | grep -v '^track\\|^browser\\|^#'"
            # GFF path: extract coords, sort, merge (already produces sorted output)
            gff_cmd = ' | cut -f1,4,5 | bedtools sort | bedtools merge' if gff else ''
            cut_cmd = f' | cut -f1-{max_col}' if max_col and not gff else ''
            sort_cmd = '' if gff else ' | sort -k1,1 -k2,2n'

            run(f'{cat_cmd} {annot_path}{strip_cmd}{gff_cmd}{cut_cmd}{sort_cmd} > {sorted_path}',
                shell=True)
            os.remove(annot_path)
            sorted_files.append(sorted_path)

    # Order files by first contig name for correct global sort
    def first_contig(path):
        with open(path) as f:
            line = f.readline()
            return line.split('\t')[0] if line else ''

    sorted_files.sort(key=first_contig)
    run(f'cat {" ".join(sorted_files)} > {full_out_path}', shell=True)
    for f in sorted_files:
        os.remove(f)

    # Clean up index file
    os.remove(index_file)

    return full_out_path


def add_local(local_path, annot_path, out_path, prefix='GRCh38#0#', max_col=None):
    """
    Combine a local reference BED with HPRC annotations.

    Reference contigs get prefixed (e.g. chr1 → GRCh38#0#chr1). Both inputs
    are already sorted, and reference prefixes sort before HPRC sample names
    (CHM13#0# < GRCh38#0# < HG...), so we just concatenate — no re-sort needed.
    Intervals with different contig prefixes cannot overlap, so bedtools merge
    is also unnecessary.
    """
    if os.path.isfile(out_path):
        sys.stderr.write(f'  {out_path} exists, skipping re-gen\n')
        return

    # Write prefixed reference annotations, then append HPRC annotations
    cut_cmd = f' | cut -f1-{max_col}' if max_col else ''
    run(f"awk -v p='{prefix}' 'BEGIN{{OFS=\"\\t\"}} {{$1=p$1; print}}' "
        f"{local_path}{cut_cmd} > {out_path}", shell=True)
    run(f'cat {annot_path} >> {out_path}', shell=True)


def main(command_line=None):
    parser = argparse.ArgumentParser(
        description='Download HPRC annotation BED files. Requires: aws, parallel, bedtools, wget, dos2unix'
    )
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of parallel downloads (default: 8)')
    parser.add_argument('--output-dir', '-o', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--skip-grch38', action='store_true',
                        help='Skip downloading GRCh38 reference annotations from UCSC')
    parser.add_argument('--skip-chm13', action='store_true',
                        help='Skip downloading CHM13 reference annotations from UCSC '
                             '(these require bigBedToBed for segdups and censat)')
    parser.add_argument('--skip-genes', action='store_true',
                        help='Skip downloading CAT genes')
    parser.add_argument('--skip-rm', action='store_true',
                        help='Skip downloading RepeatMasker')
    parser.add_argument('--skip-sd', action='store_true',
                        help='Skip downloading segmental duplications')
    parser.add_argument('--skip-censat', action='store_true',
                        help='Skip downloading CenSat (centromeric satellite) annotations')
    parser.add_argument('--test', action='store_true',
                        help='Test mode: download only one HPRC sample (both haplotypes)')

    options = parser.parse_args(command_line)

    need_bigbed = not options.skip_chm13 and (not options.skip_sd or not options.skip_censat)
    check_dependencies(need_bigbed=need_bigbed)

    # Create output directory if needed
    os.makedirs(options.output_dir, exist_ok=True)

    # Track intermediate files for cleanup
    intermediates = []

    # Download and convert reference annotations from UCSC
    hg38 = {}  # annotation type -> local BED path
    chm13 = {}

    if not options.skip_grch38:
        sys.stderr.write('\nDownloading GRCh38 reference annotations from UCSC...\n')
        if not options.skip_genes:
            raw = download_ucsc_file(UCSC_HG38_GENES, options.output_dir)
            hg38['genes'] = convert_gtf_to_bed(
                raw, os.path.join(options.output_dir, 'grch38-genes.bed'))
            intermediates += [raw, hg38['genes']]
        if not options.skip_rm:
            raw = download_ucsc_file(UCSC_HG38_RMSK, options.output_dir)
            hg38['rm'] = convert_ucsc_rmsk_to_bed(
                raw, os.path.join(options.output_dir, 'grch38-rm.bed'))
            intermediates += [raw, hg38['rm']]
        if not options.skip_sd:
            raw = download_ucsc_file(UCSC_HG38_SEGDUPS, options.output_dir)
            hg38['sd'] = convert_ucsc_segdups_to_bed(
                raw, os.path.join(options.output_dir, 'grch38-sd.bed'))
            intermediates += [raw, hg38['sd']]
        if not options.skip_censat:
            raw = download_ucsc_file(UCSC_HG38_CENTROMERES, options.output_dir)
            hg38['censat'] = convert_ucsc_centromeres_to_bed(
                raw, os.path.join(options.output_dir, 'grch38-censat.bed'))
            intermediates += [raw, hg38['censat']]

    if not options.skip_chm13:
        sys.stderr.write('\nDownloading CHM13 reference annotations from UCSC...\n')
        if not options.skip_genes:
            raw = download_ucsc_file(UCSC_HS1_GENES, options.output_dir)
            chm13['genes'] = convert_gtf_to_bed(
                raw, os.path.join(options.output_dir, 'chm13-genes.bed'))
            intermediates += [raw, chm13['genes']]
        if not options.skip_rm:
            raw = download_ucsc_file(UCSC_HS1_RMSK, options.output_dir)
            chm13['rm'] = convert_rm_out_to_bed(
                raw, os.path.join(options.output_dir, 'chm13-rm.bed'))
            intermediates += [raw, chm13['rm']]
        if not options.skip_sd:
            raw = download_ucsc_file(UCSC_HS1_SEGDUPS, options.output_dir)
            chm13['sd'] = convert_bigbed_to_bed(
                raw, os.path.join(options.output_dir, 'chm13-sd.bed'))
            intermediates += [raw, chm13['sd']]
        if not options.skip_censat:
            raw = download_ucsc_file(UCSC_HS1_CENSAT, options.output_dir)
            chm13['censat'] = convert_bigbed_to_bed(
                raw, os.path.join(options.output_dir, 'chm13-censat.bed'))
            intermediates += [raw, chm13['censat']]

    # CAT genes
    if not options.skip_genes:
        genes_path = download_from_table(
            CAT_IDX_URL, 4, options.threads, 'hprc-v2-genes.bed', options.output_dir, gff=True,
            test_sample=options.test or None
        )
        # Chain reference annotations: HPRC → +GRCh38 → +CHM13
        current = genes_path
        suffix = ''
        for prefix, ref_beds, ref_name in [('GRCh38#0#', hg38, 'grch38'),
                                            ('CHM13#0#', chm13, 'chm13')]:
            if 'genes' in ref_beds:
                intermediates.append(current)
                suffix += f'-{ref_name}'
                out = os.path.join(options.output_dir, f'hprc-v2-genes{suffix}.bed')
                add_local(ref_beds['genes'], current, out, prefix=prefix, max_col=3)
                current = out

    # RepeatMasker
    if not options.skip_rm:
        rm_path = download_from_table(
            RM_IDX_URL, 4, options.threads, 'hprc-v2-rm.bed', options.output_dir,
            test_sample=options.test or None, max_col=6
        )
        current = rm_path
        suffix = ''
        for prefix, ref_beds, ref_name in [('GRCh38#0#', hg38, 'grch38'),
                                            ('CHM13#0#', chm13, 'chm13')]:
            if 'rm' in ref_beds:
                intermediates.append(current)
                suffix += f'-{ref_name}'
                out = os.path.join(options.output_dir, f'hprc-v2-rm{suffix}.bed')
                add_local(ref_beds['rm'], current, out, prefix=prefix, max_col=6)
                current = out

    # Segmental duplications
    if not options.skip_sd:
        sd_path = download_from_table(
            SD_IDX_URL, 4, options.threads, 'hprc-v2-sd.bed', options.output_dir,
            test_sample=options.test or None, max_col=3
        )
        current = sd_path
        suffix = ''
        for prefix, ref_beds, ref_name in [('GRCh38#0#', hg38, 'grch38'),
                                            ('CHM13#0#', chm13, 'chm13')]:
            if 'sd' in ref_beds:
                intermediates.append(current)
                suffix += f'-{ref_name}'
                out = os.path.join(options.output_dir, f'hprc-v2-sd{suffix}.bed')
                add_local(ref_beds['sd'], current, out, prefix=prefix, max_col=3)
                current = out

    # CenSat (centromeric satellite)
    if not options.skip_censat:
        censat_path = download_from_table(
            CENSAT_IDX_URL, 4, options.threads, 'hprc-v2-censat.bed', options.output_dir,
            test_sample=options.test or None, max_col=3
        )
        current = censat_path
        suffix = ''
        for prefix, ref_beds, ref_name in [('GRCh38#0#', hg38, 'grch38'),
                                            ('CHM13#0#', chm13, 'chm13')]:
            if 'censat' in ref_beds:
                intermediates.append(current)
                suffix += f'-{ref_name}'
                out = os.path.join(options.output_dir, f'hprc-v2-censat{suffix}.bed')
                add_local(ref_beds['censat'], current, out, prefix=prefix, max_col=3)
                current = out

    # Clean up intermediate files
    for path in intermediates:
        if os.path.isfile(path):
            os.remove(path)
            sys.stderr.write(f'  Removed intermediate: {os.path.basename(path)}\n')

    sys.stderr.write('Done!\n')


if __name__ == '__main__':
    main()
