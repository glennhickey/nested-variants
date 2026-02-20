#!/usr/bin/env python3
"""
Download HPRC annotation BED files (RepeatMasker, Segmental Duplications, CAT genes,
CenSat) from AWS S3 buckets using index CSVs from the HPRC GitHub repository.

GRCh38 and CHM13 reference annotations are downloaded from UCSC by default;
use --skip-grch38 or --skip-chm13 to disable.
"""

import os
import sys
import re
import shlex
import subprocess
import argparse

# Absolute path to this script (embedded in parallel commands for --merge-bed-classes)
_SCRIPT_PATH = os.path.abspath(__file__)

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

# PCLAI (PC-based Local Ancestry Inference) index URLs
PCLAI_ASM_IDX_URL = 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/pclai/pclai_v0.1_asm_coord_local_hprc_r2_v1.0.index.csv'
PCLAI_CHM13_IDX_URL = 'https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/pclai/pclai_v0.1_chm13_coord_local_hprc_r2_v1.0.index.csv'

# GIAB genome stratification BED URLs (difficult regions)
GIAB_STRAT_BASE = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6'
GIAB_HG38_DIFFICULT = f'{GIAB_STRAT_BASE}/GRCh38@all/Union/GRCh38_alldifficultregions.bed.gz'
GIAB_CHM13_DIFFICULT = f'{GIAB_STRAT_BASE}/CHM13@all/Union/CHM13_alldifficultregions.bed.gz'

# Super-population centroids in PCLAI PCA space (from ai-sandbox/hprc-pclai reference metadata)
SUPERPOP_CENTROIDS = {
    'AFR': (-1.7428, 0.2066),
    'EUR': (0.4450, -1.3142),
    'EAS': (0.7508, 1.4198),
    'SAS': (0.4768, -0.5071),
    'AMR': (0.6384, 0.3940),
}

# GRCh38 standard chromosome sizes (for placeholder PCLAI BED)
GRCH38_CHROM_SIZES = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
    'chrM': 16569,
}


def classify_superpop(pc1, pc2):
    """Classify a (PC1, PC2) coordinate to the nearest super-population centroid."""
    best_pop, best_dist = 'unknown', float('inf')
    for pop, (c1, c2) in SUPERPOP_CENTROIDS.items():
        d = (pc1 - c1) ** 2 + (pc2 - c2) ** 2
        if d < best_dist:
            best_dist = d
            best_pop = pop
    return best_pop


def convert_pclai_to_ancestry_bed(input_bed, output_bed):
    """
    Convert PCLAI BED9 to 6-column ancestry BED.

    Input name field: 'HG00097/h1/chr1_w0002_(0.440,-1.405)'
    Output: chrom, start, end, window_name, score, superpopulation
    """
    if os.path.isfile(output_bed):
        sys.stderr.write(f'  {output_bed} exists, skipping conversion\n')
        return output_bed

    pc_re = re.compile(r'\(([^,]+),([^)]+)\)')
    # Extract short window name: 'HG00097/h1/chr1_w0002_(...)' -> 'chr1_w0002'
    win_re = re.compile(r'[^/]+/[^/]+/(\S+?)_\(')

    sys.stderr.write(f'  Converting {input_bed} to ancestry BED...\n')
    with open(input_bed) as fin, open(output_bed, 'w') as fout:
        for line in fin:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            chrom, start, end, name, score = fields[0], fields[1], fields[2], fields[3], fields[4]

            # Extract PC coordinates
            m = pc_re.search(name)
            if m:
                pc1, pc2 = float(m.group(1)), float(m.group(2))
                pop = classify_superpop(pc1, pc2)
            else:
                pop = 'unknown'

            # Extract short window name
            wm = win_re.search(name)
            win_name = wm.group(1) if wm else name

            fout.write(f'{chrom}\t{start}\t{end}\t{win_name}\t{score}\t{pop}\n')

    return output_bed


def generate_grch38_pclai_placeholder(out_bed):
    """Generate a placeholder PCLAI BED for GRCh38 with 'unknown' ancestry."""
    if os.path.isfile(out_bed):
        sys.stderr.write(f'  {out_bed} exists, skipping generation\n')
        return out_bed
    sys.stderr.write(f'  Generating GRCh38 PCLAI placeholder...\n')
    with open(out_bed, 'w') as f:
        for chrom in sorted(GRCH38_CHROM_SIZES.keys()):
            size = GRCH38_CHROM_SIZES[chrom]
            f.write(f'{chrom}\t0\t{size}\t{chrom}_placeholder\t0\tunknown\n')
    return out_bed


def _merge_bed_classes_stream(fin, fout, class_col):
    """
    Core per-class merge: read sorted BED lines from fin, write merged to fout.

    Input must be sorted by chrom, start. Within each class, overlapping or
    adjacent intervals are merged.  When intervals merge, the BED fields
    (name, score, etc.) from the first interval in the group are kept; only
    positions 1 and 2 (start, end) in the ``fields`` list are mutated to
    reflect the merged coordinates.

    Output within each chromosome has classes interleaved in arbitrary order
    (dict iteration order); callers must pipe through ``sort -k1,1 -k2,2n``
    if downstream needs globally sorted BED.
    """
    col_idx = class_col - 1
    # Track the current (open) interval per class.  Values are mutable lists:
    #   [chrom, start_int, end_int, fields_list]
    # where fields_list is the original split line (from the first interval
    # in the merge group).  start/end are tracked as ints separately from
    # fields to avoid repeated str↔int conversions during merging; fields[1]
    # and fields[2] are written back only at flush time.
    active = {}
    prev_chrom = None

    def flush():
        """Write all pending intervals for prev_chrom.  Order is arbitrary."""
        for cls, (chrom, start, end, fields) in active.items():
            fields[1] = str(start)
            fields[2] = str(end)
            fout.write('\t'.join(fields) + '\n')
        active.clear()

    for line in fin:
        fields = line.rstrip('\n').split('\t')
        if len(fields) <= col_idx:
            continue

        cls = fields[col_idx]
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])

        if chrom != prev_chrom:
            if prev_chrom is not None:
                flush()
            prev_chrom = chrom

        if cls in active:
            a_chrom, a_start, a_end, a_fields = active[cls]
            if start <= a_end:
                active[cls] = [chrom, a_start, max(a_end, end), a_fields]
            else:
                a_fields[1] = str(a_start)
                a_fields[2] = str(a_end)
                fout.write('\t'.join(a_fields) + '\n')
                active[cls] = [chrom, start, end, fields]
        else:
            active[cls] = [chrom, start, end, fields]

    flush()


def merge_bed_classes_stdin(class_col):
    """Pipe mode: read sorted BED from stdin, merge per class, write to stdout."""
    _merge_bed_classes_stream(sys.stdin, sys.stdout, class_col)


def merge_bed_classes_file(input_bed, output_bed, class_col):
    """
    Merge overlapping intervals per class in a sorted BED file.

    Output is re-sorted by chrom, start (the merge interleaves classes).
    """
    tmp = output_bed + '.unsorted'
    try:
        with open(input_bed) as fin, open(tmp, 'w') as fout:
            _merge_bed_classes_stream(fin, fout, class_col)
        run(f"LC_ALL=C sort -k1,1 -k2,2n '{tmp}' -o '{output_bed}'", shell=True)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


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
                        test_sample=None, max_col=None, reformat_cmd=None,
                        merge_class_col=None):
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
        reformat_cmd: Shell pipe command to reformat columns (inserted before cut/sort)
        merge_class_col: If set, merge overlapping intervals per class (1-indexed column)
                         after sorting each sample's file
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

    # Sort/convert each file individually in parallel, then concatenate in
    # contig-name order. Each file's contigs share a sample prefix (e.g.
    # HG00408#1#), so ordering files by first contig gives global sort order
    # without an expensive sort on the full concatenation.
    sorted_files = []
    cmd_file = os.path.join(out_dir, f'{out_annot_path}.sort_cmds')
    with open(index_file) as table_file, open(cmd_file, 'w') as cmds:
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
            reformat = f' | {reformat_cmd}' if reformat_cmd and not gff else ''
            cut_cmd = f' | cut -f1-{max_col}' if max_col and not gff else ''
            sort_cmd = '' if gff else ' | LC_ALL=C sort -k1,1 -k2,2n'
            if merge_class_col and not gff:
                merge_cmd = (f' | {shlex.quote(sys.executable)} {shlex.quote(_SCRIPT_PATH)}'
                             f' --merge-bed-classes {merge_class_col}'
                             f' | LC_ALL=C sort -k1,1 -k2,2n')
            else:
                merge_cmd = ''

            cmds.write(f'{cat_cmd} {annot_path}{strip_cmd}{gff_cmd}{reformat}{cut_cmd}{sort_cmd}{merge_cmd}'
                       f' > {sorted_path} && rm {annot_path}\n')
            sorted_files.append(sorted_path)

    run(f'parallel -j {threads} < {cmd_file}', shell=True)
    os.remove(cmd_file)

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
    parser.add_argument('--skip-pclai', action='store_true',
                        help='Skip downloading PCLAI (local ancestry) annotations')
    parser.add_argument('--skip-giab', action='store_true',
                        help='Skip downloading GIAB difficult regions stratification')
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
            merge_bed_classes_file(hg38['rm'], hg38['rm'] + '.merged', class_col=6)
            os.rename(hg38['rm'] + '.merged', hg38['rm'])
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
            merge_bed_classes_file(chm13['rm'], chm13['rm'] + '.merged', class_col=6)
            os.rename(chm13['rm'] + '.merged', chm13['rm'])
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
        # HPRC RM BED format: chrom, start, end, name, score, strand, class, family, ...
        # Reformat to match reference RM BEDs: chrom, start, end, name, score, class/family
        rm_reformat = r"""awk -F'\t' '{OFS="\t"; print $1,$2,$3,$4,$5,$7"/"$8}'"""
        rm_path = download_from_table(
            RM_IDX_URL, 4, options.threads, 'hprc-v2-rm.bed', options.output_dir,
            test_sample=options.test or None, reformat_cmd=rm_reformat,
            merge_class_col=6
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

        # Pre-compute class-agnostic merged BED for fast _total intersection
        total_path = os.path.splitext(current)[0] + '.total.bed'
        if not os.path.isfile(total_path):
            sys.stderr.write(f'  Pre-computing class-agnostic merged RM BED: {total_path}\n')
            tmp_total = total_path + '.tmp'
            run(f"cut -f1-3 '{current}' | bedtools merge > '{tmp_total}'", shell=True)
            os.rename(tmp_total, total_path)

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

    # PCLAI (local ancestry)
    if not options.skip_pclai:
        sys.stderr.write('\nDownloading PCLAI local ancestry annotations...\n')

        # Download HPRC sample PCLAI (asm_coord — assembly-native contigs)
        raw_pclai = download_from_table(
            PCLAI_ASM_IDX_URL, 4, options.threads, 'hprc-v2-pclai-raw.bed', options.output_dir,
            test_sample=options.test or None
        )
        pclai_path = os.path.join(options.output_dir, 'hprc-v2-pclai.bed')
        convert_pclai_to_ancestry_bed(raw_pclai, pclai_path)
        sys.stderr.write('  Merging PCLAI per-ancestry overlaps...\n')
        merge_bed_classes_file(pclai_path, pclai_path + '.merged', class_col=6)
        os.rename(pclai_path + '.merged', pclai_path)
        intermediates.append(raw_pclai)

        # CHM13 PCLAI from chm13_coord index
        if not options.skip_chm13:
            raw_chm13 = download_from_table(
                PCLAI_CHM13_IDX_URL, 4, options.threads, 'chm13-pclai-raw.bed', options.output_dir,
                test_sample='CHM13'
            )
            chm13_pclai = os.path.join(options.output_dir, 'chm13-pclai.bed')
            convert_pclai_to_ancestry_bed(raw_chm13, chm13_pclai)
            merge_bed_classes_file(chm13_pclai, chm13_pclai + '.merged', class_col=6)
            os.rename(chm13_pclai + '.merged', chm13_pclai)
            chm13['pclai'] = chm13_pclai
            intermediates += [raw_chm13, chm13_pclai]

        # GRCh38 placeholder (full chromosomes labeled 'unknown')
        if not options.skip_grch38:
            hg38['pclai'] = generate_grch38_pclai_placeholder(
                os.path.join(options.output_dir, 'grch38-pclai.bed'))
            intermediates.append(hg38['pclai'])

        # Chain: HPRC → +GRCh38 → +CHM13
        current = pclai_path
        suffix = ''
        for prefix, ref_beds, ref_name in [('GRCh38#0#', hg38, 'grch38'),
                                            ('CHM13#0#', chm13, 'chm13')]:
            if 'pclai' in ref_beds:
                intermediates.append(current)
                suffix += f'-{ref_name}'
                out = os.path.join(options.output_dir, f'hprc-v2-pclai{suffix}.bed')
                add_local(ref_beds['pclai'], current, out, prefix=prefix, max_col=6)
                current = out

        # Pre-compute class-agnostic merged BED for fast _total intersection
        total_path = os.path.splitext(current)[0] + '.total.bed'
        if not os.path.isfile(total_path):
            sys.stderr.write(f'  Pre-computing class-agnostic merged PCLAI BED: {total_path}\n')
            tmp_total = total_path + '.tmp'
            run(f"cut -f1-3 '{current}' | bedtools merge > '{tmp_total}'", shell=True)
            os.rename(tmp_total, total_path)

    # GIAB difficult regions (reference-only, no HPRC assembly-level data)
    if not options.skip_giab:
        sys.stderr.write('\nDownloading GIAB stratification BEDs...\n')
        giab_path = os.path.join(options.output_dir, 'hprc-v2-giab-difficult.bed')
        if not os.path.isfile(giab_path):
            giab_parts = []
            # Download GRCh38 (respects --skip-grch38)
            if not options.skip_grch38:
                hg38_raw = download_ucsc_file(GIAB_HG38_DIFFICULT, options.output_dir)
                hg38_giab = os.path.join(options.output_dir, 'grch38-giab-difficult.bed')
                run(f"zcat '{hg38_raw}' | cut -f1-3 | bedtools sort | bedtools merge > '{hg38_giab}'", shell=True)
                giab_parts.append(('GRCh38#0#', hg38_giab))
                intermediates += [hg38_raw, hg38_giab]
            # Download CHM13 (respects --skip-chm13)
            if not options.skip_chm13:
                chm13_raw = download_ucsc_file(GIAB_CHM13_DIFFICULT, options.output_dir)
                chm13_giab = os.path.join(options.output_dir, 'chm13-giab-difficult.bed')
                run(f"zcat '{chm13_raw}' | cut -f1-3 | bedtools sort | bedtools merge > '{chm13_giab}'", shell=True)
                giab_parts.append(('CHM13#0#', chm13_giab))
                intermediates += [chm13_raw, chm13_giab]
            # Prefix and concatenate
            first = True
            for prefix, bed in giab_parts:
                redirect = '>' if first else '>>'
                run(f"awk -v p='{prefix}' 'BEGIN{{OFS=\"\\t\"}} {{$1=p$1; print}}' '{bed}'"
                    f" {redirect} '{giab_path}'", shell=True)
                first = False
            if not giab_parts:
                # Both references skipped — create empty file
                open(giab_path, 'w').close()
        sys.stderr.write(f'  GIAB difficult regions: {giab_path}\n')

    # Clean up intermediate files
    for path in intermediates:
        if os.path.isfile(path):
            os.remove(path)
            sys.stderr.write(f'  Removed intermediate: {os.path.basename(path)}\n')

    sys.stderr.write('Done!\n')


if __name__ == '__main__':
    # Pipe mode: merge per-class overlapping intervals (used internally by parallel download)
    if '--merge-bed-classes' in sys.argv:
        idx = sys.argv.index('--merge-bed-classes')
        if idx + 1 >= len(sys.argv):
            sys.exit('Error: --merge-bed-classes requires a column number argument')
        merge_bed_classes_stdin(int(sys.argv[idx + 1]))
        sys.exit(0)
    main()
