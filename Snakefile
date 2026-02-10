############################################################################
# Snakefile — nested-variants pipeline
#
# Default config is for local chr20 testing.
# Override via config.local.yaml or command line:
#   snakemake --config ref=CHM13 vg='...' all
############################################################################

configfile: "config.yaml"

REF      = config["ref"]
AUGREF   = f"augref_{REF}"
OUT_DIR  = config["out_dir"]
OUT_NAME = config["out_name"]

# Load samples from config dict
SAMPLES = list(config.get("samples", {}).keys())

# Load samples from TSV if provided
if config.get("samples_tsv"):
    import csv
    with open(config["samples_tsv"]) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            config.setdefault("samples", {})[row["sample"]] = row["reads_index"]
    SAMPLES = list(config["samples"].keys())

def map_gbz(wildcards=None):
    """GBZ used for read mapping (giraffe). Falls back to pipeline-built GBZ."""
    return config["map_gbz"] or f"{OUT_DIR}/{OUT_NAME}.gbz"

def hapl_index(wildcards=None):
    """Haplotype index. Falls back to pipeline-built .hapl."""
    return config["hapl"] or f"{OUT_DIR}/{OUT_NAME}.hapl"

# Build deconstruct-specific option flags
def decon_opts():
    opts = []
    if config.get("cluster"):
        opts.append(f"--cluster {config['cluster']}")
    if config.get("snarls"):
        opts.append(f"--snarls {config['snarls']}")
    if config.get("star_allele"):
        opts.append("--star-allele")
    return " ".join(opts)

############################################################################
# Target rules
############################################################################

rule all:
    """Default: build graph + analysis (no genotyping)"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.offref.vcf.gz",
        f"{OUT_DIR}/{OUT_NAME}.offref.png",
        f"{OUT_DIR}/{OUT_NAME}.augref-length-hist.png",
        f"{OUT_DIR}/{OUT_NAME}.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.af-spectrum.png",

rule genotype_all:
    """Genotype all samples (vg call)"""
    input:
        expand("{out}/{s}.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.variant-types.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.size-dist.png", out=OUT_DIR, s=SAMPLES),

rule deepvariant_all:
    """Run DeepVariant on all samples"""
    input:
        expand("{out}/{s}.deepvariant.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.variant-types.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.size-dist.png", out=OUT_DIR, s=SAMPLES),

rule batch:
    """Full batch: genotype + deepvariant + merge + merged plots"""
    input:
        f"{OUT_DIR}/merged.call.vcf.gz",
        f"{OUT_DIR}/merged.deepvariant.vcf.gz",
        f"{OUT_DIR}/merged.call-offref.png",
        f"{OUT_DIR}/merged.dv-offref.png",
        f"{OUT_DIR}/{OUT_NAME}.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.af-spectrum.png",
        f"{OUT_DIR}/merged.call.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.variant-types.png",
        f"{OUT_DIR}/merged.call.size-dist.png",
        f"{OUT_DIR}/merged.call.af-spectrum.png",
        f"{OUT_DIR}/merged.dv.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.variant-types.png",
        f"{OUT_DIR}/merged.dv.size-dist.png",
        f"{OUT_DIR}/merged.dv.af-spectrum.png",

############################################################################
# Graph construction rules (run once)
############################################################################

rule paths:
    """VG → GBZ via augmented reference paths"""
    input:
        ancient(config["vg"]),
    output:
        f"{OUT_DIR}/{OUT_NAME}.gbz",
        f"{OUT_DIR}/{OUT_NAME}.gfa.gz",
        f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    threads: config.get("cpus", 8)
    resources:
        mem_mb=config.get("mem_gb", 200) * 1024,
        runtime=config.get("runtime_min", 960),
    shell:
        "scripts/paths.sh"
        " --vg '{input}'"
        " --ref {REF}"
        " --out-dir {OUT_DIR}"
        " --out-name {OUT_NAME}.gfa.gz"
        " --min-augref-len {config[min_augref_len]}"
        " --cpus {threads} --mem {config[mem_gb]}gb"
        " --local"

rule deconstruct:
    """GBZ → VCF"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
    threads: config.get("cpus", 8)
    resources:
        mem_mb=config.get("mem_gb", 200) * 1024,
        runtime=config.get("runtime_min", 960),
    shell:
        "scripts/deconstruct.sh"
        " --gbz {input}"
        " --ref {AUGREF}"
        " --out-dir {OUT_DIR}"
        " --out-name {OUT_NAME}.vcf.gz"
        " " + decon_opts() +
        " --cpus {threads} --mem {config[mem_gb]}gb"
        " --local"

rule split_vcf:
    """VCF → onref / nestedref / offref"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.onref.vcf.gz",
        f"{OUT_DIR}/{OUT_NAME}.nestedref.vcf.gz",
        f"{OUT_DIR}/{OUT_NAME}.offref.vcf.gz",
    shell:
        "cd {OUT_DIR} && {workflow.basedir}/scripts/split-ref.sh"
        " -v {workflow.basedir}/{input}"
        " -p {AUGREF}"

rule haplotypes:
    """GBZ → .hapl index"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.hapl",
    threads: config.get("cpus", 8)
    resources:
        mem_mb=config.get("mem_gb", 200) * 1024,
        runtime=config.get("runtime_min", 960),
    shell:
        "scripts/haplotypes.sh"
        " --gbz {input}"
        " --ref {REF}"
        " --out-dir {OUT_DIR}"
        " --out-name {OUT_NAME}.hapl"
        " --cpus {threads} --mem {config[mem_gb]}gb"
        " --local"

rule fasta:
    """GBZ → augmented reference FASTA"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.fa.gz",
    threads: config.get("cpus", 8)
    resources:
        mem_mb=config.get("mem_gb", 200) * 1024,
        runtime=config.get("runtime_min", 960),
    shell:
        "scripts/fasta.sh"
        " --gbz {input}"
        " --ref {AUGREF}"
        " --out-dir {OUT_DIR}"
        " --out-name {OUT_NAME}.fa.gz"
        " --cpus {threads} --mem {config[mem_gb]}gb"
        " --local"

rule length_hist:
    """Augref segments → length histogram"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{OUT_NAME}.augref-length-hist.png",
    shell:
        "Rscript scripts/offref-length-hist.R {output} {input} TRUE"

rule plots:
    """Deconstruct VCF → off-reference density ideogram"""
    input:
        vcf=f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{OUT_NAME}.offref.png",
    shell:
        "Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} Off-Reference Variant Density'"
        " 0 {config[refgaps_bed]} {config[scale_type]}"
        " --ref {REF} --offref"

############################################################################
# Per-sample rules (wildcard: {sample})
############################################################################

rule giraffe:
    """GBZ + reads → GAM"""
    input:
        gbz=map_gbz,
        hapl=hapl_index,
        reads=lambda wc: config["samples"][wc.sample],
    output:
        f"{OUT_DIR}/{{sample}}.gam",
    threads: config.get("cpus", 8)
    resources:
        mem_mb=config.get("mem_gb", 200) * 1024,
        runtime=config.get("runtime_min", 960),
    shell:
        "scripts/giraffe.sh"
        " --gbz {input.gbz}"
        " --hapl {input.hapl}"
        " --reads {input.reads}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.gam"
        " --cpus {threads} --mem {config[mem_gb]}gb"
        " --local"

rule call:
    """GAM → VCF (vg call)"""
    input:
        gam=f"{OUT_DIR}/{{sample}}.gam",
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{{sample}}.vcf.gz",
    threads: config.get("cpus", 8)
    resources:
        mem_mb=config.get("mem_gb", 200) * 1024,
        runtime=config.get("runtime_min", 960),
    shell:
        "scripts/call.sh"
        " --gbz {input.gbz}"
        " --gam {input.gam}"
        " --ref {AUGREF}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.vcf.gz"
        " --cpus {threads} --mem {config[mem_gb]}gb"
        " --local"

rule surject:
    """GAM → sorted BAM"""
    input:
        gam=f"{OUT_DIR}/{{sample}}.gam",
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{{sample}}.bam",
    threads: config.get("cpus", 8)
    resources:
        mem_mb=config.get("mem_gb", 200) * 1024,
        runtime=config.get("runtime_min", 960),
    shell:
        "scripts/surject.sh"
        " --gbz {input.gbz}"
        " --gam {input.gam}"
        " --ref {AUGREF}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.bam"
        " --cpus {threads} --mem {config[mem_gb]}gb"
        " --local"

rule deepvariant:
    """BAM + FASTA → VCF via DeepVariant Docker"""
    input:
        bam=f"{OUT_DIR}/{{sample}}.bam",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
    output:
        f"{OUT_DIR}/{{sample}}.deepvariant.vcf.gz",
    threads: config.get("cpus", 8)
    resources:
        mem_mb=config.get("mem_gb", 200) * 1024,
        runtime=config.get("runtime_min", 960),
    shell:
        "scripts/deepvariant.sh"
        " --bam {input.bam}"
        " --ref {input.ref}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.deepvariant.vcf.gz"
        " --dv-version {config[dv_version]}"
        " --cpus {threads} --mem {config[mem_gb]}gb"
        " --local"

rule call_plots:
    """Per-sample call VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/{{sample}}.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{{sample}}.call-offref.png",
    shell:
        "Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} Call Off-Reference Density ({wildcards.sample})'"
        " 0 {config[refgaps_bed]} {config[scale_type]}"
        " --ref {REF} --offref"

rule dv_plots:
    """Per-sample DeepVariant VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/{{sample}}.deepvariant.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{{sample}}.dv-offref.png",
    shell:
        "Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} DeepVariant Off-Reference Density ({wildcards.sample})'"
        " 0 {config[refgaps_bed]} {config[scale_type]}"
        " --ref {REF} --offref"

############################################################################
# Batch merge rules
############################################################################

rule merge_call_vcfs:
    """Merge per-sample call VCFs with bcftools, add AF/AC/AN tags"""
    input:
        expand("{out}/{s}.vcf.gz", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/merged.call.vcf.gz",
    shell:
        "bcftools merge {input} -Oz"
        " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
        " && tabix -p vcf {output}"

rule merge_dv_vcfs:
    """Merge per-sample DeepVariant VCFs with bcftools, add AF/AC/AN tags"""
    input:
        expand("{out}/{s}.deepvariant.vcf.gz", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/merged.deepvariant.vcf.gz",
    shell:
        "bcftools merge {input} -Oz"
        " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
        " && tabix -p vcf {output}"

rule merged_call_plots:
    """Merged call VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/merged.call.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/merged.call-offref.png",
    shell:
        "Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} Merged Call Off-Reference Density'"
        " 0 {config[refgaps_bed]} {config[scale_type]}"
        " --ref {REF} --offref"

rule merged_dv_plots:
    """Merged DeepVariant VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/merged.deepvariant.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/merged.dv-offref.png",
    shell:
        "Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} Merged DeepVariant Off-Reference Density'"
        " 0 {config[refgaps_bed]} {config[scale_type]}"
        " --ref {REF} --offref"

############################################################################
# VCF statistics rules
############################################################################

rule deconstruct_stats:
    """Deconstruct VCF → variant stats + plots (includes AF spectrum)"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.af-spectrum.png",
    shell:
        "Rscript scripts/vcf-stats.R {input} {OUT_DIR}/{OUT_NAME}"
        " --title '{REF} Deconstruct'"

rule call_stats:
    """Per-sample call VCF → variant stats + plots"""
    input:
        f"{OUT_DIR}/{{sample}}.vcf.gz",
    output:
        f"{OUT_DIR}/{{sample}}.call.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.call.variant-types.png",
        f"{OUT_DIR}/{{sample}}.call.size-dist.png",
    shell:
        "Rscript scripts/vcf-stats.R {input} {OUT_DIR}/{wildcards.sample}.call"
        " --title '{REF} Call ({wildcards.sample})'"

rule dv_stats:
    """Per-sample DeepVariant VCF → variant stats + plots"""
    input:
        f"{OUT_DIR}/{{sample}}.deepvariant.vcf.gz",
    output:
        f"{OUT_DIR}/{{sample}}.dv.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.dv.variant-types.png",
        f"{OUT_DIR}/{{sample}}.dv.size-dist.png",
    shell:
        "Rscript scripts/vcf-stats.R {input} {OUT_DIR}/{wildcards.sample}.dv"
        " --title '{REF} DeepVariant ({wildcards.sample})'"

rule merged_call_stats:
    """Merged call VCF → variant stats + plots (includes AF spectrum)"""
    input:
        f"{OUT_DIR}/merged.call.vcf.gz",
    output:
        f"{OUT_DIR}/merged.call.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.variant-types.png",
        f"{OUT_DIR}/merged.call.size-dist.png",
        f"{OUT_DIR}/merged.call.af-spectrum.png",
    shell:
        "Rscript scripts/vcf-stats.R {input} {OUT_DIR}/merged.call"
        " --title '{REF} Merged Call'"

rule merged_dv_stats:
    """Merged DeepVariant VCF → variant stats + plots (includes AF spectrum)"""
    input:
        f"{OUT_DIR}/merged.deepvariant.vcf.gz",
    output:
        f"{OUT_DIR}/merged.dv.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.variant-types.png",
        f"{OUT_DIR}/merged.dv.size-dist.png",
        f"{OUT_DIR}/merged.dv.af-spectrum.png",
    shell:
        "Rscript scripts/vcf-stats.R {input} {OUT_DIR}/merged.dv"
        " --title '{REF} Merged DeepVariant'"
