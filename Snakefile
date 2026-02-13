############################################################################
# Snakefile — nested-variants pipeline
#
# Default config is for local chr20 testing.
# Override via config.local.yaml or command line:
#   snakemake --config ref=CHM13 vg='...' graph_only
############################################################################

configfile: "config.yaml"

REF      = config["ref"]
AUGREF   = f"augref_{REF}"
OUT_DIR  = config["out_dir"]
OUT_NAME = config["out_name"]

# Expand VG input glob (supports bash extglob patterns like !(*.d*).vg)
import subprocess as _sp
_vg_result = _sp.run(
    ["bash", "-O", "extglob", "-c", f"echo {config['vg']}"],
    capture_output=True, text=True
)
VG_FILES = sorted(_vg_result.stdout.split())
if not VG_FILES or VG_FILES == [config["vg"]]:
    # No expansion happened — treat as literal path
    VG_FILES = [config["vg"]]

# Load samples from config dict
SAMPLES = list(config.get("samples", {}).keys())

# Load samples from TSV if provided
if config.get("samples_tsv"):
    import csv
    with open(config["samples_tsv"]) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            config.setdefault("samples", {})[row["sample"]] = row["reads_index"]
    SAMPLES = list(config["samples"].keys())

# Constrain {sample} wildcard to configured sample names only, preventing
# ambiguity between deconstruct ({OUT_NAME}.vcf.gz) and call ({sample}.vcf.gz)
wildcard_constraints:
    sample="|".join(SAMPLES) if SAMPLES else "$^"

# Per-rule resource helpers: look up rule-specific config, fall back to global default
def rule_cpus(rule_name, default):
    return config.get(f"{rule_name}_cpus", config.get("cpus", default))

def rule_mem_gb(rule_name, default):
    return config.get(f"{rule_name}_mem_gb", config.get("mem_gb", default))

def rule_runtime(rule_name, default=960):
    return config.get(f"{rule_name}_runtime_min", config.get("runtime_min", default))

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

# Annotation helpers
def annotation_inputs():
    """Return list of configured annotation BED files."""
    return [config[k] for k in ["annot_genes", "annot_repeats", "annot_segdups", "annot_censat"] if config.get(k, "")]

def annotation_names():
    """Return clean display names for configured annotations."""
    names = []
    for k, name in [("annot_genes", "genes"), ("annot_repeats", "repeats"), ("annot_segdups", "segdups"), ("annot_censat", "censat")]:
        if config.get(k, ""):
            names.append(name)
    return names

def annotation_outputs():
    """Return annotation output files if any annotations are configured."""
    if annotation_inputs():
        return [
            f"{OUT_DIR}/{OUT_NAME}.annot-summary.png",
            f"{OUT_DIR}/{OUT_NAME}.annot-scatter.png",
            f"{OUT_DIR}/{OUT_NAME}.annot-stats.tsv",
        ]
    return []

############################################################################
# Target rules
############################################################################

rule all:
    """Full pipeline: graph + genotype + deepvariant + merge + all plots/stats"""
    input:
        # graph_only outputs
        f"{OUT_DIR}/{OUT_NAME}.offref.png",
        f"{OUT_DIR}/{OUT_NAME}.augref-length-hist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.sites.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.af-spectrum.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variants.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.af-spectrum.png",
        *annotation_outputs(),
        # per-sample genotyping outputs
        expand("{out}/{s}.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.sites.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.sites.variant-types.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.sites.size-dist.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.variants.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.variants.variant-types.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.variants.size-dist.png", out=OUT_DIR, s=SAMPLES),
        # per-sample deepvariant outputs
        expand("{out}/{s}.deepvariant.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.sites.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.sites.variant-types.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.sites.size-dist.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.variants.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.variants.variant-types.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.variants.size-dist.png", out=OUT_DIR, s=SAMPLES),
        # merged outputs
        f"{OUT_DIR}/merged.call.vcf.gz",
        f"{OUT_DIR}/merged.deepvariant.vcf.gz",
        f"{OUT_DIR}/merged.call-offref.png",
        f"{OUT_DIR}/merged.dv-offref.png",
        f"{OUT_DIR}/merged.call.sites.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.sites.variant-types.png",
        f"{OUT_DIR}/merged.call.sites.size-dist.png",
        f"{OUT_DIR}/merged.call.sites.af-spectrum.png",
        f"{OUT_DIR}/merged.call.variants.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.variants.variant-types.png",
        f"{OUT_DIR}/merged.call.variants.size-dist.png",
        f"{OUT_DIR}/merged.call.variants.af-spectrum.png",
        f"{OUT_DIR}/merged.dv.sites.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.sites.variant-types.png",
        f"{OUT_DIR}/merged.dv.sites.size-dist.png",
        f"{OUT_DIR}/merged.dv.sites.af-spectrum.png",
        f"{OUT_DIR}/merged.dv.variants.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.variants.variant-types.png",
        f"{OUT_DIR}/merged.dv.variants.size-dist.png",
        f"{OUT_DIR}/merged.dv.variants.af-spectrum.png",

rule graph_only:
    """Graph construction + deconstruct + plots (no genotyping)"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.offref.png",
        f"{OUT_DIR}/{OUT_NAME}.augref-length-hist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.sites.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.af-spectrum.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variants.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.af-spectrum.png",
        *annotation_outputs(),

rule genotype_all:
    """Genotype all samples (vg call)"""
    input:
        expand("{out}/{s}.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.sites.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.sites.variant-types.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.sites.size-dist.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.variants.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.variants.variant-types.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.variants.size-dist.png", out=OUT_DIR, s=SAMPLES),

rule deepvariant_all:
    """Run DeepVariant on all samples"""
    input:
        expand("{out}/{s}.deepvariant.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.sites.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.sites.variant-types.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.sites.size-dist.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.variants.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.variants.variant-types.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.variants.size-dist.png", out=OUT_DIR, s=SAMPLES),

# Resolve wildcard ambiguities:
# - {OUT_NAME}.vcf.gz matches both deconstruct and call (sample={OUT_NAME})
# - {sample}.deepvariant.vcf.gz matches both deepvariant and call (sample=X.deepvariant)
ruleorder: deconstruct > call
ruleorder: deepvariant > call

############################################################################
# Graph construction rules (run once)
############################################################################

rule paths:
    """VG → GBZ via augmented reference paths"""
    input:
        ancient(VG_FILES),
    output:
        f"{OUT_DIR}/{OUT_NAME}.gbz",
        f"{OUT_DIR}/{OUT_NAME}.gfa.gz",
        f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    threads: rule_cpus("paths", 128)
    resources:
        mem_mb=rule_mem_gb("paths", 512) * 1024,
        runtime=rule_runtime("paths"),
    params:
        mem_gb=rule_mem_gb("paths", 512),
    shell:
        "scripts/paths.sh"
        " --vg '{input}'"
        " --ref {REF}"
        " --out-dir {OUT_DIR}"
        " --out-name {OUT_NAME}.gfa.gz"
        " --min-augref-len {config[min_augref_len]}"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule deconstruct:
    """GBZ → VCF"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
    threads: rule_cpus("deconstruct", 128)
    resources:
        mem_mb=rule_mem_gb("deconstruct", 512) * 1024,
        runtime=rule_runtime("deconstruct"),
    params:
        mem_gb=rule_mem_gb("deconstruct", 512),
    shell:
        "scripts/deconstruct.sh"
        " --gbz {input}"
        " --ref {AUGREF}"
        " --out-dir {OUT_DIR}"
        " --out-name {OUT_NAME}.vcf.gz"
        " " + decon_opts() +
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule split_vcf:
    """VCF → onref / nestedref / offref"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.onref.vcf.gz",
        f"{OUT_DIR}/{OUT_NAME}.nestedref.vcf.gz",
        f"{OUT_DIR}/{OUT_NAME}.offref.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
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
    threads: rule_cpus("haplotypes", 128)
    resources:
        mem_mb=rule_mem_gb("haplotypes", 512) * 1024,
        runtime=rule_runtime("haplotypes"),
    params:
        mem_gb=rule_mem_gb("haplotypes", 512),
    shell:
        "scripts/haplotypes.sh"
        " --gbz {input}"
        " --ref {REF}"
        " --out-dir {OUT_DIR}"
        " --out-name {OUT_NAME}.hapl"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule fasta:
    """GBZ → augmented reference FASTA"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.fa.gz",
    threads: rule_cpus("fasta", 128)
    resources:
        mem_mb=rule_mem_gb("fasta", 512) * 1024,
        runtime=rule_runtime("fasta"),
    params:
        mem_gb=rule_mem_gb("fasta", 512),
    shell:
        "scripts/fasta.sh"
        " --gbz {input}"
        " --ref {AUGREF}"
        " --out-dir {OUT_DIR}"
        " --out-name {OUT_NAME}.fa.gz"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule length_hist:
    """Augref segments → length histogram"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{OUT_NAME}.augref-length-hist.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/offref-length-hist.R {output} {input} TRUE"

rule plots:
    """Deconstruct VCF → off-reference density ideogram"""
    input:
        vcf=f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{OUT_NAME}.offref.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} Off-Reference Variant Density'"
        " 0 {config[refgaps_bed]} {config[scale_type]}"
        " --ref {REF} --offref"

############################################################################
# Annotation overlap rules (optional — only when annot_* keys are set)
############################################################################

rule annotation_intersect:
    """Augref segments + annotation BEDs → per-segment annotation TSV"""
    input:
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annots=annotation_inputs(),
    output:
        f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv",
    resources:
        mem_mb=512000,
        runtime=2880,
    params:
        names=" ".join(annotation_names()),
        group_arg="--group-by-column 6" if config.get("annot_repeats", "") else "",
    shell:
        "python scripts/intersect-annotations.py"
        " {input.segs} {input.annots}"
        " --per-segment --per-segment-output {output}"
        " --output-dir {OUT_DIR}"
        " --annotation-names {params.names}"
        " {params.group_arg}"

rule annotation_plots:
    """Per-segment annotation TSV → overlap plots + stats"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv",
    output:
        f"{OUT_DIR}/{OUT_NAME}.annot-summary.png",
        f"{OUT_DIR}/{OUT_NAME}.annot-scatter.png",
        f"{OUT_DIR}/{OUT_NAME}.annot-stats.tsv",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/annotation-plots.R {input} {OUT_DIR}/{OUT_NAME}"
        " --title '{REF} Annotation Overlap'"

############################################################################
# Per-sample rules (wildcard: {sample})
############################################################################

rule giraffe:
    """GBZ + reads → GAM"""
    input:
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
        hapl=f"{OUT_DIR}/{OUT_NAME}.hapl",
        reads=lambda wc: config["samples"][wc.sample],
    output:
        f"{OUT_DIR}/{{sample}}.gam",
    threads: rule_cpus("giraffe", 128)
    resources:
        mem_mb=rule_mem_gb("giraffe", 512) * 1024,
        runtime=rule_runtime("giraffe"),
    params:
        mem_gb=rule_mem_gb("giraffe", 512),
    shell:
        "scripts/giraffe.sh"
        " --gbz {input.gbz}"
        " --hapl {input.hapl}"
        " --reads {input.reads}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.gam"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule call:
    """GAM → VCF (vg call)"""
    input:
        gam=f"{OUT_DIR}/{{sample}}.gam",
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{{sample}}.vcf.gz",
    threads: rule_cpus("call", 128)
    resources:
        mem_mb=rule_mem_gb("call", 512) * 1024,
        runtime=rule_runtime("call"),
    params:
        mem_gb=rule_mem_gb("call", 512),
    shell:
        "scripts/call.sh"
        " --gbz {input.gbz}"
        " --gam {input.gam}"
        " --ref {AUGREF}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.vcf.gz"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule surject:
    """GAM → sorted BAM"""
    input:
        gam=f"{OUT_DIR}/{{sample}}.gam",
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{{sample}}.bam",
    threads: rule_cpus("surject", 128)
    resources:
        mem_mb=rule_mem_gb("surject", 512) * 1024,
        runtime=rule_runtime("surject"),
    params:
        mem_gb=rule_mem_gb("surject", 512),
    shell:
        "scripts/surject.sh"
        " --gbz {input.gbz}"
        " --gam {input.gam}"
        " --ref {AUGREF}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.bam"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule deepvariant:
    """BAM + FASTA → VCF via DeepVariant Docker"""
    input:
        bam=f"{OUT_DIR}/{{sample}}.bam",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
    output:
        f"{OUT_DIR}/{{sample}}.deepvariant.vcf.gz",
    threads: rule_cpus("deepvariant", 96)
    resources:
        mem_mb=rule_mem_gb("deepvariant", 1024) * 1024,
        runtime=rule_runtime("deepvariant"),
    params:
        mem_gb=rule_mem_gb("deepvariant", 1024),
    shell:
        "scripts/deepvariant.sh"
        " --bam {input.bam}"
        " --ref {input.ref}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.deepvariant.vcf.gz"
        " --dv-version {config[dv_version]}"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule call_plots:
    """Per-sample call VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/{{sample}}.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{{sample}}.call-offref.png",
    resources:
        mem_mb=256000,
        runtime=2880,
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
    resources:
        mem_mb=256000,
        runtime=2880,
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
    resources:
        mem_mb=256000,
        runtime=2880,
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
    resources:
        mem_mb=256000,
        runtime=2880,
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
    resources:
        mem_mb=256000,
        runtime=2880,
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
    resources:
        mem_mb=256000,
        runtime=2880,
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
        f"{OUT_DIR}/{OUT_NAME}.sites.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.sites.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.af-spectrum.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variants.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.af-spectrum.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/vcf-stats.R {input} {OUT_DIR}/{OUT_NAME}.sites"
        " --mode sites --title '{REF} Deconstruct'"
        " && Rscript scripts/vcf-stats.R {input} {OUT_DIR}/{OUT_NAME}.variants"
        " --mode variants --title '{REF} Deconstruct'"

rule call_stats:
    """Per-sample call VCF → variant stats + plots"""
    input:
        f"{OUT_DIR}/{{sample}}.vcf.gz",
    output:
        f"{OUT_DIR}/{{sample}}.call.sites.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.call.sites.variant-types.png",
        f"{OUT_DIR}/{{sample}}.call.sites.size-dist.png",
        f"{OUT_DIR}/{{sample}}.call.variants.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.call.variants.variant-types.png",
        f"{OUT_DIR}/{{sample}}.call.variants.size-dist.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/vcf-stats.R {input} {OUT_DIR}/{wildcards.sample}.call.sites"
        " --mode sites --title '{REF} Call ({wildcards.sample})'"
        " && Rscript scripts/vcf-stats.R {input} {OUT_DIR}/{wildcards.sample}.call.variants"
        " --mode variants --title '{REF} Call ({wildcards.sample})'"

rule dv_stats:
    """Per-sample DeepVariant VCF → variant stats + plots"""
    input:
        f"{OUT_DIR}/{{sample}}.deepvariant.vcf.gz",
    output:
        f"{OUT_DIR}/{{sample}}.dv.sites.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.dv.sites.variant-types.png",
        f"{OUT_DIR}/{{sample}}.dv.sites.size-dist.png",
        f"{OUT_DIR}/{{sample}}.dv.variants.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.dv.variants.variant-types.png",
        f"{OUT_DIR}/{{sample}}.dv.variants.size-dist.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/vcf-stats.R {input} {OUT_DIR}/{wildcards.sample}.dv.sites"
        " --mode sites --title '{REF} DeepVariant ({wildcards.sample})'"
        " && Rscript scripts/vcf-stats.R {input} {OUT_DIR}/{wildcards.sample}.dv.variants"
        " --mode variants --title '{REF} DeepVariant ({wildcards.sample})'"

rule merged_call_stats:
    """Merged call VCF → variant stats + plots (includes AF spectrum)"""
    input:
        f"{OUT_DIR}/merged.call.vcf.gz",
    output:
        f"{OUT_DIR}/merged.call.sites.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.sites.variant-types.png",
        f"{OUT_DIR}/merged.call.sites.size-dist.png",
        f"{OUT_DIR}/merged.call.sites.af-spectrum.png",
        f"{OUT_DIR}/merged.call.variants.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.variants.variant-types.png",
        f"{OUT_DIR}/merged.call.variants.size-dist.png",
        f"{OUT_DIR}/merged.call.variants.af-spectrum.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/vcf-stats.R {input} {OUT_DIR}/merged.call.sites"
        " --mode sites --title '{REF} Merged Call'"
        " && Rscript scripts/vcf-stats.R {input} {OUT_DIR}/merged.call.variants"
        " --mode variants --title '{REF} Merged Call'"

rule merged_dv_stats:
    """Merged DeepVariant VCF → variant stats + plots (includes AF spectrum)"""
    input:
        f"{OUT_DIR}/merged.deepvariant.vcf.gz",
    output:
        f"{OUT_DIR}/merged.dv.sites.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.sites.variant-types.png",
        f"{OUT_DIR}/merged.dv.sites.size-dist.png",
        f"{OUT_DIR}/merged.dv.sites.af-spectrum.png",
        f"{OUT_DIR}/merged.dv.variants.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.variants.variant-types.png",
        f"{OUT_DIR}/merged.dv.variants.size-dist.png",
        f"{OUT_DIR}/merged.dv.variants.af-spectrum.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/vcf-stats.R {input} {OUT_DIR}/merged.dv.sites"
        " --mode sites --title '{REF} Merged DeepVariant'"
        " && Rscript scripts/vcf-stats.R {input} {OUT_DIR}/merged.dv.variants"
        " --mode variants --title '{REF} Merged DeepVariant'"
