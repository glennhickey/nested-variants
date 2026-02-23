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
    sample="|".join(SAMPLES) if SAMPLES else "$^",
    filt="all|pass",
    mode="sites|variants"

# Per-rule resource helpers: look up rule-specific config, fall back to global default
def rule_cpus(rule_name, default):
    return config.get(f"{rule_name}_cpus", config.get("cpus", default))

def rule_mem_gb(rule_name, default):
    return config.get(f"{rule_name}_mem_gb", config.get("mem_gb", default))

def rule_runtime(rule_name, default=2880):
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
    return [config[k] for k in ["annot_genes", "annot_repeats", "annot_segdups", "annot_censat", "annot_pclai", "annot_giab"] if config.get(k, "")]

def annotation_names():
    """Return clean display names for configured annotations."""
    names = []
    for k, name in [("annot_genes", "genes"), ("annot_repeats", "repeats"), ("annot_segdups", "segdups"), ("annot_censat", "censat"), ("annot_pclai", "pclai"), ("annot_giab", "giab")]:
        if config.get(k, ""):
            names.append(name)
    return names

def annotation_outputs():
    """Return annotation output files if any annotations are configured."""
    if annotation_inputs():
        outputs = [
            f"{OUT_DIR}/{OUT_NAME}.annot-summary.png",
            f"{OUT_DIR}/{OUT_NAME}.annot-scatter.png",
            f"{OUT_DIR}/{OUT_NAME}.annot-cooccur.png",
            f"{OUT_DIR}/{OUT_NAME}.annot-stats.tsv",
        ]
        if config.get("annot_pclai", ""):
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.annot-ancestry.png")
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.annot-pclai-summary.png")
        return outputs
    return []

def annotation_snp_outputs(callers=None):
    """Return annotation SNP heatmap outputs if annotations are configured.

    callers: list of caller types to include, e.g. ["call"], ["deepvariant"],
             ["deconstruct"], or None for all.
    Deconstruct VCFs have no FILTER field, so only "all" is generated for them.
    """
    if not annotation_inputs():
        return []
    if callers is None:
        callers = ["deconstruct", "call", "deepvariant"]
    outputs = []
    for plot in ["annot-snp-counts", "annot-snp-tstv"]:
        if "deconstruct" in callers:
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.{plot}.all.png")
        if "call" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.{plot}.{filt}.png")
                outputs.append(f"{OUT_DIR}/merged.call.{plot}.{filt}.png")
        if "deepvariant" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.deepvariant.{plot}.{filt}.png")
                outputs.append(f"{OUT_DIR}/merged.deepvariant.{plot}.{filt}.png")
    if config.get("annot_pclai", ""):
        for plot in ["annot-pclai-snp-counts", "annot-pclai-snp-tstv"]:
            if "deconstruct" in callers:
                outputs.append(f"{OUT_DIR}/{OUT_NAME}.{plot}.all.png")
            if "call" in callers:
                for filt in ["all", "pass"]:
                    for s in SAMPLES:
                        outputs.append(f"{OUT_DIR}/{s}.{plot}.{filt}.png")
                    outputs.append(f"{OUT_DIR}/merged.call.{plot}.{filt}.png")
            if "deepvariant" in callers:
                for filt in ["all", "pass"]:
                    for s in SAMPLES:
                        outputs.append(f"{OUT_DIR}/{s}.deepvariant.{plot}.{filt}.png")
                    outputs.append(f"{OUT_DIR}/merged.deepvariant.{plot}.{filt}.png")
    return outputs

def annotation_stats_outputs(callers=None):
    """Return annotation-stratified VCF stats outputs when annotations are configured.

    callers: list of caller types to include, e.g. ["deconstruct"], ["call"],
             ["deepvariant"], or None for all.
    """
    if not annotation_inputs():
        return []
    if callers is None:
        callers = ["deconstruct", "call", "deepvariant"]
    outputs = []
    for suffix in ["variant-types-by-annot.png", "vcf-stats-by-annot.tsv"]:
        if "deconstruct" in callers:
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.sites.{suffix}")
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.variants.{suffix}")
        if "call" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.call.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/{s}.call.variants.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.call.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.call.variants.{filt}.{suffix}")
        if "deepvariant" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.dv.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/{s}.dv.variants.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.dv.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.dv.variants.{filt}.{suffix}")
    return outputs

def augref_annot_beds():
    """Return augref-space annotation BED files for configured annotations."""
    if not annotation_inputs():
        return []
    return [f"{OUT_DIR}/{OUT_NAME}.augref-annot-{n}.bed" for n in annotation_names()]

def polymorphism_outputs():
    """Return segment polymorphism table output."""
    return [f"{OUT_DIR}/{OUT_NAME}.segment-polymorphism.tsv"]

############################################################################
# Target rules
############################################################################

rule all:
    """Full pipeline: graph + genotype + deepvariant + merge + all plots/stats"""
    input:
        # graph_only outputs
        f"{OUT_DIR}/{OUT_NAME}.offref.png",
        f"{OUT_DIR}/{OUT_NAME}.offref-segs.png",
        f"{OUT_DIR}/{OUT_NAME}.augref-length-hist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.sites.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist-log.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.af-spectrum.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variants.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist-log.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.af-spectrum.png",
        *annotation_outputs(),
        *annotation_snp_outputs(),
        *annotation_stats_outputs(),
        *polymorphism_outputs(),
        # per-sample genotyping outputs
        expand("{out}/{s}.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.sites.{filt}.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.sites.{filt}.variant-types.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.sites.{filt}.size-dist.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.sites.{filt}.size-dist-log.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.variants.{filt}.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.variants.{filt}.variant-types.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.variants.{filt}.size-dist.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.variants.{filt}.size-dist-log.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        # per-sample deepvariant outputs
        expand("{out}/{s}.deepvariant.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.sites.{filt}.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.sites.{filt}.variant-types.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.sites.{filt}.size-dist.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.sites.{filt}.size-dist-log.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.variants.{filt}.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.variants.{filt}.variant-types.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.variants.{filt}.size-dist.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.variants.{filt}.size-dist-log.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        # merged outputs
        f"{OUT_DIR}/merged.call.vcf.gz",
        f"{OUT_DIR}/merged.deepvariant.vcf.gz",
        f"{OUT_DIR}/merged.call-offref.png",
        f"{OUT_DIR}/merged.dv-offref.png",
        f"{OUT_DIR}/merged.call.sites.all.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.sites.all.variant-types.png",
        f"{OUT_DIR}/merged.call.sites.all.size-dist.png",
        f"{OUT_DIR}/merged.call.sites.all.size-dist-log.png",
        f"{OUT_DIR}/merged.call.sites.all.af-spectrum.png",
        f"{OUT_DIR}/merged.call.sites.pass.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.sites.pass.variant-types.png",
        f"{OUT_DIR}/merged.call.sites.pass.size-dist.png",
        f"{OUT_DIR}/merged.call.sites.pass.size-dist-log.png",
        f"{OUT_DIR}/merged.call.sites.pass.af-spectrum.png",
        f"{OUT_DIR}/merged.call.variants.all.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.variants.all.variant-types.png",
        f"{OUT_DIR}/merged.call.variants.all.size-dist.png",
        f"{OUT_DIR}/merged.call.variants.all.size-dist-log.png",
        f"{OUT_DIR}/merged.call.variants.all.af-spectrum.png",
        f"{OUT_DIR}/merged.call.variants.pass.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.variants.pass.variant-types.png",
        f"{OUT_DIR}/merged.call.variants.pass.size-dist.png",
        f"{OUT_DIR}/merged.call.variants.pass.size-dist-log.png",
        f"{OUT_DIR}/merged.call.variants.pass.af-spectrum.png",
        f"{OUT_DIR}/merged.dv.sites.all.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.sites.all.variant-types.png",
        f"{OUT_DIR}/merged.dv.sites.all.size-dist.png",
        f"{OUT_DIR}/merged.dv.sites.all.size-dist-log.png",
        f"{OUT_DIR}/merged.dv.sites.all.af-spectrum.png",
        f"{OUT_DIR}/merged.dv.sites.pass.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.sites.pass.variant-types.png",
        f"{OUT_DIR}/merged.dv.sites.pass.size-dist.png",
        f"{OUT_DIR}/merged.dv.sites.pass.size-dist-log.png",
        f"{OUT_DIR}/merged.dv.sites.pass.af-spectrum.png",
        f"{OUT_DIR}/merged.dv.variants.all.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.variants.all.variant-types.png",
        f"{OUT_DIR}/merged.dv.variants.all.size-dist.png",
        f"{OUT_DIR}/merged.dv.variants.all.size-dist-log.png",
        f"{OUT_DIR}/merged.dv.variants.all.af-spectrum.png",
        f"{OUT_DIR}/merged.dv.variants.pass.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.variants.pass.variant-types.png",
        f"{OUT_DIR}/merged.dv.variants.pass.size-dist.png",
        f"{OUT_DIR}/merged.dv.variants.pass.size-dist-log.png",
        f"{OUT_DIR}/merged.dv.variants.pass.af-spectrum.png",

rule graph_only:
    """Graph construction + deconstruct + plots (no genotyping)"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.offref.png",
        f"{OUT_DIR}/{OUT_NAME}.offref-segs.png",
        f"{OUT_DIR}/{OUT_NAME}.augref-length-hist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.sites.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist-log.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.af-spectrum.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variants.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist-log.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.af-spectrum.png",
        *annotation_outputs(),
        *annotation_snp_outputs(["deconstruct"]),
        *annotation_stats_outputs(["deconstruct"]),
        *polymorphism_outputs(),

rule genotype_all:
    """Genotype all samples (vg call) + merge"""
    input:
        expand("{out}/{s}.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call.sites.{filt}.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.sites.{filt}.variant-types.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.sites.{filt}.size-dist.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.sites.{filt}.size-dist-log.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.variants.{filt}.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.variants.{filt}.variant-types.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.variants.{filt}.size-dist.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.call.variants.{filt}.size-dist-log.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        # merged call outputs
        f"{OUT_DIR}/merged.call.vcf.gz",
        f"{OUT_DIR}/merged.call-offref.png",
        f"{OUT_DIR}/merged.call.sites.all.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.sites.all.variant-types.png",
        f"{OUT_DIR}/merged.call.sites.all.size-dist.png",
        f"{OUT_DIR}/merged.call.sites.all.size-dist-log.png",
        f"{OUT_DIR}/merged.call.sites.all.af-spectrum.png",
        f"{OUT_DIR}/merged.call.sites.pass.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.sites.pass.variant-types.png",
        f"{OUT_DIR}/merged.call.sites.pass.size-dist.png",
        f"{OUT_DIR}/merged.call.sites.pass.size-dist-log.png",
        f"{OUT_DIR}/merged.call.sites.pass.af-spectrum.png",
        f"{OUT_DIR}/merged.call.variants.all.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.variants.all.variant-types.png",
        f"{OUT_DIR}/merged.call.variants.all.size-dist.png",
        f"{OUT_DIR}/merged.call.variants.all.size-dist-log.png",
        f"{OUT_DIR}/merged.call.variants.all.af-spectrum.png",
        f"{OUT_DIR}/merged.call.variants.pass.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.variants.pass.variant-types.png",
        f"{OUT_DIR}/merged.call.variants.pass.size-dist.png",
        f"{OUT_DIR}/merged.call.variants.pass.size-dist-log.png",
        f"{OUT_DIR}/merged.call.variants.pass.af-spectrum.png",
        *annotation_snp_outputs(["call"]),
        *annotation_stats_outputs(["call"]),

rule deepvariant_all:
    """Run DeepVariant on all samples"""
    input:
        expand("{out}/{s}.deepvariant.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.dv.sites.{filt}.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.sites.{filt}.variant-types.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.sites.{filt}.size-dist.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.sites.{filt}.size-dist-log.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.variants.{filt}.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.variants.{filt}.variant-types.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.variants.{filt}.size-dist.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.dv.variants.{filt}.size-dist-log.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        *annotation_snp_outputs(["deepvariant"]),
        *annotation_stats_outputs(["deepvariant"]),

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
        " --out-name {OUT_NAME}"
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

rule segment_density:
    """Augref segments → off-reference segment density ideogram"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{OUT_NAME}.offref-segs.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/chrom-density-tsv.R"
        " {input} {output}"
        " '{REF} Off-Reference Segment Density'"
        " {config[min_augref_len]} '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF}"

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
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
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
        group_arg="--group-by-column 6" if config.get("annot_repeats", "") or config.get("annot_pclai", "") else "",
    shell:
        "python scripts/intersect-annotations.py"
        " {input.segs} {input.annots}"
        " --per-segment --per-segment-output {output}"
        " --output-dir {OUT_DIR}"
        " --annotation-names {params.names}"
        " {params.group_arg}"

rule annotation_augref_beds:
    """Per-segment annotation TSV + annotation BEDs → augref-space BEDs"""
    input:
        seg_annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv",
        annots=annotation_inputs(),
    output:
        augref_annot_beds(),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        beds=lambda wc, input: ",".join(input.annots),
        names=",".join(annotation_names()),
    shell:
        "python scripts/make-augref-annotation-beds.py"
        " --per-segment-tsv {input.seg_annot}"
        " --annotation-beds {params.beds}"
        " --annotation-names {params.names}"
        " --ref {REF} --min-overlap 0.5"
        " --output-dir {OUT_DIR} --output-prefix {OUT_NAME}"

rule annotation_plots:
    """Per-segment annotation TSV → overlap plots + stats"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv",
    output:
        annotation_outputs(),
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/annotation-plots.R {input} {OUT_DIR}/{OUT_NAME}"
        " --min-overlap 0.5 --title '{REF} Annotation Overlap'"

rule norm_vcf:
    """VCF → multi-allelic split VCF via bcftools norm (shared intermediate)"""
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.normed.vcf.gz",
    threads: 4
    resources:
        mem_mb=128000,
        runtime=2880,
    shell:
        "bcftools norm -m- --threads {threads} '{input}' -Oz -o {output} 2>/dev/null"
        " && tabix -p vcf {output}"

rule biallelic_snps:
    """Normed VCF → biallelic SNP VCF (filter to true SNPs)"""
    input:
        "{prefix}.normed.vcf.gz",
    output:
        "{prefix}.biallelic-snps.{filt}.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        filt_cmd=lambda wc: "bcftools view -f PASS 2>/dev/null |" if wc.filt == "pass" else "",
    shell:
        "bcftools view '{input}' 2>/dev/null"
        " | {params.filt_cmd} bcftools view -v snps -c1"
        " -i 'STRLEN(REF)==1 && STRLEN(ALT)==1' -Oz -o {output} 2>/dev/null"
        " && tabix -p vcf {output}"

rule annotation_snp_heatmaps:
    """Annotation TSV + biallelic SNP VCF → SNP annotation heatmaps"""
    input:
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv",
        vcf=f"{OUT_DIR}/{{vcf_prefix}}.biallelic-snps.{{filt}}.vcf.gz",
    output:
        f"{OUT_DIR}/{{vcf_prefix}}.annot-snp-counts.{{filt}}.png",
        f"{OUT_DIR}/{{vcf_prefix}}.annot-snp-tstv.{{filt}}.png",
        *([f"{OUT_DIR}/{{vcf_prefix}}.annot-pclai-snp-counts.{{filt}}.png",
           f"{OUT_DIR}/{{vcf_prefix}}.annot-pclai-snp-tstv.{{filt}}.png"]
          if config.get("annot_pclai", "") else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "Rscript scripts/annotation-plots.R {input.annot}"
        " {OUT_DIR}/{wildcards.vcf_prefix}.annot-snp"
        " --vcf {input.vcf}"
        " --augref-prefix 'augref_{REF}#0#'"
        " --filter {wildcards.filt}"
        " --min-overlap 0.5"
        " --title '{REF} SNP Annotation'"

rule segment_polymorphism:
    """Deconstruct VCF → per-segment polymorphism table"""
    input:
        vcf=f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv" if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/{OUT_NAME}.segment-polymorphism.tsv",
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: f"--annot {input.annot}" if annotation_inputs() else "",
    shell:
        "Rscript scripts/segment-polymorphism.R"
        " --vcf {input.vcf}"
        " --augref-prefix '{AUGREF}#0#'"
        " {params.annot_arg}"
        " --output {output}"

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
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
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
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
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
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
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
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF} --offref"

############################################################################
# VCF statistics rules
############################################################################

rule deconstruct_sites_stats:
    """Deconstruct VCF → site-level stats + plots (includes AF spectrum)"""
    input:
        vcf=f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
        annot_beds=augref_annot_beds(),
    output:
        f"{OUT_DIR}/{OUT_NAME}.sites.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.sites.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist-log.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.af-spectrum.png",
        *([ f"{OUT_DIR}/{OUT_NAME}.sites.variant-types-by-annot.png",
            f"{OUT_DIR}/{OUT_NAME}.sites.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
    shell:
        "Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{OUT_NAME}.sites"
        " --mode sites --af-step 0.05 --title '{REF} Deconstruct'"
        " {params.annot_arg}"

rule deconstruct_variants_stats:
    """Deconstruct VCF → variant-level stats + plots (uses pre-normed VCF)"""
    input:
        vcf=f"{OUT_DIR}/{OUT_NAME}.normed.vcf.gz",
        annot_beds=augref_annot_beds(),
    output:
        f"{OUT_DIR}/{OUT_NAME}.variants.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variants.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist-log.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.af-spectrum.png",
        *([ f"{OUT_DIR}/{OUT_NAME}.variants.variant-types-by-annot.png",
            f"{OUT_DIR}/{OUT_NAME}.variants.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
    shell:
        "Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{OUT_NAME}.variants"
        " --mode variants --af-step 0.05 --title '{REF} Deconstruct'"
        " {params.annot_arg}"

rule call_stats:
    """Per-sample call VCF → variant stats + plots (one mode/filter combo)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/{wc.sample}.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/{wc.sample}.vcf.gz",
        annot_beds=augref_annot_beds(),
    output:
        f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.size-dist-log.png",
        *([ f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
    shell:
        "Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{wildcards.sample}.call.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Call ({wildcards.sample})'"
        " {params.annot_arg}"

rule dv_stats:
    """Per-sample DeepVariant VCF → variant stats + plots (one mode/filter combo)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/{wc.sample}.deepvariant.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/{wc.sample}.deepvariant.vcf.gz",
        annot_beds=augref_annot_beds(),
    output:
        f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.size-dist-log.png",
        *([ f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
    shell:
        "Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{wildcards.sample}.dv.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} DeepVariant ({wildcards.sample})'"
        " {params.annot_arg}"

rule merged_call_stats:
    """Merged call VCF → variant stats + plots (one mode/filter combo, includes AF spectrum)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.call.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.call.vcf.gz",
        annot_beds=augref_annot_beds(),
    output:
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
    shell:
        "Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.call.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged Call'"
        " {params.annot_arg}"

rule merged_dv_stats:
    """Merged DeepVariant VCF → variant stats + plots (one mode/filter combo, includes AF spectrum)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.deepvariant.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.deepvariant.vcf.gz",
        annot_beds=augref_annot_beds(),
    output:
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
    shell:
        "Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.dv.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged DeepVariant'"
        " {params.annot_arg}"
