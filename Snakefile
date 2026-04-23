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

# Long-read samples (same pipeline, separate outputs, uses -b hifi for giraffe)
LR_SAMPLES = list(config.get("longread_samples", {}).keys())
ALL_SAMPLES = SAMPLES + LR_SAMPLES
_overlap = set(SAMPLES) & set(LR_SAMPLES)
if _overlap:
    raise ValueError(f"Sample names appear in both samples and longread_samples: {_overlap}")

# Constrain {sample} wildcard to configured sample names only, preventing
# ambiguity between deconstruct ({OUT_NAME}.vcf.gz) and call ({sample}.vcf.gz)
wildcard_constraints:
    sample="|".join(ALL_SAMPLES) if ALL_SAMPLES else "$^",
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

def pangenie_enabled():
    """True when PanGenie is enabled (default true, set enable_pangenie=false to disable)."""
    val = config.get("enable_pangenie", "true")
    return str(val).lower() not in ("false", "0", "no", "")

def freebayes_longread_enabled():
    """True when FreeBayes on long reads is enabled (default false).
    FreeBayes is slow and unreliable on HiFi in complex regions (centromeres),
    so long-read samples are skipped by default. Set enable_freebayes_longread=true
    to opt in."""
    val = config.get("enable_freebayes_longread", "false")
    return str(val).lower() not in ("false", "0", "no", "")

def bcftools_enabled():
    """True when bcftools variant caller is enabled (default true)."""
    val = config.get("enable_bcftools", "true")
    return str(val).lower() not in ("false", "0", "no", "")

def bcftools_longread_enabled():
    """True when bcftools on long reads is enabled (default true)."""
    val = config.get("enable_bcftools_longread", "true")
    return str(val).lower() not in ("false", "0", "no", "")

def surject_filtering():
    """True when min_surject_len > 0 (filter contigs for surject/call)."""
    val = config.get("min_surject_len", 0)
    assert str(val).isdigit(), f"min_surject_len must be a non-negative integer, got: {val}"
    return int(val) > 0

# Annotation helpers
def annotation_inputs():
    """Return list of configured annotation BED files (excluding pclai)."""
    return [config[k] for k in ["annot_genes", "annot_repeats", "annot_segdups", "annot_censat"] if config.get(k, "")]

def annotation_names():
    """Return clean display names for configured annotations (excluding pclai)."""
    names = []
    for k, name in [("annot_genes", "genes"), ("annot_repeats", "repeats"), ("annot_segdups", "segdups"), ("annot_censat", "censat")]:
        if config.get(k, ""):
            names.append(name)
    return names

def all_annotation_inputs():
    """Return all configured annotation BED files including pclai."""
    return [config[k] for k in ["annot_genes", "annot_repeats", "annot_segdups", "annot_censat", "annot_pclai"] if config.get(k, "")]

def all_annotation_names():
    """Return all display names for configured annotations including pclai."""
    names = []
    for k, name in [("annot_genes", "genes"), ("annot_repeats", "repeats"), ("annot_segdups", "segdups"), ("annot_censat", "censat"), ("annot_pclai", "pclai")]:
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
        if config.get("annot_repeats", ""):
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.annot-repeats.png")
        if config.get("annot_pclai", ""):
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.annot-ancestry.png")
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.annot-pclai-summary.png")
        return outputs
    return []

def annotation_snp_outputs(callers=None):
    """Return annotation SNP heatmap outputs if annotations are configured.

    callers: list of caller types to include, e.g. ["call"], ["deepvariant"],
             ["freebayes"], ["pangenie"], ["deconstruct"], or None for all.
    Deconstruct VCFs have no FILTER field, so only "all" is generated for them.
    """
    if not annotation_inputs():
        return []
    if callers is None:
        callers = ["deconstruct", "call", "deepvariant", "freebayes", "pangenie"]
    outputs = []
    for plot in ["annot-snp-counts", "annot-snp-tstv"]:
        if "deconstruct" in callers:
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.{plot}.all.png")
        if "call" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.{plot}.{filt}.png")
                outputs.append(f"{OUT_DIR}/merged.call.{plot}.{filt}.png")
                for s in LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.{plot}.{filt}.png")
                if LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/merged.longread.call.{plot}.{filt}.png")
        if "deepvariant" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.deepvariant.{plot}.{filt}.png")
                outputs.append(f"{OUT_DIR}/merged.deepvariant.{plot}.{filt}.png")
                for s in LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.deepvariant.{plot}.{filt}.png")
                if LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/merged.longread.deepvariant.{plot}.{filt}.png")
        if "freebayes" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.freebayes.{plot}.{filt}.png")
                outputs.append(f"{OUT_DIR}/merged.freebayes.{plot}.{filt}.png")
                if freebayes_longread_enabled():
                    for s in LR_SAMPLES:
                        outputs.append(f"{OUT_DIR}/{s}.freebayes.{plot}.{filt}.png")
                    if LR_SAMPLES:
                        outputs.append(f"{OUT_DIR}/merged.longread.freebayes.{plot}.{filt}.png")
        if "pangenie" in callers and pangenie_enabled():
            for filt in ["all", "pass"]:
                for s in SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.pangenie.{plot}.{filt}.png")
                outputs.append(f"{OUT_DIR}/merged.pangenie.{plot}.{filt}.png")
    if config.get("annot_pclai", ""):
        for plot in ["annot-pclai-snp-counts", "annot-pclai-snp-tstv"]:
            if "deconstruct" in callers:
                outputs.append(f"{OUT_DIR}/{OUT_NAME}.{plot}.all.png")
            if "call" in callers:
                for filt in ["all", "pass"]:
                    for s in SAMPLES:
                        outputs.append(f"{OUT_DIR}/{s}.{plot}.{filt}.png")
                    outputs.append(f"{OUT_DIR}/merged.call.{plot}.{filt}.png")
                    for s in LR_SAMPLES:
                        outputs.append(f"{OUT_DIR}/{s}.{plot}.{filt}.png")
                    if LR_SAMPLES:
                        outputs.append(f"{OUT_DIR}/merged.longread.call.{plot}.{filt}.png")
            if "deepvariant" in callers:
                for filt in ["all", "pass"]:
                    for s in SAMPLES:
                        outputs.append(f"{OUT_DIR}/{s}.deepvariant.{plot}.{filt}.png")
                    outputs.append(f"{OUT_DIR}/merged.deepvariant.{plot}.{filt}.png")
                    for s in LR_SAMPLES:
                        outputs.append(f"{OUT_DIR}/{s}.deepvariant.{plot}.{filt}.png")
                    if LR_SAMPLES:
                        outputs.append(f"{OUT_DIR}/merged.longread.deepvariant.{plot}.{filt}.png")
            if "freebayes" in callers:
                for filt in ["all", "pass"]:
                    for s in SAMPLES:
                        outputs.append(f"{OUT_DIR}/{s}.freebayes.{plot}.{filt}.png")
                    outputs.append(f"{OUT_DIR}/merged.freebayes.{plot}.{filt}.png")
                    if freebayes_longread_enabled():
                        for s in LR_SAMPLES:
                            outputs.append(f"{OUT_DIR}/{s}.freebayes.{plot}.{filt}.png")
                        if LR_SAMPLES:
                            outputs.append(f"{OUT_DIR}/merged.longread.freebayes.{plot}.{filt}.png")
            if "pangenie" in callers and pangenie_enabled():
                for filt in ["all", "pass"]:
                    for s in SAMPLES:
                        outputs.append(f"{OUT_DIR}/{s}.pangenie.{plot}.{filt}.png")
                    outputs.append(f"{OUT_DIR}/merged.pangenie.{plot}.{filt}.png")
    return outputs

def annotation_stats_outputs(callers=None):
    """Return annotation-stratified VCF stats outputs when annotations are configured.

    callers: list of caller types to include, e.g. ["deconstruct"], ["call"],
             ["deepvariant"], ["freebayes"], ["pangenie"], or None for all.
    """
    if not annotation_inputs():
        return []
    if callers is None:
        callers = ["deconstruct", "call", "deepvariant", "freebayes", "pangenie"]
    outputs = []
    for suffix in ["variant-types-by-annot.png", "vcf-stats-by-annot.tsv"]:
        if "deconstruct" in callers:
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.sites.{suffix}")
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.variants.{suffix}")
        if "call" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES + LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.call.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/{s}.call.variants.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.call.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.call.variants.{filt}.{suffix}")
                if LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/merged.longread.call.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.longread.call.variants.{filt}.{suffix}")
        if "deepvariant" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES + LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.dv.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/{s}.dv.variants.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.dv.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.dv.variants.{filt}.{suffix}")
                if LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/merged.longread.dv.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.longread.dv.variants.{filt}.{suffix}")
        if "freebayes" in callers:
            for filt in ["all", "pass"]:
                fb_samples = SAMPLES + LR_SAMPLES if freebayes_longread_enabled() else SAMPLES
                for s in fb_samples:
                    outputs.append(f"{OUT_DIR}/{s}.fb.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/{s}.fb.variants.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.fb.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.fb.variants.{filt}.{suffix}")
                if LR_SAMPLES and freebayes_longread_enabled():
                    outputs.append(f"{OUT_DIR}/merged.longread.fb.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.longread.fb.variants.{filt}.{suffix}")
        if "pangenie" in callers and pangenie_enabled():
            for filt in ["all", "pass"]:
                for s in SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.pg.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/{s}.pg.variants.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.pg.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.pg.variants.{filt}.{suffix}")
    return outputs

def augref_annot_beds():
    """Return augref-space annotation BED files (excluding pclai)."""
    if not annotation_inputs():
        return []
    return [f"{OUT_DIR}/{OUT_NAME}.augref-annot-{n}.bed" for n in annotation_names()]

def all_augref_annot_beds():
    """Return all augref-space annotation BED files including pclai."""
    if not all_annotation_inputs():
        return []
    return [f"{OUT_DIR}/{OUT_NAME}.augref-annot-{n}.bed" for n in all_annotation_names()]


def polymorphism_outputs():
    """Return segment polymorphism table output."""
    return [f"{OUT_DIR}/{OUT_NAME}.segment-polymorphism.tsv"]

# GIAB genome stratification helpers
GIAB_STRAT_NAMES = ["easy", "segdup", "otherdifficult"]
GIAB_STRAT_DISPLAY = ["Easy", "Segdup", "Other_Difficult"]

def giab_strat_configured():
    return bool(config.get("giab_strat", ""))

def giab_strat_beds():
    """Return the 3 source GIAB partition BED paths from config prefix."""
    prefix = config.get("giab_strat", "")
    if not prefix:
        return []
    return [f"{prefix}-{n}.bed" for n in GIAB_STRAT_NAMES]

def augref_giab_strat_beds():
    """Return 3 augref-space GIAB partition BED paths (output of giab_strat_augref rule)."""
    if not giab_strat_configured():
        return []
    return [f"{OUT_DIR}/{OUT_NAME}.augref-giab-{n}.bed" for n in GIAB_STRAT_NAMES]

def call_giab_strat_beds():
    """Return 3 call-space GIAB partition BED paths (CHM13#0# prefix stripped).

    vg call uses the locus name (e.g. 'chr1') as CHROM, stripping the graph
    path prefix (REF#0#).  Source GIAB BEDs use graph path names
    (e.g. 'CHM13#0#chr1'), so we need a stripped copy.
    """
    if not giab_strat_configured():
        return []
    return [f"{OUT_DIR}/{OUT_NAME}.call-giab-{n}.bed" for n in GIAB_STRAT_NAMES]

def call_annot_beds():
    """Return call-space annotation BED files (augref_REF#0# prefix stripped).

    vg call uses plain chr names while augref annotation BEDs use augref path
    names (e.g. 'augref_CHM13#0#chr1').
    """
    if not annotation_inputs():
        return []
    return [f"{OUT_DIR}/{OUT_NAME}.call-annot-{n}.bed" for n in annotation_names()]

def giab_strat_stats_outputs(callers=None):
    """Return GIAB strat plot/TSV outputs for each stats rule when configured."""
    if not giab_strat_configured():
        return []
    if callers is None:
        callers = ["deconstruct", "call", "deepvariant", "freebayes", "pangenie"]
    outputs = []
    for suffix in ["giab-strat.png", "giab-strat.tsv"]:
        if "deconstruct" in callers:
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.sites.{suffix}")
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.variants.{suffix}")
        if "call" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES + LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.call.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/{s}.call.variants.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.call.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.call.variants.{filt}.{suffix}")
                if LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/merged.longread.call.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.longread.call.variants.{filt}.{suffix}")
        if "deepvariant" in callers:
            for filt in ["all", "pass"]:
                for s in SAMPLES + LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.dv.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/{s}.dv.variants.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.dv.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.dv.variants.{filt}.{suffix}")
                if LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/merged.longread.dv.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.longread.dv.variants.{filt}.{suffix}")
        if "freebayes" in callers:
            for filt in ["all", "pass"]:
                fb_samples = SAMPLES + LR_SAMPLES if freebayes_longread_enabled() else SAMPLES
                for s in fb_samples:
                    outputs.append(f"{OUT_DIR}/{s}.fb.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/{s}.fb.variants.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.fb.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.fb.variants.{filt}.{suffix}")
                if LR_SAMPLES and freebayes_longread_enabled():
                    outputs.append(f"{OUT_DIR}/merged.longread.fb.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.longread.fb.variants.{filt}.{suffix}")
        if "pangenie" in callers and pangenie_enabled():
            for filt in ["all", "pass"]:
                for s in SAMPLES:
                    outputs.append(f"{OUT_DIR}/{s}.pg.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/{s}.pg.variants.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.pg.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.pg.variants.{filt}.{suffix}")
    return outputs

def per_sample_stats_outputs(callers=None):
    """Return per-sample stats outputs for multisample VCF rules."""
    if callers is None:
        callers = ["deconstruct", "call", "deepvariant", "freebayes", "pangenie"]
    outputs = []
    for suffix in ["per-sample-types.png", "per-sample-types.tsv", "per-sample-sv-types.png"]:
        if "deconstruct" in callers:
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.sites.{suffix}")
            outputs.append(f"{OUT_DIR}/{OUT_NAME}.variants.{suffix}")
        if "call" in callers:
            for filt in ["all", "pass"]:
                outputs.append(f"{OUT_DIR}/merged.call.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.call.variants.{filt}.{suffix}")
                if LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/merged.longread.call.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.longread.call.variants.{filt}.{suffix}")
        if "deepvariant" in callers:
            for filt in ["all", "pass"]:
                outputs.append(f"{OUT_DIR}/merged.dv.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.dv.variants.{filt}.{suffix}")
                if LR_SAMPLES:
                    outputs.append(f"{OUT_DIR}/merged.longread.dv.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.longread.dv.variants.{filt}.{suffix}")
        if "freebayes" in callers:
            for filt in ["all", "pass"]:
                outputs.append(f"{OUT_DIR}/merged.fb.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.fb.variants.{filt}.{suffix}")
                if LR_SAMPLES and freebayes_longread_enabled():
                    outputs.append(f"{OUT_DIR}/merged.longread.fb.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.longread.fb.variants.{filt}.{suffix}")
        if "pangenie" in callers and pangenie_enabled():
            for filt in ["all", "pass"]:
                outputs.append(f"{OUT_DIR}/merged.pg.sites.{filt}.{suffix}")
                outputs.append(f"{OUT_DIR}/merged.pg.variants.{filt}.{suffix}")
    if giab_strat_configured():
        for suffix in ["per-sample-giab-strat.png", "per-sample-giab-strat.tsv"]:
            if "deconstruct" in callers:
                outputs.append(f"{OUT_DIR}/{OUT_NAME}.sites.{suffix}")
                outputs.append(f"{OUT_DIR}/{OUT_NAME}.variants.{suffix}")
            if "call" in callers:
                for filt in ["all", "pass"]:
                    outputs.append(f"{OUT_DIR}/merged.call.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.call.variants.{filt}.{suffix}")
                    if LR_SAMPLES:
                        outputs.append(f"{OUT_DIR}/merged.longread.call.sites.{filt}.{suffix}")
                        outputs.append(f"{OUT_DIR}/merged.longread.call.variants.{filt}.{suffix}")
            if "deepvariant" in callers:
                for filt in ["all", "pass"]:
                    outputs.append(f"{OUT_DIR}/merged.dv.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.dv.variants.{filt}.{suffix}")
                    if LR_SAMPLES:
                        outputs.append(f"{OUT_DIR}/merged.longread.dv.sites.{filt}.{suffix}")
                        outputs.append(f"{OUT_DIR}/merged.longread.dv.variants.{filt}.{suffix}")
            if "freebayes" in callers:
                for filt in ["all", "pass"]:
                    outputs.append(f"{OUT_DIR}/merged.fb.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.fb.variants.{filt}.{suffix}")
                    if LR_SAMPLES and freebayes_longread_enabled():
                        outputs.append(f"{OUT_DIR}/merged.longread.fb.sites.{filt}.{suffix}")
                        outputs.append(f"{OUT_DIR}/merged.longread.fb.variants.{filt}.{suffix}")
            if "pangenie" in callers and pangenie_enabled():
                for filt in ["all", "pass"]:
                    outputs.append(f"{OUT_DIR}/merged.pg.sites.{filt}.{suffix}")
                    outputs.append(f"{OUT_DIR}/merged.pg.variants.{filt}.{suffix}")
    return outputs

def compare_call_dv_outputs():
    """Return call-vs-DV comparison outputs when samples are configured."""
    if not SAMPLES:
        return []
    outputs = []
    for mode in ["sites", "variants"]:
        for filt in ["all", "pass"]:
            outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{mode}.{filt}.compare.png")
            outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{mode}.{filt}.compare.tsv")
    return outputs

def vcfeval_lr_compare_outputs():
    """Return vcfeval-based call-vs-DV comparison outputs for long-read samples."""
    if not LR_SAMPLES:
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-dv.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-dv.{filt}.vcfeval-compare.tsv")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-dv.{filt}.chromsplit.tsv")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-dv.{filt}.chromsplit.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-dv.{filt}.chromsplit-top.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-dv.{filt}.chromsplit-concordant.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-dv.{filt}.chromsplit-top-onref.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-dv.{filt}.chromsplit-concordant-onref.png")
        if annotation_inputs():
            outputs.append(f"{OUT_DIR}/merged.lr.call-vs-dv.{filt}.chromsplit-annot.png")
        if giab_strat_configured():
            outputs.append(f"{OUT_DIR}/merged.lr.call-vs-dv.{filt}.chromsplit-giab.png")
    return outputs

def compare_call_fb_outputs():
    """Return call-vs-FB comparison outputs when samples are configured."""
    if not SAMPLES:
        return []
    outputs = []
    for mode in ["sites", "variants"]:
        for filt in ["all", "pass"]:
            outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{mode}.{filt}.compare.png")
            outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{mode}.{filt}.compare.tsv")
    return outputs

def vcfeval_fb_compare_outputs():
    """Return vcfeval-based call-vs-FB comparison outputs when samples are configured."""
    if not SAMPLES:
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{filt}.vcfeval-compare.tsv")
        outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{filt}.chromsplit.tsv")
        outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{filt}.chromsplit.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{filt}.chromsplit-top.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{filt}.chromsplit-concordant.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{filt}.chromsplit-top-onref.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{filt}.chromsplit-concordant-onref.png")
        if annotation_inputs():
            outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{filt}.chromsplit-annot.png")
        if giab_strat_configured():
            outputs.append(f"{OUT_DIR}/merged.call-vs-fb.{filt}.chromsplit-giab.png")
    return outputs

def vcfeval_lr_fb_compare_outputs():
    """Return vcfeval-based call-vs-FB comparison outputs for long-read samples."""
    if not LR_SAMPLES or not freebayes_longread_enabled():
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-fb.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-fb.{filt}.vcfeval-compare.tsv")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-fb.{filt}.chromsplit.tsv")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-fb.{filt}.chromsplit.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-fb.{filt}.chromsplit-top.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-fb.{filt}.chromsplit-concordant.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-fb.{filt}.chromsplit-top-onref.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-fb.{filt}.chromsplit-concordant-onref.png")
        if annotation_inputs():
            outputs.append(f"{OUT_DIR}/merged.lr.call-vs-fb.{filt}.chromsplit-annot.png")
        if giab_strat_configured():
            outputs.append(f"{OUT_DIR}/merged.lr.call-vs-fb.{filt}.chromsplit-giab.png")
    return outputs

def vcfeval_dv_vs_fb_compare_outputs():
    """Return vcfeval-based DV-vs-FB comparison outputs when samples are configured."""
    if not SAMPLES:
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.dv-vs-fb.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-fb.{filt}.vcfeval-compare.tsv")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-fb.{filt}.chromsplit.tsv")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-fb.{filt}.chromsplit.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-fb.{filt}.chromsplit-top.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-fb.{filt}.chromsplit-concordant.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-fb.{filt}.chromsplit-top-onref.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-fb.{filt}.chromsplit-concordant-onref.png")
        if annotation_inputs():
            outputs.append(f"{OUT_DIR}/merged.dv-vs-fb.{filt}.chromsplit-annot.png")
        if giab_strat_configured():
            outputs.append(f"{OUT_DIR}/merged.dv-vs-fb.{filt}.chromsplit-giab.png")
    return outputs

def vcfeval_lr_dv_vs_fb_compare_outputs():
    """Return vcfeval-based DV-vs-FB comparison outputs for long-read samples."""
    if not LR_SAMPLES or not freebayes_longread_enabled():
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-fb.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-fb.{filt}.vcfeval-compare.tsv")
        outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-fb.{filt}.chromsplit.tsv")
        outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-fb.{filt}.chromsplit.png")
        outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-fb.{filt}.chromsplit-top.png")
        outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-fb.{filt}.chromsplit-concordant.png")
        outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-fb.{filt}.chromsplit-top-onref.png")
        outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-fb.{filt}.chromsplit-concordant-onref.png")
        if annotation_inputs():
            outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-fb.{filt}.chromsplit-annot.png")
        if giab_strat_configured():
            outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-fb.{filt}.chromsplit-giab.png")
    return outputs

def vcfeval_bc_compare_outputs():
    """Return vcfeval-based call-vs-bcftools comparison outputs when samples are configured."""
    if not SAMPLES or not bcftools_enabled():
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.call-vs-bc.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-bc.{filt}.vcfeval-compare.tsv")
    return outputs

def vcfeval_lr_bc_compare_outputs():
    """Return vcfeval-based call-vs-bcftools comparison outputs for long-read samples."""
    if not LR_SAMPLES or not bcftools_enabled() or not bcftools_longread_enabled():
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-bc.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.lr.call-vs-bc.{filt}.vcfeval-compare.tsv")
    return outputs

def vcfeval_dv_vs_bc_compare_outputs():
    """Return vcfeval-based DV-vs-bcftools comparison outputs when samples are configured."""
    if not SAMPLES or not bcftools_enabled():
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.dv-vs-bc.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-bc.{filt}.vcfeval-compare.tsv")
    return outputs

def vcfeval_lr_dv_vs_bc_compare_outputs():
    """Return vcfeval-based DV-vs-bcftools comparison outputs for long-read samples."""
    if not LR_SAMPLES or not bcftools_enabled() or not bcftools_longread_enabled():
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-bc.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.lr.dv-vs-bc.{filt}.vcfeval-compare.tsv")
    return outputs

def compare_call_pg_outputs():
    """Return call-vs-PG comparison outputs when samples are configured."""
    if not SAMPLES or not pangenie_enabled():
        return []
    outputs = []
    for mode in ["sites", "variants"]:
        for filt in ["all", "pass"]:
            outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{mode}.{filt}.compare.png")
            outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{mode}.{filt}.compare.tsv")
    return outputs

def vcfeval_pg_compare_outputs():
    """Return vcfeval-based call-vs-PG comparison outputs when samples are configured."""
    if not SAMPLES or not pangenie_enabled():
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{filt}.vcfeval-compare.tsv")
        outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{filt}.chromsplit.tsv")
        outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{filt}.chromsplit.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{filt}.chromsplit-top.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{filt}.chromsplit-concordant.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{filt}.chromsplit-top-onref.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{filt}.chromsplit-concordant-onref.png")
        if annotation_inputs():
            outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{filt}.chromsplit-annot.png")
        if giab_strat_configured():
            outputs.append(f"{OUT_DIR}/merged.call-vs-pg.{filt}.chromsplit-giab.png")
    return outputs

def vcfeval_lr_pg_compare_outputs():
    """PanGenie does not support long reads — always returns empty list."""
    return []

def vcfeval_dv_vs_pg_compare_outputs():
    """Return vcfeval-based DV-vs-PG comparison outputs when samples are configured."""
    if not SAMPLES or not pangenie_enabled():
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.dv-vs-pg.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-pg.{filt}.vcfeval-compare.tsv")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-pg.{filt}.chromsplit.tsv")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-pg.{filt}.chromsplit.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-pg.{filt}.chromsplit-top.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-pg.{filt}.chromsplit-concordant.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-pg.{filt}.chromsplit-top-onref.png")
        outputs.append(f"{OUT_DIR}/merged.dv-vs-pg.{filt}.chromsplit-concordant-onref.png")
        if annotation_inputs():
            outputs.append(f"{OUT_DIR}/merged.dv-vs-pg.{filt}.chromsplit-annot.png")
        if giab_strat_configured():
            outputs.append(f"{OUT_DIR}/merged.dv-vs-pg.{filt}.chromsplit-giab.png")
    return outputs

def vcfeval_lr_dv_vs_pg_compare_outputs():
    """PanGenie does not support long reads — always returns empty list."""
    return []

def vcfeval_compare_outputs():
    """Return vcfeval-based call-vs-DV comparison outputs when samples are configured."""
    if not SAMPLES:
        return []
    outputs = []
    for filt in ["all", "pass"]:
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-compare.tsv")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-squash.vcfeval-compare.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-squash.vcfeval-compare.tsv")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.chromsplit.tsv")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.chromsplit.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.chromsplit-top.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.chromsplit-concordant.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.chromsplit-top-onref.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.chromsplit-concordant-onref.png")
        if annotation_inputs():
            outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.chromsplit-annot.png")
        if giab_strat_configured():
            outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.chromsplit-giab.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-squash.chromsplit.tsv")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-squash.chromsplit.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-squash.chromsplit-top.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-squash.chromsplit-concordant.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-squash.chromsplit-top-onref.png")
        outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-squash.chromsplit-concordant-onref.png")
        if annotation_inputs():
            outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-squash.chromsplit-annot.png")
        if giab_strat_configured():
            outputs.append(f"{OUT_DIR}/merged.call-vs-dv.{filt}.vcfeval-squash.chromsplit-giab.png")
    return outputs

def merge_chromsplit_tsv(input_files, sample_names, output_file):
    """Pivot wide chromsplit TSVs to long format with sample column.

    Each per-sample chromsplit.tsv has contigs as columns and metrics as rows.
    This merges them into one TSV with columns:
      sample, contig, SNP_FP, SNP_FN, INDEL_FP, INDEL_FN, [SV cols], total_errors
    """
    with open(output_file, "w") as out:
        header_written = False
        for sample, f in zip(sample_names, input_files):
            with open(f) as fh:
                rows = [line.strip().split("\t") for line in fh]
            col_header = rows[0]  # ['type', '_total_', 'contig1', ...]
            data = {r[0]: r[1:] for r in rows[1:]}
            metrics = [r[0] for r in rows[1:] if r[0] != "ERRORS"]
            contigs = col_header[1:]
            if not header_written:
                cols = (["sample", "contig"]
                        + [m.replace("-", "_") for m in metrics]
                        + ["total_errors"])
                out.write("\t".join(cols) + "\n")
                header_written = True
            for i, contig in enumerate(contigs):
                if contig == "_total_":
                    continue
                values = [int(data.get(m, ["0"] * len(contigs))[i])
                          for m in metrics]
                total = sum(v for m, v in zip(metrics, values) if "TP" not in m)
                out.write("\t".join(
                    [sample, contig] + [str(v) for v in values]
                    + [str(total)]) + "\n")

def pantree_outputs():
    """Return pantree comparison outputs when pantree_vcf is configured."""
    if not config.get("pantree_vcf", ""):
        return []
    return [
        f"{OUT_DIR}/pantree.variant-types.png",
        f"{OUT_DIR}/pantree.density.png",
        f"{OUT_DIR}/{OUT_NAME}.pantree-types.png",
        f"{OUT_DIR}/{OUT_NAME}.pantree-types-pct.png",
        f"{OUT_DIR}/{OUT_NAME}.pantree-size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.pantree-af.png",
        f"{OUT_DIR}/{OUT_NAME}.pantree-compare.tsv",
    ]

def summary_figure_outputs():
    """Return numbered summary figure outputs based on pipeline scope."""
    outputs = [
        f"{OUT_DIR}/1.augref-summary.png",
        f"{OUT_DIR}/2.deconstruct-summary.png",
    ]
    if SAMPLES:
        outputs.append(f"{OUT_DIR}/3.call-summary.png")
        outputs.append(f"{OUT_DIR}/4.deepvariant-summary.png")
        outputs.append(f"{OUT_DIR}/4b.freebayes-summary.png")
        if pangenie_enabled():
            outputs.append(f"{OUT_DIR}/4c.pangenie-summary.png")
        if bcftools_enabled():
            outputs.append(f"{OUT_DIR}/4d.bcftools-summary.png")
        outputs.append(f"{OUT_DIR}/5.concordance-summary.png")
        outputs.append(f"{OUT_DIR}/5b.concordance-onref-summary.png")
        outputs.append(f"{OUT_DIR}/5c.coverage-summary.png")
        outputs.append(f"{OUT_DIR}/5d.mapq-summary.png")
        outputs.append(f"{OUT_DIR}/10.freebayes-concordance-summary.png")
        outputs.append(f"{OUT_DIR}/10b.freebayes-concordance-onref-summary.png")
        if pangenie_enabled():
            outputs.append(f"{OUT_DIR}/12.pangenie-concordance-summary.png")
            outputs.append(f"{OUT_DIR}/12b.pangenie-concordance-onref-summary.png")
    if LR_SAMPLES:
        outputs.append(f"{OUT_DIR}/7.call-summary-longread.png")
        outputs.append(f"{OUT_DIR}/8.deepvariant-summary-longread.png")
        if freebayes_longread_enabled():
            outputs.append(f"{OUT_DIR}/8b.freebayes-summary-longread.png")
        if bcftools_enabled() and bcftools_longread_enabled():
            outputs.append(f"{OUT_DIR}/8d.bcftools-summary-longread.png")
        outputs.append(f"{OUT_DIR}/9.concordance-summary-longread.png")
        outputs.append(f"{OUT_DIR}/9b.concordance-onref-summary-longread.png")
        outputs.append(f"{OUT_DIR}/9c.coverage-summary-longread.png")
        outputs.append(f"{OUT_DIR}/9d.mapq-summary-longread.png")
        if freebayes_longread_enabled():
            outputs.append(f"{OUT_DIR}/11.freebayes-concordance-summary-longread.png")
            outputs.append(f"{OUT_DIR}/11b.freebayes-concordance-onref-summary-longread.png")
    if config.get("pantree_vcf", ""):
        outputs.append(f"{OUT_DIR}/6.pantree-summary.png")
    return outputs

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
        *giab_strat_stats_outputs(),
        *per_sample_stats_outputs(),
        *compare_call_dv_outputs(),
        *vcfeval_compare_outputs(),
        *vcfeval_lr_compare_outputs(),
        *compare_call_fb_outputs(),
        *vcfeval_fb_compare_outputs(),
        *vcfeval_lr_fb_compare_outputs(),
        *([f"{OUT_DIR}/merged.call-vs-fb.relaxed.vcfeval-compare.png"] if SAMPLES else []),
        *([f"{OUT_DIR}/merged.lr.call-vs-fb.relaxed.vcfeval-compare.png"] if LR_SAMPLES and freebayes_longread_enabled() else []),
        *vcfeval_dv_vs_fb_compare_outputs(),
        *vcfeval_lr_dv_vs_fb_compare_outputs(),
        *vcfeval_bc_compare_outputs(),
        *vcfeval_lr_bc_compare_outputs(),
        *vcfeval_dv_vs_bc_compare_outputs(),
        *vcfeval_lr_dv_vs_bc_compare_outputs(),
        *compare_call_pg_outputs(),
        *vcfeval_pg_compare_outputs(),
        *vcfeval_lr_pg_compare_outputs(),
        *vcfeval_dv_vs_pg_compare_outputs(),
        *vcfeval_lr_dv_vs_pg_compare_outputs(),
        *polymorphism_outputs(),
        *pantree_outputs(),
        *summary_figure_outputs(),
        # per-sample genotyping outputs
        expand("{out}/{s}.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.call-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.contig-depth.png", out=OUT_DIR, s=SAMPLES),
        # long-read per-sample outputs (giraffe -b hifi + call + DV PACBIO + depth)
        expand("{out}/{s}.vcf.gz", out=OUT_DIR, s=LR_SAMPLES),
        expand("{out}/{s}.call-offref.png", out=OUT_DIR, s=LR_SAMPLES),
        expand("{out}/{s}.contig-depth.png", out=OUT_DIR, s=LR_SAMPLES),
        expand("{out}/{s}.deepvariant.vcf.gz", out=OUT_DIR, s=LR_SAMPLES),
        expand("{out}/{s}.dv-offref.png", out=OUT_DIR, s=LR_SAMPLES),
        *(expand("{out}/{s}.freebayes.vcf.gz", out=OUT_DIR, s=LR_SAMPLES)
          + expand("{out}/{s}.fb-offref.png", out=OUT_DIR, s=LR_SAMPLES)
          if freebayes_longread_enabled() else []),
        expand("{out}/{s}.call.{mode}.{filt}.vcf-stats.tsv", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"]),
        expand("{out}/{s}.call.{mode}.{filt}.variant-types.png", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"]),
        expand("{out}/{s}.call.{mode}.{filt}.size-dist.png", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"]),
        expand("{out}/{s}.call.{mode}.{filt}.size-dist-log.png", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"]),
        expand("{out}/{s}.dv.{mode}.{filt}.vcf-stats.tsv", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"]),
        expand("{out}/{s}.dv.{mode}.{filt}.variant-types.png", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"]),
        expand("{out}/{s}.dv.{mode}.{filt}.size-dist.png", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"]),
        expand("{out}/{s}.dv.{mode}.{filt}.size-dist-log.png", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"]),
        *(expand("{out}/{s}.fb.{mode}.{filt}.vcf-stats.tsv", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"])
          + expand("{out}/{s}.fb.{mode}.{filt}.variant-types.png", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"])
          + expand("{out}/{s}.fb.{mode}.{filt}.size-dist.png", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"])
          + expand("{out}/{s}.fb.{mode}.{filt}.size-dist-log.png", out=OUT_DIR, s=LR_SAMPLES, mode=["sites", "variants"], filt=["all", "pass"])
          if freebayes_longread_enabled() else []),
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
        # per-sample freebayes outputs
        expand("{out}/{s}.freebayes.vcf.gz", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.fb-offref.png", out=OUT_DIR, s=SAMPLES),
        expand("{out}/{s}.fb.sites.{filt}.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.fb.sites.{filt}.variant-types.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.fb.sites.{filt}.size-dist.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.fb.sites.{filt}.size-dist-log.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.fb.variants.{filt}.vcf-stats.tsv", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.fb.variants.{filt}.variant-types.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.fb.variants.{filt}.size-dist.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        expand("{out}/{s}.fb.variants.{filt}.size-dist-log.png", out=OUT_DIR, s=SAMPLES, filt=["all", "pass"]),
        # per-sample pangenie outputs (when enabled)
        *(expand("{out}/{s}.pangenie.vcf.gz", out=OUT_DIR, s=SAMPLES)
          + expand("{out}/{s}.pg-offref.png", out=OUT_DIR, s=SAMPLES)
          + expand("{out}/{s}.pg.{mode}.{filt}.{suffix}", out=OUT_DIR, s=SAMPLES,
                   mode=["sites", "variants"], filt=["all", "pass"],
                   suffix=["vcf-stats.tsv", "variant-types.png", "size-dist.png", "size-dist-log.png"])
          if pangenie_enabled() else []),
        # merged outputs
        f"{OUT_DIR}/merged.call.vcf.gz",
        f"{OUT_DIR}/merged.deepvariant.vcf.gz",
        f"{OUT_DIR}/merged.freebayes.vcf.gz",
        *([f"{OUT_DIR}/merged.pangenie.vcf.gz",
           f"{OUT_DIR}/merged.pg-offref.png"] if pangenie_enabled() else []),
        f"{OUT_DIR}/merged.call-offref.png",
        f"{OUT_DIR}/merged.dv-offref.png",
        f"{OUT_DIR}/merged.fb-offref.png",
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
        f"{OUT_DIR}/merged.fb.sites.all.vcf-stats.tsv",
        f"{OUT_DIR}/merged.fb.sites.all.variant-types.png",
        f"{OUT_DIR}/merged.fb.sites.all.size-dist.png",
        f"{OUT_DIR}/merged.fb.sites.all.size-dist-log.png",
        f"{OUT_DIR}/merged.fb.sites.all.af-spectrum.png",
        f"{OUT_DIR}/merged.fb.sites.pass.vcf-stats.tsv",
        f"{OUT_DIR}/merged.fb.sites.pass.variant-types.png",
        f"{OUT_DIR}/merged.fb.sites.pass.size-dist.png",
        f"{OUT_DIR}/merged.fb.sites.pass.size-dist-log.png",
        f"{OUT_DIR}/merged.fb.sites.pass.af-spectrum.png",
        f"{OUT_DIR}/merged.fb.variants.all.vcf-stats.tsv",
        f"{OUT_DIR}/merged.fb.variants.all.variant-types.png",
        f"{OUT_DIR}/merged.fb.variants.all.size-dist.png",
        f"{OUT_DIR}/merged.fb.variants.all.size-dist-log.png",
        f"{OUT_DIR}/merged.fb.variants.all.af-spectrum.png",
        f"{OUT_DIR}/merged.fb.variants.pass.vcf-stats.tsv",
        f"{OUT_DIR}/merged.fb.variants.pass.variant-types.png",
        f"{OUT_DIR}/merged.fb.variants.pass.size-dist.png",
        f"{OUT_DIR}/merged.fb.variants.pass.size-dist-log.png",
        f"{OUT_DIR}/merged.fb.variants.pass.af-spectrum.png",
        *([f"{OUT_DIR}/merged.pg.{mode}.{filt}.{suffix}"
           for mode in ["sites", "variants"]
           for filt in ["all", "pass"]
           for suffix in ["vcf-stats.tsv", "variant-types.png", "size-dist.png",
                          "size-dist-log.png", "af-spectrum.png"]]
          if pangenie_enabled() else []),
        # merged long-read outputs (when longread_samples configured)
        *([f"{OUT_DIR}/merged.longread.call.vcf.gz",
           f"{OUT_DIR}/merged.longread.deepvariant.vcf.gz"]
          + ([f"{OUT_DIR}/merged.longread.freebayes.vcf.gz"] if freebayes_longread_enabled() else [])
          + [f"{OUT_DIR}/merged.longread.{caller}.{mode}.{filt}.{suffix}"
             for caller in ["call", "dv"] + (["fb"] if freebayes_longread_enabled() else [])
             for mode in ["sites", "variants"]
             for filt in ["all", "pass"]
             for suffix in ["vcf-stats.tsv", "variant-types.png", "size-dist.png",
                            "size-dist-log.png", "af-spectrum.png"]]
          if LR_SAMPLES else []),

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
        *giab_strat_stats_outputs(["deconstruct"]),
        *per_sample_stats_outputs(["deconstruct"]),
        *polymorphism_outputs(),
        *pantree_outputs(),
        f"{OUT_DIR}/1.augref-summary.png",
        f"{OUT_DIR}/2.deconstruct-summary.png",

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
        *giab_strat_stats_outputs(["call"]),
        *per_sample_stats_outputs(["call"]),

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
        *giab_strat_stats_outputs(["deepvariant"]),
        *per_sample_stats_outputs(["deepvariant"]),

# Resolve wildcard ambiguities:
# - {OUT_NAME}.vcf.gz matches both deconstruct and call (sample={OUT_NAME})
# - {sample}.deepvariant.vcf.gz matches both deepvariant and call (sample=X.deepvariant)
ruleorder: deconstruct > call
ruleorder: deepvariant > call
ruleorder: freebayes > call
ruleorder: pangenie > call

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
        mem_mb=rule_mem_gb("paths", 1024) * 1024,
        runtime=rule_runtime("paths"),
    params:
        mem_gb=rule_mem_gb("paths", 1024),
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

rule annotate_tr:
    """Annotate indels with tandem repeat motifs"""
    input:
        vcf=f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.tr.vcf.gz",
    resources:
        mem_mb=int(rule_mem_gb("annotate_tr", 64)) * 1024,
        runtime=rule_runtime("annotate_tr", 120),
    shell:
        "python3 scripts/annotate-tr.py --vcf {input.vcf} --ref {input.ref}"
        " -o {output} && tabix -fp vcf {output}"

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
        f"{OUT_DIR}/{OUT_NAME}.fa.gz.fai",
        f"{OUT_DIR}/{OUT_NAME}.fa.gz.gzi",
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

rule filtered_paths:
    """GBZ → filtered augref path list (contigs >= min_surject_len)"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt",
    threads: 1
    resources:
        mem_mb=int(rule_mem_gb("filtered_paths", 64)) * 1024,
        runtime=30,
    shell:
        "vg paths -x {input} -S {AUGREF} -E"
        " | awk '$2 >= {config[min_surject_len]}' | cut -f1"
        " > {output}"

rule length_hist:
    """Augref segments → length histogram + log-log CCDF"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        hist=f"{OUT_DIR}/{OUT_NAME}.augref-length-hist.png",
        loglog=f"{OUT_DIR}/{OUT_NAME}.augref-length-hist-loglog.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/offref-length-hist.R {output.hist} {input}"

rule augref_depth:
    """Per-segment haplotype depth on augref _alt paths only.
    Wraps scripts/augref-depth.sh, which runs vg depth once per main
    chromosome (with trailing-underscore prefix so main refs are skipped),
    in parallel across chromosomes. Total threads / parallel = threads
    per vg depth invocation (default 120 / 10 = 12). Each vg depth load
    of an HPRC GBZ is ~13 GB; size mem / parallel accordingly."""
    input:
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.augref-depth.tsv",
    threads: rule_cpus("augref_depth", 120)
    resources:
        mem_mb=rule_mem_gb("augref_depth", 500) * 1024,
        runtime=rule_runtime("augref_depth"),
    params:
        vg=config.get("vg_bin", "vg"),
        parallel=config.get("augref_depth_parallel", 10),
    shell:
        "scripts/augref-depth.sh"
        " --gbz {input.gbz} --ref {REF} --out {output}"
        " --vg {params.vg} --threads {threads} --parallel {params.parallel}"

rule augref_frequency_plot:
    """Augref-depth TSV → population-frequency plots (histogram + length scatter).
    Filters to off-reference _alt contigs only; main chromosome paths
    saturate the depth and aren't informative for frequency."""
    input:
        f"{OUT_DIR}/{OUT_NAME}.augref-depth.tsv",
    output:
        hist=f"{OUT_DIR}/{OUT_NAME}.augref-frequency.png",
        scatter=f"{OUT_DIR}/{OUT_NAME}.augref-frequency-by-size.png",
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/augref-frequency.R"
        " {input} {OUT_DIR}/{OUT_NAME}"

rule segment_density:
    """Augref segments → off-reference segment density ideogram"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{OUT_NAME}.offref-segs.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_args=" ".join(
            [f"--censat '{config['annot_censat']}'" if config.get("annot_censat", "") else "",
             f"--segdups '{config['annot_segdups']}'" if config.get("annot_segdups", "") else "",
             f"--genes '{config['annot_genes']}'" if config.get("annot_genes", "") else ""]),
    shell:
        "ulimit -s unlimited && Rscript scripts/chrom-density-tsv.R"
        " {input} {output}"
        " '{REF} Off-Reference Segment Density'"
        " {config[min_augref_len]} '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF}"
        " {params.annot_args}"

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
    params:
        annot_args=" ".join(
            [f"--censat '{config['annot_censat']}'" if config.get("annot_censat", "") else "",
             f"--segdups '{config['annot_segdups']}'" if config.get("annot_segdups", "") else "",
             f"--genes '{config['annot_genes']}'" if config.get("annot_genes", "") else ""]),
    shell:
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} Off-Reference Variant Density'"
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF} --offref"
        " {params.annot_args}"

############################################################################
# Annotation overlap rules (optional — only when annot_* keys are set)
############################################################################

rule annotation_intersect:
    """Augref segments + annotation BEDs → per-segment annotation TSV"""
    input:
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annots=all_annotation_inputs(),
    output:
        f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv",
    resources:
        mem_mb=512000,
        runtime=2880,
    params:
        names=" ".join(all_annotation_names()),
        group_arg="--group-by-column 6" if config.get("annot_repeats", "") or config.get("annot_pclai", "") else "",
    shell:
        "python scripts/intersect-annotations.py"
        " {input.segs} {input.annots}"
        " --per-segment --per-segment-output {output}"
        " --output-dir {OUT_DIR}"
        " --annotation-names {params.names}"
        " {params.group_arg}"

rule annotation_genome_coverage:
    """Compute genome-wide annotation coverage from FAI + annotation BEDs.

    Annotation BEDs are already sorted and merged (per-class for grouped
    annotations like repeats) by the download script. We just run
    bedtools coverage against a genome BED built from the FAI.
    """
    input:
        fai=f"{OUT_DIR}/{OUT_NAME}.fa.gz.fai",
        annots=annotation_inputs(),
    output:
        f"{OUT_DIR}/{OUT_NAME}.annot-genome-coverage.tsv",
    resources:
        mem_mb=8000,
        runtime=240,
    params:
        names=annotation_names(),
        group_col=6 if config.get("annot_repeats", "") else 0,
    run:
        import os, subprocess, tempfile

        fai = str(input.fai)
        annot_files = list(input.annots)
        names = list(params.names)
        gc = int(params.group_col)

        # Build contig filter file and sum genome size from FAI (on-ref only)
        genome_bp = 0
        contigs_file = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
        with open(fai) as f:
            for line in f:
                fields = line.strip().split('\t')
                chrom, length = fields[0], int(fields[1])
                if chrom.endswith('_alt'):
                    continue
                if chrom.startswith('augref_'):
                    chrom = chrom[len('augref_'):]
                contigs_file.write(chrom + '\n')
                genome_bp += length
        contigs_file.close()

        rows = []
        for i, annot_file in enumerate(annot_files):
            name = names[i]
            is_grouped = gc > 0 and name in ('repeats', 'pclai')

            if is_grouped:
                # Class-agnostic total from .total.bed or awk on main file
                total_bed = os.path.splitext(annot_file)[0] + '.total.bed'
                if os.path.isfile(total_bed):
                    cmd = f"grep -Ff '{contigs_file.name}' '{total_bed}' | awk -F'\\t' '{{s+=$3-$2}} END{{print s+0}}'"
                else:
                    cmd = f"grep -Ff '{contigs_file.name}' '{annot_file}' | cut -f1-3 | awk -F'\\t' '{{s+=$3-$2}} END{{print s+0}}'"
                overlap = int(subprocess.check_output(cmd, shell=True, text=True).strip())
                rows.append(f'{name}\t_total\t{genome_bp}\t{overlap}\t{overlap/genome_bp:.6f}')

                # Per-class: single awk pass to sum bp per class, filtered to ref contigs
                cmd = (f"grep -Ff '{contigs_file.name}' '{annot_file}'"
                       f" | awk -F'\\t' '{{class=$" + str(gc) + "; bp[$" + str(gc) + "]+=$3-$2}"
                       f" END{{for(c in bp) print c\"\\t\"bp[c]}}'")
                result = subprocess.check_output(cmd, shell=True, text=True).strip()
                for line in result.split('\n'):
                    if line:
                        cls, bp = line.split('\t')
                        bp = int(bp)
                        rows.append(f'{name}\t{cls}\t{genome_bp}\t{bp}\t{bp/genome_bp:.6f}')
            else:
                cmd = f"grep -Ff '{contigs_file.name}' '{annot_file}' | awk -F'\\t' '{{s+=$3-$2}} END{{print s+0}}'"
                overlap = int(subprocess.check_output(cmd, shell=True, text=True).strip())
                rows.append(f'{name}\t{name}\t{genome_bp}\t{overlap}\t{overlap/genome_bp:.6f}')

        os.remove(contigs_file.name)

        # Sort rows by annotation name then class for consistent output
        rows.sort()
        with open(str(output[0]), 'w') as out:
            out.write('annotation\tannotation_class\tgenome_bp\toverlap_bp\toverlap_frac\n')
            for r in rows:
                out.write(r + '\n')

rule annotation_augref_beds:
    """Per-segment annotation TSV + annotation BEDs → augref-space BEDs"""
    input:
        seg_annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv",
        annots=all_annotation_inputs(),
    output:
        all_augref_annot_beds(),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        beds=lambda wc, input: ",".join(input.annots),
        names=",".join(all_annotation_names()),
    shell:
        "python scripts/make-augref-annotation-beds.py"
        " --per-segment-tsv {input.seg_annot}"
        " --annotation-beds {params.beds}"
        " --annotation-names {params.names}"
        " --ref {REF} --min-overlap 0.5"
        " --output-dir {OUT_DIR} --output-prefix {OUT_NAME}"

rule giab_strat_augref:
    """GIAB partition BEDs → augref-space BEDs (reference-only, add augref_ prefix)"""
    input:
        beds=giab_strat_beds(),
    output:
        augref_giab_strat_beds(),
    resources:
        mem_mb=8000,
        runtime=60,
    params:
        in_str=lambda wc, input: " ".join(input.beds),
        out_str=lambda wc, output: " ".join(output),
    shell:
        r"""
        IN=({params.in_str})
        OUT=({params.out_str})
        for i in "${{!IN[@]}}"; do
            awk 'BEGIN{{OFS="\t"}} {{$1="augref_"$1; print}}' "${{IN[$i]}}" \
              | LC_ALL=C sort -k1,1 -k2,2n > "${{OUT[$i]}}"
        done
        """

rule giab_strat_call_space:
    """GIAB partition BEDs → call-space BEDs (strip REF#0# prefix for vg call VCFs)"""
    input:
        beds=giab_strat_beds(),
    output:
        call_giab_strat_beds(),
    resources:
        mem_mb=8000,
        runtime=60,
    params:
        in_str=lambda wc, input: " ".join(input.beds),
        out_str=lambda wc, output: " ".join(output),
        prefix=f"{REF}#0#",
    shell:
        r"""
        IN=({params.in_str})
        OUT=({params.out_str})
        for i in "${{!IN[@]}}"; do
            awk -v prefix="{params.prefix}" 'BEGIN{{OFS="\t"}} {{sub("^"prefix, "", $1); print}}' "${{IN[$i]}}" \
              | LC_ALL=C sort -k1,1 -k2,2n > "${{OUT[$i]}}"
        done
        """

rule annot_beds_call_space:
    """Augref annotation BEDs → call-space BEDs (strip augref_REF#0# prefix)"""
    input:
        beds=augref_annot_beds(),
    output:
        call_annot_beds(),
    resources:
        mem_mb=8000,
        runtime=60,
    params:
        in_str=lambda wc, input: " ".join(input.beds),
        out_str=lambda wc, output: " ".join(output),
        prefix=f"augref_{REF}#0#",
    shell:
        r"""
        IN=({params.in_str})
        OUT=({params.out_str})
        for i in "${{!IN[@]}}"; do
            awk -v prefix="{params.prefix}" 'BEGIN{{OFS="\t"}} {{sub("^"prefix, "", $1); print}}' "${{IN[$i]}}" \
              | LC_ALL=C sort -k1,1 -k2,2n > "${{OUT[$i]}}"
        done
        """

rule annotation_plots:
    """Per-segment annotation TSV → overlap plots + stats"""
    input:
        per_seg=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv",
        genome_cov=f"{OUT_DIR}/{OUT_NAME}.annot-genome-coverage.tsv",
    output:
        annotation_outputs(),
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/annotation-plots.R {input.per_seg} {OUT_DIR}/{OUT_NAME}"
        " --min-overlap 0.5 --title '{REF} Annotation Overlap'"
        " --genome-coverage {input.genome_cov}"

rule pass_filter_vcf:
    """VCF → PASS-only VCF (filter before merging to avoid cross-sample filter contamination)"""
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.pass-only.vcf.gz",
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        "bcftools view -f PASS '{input}' -Oz -o {output}"
        " && tabix -p vcf {output}"

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
        "ulimit -s unlimited && Rscript scripts/annotation-plots.R {input.annot}"
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
        "ulimit -s unlimited && Rscript scripts/segment-polymorphism.R"
        " --vcf {input.vcf}"
        " --augref-prefix '{AUGREF}#0#'"
        " {params.annot_arg}"
        " --output {output}"

############################################################################
# Per-sample rules (wildcard: {sample})
############################################################################

def _get_reads(wc):
    """Look up reads index for a sample (short-read or long-read)."""
    if wc.sample in config.get("samples", {}):
        return config["samples"][wc.sample]
    return config["longread_samples"][wc.sample]

def _giraffe_preset(wc):
    """Return --preset flag for long-read samples, empty for short-read."""
    if wc.sample in config.get("longread_samples", {}):
        return "--preset hifi"
    return ""

def _giraffe_mem_gb(wc):
    """Memory for giraffe: 900 GB for long reads (HiFi index build), 512 GB for short reads."""
    default = 900 if wc.sample in config.get("longread_samples", {}) else 512
    return rule_mem_gb("giraffe", default)

rule giraffe:
    """GBZ + reads → GAM"""
    input:
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
        hapl=f"{OUT_DIR}/{OUT_NAME}.hapl",
        reads=_get_reads,
    output:
        f"{OUT_DIR}/{{sample}}.gam",
    threads: rule_cpus("giraffe", 128)
    resources:
        mem_mb=lambda wc: _giraffe_mem_gb(wc) * 1024,
        runtime=rule_runtime("giraffe"),
    params:
        mem_gb=_giraffe_mem_gb,
        preset=_giraffe_preset,
    shell:
        "scripts/giraffe.sh"
        " --gbz {input.gbz}"
        " --hapl {input.hapl}"
        " --reads {input.reads}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.gam"
        " {params.preset}"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule call:
    """GAM → VCF (vg call on all augref contigs, unfiltered)"""
    input:
        gam=f"{OUT_DIR}/{{sample}}.gam",
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        vcf=f"{OUT_DIR}/{{sample}}.vcf.gz",
        pack=f"{OUT_DIR}/{{sample}}.pack",
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

rule contig_depth:
    """Pack + GBZ → per-contig mean depth TSV"""
    input:
        pack=f"{OUT_DIR}/{{sample}}.pack",
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{{sample}}.contig-depth.tsv",
    threads: 4
    resources:
        mem_mb=256000,
        runtime=120,
    shell:
        "vg depth -k {input.pack} -b 1000000000 -P {AUGREF} -t {threads} {input.gbz}"
        " > {output}"

rule contig_depth_plot:
    """Per-sample augref contig depth scatter plot"""
    input:
        depth=f"{OUT_DIR}/{{sample}}.contig-depth.tsv",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{{sample}}.contig-depth.png",
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "ulimit -s unlimited && Rscript scripts/contig-depth-plot.R"
        " --depth {input.depth}"
        " --segs {input.segs}"
        " --sample {wildcards.sample}"
        " --output {output}"

rule gam_mapq:
    """GAM → annotated MAPQ distribution TSV (on-ref vs off-ref)"""
    input:
        gam=f"{OUT_DIR}/{{sample}}.gam",
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
    output:
        f"{OUT_DIR}/{{sample}}.gam-mapq.tsv",
    threads: rule_cpus("gam_mapq", 16)
    resources:
        mem_mb=rule_mem_gb("gam_mapq", 256) * 1024,
        runtime=rule_runtime("gam_mapq"),
    shell:
        "vg annotate -a {input.gam} -x {input.gbz} -p -m -t {threads}"
        " | vg view -aj -"
        " | python3 scripts/extract-mapq.py --mode gam --output {output}"

rule bam_mapq:
    """BAM → MAPQ distribution TSV (on-ref vs off-ref)"""
    input:
        f"{OUT_DIR}/{{sample}}.bam",
    output:
        f"{OUT_DIR}/{{sample}}.bam-mapq.tsv",
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "samtools view -F 4 {input} | cut -f3,5"
        " | python3 scripts/extract-mapq.py --mode bam --output {output}"

rule mapq_dist_plot:
    """MAPQ distribution summary: GAM vs BAM, on-ref vs off-ref"""
    input:
        gam_mapq=expand("{out}/{s}.gam-mapq.tsv", out=OUT_DIR, s=SAMPLES),
        bam_mapq=expand("{out}/{s}.bam-mapq.tsv", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/mapq-dist.png",
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        gam_arg=lambda wc, input: "--gam-mapq " + ",".join(input.gam_mapq),
        bam_arg=lambda wc, input: "--bam-mapq " + ",".join(input.bam_mapq),
    shell:
        "ulimit -s unlimited && Rscript scripts/mapq-dist-plot.R"
        " {params.gam_arg}"
        " {params.bam_arg}"
        " --output {output}"
        " --title '{REF} Mapping Quality Distribution'"

rule bam_contig_depth:
    """Surjected BAM → per-contig depth TSV (samtools coverage)"""
    input:
        f"{OUT_DIR}/{{sample}}.bam",
    output:
        f"{OUT_DIR}/{{sample}}.bam-depth.tsv",
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "samtools coverage {input} > {output}"

rule bam_contig_depth_q5:
    """Surjected BAM → per-contig depth TSV with MAPQ >= 5 filter"""
    input:
        f"{OUT_DIR}/{{sample}}.bam",
    output:
        f"{OUT_DIR}/{{sample}}.bam-depth-q5.tsv",
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        "samtools coverage -q 5 {input} > {output}"

rule contig_depth_summary:
    """Averaged pack + BAM depth across all samples → multi-panel summary"""
    input:
        pack_depths=expand("{out}/{s}.contig-depth.tsv", out=OUT_DIR, s=SAMPLES),
        bam_depths=expand("{out}/{s}.bam-depth.tsv", out=OUT_DIR, s=SAMPLES),
        bam_q5_depths=expand("{out}/{s}.bam-depth-q5.tsv", out=OUT_DIR, s=SAMPLES),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/contig-depth-summary.png",
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        pack_arg=lambda wc, input: "--pack-depths " + ",".join(input.pack_depths),
        bam_arg=lambda wc, input: "--bam-depths " + ",".join(input.bam_depths),
        bam_q5_arg=lambda wc, input: "--bam-q5-depths " + ",".join(input.bam_q5_depths),
        min_surject_len=config.get("min_surject_len", 0),
    shell:
        "ulimit -s unlimited && Rscript scripts/contig-depth-summary.R"
        " {params.pack_arg}"
        " {params.bam_arg}"
        " {params.bam_q5_arg}"
        " --segs {input.segs}"
        " --output {output}"
        " --depth-cap 60"
        " --min-surject-len {params.min_surject_len}"
        " --title '{REF} Augref Contig Read Depth'"

rule filter_call_vcf:
    """Filter call VCF to contigs >= min_surject_len (for DV comparison)"""
    input:
        vcf=f"{OUT_DIR}/{{sample}}.vcf.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt",
    output:
        f"{OUT_DIR}/{{sample}}.filtered.vcf.gz",
    resources:
        mem_mb=8000,
        runtime=60,
    params:
        strip_prefix=f"{AUGREF}#0#",
    shell:
        "awk -v prefix='{params.strip_prefix}' -v OFS='\\t'"
        " '{{sub(prefix, \"\"); print $0, 1, 2147483647}}'"
        " {input.paths} > {output}.targets.tmp"
        " && bcftools view -T {output}.targets.tmp {input.vcf}"
        " | bgzip > {output}"
        " && tabix -fp vcf {output}"
        " && rm -f {output}.targets.tmp"

rule surject:
    """GAM → sorted BAM"""
    input:
        gam=f"{OUT_DIR}/{{sample}}.gam",
        gbz=f"{OUT_DIR}/{OUT_NAME}.gbz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        f"{OUT_DIR}/{{sample}}.bam",
    threads: rule_cpus("surject", 128)
    resources:
        mem_mb=rule_mem_gb("surject", 512) * 1024,
        runtime=rule_runtime("surject"),
    params:
        mem_gb=rule_mem_gb("surject", 512),
        paths_arg=lambda wc, input: f"--paths-file {input.paths}" if surject_filtering() else "",
        longread_args=lambda wc: "--no-interleaved --read-length long" if wc.sample in config.get("longread_samples", {}) else "",
    shell:
        "scripts/surject.sh"
        " --gbz {input.gbz}"
        " --gam {input.gam}"
        " --ref {AUGREF}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.bam"
        " {params.paths_arg}"
        " {params.longread_args}"
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
        model_type=lambda wc: "PACBIO" if wc.sample in config.get("longread_samples", {}) else "WGS",
    shell:
        "scripts/deepvariant.sh"
        " --bam {input.bam}"
        " --ref {input.ref}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.deepvariant.vcf.gz"
        " --dv-version {config[dv_version]}"
        " --model-type {params.model_type}"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --tmpdir {resources.tmpdir}"
        " --local"

rule freebayes:
    """BAM + FASTA → VCF via FreeBayes (parallel by region)"""
    input:
        bam=f"{OUT_DIR}/{{sample}}.bam",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
    output:
        f"{OUT_DIR}/{{sample}}.freebayes.vcf.gz",
    threads: rule_cpus("freebayes", 96)
    resources:
        mem_mb=rule_mem_gb("freebayes", 1024) * 1024,
        runtime=rule_runtime("freebayes"),
    params:
        mem_gb=rule_mem_gb("freebayes", 1024),
        region_size=config.get("freebayes_region_size", 1000000),
        extra_args=lambda wc: (config.get("freebayes_longread_args") or "--limit-coverage 100") if wc.sample in config.get("longread_samples", {}) else (config.get("freebayes_extra_args") or ""),
        # Default to native freebayes on PATH; set freebayes_docker in config
        # to pin a container if the binary isn't available on the node.
        docker_arg=lambda wc: f"--docker {config['freebayes_docker']}" if config.get("freebayes_docker") else "",
    shell:
        "scripts/freebayes.sh"
        " --bam {input.bam}"
        " --ref {input.ref}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.freebayes.vcf.gz"
        " --region-size {params.region_size}"
        " --extra-args '{params.extra_args}'"
        " {params.docker_arg}"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule bcftools:
    """BAM + FASTA → VCF via bcftools mpileup + call (parallel by region).
    Uses the bcftools binary from PATH (must be >= 1.20 for pacbio-ccs preset)."""
    input:
        bam=f"{OUT_DIR}/{{sample}}.bam",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
    output:
        f"{OUT_DIR}/{{sample}}.bcftools.vcf.gz",
    threads: rule_cpus("bcftools", 96)
    resources:
        mem_mb=rule_mem_gb("bcftools", 256) * 1024,
        runtime=rule_runtime("bcftools"),
    params:
        mem_gb=rule_mem_gb("bcftools", 256),
        region_size=config.get("bcftools_region_size", 10000000),
        extra_args=lambda wc: config.get("bcftools_longread_args", "") if wc.sample in config.get("longread_samples", {}) else config.get("bcftools_extra_args", ""),
        long_read_flag=lambda wc: "--long-read" if wc.sample in config.get("longread_samples", {}) else "",
    shell:
        "scripts/bcftools-call.sh"
        " --bam {input.bam}"
        " --ref {input.ref}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.bcftools.vcf.gz"
        " --region-size {params.region_size}"
        " {params.long_read_flag}"
        " --extra-args '{params.extra_args}'"
        " --cpus {threads} --mem {params.mem_gb}gb"
        " --local"

rule pangenie_prepare_panel:
    """Prepare panel VCF for PanGenie: remove haploid samples, convert to diploid
    phased, filter missing, add IDs, split to biallelic, merge overlapping variants."""
    input:
        vcf=f"{OUT_DIR}/{OUT_NAME}.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.pangenie-panel.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        exclude_samples=config.get("pangenie_exclude_samples", ""),
        alt_missing=config.get("pangenie_alt_missing", 0.9),
    shell:
        "TMPDIR=$(mktemp -d \"${{TMPDIR:-.}}/pg-prep.XXXXXX\")"
        " && trap 'rm -rf \"$TMPDIR\"' EXIT"
        " && gunzip -c {input.ref} > $TMPDIR/ref.fa"
        " && samtools faidx $TMPDIR/ref.fa"
        " && if [ -n '{params.exclude_samples}' ]; then"
        "      bcftools view -s ^{params.exclude_samples} {input.vcf};"
        "    else bcftools view {input.vcf}; fi"
        "    | awk 'BEGIN{{OFS=\"\\t\"}} /^#/{{print;next}}"
        "      {{for(i=10;i<=NF;i++){{g=$i; if(g==\".\")$i=\".|.\"; else $i=g\"|\"g}} print}}'"
        "    | python3 scripts/pangenie/prepare-vcf.py --missing 0.2"
        "      --alt-missing {params.alt_missing} 2>/dev/null"
        "    | python3 scripts/pangenie/add-ids.py 2>/dev/null"
        "    | bgzip > $TMPDIR/callset.vcf.gz"
        " && tabix -fp vcf $TMPDIR/callset.vcf.gz"
        " && bcftools norm -m- $TMPDIR/callset.vcf.gz > $TMPDIR/biallelic.vcf 2>/dev/null"
        " && python3 scripts/pangenie/merge_vcfs.py merge"
        "    -vcf $TMPDIR/biallelic.vcf -r $TMPDIR/ref.fa -ploidy 2"
        "    2>/dev/null"
        "    | bgzip > {output}"
        " && tabix -fp vcf {output}"

rule pangenie:
    """FASTQ + panel VCF + ref → VCF via PanGenie Docker"""
    input:
        reads=_get_reads,
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        panel=f"{OUT_DIR}/{OUT_NAME}.pangenie-panel.vcf.gz",
    output:
        f"{OUT_DIR}/{{sample}}.pangenie.vcf.gz",
    threads: rule_cpus("pangenie", 24)
    resources:
        mem_mb=max(rule_mem_gb("pangenie", 256), 16) * 1024,
        runtime=rule_runtime("pangenie"),
    params:
        mem_gb=max(rule_mem_gb("pangenie", 256), 16),
        docker=config.get("pangenie_docker", "mgibio/pangenie:v4.2.1-bookworm"),
        jellyfish_size=config.get("pangenie_jellyfish_size", 3000000000),
    shell:
        "scripts/pangenie.sh"
        " --reads {input.reads}"
        " --ref {input.ref}"
        " --vcf {input.panel}"
        " --sample {wildcards.sample}"
        " --out-dir {OUT_DIR}"
        " --out-name {wildcards.sample}.pangenie.vcf.gz"
        " --docker {params.docker}"
        " --jellyfish-size {params.jellyfish_size}"
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
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
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
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} DeepVariant Off-Reference Density ({wildcards.sample})'"
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF} --offref"

rule fb_plots:
    """Per-sample FreeBayes VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/{{sample}}.freebayes.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{{sample}}.fb-offref.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} FreeBayes Off-Reference Density ({wildcards.sample})'"
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF} --offref"

rule bc_plots:
    """Per-sample bcftools VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/{{sample}}.bcftools.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{{sample}}.bc-offref.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} bcftools Off-Reference Density ({wildcards.sample})'"
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF} --offref"

rule pg_plots:
    """Per-sample PanGenie VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/{{sample}}.pangenie.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/{{sample}}.pg-offref.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} PanGenie Off-Reference Density ({wildcards.sample})'"
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

rule merge_longread_call_vcfs:
    """Merge long-read per-sample call VCFs"""
    input:
        expand("{out}/{s}.vcf.gz", out=OUT_DIR, s=LR_SAMPLES),
    output:
        f"{OUT_DIR}/merged.longread.call.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    run:
        if len(input) == 1:
            shell("bcftools +fill-tags {input} -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")
        else:
            shell("bcftools merge {input} -Oz"
                  " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")

rule merge_longread_dv_vcfs:
    """Merge long-read per-sample DeepVariant VCFs"""
    input:
        expand("{out}/{s}.deepvariant.vcf.gz", out=OUT_DIR, s=LR_SAMPLES),
    output:
        f"{OUT_DIR}/merged.longread.deepvariant.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    run:
        if len(input) == 1:
            shell("bcftools +fill-tags {input} -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")
        else:
            shell("bcftools merge {input} -Oz"
                  " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")

rule merge_call_pass_vcfs:
    """Merge PASS-filtered per-sample call VCFs (avoids cross-sample filter contamination)"""
    input:
        expand("{out}/{s}.pass-only.vcf.gz", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/merged.call.pass-prefiltered.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "bcftools merge {input} -Oz"
        " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
        " && tabix -p vcf {output}"

rule merge_longread_call_pass_vcfs:
    """Merge PASS-filtered long-read per-sample call VCFs"""
    input:
        expand("{out}/{s}.pass-only.vcf.gz", out=OUT_DIR, s=LR_SAMPLES),
    output:
        f"{OUT_DIR}/merged.longread.call.pass-prefiltered.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    run:
        if len(input) == 1:
            shell("bcftools +fill-tags {input} -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")
        else:
            shell("bcftools merge {input} -Oz"
                  " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")

rule merge_dv_pass_vcfs:
    """Merge PASS-filtered per-sample DeepVariant VCFs"""
    input:
        expand("{out}/{s}.deepvariant.pass-only.vcf.gz", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/merged.deepvariant.pass-prefiltered.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "bcftools merge {input} -Oz"
        " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
        " && tabix -p vcf {output}"

rule merge_longread_dv_pass_vcfs:
    """Merge PASS-filtered long-read per-sample DeepVariant VCFs"""
    input:
        expand("{out}/{s}.deepvariant.pass-only.vcf.gz", out=OUT_DIR, s=LR_SAMPLES),
    output:
        f"{OUT_DIR}/merged.longread.deepvariant.pass-prefiltered.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    run:
        if len(input) == 1:
            shell("bcftools +fill-tags {input} -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")
        else:
            shell("bcftools merge {input} -Oz"
                  " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")

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

rule merge_fb_vcfs:
    """Merge per-sample FreeBayes VCFs with bcftools, add AF/AC/AN tags"""
    input:
        expand("{out}/{s}.freebayes.vcf.gz", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/merged.freebayes.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "bcftools merge {input} -Oz"
        " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
        " && tabix -p vcf {output}"

rule merge_longread_fb_vcfs:
    """Merge long-read per-sample FreeBayes VCFs"""
    input:
        expand("{out}/{s}.freebayes.vcf.gz", out=OUT_DIR, s=LR_SAMPLES),
    output:
        f"{OUT_DIR}/merged.longread.freebayes.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    run:
        if len(input) == 1:
            shell("bcftools +fill-tags {input} -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")
        else:
            shell("bcftools merge {input} -Oz"
                  " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")

rule merge_fb_pass_vcfs:
    """Merge PASS-filtered per-sample FreeBayes VCFs"""
    input:
        expand("{out}/{s}.freebayes.pass-only.vcf.gz", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/merged.freebayes.pass-prefiltered.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "bcftools merge {input} -Oz"
        " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
        " && tabix -p vcf {output}"

rule merge_longread_fb_pass_vcfs:
    """Merge PASS-filtered long-read per-sample FreeBayes VCFs"""
    input:
        expand("{out}/{s}.freebayes.pass-only.vcf.gz", out=OUT_DIR, s=LR_SAMPLES),
    output:
        f"{OUT_DIR}/merged.longread.freebayes.pass-prefiltered.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    run:
        if len(input) == 1:
            shell("bcftools +fill-tags {input} -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")
        else:
            shell("bcftools merge {input} -Oz"
                  " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")

rule merge_bc_vcfs:
    """Merge per-sample bcftools VCFs with bcftools, add AF/AC/AN tags"""
    input:
        expand("{out}/{s}.bcftools.vcf.gz", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/merged.bcftools.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "bcftools merge {input} -Oz"
        " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
        " && tabix -p vcf {output}"

rule merge_longread_bc_vcfs:
    """Merge long-read per-sample bcftools VCFs"""
    input:
        expand("{out}/{s}.bcftools.vcf.gz", out=OUT_DIR, s=LR_SAMPLES),
    output:
        f"{OUT_DIR}/merged.longread.bcftools.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    run:
        if len(input) == 1:
            shell("bcftools +fill-tags {input} -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")
        else:
            shell("bcftools merge {input} -Oz"
                  " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")

rule merge_bc_pass_vcfs:
    """Merge PASS-filtered per-sample bcftools VCFs"""
    input:
        expand("{out}/{s}.bcftools.pass-only.vcf.gz", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/merged.bcftools.pass-prefiltered.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "bcftools merge {input} -Oz"
        " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
        " && tabix -p vcf {output}"

rule merge_longread_bc_pass_vcfs:
    """Merge PASS-filtered long-read per-sample bcftools VCFs"""
    input:
        expand("{out}/{s}.bcftools.pass-only.vcf.gz", out=OUT_DIR, s=LR_SAMPLES),
    output:
        f"{OUT_DIR}/merged.longread.bcftools.pass-prefiltered.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    run:
        if len(input) == 1:
            shell("bcftools +fill-tags {input} -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")
        else:
            shell("bcftools merge {input} -Oz"
                  " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
                  " && tabix -p vcf {output}")

rule merge_pg_vcfs:
    """Merge per-sample PanGenie VCFs with bcftools, add AF/AC/AN tags"""
    input:
        expand("{out}/{s}.pangenie.vcf.gz", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/merged.pangenie.vcf.gz",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "bcftools merge {input} -Oz"
        " | bcftools +fill-tags -Oz -o {output} -- -t AF,AC,AN"
        " && tabix -p vcf {output}"

rule merge_pg_pass_vcfs:
    """Merge PASS-filtered per-sample PanGenie VCFs"""
    input:
        expand("{out}/{s}.pangenie.pass-only.vcf.gz", out=OUT_DIR, s=SAMPLES),
    output:
        f"{OUT_DIR}/merged.pangenie.pass-prefiltered.vcf.gz",
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
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
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
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} Merged DeepVariant Off-Reference Density'"
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF} --offref"

rule merged_fb_plots:
    """Merged FreeBayes VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/merged.freebayes.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/merged.fb-offref.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} Merged FreeBayes Off-Reference Density'"
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF} --offref"

rule merged_bc_plots:
    """Merged bcftools VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/merged.bcftools.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/merged.bc-offref.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} Merged bcftools Off-Reference Density'"
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF} --offref"

rule merged_pg_plots:
    """Merged PanGenie VCF → density ideogram"""
    input:
        vcf=f"{OUT_DIR}/merged.pangenie.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/merged.pg-offref.png",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/chrom-density-segs.R"
        " {input.vcf} {input.segs} {output}"
        " '{REF} Merged PanGenie Off-Reference Density'"
        " 0 '{config[refgaps_bed]}' {config[scale_type]}"
        " --ref {REF} --offref"

############################################################################
# VCF statistics rules
############################################################################

rule deconstruct_sites_stats:
    """Deconstruct VCF → site-level stats + plots (includes AF spectrum)"""
    input:
        vcf=f"{OUT_DIR}/{OUT_NAME}.tr.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
        populations="sample-super-populations.tsv",
    threads: rule_cpus("deconstruct_stats", 32)
    output:
        f"{OUT_DIR}/{OUT_NAME}.sites.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.sites.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.size-dist-log.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.af-spectrum.png",
        *([ f"{OUT_DIR}/{OUT_NAME}.sites.variant-types-by-annot.png",
            f"{OUT_DIR}/{OUT_NAME}.sites.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/{OUT_NAME}.sites.giab-strat.png",
            f"{OUT_DIR}/{OUT_NAME}.sites.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/{OUT_NAME}.sites.per-sample-types.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.per-sample-types.tsv",
        f"{OUT_DIR}/{OUT_NAME}.sites.per-sample-types-by-pop.png",
        f"{OUT_DIR}/{OUT_NAME}.sites.per-sample-types-by-pop.tsv",
        f"{OUT_DIR}/{OUT_NAME}.sites.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/{OUT_NAME}.sites.per-sample-giab-strat.png",
            f"{OUT_DIR}/{OUT_NAME}.sites.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=int(rule_mem_gb("deconstruct_stats", 512)) * 1024,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{OUT_NAME}.sites"
        " --mode sites --af-step 0.05 --title '{REF} Deconstruct'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --per-sample"
        " --populations {input.populations} --ref-sample {REF}"
        " --threads {threads}"

rule deconstruct_variants_stats:
    """Deconstruct VCF → variant-level stats + plots (uses pre-normed VCF)"""
    input:
        vcf=f"{OUT_DIR}/{OUT_NAME}.tr.normed.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    threads: rule_cpus("deconstruct_stats", 32)
    output:
        f"{OUT_DIR}/{OUT_NAME}.variants.vcf-stats.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variants.variant-types.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.size-dist-log.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.af-spectrum.png",
        *([ f"{OUT_DIR}/{OUT_NAME}.variants.variant-types-by-annot.png",
            f"{OUT_DIR}/{OUT_NAME}.variants.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/{OUT_NAME}.variants.giab-strat.png",
            f"{OUT_DIR}/{OUT_NAME}.variants.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/{OUT_NAME}.variants.per-sample-types.png",
        f"{OUT_DIR}/{OUT_NAME}.variants.per-sample-types.tsv",
        f"{OUT_DIR}/{OUT_NAME}.variants.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/{OUT_NAME}.variants.per-sample-giab-strat.png",
            f"{OUT_DIR}/{OUT_NAME}.variants.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=int(rule_mem_gb("deconstruct_stats", 512)) * 1024,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{OUT_NAME}.variants"
        " --mode variants --af-step 0.05 --title '{REF} Deconstruct'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --per-sample --ref-sample {REF}"
        " --threads {threads}"

rule call_stats:
    """Per-sample call VCF → variant stats + plots (one mode/filter combo)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/{wc.sample}.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/{wc.sample}.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=call_annot_beds(),
        giab_beds=call_giab_strat_beds(),
    output:
        f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.size-dist-log.png",
        *([ f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/{{sample}}.call.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{wildcards.sample}.call.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Call ({wildcards.sample})'"
        " --segs {input.segs} --segs-strip-prefix '{AUGREF}#0#'"
        " {params.annot_arg} {params.giab_arg}"

rule dv_stats:
    """Per-sample DeepVariant VCF → variant stats + plots (one mode/filter combo)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/{wc.sample}.deepvariant.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/{wc.sample}.deepvariant.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    output:
        f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.size-dist-log.png",
        *([ f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/{{sample}}.dv.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{wildcards.sample}.dv.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} DeepVariant ({wildcards.sample})'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --no-sv"

rule fb_stats:
    """Per-sample FreeBayes VCF → variant stats + plots (one mode/filter combo)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/{wc.sample}.freebayes.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/{wc.sample}.freebayes.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    output:
        f"{OUT_DIR}/{{sample}}.fb.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.fb.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/{{sample}}.fb.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/{{sample}}.fb.{{mode}}.{{filt}}.size-dist-log.png",
        *([ f"{OUT_DIR}/{{sample}}.fb.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/{{sample}}.fb.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/{{sample}}.fb.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/{{sample}}.fb.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{wildcards.sample}.fb.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} FreeBayes ({wildcards.sample})'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --no-sv"

rule bc_stats:
    """Per-sample bcftools VCF → variant stats + plots (one mode/filter combo)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/{wc.sample}.bcftools.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/{wc.sample}.bcftools.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    output:
        f"{OUT_DIR}/{{sample}}.bc.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.bc.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/{{sample}}.bc.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/{{sample}}.bc.{{mode}}.{{filt}}.size-dist-log.png",
        *([ f"{OUT_DIR}/{{sample}}.bc.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/{{sample}}.bc.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/{{sample}}.bc.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/{{sample}}.bc.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{wildcards.sample}.bc.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} bcftools ({wildcards.sample})'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --no-sv"

rule pg_stats:
    """Per-sample PanGenie VCF → variant stats + plots (one mode/filter combo)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/{wc.sample}.pangenie.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/{wc.sample}.pangenie.vcf.gz",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    output:
        f"{OUT_DIR}/{{sample}}.pg.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/{{sample}}.pg.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/{{sample}}.pg.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/{{sample}}.pg.{{mode}}.{{filt}}.size-dist-log.png",
        *([ f"{OUT_DIR}/{{sample}}.pg.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/{{sample}}.pg.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/{{sample}}.pg.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/{{sample}}.pg.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/{wildcards.sample}.pg.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} PanGenie ({wildcards.sample})'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --no-sv"

rule merged_call_stats:
    """Merged call VCF → variant stats + plots (one mode/filter combo, includes AF spectrum)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.call.pass-prefiltered{'.normed' if wc.mode == 'variants' else ''}.vcf.gz" if wc.filt == "pass" else (f"{OUT_DIR}/merged.call.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.call.vcf.gz"),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=call_annot_beds(),
        giab_beds=call_giab_strat_beds(),
    threads: rule_cpus("merged_stats", 32)
    output:
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv",
            f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.annot-exclusive.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.per-sample-types.png",
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.per-sample-types.tsv",
        f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.per-sample-giab-strat.png",
            f"{OUT_DIR}/merged.call.{{mode}}.{{filt}}.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.call.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged Call'"
        " --segs {input.segs} --segs-strip-prefix '{AUGREF}#0#'"
        " {params.annot_arg} {params.giab_arg} --per-sample"
        " --threads {threads}"

rule merged_longread_call_stats:
    """Merged long-read call VCF → variant stats + plots"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.longread.call.pass-prefiltered{'.normed' if wc.mode == 'variants' else ''}.vcf.gz" if wc.filt == "pass" else (f"{OUT_DIR}/merged.longread.call.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.longread.call.vcf.gz"),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=call_annot_beds(),
        giab_beds=call_giab_strat_beds(),
    threads: rule_cpus("merged_stats", 32)
    output:
        f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv",
            f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.annot-exclusive.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.per-sample-types.png",
        f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.per-sample-types.tsv",
        f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.per-sample-giab-strat.png",
            f"{OUT_DIR}/merged.longread.call.{{mode}}.{{filt}}.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.longread.call.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged Long-Read Call'"
        " --segs {input.segs} --segs-strip-prefix '{AUGREF}#0#'"
        " {params.annot_arg} {params.giab_arg} --per-sample"
        " --threads {threads}"

rule merged_longread_dv_stats:
    """Merged long-read DeepVariant VCF → variant stats + plots"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.longread.deepvariant.pass-prefiltered{'.normed' if wc.mode == 'variants' else ''}.vcf.gz" if wc.filt == "pass" else (f"{OUT_DIR}/merged.longread.deepvariant.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.longread.deepvariant.vcf.gz"),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    threads: rule_cpus("merged_stats", 32)
    output:
        f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.per-sample-types.png",
        f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.per-sample-types.tsv",
        f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.per-sample-giab-strat.png",
            f"{OUT_DIR}/merged.longread.dv.{{mode}}.{{filt}}.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.longread.dv.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged Long-Read DeepVariant'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --per-sample --no-sv"
        " --threads {threads}"

rule merged_dv_stats:
    """Merged DeepVariant VCF → variant stats + plots (one mode/filter combo, includes AF spectrum)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.deepvariant.pass-prefiltered{'.normed' if wc.mode == 'variants' else ''}.vcf.gz" if wc.filt == "pass" else (f"{OUT_DIR}/merged.deepvariant.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.deepvariant.vcf.gz"),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    threads: rule_cpus("merged_stats", 32)
    output:
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.per-sample-types.png",
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.per-sample-types.tsv",
        f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.per-sample-giab-strat.png",
            f"{OUT_DIR}/merged.dv.{{mode}}.{{filt}}.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.dv.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged DeepVariant'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --per-sample --no-sv"
        " --threads {threads}"

rule merged_fb_stats:
    """Merged FreeBayes VCF → variant stats + plots (one mode/filter combo, includes AF spectrum)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.freebayes.pass-prefiltered{'.normed' if wc.mode == 'variants' else ''}.vcf.gz" if wc.filt == "pass" else (f"{OUT_DIR}/merged.freebayes.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.freebayes.vcf.gz"),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    threads: rule_cpus("merged_stats", 32)
    output:
        f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.per-sample-types.png",
        f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.per-sample-types.tsv",
        f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.per-sample-giab-strat.png",
            f"{OUT_DIR}/merged.fb.{{mode}}.{{filt}}.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.fb.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged FreeBayes'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --per-sample --no-sv"
        " --threads {threads}"

rule merged_longread_fb_stats:
    """Merged long-read FreeBayes VCF → variant stats + plots"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.longread.freebayes.pass-prefiltered{'.normed' if wc.mode == 'variants' else ''}.vcf.gz" if wc.filt == "pass" else (f"{OUT_DIR}/merged.longread.freebayes.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.longread.freebayes.vcf.gz"),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    threads: rule_cpus("merged_stats", 32)
    output:
        f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.per-sample-types.png",
        f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.per-sample-types.tsv",
        f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.per-sample-giab-strat.png",
            f"{OUT_DIR}/merged.longread.fb.{{mode}}.{{filt}}.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.longread.fb.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged Long-Read FreeBayes'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --per-sample --no-sv"
        " --threads {threads}"

rule merged_bc_stats:
    """Merged bcftools VCF → variant stats + plots (one mode/filter combo, includes AF spectrum)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.bcftools.pass-prefiltered{'.normed' if wc.mode == 'variants' else ''}.vcf.gz" if wc.filt == "pass" else (f"{OUT_DIR}/merged.bcftools.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.bcftools.vcf.gz"),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    threads: rule_cpus("merged_stats", 32)
    output:
        f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.per-sample-types.png",
        f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.per-sample-types.tsv",
        f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.per-sample-giab-strat.png",
            f"{OUT_DIR}/merged.bc.{{mode}}.{{filt}}.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.bc.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged bcftools'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --per-sample --no-sv"
        " --threads {threads}"

rule merged_longread_bc_stats:
    """Merged long-read bcftools VCF → variant stats + plots"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.longread.bcftools.pass-prefiltered{'.normed' if wc.mode == 'variants' else ''}.vcf.gz" if wc.filt == "pass" else (f"{OUT_DIR}/merged.longread.bcftools.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.longread.bcftools.vcf.gz"),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    threads: rule_cpus("merged_stats", 32)
    output:
        f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.per-sample-types.png",
        f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.per-sample-types.tsv",
        f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.per-sample-giab-strat.png",
            f"{OUT_DIR}/merged.longread.bc.{{mode}}.{{filt}}.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.longread.bc.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged Long-Read bcftools'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --per-sample --no-sv"
        " --threads {threads}"

rule merged_pg_stats:
    """Merged PanGenie VCF → variant stats + plots (one mode/filter combo, includes AF spectrum)"""
    input:
        vcf=lambda wc: f"{OUT_DIR}/merged.pangenie.pass-prefiltered{'.normed' if wc.mode == 'variants' else ''}.vcf.gz" if wc.filt == "pass" else (f"{OUT_DIR}/merged.pangenie.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.pangenie.vcf.gz"),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot_beds=augref_annot_beds(),
        giab_beds=augref_giab_strat_beds(),
    threads: rule_cpus("merged_stats", 32)
    output:
        f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.vcf-stats.tsv",
        f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.variant-types.png",
        f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.size-dist.png",
        f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.size-dist-log.png",
        f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.af-spectrum.png",
        *([ f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.variant-types-by-annot.png",
            f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.vcf-stats-by-annot.tsv"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.giab-strat.png",
            f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.giab-strat.tsv"]
          if giab_strat_configured() else []),
        f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.per-sample-types.png",
        f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.per-sample-types.tsv",
        f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.per-sample-sv-types.png",
        *([ f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.per-sample-giab-strat.png",
            f"{OUT_DIR}/merged.pg.{{mode}}.{{filt}}.per-sample-giab-strat.tsv"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=256000,
        runtime=2880,
    params:
        annot_arg=lambda wc, input: (
            f"--annot-beds {','.join(input.annot_beds)} --annot-names {','.join(annotation_names())}"
            if annotation_inputs() else ""),
        giab_arg=lambda wc, input: (
            f"--giab-strat-beds {','.join(input.giab_beds)}"
            f" --giab-strat-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input.vcf} {OUT_DIR}/merged.pg.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt} --title '{REF} Merged PanGenie'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg} --per-sample --no-sv"
        " --threads {threads}"

############################################################################
# Call vs DeepVariant comparison
############################################################################

rule compare_call_dv:
    """Compare merged call and merged DeepVariant VCFs per-sample"""
    input:
        call_vcf=lambda wc: f"{OUT_DIR}/merged.call.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.call.vcf.gz",
        dv_vcf=lambda wc: f"{OUT_DIR}/merged.deepvariant.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.deepvariant.vcf.gz",
    output:
        f"{OUT_DIR}/merged.call-vs-dv.{{mode}}.{{filt}}.compare.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{mode}}.{{filt}}.compare.tsv",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare.R {input.call_vcf} {input.dv_vcf}"
        " {OUT_DIR}/merged.call-vs-dv.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt}"
        " --label-a Call --label-b DeepVariant"
        " --title '{REF} Call vs DeepVariant'"
        " --strip-prefix '{AUGREF}#0#'"
        " --no-sv"

############################################################################
# VCF comparison: Call vs DeepVariant (vcfeval or aardvark)
############################################################################

rule vcfeval_per_sample:
    """Run VCF comparison per sample: call VCF (truth) vs DeepVariant VCF (calls)

    Dispatches to vcfeval or aardvark based on config['eval_tool'].
    When filt=pass, both VCFs are pre-filtered to PASS before comparison
    (aardvark/vcfeval strip FILTER, so post-hoc filtering doesn't work).

    vg call emits plain locus names (e.g. 'chr1') as CHROM while DeepVariant
    uses the full augref path (e.g. 'augref_CHM13#0#chr1').  The rename step
    restores the prefix on the call VCF so all inputs share the same namespace.
    It is idempotent: only renames CHROMs that lack the prefix.
    """
    input:
        call_vcf=f"{OUT_DIR}/{{sample}}.filtered.vcf.gz" if surject_filtering() else f"{OUT_DIR}/{{sample}}.vcf.gz",
        dv_vcf=f"{OUT_DIR}/{{sample}}.deepvariant.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        tp_baseline=f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fn.vcf.gz",
    threads: rule_cpus("vcfeval", 64)
    resources:
        mem_mb=rule_mem_gb("vcfeval", 128) * 1024,
        runtime=rule_runtime("vcfeval"),
    params:
        out_dir=f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}",
        eval_tool=config.get("eval_tool", "aardvark"),
        docker_arg=lambda wc: f"--docker {config['vcfeval_docker']}" if config.get("vcfeval_docker") else "",
        no_docker="" if config.get("vcfeval_docker") else "--no-docker",
        augref_prefix=f"{AUGREF}#0#",
        min_vcfeval_len=config.get("min_vcfeval_len", 0),
    shell:
        # Build contig rename map: only rename CHROMs missing the augref prefix
        # (idempotent: works with both old vg [plain names] and fixed vg [augref names])
        "export RTG_MEM=$(({resources.mem_mb} / 1024))g"
        " && mkdir -p {params.out_dir}"
        " && bcftools query -f '%CHROM\\n' {input.call_vcf} | sort -u"
        "    | sed -n '/^{params.augref_prefix}/!s/^\\(.*\\)/\\1\\t{params.augref_prefix}\\1/p'"
        "    > {params.out_dir}/rename-chrs.txt"
        # Rename chroms (no-op if rename file is empty); optionally pre-filter to PASS
        " && if [ -s {params.out_dir}/rename-chrs.txt ]; then"
        "      bcftools annotate --rename-chrs {params.out_dir}/rename-chrs.txt"
        "        {input.call_vcf}"
        "        | awk '/^##contig=/{{id=$0; sub(/.*ID=/, \"\", id); sub(/[,>].*/, \"\", id);"
        "                if(seen[id]++) next}} {{print}}';"
        "    else"
        "      bcftools view {input.call_vcf};"
        "    fi"
        "    | if [ '{wildcards.filt}' = 'pass' ]; then"
        "        bcftools view -f PASS 2>/dev/null;"
        "      else cat; fi"
        "    | bgzip > {params.out_dir}/call.renamed.vcf.gz"
        " && tabix -fp vcf {params.out_dir}/call.renamed.vcf.gz"
        # Pre-filter DV VCF to PASS if needed
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      bcftools view -f PASS {input.dv_vcf} 2>/dev/null"
        "        | bgzip > {params.out_dir}/dv.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/dv.pass.vcf.gz;"
        "    fi"
        # Stage inputs to node-local scratch for fast I/O
        " && WORK_TMPDIR=$(mktemp -d \"${{TMPDIR:-{params.out_dir}}}/vcfeval.XXXXXX\")"
        " && trap 'rm -rf \"$WORK_TMPDIR\"' EXIT"
        " && echo \"Staging inputs to $WORK_TMPDIR\""
        " && cp {params.out_dir}/call.renamed.vcf.gz"
        "       {params.out_dir}/call.renamed.vcf.gz.tbi"
        "       {input.ref} {input.ref}.fai"
        "       \"$WORK_TMPDIR/\""
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      cp {params.out_dir}/dv.pass.vcf.gz"
        "         {params.out_dir}/dv.pass.vcf.gz.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      DV_VCF=$WORK_TMPDIR/dv.pass.vcf.gz;"
        "    else"
        "      cp {input.dv_vcf} {input.dv_vcf}.tbi \"$WORK_TMPDIR/\";"
        "      DV_VCF=$WORK_TMPDIR/$(basename {input.dv_vcf});"
        "    fi"
        " && {{ [ -f {input.ref}.gzi ]"
        "       && cp {input.ref}.gzi \"$WORK_TMPDIR/\" || true; }}"
        " && python3 scripts/vcfcomp.py {params.eval_tool}"
        "    --truth $WORK_TMPDIR/call.renamed.vcf.gz"
        "    --calls $DV_VCF"
        "    --ref $WORK_TMPDIR/$(basename {input.ref})"
        "    --out-dir $WORK_TMPDIR"
        "    --threads {threads}"
        "    --min-contig-len {params.min_vcfeval_len}"
        "    {params.docker_arg} {params.no_docker}"
        # Copy results back from local scratch
        " && for f in tp-baseline.vcf.gz tp-baseline.vcf.gz.tbi"
        "          fp.vcf.gz fp.vcf.gz.tbi fn.vcf.gz fn.vcf.gz.tbi"
        "          summary.txt snp_roc.tsv.gz non_snp_roc.tsv.gz weighted_roc.tsv.gz"
        "          phasing.txt vcfeval.log progress"
        "          query.vcf.gz query.vcf.gz.tbi truth.vcf.gz truth.vcf.gz.tbi; do"
        "    [ -f \"$WORK_TMPDIR/$f\" ] && cp \"$WORK_TMPDIR/$f\" {params.out_dir}/;"
        "  done"
        # Clean up staged files on shared storage
        " && rm -f {params.out_dir}/call.renamed.vcf.gz"
        "    {params.out_dir}/call.renamed.vcf.gz.tbi"
        "    {params.out_dir}/rename-chrs.txt"
        "    {params.out_dir}/dv.pass.vcf.gz"
        "    {params.out_dir}/dv.pass.vcf.gz.tbi"

rule vcfeval_compare_plot:
    """Aggregate per-sample vcfeval results into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval/{wc.filt}/{s}" for s in SAMPLES]),
        sample_names=",".join(SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.call-vs-dv.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a Call --label-b DeepVariant"
        " --title '{REF} Call vs DeepVariant (vcfeval)'"
        " --no-sv"

rule vcfeval_per_sample_squash:
    """Run VCF comparison per sample with genotypes squashed (het→hom).

    For vcfeval: uses --squash-ploidy.
    For aardvark: preprocesses VCFs to set all non-ref GTs to 1/1.
    When filt=pass, both VCFs are pre-filtered to PASS before comparison.
    """
    input:
        call_vcf=f"{OUT_DIR}/{{sample}}.filtered.vcf.gz" if surject_filtering() else f"{OUT_DIR}/{{sample}}.vcf.gz",
        dv_vcf=f"{OUT_DIR}/{{sample}}.deepvariant.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        tp_baseline=f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/fn.vcf.gz",
    threads: rule_cpus("vcfeval", 64)
    resources:
        mem_mb=rule_mem_gb("vcfeval", 128) * 1024,
        runtime=rule_runtime("vcfeval"),
    params:
        out_dir=f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}",
        eval_tool=config.get("eval_tool", "aardvark"),
        docker_arg=lambda wc: f"--docker {config['vcfeval_docker']}" if config.get("vcfeval_docker") else "",
        no_docker="" if config.get("vcfeval_docker") else "--no-docker",
        augref_prefix=f"{AUGREF}#0#",
        min_vcfeval_len=config.get("min_vcfeval_len", 0),
    shell:
        # Build contig rename map: only rename CHROMs missing the augref prefix
        # (idempotent: works with both old vg [plain names] and fixed vg [augref names])
        "export RTG_MEM=$(({resources.mem_mb} / 1024))g"
        " && mkdir -p {params.out_dir}"
        " && bcftools query -f '%CHROM\\n' {input.call_vcf} | sort -u"
        "    | sed -n '/^{params.augref_prefix}/!s/^\\(.*\\)/\\1\\t{params.augref_prefix}\\1/p'"
        "    > {params.out_dir}/rename-chrs.txt"
        # Rename chroms (no-op if rename file is empty); optionally pre-filter to PASS
        " && if [ -s {params.out_dir}/rename-chrs.txt ]; then"
        "      bcftools annotate --rename-chrs {params.out_dir}/rename-chrs.txt"
        "        {input.call_vcf}"
        "        | awk '/^##contig=/{{id=$0; sub(/.*ID=/, \"\", id); sub(/[,>].*/, \"\", id);"
        "                if(seen[id]++) next}} {{print}}';"
        "    else"
        "      bcftools view {input.call_vcf};"
        "    fi"
        "    | if [ '{wildcards.filt}' = 'pass' ]; then"
        "        bcftools view -f PASS 2>/dev/null;"
        "      else cat; fi"
        "    | bgzip > {params.out_dir}/call.renamed.vcf.gz"
        " && tabix -fp vcf {params.out_dir}/call.renamed.vcf.gz"
        # Pre-filter DV VCF to PASS if needed
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      bcftools view -f PASS {input.dv_vcf} 2>/dev/null"
        "        | bgzip > {params.out_dir}/dv.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/dv.pass.vcf.gz;"
        "    fi"
        # Stage inputs to node-local scratch for fast I/O
        " && WORK_TMPDIR=$(mktemp -d \"${{TMPDIR:-{params.out_dir}}}/vcfeval.XXXXXX\")"
        " && trap 'rm -rf \"$WORK_TMPDIR\"' EXIT"
        " && echo \"Staging inputs to $WORK_TMPDIR\""
        " && cp {params.out_dir}/call.renamed.vcf.gz"
        "       {params.out_dir}/call.renamed.vcf.gz.tbi"
        "       {input.ref} {input.ref}.fai"
        "       \"$WORK_TMPDIR/\""
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      cp {params.out_dir}/dv.pass.vcf.gz"
        "         {params.out_dir}/dv.pass.vcf.gz.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      DV_VCF=$WORK_TMPDIR/dv.pass.vcf.gz;"
        "    else"
        "      cp {input.dv_vcf} {input.dv_vcf}.tbi \"$WORK_TMPDIR/\";"
        "      DV_VCF=$WORK_TMPDIR/$(basename {input.dv_vcf});"
        "    fi"
        " && {{ [ -f {input.ref}.gzi ]"
        "       && cp {input.ref}.gzi \"$WORK_TMPDIR/\" || true; }}"
        # Squash genotypes: for vcfeval use --squash-ploidy; for aardvark preprocess VCFs
        " && if [ '{params.eval_tool}' = 'vcfeval' ]; then"
        "      python3 scripts/vcfcomp.py {params.eval_tool}"
        "        --truth $WORK_TMPDIR/call.renamed.vcf.gz"
        "        --calls $DV_VCF"
        "        --ref $WORK_TMPDIR/$(basename {input.ref})"
        "        --out-dir $WORK_TMPDIR"
        "        --threads {threads}"
        "        --min-contig-len {params.min_vcfeval_len}"
        "        --options '--decompose --ref-overlap --squash-ploidy'"
        "        {params.docker_arg} {params.no_docker};"
        "    else"
        "      bcftools +setGT $WORK_TMPDIR/call.renamed.vcf.gz"
        "        -- -t q -n c:'1/1' -i 'GT=\"het\"'"
        "        | bgzip > $WORK_TMPDIR/call.squash.vcf.gz"
        "      && tabix -fp vcf $WORK_TMPDIR/call.squash.vcf.gz"
        "      && bcftools +setGT $DV_VCF"
        "        -- -t q -n c:'1/1' -i 'GT=\"het\"'"
        "        | bgzip > $WORK_TMPDIR/dv.squash.vcf.gz"
        "      && tabix -fp vcf $WORK_TMPDIR/dv.squash.vcf.gz"
        "      && python3 scripts/vcfcomp.py {params.eval_tool}"
        "        --truth $WORK_TMPDIR/call.squash.vcf.gz"
        "        --calls $WORK_TMPDIR/dv.squash.vcf.gz"
        "        --ref $WORK_TMPDIR/$(basename {input.ref})"
        "        --out-dir $WORK_TMPDIR"
        "        --threads {threads}"
        "        --min-contig-len {params.min_vcfeval_len}"
        "        {params.docker_arg} {params.no_docker};"
        "    fi"
        # Copy results back from local scratch
        " && for f in tp-baseline.vcf.gz tp-baseline.vcf.gz.tbi"
        "          fp.vcf.gz fp.vcf.gz.tbi fn.vcf.gz fn.vcf.gz.tbi"
        "          summary.txt snp_roc.tsv.gz non_snp_roc.tsv.gz weighted_roc.tsv.gz"
        "          phasing.txt vcfeval.log progress"
        "          query.vcf.gz query.vcf.gz.tbi truth.vcf.gz truth.vcf.gz.tbi; do"
        "    [ -f \"$WORK_TMPDIR/$f\" ] && cp \"$WORK_TMPDIR/$f\" {params.out_dir}/;"
        "  done"
        " && rm -f {params.out_dir}/call.renamed.vcf.gz"
        "    {params.out_dir}/call.renamed.vcf.gz.tbi"
        "    {params.out_dir}/rename-chrs.txt"
        "    {params.out_dir}/dv.pass.vcf.gz"
        "    {params.out_dir}/dv.pass.vcf.gz.tbi"

rule vcfeval_compare_plot_squash:
    """Aggregate per-sample squashed vcfeval results into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.vcfeval-compare.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-squash/{wc.filt}/{s}" for s in SAMPLES]),
        sample_names=",".join(SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.call-vs-dv.{wildcards.filt}.vcfeval-squash"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a Call --label-b DeepVariant"
        " --title '{REF} Call vs DeepVariant (vcfeval, squash-ploidy)'"
        " --no-sv"

rule vcfeval_chromsplit:
    """Per-contig FP/FN/TP breakdown from vcfeval/aardvark output"""
    input:
        fp=f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fn.vcf.gz",
        tp=f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
    output:
        f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        out_dir=f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}",
        subcommand="aardvark-breakdown" if config.get("eval_tool", "aardvark") == "aardvark" else "vcfeval-breakdown",
    shell:
        "python3 scripts/vcfcomp.py {params.subcommand}"
        " --dir {params.out_dir} > {output}"

rule vcfeval_chromsplit_squash:
    """Per-contig FP/FN/TP breakdown from squashed vcfeval/aardvark output"""
    input:
        fp=f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/fn.vcf.gz",
        tp=f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
    output:
        f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        out_dir=f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}",
        subcommand="aardvark-breakdown" if config.get("eval_tool", "aardvark") == "aardvark" else "vcfeval-breakdown",
    shell:
        "python3 scripts/vcfcomp.py {params.subcommand}"
        " --dir {params.out_dir} > {output}"

rule vcfeval_chromsplit_merge:
    """Merge per-sample chromsplit breakdowns into a single long-format TSV"""
    input:
        expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/chromsplit.tsv",
               sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=30,
    run:
        merge_chromsplit_tsv(input, SAMPLES, output[0])

rule vcfeval_chromsplit_squash_merge:
    """Merge per-sample squashed chromsplit breakdowns into a single long-format TSV"""
    input:
        expand(f"{OUT_DIR}/vcfeval-squash/{{filt}}/{{sample}}/chromsplit.tsv",
               sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=30,
    run:
        merge_chromsplit_tsv(input, SAMPLES, output[0])

rule vcfeval_chromsplit_plot:
    """Per-contig FP/FN scatter and top-discordant bar chart"""
    input:
        tsv=f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.chromsplit.tsv",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv" if annotation_inputs() else [],
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.chromsplit.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.chromsplit-top.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.chromsplit-concordant.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.chromsplit-top-onref.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.chromsplit-concordant-onref.png",
        *([ f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.chromsplit-annot.png"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.chromsplit-giab.png"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        strip_prefix=f"{AUGREF}#0#",
        annot_arg=lambda wc, input: f"--annot {input.annot}" if annotation_inputs() else "",
        giab_arg=lambda wc, input: (
            f"--giab-beds {','.join(input.giab_beds)} --giab-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-chromsplit-plot.R {input.tsv}"
        " {OUT_DIR}/merged.call-vs-dv.{wildcards.filt}"
        " --title '{REF} Call vs DeepVariant Per-Contig'"
        " --strip-prefix '{params.strip_prefix}'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg}"

rule vcfeval_chromsplit_squash_plot:
    """Per-contig FP/FN/TP scatter and top bar charts (squash-ploidy)"""
    input:
        tsv=f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.chromsplit.tsv",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv" if annotation_inputs() else [],
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.chromsplit.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.chromsplit-top.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.chromsplit-concordant.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.chromsplit-top-onref.png",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.chromsplit-concordant-onref.png",
        *([ f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.chromsplit-annot.png"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.vcfeval-squash.chromsplit-giab.png"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        strip_prefix=f"{AUGREF}#0#",
        annot_arg=lambda wc, input: f"--annot {input.annot}" if annotation_inputs() else "",
        giab_arg=lambda wc, input: (
            f"--giab-beds {','.join(input.giab_beds)} --giab-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-chromsplit-plot.R {input.tsv}"
        " {OUT_DIR}/merged.call-vs-dv.{wildcards.filt}.vcfeval-squash"
        " --title '{REF} Call vs DeepVariant Per-Contig (squash-ploidy)'"
        " --strip-prefix '{params.strip_prefix}'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg}"

############################################################################
# Long-read vcfeval comparison rules (call vs DV, no squash-ploidy)
############################################################################

rule vcfeval_lr_compare_plot:
    """Aggregate long-read per-sample vcfeval results into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fp.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fn.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval/{wc.filt}/{s}" for s in LR_SAMPLES]),
        sample_names=",".join(LR_SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.lr.call-vs-dv.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a Call --label-b DeepVariant"
        " --title '{REF} Long-Read Call vs DeepVariant (vcfeval)'"
        " --no-sv"

rule vcfeval_lr_chromsplit_merge:
    """Merge long-read per-sample chromsplit breakdowns into a single long-format TSV"""
    input:
        expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/chromsplit.tsv",
               sample=LR_SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=30,
    run:
        merge_chromsplit_tsv(input, LR_SAMPLES, output[0])

rule vcfeval_lr_chromsplit_plot:
    """Long-read per-contig FP/FN scatter and top-discordant bar chart"""
    input:
        tsv=f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.chromsplit.tsv",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv" if annotation_inputs() else [],
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.chromsplit.png",
        f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.chromsplit-top.png",
        f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.chromsplit-concordant.png",
        f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.chromsplit-top-onref.png",
        f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.chromsplit-concordant-onref.png",
        *([ f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.chromsplit-annot.png"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.chromsplit-giab.png"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        strip_prefix=f"{AUGREF}#0#",
        annot_arg=lambda wc, input: f"--annot {input.annot}" if annotation_inputs() else "",
        giab_arg=lambda wc, input: (
            f"--giab-beds {','.join(input.giab_beds)} --giab-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-chromsplit-plot.R {input.tsv}"
        " {OUT_DIR}/merged.lr.call-vs-dv.{wildcards.filt}"
        " --title '{REF} Long-Read Call vs DeepVariant Per-Contig'"
        " --strip-prefix '{params.strip_prefix}'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg}"

############################################################################
# Call vs FreeBayes comparison
############################################################################

rule compare_call_fb:
    """Compare merged call and merged FreeBayes VCFs per-sample"""
    input:
        call_vcf=lambda wc: f"{OUT_DIR}/merged.call.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.call.vcf.gz",
        fb_vcf=lambda wc: f"{OUT_DIR}/merged.freebayes.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.freebayes.vcf.gz",
    output:
        f"{OUT_DIR}/merged.call-vs-fb.{{mode}}.{{filt}}.compare.png",
        f"{OUT_DIR}/merged.call-vs-fb.{{mode}}.{{filt}}.compare.tsv",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare.R {input.call_vcf} {input.fb_vcf}"
        " {OUT_DIR}/merged.call-vs-fb.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt}"
        " --label-a Call --label-b FreeBayes"
        " --title '{REF} Call vs FreeBayes'"
        " --strip-prefix '{AUGREF}#0#'"
        " --no-sv"

rule vcfeval_fb_per_sample:
    """Run VCF comparison per sample: call VCF (truth) vs FreeBayes VCF (calls)

    Dispatches to vcfeval or aardvark based on config['eval_tool'].
    When filt=pass, both VCFs are pre-filtered to PASS before comparison
    (aardvark/vcfeval strip FILTER, so post-hoc filtering doesn't work).

    vg call emits plain locus names (e.g. 'chr1') as CHROM while FreeBayes
    uses the full augref path (e.g. 'augref_CHM13#0#chr1').  The rename step
    restores the prefix on the call VCF so all inputs share the same namespace.
    It is idempotent: only renames CHROMs that lack the prefix.
    """
    input:
        call_vcf=f"{OUT_DIR}/{{sample}}.filtered.vcf.gz" if surject_filtering() else f"{OUT_DIR}/{{sample}}.vcf.gz",
        fb_vcf=f"{OUT_DIR}/{{sample}}.freebayes.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        tp_baseline=f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fn.vcf.gz",
    threads: rule_cpus("vcfeval", 64)
    resources:
        mem_mb=rule_mem_gb("vcfeval", 128) * 1024,
        runtime=rule_runtime("vcfeval"),
    params:
        out_dir=f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}",
        eval_tool=config.get("eval_tool", "aardvark"),
        docker_arg=lambda wc: f"--docker {config['vcfeval_docker']}" if config.get("vcfeval_docker") else "",
        no_docker="" if config.get("vcfeval_docker") else "--no-docker",
        augref_prefix=f"{AUGREF}#0#",
        min_vcfeval_len=config.get("min_vcfeval_len", 0),
        fb_min_qual=config.get("freebayes_min_qual", 20),
    shell:
        # Build contig rename map: only rename CHROMs missing the augref prefix
        # (idempotent: works with both old vg [plain names] and fixed vg [augref names])
        "export RTG_MEM=$(({resources.mem_mb} / 1024))g"
        " && mkdir -p {params.out_dir}"
        " && bcftools query -f '%CHROM\\n' {input.call_vcf} | sort -u"
        "    | sed -n '/^{params.augref_prefix}/!s/^\\(.*\\)/\\1\\t{params.augref_prefix}\\1/p'"
        "    > {params.out_dir}/rename-chrs.txt"
        # Rename chroms (no-op if rename file is empty); optionally pre-filter to PASS
        " && if [ -s {params.out_dir}/rename-chrs.txt ]; then"
        "      bcftools annotate --rename-chrs {params.out_dir}/rename-chrs.txt"
        "        {input.call_vcf}"
        "        | awk '/^##contig=/{{id=$0; sub(/.*ID=/, \"\", id); sub(/[,>].*/, \"\", id);"
        "                if(seen[id]++) next}} {{print}}';"
        "    else"
        "      bcftools view {input.call_vcf};"
        "    fi"
        "    | if [ '{wildcards.filt}' = 'pass' ]; then"
        "        bcftools view -f PASS 2>/dev/null;"
        "      else cat; fi"
        "    | bgzip > {params.out_dir}/call.renamed.vcf.gz"
        " && tabix -fp vcf {params.out_dir}/call.renamed.vcf.gz"
        # Pre-filter FB VCF: apply QUAL threshold then PASS filter
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      bcftools filter -e 'QUAL<{params.fb_min_qual}' -s LowQual {input.fb_vcf} 2>/dev/null"
        "        | bcftools view -f PASS 2>/dev/null"
        "        | bgzip > {params.out_dir}/fb.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/fb.pass.vcf.gz;"
        "    fi"
        # Stage inputs to node-local scratch for fast I/O
        " && WORK_TMPDIR=$(mktemp -d \"${{TMPDIR:-{params.out_dir}}}/vcfeval.XXXXXX\")"
        " && trap 'rm -rf \"$WORK_TMPDIR\"' EXIT"
        " && echo \"Staging inputs to $WORK_TMPDIR\""
        " && cp {params.out_dir}/call.renamed.vcf.gz"
        "       {params.out_dir}/call.renamed.vcf.gz.tbi"
        "       {input.ref} {input.ref}.fai"
        "       \"$WORK_TMPDIR/\""
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      cp {params.out_dir}/fb.pass.vcf.gz"
        "         {params.out_dir}/fb.pass.vcf.gz.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      FB_VCF=$WORK_TMPDIR/fb.pass.vcf.gz;"
        "    else"
        "      cp {input.fb_vcf} {input.fb_vcf}.tbi \"$WORK_TMPDIR/\";"
        "      FB_VCF=$WORK_TMPDIR/$(basename {input.fb_vcf});"
        "    fi"
        " && {{ [ -f {input.ref}.gzi ]"
        "       && cp {input.ref}.gzi \"$WORK_TMPDIR/\" || true; }}"
        " && python3 scripts/vcfcomp.py {params.eval_tool}"
        "    --truth $WORK_TMPDIR/call.renamed.vcf.gz"
        "    --calls $FB_VCF"
        "    --ref $WORK_TMPDIR/$(basename {input.ref})"
        "    --out-dir $WORK_TMPDIR"
        "    --threads {threads}"
        "    --min-contig-len {params.min_vcfeval_len}"
        "    {params.docker_arg} {params.no_docker}"
        # Copy results back from local scratch
        " && for f in tp-baseline.vcf.gz tp-baseline.vcf.gz.tbi"
        "          fp.vcf.gz fp.vcf.gz.tbi fn.vcf.gz fn.vcf.gz.tbi"
        "          summary.txt snp_roc.tsv.gz non_snp_roc.tsv.gz weighted_roc.tsv.gz"
        "          phasing.txt vcfeval.log progress"
        "          query.vcf.gz query.vcf.gz.tbi truth.vcf.gz truth.vcf.gz.tbi; do"
        "    [ -f \"$WORK_TMPDIR/$f\" ] && cp \"$WORK_TMPDIR/$f\" {params.out_dir}/;"
        "  done"
        # Clean up staged files on shared storage
        " && rm -f {params.out_dir}/call.renamed.vcf.gz"
        "    {params.out_dir}/call.renamed.vcf.gz.tbi"
        "    {params.out_dir}/rename-chrs.txt"
        "    {params.out_dir}/fb.pass.vcf.gz"
        "    {params.out_dir}/fb.pass.vcf.gz.tbi"

rule vcfeval_fb_relaxed_per_sample:
    """Relaxed call-vs-FB comparison: unfiltered call (truth) vs QUAL-filtered FB (calls).

    TP = FB variant confirmed by any call record (including filtered).
    FP = FB variant absent from call entirely.
    FN = PASS call variant not in FB.
    """
    input:
        call_vcf=f"{OUT_DIR}/{{sample}}.filtered.vcf.gz" if surject_filtering() else f"{OUT_DIR}/{{sample}}.vcf.gz",
        fb_vcf=f"{OUT_DIR}/{{sample}}.freebayes.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        tp_baseline=f"{OUT_DIR}/vcfeval-fb-relaxed/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval-fb-relaxed/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-fb-relaxed/{{sample}}/fn.vcf.gz",
    threads: rule_cpus("vcfeval", 64)
    resources:
        mem_mb=rule_mem_gb("vcfeval", 128) * 1024,
        runtime=rule_runtime("vcfeval"),
    params:
        out_dir=f"{OUT_DIR}/vcfeval-fb-relaxed/{{sample}}",
        eval_tool=config.get("eval_tool", "aardvark"),
        docker_arg=lambda wc: f"--docker {config['vcfeval_docker']}" if config.get("vcfeval_docker") else "",
        no_docker="" if config.get("vcfeval_docker") else "--no-docker",
        augref_prefix=f"{AUGREF}#0#",
        min_vcfeval_len=config.get("min_vcfeval_len", 0),
        fb_min_qual=config.get("freebayes_min_qual", 20),
    shell:
        # Rename call CHROMs (unfiltered — keep all variants including lowad/lowdepth)
        "export RTG_MEM=$(({resources.mem_mb} / 1024))g"
        " && mkdir -p {params.out_dir}"
        " && bcftools query -f '%CHROM\\n' {input.call_vcf} | sort -u"
        "    | sed -n '/^{params.augref_prefix}/!s/^\\(.*\\)/\\1\\t{params.augref_prefix}\\1/p'"
        "    > {params.out_dir}/rename-chrs.txt"
        " && if [ -s {params.out_dir}/rename-chrs.txt ]; then"
        "      bcftools annotate --rename-chrs {params.out_dir}/rename-chrs.txt"
        "        {input.call_vcf}"
        "        | awk '/^##contig=/{{id=$0; sub(/.*ID=/, \"\", id); sub(/[,>].*/, \"\", id);"
        "                if(seen[id]++) next}} {{print}}';"
        "    else"
        "      bcftools view {input.call_vcf};"
        "    fi"
        "    | bgzip > {params.out_dir}/call.renamed.vcf.gz"
        " && tabix -fp vcf {params.out_dir}/call.renamed.vcf.gz"
        # QUAL-filter FB VCF, then PASS-filter
        " && bcftools filter -e 'QUAL<{params.fb_min_qual}' -s LowQual {input.fb_vcf} 2>/dev/null"
        "    | bcftools view -f PASS 2>/dev/null"
        "    | bgzip > {params.out_dir}/fb.qfilt.vcf.gz"
        " && tabix -fp vcf {params.out_dir}/fb.qfilt.vcf.gz"
        # Stage to scratch
        " && WORK_TMPDIR=$(mktemp -d \"${{TMPDIR:-{params.out_dir}}}/vcfeval.XXXXXX\")"
        " && trap 'rm -rf \"$WORK_TMPDIR\"' EXIT"
        " && cp {params.out_dir}/call.renamed.vcf.gz"
        "       {params.out_dir}/call.renamed.vcf.gz.tbi"
        "       {params.out_dir}/fb.qfilt.vcf.gz"
        "       {params.out_dir}/fb.qfilt.vcf.gz.tbi"
        "       {input.ref} {input.ref}.fai"
        "       \"$WORK_TMPDIR/\""
        " && {{ [ -f {input.ref}.gzi ]"
        "       && cp {input.ref}.gzi \"$WORK_TMPDIR/\" || true; }}"
        " && python3 scripts/vcfcomp.py {params.eval_tool}"
        "    --truth $WORK_TMPDIR/call.renamed.vcf.gz"
        "    --calls $WORK_TMPDIR/fb.qfilt.vcf.gz"
        "    --ref $WORK_TMPDIR/$(basename {input.ref})"
        "    --out-dir $WORK_TMPDIR"
        "    --threads {threads}"
        "    --min-contig-len {params.min_vcfeval_len}"
        "    {params.docker_arg} {params.no_docker}"
        " && for f in tp-baseline.vcf.gz tp-baseline.vcf.gz.tbi"
        "          fp.vcf.gz fp.vcf.gz.tbi fn.vcf.gz fn.vcf.gz.tbi"
        "          summary.txt snp_roc.tsv.gz non_snp_roc.tsv.gz weighted_roc.tsv.gz"
        "          phasing.txt vcfeval.log progress"
        "          query.vcf.gz query.vcf.gz.tbi truth.vcf.gz truth.vcf.gz.tbi; do"
        "    [ -f \"$WORK_TMPDIR/$f\" ] && cp \"$WORK_TMPDIR/$f\" {params.out_dir}/;"
        "  done"
        " && rm -f {params.out_dir}/call.renamed.vcf.gz"
        "    {params.out_dir}/call.renamed.vcf.gz.tbi"
        "    {params.out_dir}/rename-chrs.txt"
        "    {params.out_dir}/fb.qfilt.vcf.gz"
        "    {params.out_dir}/fb.qfilt.vcf.gz.tbi"

rule vcfeval_fb_recall_per_sample:
    """Relaxed recall: PASS call (truth) vs ALL FB (calls, no QUAL filter).

    Measures what fraction of PASS call variants FB detects at any confidence.
    FN from this run = PASS call variants that FB doesn't see at all.
    """
    input:
        call_vcf=f"{OUT_DIR}/{{sample}}.filtered.vcf.gz" if surject_filtering() else f"{OUT_DIR}/{{sample}}.vcf.gz",
        fb_vcf=f"{OUT_DIR}/{{sample}}.freebayes.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        tp_baseline=f"{OUT_DIR}/vcfeval-fb-recall/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval-fb-recall/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-fb-recall/{{sample}}/fn.vcf.gz",
    threads: rule_cpus("vcfeval", 64)
    resources:
        mem_mb=rule_mem_gb("vcfeval", 128) * 1024,
        runtime=rule_runtime("vcfeval"),
    params:
        out_dir=f"{OUT_DIR}/vcfeval-fb-recall/{{sample}}",
        eval_tool=config.get("eval_tool", "aardvark"),
        docker_arg=lambda wc: f"--docker {config['vcfeval_docker']}" if config.get("vcfeval_docker") else "",
        no_docker="" if config.get("vcfeval_docker") else "--no-docker",
        augref_prefix=f"{AUGREF}#0#",
        min_vcfeval_len=config.get("min_vcfeval_len", 0),
    shell:
        # Rename call CHROMs and PASS-filter (same as standard comparison)
        "export RTG_MEM=$(({resources.mem_mb} / 1024))g"
        " && mkdir -p {params.out_dir}"
        " && bcftools query -f '%CHROM\\n' {input.call_vcf} | sort -u"
        "    | sed -n '/^{params.augref_prefix}/!s/^\\(.*\\)/\\1\\t{params.augref_prefix}\\1/p'"
        "    > {params.out_dir}/rename-chrs.txt"
        " && if [ -s {params.out_dir}/rename-chrs.txt ]; then"
        "      bcftools annotate --rename-chrs {params.out_dir}/rename-chrs.txt"
        "        {input.call_vcf}"
        "        | awk '/^##contig=/{{id=$0; sub(/.*ID=/, \"\", id); sub(/[,>].*/, \"\", id);"
        "                if(seen[id]++) next}} {{print}}';"
        "    else"
        "      bcftools view {input.call_vcf};"
        "    fi"
        "    | bcftools view -f PASS 2>/dev/null"
        "    | bgzip > {params.out_dir}/call.renamed.vcf.gz"
        " && tabix -fp vcf {params.out_dir}/call.renamed.vcf.gz"
        # FB VCF: NO QUAL filter, just use raw (FILTER already set to PASS by freebayes.sh)
        # Stage to scratch
        " && WORK_TMPDIR=$(mktemp -d \"${{TMPDIR:-{params.out_dir}}}/vcfeval.XXXXXX\")"
        " && trap 'rm -rf \"$WORK_TMPDIR\"' EXIT"
        " && cp {params.out_dir}/call.renamed.vcf.gz"
        "       {params.out_dir}/call.renamed.vcf.gz.tbi"
        "       {input.fb_vcf} {input.fb_vcf}.tbi"
        "       {input.ref} {input.ref}.fai"
        "       \"$WORK_TMPDIR/\""
        " && {{ [ -f {input.ref}.gzi ]"
        "       && cp {input.ref}.gzi \"$WORK_TMPDIR/\" || true; }}"
        " && python3 scripts/vcfcomp.py {params.eval_tool}"
        "    --truth $WORK_TMPDIR/call.renamed.vcf.gz"
        "    --calls $WORK_TMPDIR/$(basename {input.fb_vcf})"
        "    --ref $WORK_TMPDIR/$(basename {input.ref})"
        "    --out-dir $WORK_TMPDIR"
        "    --threads {threads}"
        "    --min-contig-len {params.min_vcfeval_len}"
        "    {params.docker_arg} {params.no_docker}"
        " && for f in tp-baseline.vcf.gz tp-baseline.vcf.gz.tbi"
        "          fp.vcf.gz fp.vcf.gz.tbi fn.vcf.gz fn.vcf.gz.tbi"
        "          summary.txt snp_roc.tsv.gz non_snp_roc.tsv.gz weighted_roc.tsv.gz"
        "          phasing.txt vcfeval.log progress"
        "          query.vcf.gz query.vcf.gz.tbi truth.vcf.gz truth.vcf.gz.tbi; do"
        "    [ -f \"$WORK_TMPDIR/$f\" ] && cp \"$WORK_TMPDIR/$f\" {params.out_dir}/;"
        "  done"
        " && rm -f {params.out_dir}/call.renamed.vcf.gz"
        "    {params.out_dir}/call.renamed.vcf.gz.tbi"
        "    {params.out_dir}/rename-chrs.txt"

rule vcfeval_fb_relaxed_combine:
    """Combine PASS precision (TP/FP) with relaxed recall (FN) per sample.

    Precision from standard PASS comparison: TP and FP from vcfeval-fb/pass/
    Recall from recall comparison: FN from vcfeval-fb-recall/ (PASS call vs ALL FB)
    Combined into vcfeval-fb-combined/ for plotting.
    """
    input:
        pass_tp=f"{OUT_DIR}/vcfeval-fb/pass/{{sample}}/tp-baseline.vcf.gz",
        pass_fp=f"{OUT_DIR}/vcfeval-fb/pass/{{sample}}/fp.vcf.gz",
        recall_fn=f"{OUT_DIR}/vcfeval-fb-recall/{{sample}}/fn.vcf.gz",
    output:
        tp=f"{OUT_DIR}/vcfeval-fb-combined/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval-fb-combined/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-fb-combined/{{sample}}/fn.vcf.gz",
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "mkdir -p {OUT_DIR}/vcfeval-fb-combined/{wildcards.sample}"
        " && cp {input.pass_tp} {output.tp}"
        " && cp {input.pass_fp} {output.fp}"
        " && cp {input.recall_fn} {output.fn}"

rule vcfeval_fb_relaxed_compare_plot:
    """Combined comparison plot: PASS precision + relaxed recall"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-fb-combined/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES),
        fp=expand(f"{OUT_DIR}/vcfeval-fb-combined/{{sample}}/fp.vcf.gz", sample=SAMPLES),
        fn=expand(f"{OUT_DIR}/vcfeval-fb-combined/{{sample}}/fn.vcf.gz", sample=SAMPLES),
    output:
        f"{OUT_DIR}/merged.call-vs-fb.relaxed.vcfeval-compare.png",
        f"{OUT_DIR}/merged.call-vs-fb.relaxed.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-fb-combined/{s}" for s in SAMPLES]),
        sample_names=",".join(SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.call-vs-fb.relaxed"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a 'Call (PASS)' --label-b 'FreeBayes (QUAL>={config[freebayes_min_qual]})'"
        " --title '{REF} Call vs FreeBayes — relaxed recall (vcfeval)'"
        " --no-sv"

rule vcfeval_fb_relaxed_lr_compare_plot:
    """Combined comparison plot for long-read: PASS precision + relaxed recall"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-fb-combined/{{sample}}/tp-baseline.vcf.gz", sample=LR_SAMPLES),
        fp=expand(f"{OUT_DIR}/vcfeval-fb-combined/{{sample}}/fp.vcf.gz", sample=LR_SAMPLES),
        fn=expand(f"{OUT_DIR}/vcfeval-fb-combined/{{sample}}/fn.vcf.gz", sample=LR_SAMPLES),
    output:
        f"{OUT_DIR}/merged.lr.call-vs-fb.relaxed.vcfeval-compare.png",
        f"{OUT_DIR}/merged.lr.call-vs-fb.relaxed.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-fb-combined/{s}" for s in LR_SAMPLES]),
        sample_names=",".join(LR_SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.lr.call-vs-fb.relaxed"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a 'Call (PASS)' --label-b 'FreeBayes (QUAL>={config[freebayes_min_qual]})'"
        " --title '{REF} Long-Read Call vs FreeBayes — relaxed recall (vcfeval)'"
        " --no-sv"

rule vcfeval_fb_compare_plot:
    """Aggregate per-sample vcfeval results (FreeBayes) into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-fb/{wc.filt}/{s}" for s in SAMPLES]),
        sample_names=",".join(SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.call-vs-fb.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a Call --label-b FreeBayes"
        " --title '{REF} Call vs FreeBayes (vcfeval)'"
        " --no-sv"

rule vcfeval_fb_detailed_plot:
    """FreeBayes-vs-vg-call vcfeval plot with FB-only split by graph membership
    and stacks by GIAB region; Ts/Tv annotated on SNP bars. Sibling of
    vcfeval_fb_compare_plot — kept alongside it in figure 4b."""
    input:
        # tp-baseline being present implies truth.vcf.gz(.tbi) are in the same
        # vcfeval-fb dir (they're side-effect outputs of vcfeval_fb_per_sample).
        # The R script reads truth.vcf.gz via bcftools view -R.
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
        giab_beds=augref_giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.vcfeval-detailed.png",
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.vcfeval-detailed.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-fb/{wc.filt}/{s}" for s in SAMPLES]),
        sample_names=",".join(SAMPLES),
    resources:
        mem_mb=32000,
        runtime=240,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcfeval-fb-detailed.R"
        " {OUT_DIR}/merged.call-vs-fb.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --giab-beds {input.giab_beds[0]},{input.giab_beds[1]},{input.giab_beds[2]}"
        " --giab-names " + ",".join(GIAB_STRAT_DISPLAY) +
        " --filter-label '{wildcards.filt} only'"
        " --title '{REF} Call vs FreeBayes — detailed (vcfeval)'"

rule vcfeval_fb_chromsplit:
    """Per-contig FP/FN/TP breakdown from vcfeval/aardvark output (FreeBayes)"""
    input:
        fp=f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fn.vcf.gz",
        tp=f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
    output:
        f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        out_dir=f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}",
        subcommand="aardvark-breakdown" if config.get("eval_tool", "aardvark") == "aardvark" else "vcfeval-breakdown",
    shell:
        "python3 scripts/vcfcomp.py {params.subcommand}"
        " --dir {params.out_dir} > {output}"

rule vcfeval_fb_chromsplit_merge:
    """Merge per-sample FreeBayes chromsplit breakdowns into a single long-format TSV"""
    input:
        expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/chromsplit.tsv",
               sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=30,
    run:
        merge_chromsplit_tsv(input, SAMPLES, output[0])

rule vcfeval_fb_chromsplit_plot:
    """Per-contig FP/FN scatter and top-discordant bar chart (FreeBayes)"""
    input:
        tsv=f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.chromsplit.tsv",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv" if annotation_inputs() else [],
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.chromsplit.png",
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.chromsplit-top.png",
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.chromsplit-concordant.png",
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.chromsplit-top-onref.png",
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.chromsplit-concordant-onref.png",
        *([ f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.chromsplit-annot.png"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.chromsplit-giab.png"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        strip_prefix=f"{AUGREF}#0#",
        annot_arg=lambda wc, input: f"--annot {input.annot}" if annotation_inputs() else "",
        giab_arg=lambda wc, input: (
            f"--giab-beds {','.join(input.giab_beds)} --giab-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-chromsplit-plot.R {input.tsv}"
        " {OUT_DIR}/merged.call-vs-fb.{wildcards.filt}"
        " --title '{REF} Call vs FreeBayes Per-Contig'"
        " --strip-prefix '{params.strip_prefix}'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg}"

rule vcfeval_fb_lr_compare_plot:
    """Aggregate long-read per-sample vcfeval results (FreeBayes) into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fp.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fn.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-fb/{wc.filt}/{s}" for s in LR_SAMPLES]),
        sample_names=",".join(LR_SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.lr.call-vs-fb.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a Call --label-b FreeBayes"
        " --title '{REF} Long-Read Call vs FreeBayes (vcfeval)'"
        " --no-sv"

rule vcfeval_fb_lr_chromsplit_merge:
    """Merge long-read per-sample FreeBayes chromsplit breakdowns into a single long-format TSV"""
    input:
        expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/chromsplit.tsv",
               sample=LR_SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=30,
    run:
        merge_chromsplit_tsv(input, LR_SAMPLES, output[0])

rule vcfeval_fb_lr_chromsplit_plot:
    """Long-read per-contig FP/FN scatter and top-discordant bar chart (FreeBayes)"""
    input:
        tsv=f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.chromsplit.tsv",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv" if annotation_inputs() else [],
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.chromsplit.png",
        f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.chromsplit-top.png",
        f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.chromsplit-concordant.png",
        f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.chromsplit-top-onref.png",
        f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.chromsplit-concordant-onref.png",
        *([ f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.chromsplit-annot.png"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.chromsplit-giab.png"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        strip_prefix=f"{AUGREF}#0#",
        annot_arg=lambda wc, input: f"--annot {input.annot}" if annotation_inputs() else "",
        giab_arg=lambda wc, input: (
            f"--giab-beds {','.join(input.giab_beds)} --giab-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-chromsplit-plot.R {input.tsv}"
        " {OUT_DIR}/merged.lr.call-vs-fb.{wildcards.filt}"
        " --title '{REF} Long-Read Call vs FreeBayes Per-Contig'"
        " --strip-prefix '{params.strip_prefix}'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg}"

############################################################################
# DV vs FreeBayes comparison (vcfeval/aardvark)
############################################################################

rule vcfeval_dv_vs_fb_per_sample:
    """Run VCF comparison per sample: DeepVariant (truth) vs FreeBayes (calls)"""
    input:
        dv_vcf=f"{OUT_DIR}/{{sample}}.deepvariant.vcf.gz",
        fb_vcf=f"{OUT_DIR}/{{sample}}.freebayes.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        tp_baseline=f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/fn.vcf.gz",
    threads: rule_cpus("vcfeval", 64)
    resources:
        mem_mb=rule_mem_gb("vcfeval", 128) * 1024,
        runtime=rule_runtime("vcfeval"),
    params:
        out_dir=f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}",
        eval_tool=config.get("eval_tool", "aardvark"),
        docker_arg=lambda wc: f"--docker {config['vcfeval_docker']}" if config.get("vcfeval_docker") else "",
        no_docker="" if config.get("vcfeval_docker") else "--no-docker",
        augref_prefix=f"{AUGREF}#0#",
        min_vcfeval_len=config.get("min_vcfeval_len", 0),
        fb_min_qual=config.get("freebayes_min_qual", 20),
    shell:
        # Both DV and FB VCFs use augref CHROM names — no rename needed
        "export RTG_MEM=$(({resources.mem_mb} / 1024))g"
        " && mkdir -p {params.out_dir}"
        # Pre-filter to PASS if needed (FB gets QUAL filter first)
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      bcftools view -f PASS {input.dv_vcf} 2>/dev/null"
        "        | bgzip > {params.out_dir}/dv.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/dv.pass.vcf.gz;"
        "      bcftools filter -e 'QUAL<{params.fb_min_qual}' -s LowQual {input.fb_vcf} 2>/dev/null"
        "        | bcftools view -f PASS 2>/dev/null"
        "        | bgzip > {params.out_dir}/fb.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/fb.pass.vcf.gz;"
        "    fi"
        # Stage inputs to node-local scratch
        " && WORK_TMPDIR=$(mktemp -d \"${{TMPDIR:-{params.out_dir}}}/vcfeval.XXXXXX\")"
        " && trap 'rm -rf \"$WORK_TMPDIR\"' EXIT"
        " && cp {input.ref} {input.ref}.fai \"$WORK_TMPDIR/\""
        " && {{ [ -f {input.ref}.gzi ]"
        "       && cp {input.ref}.gzi \"$WORK_TMPDIR/\" || true; }}"
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      cp {params.out_dir}/dv.pass.vcf.gz"
        "         {params.out_dir}/dv.pass.vcf.gz.tbi"
        "         {params.out_dir}/fb.pass.vcf.gz"
        "         {params.out_dir}/fb.pass.vcf.gz.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      TRUTH_VCF=$WORK_TMPDIR/dv.pass.vcf.gz;"
        "      CALLS_VCF=$WORK_TMPDIR/fb.pass.vcf.gz;"
        "    else"
        "      cp {input.dv_vcf} {input.dv_vcf}.tbi"
        "         {input.fb_vcf} {input.fb_vcf}.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      TRUTH_VCF=$WORK_TMPDIR/$(basename {input.dv_vcf});"
        "      CALLS_VCF=$WORK_TMPDIR/$(basename {input.fb_vcf});"
        "    fi"
        " && python3 scripts/vcfcomp.py {params.eval_tool}"
        "      --truth \"$TRUTH_VCF\""
        "      --calls \"$CALLS_VCF\""
        "      --ref $WORK_TMPDIR/$(basename {input.ref})"
        "      --out-dir $WORK_TMPDIR/eval_out"
        "      --threads {threads}"
        "      --min-contig-len {params.min_vcfeval_len}"
        "      {params.docker_arg} {params.no_docker}"
        " && for f in tp.vcf.gz tp.vcf.gz.tbi tp-baseline.vcf.gz tp-baseline.vcf.gz.tbi"
        "          fp.vcf.gz fp.vcf.gz.tbi fn.vcf.gz fn.vcf.gz.tbi"
        "          summary.txt non_snp_roc.tsv.gz snp_roc.tsv.gz weighted_roc.tsv.gz"
        "          phasing.txt vcfeval.log progress"
        "          query.vcf.gz query.vcf.gz.tbi truth.vcf.gz truth.vcf.gz.tbi; do"
        "    [ -f \"$WORK_TMPDIR/eval_out/$f\" ] && cp \"$WORK_TMPDIR/eval_out/$f\" {params.out_dir}/;"
        "  done"
        " && rm -f {params.out_dir}/dv.pass.vcf.gz"
        "    {params.out_dir}/dv.pass.vcf.gz.tbi"
        "    {params.out_dir}/fb.pass.vcf.gz"
        "    {params.out_dir}/fb.pass.vcf.gz.tbi"

rule vcfeval_dv_vs_fb_compare_plot:
    """Aggregate per-sample DV-vs-FB vcfeval results into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-dv-vs-fb/{wc.filt}/{s}" for s in SAMPLES]),
        sample_names=",".join(SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.dv-vs-fb.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a DeepVariant --label-b FreeBayes"
        " --title '{REF} DeepVariant vs FreeBayes (vcfeval)'"
        " --no-sv"

rule vcfeval_dv_vs_fb_chromsplit:
    """Per-contig FP/FN/TP breakdown from DV-vs-FB vcfeval output"""
    input:
        fp=f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/fn.vcf.gz",
        tp=f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
    output:
        f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        out_dir=f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}",
        subcommand="aardvark-breakdown" if config.get("eval_tool", "aardvark") == "aardvark" else "vcfeval-breakdown",
    shell:
        "python3 scripts/vcfcomp.py {params.subcommand}"
        " --dir {params.out_dir} > {output}"

rule vcfeval_dv_vs_fb_chromsplit_merge:
    """Merge per-sample DV-vs-FB chromsplit breakdowns"""
    input:
        expand(f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/chromsplit.tsv",
               sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=30,
    run:
        merge_chromsplit_tsv(input, SAMPLES, output[0])

rule vcfeval_dv_vs_fb_chromsplit_plot:
    """DV-vs-FB per-contig concordance plots"""
    input:
        tsv=f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.chromsplit.tsv",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv" if annotation_inputs() else [],
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.chromsplit.png",
        f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.chromsplit-top.png",
        f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.chromsplit-concordant.png",
        f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.chromsplit-top-onref.png",
        f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.chromsplit-concordant-onref.png",
        *([ f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.chromsplit-annot.png"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.dv-vs-fb.{{filt}}.chromsplit-giab.png"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        strip_prefix=f"{AUGREF}#0#",
        annot_arg=lambda wc, input: f"--annot {input.annot}" if annotation_inputs() else "",
        giab_arg=lambda wc, input: (
            f"--giab-beds {','.join(input.giab_beds)} --giab-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-chromsplit-plot.R {input.tsv}"
        " {OUT_DIR}/merged.dv-vs-fb.{wildcards.filt}"
        " --title '{REF} DeepVariant vs FreeBayes Per-Contig'"
        " --strip-prefix '{params.strip_prefix}'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg}"

# Long-read DV-vs-FB comparison
rule vcfeval_dv_vs_fb_lr_compare_plot:
    """Aggregate long-read per-sample DV-vs-FB vcfeval results"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/fp.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/fn.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-dv-vs-fb/{wc.filt}/{s}" for s in LR_SAMPLES]),
        sample_names=",".join(LR_SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.lr.dv-vs-fb.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a DeepVariant --label-b FreeBayes"
        " --title '{REF} Long-Read DeepVariant vs FreeBayes (vcfeval)'"
        " --no-sv"

rule vcfeval_dv_vs_fb_lr_chromsplit_merge:
    """Merge long-read per-sample DV-vs-FB chromsplit breakdowns"""
    input:
        expand(f"{OUT_DIR}/vcfeval-dv-vs-fb/{{filt}}/{{sample}}/chromsplit.tsv",
               sample=LR_SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=30,
    run:
        merge_chromsplit_tsv(input, LR_SAMPLES, output[0])

rule vcfeval_dv_vs_fb_lr_chromsplit_plot:
    """Long-read DV-vs-FB per-contig concordance plots"""
    input:
        tsv=f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.chromsplit.tsv",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv" if annotation_inputs() else [],
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.chromsplit.png",
        f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.chromsplit-top.png",
        f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.chromsplit-concordant.png",
        f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.chromsplit-top-onref.png",
        f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.chromsplit-concordant-onref.png",
        *([ f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.chromsplit-annot.png"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.lr.dv-vs-fb.{{filt}}.chromsplit-giab.png"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        strip_prefix=f"{AUGREF}#0#",
        annot_arg=lambda wc, input: f"--annot {input.annot}" if annotation_inputs() else "",
        giab_arg=lambda wc, input: (
            f"--giab-beds {','.join(input.giab_beds)} --giab-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-chromsplit-plot.R {input.tsv}"
        " {OUT_DIR}/merged.lr.dv-vs-fb.{wildcards.filt}"
        " --title '{REF} Long-Read DeepVariant vs FreeBayes Per-Contig'"
        " --strip-prefix '{params.strip_prefix}'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg}"

############################################################################
# bcftools comparison rules (mirror of FB rules, using .bcftools.vcf.gz)
############################################################################

rule vcfeval_bc_per_sample:
    """Run VCF comparison per sample: call VCF (truth) vs bcftools VCF (calls).

    Mirrors vcfeval_fb_per_sample. vg call emits plain locus names as CHROM
    while bcftools uses the full augref path; the rename step restores the
    prefix on the call VCF so all inputs share the same namespace (idempotent).
    """
    input:
        call_vcf=f"{OUT_DIR}/{{sample}}.filtered.vcf.gz" if surject_filtering() else f"{OUT_DIR}/{{sample}}.vcf.gz",
        bc_vcf=f"{OUT_DIR}/{{sample}}.bcftools.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        tp_baseline=f"{OUT_DIR}/vcfeval-bc/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval-bc/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-bc/{{filt}}/{{sample}}/fn.vcf.gz",
    threads: rule_cpus("vcfeval", 64)
    resources:
        mem_mb=rule_mem_gb("vcfeval", 128) * 1024,
        runtime=rule_runtime("vcfeval"),
    params:
        out_dir=f"{OUT_DIR}/vcfeval-bc/{{filt}}/{{sample}}",
        eval_tool=config.get("eval_tool", "aardvark"),
        docker_arg=lambda wc: f"--docker {config['vcfeval_docker']}" if config.get("vcfeval_docker") else "",
        no_docker="" if config.get("vcfeval_docker") else "--no-docker",
        augref_prefix=f"{AUGREF}#0#",
        min_vcfeval_len=config.get("min_vcfeval_len", 0),
        bc_min_qual=config.get("bcftools_min_qual", 20),
    shell:
        "export RTG_MEM=$(({resources.mem_mb} / 1024))g"
        " && mkdir -p {params.out_dir}"
        " && bcftools query -f '%CHROM\\n' {input.call_vcf} | sort -u"
        "    | sed -n '/^{params.augref_prefix}/!s/^\\(.*\\)/\\1\\t{params.augref_prefix}\\1/p'"
        "    > {params.out_dir}/rename-chrs.txt"
        " && if [ -s {params.out_dir}/rename-chrs.txt ]; then"
        "      bcftools annotate --rename-chrs {params.out_dir}/rename-chrs.txt"
        "        {input.call_vcf}"
        "        | awk '/^##contig=/{{id=$0; sub(/.*ID=/, \"\", id); sub(/[,>].*/, \"\", id);"
        "                if(seen[id]++) next}} {{print}}';"
        "    else"
        "      bcftools view {input.call_vcf};"
        "    fi"
        "    | if [ '{wildcards.filt}' = 'pass' ]; then"
        "        bcftools view -f PASS 2>/dev/null;"
        "      else cat; fi"
        "    | bgzip > {params.out_dir}/call.renamed.vcf.gz"
        " && tabix -fp vcf {params.out_dir}/call.renamed.vcf.gz"
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      bcftools filter -e 'QUAL<{params.bc_min_qual}' -s LowQual {input.bc_vcf} 2>/dev/null"
        "        | bcftools view -f PASS 2>/dev/null"
        "        | bgzip > {params.out_dir}/bc.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/bc.pass.vcf.gz;"
        "    fi"
        " && WORK_TMPDIR=$(mktemp -d \"${{TMPDIR:-{params.out_dir}}}/vcfeval.XXXXXX\")"
        " && trap 'rm -rf \"$WORK_TMPDIR\"' EXIT"
        " && echo \"Staging inputs to $WORK_TMPDIR\""
        " && cp {params.out_dir}/call.renamed.vcf.gz"
        "       {params.out_dir}/call.renamed.vcf.gz.tbi"
        "       {input.ref} {input.ref}.fai"
        "       \"$WORK_TMPDIR/\""
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      cp {params.out_dir}/bc.pass.vcf.gz"
        "         {params.out_dir}/bc.pass.vcf.gz.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      BC_VCF=$WORK_TMPDIR/bc.pass.vcf.gz;"
        "    else"
        "      cp {input.bc_vcf} {input.bc_vcf}.tbi \"$WORK_TMPDIR/\";"
        "      BC_VCF=$WORK_TMPDIR/$(basename {input.bc_vcf});"
        "    fi"
        " && {{ [ -f {input.ref}.gzi ]"
        "       && cp {input.ref}.gzi \"$WORK_TMPDIR/\" || true; }}"
        " && python3 scripts/vcfcomp.py {params.eval_tool}"
        "    --truth $WORK_TMPDIR/call.renamed.vcf.gz"
        "    --calls $BC_VCF"
        "    --ref $WORK_TMPDIR/$(basename {input.ref})"
        "    --out-dir $WORK_TMPDIR"
        "    --threads {threads}"
        "    --min-contig-len {params.min_vcfeval_len}"
        "    {params.docker_arg} {params.no_docker}"
        " && for f in tp-baseline.vcf.gz tp-baseline.vcf.gz.tbi"
        "          fp.vcf.gz fp.vcf.gz.tbi fn.vcf.gz fn.vcf.gz.tbi"
        "          summary.txt snp_roc.tsv.gz non_snp_roc.tsv.gz weighted_roc.tsv.gz"
        "          phasing.txt vcfeval.log progress"
        "          query.vcf.gz query.vcf.gz.tbi truth.vcf.gz truth.vcf.gz.tbi; do"
        "    [ -f \"$WORK_TMPDIR/$f\" ] && cp \"$WORK_TMPDIR/$f\" {params.out_dir}/;"
        "  done"
        " && rm -f {params.out_dir}/call.renamed.vcf.gz"
        "    {params.out_dir}/call.renamed.vcf.gz.tbi"
        "    {params.out_dir}/rename-chrs.txt"
        "    {params.out_dir}/bc.pass.vcf.gz"
        "    {params.out_dir}/bc.pass.vcf.gz.tbi"

rule vcfeval_bc_compare_plot:
    """Aggregate per-sample vcfeval results (bcftools) into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-bc/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-bc/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-bc/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.call-vs-bc.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.call-vs-bc.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-bc/{wc.filt}/{s}" for s in SAMPLES]),
        sample_names=",".join(SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.call-vs-bc.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a Call --label-b bcftools"
        " --title '{REF} Call vs bcftools (vcfeval)'"
        " --no-sv"

rule vcfeval_bc_lr_compare_plot:
    """Aggregate long-read per-sample vcfeval results (bcftools) into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-bc/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-bc/{{filt}}/{{sample}}/fp.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-bc/{{filt}}/{{sample}}/fn.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.lr.call-vs-bc.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.lr.call-vs-bc.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-bc/{wc.filt}/{s}" for s in LR_SAMPLES]),
        sample_names=",".join(LR_SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.lr.call-vs-bc.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a Call --label-b bcftools"
        " --title '{REF} Long-Read Call vs bcftools (vcfeval)'"
        " --no-sv"

rule vcfeval_dv_vs_bc_per_sample:
    """Run VCF comparison per sample: DeepVariant (truth) vs bcftools (calls)"""
    input:
        dv_vcf=f"{OUT_DIR}/{{sample}}.deepvariant.vcf.gz",
        bc_vcf=f"{OUT_DIR}/{{sample}}.bcftools.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        tp_baseline=f"{OUT_DIR}/vcfeval-dv-vs-bc/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval-dv-vs-bc/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-dv-vs-bc/{{filt}}/{{sample}}/fn.vcf.gz",
    threads: rule_cpus("vcfeval", 64)
    resources:
        mem_mb=rule_mem_gb("vcfeval", 128) * 1024,
        runtime=rule_runtime("vcfeval"),
    params:
        out_dir=f"{OUT_DIR}/vcfeval-dv-vs-bc/{{filt}}/{{sample}}",
        eval_tool=config.get("eval_tool", "aardvark"),
        docker_arg=lambda wc: f"--docker {config['vcfeval_docker']}" if config.get("vcfeval_docker") else "",
        no_docker="" if config.get("vcfeval_docker") else "--no-docker",
        augref_prefix=f"{AUGREF}#0#",
        min_vcfeval_len=config.get("min_vcfeval_len", 0),
        bc_min_qual=config.get("bcftools_min_qual", 20),
    shell:
        "export RTG_MEM=$(({resources.mem_mb} / 1024))g"
        " && mkdir -p {params.out_dir}"
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      bcftools view -f PASS {input.dv_vcf} 2>/dev/null"
        "        | bgzip > {params.out_dir}/dv.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/dv.pass.vcf.gz;"
        "      bcftools filter -e 'QUAL<{params.bc_min_qual}' -s LowQual {input.bc_vcf} 2>/dev/null"
        "        | bcftools view -f PASS 2>/dev/null"
        "        | bgzip > {params.out_dir}/bc.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/bc.pass.vcf.gz;"
        "    fi"
        " && WORK_TMPDIR=$(mktemp -d \"${{TMPDIR:-{params.out_dir}}}/vcfeval.XXXXXX\")"
        " && trap 'rm -rf \"$WORK_TMPDIR\"' EXIT"
        " && cp {input.ref} {input.ref}.fai \"$WORK_TMPDIR/\""
        " && {{ [ -f {input.ref}.gzi ]"
        "       && cp {input.ref}.gzi \"$WORK_TMPDIR/\" || true; }}"
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      cp {params.out_dir}/dv.pass.vcf.gz"
        "         {params.out_dir}/dv.pass.vcf.gz.tbi"
        "         {params.out_dir}/bc.pass.vcf.gz"
        "         {params.out_dir}/bc.pass.vcf.gz.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      TRUTH_VCF=$WORK_TMPDIR/dv.pass.vcf.gz;"
        "      CALLS_VCF=$WORK_TMPDIR/bc.pass.vcf.gz;"
        "    else"
        "      cp {input.dv_vcf} {input.dv_vcf}.tbi"
        "         {input.bc_vcf} {input.bc_vcf}.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      TRUTH_VCF=$WORK_TMPDIR/$(basename {input.dv_vcf});"
        "      CALLS_VCF=$WORK_TMPDIR/$(basename {input.bc_vcf});"
        "    fi"
        " && python3 scripts/vcfcomp.py {params.eval_tool}"
        "      --truth \"$TRUTH_VCF\""
        "      --calls \"$CALLS_VCF\""
        "      --ref $WORK_TMPDIR/$(basename {input.ref})"
        "      --out-dir $WORK_TMPDIR/eval_out"
        "      --threads {threads}"
        "      --min-contig-len {params.min_vcfeval_len}"
        "      {params.docker_arg} {params.no_docker}"
        " && for f in tp.vcf.gz tp.vcf.gz.tbi tp-baseline.vcf.gz tp-baseline.vcf.gz.tbi"
        "          fp.vcf.gz fp.vcf.gz.tbi fn.vcf.gz fn.vcf.gz.tbi"
        "          summary.txt non_snp_roc.tsv.gz snp_roc.tsv.gz weighted_roc.tsv.gz"
        "          phasing.txt vcfeval.log progress"
        "          query.vcf.gz query.vcf.gz.tbi truth.vcf.gz truth.vcf.gz.tbi; do"
        "    [ -f \"$WORK_TMPDIR/eval_out/$f\" ] && cp \"$WORK_TMPDIR/eval_out/$f\" {params.out_dir}/;"
        "  done"
        " && rm -f {params.out_dir}/dv.pass.vcf.gz"
        "    {params.out_dir}/dv.pass.vcf.gz.tbi"
        "    {params.out_dir}/bc.pass.vcf.gz"
        "    {params.out_dir}/bc.pass.vcf.gz.tbi"

rule vcfeval_dv_vs_bc_compare_plot:
    """Aggregate per-sample DV-vs-bcftools vcfeval results into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-dv-vs-bc/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-dv-vs-bc/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-dv-vs-bc/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.dv-vs-bc.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.dv-vs-bc.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-dv-vs-bc/{wc.filt}/{s}" for s in SAMPLES]),
        sample_names=",".join(SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.dv-vs-bc.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a DeepVariant --label-b bcftools"
        " --title '{REF} DeepVariant vs bcftools (vcfeval)'"
        " --no-sv"

rule vcfeval_dv_vs_bc_lr_compare_plot:
    """Aggregate long-read per-sample DV-vs-bcftools vcfeval results"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-dv-vs-bc/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-dv-vs-bc/{{filt}}/{{sample}}/fp.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-dv-vs-bc/{{filt}}/{{sample}}/fn.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.lr.dv-vs-bc.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.lr.dv-vs-bc.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-dv-vs-bc/{wc.filt}/{s}" for s in LR_SAMPLES]),
        sample_names=",".join(LR_SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.lr.dv-vs-bc.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a DeepVariant --label-b bcftools"
        " --title '{REF} Long-Read DeepVariant vs bcftools (vcfeval)'"
        " --no-sv"

############################################################################
# On-reference GIAB stratification for vcfeval comparisons
# (direct intersection of on-ref variants with GIAB BEDs)
############################################################################

rule vcfeval_onref_giab:
    """On-ref GIAB stratification for call-vs-DV comparison"""
    input:
        tp=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.onref-giab.tsv",
        f"{OUT_DIR}/merged.call-vs-dv.{{filt}}.onref-giab.png",
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        tp_vcfs=lambda wc, input: ",".join(input.tp),
        fp_vcfs=lambda wc, input: ",".join(input.fp),
        fn_vcfs=lambda wc, input: ",".join(input.fn),
        giab_beds=lambda wc, input: ",".join(input.giab_beds),
        giab_names=",".join(GIAB_STRAT_DISPLAY),
    shell:
        "bash scripts/vcfeval-onref-giab.sh"
        " {params.tp_vcfs} {params.fp_vcfs} {params.fn_vcfs}"
        " {params.giab_beds} {params.giab_names}"
        " {OUT_DIR}/merged.call-vs-dv.{wildcards.filt}.onref-giab.tsv"
        " && ulimit -s unlimited && Rscript scripts/vcfeval-onref-giab-plot.R"
        " {OUT_DIR}/merged.call-vs-dv.{wildcards.filt}.onref-giab.tsv"
        " {OUT_DIR}/merged.call-vs-dv.{wildcards.filt}.onref-giab.png"
        " --title '{REF} Call vs DeepVariant On-Reference'"

rule vcfeval_lr_onref_giab:
    """On-ref GIAB stratification for long-read call-vs-DV comparison"""
    input:
        tp=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fp.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval/{{filt}}/{{sample}}/fn.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.onref-giab.tsv",
        f"{OUT_DIR}/merged.lr.call-vs-dv.{{filt}}.onref-giab.png",
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        tp_vcfs=lambda wc, input: ",".join(input.tp),
        fp_vcfs=lambda wc, input: ",".join(input.fp),
        fn_vcfs=lambda wc, input: ",".join(input.fn),
        giab_beds=lambda wc, input: ",".join(input.giab_beds),
        giab_names=",".join(GIAB_STRAT_DISPLAY),
    shell:
        "bash scripts/vcfeval-onref-giab.sh"
        " {params.tp_vcfs} {params.fp_vcfs} {params.fn_vcfs}"
        " {params.giab_beds} {params.giab_names}"
        " {OUT_DIR}/merged.lr.call-vs-dv.{wildcards.filt}.onref-giab.tsv"
        " && ulimit -s unlimited && Rscript scripts/vcfeval-onref-giab-plot.R"
        " {OUT_DIR}/merged.lr.call-vs-dv.{wildcards.filt}.onref-giab.tsv"
        " {OUT_DIR}/merged.lr.call-vs-dv.{wildcards.filt}.onref-giab.png"
        " --title '{REF} Long-Read Call vs DeepVariant On-Reference'"

rule vcfeval_fb_onref_giab:
    """On-ref GIAB stratification for call-vs-FB comparison"""
    input:
        tp=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.onref-giab.tsv",
        f"{OUT_DIR}/merged.call-vs-fb.{{filt}}.onref-giab.png",
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        tp_vcfs=lambda wc, input: ",".join(input.tp),
        fp_vcfs=lambda wc, input: ",".join(input.fp),
        fn_vcfs=lambda wc, input: ",".join(input.fn),
        giab_beds=lambda wc, input: ",".join(input.giab_beds),
        giab_names=",".join(GIAB_STRAT_DISPLAY),
    shell:
        "bash scripts/vcfeval-onref-giab.sh"
        " {params.tp_vcfs} {params.fp_vcfs} {params.fn_vcfs}"
        " {params.giab_beds} {params.giab_names}"
        " {OUT_DIR}/merged.call-vs-fb.{wildcards.filt}.onref-giab.tsv"
        " && ulimit -s unlimited && Rscript scripts/vcfeval-onref-giab-plot.R"
        " {OUT_DIR}/merged.call-vs-fb.{wildcards.filt}.onref-giab.tsv"
        " {OUT_DIR}/merged.call-vs-fb.{wildcards.filt}.onref-giab.png"
        " --title '{REF} Call vs FreeBayes On-Reference'"

rule vcfeval_fb_lr_onref_giab:
    """On-ref GIAB stratification for long-read call-vs-FB comparison"""
    input:
        tp=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fp.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-fb/{{filt}}/{{sample}}/fn.vcf.gz", sample=LR_SAMPLES, allow_missing=True),
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.onref-giab.tsv",
        f"{OUT_DIR}/merged.lr.call-vs-fb.{{filt}}.onref-giab.png",
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        tp_vcfs=lambda wc, input: ",".join(input.tp),
        fp_vcfs=lambda wc, input: ",".join(input.fp),
        fn_vcfs=lambda wc, input: ",".join(input.fn),
        giab_beds=lambda wc, input: ",".join(input.giab_beds),
        giab_names=",".join(GIAB_STRAT_DISPLAY),
    shell:
        "bash scripts/vcfeval-onref-giab.sh"
        " {params.tp_vcfs} {params.fp_vcfs} {params.fn_vcfs}"
        " {params.giab_beds} {params.giab_names}"
        " {OUT_DIR}/merged.lr.call-vs-fb.{wildcards.filt}.onref-giab.tsv"
        " && ulimit -s unlimited && Rscript scripts/vcfeval-onref-giab-plot.R"
        " {OUT_DIR}/merged.lr.call-vs-fb.{wildcards.filt}.onref-giab.tsv"
        " {OUT_DIR}/merged.lr.call-vs-fb.{wildcards.filt}.onref-giab.png"
        " --title '{REF} Long-Read Call vs FreeBayes On-Reference'"

############################################################################
# Pantree comparison rules (optional — only when pantree_vcf is set)
############################################################################

rule pantree_extract:
    """Pantree VCF → standardized records TSV"""
    input:
        vcf=config["pantree_vcf"] if config.get("pantree_vcf", "") else [],
    output:
        f"{OUT_DIR}/pantree.records.tsv",
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "python3 scripts/pantree-extract.py"
        " --vcf {input.vcf} --output {output}"

rule pantree_stats:
    """Pantree records TSV → standalone variant-type/size/AF plots"""
    input:
        f"{OUT_DIR}/pantree.records.tsv",
    output:
        f"{OUT_DIR}/pantree.vcf-stats.tsv",
        f"{OUT_DIR}/pantree.variant-types.png",
        f"{OUT_DIR}/pantree.size-dist.png",
        f"{OUT_DIR}/pantree.size-dist-log.png",
        f"{OUT_DIR}/pantree.af-spectrum.png",
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input} {OUT_DIR}/pantree"
        " --tsv --title 'Pantree'"

rule deconstruct_records:
    """Export deconstruct VCF as records TSV for comparison"""
    input:
        f"{OUT_DIR}/{OUT_NAME}.tr.vcf.gz",
    output:
        f"{OUT_DIR}/{OUT_NAME}.records.tsv",
    resources:
        mem_mb=int(rule_mem_gb("deconstruct_stats", 512)) * 1024,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-stats.R {input} {OUT_DIR}/{OUT_NAME}"
        " --mode sites --dump-records --records-only"
        " --title '{REF} Deconstruct'"

rule pantree_compare:
    """Side-by-side comparison of our deconstruct vs pantree"""
    input:
        ours=f"{OUT_DIR}/{OUT_NAME}.records.tsv",
        pantree=f"{OUT_DIR}/pantree.records.tsv",
    output:
        f"{OUT_DIR}/{OUT_NAME}.pantree-types.png",
        f"{OUT_DIR}/{OUT_NAME}.pantree-types-pct.png",
        f"{OUT_DIR}/{OUT_NAME}.pantree-size-dist.png",
        f"{OUT_DIR}/{OUT_NAME}.pantree-af.png",
        f"{OUT_DIR}/{OUT_NAME}.pantree-compare.tsv",
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/pantree-compare.R"
        " --ours {input.ours} --pantree {input.pantree}"
        " --prefix {OUT_DIR}/{OUT_NAME}"
        " --title 'Deconstruct vs Pantree'"

rule pantree_density:
    """Pantree records TSV → chromosome density ideogram"""
    input:
        f"{OUT_DIR}/pantree.records.tsv",
    output:
        f"{OUT_DIR}/pantree.density.png",
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        annot_args=" ".join(
            [f"--censat '{config['annot_censat']}'" if config.get("annot_censat", "") else "",
             f"--segdups '{config['annot_segdups']}'" if config.get("annot_segdups", "") else "",
             f"--genes '{config['annot_genes']}'" if config.get("annot_genes", "") else ""]),
    shell:
        "ulimit -s unlimited && Rscript scripts/pantree-density.R"
        " {input} {output}"
        " 'Pantree Off-Reference Variant Density'"
        " --ref {REF}"
        " --bed '{config[refgaps_bed]}'"
        " {params.annot_args}"

############################################################################
# Summary figures — numbered multi-panel composites for quick overview
############################################################################

rule summary_augref:
    """Compose augref segment summary figure"""
    input:
        length_hist=f"{OUT_DIR}/{OUT_NAME}.augref-length-hist.png",
        length_loglog=f"{OUT_DIR}/{OUT_NAME}.augref-length-hist-loglog.png",
        ideogram=f"{OUT_DIR}/{OUT_NAME}.offref-segs.png",
        frequency=f"{OUT_DIR}/{OUT_NAME}.augref-frequency.png",
        frequency_by_size=f"{OUT_DIR}/{OUT_NAME}.augref-frequency-by-size.png",
        annot_summary=[f"{OUT_DIR}/{OUT_NAME}.annot-summary.png"] if annotation_inputs() else [],
        annot_repeats=[f"{OUT_DIR}/{OUT_NAME}.annot-repeats.png"] if config.get("annot_repeats", "") else [],
    output:
        f"{OUT_DIR}/1.augref-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Segment Lengths:{input.length_hist}'",
             f"'Segment Lengths (log-log):{input.length_loglog}'",
             f"'Density Ideogram:{input.ideogram}'",
             f"'Population Frequency:{input.frequency}'",
             f"'Frequency vs Length:{input.frequency_by_size}'"]
            + ([f"'Annotation Overlap:{input.annot_summary[0]}'"] if input.annot_summary else [])
            + ([f"'Repeat Classes:{input.annot_repeats[0]}'"] if input.annot_repeats else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Augref Segment Summary'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_deconstruct:
    """Compose deconstruct variant catalog summary figure"""
    input:
        size_dist=f"{OUT_DIR}/{OUT_NAME}.sites.size-dist-log.png",
        af_spectrum=f"{OUT_DIR}/{OUT_NAME}.sites.af-spectrum.png",
        per_sample=f"{OUT_DIR}/{OUT_NAME}.sites.per-sample-types.png",
        per_sample_pop=f"{OUT_DIR}/{OUT_NAME}.sites.per-sample-types-by-pop.png",
        annot_snp=[f"{OUT_DIR}/{OUT_NAME}.annot-snp-tstv.all.png"] if annotation_inputs() else [],
        giab_strat=[f"{OUT_DIR}/{OUT_NAME}.sites.giab-strat.png"] if giab_strat_configured() else [],
    output:
        f"{OUT_DIR}/2.deconstruct-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Size Distribution:{input.size_dist}'",
             f"'AF Spectrum:{input.af_spectrum}'",
             f"'Per-Sample Types:{input.per_sample}'",
             f"'Per-Sample Types by Super-Population:{input.per_sample_pop}'"]
            + ([f"'SNP Ts/Tv by Annotation:{input.annot_snp[0]}'"] if input.annot_snp else [])
            + ([f"'GIAB Stratification:{input.giab_strat[0]}'"] if input.giab_strat else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Deconstruct Variant Catalog'"
        " --cols 2"
        " --panels {params.panels}"

rule call_summary_panel:
    """Per-sample call summary: annotation-stacked bars with Ts/Tv"""
    input:
        per_sample=f"{OUT_DIR}/merged.call.sites.pass.per-sample-types.tsv",
        vcf=f"{OUT_DIR}/merged.call.vcf.gz",
        annot=[f"{OUT_DIR}/merged.call.sites.pass.annot-exclusive.tsv"] if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/merged.call.sites.pass.call-summary-panel.png",
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        annot_arg=lambda wc, input: f"--annot {input.annot[0]}" if input.annot else "",
        min_sv_size=config.get("min_augref_len", 50),
    shell:
        "ulimit -s unlimited && Rscript scripts/call-summary-panel.R"
        " --per-sample {input.per_sample}"
        " --vcf {input.vcf}"
        " --min-sv-size {params.min_sv_size}"
        " {params.annot_arg}"
        " --output {output}"
        " --title 'Off-reference genotyping'"

rule summary_call:
    """Compose vg call genotyping summary figure"""
    input:
        per_sample=f"{OUT_DIR}/merged.call.sites.pass.per-sample-types.png",
        call_panel=[f"{OUT_DIR}/merged.call.sites.pass.call-summary-panel.png"] if SAMPLES else [],
        giab_per_sample=[f"{OUT_DIR}/merged.call.sites.pass.per-sample-giab-strat.png"] if giab_strat_configured() else [],
        annot_snp=[f"{OUT_DIR}/merged.call.annot-snp-tstv.pass.png"] if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/3.call-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Per-Sample Types (PASS):{input.per_sample}'"]
            + ([f"'Call Summary:{input.call_panel[0]}'"] if input.call_panel else [])
            + ([f"'Per-Sample GIAB:{input.giab_per_sample[0]}'"] if input.giab_per_sample else [])
            + ([f"'SNP Ts/Tv by Annotation:{input.annot_snp[0]}'"] if input.annot_snp else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'vg call Genotyping (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_deepvariant:
    """Compose deepvariant + comparison summary figure"""
    input:
        dv_per_sample=f"{OUT_DIR}/merged.dv.sites.pass.per-sample-types.png",
        vcfeval=f"{OUT_DIR}/merged.call-vs-dv.pass.vcfeval-compare.png",
        annot_snp=[f"{OUT_DIR}/merged.deepvariant.annot-snp-tstv.pass.png"] if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/4.deepvariant-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'DV Per-Sample Types (PASS):{input.dv_per_sample}'",
             f"'Call vs DV (vcfeval):{input.vcfeval}'"]
            + ([f"'DV SNP Ts/Tv by Annotation:{input.annot_snp[0]}'"] if input.annot_snp else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'DeepVariant + Comparison (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_freebayes:
    """Compose FreeBayes summary figure"""
    input:
        fb_types=f"{OUT_DIR}/merged.fb.sites.pass.variant-types.png",
        fb_per_sample=f"{OUT_DIR}/merged.fb.sites.pass.per-sample-types.png",
        vcfeval=f"{OUT_DIR}/merged.call-vs-fb.pass.vcfeval-compare.png",
        vcfeval_detailed=([f"{OUT_DIR}/merged.call-vs-fb.pass.vcfeval-detailed.png"]
                          if giab_strat_configured() else []),
        fb_giab=[f"{OUT_DIR}/merged.fb.sites.pass.giab-strat.png"] if giab_strat_configured() else [],
        annot_snp=[f"{OUT_DIR}/merged.freebayes.annot-snp-tstv.pass.png"] if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/4b.freebayes-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'FB Variant Types (PASS):{input.fb_types}'",
             f"'FB Per-Sample Types (PASS):{input.fb_per_sample}'",
             f"'Call vs FB (vcfeval):{input.vcfeval}'"]
            + ([f"'Call vs FB — detailed:{input.vcfeval_detailed[0]}'"] if input.vcfeval_detailed else [])
            + ([f"'FB GIAB Stratification:{input.fb_giab[0]}'"] if input.fb_giab else [])
            + ([f"'FB SNP Ts/Tv by Annotation:{input.annot_snp[0]}'"] if input.annot_snp else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'FreeBayes (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_bcftools:
    """Compose bcftools summary figure"""
    input:
        bc_types=f"{OUT_DIR}/merged.bc.sites.pass.variant-types.png",
        bc_per_sample=f"{OUT_DIR}/merged.bc.sites.pass.per-sample-types.png",
        vcfeval=f"{OUT_DIR}/merged.call-vs-bc.pass.vcfeval-compare.png",
        bc_giab=[f"{OUT_DIR}/merged.bc.sites.pass.giab-strat.png"] if giab_strat_configured() else [],
        annot_snp=[f"{OUT_DIR}/merged.bcftools.annot-snp-tstv.pass.png"] if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/4d.bcftools-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'BC Variant Types (PASS):{input.bc_types}'",
             f"'BC Per-Sample Types (PASS):{input.bc_per_sample}'",
             f"'Call vs BC (vcfeval):{input.vcfeval}'"]
            + ([f"'BC GIAB Stratification:{input.bc_giab[0]}'"] if input.bc_giab else [])
            + ([f"'BC SNP Ts/Tv by Annotation:{input.annot_snp[0]}'"] if input.annot_snp else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'bcftools (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_concordance:
    """Compose call-vs-DV concordance stratification summary figure"""
    input:
        discordant=f"{OUT_DIR}/merged.call-vs-dv.pass.chromsplit-top.png",
        concordant=f"{OUT_DIR}/merged.call-vs-dv.pass.chromsplit-concordant.png",
        annot=f"{OUT_DIR}/merged.call-vs-dv.pass.chromsplit-annot.png" if annotation_inputs() else [],
        giab=f"{OUT_DIR}/merged.call-vs-dv.pass.chromsplit-giab.png" if giab_strat_configured() else [],
    output:
        f"{OUT_DIR}/5.concordance-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Top Discordant Off-Ref:{input.discordant}'",
             f"'Top Concordant Off-Ref:{input.concordant}'"]
            + ([f"'TP/FP/FN by Annotation:{input.annot}'"] if input.annot else [])
            + ([f"'TP/FP/FN by GIAB Region:{input.giab}'"] if input.giab else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Call vs DeepVariant Concordance — Off-Ref (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_concordance_onref:
    """Compose on-reference concordance summary figure"""
    input:
        discordant=f"{OUT_DIR}/merged.call-vs-dv.pass.chromsplit-top-onref.png",
        concordant=f"{OUT_DIR}/merged.call-vs-dv.pass.chromsplit-concordant-onref.png",
        onref_giab=[f"{OUT_DIR}/merged.call-vs-dv.pass.onref-giab.png"] if giab_strat_configured() else [],
    output:
        f"{OUT_DIR}/5b.concordance-onref-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Top Discordant On-Ref:{input.discordant}'",
             f"'Top Concordant On-Ref:{input.concordant}'"]
            + ([f"'On-Ref GIAB Stratification:{input.onref_giab[0]}'"] if input.onref_giab else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Call vs DeepVariant Concordance — On-Ref (PASS)'"
        " --cols 2"
        " --panels {params.panels}"
        " 'Top Discordant On-Ref:{input.discordant}'"
        " 'Top Concordant On-Ref:{input.concordant}'"

rule summary_coverage:
    """Averaged contig depth summary: pack vs BAM"""
    input:
        f"{OUT_DIR}/contig-depth-summary.png",
    output:
        f"{OUT_DIR}/5c.coverage-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "cp {input} {output}"

rule summary_mapq:
    """MAPQ distribution summary: GAM vs BAM"""
    input:
        f"{OUT_DIR}/mapq-dist.png",
    output:
        f"{OUT_DIR}/5d.mapq-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "cp {input} {output}"

############################################################################
# Long-read merged rules and summary figures
############################################################################

rule longread_call_summary_panel:
    """Long-read per-sample call summary: annotation-stacked bars with Ts/Tv"""
    input:
        per_sample=f"{OUT_DIR}/merged.longread.call.sites.pass.per-sample-types.tsv",
        vcf=f"{OUT_DIR}/merged.longread.call.vcf.gz",
        annot=[f"{OUT_DIR}/merged.longread.call.sites.pass.annot-exclusive.tsv"] if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/merged.longread.call.sites.pass.call-summary-panel.png",
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        annot_arg=lambda wc, input: f"--annot {input.annot[0]}" if input.annot else "",
        min_sv_size=config.get("min_augref_len", 50),
    shell:
        "ulimit -s unlimited && Rscript scripts/call-summary-panel.R"
        " --per-sample {input.per_sample}"
        " --vcf {input.vcf}"
        " --min-sv-size {params.min_sv_size}"
        " {params.annot_arg}"
        " --output {output}"
        " --title 'Off-reference genotyping (long reads)'"

rule longread_contig_depth_summary:
    """Averaged long-read pack + BAM depth → multi-panel summary"""
    input:
        pack_depths=expand("{out}/{s}.contig-depth.tsv", out=OUT_DIR, s=LR_SAMPLES),
        bam_depths=expand("{out}/{s}.bam-depth.tsv", out=OUT_DIR, s=LR_SAMPLES),
        bam_q5_depths=expand("{out}/{s}.bam-depth-q5.tsv", out=OUT_DIR, s=LR_SAMPLES),
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
    output:
        f"{OUT_DIR}/contig-depth-summary.lr.png",
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        pack_arg=lambda wc, input: "--pack-depths " + ",".join(input.pack_depths),
        bam_arg=lambda wc, input: "--bam-depths " + ",".join(input.bam_depths),
        bam_q5_arg=lambda wc, input: "--bam-q5-depths " + ",".join(input.bam_q5_depths),
        min_surject_len=config.get("min_surject_len", 0),
    shell:
        "ulimit -s unlimited && Rscript scripts/contig-depth-summary.R"
        " {params.pack_arg}"
        " {params.bam_arg}"
        " {params.bam_q5_arg}"
        " --segs {input.segs}"
        " --output {output}"
        " --depth-cap 60"
        " --min-surject-len {params.min_surject_len}"
        " --title '{REF} Long-Read Contig Depth'"

rule longread_mapq_dist_plot:
    """Long-read MAPQ distribution: GAM vs BAM"""
    input:
        gam_mapq=expand("{out}/{s}.gam-mapq.tsv", out=OUT_DIR, s=LR_SAMPLES),
        bam_mapq=expand("{out}/{s}.bam-mapq.tsv", out=OUT_DIR, s=LR_SAMPLES),
    output:
        f"{OUT_DIR}/mapq-dist.lr.png",
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        gam_arg=lambda wc, input: "--gam-mapq " + ",".join(input.gam_mapq),
        bam_arg=lambda wc, input: "--bam-mapq " + ",".join(input.bam_mapq),
    shell:
        "ulimit -s unlimited && Rscript scripts/mapq-dist-plot.R"
        " {params.gam_arg}"
        " {params.bam_arg}"
        " --output {output}"
        " --title '{REF} Long-Read MAPQ Distribution'"

rule summary_longread_call:
    """Compose long-read call summary figure"""
    input:
        call_panel=f"{OUT_DIR}/merged.longread.call.sites.pass.call-summary-panel.png",
        giab_per_sample=[f"{OUT_DIR}/merged.longread.call.sites.pass.per-sample-giab-strat.png"] if giab_strat_configured() else [],
    output:
        f"{OUT_DIR}/7.call-summary-longread.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Call Summary:{input.call_panel}'"]
            + ([f"'Per-Sample GIAB:{input.giab_per_sample[0]}'"] if input.giab_per_sample else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Long-Read vg call Genotyping (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_longread_deepvariant:
    """Compose long-read DeepVariant + comparison summary figure"""
    input:
        dv_per_sample=f"{OUT_DIR}/merged.longread.dv.sites.pass.per-sample-types.png",
        vcfeval=f"{OUT_DIR}/merged.lr.call-vs-dv.pass.vcfeval-compare.png",
        dv_per_sample_giab=[f"{OUT_DIR}/merged.longread.dv.sites.pass.per-sample-giab-strat.png"] if giab_strat_configured() else [],
        annot_snp=[f"{OUT_DIR}/merged.longread.dv.sites.pass.variant-types-by-annot.png"] if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/8.deepvariant-summary-longread.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'DV Per-Sample Types (PASS):{input.dv_per_sample}'",
             f"'Call vs DV (vcfeval):{input.vcfeval}'"]
            + ([f"'DV Per-Sample GIAB:{input.dv_per_sample_giab[0]}'"] if input.dv_per_sample_giab else [])
            + ([f"'DV Types by Annotation:{input.annot_snp[0]}'"] if input.annot_snp else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Long-Read DeepVariant + Comparison (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_longread_freebayes:
    """Compose long-read FreeBayes summary figure"""
    input:
        fb_types=f"{OUT_DIR}/merged.longread.fb.sites.pass.variant-types.png",
        fb_per_sample=f"{OUT_DIR}/merged.longread.fb.sites.pass.per-sample-types.png",
        vcfeval=f"{OUT_DIR}/merged.lr.call-vs-fb.pass.vcfeval-compare.png",
        fb_giab=[f"{OUT_DIR}/merged.longread.fb.sites.pass.giab-strat.png"] if giab_strat_configured() else [],
        fb_per_sample_giab=[f"{OUT_DIR}/merged.longread.fb.sites.pass.per-sample-giab-strat.png"] if giab_strat_configured() else [],
        annot_snp=[f"{OUT_DIR}/merged.longread.fb.sites.pass.variant-types-by-annot.png"] if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/8b.freebayes-summary-longread.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'FB Variant Types (PASS):{input.fb_types}'",
             f"'FB Per-Sample Types (PASS):{input.fb_per_sample}'",
             f"'Call vs FB (vcfeval):{input.vcfeval}'"]
            + ([f"'FB GIAB Stratification:{input.fb_giab[0]}'"] if input.fb_giab else [])
            + ([f"'FB Per-Sample GIAB:{input.fb_per_sample_giab[0]}'"] if input.fb_per_sample_giab else [])
            + ([f"'FB Types by Annotation:{input.annot_snp[0]}'"] if input.annot_snp else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Long-Read FreeBayes (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_longread_bcftools:
    """Compose long-read bcftools summary figure"""
    input:
        bc_types=f"{OUT_DIR}/merged.longread.bc.sites.pass.variant-types.png",
        bc_per_sample=f"{OUT_DIR}/merged.longread.bc.sites.pass.per-sample-types.png",
        vcfeval=f"{OUT_DIR}/merged.lr.call-vs-bc.pass.vcfeval-compare.png",
        bc_giab=[f"{OUT_DIR}/merged.longread.bc.sites.pass.giab-strat.png"] if giab_strat_configured() else [],
        bc_per_sample_giab=[f"{OUT_DIR}/merged.longread.bc.sites.pass.per-sample-giab-strat.png"] if giab_strat_configured() else [],
        annot_snp=[f"{OUT_DIR}/merged.longread.bc.sites.pass.variant-types-by-annot.png"] if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/8d.bcftools-summary-longread.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'BC Variant Types (PASS):{input.bc_types}'",
             f"'BC Per-Sample Types (PASS):{input.bc_per_sample}'",
             f"'Call vs BC (vcfeval):{input.vcfeval}'"]
            + ([f"'BC GIAB Stratification:{input.bc_giab[0]}'"] if input.bc_giab else [])
            + ([f"'BC Per-Sample GIAB:{input.bc_per_sample_giab[0]}'"] if input.bc_per_sample_giab else [])
            + ([f"'BC Types by Annotation:{input.annot_snp[0]}'"] if input.annot_snp else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Long-Read bcftools (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_longread_concordance:
    """Compose long-read call-vs-DV concordance summary figure (off-ref)"""
    input:
        discordant=f"{OUT_DIR}/merged.lr.call-vs-dv.pass.chromsplit-top.png",
        concordant=f"{OUT_DIR}/merged.lr.call-vs-dv.pass.chromsplit-concordant.png",
        annot=f"{OUT_DIR}/merged.lr.call-vs-dv.pass.chromsplit-annot.png" if annotation_inputs() else [],
        giab=f"{OUT_DIR}/merged.lr.call-vs-dv.pass.chromsplit-giab.png" if giab_strat_configured() else [],
    output:
        f"{OUT_DIR}/9.concordance-summary-longread.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Top Discordant Off-Ref:{input.discordant}'",
             f"'Top Concordant Off-Ref:{input.concordant}'"]
            + ([f"'TP/FP/FN by Annotation:{input.annot}'"] if input.annot else [])
            + ([f"'TP/FP/FN by GIAB Region:{input.giab}'"] if input.giab else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Long-Read Call vs DeepVariant Concordance — Off-Ref (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_longread_concordance_onref:
    """Compose long-read on-reference concordance summary figure"""
    input:
        discordant=f"{OUT_DIR}/merged.lr.call-vs-dv.pass.chromsplit-top-onref.png",
        concordant=f"{OUT_DIR}/merged.lr.call-vs-dv.pass.chromsplit-concordant-onref.png",
        onref_giab=[f"{OUT_DIR}/merged.lr.call-vs-dv.pass.onref-giab.png"] if giab_strat_configured() else [],
    output:
        f"{OUT_DIR}/9b.concordance-onref-summary-longread.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Top Discordant On-Ref:{input.discordant}'",
             f"'Top Concordant On-Ref:{input.concordant}'"]
            + ([f"'On-Ref GIAB Stratification:{input.onref_giab[0]}'"] if input.onref_giab else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Long-Read Call vs DeepVariant Concordance — On-Ref (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_longread_coverage:
    """Long-read contig depth summary"""
    input:
        f"{OUT_DIR}/contig-depth-summary.lr.png",
    output:
        f"{OUT_DIR}/9c.coverage-summary-longread.png",
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "cp {input} {output}"

rule summary_longread_mapq:
    """Long-read MAPQ distribution summary"""
    input:
        f"{OUT_DIR}/mapq-dist.lr.png",
    output:
        f"{OUT_DIR}/9d.mapq-summary-longread.png",
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "cp {input} {output}"

rule summary_pantree:
    """Compose pantree comparison summary figure"""
    input:
        variant_types=f"{OUT_DIR}/pantree.variant-types.png",
        pantree_types=f"{OUT_DIR}/{OUT_NAME}.pantree-types.png",
        size_dist=f"{OUT_DIR}/{OUT_NAME}.pantree-size-dist.png",
        af=f"{OUT_DIR}/{OUT_NAME}.pantree-af.png",
        density=f"{OUT_DIR}/pantree.density.png",
    output:
        f"{OUT_DIR}/6.pantree-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Pantree Comparison'"
        " --cols 2"
        " --panels"
        " 'Pantree Variant Types:{input.variant_types}'"
        " 'Type Comparison:{input.pantree_types}'"
        " 'Size Distribution:{input.size_dist}'"
        " 'AF Spectrum:{input.af}'"
        " 'Density Ideogram:{input.density}'"

############################################################################
# FreeBayes concordance summary panels
############################################################################

rule summary_freebayes_concordance:
    """Compose FB concordance summary: call-vs-FB and DV-vs-FB (off-ref)"""
    input:
        call_disc=f"{OUT_DIR}/merged.call-vs-fb.pass.chromsplit-top.png",
        call_conc=f"{OUT_DIR}/merged.call-vs-fb.pass.chromsplit-concordant.png",
        dv_disc=f"{OUT_DIR}/merged.dv-vs-fb.pass.chromsplit-top.png",
        dv_conc=f"{OUT_DIR}/merged.dv-vs-fb.pass.chromsplit-concordant.png",
        annot=f"{OUT_DIR}/merged.call-vs-fb.pass.chromsplit-annot.png" if annotation_inputs() else [],
        giab=f"{OUT_DIR}/merged.call-vs-fb.pass.chromsplit-giab.png" if giab_strat_configured() else [],
    output:
        f"{OUT_DIR}/10.freebayes-concordance-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Call vs FB Discordant:{input.call_disc}'",
             f"'DV vs FB Discordant:{input.dv_disc}'",
             f"'Call vs FB Concordant:{input.call_conc}'",
             f"'DV vs FB Concordant:{input.dv_conc}'"]
            + ([f"'Call vs FB by Annotation:{input.annot}'"] if input.annot else [])
            + ([f"'Call vs FB by GIAB:{input.giab}'"] if input.giab else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'FreeBayes Concordance — Off-Ref (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_freebayes_concordance_onref:
    """Compose FB concordance summary: call-vs-FB and DV-vs-FB (on-ref)"""
    input:
        call_disc=f"{OUT_DIR}/merged.call-vs-fb.pass.chromsplit-top-onref.png",
        call_conc=f"{OUT_DIR}/merged.call-vs-fb.pass.chromsplit-concordant-onref.png",
        dv_disc=f"{OUT_DIR}/merged.dv-vs-fb.pass.chromsplit-top-onref.png",
        dv_conc=f"{OUT_DIR}/merged.dv-vs-fb.pass.chromsplit-concordant-onref.png",
    output:
        f"{OUT_DIR}/10b.freebayes-concordance-onref-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'FreeBayes Concordance — On-Ref (PASS)'"
        " --cols 2"
        " --panels"
        " 'Call vs FB Discordant:{input.call_disc}'"
        " 'DV vs FB Discordant:{input.dv_disc}'"
        " 'Call vs FB Concordant:{input.call_conc}'"
        " 'DV vs FB Concordant:{input.dv_conc}'"

rule summary_longread_freebayes_concordance:
    """Compose LR FB concordance summary: call-vs-FB and DV-vs-FB (off-ref)"""
    input:
        call_disc=f"{OUT_DIR}/merged.lr.call-vs-fb.pass.chromsplit-top.png",
        call_conc=f"{OUT_DIR}/merged.lr.call-vs-fb.pass.chromsplit-concordant.png",
        dv_disc=f"{OUT_DIR}/merged.lr.dv-vs-fb.pass.chromsplit-top.png",
        dv_conc=f"{OUT_DIR}/merged.lr.dv-vs-fb.pass.chromsplit-concordant.png",
        annot=f"{OUT_DIR}/merged.lr.call-vs-fb.pass.chromsplit-annot.png" if annotation_inputs() else [],
        giab=f"{OUT_DIR}/merged.lr.call-vs-fb.pass.chromsplit-giab.png" if giab_strat_configured() else [],
    output:
        f"{OUT_DIR}/11.freebayes-concordance-summary-longread.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Call vs FB Discordant:{input.call_disc}'",
             f"'DV vs FB Discordant:{input.dv_disc}'",
             f"'Call vs FB Concordant:{input.call_conc}'",
             f"'DV vs FB Concordant:{input.dv_conc}'"]
            + ([f"'Call vs FB by Annotation:{input.annot}'"] if input.annot else [])
            + ([f"'Call vs FB by GIAB:{input.giab}'"] if input.giab else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Long-Read FreeBayes Concordance — Off-Ref (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_longread_freebayes_concordance_onref:
    """Compose LR FB concordance summary: call-vs-FB and DV-vs-FB (on-ref)"""
    input:
        call_disc=f"{OUT_DIR}/merged.lr.call-vs-fb.pass.chromsplit-top-onref.png",
        call_conc=f"{OUT_DIR}/merged.lr.call-vs-fb.pass.chromsplit-concordant-onref.png",
        dv_disc=f"{OUT_DIR}/merged.lr.dv-vs-fb.pass.chromsplit-top-onref.png",
        dv_conc=f"{OUT_DIR}/merged.lr.dv-vs-fb.pass.chromsplit-concordant-onref.png",
    output:
        f"{OUT_DIR}/11b.freebayes-concordance-onref-summary-longread.png",
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'Long-Read FreeBayes Concordance — On-Ref (PASS)'"
        " --cols 2"
        " --panels"
        " 'Call vs FB Discordant:{input.call_disc}'"
        " 'DV vs FB Discordant:{input.dv_disc}'"
        " 'Call vs FB Concordant:{input.call_conc}'"
        " 'DV vs FB Concordant:{input.dv_conc}'"

############################################################################
# Call vs PanGenie comparison
############################################################################

rule compare_call_pg:
    """Compare merged call and merged PanGenie VCFs per-sample"""
    input:
        call_vcf=lambda wc: f"{OUT_DIR}/merged.call.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.call.vcf.gz",
        pg_vcf=lambda wc: f"{OUT_DIR}/merged.pangenie.normed.vcf.gz" if wc.mode == "variants" else f"{OUT_DIR}/merged.pangenie.vcf.gz",
    output:
        f"{OUT_DIR}/merged.call-vs-pg.{{mode}}.{{filt}}.compare.png",
        f"{OUT_DIR}/merged.call-vs-pg.{{mode}}.{{filt}}.compare.tsv",
    resources:
        mem_mb=256000,
        runtime=2880,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare.R {input.call_vcf} {input.pg_vcf}"
        " {OUT_DIR}/merged.call-vs-pg.{wildcards.mode}.{wildcards.filt}"
        " --mode {wildcards.mode} --filter {wildcards.filt}"
        " --label-a Call --label-b PanGenie"
        " --title '{REF} Call vs PanGenie'"
        " --strip-prefix '{AUGREF}#0#'"
        " --no-sv"

############################################################################
# VCF comparison: Call vs PanGenie (vcfeval or aardvark)
############################################################################

rule vcfeval_pg_per_sample:
    """Run VCF comparison per sample: call VCF (truth) vs PanGenie VCF (calls)

    Dispatches to vcfeval or aardvark based on config['eval_tool'].
    When filt=pass, both VCFs are pre-filtered to PASS before comparison
    (aardvark/vcfeval strip FILTER, so post-hoc filtering doesn't work).

    vg call emits plain locus names (e.g. 'chr1') as CHROM while PanGenie
    uses the full augref path (e.g. 'augref_CHM13#0#chr1').  The rename step
    restores the prefix on the call VCF so all inputs share the same namespace.
    It is idempotent: only renames CHROMs that lack the prefix.
    """
    input:
        call_vcf=f"{OUT_DIR}/{{sample}}.filtered.vcf.gz" if surject_filtering() else f"{OUT_DIR}/{{sample}}.vcf.gz",
        pg_vcf=f"{OUT_DIR}/{{sample}}.pangenie.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        tp_baseline=f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/fn.vcf.gz",
    threads: rule_cpus("vcfeval", 64)
    resources:
        mem_mb=rule_mem_gb("vcfeval", 128) * 1024,
        runtime=rule_runtime("vcfeval"),
    params:
        out_dir=f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}",
        eval_tool=config.get("eval_tool", "aardvark"),
        docker_arg=lambda wc: f"--docker {config['vcfeval_docker']}" if config.get("vcfeval_docker") else "",
        no_docker="" if config.get("vcfeval_docker") else "--no-docker",
        augref_prefix=f"{AUGREF}#0#",
        min_vcfeval_len=config.get("min_vcfeval_len", 0),
    shell:
        # Build contig rename map: only rename CHROMs missing the augref prefix
        # (idempotent: works with both old vg [plain names] and fixed vg [augref names])
        "export RTG_MEM=$(({resources.mem_mb} / 1024))g"
        " && mkdir -p {params.out_dir}"
        " && bcftools query -f '%CHROM\\n' {input.call_vcf} | sort -u"
        "    | sed -n '/^{params.augref_prefix}/!s/^\\(.*\\)/\\1\\t{params.augref_prefix}\\1/p'"
        "    > {params.out_dir}/rename-chrs.txt"
        # Rename chroms (no-op if rename file is empty); optionally pre-filter to PASS
        " && if [ -s {params.out_dir}/rename-chrs.txt ]; then"
        "      bcftools annotate --rename-chrs {params.out_dir}/rename-chrs.txt"
        "        {input.call_vcf}"
        "        | awk '/^##contig=/{{id=$0; sub(/.*ID=/, \"\", id); sub(/[,>].*/, \"\", id);"
        "                if(seen[id]++) next}} {{print}}';"
        "    else"
        "      bcftools view {input.call_vcf};"
        "    fi"
        "    | if [ '{wildcards.filt}' = 'pass' ]; then"
        "        bcftools view -f PASS 2>/dev/null;"
        "      else cat; fi"
        "    | bgzip > {params.out_dir}/call.renamed.vcf.gz"
        " && tabix -fp vcf {params.out_dir}/call.renamed.vcf.gz"
        # Pre-filter PG VCF to PASS if needed
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      bcftools view -f PASS {input.pg_vcf} 2>/dev/null"
        "        | bgzip > {params.out_dir}/pg.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/pg.pass.vcf.gz;"
        "    fi"
        # Stage inputs to node-local scratch for fast I/O
        " && WORK_TMPDIR=$(mktemp -d \"${{TMPDIR:-{params.out_dir}}}/vcfeval.XXXXXX\")"
        " && trap 'rm -rf \"$WORK_TMPDIR\"' EXIT"
        " && echo \"Staging inputs to $WORK_TMPDIR\""
        " && cp {params.out_dir}/call.renamed.vcf.gz"
        "       {params.out_dir}/call.renamed.vcf.gz.tbi"
        "       {input.ref} {input.ref}.fai"
        "       \"$WORK_TMPDIR/\""
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      cp {params.out_dir}/pg.pass.vcf.gz"
        "         {params.out_dir}/pg.pass.vcf.gz.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      PG_VCF=$WORK_TMPDIR/pg.pass.vcf.gz;"
        "    else"
        "      cp {input.pg_vcf} {input.pg_vcf}.tbi \"$WORK_TMPDIR/\";"
        "      PG_VCF=$WORK_TMPDIR/$(basename {input.pg_vcf});"
        "    fi"
        " && {{ [ -f {input.ref}.gzi ]"
        "       && cp {input.ref}.gzi \"$WORK_TMPDIR/\" || true; }}"
        " && python3 scripts/vcfcomp.py {params.eval_tool}"
        "    --truth $WORK_TMPDIR/call.renamed.vcf.gz"
        "    --calls $PG_VCF"
        "    --ref $WORK_TMPDIR/$(basename {input.ref})"
        "    --out-dir $WORK_TMPDIR"
        "    --threads {threads}"
        "    --min-contig-len {params.min_vcfeval_len}"
        "    {params.docker_arg} {params.no_docker}"
        # Copy results back from local scratch
        " && for f in tp-baseline.vcf.gz tp-baseline.vcf.gz.tbi"
        "          fp.vcf.gz fp.vcf.gz.tbi fn.vcf.gz fn.vcf.gz.tbi"
        "          summary.txt snp_roc.tsv.gz non_snp_roc.tsv.gz weighted_roc.tsv.gz"
        "          phasing.txt vcfeval.log progress"
        "          query.vcf.gz query.vcf.gz.tbi truth.vcf.gz truth.vcf.gz.tbi; do"
        "    [ -f \"$WORK_TMPDIR/$f\" ] && cp \"$WORK_TMPDIR/$f\" {params.out_dir}/;"
        "  done"
        # Clean up staged files on shared storage
        " && rm -f {params.out_dir}/call.renamed.vcf.gz"
        "    {params.out_dir}/call.renamed.vcf.gz.tbi"
        "    {params.out_dir}/rename-chrs.txt"
        "    {params.out_dir}/pg.pass.vcf.gz"
        "    {params.out_dir}/pg.pass.vcf.gz.tbi"

rule vcfeval_pg_compare_plot:
    """Aggregate per-sample vcfeval results (PanGenie) into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-pg/{wc.filt}/{s}" for s in SAMPLES]),
        sample_names=",".join(SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.call-vs-pg.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a Call --label-b PanGenie"
        " --title '{REF} Call vs PanGenie (vcfeval)'"
        " --no-sv"

rule vcfeval_pg_chromsplit:
    """Per-contig FP/FN/TP breakdown from vcfeval/aardvark output (PanGenie)"""
    input:
        fp=f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/fn.vcf.gz",
        tp=f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
    output:
        f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        out_dir=f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}",
        subcommand="aardvark-breakdown" if config.get("eval_tool", "aardvark") == "aardvark" else "vcfeval-breakdown",
    shell:
        "python3 scripts/vcfcomp.py {params.subcommand}"
        " --dir {params.out_dir} > {output}"

rule vcfeval_pg_chromsplit_merge:
    """Merge per-sample PanGenie chromsplit breakdowns into a single long-format TSV"""
    input:
        expand(f"{OUT_DIR}/vcfeval-pg/{{filt}}/{{sample}}/chromsplit.tsv",
               sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=30,
    run:
        merge_chromsplit_tsv(input, SAMPLES, output[0])

rule vcfeval_pg_chromsplit_plot:
    """Per-contig FP/FN scatter and top-discordant bar chart (PanGenie)"""
    input:
        tsv=f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.chromsplit.tsv",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv" if annotation_inputs() else [],
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.chromsplit.png",
        f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.chromsplit-top.png",
        f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.chromsplit-concordant.png",
        f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.chromsplit-top-onref.png",
        f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.chromsplit-concordant-onref.png",
        *([ f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.chromsplit-annot.png"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.call-vs-pg.{{filt}}.chromsplit-giab.png"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        strip_prefix=f"{AUGREF}#0#",
        annot_arg=lambda wc, input: f"--annot {input.annot}" if annotation_inputs() else "",
        giab_arg=lambda wc, input: (
            f"--giab-beds {','.join(input.giab_beds)} --giab-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-chromsplit-plot.R {input.tsv}"
        " {OUT_DIR}/merged.call-vs-pg.{wildcards.filt}"
        " --title '{REF} Call vs PanGenie Per-Contig'"
        " --strip-prefix '{params.strip_prefix}'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg}"

############################################################################
# DV vs PanGenie comparison (vcfeval/aardvark)
############################################################################

rule vcfeval_dv_vs_pg_per_sample:
    """Run VCF comparison per sample: DeepVariant (truth) vs PanGenie (calls)"""
    input:
        dv_vcf=f"{OUT_DIR}/{{sample}}.deepvariant.vcf.gz",
        pg_vcf=f"{OUT_DIR}/{{sample}}.pangenie.vcf.gz",
        ref=f"{OUT_DIR}/{OUT_NAME}.fa.gz",
        paths=f"{OUT_DIR}/{OUT_NAME}.filtered-paths.txt" if surject_filtering() else [],
    output:
        tp_baseline=f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
        fp=f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/fn.vcf.gz",
    threads: rule_cpus("vcfeval", 64)
    resources:
        mem_mb=rule_mem_gb("vcfeval", 128) * 1024,
        runtime=rule_runtime("vcfeval"),
    params:
        out_dir=f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}",
        eval_tool=config.get("eval_tool", "aardvark"),
        docker_arg=lambda wc: f"--docker {config['vcfeval_docker']}" if config.get("vcfeval_docker") else "",
        no_docker="" if config.get("vcfeval_docker") else "--no-docker",
        augref_prefix=f"{AUGREF}#0#",
        min_vcfeval_len=config.get("min_vcfeval_len", 0),
    shell:
        # Both DV and PG VCFs use augref CHROM names — no rename needed
        "export RTG_MEM=$(({resources.mem_mb} / 1024))g"
        " && mkdir -p {params.out_dir}"
        # Pre-filter to PASS if needed
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      bcftools view -f PASS {input.dv_vcf} 2>/dev/null"
        "        | bgzip > {params.out_dir}/dv.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/dv.pass.vcf.gz;"
        "      bcftools view -f PASS {input.pg_vcf} 2>/dev/null"
        "        | bgzip > {params.out_dir}/pg.pass.vcf.gz"
        "      && tabix -fp vcf {params.out_dir}/pg.pass.vcf.gz;"
        "    fi"
        # Stage inputs to node-local scratch
        " && WORK_TMPDIR=$(mktemp -d \"${{TMPDIR:-{params.out_dir}}}/vcfeval.XXXXXX\")"
        " && trap 'rm -rf \"$WORK_TMPDIR\"' EXIT"
        " && cp {input.ref} {input.ref}.fai \"$WORK_TMPDIR/\""
        " && {{ [ -f {input.ref}.gzi ]"
        "       && cp {input.ref}.gzi \"$WORK_TMPDIR/\" || true; }}"
        " && if [ '{wildcards.filt}' = 'pass' ]; then"
        "      cp {params.out_dir}/dv.pass.vcf.gz"
        "         {params.out_dir}/dv.pass.vcf.gz.tbi"
        "         {params.out_dir}/pg.pass.vcf.gz"
        "         {params.out_dir}/pg.pass.vcf.gz.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      TRUTH_VCF=$WORK_TMPDIR/dv.pass.vcf.gz;"
        "      CALLS_VCF=$WORK_TMPDIR/pg.pass.vcf.gz;"
        "    else"
        "      cp {input.dv_vcf} {input.dv_vcf}.tbi"
        "         {input.pg_vcf} {input.pg_vcf}.tbi"
        "         \"$WORK_TMPDIR/\";"
        "      TRUTH_VCF=$WORK_TMPDIR/$(basename {input.dv_vcf});"
        "      CALLS_VCF=$WORK_TMPDIR/$(basename {input.pg_vcf});"
        "    fi"
        " && python3 scripts/vcfcomp.py {params.eval_tool}"
        "      --truth \"$TRUTH_VCF\""
        "      --calls \"$CALLS_VCF\""
        "      --ref $WORK_TMPDIR/$(basename {input.ref})"
        "      --out-dir $WORK_TMPDIR/eval_out"
        "      --threads {threads}"
        "      --min-contig-len {params.min_vcfeval_len}"
        "      {params.docker_arg} {params.no_docker}"
        " && for f in tp.vcf.gz tp.vcf.gz.tbi tp-baseline.vcf.gz tp-baseline.vcf.gz.tbi"
        "          fp.vcf.gz fp.vcf.gz.tbi fn.vcf.gz fn.vcf.gz.tbi"
        "          summary.txt non_snp_roc.tsv.gz snp_roc.tsv.gz weighted_roc.tsv.gz"
        "          phasing.txt vcfeval.log progress"
        "          query.vcf.gz query.vcf.gz.tbi truth.vcf.gz truth.vcf.gz.tbi; do"
        "    [ -f \"$WORK_TMPDIR/eval_out/$f\" ] && cp \"$WORK_TMPDIR/eval_out/$f\" {params.out_dir}/;"
        "  done"
        " && rm -f {params.out_dir}/dv.pass.vcf.gz"
        "    {params.out_dir}/dv.pass.vcf.gz.tbi"
        "    {params.out_dir}/pg.pass.vcf.gz"
        "    {params.out_dir}/pg.pass.vcf.gz.tbi"

rule vcfeval_dv_vs_pg_compare_plot:
    """Aggregate per-sample DV-vs-PG vcfeval results into comparison plot"""
    input:
        tp_baseline=expand(f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/tp-baseline.vcf.gz", sample=SAMPLES, allow_missing=True),
        fp=expand(f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/fp.vcf.gz", sample=SAMPLES, allow_missing=True),
        fn=expand(f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/fn.vcf.gz", sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.vcfeval-compare.png",
        f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.vcfeval-compare.tsv",
    params:
        vcfeval_dirs=lambda wc, input: ",".join(
            [f"{OUT_DIR}/vcfeval-dv-vs-pg/{wc.filt}/{s}" for s in SAMPLES]),
        sample_names=",".join(SAMPLES),
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-compare-vcfeval.R"
        " {OUT_DIR}/merged.dv-vs-pg.{wildcards.filt}"
        " --vcfeval-dirs {params.vcfeval_dirs}"
        " --samples {params.sample_names}"
        " --label-a DeepVariant --label-b PanGenie"
        " --title '{REF} DeepVariant vs PanGenie (vcfeval)'"
        " --no-sv"

rule vcfeval_dv_vs_pg_chromsplit:
    """Per-contig FP/FN/TP breakdown from DV-vs-PG vcfeval output"""
    input:
        fp=f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/fp.vcf.gz",
        fn=f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/fn.vcf.gz",
        tp=f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/tp-baseline.vcf.gz",
    output:
        f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        out_dir=f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}",
        subcommand="aardvark-breakdown" if config.get("eval_tool", "aardvark") == "aardvark" else "vcfeval-breakdown",
    shell:
        "python3 scripts/vcfcomp.py {params.subcommand}"
        " --dir {params.out_dir} > {output}"

rule vcfeval_dv_vs_pg_chromsplit_merge:
    """Merge per-sample DV-vs-PG chromsplit breakdowns"""
    input:
        expand(f"{OUT_DIR}/vcfeval-dv-vs-pg/{{filt}}/{{sample}}/chromsplit.tsv",
               sample=SAMPLES, allow_missing=True),
    output:
        f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.chromsplit.tsv",
    resources:
        mem_mb=8000,
        runtime=30,
    run:
        merge_chromsplit_tsv(input, SAMPLES, output[0])

rule vcfeval_dv_vs_pg_chromsplit_plot:
    """DV-vs-PG per-contig concordance plots"""
    input:
        tsv=f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.chromsplit.tsv",
        segs=f"{OUT_DIR}/{OUT_NAME}.augref-segs.tsv",
        annot=f"{OUT_DIR}/{OUT_NAME}.annot-per-segment.tsv" if annotation_inputs() else [],
        giab_beds=giab_strat_beds(),
    output:
        f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.chromsplit.png",
        f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.chromsplit-top.png",
        f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.chromsplit-concordant.png",
        f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.chromsplit-top-onref.png",
        f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.chromsplit-concordant-onref.png",
        *([ f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.chromsplit-annot.png"]
          if annotation_inputs() else []),
        *([ f"{OUT_DIR}/merged.dv-vs-pg.{{filt}}.chromsplit-giab.png"]
          if giab_strat_configured() else []),
    resources:
        mem_mb=32000,
        runtime=120,
    params:
        strip_prefix=f"{AUGREF}#0#",
        annot_arg=lambda wc, input: f"--annot {input.annot}" if annotation_inputs() else "",
        giab_arg=lambda wc, input: (
            f"--giab-beds {','.join(input.giab_beds)} --giab-names {','.join(GIAB_STRAT_DISPLAY)}"
            if giab_strat_configured() else ""),
    shell:
        "ulimit -s unlimited && Rscript scripts/vcf-chromsplit-plot.R {input.tsv}"
        " {OUT_DIR}/merged.dv-vs-pg.{wildcards.filt}"
        " --title '{REF} DeepVariant vs PanGenie Per-Contig'"
        " --strip-prefix '{params.strip_prefix}'"
        " --segs {input.segs}"
        " {params.annot_arg} {params.giab_arg}"

############################################################################
# PanGenie summary panels
############################################################################

rule summary_pangenie:
    """Compose PanGenie summary figure"""
    input:
        pg_types=f"{OUT_DIR}/merged.pg.sites.pass.variant-types.png",
        pg_per_sample=f"{OUT_DIR}/merged.pg.sites.pass.per-sample-types.png",
        vcfeval=f"{OUT_DIR}/merged.call-vs-pg.pass.vcfeval-compare.png",
        pg_giab=[f"{OUT_DIR}/merged.pg.sites.pass.giab-strat.png"] if giab_strat_configured() else [],
        annot_snp=[f"{OUT_DIR}/merged.pangenie.annot-snp-tstv.pass.png"] if annotation_inputs() else [],
    output:
        f"{OUT_DIR}/4c.pangenie-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'PG Variant Types (PASS):{input.pg_types}'",
             f"'PG Per-Sample Types (PASS):{input.pg_per_sample}'",
             f"'Call vs PG (vcfeval):{input.vcfeval}'"]
            + ([f"'PG GIAB Stratification:{input.pg_giab[0]}'"] if input.pg_giab else [])
            + ([f"'PG SNP Ts/Tv by Annotation:{input.annot_snp[0]}'"] if input.annot_snp else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'PanGenie (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

############################################################################
# PanGenie concordance summary panels
############################################################################

rule summary_pangenie_concordance:
    """Compose PG concordance summary: call-vs-PG and DV-vs-PG (off-ref)"""
    input:
        call_disc=f"{OUT_DIR}/merged.call-vs-pg.pass.chromsplit-top.png",
        call_conc=f"{OUT_DIR}/merged.call-vs-pg.pass.chromsplit-concordant.png",
        dv_disc=f"{OUT_DIR}/merged.dv-vs-pg.pass.chromsplit-top.png",
        dv_conc=f"{OUT_DIR}/merged.dv-vs-pg.pass.chromsplit-concordant.png",
        annot=f"{OUT_DIR}/merged.call-vs-pg.pass.chromsplit-annot.png" if annotation_inputs() else [],
        giab=f"{OUT_DIR}/merged.call-vs-pg.pass.chromsplit-giab.png" if giab_strat_configured() else [],
    output:
        f"{OUT_DIR}/12.pangenie-concordance-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    params:
        panels=lambda wc, input: " ".join(
            [f"'Call vs PG Discordant:{input.call_disc}'",
             f"'DV vs PG Discordant:{input.dv_disc}'",
             f"'Call vs PG Concordant:{input.call_conc}'",
             f"'DV vs PG Concordant:{input.dv_conc}'"]
            + ([f"'Call vs PG by Annotation:{input.annot}'"] if input.annot else [])
            + ([f"'Call vs PG by GIAB:{input.giab}'"] if input.giab else [])
        ),
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'PanGenie Concordance — Off-Ref (PASS)'"
        " --cols 2"
        " --panels {params.panels}"

rule summary_pangenie_concordance_onref:
    """Compose PG concordance summary: call-vs-PG and DV-vs-PG (on-ref)"""
    input:
        call_disc=f"{OUT_DIR}/merged.call-vs-pg.pass.chromsplit-top-onref.png",
        call_conc=f"{OUT_DIR}/merged.call-vs-pg.pass.chromsplit-concordant-onref.png",
        dv_disc=f"{OUT_DIR}/merged.dv-vs-pg.pass.chromsplit-top-onref.png",
        dv_conc=f"{OUT_DIR}/merged.dv-vs-pg.pass.chromsplit-concordant-onref.png",
    output:
        f"{OUT_DIR}/12b.pangenie-concordance-onref-summary.png",
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "python3 scripts/compose-summary.py"
        " --output {output}"
        " --title 'PanGenie Concordance — On-Ref (PASS)'"
        " --cols 2"
        " --panels"
        " 'Call vs PG Discordant:{input.call_disc}'"
        " 'DV vs PG Discordant:{input.dv_disc}'"
        " 'Call vs PG Concordant:{input.call_conc}'"
        " 'DV vs PG Concordant:{input.dv_conc}'"

