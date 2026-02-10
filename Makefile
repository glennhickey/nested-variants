############################################################################
# Makefile — nested-variants pipeline
#
# Default config is for local chr20 testing.
# Override via config.local.mk or command line:
#   make EXEC_MODE=slurm REF=CHM13 VG='...' all
############################################################################

# Load default config, then user overrides (if present)
include config.mk
-include config.local.mk

# Compute execution flag from mode
ifeq ($(EXEC_MODE),local)
  EXEC_FLAG := --local
else
  EXEC_FLAG :=
endif

# Augmented reference sample name (used by deconstruct, call, split-vcf)
AUGREF := augref_$(REF)

# GBZ used for read mapping (giraffe).  Falls back to the pipeline-built GBZ.
GIRAFFE_GBZ := $(if $(MAP_GBZ),$(MAP_GBZ),$(OUT_DIR)/$(OUT_NAME).gbz)

# Derived file paths
GBZ       := $(OUT_DIR)/$(OUT_NAME).gbz
VCF       := $(OUT_DIR)/$(OUT_NAME).vcf.gz
ONREF     := $(OUT_DIR)/$(OUT_NAME).onref.vcf.gz
NESTEDREF := $(OUT_DIR)/$(OUT_NAME).nestedref.vcf.gz
OFFREF    := $(OUT_DIR)/$(OUT_NAME).offref.vcf.gz
PLOT      := $(OUT_DIR)/$(OUT_NAME).offref.png
BAM       := $(OUT_DIR)/$(SAMPLE).bam
FASTA     := $(OUT_DIR)/$(OUT_NAME).fa.gz
HAPL_INDEX := $(OUT_DIR)/$(OUT_NAME).hapl
CALL_VCF  := $(OUT_DIR)/$(SAMPLE).vcf.gz
DV_VCF    := $(OUT_DIR)/$(SAMPLE).deepvariant.vcf.gz
AUGREF_SEGS := $(OUT_DIR)/$(OUT_NAME).augref-segs.tsv
LENGTH_HIST := $(OUT_DIR)/$(OUT_NAME).augref-length-hist.png
CALL_PLOT  := $(OUT_DIR)/$(SAMPLE).call-offref.png
DV_PLOT    := $(OUT_DIR)/$(SAMPLE).dv-offref.png

# Build SLURM option flags
SLURM_OPTS := --cpus $(CPUS) --mem $(MEM) --time $(TIME) --partition $(PARTITION)

# Build deconstruct-specific option flags
DECON_OPTS :=
ifneq ($(CLUSTER),)
  DECON_OPTS += --cluster $(CLUSTER)
endif
ifneq ($(SNARLS),)
  DECON_OPTS += --snarls $(SNARLS)
endif
ifneq ($(STAR_ALLELE),)
  DECON_OPTS += --star-allele
endif

# Build BED overlay flag for plots
BED_FLAG :=
ifneq ($(REFGAPS_BED),)
  BED_FLAG := $(REFGAPS_BED)
endif

SCRIPTS := $(CURDIR)/scripts

############################################################################
# Top-level targets
############################################################################

.PHONY: all analysis paths deconstruct genotype surject haplotypes fasta deepvariant split-vcf length-hist plots call-plots dv-plots test clean help

all: analysis

analysis: split-vcf plots

## paths: VG file(s) → GBZ via augmented reference paths
paths: $(GBZ)

$(GBZ): scripts/paths.sh
	$(SCRIPTS)/paths.sh \
		--vg '$(VG)' \
		--ref $(REF) \
		--out-dir $(OUT_DIR) \
		--out-name $(OUT_NAME).gfa.gz \
		--min-augref-len $(MIN_AUGREF_LEN) \
		$(SLURM_OPTS) $(EXEC_FLAG)

## deconstruct: GBZ → VCF via vg deconstruct
deconstruct: $(VCF)

$(VCF): $(GBZ) scripts/deconstruct.sh
	$(SCRIPTS)/deconstruct.sh \
		--gbz $(GBZ) \
		--ref $(AUGREF) \
		--out-dir $(OUT_DIR) \
		--out-name $(OUT_NAME).vcf.gz \
		$(DECON_OPTS) $(SLURM_OPTS) $(EXEC_FLAG)

## genotype: GBZ + reads → GAM → sample VCF (optional)
genotype: $(OUT_DIR)/$(SAMPLE).vcf.gz

$(OUT_DIR)/$(SAMPLE).gam: scripts/giraffe.sh
	$(SCRIPTS)/giraffe.sh \
		--gbz $(GIRAFFE_GBZ) \
		--hapl $(HAPL) \
		--reads $(READS) \
		--sample $(SAMPLE) \
		--out-dir $(OUT_DIR) \
		--out-name $(SAMPLE).gam \
		$(SLURM_OPTS) $(EXEC_FLAG)

$(OUT_DIR)/$(SAMPLE).vcf.gz: $(OUT_DIR)/$(SAMPLE).gam $(GBZ) scripts/call.sh
	$(SCRIPTS)/call.sh \
		--gbz $(GBZ) \
		--gam $(OUT_DIR)/$(SAMPLE).gam \
		--ref $(AUGREF) \
		--sample $(SAMPLE) \
		--out-dir $(OUT_DIR) \
		--out-name $(SAMPLE).vcf.gz \
		$(SLURM_OPTS) $(EXEC_FLAG)

## surject: GAM → sorted BAM via vg surject (onto augmented reference paths)
surject: $(BAM)

$(BAM): $(OUT_DIR)/$(SAMPLE).gam $(GBZ) scripts/surject.sh
	$(SCRIPTS)/surject.sh \
		--gbz $(GBZ) \
		--gam $(OUT_DIR)/$(SAMPLE).gam \
		--ref $(AUGREF) \
		--sample $(SAMPLE) \
		--out-dir $(OUT_DIR) \
		--out-name $(SAMPLE).bam \
		$(SLURM_OPTS) $(EXEC_FLAG)

## haplotypes: GBZ → .hapl index for giraffe haplotype-aware mapping
haplotypes: $(HAPL_INDEX)

$(HAPL_INDEX): $(GBZ) scripts/haplotypes.sh
	$(SCRIPTS)/haplotypes.sh \
		--gbz $(GBZ) \
		--ref $(REF) \
		--out-dir $(OUT_DIR) \
		--out-name $(OUT_NAME).hapl \
		$(SLURM_OPTS) $(EXEC_FLAG)

## fasta: GBZ → augmented reference FASTA
fasta: $(FASTA)

$(FASTA): $(GBZ) scripts/fasta.sh
	$(SCRIPTS)/fasta.sh \
		--gbz $(GBZ) \
		--ref $(AUGREF) \
		--out-dir $(OUT_DIR) \
		--out-name $(OUT_NAME).fa.gz \
		$(SLURM_OPTS) $(EXEC_FLAG)

## deepvariant: BAM + FASTA → VCF via DeepVariant Docker
deepvariant: $(DV_VCF)

$(DV_VCF): $(BAM) $(FASTA) scripts/deepvariant.sh
	$(SCRIPTS)/deepvariant.sh \
		--bam $(BAM) \
		--ref $(FASTA) \
		--sample $(SAMPLE) \
		--out-dir $(OUT_DIR) \
		--out-name $(SAMPLE).deepvariant.vcf.gz \
		--dv-version $(DV_VERSION) \
		$(SLURM_OPTS) $(EXEC_FLAG)

## split-vcf: VCF → onref / nestedref / offref VCFs
split-vcf: $(OFFREF)

$(OFFREF): $(VCF) scripts/split-ref.sh
	cd $(OUT_DIR) && $(SCRIPTS)/split-ref.sh -v $(CURDIR)/$(VCF) -p $(AUGREF)

## length-hist: augref segment table → length histogram PNG
length-hist: $(LENGTH_HIST)

$(LENGTH_HIST): $(AUGREF_SEGS) scripts/offref-length-hist.R
	Rscript $(SCRIPTS)/offref-length-hist.R \
		$(LENGTH_HIST) $(AUGREF_SEGS) TRUE

## plots: deconstruct VCF → off-reference density ideogram PNG
plots: $(PLOT)

$(PLOT): $(VCF) $(AUGREF_SEGS) scripts/chrom-density-segs.R
	Rscript $(SCRIPTS)/chrom-density-segs.R \
		$(VCF) $(AUGREF_SEGS) $(PLOT) \
		"$(REF) Off-Reference Variant Density" \
		0 $(BED_FLAG) $(SCALE_TYPE) --ref $(REF) --offref

## call-plots: genotyped VCF → off-reference density ideogram PNG
call-plots: $(CALL_PLOT)

$(CALL_PLOT): $(AUGREF_SEGS) scripts/chrom-density-segs.R
	Rscript $(SCRIPTS)/chrom-density-segs.R \
		$(CALL_VCF) $(AUGREF_SEGS) $(CALL_PLOT) \
		"$(REF) Call Off-Reference Density ($(SAMPLE))" \
		0 $(BED_FLAG) $(SCALE_TYPE) --ref $(REF) --offref

## dv-plots: DeepVariant VCF → off-reference density ideogram PNG
dv-plots: $(DV_PLOT)

$(DV_PLOT): $(AUGREF_SEGS) scripts/chrom-density-segs.R
	Rscript $(SCRIPTS)/chrom-density-segs.R \
		$(DV_VCF) $(AUGREF_SEGS) $(DV_PLOT) \
		"$(REF) DeepVariant Off-Reference Density ($(SAMPLE))" \
		0 $(BED_FLAG) $(SCALE_TYPE) --ref $(REF) --offref

## test: run shellcheck and help-flag tests
test: test/test-pipeline.sh
	bash test/test-pipeline.sh

## clean: remove all generated outputs
clean:
	rm -rf $(OUT_DIR)

## help: list available targets
help:
	@echo "nested-variants pipeline"
	@echo ""
	@echo "Targets:"
	@echo "  all          Build everything (default: analysis)"
	@echo "  paths        VG → GBZ (augmented reference paths)"
	@echo "  deconstruct  GBZ → VCF (vg deconstruct)"
	@echo "  genotype     GBZ + reads → GAM → sample VCF (requires READS, HAPL, SAMPLE, MAP_GBZ)"
	@echo "  surject      GAM → sorted BAM (requires SAMPLE; uses augmented GBZ)"
	@echo "  haplotypes   GBZ → .hapl index for giraffe haplotype-aware mapping"
	@echo "  fasta        GBZ → augmented reference FASTA (requires paths output)"
	@echo "  deepvariant  BAM + FASTA → VCF via DeepVariant Docker (requires surject, fasta)"
	@echo "  split-vcf    VCF → onref / nestedref / offref VCFs"
	@echo "  length-hist        augref segments → length histogram"
	@echo "  plots              offref VCF → chromosome density ideogram"
	@echo "  call-plots         genotyped VCF → off-reference density ideogram (requires SAMPLE)"
	@echo "  dv-plots           DeepVariant VCF → off-reference density ideogram (requires SAMPLE)"
	@echo "  analysis     split-vcf + plots"
	@echo "  test         Run shellcheck and help-flag tests"
	@echo "  clean        Remove output directory"
	@echo ""
	@echo "Configuration:"
	@echo "  Edit config.mk or create config.local.mk (see config.local.mk.example)"
	@echo "  Override on command line:  make REF=CHM13 EXEC_MODE=slurm all"
	@echo ""
	@echo "Current settings:"
	@echo "  EXEC_MODE=$(EXEC_MODE)  REF=$(REF)  AUGREF=$(AUGREF)"
	@echo "  VG=$(VG)"
	@echo "  OUT_DIR=$(OUT_DIR)  OUT_NAME=$(OUT_NAME)"
	@echo "  MAP_GBZ=$(MAP_GBZ)  GIRAFFE_GBZ=$(GIRAFFE_GBZ)"
