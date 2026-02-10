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
CALL_VCF  := $(OUT_DIR)/$(SAMPLE).vcf.gz
CALL_PLOT := $(OUT_DIR)/$(SAMPLE).call-density.png

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

.PHONY: all analysis paths deconstruct genotype surject split-vcf plots call-plots test clean help

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

## split-vcf: VCF → onref / nestedref / offref VCFs
split-vcf: $(OFFREF)

$(OFFREF): $(VCF) scripts/split-ref.sh
	cd $(OUT_DIR) && $(SCRIPTS)/split-ref.sh -v $(CURDIR)/$(VCF) -p $(AUGREF)

## plots: offref VCF → density ideogram PNG
plots: $(PLOT)

$(PLOT): $(OFFREF) scripts/chrom-density-vcf.R
	Rscript $(SCRIPTS)/chrom-density-vcf.R \
		$(OFFREF) $(PLOT) \
		"$(REF) Off-Reference Variant Density (>=$(MIN_LENGTH)bp)" \
		$(MIN_LENGTH) $(BED_FLAG) $(SCALE_TYPE)

## call-plots: genotyped VCF → density ideogram PNG
call-plots: $(CALL_PLOT)

$(CALL_PLOT): scripts/chrom-density-call.R
	Rscript $(SCRIPTS)/chrom-density-call.R \
		$(CALL_VCF) $(CALL_PLOT) \
		"$(REF) Genotyped Variant Density ($(SAMPLE))" \
		0 $(BED_FLAG) $(SCALE_TYPE) --ref $(REF)

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
	@echo "  split-vcf    VCF → onref / nestedref / offref VCFs"
	@echo "  plots        offref VCF → chromosome density ideogram"
	@echo "  call-plots   genotyped VCF → chromosome density ideogram (requires SAMPLE)"
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
