# config.mk — Default configuration for local chr20 testing
# Override any variable on the command line or in config.local.mk

# Execution mode: "local" or "slurm"
EXEC_MODE   ?= local

# Reference name (used by vg paths -Q; deconstruct/call use augref_REF automatically)
REF         ?= GRCh38

# Input VG file(s) — glob patterns accepted (quoted)
VG          ?= data/chr20.vg

# Output directory and naming prefix
OUT_DIR     ?= output
OUT_NAME    ?= chr20.nested

# Minimum augmented-reference fragment length
MIN_AUGREF_LEN ?= 50

# Identity threshold for deconstruct clustering
CLUSTER     ?=

# Snarls file for deconstruct (optional)
SNARLS      ?=

# Star allele mode for deconstruct (set to 1 to enable)
STAR_ALLELE ?=

# Reads and sample for genotyping (optional)
READS       ?=
HAPL        ?=
SAMPLE      ?=
GAM         ?= data/chr20.sim.10.gam

# SLURM resource defaults (ignored when EXEC_MODE=local)
CPUS        ?= 8
MEM         ?= 200gb
TIME        ?= 16:00:00
PARTITION   ?= long

# Bin size for density plots (bp)
BIN_SIZE    ?= 1000000

# Minimum variant length for density plots
MIN_LENGTH  ?= 50

# BED file for reference gap overlay on plots (optional)
REFGAPS_BED ?=

# Scale type for density plots: log1p, sqrt, log, identity
SCALE_TYPE  ?= log1p
