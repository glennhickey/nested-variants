# Convert GFA path names to PanSN format matching chr12.vg.
# Only modifies P-line path names (column 2).
#
# Reference path CHM13.0 is left as-is (generic name, not PanSN)
# so vg convert won't add a #0 fragment.
#
# Haplotype paths: SAMPLE.HAP -> SAMPLE#HAP#SAMPLE.HAP#0

/^P\t[^C]/ s/^P\t\([^.]*\)\.\([^\t]*\)\t/P\t\1#\2#\1.\2#0\t/
