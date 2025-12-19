#!/bin/bash

# Generate all chromosome density plots from VCF files only (using reference coordinates from INFO field)

echo "Starting VCF plot generation..."

./chrom-density-vcf.R off-ref-calls/hprc-v2.0-mc-chm13.nested.95.offref.vcf.gz chrom-density-vcf-v2-chm13-50bp.png "CHM13 v2.0 Variant Density (>=50bp)" 50 hprc-v2.0-mc-chm13.refgaps.bed
echo "Completed 1/16"

./chrom-density-vcf.R off-ref-calls/hprc-v2.0-mc-chm13.nested.95.offref.vcf.gz chrom-density-vcf-v2-chm13-1000bp.png "CHM13 v2.0 Variant Density (>=1000bp)" 1000 hprc-v2.0-mc-chm13.refgaps.bed
echo "Completed 2/16"

./chrom-density-vcf.R off-ref-calls/hprc-v2.0-mc-chm13.nested.95.offref.vcf.gz chrom-density-vcf-v2-chm13-10000bp.png "CHM13 v2.0 Variant Density (>=10000bp)" 10000 hprc-v2.0-mc-chm13.refgaps.bed
echo "Completed 3/16"

./chrom-density-vcf.R off-ref-calls/hprc-v2.0-mc-chm13.nested.95.offref.vcf.gz chrom-density-vcf-v2-chm13-100000bp.png "CHM13 v2.0 Variant Density (>=100000bp)" 100000 hprc-v2.0-mc-chm13.refgaps.bed
echo "Completed 4/16"

./chrom-density-vcf.R off-ref-calls/hprc-v2.0-mc-grch38.nested.95.offref.vcf.gz chrom-density-vcf-v2-grch38-50bp.png "GRCh38 v2.0 Variant Density (>=50bp)" 50 hprc-v2.0-mc-grch38.refgaps.bed
echo "Completed 5/16"

./chrom-density-vcf.R off-ref-calls/hprc-v2.0-mc-grch38.nested.95.offref.vcf.gz chrom-density-vcf-v2-grch38-1000bp.png "GRCh38 v2.0 Variant Density (>=1000bp)" 1000 hprc-v2.0-mc-grch38.refgaps.bed
echo "Completed 6/16"

./chrom-density-vcf.R off-ref-calls/hprc-v2.0-mc-grch38.nested.95.offref.vcf.gz chrom-density-vcf-v2-grch38-10000bp.png "GRCh38 v2.0 Variant Density (>=10000bp)" 10000 hprc-v2.0-mc-grch38.refgaps.bed
echo "Completed 7/16"

./chrom-density-vcf.R off-ref-calls/hprc-v2.0-mc-grch38.nested.95.offref.vcf.gz chrom-density-vcf-v2-grch38-100000bp.png "GRCh38 v2.0 Variant Density (>=100000bp)" 100000 hprc-v2.0-mc-grch38.refgaps.bed
echo "Completed 8/16"

./chrom-density-vcf.R off-ref-calls/hprc-v1.1-mc-chm13.nested.95.offref.vcf.gz chrom-density-vcf-v1-chm13-50bp.png "CHM13 v1.1 Variant Density (>=50bp)" 50 hprc-v2.0-mc-chm13.refgaps.bed
echo "Completed 9/16"

./chrom-density-vcf.R off-ref-calls/hprc-v1.1-mc-chm13.nested.95.offref.vcf.gz chrom-density-vcf-v1-chm13-1000bp.png "CHM13 v1.1 Variant Density (>=1000bp)" 1000 hprc-v2.0-mc-chm13.refgaps.bed
echo "Completed 10/16"

./chrom-density-vcf.R off-ref-calls/hprc-v1.1-mc-chm13.nested.95.offref.vcf.gz chrom-density-vcf-v1-chm13-10000bp.png "CHM13 v1.1 Variant Density (>=10000bp)" 10000 hprc-v2.0-mc-chm13.refgaps.bed
echo "Completed 11/16"

./chrom-density-vcf.R off-ref-calls/hprc-v1.1-mc-chm13.nested.95.offref.vcf.gz chrom-density-vcf-v1-chm13-100000bp.png "CHM13 v1.1 Variant Density (>=100000bp)" 100000 hprc-v2.0-mc-chm13.refgaps.bed
echo "Completed 12/16"

./chrom-density-vcf.R off-ref-calls/hprc-v1.1-mc-grch38.nested.95.offref.vcf.gz chrom-density-vcf-v1-grch38-50bp.png "GRCh38 v1.1 Variant Density (>=50bp)" 50 hprc-v2.0-mc-grch38.refgaps.bed
echo "Completed 13/16"

./chrom-density-vcf.R off-ref-calls/hprc-v1.1-mc-grch38.nested.95.offref.vcf.gz chrom-density-vcf-v1-grch38-1000bp.png "GRCh38 v1.1 Variant Density (>=1000bp)" 1000 hprc-v2.0-mc-grch38.refgaps.bed
echo "Completed 14/16"

./chrom-density-vcf.R off-ref-calls/hprc-v1.1-mc-grch38.nested.95.offref.vcf.gz chrom-density-vcf-v1-grch38-10000bp.png "GRCh38 v1.1 Variant Density (>=10000bp)" 10000 hprc-v2.0-mc-grch38.refgaps.bed
echo "Completed 15/16"

./chrom-density-vcf.R off-ref-calls/hprc-v1.1-mc-grch38.nested.95.offref.vcf.gz chrom-density-vcf-v1-grch38-100000bp.png "GRCh38 v1.1 Variant Density (>=100000bp)" 100000 hprc-v2.0-mc-grch38.refgaps.bed
echo "Completed 16/16"

echo "All VCF plots generated successfully!"
