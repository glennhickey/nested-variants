# Pantree Comparison & Integration Notes

Reference: Salehi Nowbandegani et al., "Defining and cataloging variants in pangenome graphs" (bioRxiv 2025, doi:10.1101/2025.08.04.668502)

## 1. Direct Comparisons

### What's comparable and what isn't

Pantree analyzed the **HPRC v1.1 MC GRCh38** graph. We're running on **HPRC v2.1 MC CHM13**. So there are two confounds in any direct comparison: graph version (v1.1 vs v2.1) and reference initialization (GRCh38 vs CHM13). The paper itself notes that the CHM13 graph has 14.7% more variants, mostly in centromeres/subtelomeres/acrocentric arms. That said, there are still meaningful comparisons:

### a) Variant type distributions

Pantree reports a detailed breakdown (Fig 2a/b): 21.7M SNPs, 3.2M insertions, 3.5M deletions, 660K MNPs, etc., split by small (<50bp) vs large (>=50bp). Our pipeline produces `variant-types.png` and `vcf-stats.tsv` from `vg deconstruct`. We could compare the proportions even if absolute counts differ due to graph version. The key question: does the superbubble approach (`vg deconstruct`) produce a similar type distribution to the reference tree approach, or does it systematically miss certain categories?

### b) SNPs in segmental duplications

This is pantree's strongest claim against the superbubble approach: they find **46.1% more SNPs** in segdups (2.3M vs 1.6M) and 43.4% more in other difficult regions. We already have segdup annotation overlap. We could directly count SNPs in segdup-overlapping off-reference segments and compare with pantree's numbers for the v1.1 graph.

### c) Off-reference vs non-GRCh38 variants

Pantree identifies 3.5M non-GRCh38 variants (11.7%) -- variants whose reference allele is not on GRCh38. Our augmented reference approach explicitly models off-reference segments. These concepts are related but not identical: our off-reference variants are variants *called on* off-reference segments, while pantree's non-GRCh38 variants are those whose *reference tree allele* is off GRCh38. Still, we could compare the fraction of total variants that are "off-reference" in each approach.

### d) AF spectrum comparison

Pantree shows allele frequency distributions for non-GRCh38 vs GRCh38 variants (Fig 2e). Our pipeline already produces AF spectra. We could split by on-ref vs off-ref to see if the enrichment patterns match.

### Practical comparison steps

1. **Download pantree's VCF** from Zenodo (https://zenodo.org/records/15374896) -- they provide per-chromosome VCFs for both GRCh38 and CHM13 graphs
2. Run our `vcf-stats.R` on the pantree VCF to get comparable type/size/AF plots
3. Intersect pantree variants with our same annotation BEDs (segdups, repeats, etc.) to stratify by region
4. Compare pantree's CHM13 graph results directly with our v2.1 CHM13 results (same reference, different graph version)

## 2. Ideas to Incorporate

### High-value additions

#### a) Easy / Segdup / Hard region stratification

Pantree uses GIAB genomic stratification BEDs to partition variants into Easy, Segdup, and Hard regions. We already have segdup annotation; adding the GIAB "easy" regions BED would let us partition variants the same way. This is a lightweight addition -- just another `annot_*` config key pointing to the GIAB notinalldifficultregions BED. Then all our existing plots (summary, scatter, co-occurrence, SNP heatmaps) would automatically include it.

#### b) Non-reference variant flagging

Pantree's NR (non-reference) INFO field identifies variants whose reference allele is off the linear reference. Our `vg deconstruct` VCF already encodes this implicitly -- variants on off-reference contigs (`_alt` CHROMs) have non-reference alleles by definition. We could add a column to `segment-polymorphism.tsv` counting non-reference variants per segment, giving a direct analog to pantree's non-GRCh38 counts. This would let us answer: "of the variants in our off-reference segments, how many are entirely off-reference (both alleles) vs straddling on/off?"

#### c) Tandem repeat / homopolymer annotation for indels

Pantree annotates which indels are tandem repeat expansions/contractions (the hatched bars in Fig 2a/b). They find most small indels are repeat expansions, mostly in short tandem repeats and homopolymers. We could do something similar by cross-referencing our indels with the RepeatMasker classes we already have. Our repeat class breakdown plot already shows overlap fractions -- extending it to variant-type-specific breakdown would be natural.

### Medium-value additions

#### d) Multiallelic superbubble characterization

Pantree's analysis of triallelic superbubbles (Fig 3e-f: properly triallelic, overlapping, nested, interlocking) is interesting because it characterizes the structural complexity that `vg deconstruct` has to deal with. We could use the `.snarls` file (which our pipeline can produce via `vg snarls`) to classify superbubbles by allele count and nesting structure, then correlate with our annotation overlaps. This would answer: "are the off-reference segments with the most annotation complexity also the ones with the most multiallelic superbubbles?"

#### e) Per-variant region-of-origin annotation

Pantree annotates each variant with DR (distance from reference), telling you how far each node in a variant edge is from GRCh38. Our augmented reference segments already encode this implicitly (the segment's on-reference coordinate span). We could add a per-variant distance-from-reference metric to our polymorphism table.

### Lower-priority but interesting

#### f) vcfwave comparison

Pantree compares with vcfwave and finds 1.9M SNPs missed by vcfwave (mostly non-GRCh38). If we were to run vcfwave on our deconstruct VCF, we could see if the nested variants recovered by vcfwave's re-alignment match what we find via our off-reference segment analysis.

#### g) Complex region case studies

Pantree's HLA-A and RHD analyses (Fig 4) are compelling case studies. We could do analogous case studies for complex regions in the v2.1 CHM13 graph, using our annotation overlap data to identify candidate regions (segments with high repeat/segdup overlap and many variants).

## Priority Summary

- **Highest-impact comparison**: Download the pantree CHM13 VCF and run it through our existing stats pipeline for a side-by-side.
- **Highest-impact incorporation**: Add the GIAB Easy/Hard stratification -- essentially free with our existing annotation infrastructure and gives us the same regional breakdown pantree uses.
