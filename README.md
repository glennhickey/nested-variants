# nested-variants

## Finding the off-reference sequence

(all commands run in `./construction')

Right now, the off-reference sequence (ie rGFA cover) is computed within `vg deconstruct`.  It is output in three inter-related files:

* Nested VCF (contains variant calls both on and off the chosen reference)
* Off-reference FASTA (contains all off-reference contigs used by the VCF)
* Off-reference TSV (essentially a BED-formatted version of the Fasta, with some additional fields linking back to reference intervals)

To create these files, `vg deconstruct` must be run on `.vg` (and not `.gbz`) files.  For the HPRC graphs, these are normally found in the `.chroms` subdirectory alongside the main output.  This repo contains a script to help with this:

```
./slurm-deconstruct.sh --vg /private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.chroms/*.vg --ref CHM13 -L ${L} --out-dir $(pwd) --out-name hprc-v2.0-mc-chm13 --cpus 8 &
./slurm-deconstruct.sh --vg /private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.chroms/*.vg --ref CHM13 -L ${L} --out-dir $(pwd) --out-name hprc-v2.0-mc-chm13 --cpus 8 &
./slurm-deconstruct.sh --vg /private/groups/cgl/hprc-graphs/hprc-v1.1-jul4/hprc-v2.0-mc-chm13/hprc-v1.1-mc-chm13.chroms/*.vg --ref CHM13 -L ${L} --out-dir $(pwd) --out-name hprc-v1.1-mc-chm13 --cpus 8 &
./slurm-deconstruct.sh --vg /private/groups/cgl/hprc-graphs/hprc-v1.1-jul4/hprc-v2.0-mc-chm13/hprc-v1.1-mc-chm13.chroms/*.vg --ref CHM13 -L ${L} --out-dir $(pwd) --out-name hprc-v1.1-mc-chm13 --cpus 8 &
wait
```

Minigraph can produce something similar to the TSV and FASTA output above, so we get the BEDs for comparison.  The `grep` commands filter out the reference contigs (to be consistent with above) and the `sed` mess is to accommodate the v1.1 data which has the native cactus prefixes in the minigraph files, as opposed to PANSN.

```
gfatools gfa2bed -s /private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.sv.gfa.gz | grep -v ^CHM13 > hprc-v2.0-mc-chm13.sv.offref.bed
gfatools gfa2bed -s /private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-grch38/hprc-v2.0-mc-grch38.sv.gfa.gz | grep -v ^GRCh38 > hprc-v2.0-mc-grch38.sv.offref.bed
gfatools gfa2bed -s /private/groups/cgl/hprc-graphs/hprc-v1.1-jul4/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.sv.gfa.gz | sed -E 's/id=([[:alnum:]]+)\.([0-9])\|/\1#\2#/g; s/id=([[:alnum:]]+)\|/\1#0#/g' | grep -v ^CHM13 > hprc-v1.1-mc-chm13.sv.offref.bed
gfatools gfa2bed -s /private/groups/cgl/hprc-graphs/hprc-v1.1-jul4/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.sv.gfa.gz | sed -E 's/id=([[:alnum:]]+)\.([0-9])\|/\1#\2#/g; s/id=([[:alnum:]]+)\|/\1#0#/g' | grep -v ^GRCh38 > hprc-v1.1-mc-grch38.sv.offref.bed

```

### Variant Identity Threshold

The `-L` option sets a threshold for merging similar SV alt alleles.  This helps simplify the output VCF.  I've been using `0.95` but it could be interesting to compare other values (I've also tried 75,90,99,100).

Todo: script to summarize results here

```
for L in 0.75 0.90 0.99 1.00; do
./slurm-deconstruct.sh --vg /private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.chroms/*.vg --ref CHM13 -L ${L} --out-dir $(pwd) --out-name hprc-v2.0-mc-chm13 --cpus 8 &
./slurm-deconstruct.sh --vg /private/groups/cgl/hprc-graphs/hprc-v2.0-feb28/hprc-v2.0-mc-chm13/hprc-v2.0-mc-chm13.chroms/*.vg --ref CHM13 -L ${L} --out-dir $(pwd) --out-name hprc-v2.0-mc-chm13 --cpus 8 &
./slurm-deconstruct.sh --vg /private/groups/cgl/hprc-graphs/hprc-v1.1-jul4/hprc-v2.0-mc-chm13/hprc-v1.1-mc-chm13.chroms/*.vg --ref CHM13 -L ${L} --out-dir $(pwd) --out-name hprc-v1.1-mc-chm13 --cpus 8 &
./slurm-deconstruct.sh --vg /private/groups/cgl/hprc-graphs/hprc-v1.1-jul4/hprc-v2.0-mc-chm13/hprc-v1.1-mc-chm13.chroms/*.vg --ref CHM13 -L ${L} --out-dir $(pwd) --out-name hprc-v1.1-mc-chm13 --cpus 8 &
wait
```

## Off-reference Variant Stats

(all commands run in `./variants')

### Size Distribution

The size distribution of the off-reference variants can be computed from the tsv files

```
./offref-length-hist.R grch38-offref-lengths.png ../construction/hprc-v1.1-mc-grch38.nested.95.fa.nesting.tsv ../construction/hprc-v2.0-mc-grch38.nested.95.fa.nesting.tsv ../construction/hprc-v1.1-mc-grch38.sv.offref.bed ../construction/hprc-v2.0-mc-grch38.sv.offref.bed  50 TRUE

./offref-length-hist.R chm13-offref-lengths.png ../construction/hprc-v1.1-mc-chm13.nested.95.fa.nesting.tsv ../construction/hprc-v2.0-mc-chm13.nested.95.fa.nesting.tsv ../construction/hprc-v1.1-mc-chm13.sv.offref.bed ../construction/hprc-v2.0-mc-chm13.sv.offref.bed  50 TRUE
```

