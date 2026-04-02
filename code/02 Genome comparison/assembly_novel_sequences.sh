#!/bin/bash

cat Gal4.fa Gal5.fa GRCg7b.fa GRCg7w.fa GGsu.fa > others_combined.fa
minimap2 -cx asm10 -t 200 ZJ.fa others_combined.fa > ZJ.paf
cut -f6,8,9 ZJ.paf | sort -k1,1 -k2,2n -t$'\t' > ZJ.covered.bed
bedtools merge -i ZJ.covered.bed > ZJ.covered_merged.bed
samtools faidx ZJ.fa
cut -f1,2 ZJ.fa.fai > ZJ.fa.sizes
bedtools sort -i ZJ.covered_merged.bed -g ZJ.fa.sizes > ZJ.covered_merged.sorted.bed
bedtools complement -i ZJ.covered_merged.sorted.bed -g ZJ.fa.sizes > ZJ.unique_regions.bed
awk '($3-$2) > 100' ZJ.unique_regions.bed > ZJ.unique_filtered.bed
bedtools getfasta -fi ZJ.fa -bed ZJ.unique_filtered.bed -fo ZJ.unique_sequences.fasta
grep -v "^>" ZJ.unique_sequences.fasta | tr -d '\n' | wc -c
