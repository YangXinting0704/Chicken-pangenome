###incorporate the GRCg6a assembly into pangenome graph
cactus-pangenome ./total/js ./total/seqfile.txt --outDir ./total --outName BYpan_GRCg6a --reference BY --gbz --gfa --vcf --mapCores 192 --indexCores 192 --maxCores 192

### View the genomic composition in the pangenome graph
halStats BYpan_GRCg6a.full.hal

### projection of functional annotation
halLiftover BYpan_GRCg6a.full.hal GRCg6a peak_GRCg6a.bed BY BY_peak_GRCg6a.bed
halLiftover BYpan_GRCg6a.full.hal GRCg7b peak_GRCg7b.bed BY BY_peak_GRCg7b.bed
