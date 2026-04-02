#!/usr/bin/env bash

cactus-pangenome ./BY_WC/js ./BY_WC/seqfile.txt --outDir ./BY_WC --outName BYpan_WC --reference BY --gbz --gfa --vcf --mapCores 292 --indexCores 292 --maxCores 292

halStats BYpan_WC.full.hal | grep Genomes

bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' WChang.vcf.gz > WChang.SV.bed

halLiftover BYpan_WC.full.hal WChang WChang.SV.bed BY WChang.SV.toBY.bed

sed -i 's/\r//' merge_liftover_sv.sh
chmod +x merge_liftover_sv.sh
./merge_liftover_sv.sh

sed -i 's/\r//' vcf_list.txt
SURVIVOR merge vcf_list.txt 0.5 1 1 0 0 50 merged.vcf
plink --chr-set 40 --allow-extra-chr --vcf merged.vcf --recode vcf-iid --out plink.merged
