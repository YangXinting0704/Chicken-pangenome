#################################################################Excess heterozygosity filtering
plink --chr-set 50 --allow-extra-chr --vcf all_samples_vg_call.renamedID.vcf.gz --double-id --recode vcf-iid --out white --maf 0.01 --geno 0.1
bcftools view white.vcf -Oz -o white.vcf.gz
bcftools index white.vcf.gz

bcftools query -f '%CHROM\t%POS\t%ID\t[%GT\t]\n' white.vcf.gz | \
awk '{
  chrom=$1; pos=$2; id=$3;
  ref=0; het=0; hom=0; miss=0;
  for(i=4;i<=NF;i++){
    g=$i;
    if(g=="./."||g==".") miss++;
    else if(g=="0/0"||g=="0|0") ref++;
    else if(g=="1/1"||g=="1|1") hom++;
    else het++;
  }
  called = ref + het + hom;
  print chrom"\t"pos"\t"id"\t"called"\t"ref"\t"het"\t"hom"\t"miss
}' > white_per_site_genocounts.tsv

Rscript excessHet_test.R white_per_site_genocounts.tsv white_per_site_excessHet_test.tsv



#################################################################SV genotyping calibration versus long-read truth
python sv_matrix_qc_with_extra_metrics.py \
  --tgs TGS.vcf \
  --ngs NGS.vcf \
  --out per_variant_TGS_NGS.tsv



#################################################################SV genotyping accuracy based on pedigree-based concordance analysis
Rscript sv_family_concordance.R






















