################################################################# Disassemble multiple alleles and left-align them (reference genome ref.fa)
bcftools index SNP.vcf.gz
bcftools norm -m -any -f BY.fa -Oz -o norm.SNP.vcf.gz SNP.vcf.gz -W
bcftools norm -m -any -f BY.fa -Oz -o norm.indel.vcf.gz indel.vcf.gz -W


################################################################# site-level hard filter with bcftools
bcftools filter -e 'INFO/QD < 2.0 || INFO/FS > 60 || INFO/MQ < 40 || INFO/MQRankSum < -12.5 || INFO/ReadPosRankSum < -8' -Oz -o hardfiltered.norm.SNP.vcf.gz norm.SNP.vcf.gz
bcftools filter -e 'INFO/QD < 2.0 || INFO/FS > 200 || INFO/ReadPosRankSum < -20' -Oz -o hardfiltered.norm.indel.vcf.gz norm.indel.vcf.gz



################################################################# MAF and CallRate filtering
bcftools +fill-tags hardfiltered.norm.SNP.vcf.gz -Oz -o tag.hardfiltered.norm.SNP.vcf.gz -- -t MAF,F_MISSING
bcftools view -i 'REF!="*" && ALT!="*" && MAF>=0.01 && F_MISSING<=0.1' -Oz -o final.SNP.vcf.gz tag.hardfiltered.norm.SNP.vcf.gz

bcftools +fill-tags hardfiltered.norm.indel.vcf.gz -Oz -o tag.hardfiltered.norm.indel.vcf.gz -- -t MAF,F_MISSING
bcftools view -i 'REF!="*" && ALT!="*" && MAF>=0.01 && F_MISSING<=0.1' -Oz -o final.indel.vcf.gz tag.hardfiltered.norm.indel.vcf.gz



################################################################# excess heterozygosity filtering
Rscript excessHet_genomewide.R \
  final_INDELSNP_per_site_genocounts.tsv \
  final_INDELSNP_excessHet.signif.tsv
