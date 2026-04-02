#!/bin/bash
set -euo pipefail

REF=/path/to/BY.fa
VCF_DIR=vcfs
WORK=work
mkdir -p ${WORK}/{siteonly,norm,prefixed}
# 1) Compress all the inputs
for f in ${VCF_DIR}/*.vcf; do
  bgzip -c "$f" > "${f}.gz"
  tabix -p vcf "${f}.gz"
done || true

# 2) List all.vcf.gz inputs
ls ${VCF_DIR}/*.vcf.gz > ${WORK}/all_vcfs.list

# 3) One by one: Remove GT (site-only) -> Normalize (left-align & split multi-allelic) -> Filter length >=50
while read v; do
  base=$(basename "$v" .vcf.gz)
  echo "Processing $base"
  # remove genotypes
  bcftools view -G -Oz -o ${WORK}/siteonly/${base}.siteonly.vcf.gz "$v"
  tabix -p vcf ${WORK}/siteonly/${base}.siteonly.vcf.gz

  # normalize: left-align, split multiallelic
  bcftools norm -f ${REF} -m -any -Oz -o ${WORK}/norm/${base}.norm.vcf.gz ${WORK}/siteonly/${base}.siteonly.vcf.gz
  tabix -p vcf ${WORK}/norm/${base}.norm.vcf.gz

  # filter: keep only records where REF length>=50 OR ALT length>=50
  zcat ${WORK}/norm/${base}.norm.vcf.gz \
    | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {ref_len=length($4); alt_len=length($5); if(ref_len>=50 || alt_len>=50) print}' \
    | bgzip -c > ${WORK}/norm/${base}.norm.long50.vcf.gz
  tabix -p vcf ${WORK}/norm/${base}.norm.long50.vcf.gz

done < ${WORK}/all_vcfs.list