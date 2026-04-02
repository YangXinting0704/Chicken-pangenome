######################################################################################(1) compute top three principal components (PCs)
plink --chr-set 50 --allow-extra-chr --bfile ALLB --indep-pairwise 50 5 0.2 --out ALLB
plink --chr-set 50 --allow-extra-chr --bfile ALLB --extract ALLB.prune.in --pca 3 --out ALLB



######################################################################################(2) estimate the effective number of independent markers (M_eff)

plink --chr-set 50 --allow-extra-chr --bfile order --recode vcf --out order
bcftools view order.vcf -Oz -o order.vcf.gz
bcftools index order.vcf.gz

R
library(SNPRelate)
snpgdsVCF2GDS("order.vcf.gz", "geno.gds", method="biallelic.only")
snpgdsOpen("geno.gds")

library(SNPRelate)
gds <- snpgdsOpen("geno.gds")

pruned <- snpgdsLDpruning(gds,
                           method="corr",
                           ld.threshold=0.2,
                           num.thread=80)

M_eff <- length(unlist(pruned))
M_eff



######################################################################################(3) Construction of LOCO kinship
sort -k1,1 ALL.bim | awk '{print > ($1 ".txt")}'

for chr in {1..38} W Z
do
    plink --chr-set 50 --allow-extra-chr \
          --bfile ALLB \
          --extract ${chr}.txt \
          --make-bed \
          --out chr${chr}

    plink --chr-set 50 --allow-extra-chr \
          --bfile ALLB \
          --exclude ${chr}.txt \
          --make-bed \
          --out nochr${chr}

    gemma -bfile nochr${chr} -gk 1 -o kinship_nochr${chr}
done



######################################################################################(4) Perform chromosome-specific association and then merge them
for chr in {1..38}; do
    gemma -bfile chr${chr} -k output/kinship_nochr${chr}.cXX.txt -lmm 1 -miss 0.1 -maf 0.01 -c COV.SEX.txt -o chr${chr}.trait01
done

for chr in W Z; do
    gemma -bfile chr${chr} -k output/kinship_nochr${chr}.cXX.txt -lmm 1 -miss 0.1 -maf 0.01 -c COV.SEX.txt -o chr${chr}.trait01
done

### Merge chromosome-specific association result
sed -i 's/\r//' Merging.suffix.sh
chmod +x Merging.suffix.sh
./Merging.suffix.sh



######################################################################################(5) genomic control factor (λ_GC) + Quantile–quantile (QQ) plots
ls *.assoc.txt | parallel -j 10 Rscript qq_lambda_one.R {}

head -n 1 FINAL.*.lambdaGC.txt | head -1 > all.lambdaGC.txt
cat FINAL.*.lambdaGC.txt | grep -v trait >> all.lambdaGC.txt



######################################################################################(6) LD-clumping
sed -i 's/\r//' run_clump_parallel.sh
chmod +x run_clump_parallel.sh
./run_clump_parallel.sh



######################################################################################(7) cross-trait FDR with Benjamini-Hochberg (BH)
R
dat <- read.table("ALL.leadQTL.txt", header = TRUE, stringsAsFactors = FALSE)
p <- dat$P
keep <- !is.na(p) & p > 0 & p < 1
dat <- dat[keep, ]
# BH-FDR
dat$FDR_BH <- p.adjust(dat$P, method = "BH")
# Bonferroni
dat$FDR_Bonf <- p.adjust(dat$P, method = "bonferroni")
# OUTPUT
write.table(dat,
            "ALL.leadQTL.withFDR.txt",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

