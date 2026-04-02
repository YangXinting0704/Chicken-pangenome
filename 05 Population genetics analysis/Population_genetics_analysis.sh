#################################################################################################################################PCA
plink --chr-set 50 --allow-extra-chr --bfile selection.SV.final --indep-pairwise 500 50 0.2 --out selection.SV
plink --chr-set 50 --allow-extra-chr --bfile selection.SV.final --extract selection.SV.prune.in --pca 500 --out selection.SV
awk '{$1=$1}1' OFS='\t' selection.SV.eigenvec > selection.SV.eigenvec.tsv



#################################################################################################################################ADMIXTURE
plink --chr-set 50 --allow-extra-chr --bfile selection.SV.final --extract selection.SV.prune.in --make-bed --out selection.SV

for K in {1..10}; do
    ./admixture --cv selection.SV.bed $K | tee selection.SV.log${K}.out
done

grep -h "CV error" selection.SV.log*.out



#################################################################################################################################phylogenetic relationships
plink --chr-set 50 --allow-extra-chr --bfile selection.SV --recode vcf-iid --out selection.SV
python vcf2phylip.py -i selection.SV.vcf --fasta -o selection.SV.min4
FastTree -nt -gtr selection.SV.min4.fasta > selection.SV.tree



#################################################################################################################################Selection scan analysis
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop A.txt --weir-fst-pop B1.txt --out SV_A_B1
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop A.txt --weir-fst-pop B2.txt --out SV_A_B2
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop A.txt --weir-fst-pop B3.txt --out SV_A_B3
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop A.txt --weir-fst-pop B4.txt --out SV_A_B4
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop A.txt --weir-fst-pop B5.txt --out SV_A_B5
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop A.txt --weir-fst-pop B6.txt --out SV_A_B6
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop A.txt --weir-fst-pop B7.txt --out SV_A_B7
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop A.txt --weir-fst-pop C.txt --out SV_A_C
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop B1.txt --weir-fst-pop C.txt --out SV_B1_C
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop B2.txt --weir-fst-pop C.txt --out SV_B2_C
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop B3.txt --weir-fst-pop C.txt --out SV_B3_C
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop B4.txt --weir-fst-pop C.txt --out SV_B4_C
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop B5.txt --weir-fst-pop C.txt --out SV_B5_C
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop B6.txt --weir-fst-pop C.txt --out SV_B6_C
vcftools --vcf selection.SV.final.vcf --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop B7.txt --weir-fst-pop C.txt --out SV_B7_C

vcftools --vcf selection.SV.final.vcf --keep A.txt --out SV_A --window-pi 40000 --window-pi-step 20000
vcftools --vcf selection.SV.final.vcf --keep B1.txt --out SV_B1 --window-pi 40000 --window-pi-step 20000
vcftools --vcf selection.SV.final.vcf --keep B2.txt --out SV_B2 --window-pi 40000 --window-pi-step 20000
vcftools --vcf selection.SV.final.vcf --keep B3.txt --out SV_B3 --window-pi 40000 --window-pi-step 20000
vcftools --vcf selection.SV.final.vcf --keep B4.txt --out SV_B4 --window-pi 40000 --window-pi-step 20000
vcftools --vcf selection.SV.final.vcf --keep B5.txt --out SV_B5 --window-pi 40000 --window-pi-step 20000
vcftools --vcf selection.SV.final.vcf --keep B6.txt --out SV_B6 --window-pi 40000 --window-pi-step 20000
vcftools --vcf selection.SV.final.vcf --keep B7.txt --out SV_B7 --window-pi 40000 --window-pi-step 20000
vcftools --vcf selection.SV.final.vcf --keep C.txt --out SV_C --window-pi 40000 --window-pi-step 20000
