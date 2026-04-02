#!/bin/bash
set -e

tissues=("breast" "liver" "whole blood")
analyses=("eqtl" "sqtl")

for tissue in "${tissues[@]}"; do
    for analysis in "${analyses[@]}"; do
        
        if [ "$analysis" == "eqtl" ]; then
            # cis-eQTL analysis
            omiga --mode cis --genotype geno.imputed --phenotype "${tissue}.expr.bed" --covariates "${tissue}.cov.txt" --dcovar-name gender sex --calcu-variant-threshold --prefix "${tissue}" --output-dir "${tissue}_eqtl"
            
            # independent cis-eQTL analysis
            omiga --mode cis_independent --genotype geno.imputed --phenotype "${tissue}.expr.bed" --covariates "${tissue}.cov.txt" --dcovar-name gender sex --prefix "${tissue}" --cis-file "${tissue}.TPM.cis_qtl.txt.gz" --output-dir "${tissue}_eqtl"
        
        else
            # cis-sQTL analysis
            omiga --mode cis --genotype geno.imputed --phenotype "${tissue}.expr.bed" --covariates "${tissue}.cov.txt" --dcovar-name gender sex --pheno-group "${tissue}.splicing.group.txt" --calcu-variant-threshold --prefix "${tissue}" --output-dir "${tissue}_sqtl"
            
            # independent cis-sQTL analysis
            omiga --mode cis_independent --genotype geno.imputed --phenotype "${tissue}.expr.bed" --covariates "${tissue}.cov.txt" --dcovar-name gender sex --pheno-group "${tissue}.splicing.group.txt" --prefix "${tissue}" --cis-file "${tissue}.PSI.cis_qtl.txt.gz" --output-dir "${tissue}_sqtl"
        fi
    
    done
done