#Derivation of PGC and Splicing Expression Levels from Transcriptomic Raw Data
sh 1. RNA_seq_quant.sh

#Mapping Alternative Splicing Events to Their Respective Genes
Rscript 2.1 Splicing.map_clusers_to_genes.R

#Normalization of Splicing Expression and Conversion to Percent Spliced In (PSI)
sh 2.1 Splicing.Normalize_and_convert_to_PSI.sh

#Merge the results of each chromosome (tissue_perind.counts.gz.qqnorm_chrN)
Rscript 2.1 Splicing.prepare_PSI_bed.R

#Normalization of Raw PGC Expression Using the TMM Method
Rscript 2.2 PGCs.calculate_TMM.R

#Calculation of PEER Latent Factors as Covariates for Subsequent QTL Analysis
Rscript 2.3 Covariate.PEER_factor_analysis.R

#Analysis of cis-molQTL and Independent cis-molQTL
sh 2.4 Cis-molQTL_analysis.sh