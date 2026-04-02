#!/bin/bash
#Making colorful identity heatmaps of genomic sequence
#mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy 
#github:https://github.com/mrvollger/StainedGlass https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/mrvollger/StainedGlass.html

samtools faidx $chr.fa $chr_region > $chr_region.fa

samtools faidx  $chr_region

snakemake -s /home/ypzhang/reference/6.StainedGlass/project-workdir/workflow/Snakefile --configfile /home/ypzhang/reference/6.StainedGlass/project-workdir/config/config.yaml --config sample=$chr_region fasta=/home/ypzhang/reference/6.StainedGlass/1.BY/$chr_region.fa --cores 24  make_figures -p

