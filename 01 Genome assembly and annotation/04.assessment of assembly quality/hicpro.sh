#!/bin/bash
##########################download hicpro#####################################################################
#wget -c https://github.com/nservant/HiC-Pro/archive/refs/tags/v3.1.0.tar.gz
#tar -zxvf v3.1.0.tar.gz
#######################install hicpro####################################################################
#conda env create -f /home/ypzhang/software/HiC-Pro-3.1.0/environment.yml -p /home/ypzhang/miniconda3/envs/hic
#PREFIX =/home/ypzhang/software/Hicpro
#BOWTIE2_PATH =/home/ypzhang/miniconda3/envs/hic/bin/bowtie2
#SAMTOOLS_PATH =/home/ypzhang/miniconda3/envs/hic/bin/samtools
#R_PATH =/home/ypzhang/miniconda3/envs/hic/bin/R
#PYTHON_PATH =/home/ypzhang/miniconda3/envs/hic/bin/python3.8
#CLUSTER_SYS =TORQUE
#make configure
#make install

bowtie2-build -f $breed.fa $breed
python /share/home/hnaas_liuj/software/HiC-Pro-3.1.0/bin/utils/digest_genome.py -r ^GATC -o $breed.bed  $breed.fa
samtools faidx $breed.fa
awk -v OFS='\t' '{print $1,$2}' $breed.fa.fai > $breed.fa.sizes
bwa index $breed.fa
/share/org/HNAASZHANGYP/hnaas_liuj/software/Hicpro/HiC-Pro_3.1.0/bin/HiC-Pro -i 02.reads/ -o  ./03.hicpro_out -c config-hicpro.txt 




