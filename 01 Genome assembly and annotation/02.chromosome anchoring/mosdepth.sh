#!/bin/bash
#Align Illumina read pairs and calculate read depth
#conda install bioconda::mosdepth 
#conda install bioconda::samtools
#conda install bioconda::hifiasm 
#conda install bioconda::bwa
#ERR4351497
#github:https://github.com/brentp/mosdepth

bwa index -a bwtsw contig.fa
bwa mem -k 32 -w 10 -B 3 -O 11 -E 4 -t 60 contig.fa ERR4351497_1.fastq.gz ERR4351497_2.fastq.gz |samtools view -b -S - > ERR4351497_1.fastq.gz.bam
samtools sort -m 2000000000 -@ 60 -T tmp.bam -o ERR4351497_1.fastq.sort.bam ERR4351497_1.fastq.gz.bam
samtools index -@ 50 ERR4351497_1.fastq.sort.bam 
mosdepth -t 10 --fast-mode --by 500 sample-output ERR4351497_1.fastq.sort.bam
plot-dist.py *global.dist.txt
