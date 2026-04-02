#!/bin/bash
#Remove redundant sequences
#conda install bioconda::purge_haplotigs 
#github:https://github.com/skingan/purge_haplotigs_multiBAM

minimap2 -ax map-pb  -t 20 $breed.fa ../*bam.fasta.gz | samtools view -hF 256 - | samtools sort -@ 20 -m 1G -o aligned.bam -T tmp.ali
purge_haplotigs  readhist  -b aligned.bam  -g $breed.fa  -t 20
purge_haplotigs  contigcov  -i aligned.bam.gencov  -l 3  -m 150  -h 590
purge_haplotigs purge  -g $breed.fa  -c coverage_stats.csv -a 80

