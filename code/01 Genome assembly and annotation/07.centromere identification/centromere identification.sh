#!/bin/bash
#Identification of centromeres
#conda install bioconda::bwa 
#conda install bioconda::sambamba
#conda install bioconda::bedtools
#github:https://github.com/arq5x/bedtools2

ref=/home/ypzhang/1.chip/09.ref
raw=/home/ypzhang/1.chip/01.raw_data
res=/home/ypzhang/1.chip/06.chicken
TMPDIR=/home/ypzhang/1.chip/06.chicken

bwa mem -t 12  -k 50 -c 1000000 $ref/$breed.fa <(zcat $raw/DRR018430_1.fastq.gz) <(zcat $raw/DRR018430_2.fastq.gz) |  samtools sort -@ 8 -o $res/$breed.bwa.bam  - -O BAM

sambamba markdup --overflow-list-size 600000  --tmpdir=$TMPDIR -r $breed.bwa.bam $breed.bwa.rmdup.bam

samtools view -@ 12 -q 30 -F 2308 $breed.bwa.rmdup.bam -O BAM -o $breed.bwa.rmdup.2308.bam

samtools faidx $breed.fasta

bedtools makewindows -g $breed.fa.fai -w 10000 > $breed_genome.10k.bed

bedtools genomecov -ibam $breed.bwa.rmdup.2308.bam -d |  awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | bedtools map -c 4 -o sum -b - -a  $breed_genome.10k.bed > $breed.bwa.rmdup.2308.bam.10k




