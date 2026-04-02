#!/bin/bash
#Methylation analysis
#conda install bioconda::nanopolish 
#github:https://github.com/jts/nanopolish


raw=/home/ypzhang/reference/5.genome_assembly/2.ONT/1.BY/01.raw_data
ref=/home/ypzhang/reference/2.genome
bam=/home/ypzhang/reference/5.genome_assembly/2.ONT/1.BY/02.bam

minimap2  -K 200M  --secondary=no -t 20 -ax  map-ont ${ref}/$breed.fa ${raw}/$breed.ont.fastq.gz | samtools view --threads 10 -T ${ref}/$breed.fa -bS | samtools sort --threads 8 -m 1G -o ${bam}/$breed.ont.sort.bam

samtools index  -@ 10  ${bam}/$breed.ont.sort.bam

#methylation calling

nanopolish call-methylation -t 20 -r  ${raw}/$breed.ont.fastq.gz  -b ${bam}/$breed.ont.sort.bam -g ${ref}/$breed.fa   > $breed.ont.methylation_calls.tsv

calculate_methylation_frequency.py $breed.ont.methylation_calls.tsv > $breed.ont.methylation_frequency.tsv

# merge methylation frequency from multiple Nanopore cells

cat *methylation_frequency.tsv   |    grep -v ^chro | awk '{a[$1"__"$2]+=$5;b[$1"__"$2]+=$6}END{for(i in a){print i"\t"a[i]"\t"b[i]}}' | sed 's/__/ /' | awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4}'  > merged.methy

# get methylation level in 500 bp windows

bedtools map -c 4,5 -o sum -a <(bedtools sort -i $breed.fa.fai.500) -b $breed.merged.methy > $breed.merged.methy.sort.500-sum 
