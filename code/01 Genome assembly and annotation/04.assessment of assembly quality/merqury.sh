#!/bin/bash
#Measure the completeness of the genome and the quality of the assembly
#conda install bioconda::merqury 
#github:https://github.com/marbl/merqury

#Evaluate the optimal K value

best_k.sh $genome_size

#Build meryl dbs

meryl k=20 count output read.meryl  $breed.hifi_reads.fastq.gz

merqury.sh ./read.meryl $breed.fa out_prefix







