#!/bin/bash
#Assess the integrity of the genome
#conda install bioconda::busco 
#github:https://github.com/metashot/busco

busco -i $breed.fa -c 8 -o busco -m geno -l home/ypzhang/reference/9.busco_data/aves_odb10 --offline



