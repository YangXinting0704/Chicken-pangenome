#!/bin/bash
#Statistical assembly results
#conda install bioconda::assembly-stats
#github:https://github.com/rjchallis/assembly-stats

assembly-stats $breed.fa > ./genome_N50/$breed.N50.stat

