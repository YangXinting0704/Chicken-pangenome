#!/bin/bash
#mapping contigs to the assembly
#conda install bioconda::mummer 
#github:https://github.com/mummer4/mummer

nucmer -t 40 --mum -c 200 -l 100 ref.fa contig.fa  -p scf
delta-filter -i 85 -l 100 -1 contig.delta > contig.flt.delta
show-coords -cdlqoTH contig.flt.delta > contig.flt.coords
mummerplot --fat --large --png -p contig.fa contig.flt.delta
