#!/bin/bash
#W chromosome anchoring
#conda install hcc::juicer
#conda install bioconda::3d-dna
#github:https://github.com/aidenlab/juicer
#github:https://github.com/aidenlab/3d-dna

path=`pwd`
awk -f wrap-fasta-sequence.awk p.fa > wrapped.fa
bwa index wrapped.fa
python generate_site_positions.py MboI wrapped  wrapped.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' wrapped_MboI.txt > wrapped.chrom.sizes

juicer.sh -t 100   -z ${path}/wrapped.fa -p ${path}/wrapped.chrom.sizes -y ${path}/wrapped_MboI.txt -s MboI -D /PUBLIC/software/bin/script/ok/juicer/juicer
run-asm-pipeline.sh -g 100 -m haploid --build-gapped-map --editor-repeat-coverage 4   -i 1000000  ${path}/wrapped.fa  ${path}/aligned/merged_nodups.txt
run-asm-pipeline-post-review.sh -r wrapped.review.assembly  -g 100  --build-gapped-map --sort-output wrapped.fasta aligned/merged_nodups.txt &> 3d.log &

