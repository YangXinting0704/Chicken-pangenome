#!/bin/bash
#hifiasm assembly
#conda install bioconda::hifiasm 
#github:https://github.com/chhylp123/hifiasm

hifiasm  --b-cov 1   --ul-cut 40000  --path-min 0.2 --path-max 0.9  -y 0.2 -x 0.9 --ul ont*.fa.gz --h1 hic*_1.fq.gz  --h2 hic*_2.fq.gz -t 80  hifi*.fa.gz
