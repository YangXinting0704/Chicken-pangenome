#!/bin/bash

vg autoindex --workflow lr-giraffe --prefix pangenome_LONG --gbz pangenome.gbz
vg giraffe -b hifi -x pangenome.xg -Z pangenome.gbz -z pangenome_LONG.longread.zipcodes -m pangenome_LONG.longread.withzip.min -d pangenome_LONG.dist -f ./HiFi/A79455.fastq.gz -N A79455 -R A79455 -p >hifi.A79455.gam
vg pack -x pangenome.xg -g hifi.A79455.gam -Q 10 -o hifi.A79455.pack
vg call pangenome.xg -k hifi.A79455.pack -r BYpanSV.snarls -s A79455 -a > hifi.A79455.vcf
