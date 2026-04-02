#!/bin/bash
#verkko assembly
#conda install bioconda::verkko
#github:https://github.com/marbl/verkko

verkko -d verko --min-ont-length 50000  --hifi hifi*fa.gz  --nano ont*.fa.gz --threads 70 --unitig-abundance 4 --hic1 hic*_1.fq.gz --hic2 hic*_2.fq.gz  --haplo-divergence 0.12
