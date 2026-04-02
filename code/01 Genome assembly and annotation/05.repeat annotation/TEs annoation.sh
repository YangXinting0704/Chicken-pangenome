#!/bin/bash
#TE annotation
#conda install bioconda::repeatmodeler
#conda install bioconda::repeatmasker
#github:https://github.com/Dfam-consortium/FamDB https://github.com/Dfam-consortium/RepeatModeler https://github.com/Dfam-consortium/RepeatMasker/

sample=genome
mkdir result
famdb.py -i Libraries/RepeatMaskerLib.h5 families -f embl  -a -d Aves  > Aves.embl
buildRMLibFromEMBL.pl ~/Aves.embl > Aves.embl.fa
BuildDatabase -name ${sample}db -engine rmblast ${sample}.fa
RepeatModeler -database ${sample}db -threads 20 -LTRStruct
cat ${sample}db-families.fa Aves.embl.fa >TE.fa
RepeatMasker -nolow -e rmblast -pa 20 -norna -gff -lib TE.fa genome.fa -dir ./result  