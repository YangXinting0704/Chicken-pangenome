#!/bin/bash
#Fill the gaps
#conda install conda-forge::quartet
#github:https://github.com/aaranyue/quarTeT

python quartet_gapfiller.py -d assembly.fa -g contigs.fa -f 5000 -l 1000 -m 1000000 -i 40 -p fill -t 20 --fillonly --overwrite 