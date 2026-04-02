#!/bin/bash
#validate anchor contigs
#conda install conda-forge::quartet
#github:https://github.com/aaranyue/quarTeT

python quartet_assemblymapper.py -r GGswu.fa -q contig.fa -c 50000 -p quarTeT -t 10  --plot