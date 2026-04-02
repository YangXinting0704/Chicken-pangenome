#!/bin/bash
#Identify the telomeres
#conda install conda-forge::quartet
#github:https://github.com/aaranyue/quarTeT

python quartet_teloexplorer.py -i $breed.fa -c animal -m 100 -p $breed