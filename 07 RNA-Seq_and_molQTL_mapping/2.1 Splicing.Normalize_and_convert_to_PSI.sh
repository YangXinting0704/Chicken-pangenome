#!/bin/bash
set -e

tissues=("breast" "liver" "whole blood")

for tissue in "${tissues[@]}"; do
    python prepare_phenotype_table.py "${tissue}_perind.counts.gz" -p 0
done