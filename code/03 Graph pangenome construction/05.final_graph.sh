#!/bin/bash

tabix -p vcf giraffe.vcf.gz
vg construct -r BY.fa -v giraffe.vcf.gz -a > BYpanSV.vg
