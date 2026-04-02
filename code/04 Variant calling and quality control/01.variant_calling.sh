#######################################################################################################Build an index
vg index BYpanSV.vg -x BYpanSV.xg -L
vg snarls BYpanSV.xg > BYpanSV.snarls
vg gbwt -x BYpanSV.xg -v giraffe.vcf.gz --force-phasing --ignore-missing -o BYpanSV.haps.gbwt
vg gbwt -x BYpanSV.xg -E -o BYpanSV.paths.gbwt
vg gbwt -m BYpanSV.haps.gbwt BYpanSV.paths.gbwt -o BYpanSV.combined.gbwt
vg gbwt -x BYpanSV.xg BYpanSV.combined.gbwt --gbz-format -g BYpanSV.giraffe.gbz
vg index -j BYpanSV.dist BYpanSV.xg
vg minimizer -d BYpanSV.dist -o BYpanSV.shortread.withzip.min -z BYpanSV.shortread.zipcodes BYpanSV.giraffe.gbz



#####################################################################################################SV calling based on short-reads
vg giraffe \
  -x BYpanSV.xg -Z BYpanSV.giraffe.gbz -z BYpanSV.shortread.zipcodes -m BYpanSV.shortread.withzip.min -d BYpanSV.dist\
  -f ./$ID_R1.fastq.gz -f ./$ID_R2.fastq.gz \
  --parameter-preset default \
  -N $ID \
  > ./V02/$ID.aln.gam

vg pack -x BYpanSV.xg -g ./V02/$ID.aln.gam -Q 10 -o ./V02/$ID.aln.pack

vg call BYpanSV.xg -k ./V02/$ID.aln.pack -v giraffe.vcf.gz -r BYpanSV.snarls -m 3 -s $ID > ./V02/$ID.genotypes.vcf



#####################################################################################################SNP and INDEL calling based on short-reads
# Surject graph alignments to linear reference
vg surject -x BYpanSV.xg -b aln.gam > aln.bam

# Sort and index
samtools sort -@ 48 -o aln.sorted.bam aln.bam
samtools index aln.sorted.bam

# Step 1: Deduplication
sentieon driver -t 48 -i aln.sorted.bam \
  --algo LocusCollector --fun score_info score.txt

sentieon driver -t 48 -i aln.sorted.bam \
  --algo Dedup --score_info score.txt \
  --metrics dedup_metrics.txt aln.dedup.bam

# Step 2: Base Quality Score Recalibration (BQSR)
sentieon driver -t 48 -r BY.fa -i aln.dedup.bam \
  --algo QualCal recal_data.table

sentieon driver -t 48 -r BY.fa -i aln.dedup.bam \
  --algo ApplyBQSR -q recal_data.table aln.recal.bam

# Step 3: Variant Calling (gVCF mode)
sentieon driver -t 48 -r BY.fa -i aln.recal.bam \
  --algo Haplotyper \
  --emit_mode gvcf \
  sample.g.vcf.gz



#####################################################################################################SV calling based on long-reads
vg autoindex --workflow lr-giraffe --prefix BYpanSV_LONG --gbz BYpanSV.giraffe.gbz
vg giraffe -b hifi -x BYpanSV.xg -Z BYpanSV.giraffe.gbz -z BYpanSV_LONG.longread.zipcodes -m BYpanSV_LONG.longread.withzip.min -d BYpanSV_LONG.dist -f ./HiFi/A79455.fastq.gz -N A79455 -R A79455 -p >hifi.A79455.gam
vg pack -x BYpanSV.xg -g hifi.A79455.gam -Q 10 -o hifi.A79455.pack
vg call BYpanSV.xg -k hifi.A79455.pack -v giraffe.vcf.gz -r BYpanSV.snarls -m 3 -s A79455 > hifi.A79455.vcf


