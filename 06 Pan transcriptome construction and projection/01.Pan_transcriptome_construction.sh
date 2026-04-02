#########Construction of linear pan-genome (PSVCP)
python3 ./psvcp_v1.01/1Genome_construct_Pangenome.py ./pantranscriptome ./psvcp/genome_list 



##########Extract the gene & transcript sequence from GFF3 + FASTA of each genome (to prepare for mapping)

GENOMES_DIR=./pantranscriptome
OUTSEQ_DIR=./mapping/sequence
mkdir -p $OUTSEQ_DIR

# loop over samples
for gff in $GENOMES_DIR/*.gff; do
  prefix=$(basename "$gff" .gff)
  fa="$GENOMES_DIR/${prefix}.fa"
  if [ ! -f "$fa" ]; then
    echo "Fasta not found for $prefix, skipping"
    continue
  fi
  echo "Processing $prefix"

  # 1) transcript cDNA (spliced) sequences using gffread
  # Requires gffread (from cufflinks/stringtie)
  mkdir -p $OUTSEQ_DIR/$prefix
  gffread "$gff" -g "$fa" -w "$OUTSEQ_DIR/$prefix/${prefix}_transcripts.cDNA.fa"
  # output: spliced transcript sequences (one record per transcript)

  # 2) gene genomic sequences (gene locus including introns)
  # Extract gene coordinates from GFF3, then bedtools getfasta
  # Produce BED of genes: chrom  start  end  gene_id  .  strand
  awk '$3=="gene"{ gene="."; 
       # try find ID in attributes
       match($0,/ID=([^;]+)/,a); if(a[1]!="") gene=a[1];
       print $1"\t"($4-1)"\t"$5"\t"gene"\t.\t"$7 }' "$gff" > $OUTSEQ_DIR/$prefix/${prefix}_genes.bed

  bedtools getfasta -fi "$fa" -bed $OUTSEQ_DIR/$prefix/${prefix}_genes.bed -name -s -fo $OUTSEQ_DIR/$prefix/${prefix}_genes.genomic.fa

  # 3) Optional: sanity check counts
  echo -n "  transcripts: "
  grep -c "^>" $OUTSEQ_DIR/$prefix/${prefix}_transcripts.cDNA.fa || true
  echo -n "  genes fasta seqs: "
  grep -c "^>" $OUTSEQ_DIR/$prefix/${prefix}_genes.genomic.fa || true
done



###############Map all gene (including intron) sequences to the pan-genome (identify gene locus)
OUTSEQ_DIR=./mapping/sequence
PAN_FA=./psvcp/pan_genome.fa
MAP_OUT_DIR=./mapping/gene_level
mkdir -p $MAP_OUT_DIR

# loop per sample (parallelize as desired)
for d in $OUTSEQ_DIR/*; do
  prefix=$(basename "$d")
  gene_fa="$d/${prefix}_genes.genomic.fa"
  if [ ! -f "$gene_fa" ]; then continue; fi

  outsam=$MAP_OUT_DIR/${prefix}_genes_to_pan.sam
  minimap2 -t 300 -G 500000 -L --secondary=no --MD -ax splice:hq -uf "$PAN_FA" "$gene_fa" > "$outsam"
  # convert to sorted bam
  prefix_noext=${outsam%.*}
  samtools view -Sb -F 2048 "$outsam" > ${prefix_noext}.bam
  samtools sort -@4 ${prefix_noext}.bam -o ${prefix_noext}_sorted.bam
  samtools index ${prefix_noext}_sorted.bam
done



###############Precisely map transcripts (cDNA) to the local region of each gene (to avoid cross-site mismatch)
MAP_OUT_DIR=./mapping/gene_level
MAP_TRANS_DIR=./mapping/transcript_level
PADDING=500      # Expand both ends of the locus by 500bp to accommodate the UTR variation
mkdir -p $MAP_TRANS_DIR

# example for one sample prefix
prefix=D
gene_bam=$MAP_OUT_DIR/${prefix}_genes_to_pan_sorted.bam
trans_fa=$OUTSEQ_DIR/$prefix/${prefix}_transcripts.cDNA.fa

# 1)  BAM->bed12
bedtools bamtobed -bed12 -i "$MAP_OUT_DIR/${prefix}_genes_to_pan_sorted.bam" > "$MAP_TRANS_DIR/${prefix}_genes_to_pan.bed12"

# 2) Aggregate by gene id（min start, max end），output chr start end gene_id 0 strand
awk 'BEGIN{OFS="\t"}{
  chr=$1; start=$2; end=$3; name=$4; strand=$6;
  if(!(name in min) || start < min[name]) min[name]=start;
  if(!(name in max) || end > max[name]) max[name]=end;
  chrname[name]=chr; strandname[name]=strand;
}
END{
  for(name in min) print chrname[name], min[name], max[name], name, 0, strandname[name];
}' "$MAP_TRANS_DIR/${prefix}_genes_to_pan.bed12" | sort -k1,1 -k2,2n > "$MAP_TRANS_DIR/${prefix}_genes_to_pan.coords.bed"

# 3) pad loci safely
samtools faidx "$PAN_FA"
bedtools slop -i "$MAP_TRANS_DIR/${prefix}_genes_to_pan.coords.bed" -g "${PAN_FA}.fai" -b "$PADDING" > "$MAP_TRANS_DIR/${prefix}_genes_to_pan.padded.bed"

# 4) extract loci fasta (strand-aware name and sequence)
bedtools getfasta -fi "$PAN_FA" -bed "$MAP_TRANS_DIR/${prefix}_genes_to_pan.padded.bed" -name -s -fo "$MAP_TRANS_DIR/${prefix}_pan_loci.fa"

# 5) map transcripts to loci
minimap2 -t 256 -ax splice:hq -uf --secondary=no --MD "$MAP_TRANS_DIR/${prefix}_pan_loci.fa" "$OUTSEQ_DIR/$prefix/${prefix}_transcripts.cDNA.fa" > "$MAP_TRANS_DIR/${prefix}_trans_to_loci.sam"

# 6) convert, sort, index
samtools view -Sb -F 2048 "$MAP_TRANS_DIR/${prefix}_trans_to_loci.sam" > "$MAP_TRANS_DIR/${prefix}_trans_to_loci.bam"
samtools sort -@4 "$MAP_TRANS_DIR/${prefix}_trans_to_loci.bam" -o "$MAP_TRANS_DIR/${prefix}_trans_to_loci_sorted.bam"
samtools index "$MAP_TRANS_DIR/${prefix}_trans_to_loci_sorted.bam"



##########################Generate a unified GTF from all transcript→pan mapping results (merge all samples)

MAP_TRANS_DIR=./mapping/transcript_level
MERGE_MAPPING_DIR=./mapping/merged
mkdir -p $MERGE_MAPPING_DIR

# create a combined bed12 of all transcripts mapped (one file)
for bam in $MAP_TRANS_DIR/*_trans_to_loci_sorted.bam; do
  prefix=$(basename "$bam" _trans_to_loci_sorted.bam)
  bedtools bamtobed -bed12 -i $bam > $MERGE_MAPPING_DIR/${prefix}_trans.bed
done

# concatenate all sample beds into a single bed
cat $MERGE_MAPPING_DIR/*_trans.bed > $MERGE_MAPPING_DIR/all_samples_trans.bed

# convert bed12 -> genePred -> gtf using UCSC tools
bedToGenePred $MERGE_MAPPING_DIR/all_samples_trans.bed $MERGE_MAPPING_DIR/all_samples_trans.genepred
genePredToGtf file $MERGE_MAPPING_DIR/all_samples_trans.genepred $MERGE_MAPPING_DIR/all_samples_trans.gtf

python3 normalize_transcripts_to_pan.py all_samples_trans.gtf > transcripts_to_pan.normalized.gtf



#########################The core principle of generating PanBaRT (merging transcripts and creating a look-up table) using R is "structure priority, position priority, and rule merging": Transcripts with the same set of introns (representing the same splicing pattern) are prioritized for merging; single-exon transcripts are merged based on position overlap; special determination and handling are performed for small overlaps or nested cases to avoid incorrect merging.

Rscript PanBaRT_build_from_gtf.R transcripts_to_pan.normalized.gtf ChickenPan_ 0.05

