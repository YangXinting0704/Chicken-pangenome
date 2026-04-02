#!/usr/bin/env bash

BED=GGsu_4.SV.toBY.bed
VCF=sorted.long50.norm.GGsu_4.vcf.gz
FASTA=BY.fa
OUT=GGsu_4.SV.toBY.vcf

# temporary files
ALT_TABLE=ALT_table.tsv

# ---------- 1. Extract the ID->ALT mapping table ----------
if [ ! -f "$ALT_TABLE" ]; then
    echo "Generating ID->ALT table..."
    bcftools query -f '%ID\t%ALT\n' $VCF > $ALT_TABLE
fi

# ---------- 2. Output the VCF header ----------
bcftools view -h $VCF | sed 's/WChang/BY/' > $OUT

# ---------- 3. Parallel processing BED ----------
echo "Processing BED in parallel..."

cat $BED | parallel -j 80 --colsep '\t' '
    CHR={1}
    START={2}
    END={3}
    ID={4}
    POS=$((START + 1))
    LEN=$((END - START))

    # extract REF
    REF=$(samtools faidx '"$FASTA"' ${CHR}:${POS}-$((${POS}+LEN-1)) | tail -n +2 | tr -d "\n")

    # Look up ALT from hash-table (avoid calling bcftools)
    ALT=$(grep -P "^${ID}\t" '"$ALT_TABLE"' | cut -f2)

    # Output the VCF line
    printf "%s\t%d\t%s\t%s\t%s\t.\tPASS\tSVTYPE=NA\n" "$CHR" "$POS" "$ID" "$REF" "$ALT"
' >> $OUT

echo "Done! Output: $OUT"
