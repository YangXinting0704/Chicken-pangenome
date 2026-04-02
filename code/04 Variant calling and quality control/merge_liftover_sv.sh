#!/usr/bin/env bash

BED=GGsu_4.SV.toBY.bed
VCF=sorted.long50.norm.GGsu_4.vcf.gz
FASTA=BY.fa
OUT=GGsu_4.SV.toBY.vcf

# 临时文件
ALT_TABLE=ALT_table.tsv

# ---------- 1. 提取 ID->ALT 映射表 ----------
if [ ! -f "$ALT_TABLE" ]; then
    echo "Generating ID->ALT table..."
    bcftools query -f '%ID\t%ALT\n' $VCF > $ALT_TABLE
fi

# ---------- 2. 输出 VCF header ----------
bcftools view -h $VCF | sed 's/WChang/BY/' > $OUT

# ---------- 3. 并行处理 BED ----------
echo "Processing BED in parallel..."

cat $BED | parallel -j 80 --colsep '\t' '
    CHR={1}
    START={2}
    END={3}
    ID={4}
    POS=$((START + 1))
    LEN=$((END - START))

    # 取 REF
    REF=$(samtools faidx '"$FASTA"' ${CHR}:${POS}-$((${POS}+LEN-1)) | tail -n +2 | tr -d "\n")

    # 从 hash-table 查 ALT（避免 bcftools 调用）
    ALT=$(grep -P "^${ID}\t" '"$ALT_TABLE"' | cut -f2)

    # 输出 VCF 行
    printf "%s\t%d\t%s\t%s\t%s\t.\tPASS\tSVTYPE=NA\n" "$CHR" "$POS" "$ID" "$REF" "$ALT"
' >> $OUT

echo "Done! Output: $OUT"
