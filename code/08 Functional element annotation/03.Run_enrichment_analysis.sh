#!/bin/bash

# input file
QTL="eQTL.SV.bed"
BG="background.SV.bed"
ANNOT="annotation.bed"

# output file
OUT="annotation.eQTL.SV.txt"
TMP="annotation.tmp.txt"
echo -e "annotation\tp-value\todds ratio\tFold Enrichment\tconfidence interval_1\tconfidence interval_2" > $OUT
> $TMP  # Clear the temporary files

# Get all annotation types
TYPES=$(cut -f4 $ANNOT | sort | uniq)

# Make statistics for each type of annotation separately
for TYPE in $TYPES; do
    echo "Processing annotation: $TYPE"

    # Extract the annotation area of this category
    awk -v t="$TYPE" '$4==t' $ANNOT > temp.annot.bed

    # A: QTL in annotation
    A=$(bedtools intersect -a $QTL -b temp.annot.bed -u | wc -l)

    # B: QTL not in annotation
    B=$(bedtools intersect -a $QTL -b temp.annot.bed -v | wc -l)

    # C: BG in annotation
    C=$(bedtools intersect -a $BG -b temp.annot.bed -u | wc -l)

    # D: BG not in annotation
    D=$(bedtools intersect -a $BG -b temp.annot.bed -v | wc -l)

    # Generate R scripts
    cat <<EOF > temp_fisher.R
a <- $A
b <- $B
c <- $C
d <- $D
mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
              dimnames = list(c("QTL", "Background"), c("In_Annot", "Not_in_Annot")))
fisher_result <- fisher.test(mat)
fe <- (a / (a + b)) / (c / (c + d))
out <- data.frame(
    annotation = "$TYPE",
    p_value = fisher_result\$p.value,
    odds_ratio = fisher_result\$estimate,
    fold_enrichment = fe,
    ci1 = fisher_result\$conf.int[1],
    ci2 = fisher_result\$conf.int[2]
)
write.table(out, file="$TMP", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
EOF

    Rscript temp_fisher.R > /dev/null

done

# Add FDR correction and write the final result
cat <<EOF > temp_fdr.R
df <- read.table("$TMP", header=FALSE, sep="\t",
    col.names=c("annotation", "p_value", "odds_ratio", "Fold_Enrichment", "confidence_interval_1", "confidence_interval_2"))
df\$FDR <- p.adjust(df\$p_value, method = "BH")
write.table(df, file="$OUT", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
EOF

Rscript temp_fdr.R > /dev/null

# Clear temporary documents
rm temp.annot.bed temp_fisher.R temp_fdr.R "$TMP"
