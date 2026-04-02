#!/usr/bin/env bash
set -euo pipefail

############################################

BFILE="final.imputed.SVINDELSNP"

# clump parameter
P1="4.11957e-7"
P2="8.23913e-6"
R2="0.2"
KB="1000"

N_JOBS=150

############################################

run_one_clump () {
    infile="$1"
    base=$(basename "$infile" .assoc.txt)

    tmp="temporary.${base}.assoc.txt"
    out="clumped.${base}"

    echo "[INFO] Processing $infile"

    # 1. Generate the Plink-Clump standard input file
    awk '
    BEGIN {print "SNP\tCHR\tBP\tP"}
    NR>1 {print $2, $1, $3, $12}
    ' "$infile" > "$tmp"

    # 2. Run plink clump
    plink \
      --bfile "$BFILE" \
      --chr-set 50 \
      --allow-extra-chr \
      --clump "$tmp" \
      --clump-p1 "$P1" \
      --clump-p2 "$P2" \
      --clump-r2 "$R2" \
      --clump-kb "$KB" \
      --out "$out"

    # 3. Delete temporary files
    rm -f "$tmp"

    echo "[INFO] Finished $infile"
}

export -f run_one_clump
export BFILE P1 P2 R2 KB

############################################
# Run all FINAL.*.assoc.txt
############################################

ls FINAL.*.assoc.txt | \
    parallel -j "$N_JOBS" run_one_clump {}
