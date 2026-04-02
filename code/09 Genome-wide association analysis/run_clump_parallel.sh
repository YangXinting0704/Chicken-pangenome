#!/usr/bin/env bash
set -euo pipefail

############################################
# 需要你确认 / 修改的参数
############################################

# PLINK 二进制基因型前缀（不要加 .bed）
BFILE="final.imputed.SVINDELSNP"

# clump 参数（与你论文里保持一致即可）
P1="4.11957e-7"
P2="8.23913e-6"
R2="0.2"
KB="1000"

# 并行任务数（根据机器核数改）
N_JOBS=150

############################################
# 不需要再改的部分
############################################

run_one_clump () {
    infile="$1"
    base=$(basename "$infile" .assoc.txt)

    tmp="temporary.${base}.assoc.txt"
    out="clumped.${base}"

    echo "[INFO] Processing $infile"

    # 1. 生成 plink --clump 标准输入文件
    #    必须包含 header：SNP CHR BP P
    awk '
    BEGIN {print "SNP\tCHR\tBP\tP"}
    NR>1 {print $2, $1, $3, $12}
    ' "$infile" > "$tmp"

    # 2. 运行 plink clump
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

    # 3. 删除临时文件
    rm -f "$tmp"

    echo "[INFO] Finished $infile"
}

export -f run_one_clump
export BFILE P1 P2 R2 KB

############################################
# 并行运行所有 FINAL.*.assoc.txt
############################################

ls FINAL.*.assoc.txt | \
    parallel -j "$N_JOBS" run_one_clump {}
