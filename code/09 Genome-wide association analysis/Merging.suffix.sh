#!/bin/bash

# 获取所有后缀
find . -maxdepth 1 -type f -name "chr*" \
    | sed 's|^\./||' \
    | sed 's/^chr[0-9A-Z]*\.//' \
    | sort -u > suffix.list

# 并行处理每个 suffix
cat suffix.list | parallel --jobs 100 '
    suf={}

    echo "Processing suffix: $suf"
    out="indelSMALL.$suf"

    # 提取表头（找到的第一个文件）
    find . -maxdepth 1 -type f -name "chr*.$suf" \
        | head -n 1 \
        | xargs head -n 1 > "$out"

    # 追加内容
    find . -maxdepth 1 -type f -name "chr*.$suf" \
        | while read f; do
            tail -n +2 "$f"
        done >> "$out"

    echo "Done: $out"
'
