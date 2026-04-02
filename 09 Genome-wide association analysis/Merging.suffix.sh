#!/bin/bash

# Get all suffixes
find . -maxdepth 1 -type f -name "chr*" \
    | sed 's|^\./||' \
    | sed 's/^chr[0-9A-Z]*\.//' \
    | sort -u > suffix.list

# Get all suffixes
cat suffix.list | parallel --jobs 100 '
    suf={}

    echo "Processing suffix: $suf"
    out="indelSMALL.$suf"

    # Extract the table header (the first file found)
    find . -maxdepth 1 -type f -name "chr*.$suf" \
        | head -n 1 \
        | xargs head -n 1 > "$out"

    find . -maxdepth 1 -type f -name "chr*.$suf" \
        | while read f; do
            tail -n +2 "$f"
        done >> "$out"

    echo "Done: $out"
'
