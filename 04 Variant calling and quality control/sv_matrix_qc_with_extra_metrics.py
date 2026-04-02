#!/usr/bin/env python3
"""
sv_matrix_qc_with_extra_metrics.py

Usage:
python sv_matrix_qc_with_extra_metrics.py \
  --tgs IDextract.TGS01.vcf \
  --ngs IDextract.NGS02.vcf \
  --out per_variant_new_TGS01_NGS02.tsv
"""

import argparse
import sys
import gzip
import csv
from collections import OrderedDict


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    else:
        return open(path, "r")


def normalize_gt(raw):
    """
    Normalize genotype string to one of:
      '0/0', '0/1', '1/1', './.' or 'other'

    - Accepts '|' or '/' separators
    - Accepts leading fields with ':' (like '0/1:...'), only keep first part
    - Treat '.', '.|.', './.' as missing -> './.'
    """
    if raw is None:
        return "./."
    raw = raw.strip()
    if raw == "":
        return "./."

    # if colon present, take first field
    if ":" in raw:
        raw = raw.split(":", 1)[0]

    # missing forms
    if raw in (".", "./.", ".|."):
        return "./."

    # unify pipes to slashes
    raw = raw.replace("|", "/")

    # expected forms
    if raw in ("0/0", "0/1", "1/0", "1/1"):
        if raw == "1/0":
            return "0/1"
        return raw

    return "other"


def read_matrix(path):
    """
    Read a tab-delimited matrix file where first column is ID and rest are samples.
    Returns:
      sample_names (list),
      variants OrderedDict {variant_id: [genotypes aligned with sample_names]}
    Ignores comment lines starting with '#'
    """
    with open_maybe_gz(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = None
        variants = OrderedDict()

        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue

            if header is None:
                header = row
                if len(header) < 2:
                    raise ValueError(
                        f"Header in {path} must have at least 2 columns (ID + samples). Found: {header}"
                    )
                sample_names = header[1:]
                continue

            variant_id = row[0]
            genos = row[1:]

            # pad / trim to match header sample count
            if len(genos) < len(sample_names):
                genos = genos + [""] * (len(sample_names) - len(genos))
            elif len(genos) > len(sample_names):
                genos = genos[: len(sample_names)]

            variants[variant_id] = genos

    return sample_names, variants


def main():
    p = argparse.ArgumentParser(
        description="Per-variant QC between two genotype-matrix files (TGS vs NGS)."
    )
    p.add_argument(
        "--tgs",
        required=True,
        help="TGS matrix file (tab-delimited). First column variant ID, rest sample genotypes.",
    )
    p.add_argument(
        "--ngs",
        required=True,
        help="NGS matrix file (same format as TGS).",
    )
    p.add_argument("--out", required=True, help="Output TSV file path.")
    args = p.parse_args()

    # read files
    try:
        tgs_samples, tgs_vars = read_matrix(args.tgs)
    except Exception as e:
        print(f"ERROR reading TGS file {args.tgs}: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        ngs_samples, ngs_vars = read_matrix(args.ngs)
    except Exception as e:
        print(f"ERROR reading NGS file {args.ngs}: {e}", file=sys.stderr)
        sys.exit(1)

    # align sample sets by intersection (keep TGS order)
    set_tgs = list(tgs_samples)
    set_ngs = list(ngs_samples)
    inter = [s for s in set_tgs if s in set_ngs]

    if len(inter) == 0:
        print("ERROR: No shared samples between TGS and NGS headers.", file=sys.stderr)
        sys.exit(1)

    if len(inter) != len(set_tgs) or len(inter) != len(set_ngs):
        print(
            f"Warning: sample sets differ. "
            f"TGS samples: {len(set_tgs)}, NGS samples: {len(set_ngs)}, intersection used: {len(inter)}. "
            f"Samples not in intersection will be ignored.",
            file=sys.stderr,
        )

    # build index maps
    tgs_idx = {s: i for i, s in enumerate(tgs_samples)}
    ngs_idx = {s: i for i, s in enumerate(ngs_samples)}
    use_idx_tgs = [tgs_idx[s] for s in inter]
    use_idx_ngs = [ngs_idx[s] for s in inter]

    # union of variant IDs (keep order)
    all_variants = list(
        OrderedDict.fromkeys(list(tgs_vars.keys()) + list(ngs_vars.keys()))
    )

    with open(args.out, "w") as outfh:
        header = [
            "variant_id",
            "n_samples",
            "tgs_alt_called",
            "ngs_alt_called",
            "both_alt_called",
            "concordant_samples",
            "concordance_rate",
            "non_missing_samples",
            "concordant_non_missing",
            "concordance_rate_non_missing",
            # NEW: 4x4 checked matrix proportion (denominator = all samples)
            "checked_count_custom",
            "checked_rate_custom",
        ]
        outfh.write("\t".join(header) + "\n")

        for vid in all_variants:
            # get genotype arrays aligned to intersection samples
            if vid in tgs_vars:
                tgs_row = tgs_vars[vid]
                tgs_gts = [normalize_gt(tgs_row[i]) for i in use_idx_tgs]
            else:
                tgs_gts = ["./."] * len(inter)

            if vid in ngs_vars:
                ngs_row = ngs_vars[vid]
                ngs_gts = [normalize_gt(ngs_row[i]) for i in use_idx_ngs]
            else:
                ngs_gts = ["./."] * len(inter)

            n_samples = len(inter)

            tgs_alt = 0
            ngs_alt = 0
            both_alt = 0

            concordant = 0  # exact match among (0/0,0/1,1/1,./.)
            non_missing_samples = 0  # both not missing
            concordant_non_missing = 0  # both not missing and exact match among (0/0,0/1,1/1)

            # NEW: your 4x4 checked matrix count (denominator = all samples)
            checked_count_custom = 0

            for gt_t, gt_n in zip(tgs_gts, ngs_gts):
                # alt-called counts
                if gt_t in ("0/1", "1/1"):
                    tgs_alt += 1
                if gt_n in ("0/1", "1/1"):
                    ngs_alt += 1
                if (gt_t in ("0/1", "1/1")) and (gt_n in ("0/1", "1/1")):
                    both_alt += 1

                # overall concordance: exact match among (0/0,0/1,1/1,./.)
                if gt_t == gt_n and gt_t in ("0/0", "0/1", "1/1", "./."):
                    concordant += 1

                # non-missing concordance: only when both != ./.
                if gt_t != "./." and gt_n != "./.":
                    non_missing_samples += 1
                    if gt_t == gt_n and gt_t in ("0/0", "0/1", "1/1"):
                        concordant_non_missing += 1

                # NEW: 4x4 checked matrix (as you provided)
                # columns: 0/0 0/1 1/1 ./.
                # rows:    0/0 0/1 1/1 ./.
                #
                # allowed (√):
                # 0/0 vs 0/0, 0/0 vs ./.
                # 0/1 vs 0/1, 0/1 vs 1/1
                # 1/1 vs 0/1, 1/1 vs 1/1
                # ./. vs 0/0, ./. vs ./.
                if (
                    (gt_t == "0/0" and gt_n in ("0/0", "./.")) or
                    (gt_t == "0/1" and gt_n in ("0/1", "1/1")) or
                    (gt_t == "1/1" and gt_n in ("0/1", "1/1")) or
                    (gt_t == "./." and gt_n in ("0/0", "./."))
                ):
                    checked_count_custom += 1

            concordance_rate = concordant / n_samples if n_samples > 0 else 0.0
            concordance_rate_non_missing = (
                concordant_non_missing / non_missing_samples
                if non_missing_samples > 0
                else 0.0
            )

            checked_rate_custom = (
                checked_count_custom / n_samples if n_samples > 0 else 0.0
            )

            outfh.write(
                "\t".join(
                    [
                        vid,
                        str(n_samples),
                        str(tgs_alt),
                        str(ngs_alt),
                        str(both_alt),
                        str(concordant),
                        f"{concordance_rate:.4f}",
                        str(non_missing_samples),
                        str(concordant_non_missing),
                        f"{concordance_rate_non_missing:.4f}",
                        str(checked_count_custom),
                        f"{checked_rate_custom:.4f}",
                    ]
                )
                + "\n"
            )

    print(f"Done. Results written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
