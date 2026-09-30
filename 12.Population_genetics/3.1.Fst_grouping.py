#!/usr/bin/env python3
"""
fst_grouping.py

Group FST windows into deciles (10 equal-sized bins) based on WEIGHTED_FST,
and output BED and detailed tables for downstream analysis.

The deciles are ranked by FST:
  - Group10 contains the top 10% of windows with the highest WEIGHTED_FST.
  - Group9  contains the next 10%.
  - ...
  - Group1  contains the bottom 10% with the lowest WEIGHTED_FST.

Each group (except possibly the last) contains approximately the same
number of windows.

"""

import argparse
import os
import sys

import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Group FST windows into deciles based on WEIGHTED_FST."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input FST window file (e.g., .windowed.weir.fst or .filter)."
    )
    parser.add_argument(
        "-o", "--outdir", required=True,
        help="Output directory (will be created if missing)."
    )
    parser.add_argument(
        "--prefix", default=None,
        help="Prefix for output files. Default: basename of input file."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input_file = args.input
    outdir = args.outdir
    prefix = args.prefix or os.path.splitext(os.path.basename(input_file))[0]

    os.makedirs(outdir, exist_ok=True)

    print(f"[INFO] Reading: {input_file}")
    try:
        df = pd.read_csv(input_file, sep=r'\s+')
    except Exception as e:
        print(f"[ERROR] Failed to read input file: {e}", file=sys.stderr)
        sys.exit(1)

    required_cols = ["CHROM", "BIN_START", "BIN_END", "WEIGHTED_FST"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"[ERROR] Missing required columns: {missing}", file=sys.stderr)
        sys.exit(1)

    # Sort by WEIGHTED_FST descending
    df_sorted = df.sort_values("WEIGHTED_FST", ascending=False).reset_index(drop=True)

    n_windows = len(df_sorted)
    if n_windows == 0:
        print("[ERROR] No windows found in input file.", file=sys.stderr)
        sys.exit(1)

    # Divide into deciles: top 10% = Group10, bottom 10% = Group1
    windows_per_group = n_windows // 10
    df_sorted["Group"] = ""

    for i in range(10):
        start = i * windows_per_group
        if i < 9:
            end = (i + 1) * windows_per_group
        else:
            end = n_windows
        # Group10 = highest FST (first chunk), Group1 = lowest FST (last chunk)
        group_name = f"Group{10 - i}"
        df_sorted.loc[start:end - 1, "Group"] = group_name

    # Ensure CHROM has 'chr' prefix for BED
    df_sorted["CHROM"] = df_sorted["CHROM"].astype(str)
    df_sorted["CHROM"] = df_sorted["CHROM"].apply(
        lambda x: x if x.startswith("chr") else f"chr{x}"
    )

    # BED output
    bed_df = pd.DataFrame({
        "CHROM": df_sorted["CHROM"],
        "START": df_sorted["BIN_START"].astype(int),
        "END": df_sorted["BIN_END"].astype(int),
        "GROUP": df_sorted["Group"],
    })
    bed_df = bed_df.sort_values(["CHROM", "START"]).reset_index(drop=True)

    bed_file = os.path.join(outdir, f"{prefix}_fst_groups.bed")
    bed_df.to_csv(bed_file, sep="\t", header=False, index=False)
    print(f"[INFO] BED file written: {bed_file}")

    # Detailed output
    detailed_file = os.path.join(outdir, f"{prefix}_fst_groups_detailed.txt")
    df_sorted.to_csv(detailed_file, sep="\t", index=False)
    print(f"[INFO] Detailed file written: {detailed_file}")

    # Summary
    print("\n=== Decile summary (Group10 = highest FST, Group1 = lowest FST) ===")
    for gnum in range(1, 11):
        gname = f"Group{gnum}"
        sub = df_sorted[df_sorted["Group"] == gname]
        if len(sub) == 0:
            continue
        print(f"{gname}:")
        print(f"  Windows : {len(sub)}")
        print(f"  FST min : {sub['WEIGHTED_FST'].min():.6f}")
        print(f"  FST max : {sub['WEIGHTED_FST'].max():.6f}")
        print(f"  FST mean: {sub['WEIGHTED_FST'].mean():.6f}")
        chroms = sub["CHROM"].unique()
        print(f"  Chromosomes: {', '.join(map(str, chroms[:5]))}"
              + (" ..." if len(chroms) > 5 else ""))
        print()


if __name__ == "__main__":
    main()