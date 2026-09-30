#!/usr/bin/env python3
"""Create gene-centered TSS windows from a GTF.

Output columns are: gene_id, chromosome, zero-based start, one-based end.
"""

import argparse
import re


def get_attr(attributes: str, key: str):
    match = re.search(rf'(?:^|;\s*){re.escape(key)}\s+["\']?([^"\';]+)', attributes)
    return match.group(1).strip() if match else None


parser = argparse.ArgumentParser()
parser.add_argument("--gtf", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--window", type=int, default=1_000_000)
args = parser.parse_args()

seen = set()
with open(args.gtf, encoding="utf-8") as src, open(args.output, "w", encoding="utf-8") as out:
    out.write("gene_id\tchrom\tstart\tend\n")
    for line in src:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        gene = get_attr(fields[8], "gene_id")
        if not gene or gene in seen:
            continue
        seen.add(gene)
        chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
        tss = start if strand == "+" else end
        window_start = max(0, tss - 1 - args.window)
        window_end = tss + args.window
        out.write(f"{gene}\t{chrom}\t{window_start}\t{window_end}\n")

