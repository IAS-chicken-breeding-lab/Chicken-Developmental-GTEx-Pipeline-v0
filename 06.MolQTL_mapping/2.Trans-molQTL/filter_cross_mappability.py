#!/usr/bin/env python3
"""Remove likely trans-eQTL artifacts caused by cross-mapping.

A pair is removed when its phenotype gene has a cross-mappable partner and the
variant falls inside that partner's supplied cis/TSS window.
"""

import argparse
import csv
import gzip
import re
from collections import defaultdict


def open_text(path, mode="rt"):
    return gzip.open(path, mode, newline="") if path.endswith(".gz") else open(path, mode, newline="")


def norm_chrom(value):
    return re.sub(r"^chr", "", str(value), flags=re.IGNORECASE)


def parse_variant(value):
    match = re.match(r"^(?:chr)?([^_:]+)[_:](\d+)", value)
    if not match:
        raise ValueError(f"Cannot parse variant_id as chromosome/position: {value}")
    return norm_chrom(match.group(1)), int(match.group(2))


parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--crossmap", required=True)
parser.add_argument("--tss-windows", required=True)
parser.add_argument("--tissue", required=True)
parser.add_argument("--stats", required=True)
args = parser.parse_args()

partners = defaultdict(set)
with open_text(args.crossmap) as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2 or fields[0].lower() in {"gene", "gene1", "gene_id"}:
            continue
        a, b = fields[0], fields[1]
        partners[a].add(b)
        partners[b].add(a)

windows = {}
with open_text(args.tss_windows) as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4 or fields[0] == "gene_id":
            continue
        windows[fields[0]] = (norm_chrom(fields[1]), int(fields[2]), int(fields[3]))

total = filtered = kept = 0
with open_text(args.input) as src, open_text(args.output, "wt") as dst:
    reader = csv.DictReader(src, delimiter="\t")
    if not reader.fieldnames or not {"pheno_id", "variant_id"}.issubset(reader.fieldnames):
        raise SystemExit("Input must contain pheno_id and variant_id columns")
    fields = list(reader.fieldnames)
    if "tissue" not in fields:
        fields.insert(0, "tissue")
    writer = csv.DictWriter(dst, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in reader:
        total += 1
        chrom, pos = parse_variant(row["variant_id"])
        artifact = False
        for partner in partners.get(row["pheno_id"], ()):
            if partner not in windows:
                continue
            partner_chrom, start, end = windows[partner]
            # BED interval is zero-based, half-open; variant position is one-based.
            if chrom == partner_chrom and start < pos <= end:
                artifact = True
                break
        if artifact:
            filtered += 1
            continue
        row["tissue"] = row.get("tissue") or args.tissue
        writer.writerow(row)
        kept += 1

with open(args.stats, "w", encoding="utf-8") as out:
    out.write("Tissue\tInput_pairs\tFiltered_pairs\tRetained_pairs\tFiltered_fraction\n")
    out.write(f"{args.tissue}\t{total}\t{filtered}\t{kept}\t{filtered / total if total else 0:.8g}\n")
