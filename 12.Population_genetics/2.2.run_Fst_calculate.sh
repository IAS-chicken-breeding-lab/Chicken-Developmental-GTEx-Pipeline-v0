#!/usr/bin/env bash
#
# Wrapper for fst_from_merged_vcf.sh with preset parameters.

set -euo pipefail

VCF="merged.vcf.gz"
POP1="pop1.samples.txt"
POP2="pop2.samples.txt"
OUTDIR="results"
WIN=50000
STEP=10000
MAF=0.01
MAX_MISS=0.8
MIN_SNPS=10
THREADS=1

bash 2.1.Fst_calculate.sh \
  --vcf "$VCF" \
  --pop1 "$POP1" \
  --pop2 "$POP2" \
  --outdir "$OUTDIR" \
  --win "$WIN" \
  --step "$STEP" \
  --maf "$MAF" \
  --max-miss "$MAX_MISS" \
  --min-snps "$MIN_SNPS" \
  --threads "$THREADS"