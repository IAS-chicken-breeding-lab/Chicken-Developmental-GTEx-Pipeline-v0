#!/usr/bin/env bash
set -euo pipefail

INPUT="RJF_Fst.filter"
OUTDIR="RJF-Group"
PREFIX="RJF"

python3 3.1.Fst_grouping.py \
  --input "$INPUT" \
  --outdir "$OUTDIR" \
  --prefix "$PREFIX"