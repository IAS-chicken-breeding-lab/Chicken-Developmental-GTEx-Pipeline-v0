#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

bash 4.1.extract_top_fst.sh \
  --input Layer_Fst.filter \
  --outdir Layer-Group \
  --fst-col 5 \
  --chr-col 1 \
  --start-col 2 \
  --percentages "1 2 3 4 5 6 7 8 9 10"