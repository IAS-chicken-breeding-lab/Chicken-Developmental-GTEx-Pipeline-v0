#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 1 ]] || { echo "Usage: $0 CONFIG" >&2; exit 1; }
CONFIG=$1
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
[[ -f "$CONFIG" ]] || { echo "Missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

for tissue in "${TISSUES[@]}"; do
  out="$ANALYSIS_DIR/$tissue"
  result_dir="$out/results/omiga/trans_lmm"
  inputs=()
  for ((task=1; task<=N_TASKS; task++)); do
    file="$result_dir/$tissue.trans_qtl_pairs.task_${task}.crossmap.txt.gz"
    [[ -s "$file" ]] || { echo "Missing filtered task: $file" >&2; exit 1; }
    inputs+=("$file")
  done

  merged="$result_dir/$tissue.trans_qtl_pairs.crossmap.txt.gz"
  python3 "$SCRIPT_DIR/merge_gzip_tables.py" --output "$merged" "${inputs[@]}"
  prefix="$result_dir/$tissue.trans_qtl_pairs.fdr${FDR_THRESHOLD}"
  Rscript "$SCRIPT_DIR/07_call_trans_egenes.R" \
    "$merged" "$prefix" "$FDR_THRESHOLD" "$EFFECTIVE_TRANS_TESTS" \
    >"$out/logs/05_trans_eGene.log" 2>&1
done

echo "Merged task results and called trans-eGenes"

