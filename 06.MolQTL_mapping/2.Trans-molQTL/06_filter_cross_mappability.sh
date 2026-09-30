#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 1 ]] || { echo "Usage: $0 CONFIG" >&2; exit 1; }
CONFIG=$1
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
[[ -f "$CONFIG" ]] || { echo "Missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

crossmap="$CROSSMAP_DIR/cross_mappability.tsv"
[[ -s "$crossmap" ]] || { echo "Missing cross-mappability file: $crossmap" >&2; exit 1; }
[[ -s "$GENE_TSS_WINDOW_BED" ]] || { echo "Missing gene TSS windows: $GENE_TSS_WINDOW_BED" >&2; exit 1; }

for tissue in "${TISSUES[@]}"; do
  out="$ANALYSIS_DIR/$tissue"
  result_dir="$out/results/omiga/trans_lmm"
  mkdir -p "$result_dir/filter_stats"
  for ((task=1; task<=N_TASKS; task++)); do
    input="$result_dir/$tissue.trans_qtl_pairs.task_${task}.txt.gz"
    output="$result_dir/$tissue.trans_qtl_pairs.task_${task}.crossmap.txt.gz"
    [[ -s "$input" ]] || { echo "Missing task result: $input" >&2; exit 1; }
    python3 "$SCRIPT_DIR/filter_cross_mappability.py" \
      --input "$input" --output "$output" --crossmap "$crossmap" \
      --tss-windows "$GENE_TSS_WINDOW_BED" --tissue "$tissue" \
      --stats "$result_dir/filter_stats/task_${task}.tsv" \
      >"$out/logs/04_crossmap_filter_${task}.log" 2>&1
  done
done

echo "Cross-mappability filtering completed"

