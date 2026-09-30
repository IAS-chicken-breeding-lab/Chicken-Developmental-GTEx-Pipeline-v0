#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 1 ]] || { echo "Usage: $0 CONFIG" >&2; exit 1; }
CONFIG=$1
[[ -f "$CONFIG" ]] || { echo "Missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
[[ -x "$OMIGA_BIN" ]] || { echo "OmiGA is not executable: $OMIGA_BIN" >&2; exit 1; }

run_task() {
  local tissue=$1 task=$2 out result_dir
  out="$ANALYSIS_DIR/$tissue"
  result_dir="$out/results/omiga/trans_lmm"
  mkdir -p "$result_dir" "$out/logs"
  "$OMIGA_BIN" --mode trans --verbose --threads "$THREADS_PER_TASK" \
    --multi-task "$N_TASKS" "$task" \
    --covariates "$out/covar/${tissue}_combined_covar.txt" \
    --dcovar-name Sex Period \
    --genotype "$out/geno/$tissue" --phenotype "$out/pheno/$tissue.bed" \
    --prefix "$tissue" --rm-collinear-covar 0.95 --output-dir "$result_dir" \
    >"$out/logs/03_trans_${task}.log" 2>&1
}

active=0
for tissue in "${TISSUES[@]}"; do
  for ((task=1; task<=N_TASKS; task++)); do
    run_task "$tissue" "$task" &
    active=$((active + 1))
    if (( active >= MAX_PARALLEL_TASKS )); then
      wait -n
      active=$((active - 1))
    fi
  done
done
wait
echo "All OmiGA trans-eQTL tasks completed"

