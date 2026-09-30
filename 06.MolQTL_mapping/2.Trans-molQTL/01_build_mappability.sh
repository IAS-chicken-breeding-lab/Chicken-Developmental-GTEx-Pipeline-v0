#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 CONFIG [N_PARALLEL_CHROMOSOMES]" >&2
  exit 1
}

[[ $# -ge 1 && $# -le 2 ]] || usage
CONFIG=$1
JOBS=${2:-1}
[[ -f "$CONFIG" ]] || { echo "Missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

for tool in samtools genmap awk bedtools; do
  command -v "$tool" >/dev/null || { echo "Missing executable: $tool" >&2; exit 1; }
done
[[ -f "$REFERENCE_FASTA" ]] || { echo "Missing FASTA: $REFERENCE_FASTA" >&2; exit 1; }
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || { echo "N_PARALLEL_CHROMOSOMES must be a positive integer" >&2; exit 1; }

mkdir -p "$MAPPABILITY_DIR/chromosomes" "$MAPPABILITY_DIR/logs"

run_chromosome() {
  local chrom=$1
  local work="$MAPPABILITY_DIR/chromosomes/$chrom"
  mkdir -p "$work"
  samtools faidx "$REFERENCE_FASTA" "$chrom" > "$work/$chrom.fa"
  ln -sfn "$work/$chrom.fa" "$MAPPABILITY_DIR/chromosomes/$chrom.fa"
  genmap index -F "$work/$chrom.fa" -I "$work/${chrom}_index"
  genmap map -K "$UTR_K" -E "$MISMATCHES" -I "$work/${chrom}_index" \
    -O "$work/${chrom}_${UTR_K}_${MISMATCHES}" -t -w -bg
  genmap map -K "$EXON_K" -E "$MISMATCHES" -I "$work/${chrom}_index" \
    -O "$work/${chrom}_${EXON_K}_${MISMATCHES}" -t -w -bg
}

active=0
for chrom in "${CHROMOSOMES[@]}"; do
  run_chromosome "$chrom" >"$MAPPABILITY_DIR/logs/${chrom}.log" 2>&1 &
  active=$((active + 1))
  if (( active >= JOBS )); then
    wait -n
    active=$((active - 1))
  fi
done
wait

for k in "$UTR_K" "$EXON_K"; do
  merged="$MAPPABILITY_DIR/chicken_${k}_${MISMATCHES}.bed"
  : > "$merged"
  for chrom in "${CHROMOSOMES[@]}"; do
    f="$MAPPABILITY_DIR/chromosomes/$chrom/${chrom}_${k}_${MISMATCHES}.bedgraph"
    [[ -s "$f" ]] || { echo "Missing chromosome output: $f" >&2; exit 1; }
    awk -v c="$chrom" 'BEGIN{OFS="\t"} {$1=c; print}' "$f" >> "$merged"
  done
  awk '$4 < 0.5' "$merged" | bedtools sort -i - | bedtools merge -i - \
    > "$MAPPABILITY_DIR/chicken_${k}_${MISMATCHES}_lt0.5.bed"
done

selected="$MAPPABILITY_DIR/chicken_${SNP_MAPPABILITY_K}_${MISMATCHES}.bed"
[[ -s "$selected" ]] || { echo "Missing selected SNP mappability track: $selected" >&2; exit 1; }
awk '$4 == 1' "$selected" | bedtools sort -i - | bedtools merge -i - \
  > "$MAPPABILITY_DIR/snp_mappability1.bed"

echo "Mappability tracks written to $MAPPABILITY_DIR"
