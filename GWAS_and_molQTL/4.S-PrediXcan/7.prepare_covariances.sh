#!/usr/bin/env bash
#
# prepare_covariances.sh
# Merge and deduplicate covariance files for S-PrediXcan.
#
# USAGE
#   bash prepare_covariances.sh \
#     --src-dir <dir> \
#     --out-dir <dir> \
#     [--tissues <list or file>] \
#     [--tissue-pattern <glob>] \
#     [--dedup] \
#     [--header <string>]
#
# REQUIRED ARGUMENTS
#   --src-dir     Directory containing per-tissue subdirs with
#                 *_chr*_covariances.txt files
#   --out-dir     Output directory for merged and deduplicated files
#
# OPTIONAL ARGUMENTS
#   --tissues         Comma-separated tissue list, or a file with one tissue
#                     per line. If omitted, all subdirs in --src-dir are used.
#   --tissue-pattern  Glob pattern for tissue directories (default: "*")
#   --dedup           Deduplicate merged files by GENE + RSID1 + RSID2
#   --header          Header line for merged file
#                     (default: "GENE RSID1 RSID2 VALUE")
#   -h, --help        Show this help
#


set -euo pipefail

#  Defaults 
SRC_DIR=""
OUT_DIR=""
TISSUES=""
TISSUE_PATTERN="*"
DEDUP=false
HEADER="GENE RSID1 RSID2 VALUE"

usage() {
  sed -n '3,50p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

#  Parse args 
while [[ $# -gt 0 ]]; do
  case $1 in
    --src-dir)        SRC_DIR="$2";        shift 2 ;;
    --out-dir)        OUT_DIR="$2";        shift 2 ;;
    --tissues)        TISSUES="$2";        shift 2 ;;
    --tissue-pattern) TISSUE_PATTERN="$2"; shift 2 ;;
    --dedup)          DEDUP=true;          shift ;;
    --header)         HEADER="$2";         shift 2 ;;
    -h|--help)        usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$SRC_DIR" || -z "$OUT_DIR" ]]; then
  echo "Error: --src-dir and --out-dir are required."
  usage
fi

if [[ ! -d "$SRC_DIR" ]]; then
  echo "Error: source directory not found: $SRC_DIR"
  exit 1
fi

mkdir -p "$OUT_DIR"

#  Tissue list 
if [[ -n "$TISSUES" ]]; then
  if [[ -f "$TISSUES" ]]; then
    mapfile -t tissue_array < "$TISSUES"
  else
    IFS=',' read -r -a tissue_array <<< "$TISSUES"
  fi
else
  shopt -s nullglob
  tissue_array=()
  for d in "$SRC_DIR"/$TISSUE_PATTERN; do
    [[ -d "$d" ]] && tissue_array+=("$(basename "$d")")
  done
  shopt -u nullglob
fi

if [[ ${#tissue_array[@]} -eq 0 ]]; then
  echo "Error: no tissue directories found in $SRC_DIR"
  exit 1
fi

echo "Prepare covariance files for S-PrediXcan"
echo "  Source dir : $SRC_DIR"
echo "  Output dir : $OUT_DIR"
echo "  Tissues    : ${tissue_array[*]}"
echo "  Dedup      : $DEDUP"

#  Process each tissue 
merge_and_dedup() {
  local tissue="$1"
  local src_tissue_dir="$SRC_DIR/$tissue"
  local out_tissue_dir="$OUT_DIR/$tissue"
  local out_file="$out_tissue_dir/${tissue}_Model_training_covariances.txt.gz"
  local bak_file="$out_tissue_dir/${tissue}_Model_training_covariancesold.txt.gz"

  echo ">>> Processing: $tissue"

  if [[ ! -d "$src_tissue_dir" ]]; then
    echo "  [WARN] source dir not found: $src_tissue_dir"
    return
  fi

  shopt -s nullglob
  local cov_files=("$src_tissue_dir"/*_chr*_covariances.txt)
  shopt -u nullglob

  if [[ ${#cov_files[@]} -eq 0 ]]; then
    echo "  [INFO] no covariance files in $src_tissue_dir"
    return
  fi

  mkdir -p "$out_tissue_dir"
  local tmpfile
  tmpfile=$(mktemp)

  # Merge all chromosome files, drop original headers
  for f in "${cov_files[@]}"; do
    awk 'NR>1' "$f" >> "$tmpfile"
  done

  # Add unified header
  sed -i "1i${HEADER}" "$tmpfile"

  # Deduplicate by GENE + RSID1 + RSID2 if requested
  if [[ "$DEDUP" == true ]]; then
    local tmp_dedup
    tmp_dedup=$(mktemp)
    awk '!seen[$1 FS $2 FS $3]++' "$tmpfile" > "$tmp_dedup"

    if [[ -f "$out_file" ]]; then
      mv "$out_file" "$bak_file"
      echo "  Backup: $(basename "$bak_file")"
    fi

    gzip -c "$tmp_dedup" > "$out_file"
    rm -f "$tmp_dedup"
  else
    gzip -c "$tmpfile" > "$out_file"
  fi

  rm -f "$tmpfile"
  echo "  Written: $(basename "$out_file")"
}

for tissue in "${tissue_array[@]}"; do
  merge_and_dedup "$tissue"
done

echo "All done at $(date)"