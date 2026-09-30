#!/usr/bin/env bash
#
# Convert FastQTL nominal eQTL summary files to SMR BESD format.
#
# USAGE
#   bash fastqtl_to_besd.sh \
#     --input-dir <dir> \
#     --output-dir <dir> \
#     --smr <path> \
#     [--pattern <glob>] \
#     [--jobs <n>]
#
# OPTIONAL ARGUMENTS
#   --pattern      Glob pattern for input files (default: *_eGene_ciseQTL.fastqtl.txt.gz)
#   --jobs         Number of parallel jobs (default: 1)
#   -h, --help     Show this help
#
# EXAMPLE
#   bash fastqtl_to_besd.sh \
#     --input-dir /path/to/fastqtl \
#     --output-dir /path/to/besd \
#     --smr /path/to/smr \
#     --jobs 4
#

set -euo pipefail

#  Defaults 
INPUT_DIR=""
OUTPUT_DIR=""
SMR=""
PATTERN="*_eGene_ciseQTL.fastqtl.txt.gz"
JOBS=1

usage() {
  sed -n '3,45p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

#  Parse arguments 
while [[ $# -gt 0 ]]; do
  case $1 in
    --input-dir)   INPUT_DIR="$2";  shift 2 ;;
    --output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
    --smr)         SMR="$2";        shift 2 ;;
    --pattern)     PATTERN="$2";    shift 2 ;;
    --jobs)        JOBS="$2";       shift 2 ;;
    -h|--help)     usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$SMR" ]]; then
  echo "Error: --input-dir, --output-dir and --smr are required."
  usage
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: input directory not found: $INPUT_DIR"
  exit 1
fi

if [[ ! -x "$SMR" ]]; then
  echo "Error: SMR executable not found or not executable: $SMR"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

#  Collect files 
shopt -s nullglob
files=("$INPUT_DIR"/$PATTERN)
shopt -u nullglob

if [[ ${#files[@]} -eq 0 ]]; then
  echo "Error: no files matching '$PATTERN' in $INPUT_DIR"
  exit 1
fi

echo "Found ${#files[@]} file(s) to process."

#  Worker function 
process_one() {
  local input_file="$1"
  local filename
  filename=$(basename "$input_file")

  # Strip suffix to get "Period_Tissue"
  local base
  base=$(echo "$filename" | sed 's/_eGene_ciseQTL\.fastqtl\.txt\.gz$//')

  local period tissue
  period=$(echo "$base" | cut -d'_' -f1)
  tissue=$(echo "$base" | cut -d'_' -f2-)

  if [[ -z "$period" || -z "$tissue" ]]; then
    echo "[WARN] Cannot parse Period/Tissue from $filename, skipping."
    return
  fi

  local temp_file="${OUTPUT_DIR}/${base}_temp.txt"
  local out_prefix="${OUTPUT_DIR}/${base}_eGene_ciseQTL"

  echo "[PROCESS] $filename -> ${base}_eGene_ciseQTL"

  if ! gunzip -c "$input_file" > "$temp_file"; then
    echo "[ERROR] Failed to unzip $filename"
    return
  fi

  if "$SMR" --eqtl-summary "$temp_file" \
            --fastqtl-nominal-format \
            --make-besd \
            --out "$out_prefix"; then
    echo "[SUCCESS] $base"
    rm -f "$temp_file"
  else
    echo "[ERROR] SMR failed for $filename"
    echo "        Temporary file kept: $temp_file"
  fi
}

export -f process_one
export OUTPUT_DIR SMR

#  Run 
if [[ "$JOBS" -gt 1 ]] && command -v parallel >/dev/null 2>&1; then
  printf "%s\n" "${files[@]}" | parallel -j "$JOBS" process_one {}
else
  for f in "${files[@]}"; do
    process_one "$f"
  done
fi

echo "All done. BESD files are in: $OUTPUT_DIR"