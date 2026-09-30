#!/usr/bin/env bash
#
# extract_top_fst.sh
# Extract top N% FST windows from a windowed FST file.

set -euo pipefail

#  Defaults
INPUT_FILE=""
OUT_DIR=""
FST_COL=5
CHR_COL=1
START_COL=2
PERCENTAGES="1 2 3 4 5 6 7 8 9 10"

usage() {
  cat <<EOF
Usage: $0 --input <fst_file> --outdir <dir> [options]

Required:
  --input      Input FST window file (with header)
  --outdir     Output directory (created if missing)

Optional:
  --fst-col    Column index of WEIGHTED_FST (default: 5)
  --chr-col    Column index of CHROM       (default: 1)
  --start-col  Column index of BIN_START   (default: 2)
  --percentages  Space-separated list of percentages to extract
                 (default: "1 2 3 4 5 6 7 8 9 10")
  -h, --help   Show this help

Example:
  bash $0 --input Layer_Fst.filter --outdir Layer-Group
EOF
  exit 0
}

#  Parse args 
while [[ $# -gt 0 ]]; do
  case $1 in
    --input)       INPUT_FILE="$2"; shift 2 ;;
    --outdir)      OUT_DIR="$2"; shift 2 ;;
    --fst-col)     FST_COL="$2"; shift 2 ;;
    --chr-col)     CHR_COL="$2"; shift 2 ;;
    --start-col)   START_COL="$2"; shift 2 ;;
    --percentages) PERCENTAGES="$2"; shift 2 ;;
    -h|--help)     usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$INPUT_FILE" || -z "$OUT_DIR" ]]; then
  echo "Error: --input and --outdir are required."
  usage
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: input file not found: $INPUT_FILE"
  exit 1
fi

mkdir -p "$OUT_DIR"

#  Count windows 
TOTAL=$(tail -n +2 "$INPUT_FILE" | wc -l)
echo "Total windows: $TOTAL"

# Header
HEADER=$(head -n 1 "$INPUT_FILE")

#  Function to extract top N% 
extract_top() {
  local pct="$1"
  local num_lines="$2"
  local outfile="${OUT_DIR}/Top${pct}percent_FST_windows.txt"

  echo "Writing ${outfile} ..."

  {
    echo "$HEADER"
    tail -n +2 "$INPUT_FILE" \
      | sort -k${FST_COL},${FST_COL}gr \
      | head -n "$num_lines" \
      | sort -k${CHR_COL},${CHR_COL}n -k${START_COL},${START_COL}n
  } > "$outfile"
}

#  Extract for each percentage 
for pct in $PERCENTAGES; do
  num=$(( TOTAL * pct / 100 ))
  [[ $num -lt 1 ]] && num=1
  echo "Top ${pct}% -> ${num} windows"
  extract_top "$pct" "$num"
done

#  Summary 
echo
echo "Done. Results in: ${OUT_DIR}"
ls -lh "${OUT_DIR}/Top"*"_FST_windows.txt"

echo -e "\nSummary:"
for file in "${OUT_DIR}"/Top*_FST_windows.txt; do
  if [[ -f "$file" ]]; then
    lines=$(wc -l < "$file")
    data_lines=$((lines - 1))
    if [[ $data_lines -gt 0 ]]; then
      min_fst=$(tail -n +2 "$file" | awk -v col="$FST_COL" '{print $col}' | sort -n | head -1)
      max_fst=$(tail -n +2 "$file" | awk -v col="$FST_COL" '{print $col}' | sort -nr | head -1)
      echo "$(basename "$file"): $data_lines windows, FST range: $min_fst - $max_fst"
    else
      echo "$(basename "$file"): empty"
    fi
  fi
done