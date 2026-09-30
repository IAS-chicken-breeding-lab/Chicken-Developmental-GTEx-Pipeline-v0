#!/usr/bin/env bash
#
# Batch-run update_besd.R for multiple tissues, optionally inside Docker.
#
# Usage:
#   bash run_update_besd.sh \
#     --tss <file> \
#     --gtf <file> \
#     --frq-dir <dir> \
#     --besd-dir <dir> \
#     --output-dir <dir> \
#     --tissues <list or file> \
#     [--r-script <path>] \
#     [--docker-image <image>] \
#     [--jobs <n>]
#
#

set -euo pipefail

#  Defaults 
TSS=""
GTF=""
FRQ_DIR=""
BESD_DIR=""
OUTPUT_DIR=""
TISSUES=""
R_SCRIPT="update_besd.R"
DOCKER_IMAGE=""
JOBS=1

usage() {
  sed -n '3,30p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

#  Parse arguments 
while [[ $# -gt 0 ]]; do
  case $1 in
    --tss)          TSS="$2"; shift 2 ;;
    --gtf)          GTF="$2"; shift 2 ;;
    --frq-dir)      FRQ_DIR="$2"; shift 2 ;;
    --besd-dir)     BESD_DIR="$2"; shift 2 ;;
    --output-dir)   OUTPUT_DIR="$2"; shift 2 ;;
    --tissues)      TISSUES="$2"; shift 2 ;;
    --r-script)     R_SCRIPT="$2"; shift 2 ;;
    --docker-image) DOCKER_IMAGE="$2"; shift 2 ;;
    --jobs)         JOBS="$2"; shift 2 ;;
    -h|--help)      usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$TSS" || -z "$GTF" || -z "$FRQ_DIR" || -z "$BESD_DIR" || -z "$OUTPUT_DIR" || -z "$TISSUES" ]]; then
  echo "Error: missing required arguments."
  usage
fi

mkdir -p "$OUTPUT_DIR"

#  Resolve tissue list 
if [[ -f "$TISSUES" ]]; then
  mapfile -t tissue_array < "$TISSUES"
else
  IFS=',' read -r -a tissue_array <<< "$TISSUES"
fi

if [[ ${#tissue_array[@]} -eq 0 ]]; then
  echo "Error: no tissues provided."
  exit 1
fi

#  Worker function 
process_tissue() {
  local tis="$1"

  if [[ -n "$DOCKER_IMAGE" ]]; then
    echo "[DOCKER] Processing $tis ..."
    docker run --rm \
      -v /ossfs:/ossfs \
      -v /data:/data \
      -v /juice_data:/juice_data \
      -w /data \
      "$DOCKER_IMAGE" \
      Rscript "$R_SCRIPT" "$TSS" "$GTF" "$FRQ_DIR" "$BESD_DIR" "$OUTPUT_DIR" "$tis"
  else
    echo "[LOCAL] Processing $tis ..."
    Rscript "$R_SCRIPT" "$TSS" "$GTF" "$FRQ_DIR" "$BESD_DIR" "$OUTPUT_DIR" "$tis"
  fi
}

export -f process_tissue
export TSS GTF FRQ_DIR BESD_DIR OUTPUT_DIR R_SCRIPT DOCKER_IMAGE

#  Run 
if [[ "$JOBS" -gt 1 ]] && command -v parallel >/dev/null 2>&1; then
  printf "%s\n" "${tissue_array[@]}" | parallel -j "$JOBS" process_tissue {}
else
  for tis in "${tissue_array[@]}"; do
    process_tissue "$tis"
  done
fi

echo "All tissues processed."