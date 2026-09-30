#!/usr/bin/env bash
#
# Convert cis-eQTL result files to FastQTL input format.
#
# USAGE
#   bash convert_to_fastqtl.sh \
#     --input-dir <dir> \
#     --output-dir <dir> \
#     --tissues <list> \
#     [--periods <list>] \
#     [--threads <n>] \
#     [--col-gene <n>] \
#     [--col-snp <n>] \
#     [--col-dist <n>] \
#     [--col-pval <n>] \
#     [--col-beta <n>]
#
# OPTIONAL ARGUMENTS
#   --periods      Comma-separated list of periods (default: A,B,C)
#   --threads      Number of parallel jobs (default: 4)
#   --col-gene     Column index for gene ID    (default: 1)
#   --col-snp      Column index for SNP ID     (default: 2)
#   --col-dist     Column index for distance   (default: 3)
#   --col-pval     Column index for p-value    (default: 7)
#   --col-beta     Column index for beta       (default: 5)
#   -h, --help     Show this help message
#
# EXAMPLE
#   bash convert_to_fastqtl.sh \
#     --input-dir /path/to/input \
#     --output-dir /path/to/output \
#     --tissues AF,PZ,CT,SWM,FZ,SZ,GZ,XJ,HC,XQN,JG,XW,KC,XZ,MC,ZC \
#     --periods A,B,C \
#     --threads 8
# 

set -euo pipefail

#  Defaults 
INPUT_DIR=""
OUTPUT_DIR=""
TISSUES=""
PERIODS="A,B,C"
THREADS=4
COL_GENE=1
COL_SNP=2
COL_DIST=3
COL_PVAL=7
COL_BETA=5
COMPRESS_CMD=""

usage() {
  sed -n '2,60p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

#  Parse arguments 
while [[ $# -gt 0 ]]; do
  case $1 in
    --input-dir)   INPUT_DIR="$2";   shift 2 ;;
    --output-dir)  OUTPUT_DIR="$2";  shift 2 ;;
    --tissues)     TISSUES="$2";     shift 2 ;;
    --periods)     PERIODS="$2";     shift 2 ;;
    --threads)     THREADS="$2";     shift 2 ;;
    --col-gene)    COL_GENE="$2";    shift 2 ;;
    --col-snp)     COL_SNP="$2";     shift 2 ;;
    --col-dist)    COL_DIST="$2";    shift 2 ;;
    --col-pval)    COL_PVAL="$2";    shift 2 ;;
    --col-beta)    COL_BETA="$2";    shift 2 ;;
    -h|--help)     usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$TISSUES" ]]; then
  echo "Error: --input-dir, --output-dir and --tissues are required."
  usage
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: input directory not found: $INPUT_DIR"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Convert comma-separated lists to arrays
IFS=',' read -r -a tissue_array <<< "$TISSUES"
IFS=',' read -r -a period_array <<< "$PERIODS"

#  Compression tool 
if command -v pigz >/dev/null 2>&1; then
  COMPRESS_CMD="pigz -c"
else
  echo "[WARN] pigz not found, falling back to gzip"
  COMPRESS_CMD="gzip -c"
fi

#  Worker function 
process_tissue() {
  local tissue="$1"
  local periods=("${!2}")

  echo "[INFO] Processing tissue: $tissue"

  for period in "${periods[@]}"; do
    local input_file="${INPUT_DIR}/${period}_${tissue}_eGene_ciseQTL.txt.gz"
    local output_file="${OUTPUT_DIR}/${period}_${tissue}_eGene_ciseQTL.fastqtl.txt.gz"

    if [[ ! -f "$input_file" ]]; then
      echo "[WARN] File not found: $input_file"
      continue
    fi

    echo "[PROCESS] $period$tissue -> $(basename "$output_file")"

    zcat "$input_file" | \
      awk -F'\t' -v OFS='\t' \
          -v c_gene="$COL_GENE" \
          -v c_snp="$COL_SNP" \
          -v c_dist="$COL_DIST" \
          -v c_pval="$COL_PVAL" \
          -v c_beta="$COL_BETA" '
        FNR==1 && $1=="pheno_id" && $2=="variant_id" { next }
        { print $c_gene, $c_snp, $c_dist, $c_pval, $c_beta }
      ' | \
      $COMPRESS_CMD > "$output_file"

    if [[ -s "$output_file" ]]; then
      local col_count
      col_count=$(zcat "$output_file" | head -1 | awk -F'\t' '{print NF}')
      local line_count
      line_count=$(zcat "$output_file" | wc -l)
      echo "[SUCCESS] $(basename "$output_file") (${line_count} lines, ${col_count} columns)"
    else
      echo "[ERROR] Failed to create: $output_file"
    fi
  done

  echo "[DONE] Tissue: $tissue"
}

export -f process_tissue
export INPUT_DIR OUTPUT_DIR COMPRESS_CMD
export COL_GENE COL_SNP COL_DIST COL_PVAL COL_BETA

#  Main loop 
echo "[START] Converting eQTL files to FastQTL format"
echo "[INFO] Input directory : $INPUT_DIR"
echo "[INFO] Output directory: $OUTPUT_DIR"
echo "[INFO] Tissues         : ${tissue_array[*]}"
echo "[INFO] Periods         : ${period_array[*]}"
echo "[INFO] Threads         : $THREADS"

if command -v parallel >/dev/null 2>&1; then
  printf "%s\n" "${tissue_array[@]}" | \
    parallel -j "$THREADS" process_tissue {} period_array
else
  echo "[WARN] GNU parallel not found, using background jobs"
  pids=()
  for tissue in "${tissue_array[@]}"; do
    process_tissue "$tissue" period_array &
    pids+=($!)
    while [[ ${#pids[@]} -ge $THREADS ]]; do
      wait -n
      new_pids=()
      for pid in "${pids[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
          new_pids+=("$pid")
        fi
      done
      pids=("${new_pids[@]}")
    done
  done
  wait
fi

#  Final summary 
echo "[SUMMARY] Conversion finished."
count=$(find "$OUTPUT_DIR" -name "*.fastqtl.txt.gz" | wc -l)
echo "[INFO] Total FastQTL files: $count"