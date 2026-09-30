#!/usr/bin/env bash
#
# run_spredixcan.sh
# Batch-run S-PrediXcan (TWAS) for multiple eQTL models and GWAS traits
#
# USAGE
#   bash run_spredixcan.sh \
#     --spredixcan <path> \
#     --model-dir <dir> \
#     --cov-dir <dir> \
#     --gwas-dir <dir> \
#     --out-dir <dir> \
#     [--gwas-pattern <glob>] \
#     [--model-prefix <string>] \
#     [--sep <string>] \
#     [--jobs <n>] \
#     [--model-suffix <string>] \
#     [--cov-suffix <string>]
#
# REQUIRED ARGUMENTS
#   --spredixcan     Path to SPrediXcan.py
#   --model-dir      Directory containing eQTL model .db files
#   --cov-dir        Directory containing covariance files
#   --gwas-dir       Directory with GWAS summary files
#   --out-dir        Output directory (created if missing)
#
# OPTIONAL ARGUMENTS
#   --gwas-pattern   Glob pattern for GWAS files (default: "A_*_B_*.summary.txt")
#   --model-prefix   Prefix of model db files (default: "Chicken_")
#   --model-suffix   Suffix of model db files (default: "_ElasticNet_models_filtered_signif.db")
#   --cov-suffix     Suffix of covariance files (default: "_Model_training_covariances.txt.gz")
#   --sep            Separator between model tag and trait tag in the GWAS
#                    filename (default: "_B_")
#   --jobs           Max parallel jobs (default: 8)
#   -h, --help       Show this help
#
# EXAMPLE
#   bash run_spredixcan.sh \
#     --spredixcan /path/to/SPrediXcan.py \
#     --model-dir /path/to/ModelDB \
#     --cov-dir /path/to/Covariances \
#     --gwas-dir /path/to/GWAS \
#     --out-dir /path/to/Result \
#     --jobs 8
#

set -euo pipefail

#  Defaults 
SPREDIXCAN=""
MODEL_DIR=""
COV_DIR=""
GWAS_DIR=""
OUT_DIR=""
GWAS_PATTERN="A_*_B_*.summary.txt"
MODEL_PREFIX="Chicken_"
MODEL_SUFFIX="_ElasticNet_models_filtered_signif.db"
COV_SUFFIX="_Model_training_covariances.txt.gz"
SEP="_B_"
JOBS=8

# GWAS column names (adjust if your GWAS format differs)
SNP_COL="variant_id"
EA_COL="effect_allele"
NEA_COL="non_effect_allele"
BETA_COL="effect_size"
PVAL_COL="pvalue"

usage() {
  sed -n '3,70p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

#  Parse arguments 
while [[ $# -gt 0 ]]; do
  case $1 in
    --spredixcan)    SPREDIXCAN="$2";    shift 2 ;;
    --model-dir)     MODEL_DIR="$2";     shift 2 ;;
    --cov-dir)       COV_DIR="$2";       shift 2 ;;
    --gwas-dir)      GWAS_DIR="$2";      shift 2 ;;
    --out-dir)       OUT_DIR="$2";       shift 2 ;;
    --gwas-pattern)  GWAS_PATTERN="$2";  shift 2 ;;
    --model-prefix)  MODEL_PREFIX="$2";  shift 2 ;;
    --model-suffix)  MODEL_SUFFIX="$2";  shift 2 ;;
    --cov-suffix)    COV_SUFFIX="$2";    shift 2 ;;
    --sep)           SEP="$2";           shift 2 ;;
    --jobs)          JOBS="$2";          shift 2 ;;
    -h|--help)       usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

#  Validate 
if [[ -z "$SPREDIXCAN" || -z "$MODEL_DIR" || -z "$COV_DIR" || -z "$GWAS_DIR" || -z "$OUT_DIR" ]]; then
  echo "Error: --spredixcan, --model-dir, --cov-dir, --gwas-dir and --out-dir are required."
  usage
fi

for f in "$SPREDIXCAN"; do
  [[ -f "$f" ]] || { echo "Error: file not found: $f"; exit 1; }
done

for d in "$MODEL_DIR" "$COV_DIR" "$GWAS_DIR"; do
  [[ -d "$d" ]] || { echo "Error: directory not found: $d"; exit 1; }
done

mkdir -p "$OUT_DIR"

echo "Batch S-PrediXcan (TWAS)"
echo "  SPrediXcan    : $SPREDIXCAN"
echo "  Model dir     : $MODEL_DIR"
echo "  Covariance dir: $COV_DIR"
echo "  GWAS dir      : $GWAS_DIR"
echo "  Output dir    : $OUT_DIR"
echo "  GWAS pattern  : $GWAS_PATTERN"
echo "  Max jobs      : $JOBS"
echo "  Start time    : $(date)"

#  Collect GWAS files 
shopt -s nullglob
GWAS_FILES=("$GWAS_DIR"/$GWAS_PATTERN)
shopt -u nullglob

if [[ ${#GWAS_FILES[@]} -eq 0 ]]; then
  echo "Error: no GWAS files matching '$GWAS_PATTERN' found in $GWAS_DIR"
  exit 1
fi

echo "Found ${#GWAS_FILES[@]} GWAS file(s)."
echo

#  Job control 
wait_for_jobs() {
  while [[ "$(jobs -rp | wc -l)" -ge "$JOBS" ]]; do
    sleep 1
  done
}

#  Main loop 
for GWAS in "${GWAS_FILES[@]}"; do
  base=$(basename "$GWAS" .summary.txt)

  # Split into model_tag and trait_tag using SEP
  model_tag="${base%%${SEP}*}"
  trait_tag="${base#*${SEP}}"

  model_db="$MODEL_DIR/${MODEL_PREFIX}${model_tag}${MODEL_SUFFIX}"
  cov_file="$COV_DIR/${model_tag}/${model_tag}${COV_SUFFIX}"
  out_file="$OUT_DIR/${model_tag}${SEP}${trait_tag}.csv"

  if [[ ! -f "$model_db" || ! -f "$cov_file" ]]; then
    continue
  fi

  wait_for_jobs

  (
    echo "[$(date '+%H:%M:%S')] start ${model_tag} x ${trait_tag}"
    python3 "$SPREDIXCAN" \
      --model_db_path "$model_db" \
      --covariance "$cov_file" \
      --gwas_file "$GWAS" \
      --snp_column "$SNP_COL" \
      --effect_allele_column "$EA_COL" \
      --non_effect_allele_column "$NEA_COL" \
      --beta_column "$BETA_COL" \
      --pvalue_column "$PVAL_COL" \
      --keep_non_rsid \
      --overwrite \
      --output_file "$out_file"
    echo "[$(date '+%H:%M:%S')] done  ${model_tag} x ${trait_tag}"
  ) &

done

wait
echo "All jobs finished at $(date)"