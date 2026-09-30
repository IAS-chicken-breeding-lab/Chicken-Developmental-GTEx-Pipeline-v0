#!/usr/bin/env bash
#
# combine_twas_models.sh
# Combine per-gene TWAS model outputs (summary, covariances, weights)
# into per-tissue, per-chromosome files.
#
# USAGE
#   bash combine_twas_models.sh \
#     --src-dir <dir> \
#     --dst-dir <dir> \
#     [--chromosomes <list>] \
#     [--categories <list>]
#
# REQUIRED ARGUMENTS
#   --src-dir     Source directory containing per-tissue, per-chromosome
#                 subdirectories under: summary/, covariances/, weights/
#   --dst-dir     Output directory for combined files
#
# OPTIONAL ARGUMENTS
#   --chromosomes  Space-separated chromosome list (default: 1..39)
#   --categories   Space-separated categories (default: "summary covariances weights")
#   -h, --help     Show this help
#
# INPUT STRUCTURE
#   <src-dir>/
#     ├── summary/<tissue>/<chr>/<tissue>_Model_training_chr<chr>.*.model_summaries.txt
#     ├── covariances/<tissue>/<chr>/<tissue>_Model_training_chr<chr>.*.covariances.txt
#     └── weights/<tissue>/<chr>/<tissue>_Model_training_chr<chr>.*.weights.txt
#

set -euo pipefail

#  Defaults 
SRC_BASE=""
DST_BASE=""
CHROMOSOMES=""
CATEGORIES="summary covariances weights"

usage() {
  sed -n '3,60p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

#  Parse arguments 
while [[ $# -gt 0 ]]; do
  case $1 in
    --src-dir)      SRC_BASE="$2";     shift 2 ;;
    --dst-dir)      DST_BASE="$2";     shift 2 ;;
    --chromosomes)  CHROMOSOMES="$2";  shift 2 ;;
    --categories)   CATEGORIES="$2";   shift 2 ;;
    -h|--help)      usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

#  Validate 
if [[ -z "$SRC_BASE" || -z "$DST_BASE" ]]; then
  echo "Error: --src-dir and --dst-dir are required."
  usage
fi

if [[ ! -d "$SRC_BASE" ]]; then
  echo "Error: source directory not found: $SRC_BASE"
  exit 1
fi

mkdir -p "$DST_BASE"

#  Chromosome list 
if [[ -n "$CHROMOSOMES" ]]; then
  read -r -a chr_array <<< "$CHROMOSOMES"
else
  chr_array=($(seq 1 39))
fi

#  Category list 
read -r -a cat_array <<< "$CATEGORIES"

echo "Combine TWAS model outputs"
echo "  Source dir    : $SRC_BASE"
echo "  Output dir    : $DST_BASE"
echo "  Chromosomes   : ${chr_array[*]}"
echo "  Categories    : ${cat_array[*]}"
echo "  Start time    : $(date)"

#  Discover tissues 
# Tissues are inferred from the first category that exists in src-dir.
tissue_array=()
for cat in "${cat_array[@]}"; do
  if [[ -d "$SRC_BASE/$cat" ]]; then
    for tdir in "$SRC_BASE/$cat"/*/; do
      [[ -d "$tdir" ]] || continue
      tissue_array+=("$(basename "$tdir")")
    done
    break
  fi
done

if [[ ${#tissue_array[@]} -eq 0 ]]; then
  echo "Error: no tissues found under $SRC_BASE"
  exit 1
fi

# Deduplicate
readarray -t tissue_array < <(printf '%s\n' "${tissue_array[@]}" | sort -u)

echo "Found ${#tissue_array[@]} tissue(s): ${tissue_array[*]}"
echo

#  Worker 
process_tissue() {
  local tissue="$1"
  echo ">>> Processing tissue: $tissue"
  local start_time=$(date +%s)

  # Create output directories
  for cat in "${cat_array[@]}"; do
    mkdir -p "$DST_BASE/$cat/$tissue"
  done

  for chr in "${chr_array[@]}"; do

    ## ----- 1. summary -----
    if [[ " ${cat_array[*]} " == *" summary "* ]]; then
      local summary_src_dir="$SRC_BASE/summary/$tissue/$chr"
      if [[ -d "$summary_src_dir" ]]; then
        local model_files=("$summary_src_dir"/${tissue}_Model_training_chr${chr}.*.model_summaries.txt)
        if [[ ${#model_files[@]} -gt 0 && -f "${model_files[0]}" ]]; then
          local summary_out="$DST_BASE/summary/$tissue/${tissue}_Model_training_chr${chr}_model_summaries.txt"
          echo "    [summary] chr$chr -> $(basename "$summary_out") (${#model_files[@]} files)"
          {
            printf "gene_id\tgene_name\tgene_type\talpha\tn_snps_in_window\tn_snps_in_model\tlambda_min_mse\ttest_R2_avg\ttest_R2_sd\tcv_R2_avg\tcv_R2_sd\tin_sample_R2\tnested_cv_fisher_pval\trho_avg\trho_se\trho_zscore\trho_avg_squared\tzscore_pval\tcv_rho_avg\tcv_rho_se\tcv_rho_avg_squared\tcv_zscore_est\tcv_zscore_pval\tcv_pval_est\n"
            for f in "${model_files[@]}"; do
              tail -n +2 "$f"
            done
          } > "$summary_out"
        fi

        local summary_txt="$summary_src_dir/${tissue}_Model_training_chr${chr}_summary.txt"
        if [[ -f "$summary_txt" ]]; then
          cp -f "$summary_txt" "$DST_BASE/summary/$tissue/"
        fi
      fi
    fi

    ## 2. covariances
    if [[ " ${cat_array[*]} " == *" covariances "* ]]; then
      local cov_src_dir="$SRC_BASE/covariances/$tissue/$chr"
      if [[ -d "$cov_src_dir" ]]; then
        local cov_files=("$cov_src_dir"/${tissue}_Model_training_chr${chr}.*.covariances.txt)
        if [[ ${#cov_files[@]} -gt 0 && -f "${cov_files[0]}" ]]; then
          local cov_out="$DST_BASE/covariances/$tissue/${tissue}_Model_training_chr${chr}_covariances.txt"
          echo "    [covariances] chr$chr -> $(basename "$cov_out") (${#cov_files[@]} files)"
          {
            printf "GENE RSID1 RSID2 VALUE\n"
            for f in "${cov_files[@]}"; do
              tail -n +2 "$f"
            done
          } > "$cov_out"
        fi
      fi
    fi

    ## 3. weights
    if [[ " ${cat_array[*]} " == *" weights "* ]]; then
      local weights_src_dir="$SRC_BASE/weights/$tissue/$chr"
      if [[ -d "$weights_src_dir" ]]; then
        local weights_files=("$weights_src_dir"/${tissue}_Model_training_chr${chr}.*.weights.txt)
        if [[ ${#weights_files[@]} -gt 0 && -f "${weights_files[0]}" ]]; then
          local weights_out="$DST_BASE/weights/$tissue/${tissue}_Model_training_chr${chr}_weights.txt"
          echo "    [weights] chr$chr -> $(basename "$weights_out") (${#weights_files[@]} files)"
          {
            printf "gene_id\trsid\tvarID\tref\talt\tbeta\n"
            for f in "${weights_files[@]}"; do
              tail -n +2 "$f"
            done
          } > "$weights_out"
        fi
      fi
    fi

  done

  local end_time=$(date +%s)
  echo "<<< Finished tissue: $tissue (${end_time}-${start_time}s)"
}

#  Main loop 
for tissue in "${tissue_array[@]}"; do
  process_tissue "$tissue"
done

#  Summary 
echo "All done!"
for cat in "${cat_array[@]}"; do
  if [[ -d "$DST_BASE/$cat" ]]; then
    count=$(find "$DST_BASE/$cat" -name "*.txt" | wc -l)
    echo "  $cat: $count files"
  fi
done