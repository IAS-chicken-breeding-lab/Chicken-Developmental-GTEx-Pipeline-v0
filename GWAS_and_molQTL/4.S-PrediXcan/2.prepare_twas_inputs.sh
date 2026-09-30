#!/usr/bin/env bash
#
# prepare_twas_inputs.sh
# Prepare SNP annotation and transformed expression matrices for TWAS
# prediction models, then submit per-chromosome jobs for each tissue.
#
# USAGE
#   bash prepare_twas_inputs.sh \
#     --geno-dir <dir> \
#     --expr-dir <dir> \
#     --out-dir <dir> \
#     --chr-job-script <path> \
#     [--geno-pattern <glob>] \
#     [--expr-pattern <glob>] \
#     [--chromosomes <list>] \
#     [--chr-set <n>] \
#     [--type <string>]
#
# REQUIRED ARGUMENTS
#   --geno-dir       Directory with PLINK binary genotype files.
#                    Expected per-tissue files: <tissue>_qc.bed/.bim/.fam
#   --expr-dir       Directory with expression phenotype files.
#                    Expected per-tissue file: <tissue>_pheno_autosome.bed
#   --out-dir        Output directory (created if missing)
#   --chr-job-script Path to the per-chromosome job script
#                    (called as: bash <script> <chr> <tissue>)
#
# OPTIONAL ARGUMENTS
#   --geno-pattern   Glob pattern for genotype files (default: *_qc.bed)
#   --expr-pattern   Expression file pattern (default: _pheno_autosome.bed)
#   --chromosomes    Space-separated chromosome list (default: 1..39)
#   --chr-set        PLINK --chr-set value (default: 39)
#   --type           Label used in output filenames (default: eQTL)
#   -h, --help       Show this help
#


set -euo pipefail

#  Defaults 
GENO_DIR=""
EXPR_DIR=""
OUT_DIR=""
CHR_JOB_SCRIPT=""
GENO_PATTERN="*_qc.bed"
EXPR_PATTERN="_pheno_autosome.bed"
CHROMOSOMES=""
CHR_SET=39
TYPE="eQTL"

usage() {
  sed -n '3,60p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

#  Parse arguments 
while [[ $# -gt 0 ]]; do
  case $1 in
    --geno-dir)       GENO_DIR="$2";       shift 2 ;;
    --expr-dir)       EXPR_DIR="$2";       shift 2 ;;
    --out-dir)        OUT_DIR="$2";        shift 2 ;;
    --chr-job-script) CHR_JOB_SCRIPT="$2"; shift 2 ;;
    --geno-pattern)   GENO_PATTERN="$2";   shift 2 ;;
    --expr-pattern)   EXPR_PATTERN="$2";   shift 2 ;;
    --chromosomes)    CHROMOSOMES="$2";    shift 2 ;;
    --chr-set)        CHR_SET="$2";        shift 2 ;;
    --type)           TYPE="$2";           shift 2 ;;
    -h|--help)        usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

#  Validate 
if [[ -z "$GENO_DIR" || -z "$EXPR_DIR" || -z "$OUT_DIR" || -z "$CHR_JOB_SCRIPT" ]]; then
  echo "Error: --geno-dir, --expr-dir, --out-dir and --chr-job-script are required."
  usage
fi

for d in "$GENO_DIR" "$EXPR_DIR"; do
  [[ -d "$d" ]] || { echo "Error: directory not found: $d"; exit 1; }
done

[[ -f "$CHR_JOB_SCRIPT" ]] || { echo "Error: chromosome job script not found: $CHR_JOB_SCRIPT"; exit 1; }

command -v plink >/dev/null 2>&1 || { echo "Error: plink not found in PATH."; exit 1; }
command -v awk   >/dev/null 2>&1 || { echo "Error: awk not found."; exit 1; }

mkdir -p "$OUT_DIR" "$OUT_DIR/logs" "$OUT_DIR/tissue_snp_annotation" "$OUT_DIR/tissue_expression"
cd "$OUT_DIR"

#  Chromosome list 
if [[ -n "$CHROMOSOMES" ]]; then
  read -r -a chr_array <<< "$CHROMOSOMES"
else
  chr_array=($(seq 1 "$CHR_SET"))
fi

echo "=============================================="
echo "TWAS input preparation"
echo "  Genotype dir    : $GENO_DIR"
echo "  Expression dir  : $EXPR_DIR"
echo "  Output dir      : $OUT_DIR"
echo "  Chr job script  : $CHR_JOB_SCRIPT"
echo "  Geno pattern    : $GENO_PATTERN"
echo "  Expr pattern    : $EXPR_PATTERN"
echo "  Chromosomes     : ${chr_array[*]}"
echo "  chr-set         : $CHR_SET"
echo "  Type label      : $TYPE"
echo "=============================================="

#  Discover tissues 
shopt -s nullglob
geno_files=("$GENO_DIR"/$GENO_PATTERN)
shopt -u nullglob

if [[ ${#geno_files[@]} -eq 0 ]]; then
  echo "Error: no genotype files matching '$GENO_PATTERN' in $GENO_DIR"
  exit 1
fi

echo "Found ${#geno_files[@]} genotype file(s)."

#  Main loop 
for bed in "${geno_files[@]}"; do
  tissue=$(basename "$bed" _qc.bed)
  echo "------ Processing tissue: $tissue ------"

  #  1. SNP annotation 
  snp_dir="$OUT_DIR/tissue_snp_annotation/$tissue"
  mkdir -p "$snp_dir"

  plink --allow-extra-chr --chr-set "$CHR_SET" \
        --bfile "$GENO_DIR/${tissue}_qc" \
        --make-just-bim \
        --out "$snp_dir/$tissue"

  awk '{print $1"\t"$4"\t"$2"\t"$5"\t"$6}' \
      "$snp_dir/${tissue}.bim" \
      > "$snp_dir/${tissue}.snp_annotation.txt"

  #  2. Expression matrix 
  expr_file="$EXPR_DIR/${tissue}${EXPR_PATTERN}"
  if [[ ! -f "$expr_file" ]]; then
    echo "  [WARN] Expression file not found: $expr_file, skipping."
    continue
  fi

  expr_out_dir="$OUT_DIR/tissue_expression"
  mkdir -p "$expr_out_dir"

  tmp_expr="$expr_out_dir/${tissue}.expr.tmp1.txt"
  tmp_samples="$expr_out_dir/${tissue}.sample_ids.txt"
  final_expr="$expr_out_dir/${tissue}.${TYPE}.transformed_expression.txt"

  # gene x sample (skip header)
  tail -n +2 "$expr_file" | awk '{printf "%s", $4; for(i=5;i<=NF;i++){printf "\t%s",$i}; printf "\n"}' > "$tmp_expr"

  # sample IDs from header
  head -n 1 "$expr_file" | cut -f5- | tr '\t' '\n' > "$tmp_samples"

  # transpose to sample x gene
  awk '
  BEGIN {
      sample_count = 0
      while ((getline < ARGV[1]) > 0) {
          sample_count++
          samples[sample_count] = $0
      }
      close(ARGV[1])

      gene_count = 0
      while ((getline < ARGV[2]) > 0) {
          gene_count++
          split($0, fields, "\t")
          gene_ids[gene_count] = fields[1]
          for (i = 2; i <= length(fields); i++) {
              expr_matrix[gene_count, i-1] = fields[i]
          }
      }
      close(ARGV[2])

      printf "id"
      for (g = 1; g <= gene_count; g++) printf "\t%s", gene_ids[g]
      printf "\n"

      for (s = 1; s <= sample_count; s++) {
          printf "%s", samples[s]
          for (g = 1; g <= gene_count; g++) printf "\t%s", expr_matrix[g, s]
          printf "\n"
      }
  }' "$tmp_samples" "$tmp_expr" > "$final_expr"

  rm -f "$tmp_expr" "$tmp_samples"

  echo "  Expression matrix written: $final_expr"

  #  3. Submit per-chromosome jobs 
  for chr in "${chr_array[@]}"; do
    bash "$CHR_JOB_SCRIPT" "$chr" "$tissue" \
      > "$OUT_DIR/logs/${tissue}_chr${chr}.log" 2>&1 &
  done

  echo "  Submitted ${#chr_array[@]} chromosome jobs for $tissue"
done

wait
echo "====== All jobs submitted and completed ======"