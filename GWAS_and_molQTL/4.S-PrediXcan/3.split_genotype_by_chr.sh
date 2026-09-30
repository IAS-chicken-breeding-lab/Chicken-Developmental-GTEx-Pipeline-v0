#!/usr/bin/env bash
#
# split_genotype_by_chr.sh
# Split genotype data by chromosome, recode to additive format,
# and transpose to a SNP x sample matrix for TWAS prediction models.
#
# USAGE
#   bash split_genotype_by_chr.sh \
#     --geno-dir <dir> \
#     --snp-annot-dir <dir> \
#     --out-dir <dir> \
#     [--tissues <list or file>] \
#     [--chromosomes <list>] \
#     [--chr-set <n>] \
#     [--jobs <n>]
#
# REQUIRED ARGUMENTS
#   --geno-dir       Directory with PLINK binary genotype files.
#                    Expected per-tissue file: <tissue>_qc.bed/.bim/.fam
#   --snp-annot-dir  Directory containing per-tissue SNP annotation files,
#                    expected at: <snp-annot-dir>/<tissue>/<tissue>.snp_annotation.txt
#   --out-dir        Output directory (created if missing)
#
# OPTIONAL ARGUMENTS
#   --tissues        Comma-separated tissue list, or path to a file with one
#                    tissue per line. If omitted, all tissues matching
#                    <geno-dir>/*_qc.bed are processed.
#   --chromosomes    Space-separated chromosome list (default: 1..39)
#   --chr-set        PLINK --chr-set value (default: 39)
#   --jobs           Max parallel chromosome jobs per tissue (default: 10)
#   -h, --help       Show this help
#
#


set -euo pipefail

#  Defaults 
GENO_DIR=""
SNP_ANNOT_DIR=""
OUT_DIR=""
TISSUES=""
CHROMOSOMES=""
CHR_SET=39
JOBS=10

usage() {
  sed -n '3,60p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

#  Parse arguments 
while [[ $# -gt 0 ]]; do
  case $1 in
    --geno-dir)       GENO_DIR="$2";       shift 2 ;;
    --snp-annot-dir)  SNP_ANNOT_DIR="$2";  shift 2 ;;
    --out-dir)        OUT_DIR="$2";        shift 2 ;;
    --tissues)        TISSUES="$2";        shift 2 ;;
    --chromosomes)    CHROMOSOMES="$2";    shift 2 ;;
    --chr-set)        CHR_SET="$2";        shift 2 ;;
    --jobs)           JOBS="$2";           shift 2 ;;
    -h|--help)        usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

#  Validate 
if [[ -z "$GENO_DIR" || -z "$SNP_ANNOT_DIR" || -z "$OUT_DIR" ]]; then
  echo "Error: --geno-dir, --snp-annot-dir and --out-dir are required."
  usage
fi

for d in "$GENO_DIR" "$SNP_ANNOT_DIR"; do
  [[ -d "$d" ]] || { echo "Error: directory not found: $d"; exit 1; }
done

command -v plink >/dev/null 2>&1 || { echo "Error: plink not found in PATH."; exit 1; }
command -v awk   >/dev/null 2>&1 || { echo "Error: awk not found."; exit 1; }

mkdir -p "$OUT_DIR"
OUT_DIR=$(realpath "$OUT_DIR")

# Output subdirectories
SNP_OUT_DIR="$OUT_DIR/tissue_snp_annotation"
RAW_DIR="$OUT_DIR/raw_tissue_gene_data/tissue_chr_data"
MATRIX_DIR="$OUT_DIR/tissue_gene_matrix"

mkdir -p "$SNP_OUT_DIR" "$RAW_DIR" "$MATRIX_DIR"

#  Chromosome list 
if [[ -n "$CHROMOSOMES" ]]; then
  read -r -a chr_array <<< "$CHROMOSOMES"
else
  chr_array=($(seq 1 "$CHR_SET"))
fi

#  Tissue list 
if [[ -n "$TISSUES" ]]; then
  if [[ -f "$TISSUES" ]]; then
    mapfile -t tissue_array < "$TISSUES"
  else
    IFS=',' read -r -a tissue_array <<< "$TISSUES"
  fi
else
  shopt -s nullglob
  geno_files=("$GENO_DIR"/*_qc.bed)
  shopt -u nullglob
  tissue_array=()
  for bed in "${geno_files[@]}"; do
    tissue_array+=("$(basename "$bed" _qc.bed)")
  done
fi

if [[ ${#tissue_array[@]} -eq 0 ]]; then
  echo "Error: no tissues to process."
  exit 1
fi

echo "Split genotype data by chromosome"
echo "  Genotype dir   : $GENO_DIR"
echo "  SNP annot dir  : $SNP_ANNOT_DIR"
echo "  Output dir     : $OUT_DIR"
echo "  Tissues        : ${tissue_array[*]}"
echo "  Chromosomes    : ${chr_array[*]}"
echo "  chr-set        : $CHR_SET"
echo "  Parallel jobs  : $JOBS"

#  Worker: one chromosome 
process_chr() {
  local tissue="$1"
  local chr="$2"

  local snp_annot_file="$SNP_ANNOT_DIR/$tissue/$tissue.snp_annotation.txt"
  local snp_chr_file="$SNP_OUT_DIR/$tissue/$tissue.snp_annot.chr${chr}.txt"
  local plink_prefix="$RAW_DIR/$tissue/$tissue.ref.chr${chr}"
  local raw_file="${plink_prefix}.raw"
  local matrix_file="$MATRIX_DIR/$tissue/$tissue.genotype.chr${chr}.txt"
  local tmp_matrix_file="$MATRIX_DIR/$tissue/${tissue}.genotype_chr${chr}.txt"

  # 1. Split SNP annotation by chromosome
  awk -v chr="$chr" '$1==chr {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' \
      "$snp_annot_file" > "$snp_chr_file"
  sed -i '1i\chromosome\tpos\tvarID\tref_vcf\talt_vcf' "$snp_chr_file"

  # 2. Extract genotypes for this chromosome with PLINK
  plink --allow-extra-chr --chr-set "$CHR_SET" \
        --bfile "$GENO_DIR/${tissue}_qc" \
        --chr "$chr" \
        --recode A \
        --out "$plink_prefix" > /dev/null 2>&1

  if [[ ! -f "$raw_file" ]]; then
    echo "[WARN] PLINK failed for $tissue chr$chr"
    return
  fi

  # 3. Keep IID and genotype columns only
  cut -d" " -f2,7- "$raw_file" > "$tmp_matrix_file"

  # 4. Clean column names: remove allele suffix and rename IID -> Id
  sed -i '1s/_[A-Z] / /g'  "$tmp_matrix_file"
  sed -i '1s/_[A-Z]$//'    "$tmp_matrix_file"
  sed -i '1s/IID/Id/'      "$tmp_matrix_file"

  # 5. Transpose to SNP x sample
  awk '{
    for (i = 1; i <= NF; i++) col[i] = col[i] $i "\t"
  } END {
    for (i = 1; i <= NF; i++) { sub(/\t$/, "", col[i]); print col[i] }
  }' "$tmp_matrix_file" | sed 's/[ \t]*$//g' > "$matrix_file"

  rm -f "$tmp_matrix_file"

  echo "[DONE] $tissue chr$chr"
}

export -f process_chr
export GENO_DIR SNP_ANNOT_DIR SNP_OUT_DIR RAW_DIR MATRIX_DIR CHR_SET

#  Main loop 
for tissue in "${tissue_array[@]}"; do
  echo "Processing tissue: $tissue "

  # Check input files
  if [[ ! -f "$GENO_DIR/${tissue}_qc.bed" ]]; then
    echo "  [WARN] Genotype file not found: $GENO_DIR/${tissue}_qc.bed, skipping."
    continue
  fi

  snp_annot_file="$SNP_ANNOT_DIR/$tissue/$tissue.snp_annotation.txt"
  if [[ ! -f "$snp_annot_file" ]]; then
    echo "  [WARN] SNP annotation not found: $snp_annot_file, skipping."
    continue
  fi

  mkdir -p "$SNP_OUT_DIR/$tissue" "$RAW_DIR/$tissue" "$MATRIX_DIR/$tissue"

  # Parallel jobs with a cap of $JOBS
  running=0
  for chr in "${chr_array[@]}"; do
    process_chr "$tissue" "$chr" &
    running=$((running + 1))

    if (( running >= JOBS )); then
      wait -n
      running=$((running - 1))
    fi
  done

  wait
  echo " Finished tissue: $tissue "
done

echo " All tissues processed "