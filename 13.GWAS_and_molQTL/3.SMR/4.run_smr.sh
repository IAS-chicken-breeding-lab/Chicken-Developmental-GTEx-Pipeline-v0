#!/usr/bin/env bash
#
# Batch SMR analysis: run SMR for each BESD file against every GWAS file
# 
# USAGE
#   bash run_smr.sh \
#     --smr <path> \
#     --besd-dir <dir> \
#     --gwas-dir <dir> \
#     --geno-dir <dir> \
#     --out-dir <dir> \
#     [--gwas-pattern <glob>] \
#     [--besd-pattern <glob>] \
#     [--threads <n>] \
#     [--peqtl <float>] \
#     [--diff-freq <float>] \
#     [--diff-freq-prop <float>] \
#     [--compress]
#
# OPTIONAL ARGUMENTS
#   --gwas-pattern    Glob pattern for GWAS files (default: "*.ma")
#   --besd-pattern    Glob pattern for BESD files (default: "*.besd")
#   --threads         Number of threads for SMR (default: 8)
#   --peqtl           eQTL p-value threshold for SMR (default: 1e-5)
#   --diff-freq       Allele frequency difference threshold (default: 0.9)
#   --diff-freq-prop  Allele frequency difference proportion (default: 0.1)
#   --compress        Compress .smr output with pigz if available
#   -h, --help        Show this help
#
# EXAMPLE
#   bash run_smr_batch.sh \
#     --smr /path/to/smr \
#     --besd-dir /path/to/besd \
#     --gwas-dir /path/to/gwas \
#     --geno-dir /path/to/geno \
#     --out-dir /path/to/results \
#     --threads 8 \
#     --compress
#

set -euo pipefail

#  Defaults 
SMR_BIN=""
BESD_DIR=""
GWAS_DIR=""
GENO_DIR=""
OUT_DIR=""
GWAS_PATTERN="*.ma"
BESD_PATTERN="*.besd"
THREADS=8
P_THRESHOLD=1e-5
DIFF_FREQ=0.9
DIFF_FREQ_PROP=0.1
COMPRESS=false

usage() {
  sed -n '3,70p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

#  Parse arguments 
while [[ $# -gt 0 ]]; do
  case $1 in
    --smr)             SMR_BIN="$2"; shift 2 ;;
    --besd-dir)        BESD_DIR="$2"; shift 2 ;;
    --gwas-dir)        GWAS_DIR="$2"; shift 2 ;;
    --geno-dir)        GENO_DIR="$2"; shift 2 ;;
    --out-dir)         OUT_DIR="$2"; shift 2 ;;
    --gwas-pattern)    GWAS_PATTERN="$2"; shift 2 ;;
    --besd-pattern)    BESD_PATTERN="$2"; shift 2 ;;
    --threads)         THREADS="$2"; shift 2 ;;
    --peqtl)           P_THRESHOLD="$2"; shift 2 ;;
    --diff-freq)       DIFF_FREQ="$2"; shift 2 ;;
    --diff-freq-prop)  DIFF_FREQ_PROP="$2"; shift 2 ;;
    --compress)        COMPRESS=true; shift ;;
    -h|--help)         usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

#  Validate 
if [[ -z "$SMR_BIN" || -z "$BESD_DIR" || -z "$GWAS_DIR" || -z "$GENO_DIR" || -z "$OUT_DIR" ]]; then
  echo "Error: --smr, --besd-dir, --gwas-dir, --geno-dir and --out-dir are required."
  usage
fi

if [[ ! -x "$SMR_BIN" ]]; then
  echo "Error: SMR executable not found or not executable: $SMR_BIN"
  exit 1
fi

for d in "$BESD_DIR" "$GWAS_DIR" "$GENO_DIR"; do
  if [[ ! -d "$d" ]]; then
    echo "Error: directory not found: $d"
    exit 1
  fi
done

mkdir -p "$OUT_DIR"

#  Collect files 
shopt -s nullglob
besd_files=("$BESD_DIR"/$BESD_PATTERN)
gwas_files=("$GWAS_DIR"/$GWAS_PATTERN)
shopt -u nullglob

if [[ ${#besd_files[@]} -eq 0 ]]; then
  echo "Error: no BESD files matching '$BESD_PATTERN' in $BESD_DIR"
  exit 1
fi
if [[ ${#gwas_files[@]} -eq 0 ]]; then
  echo "Error: no GWAS files matching '$GWAS_PATTERN' in $GWAS_DIR"
  exit 1
fi

echo "Found ${#besd_files[@]} BESD file(s) and ${#gwas_files[@]} GWAS file(s)."
echo

#  Helper: derive genotype prefix 
derive_geno_prefix() {
  local prefix="$1"
  # Replace "_eGene_ciseQTL" with "_qc" to get the genotype prefix.
  # Users with different conventions should edit this line.
  echo "${prefix}" | sed 's/_eGene_ciseQTL/_qc/'
}

#  Main loop 
for besd_file in "${besd_files[@]}"; do
  prefix=$(basename "$besd_file" .besd)
  echo "Processing BESD: ${prefix}"

  geno_prefix=$(derive_geno_prefix "$prefix")
  geno_path="${GENO_DIR}/${geno_prefix}"

  if [[ ! -f "${geno_path}.bed" ]]; then
    echo "  [WARN] No genotype file found for ${prefix} (expected ${geno_path}.bed), skipping."
    continue
  fi

  mkdir -p "${OUT_DIR}/${prefix}"

  for gwas_file in "${gwas_files[@]}"; do
    trait=$(basename "$gwas_file")
    trait="${trait%.ma}"     # strip .ma
    out_prefix="${OUT_DIR}/${prefix}/${prefix}_${trait}"

    echo "  -> SMR for trait: ${trait}"

    "$SMR_BIN" \
      --bfile "${geno_path}" \
      --gwas-summary "${gwas_file}" \
      --beqtl-summary "${BESD_DIR}/${prefix}" \
      --peqtl-smr "${P_THRESHOLD}" \
      --diff-freq "${DIFF_FREQ}" \
      --diff-freq-prop "${DIFF_FREQ_PROP}" \
      --out "${out_prefix}" \
      --thread-num "${THREADS}"

    if [[ "$COMPRESS" == true && -f "${out_prefix}.smr" ]]; then
      if command -v pigz >/dev/null 2>&1; then
        pigz -f "${out_prefix}.smr"
      else
        gzip -f "${out_prefix}.smr"
      fi
    fi
  done

  echo "Finished BESD: ${prefix}"
done

echo "All SMR analyses completed. Results saved to: ${OUT_DIR}"