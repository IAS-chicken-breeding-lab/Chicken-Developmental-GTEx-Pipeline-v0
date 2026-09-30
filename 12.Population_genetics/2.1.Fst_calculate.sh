#!/usr/bin/env bash
#
# Compute windowed Weir & Cockerham FST between two populations from a merged VCF.

set -euo pipefail

# Defaults 
MERGED_VCF=""
POP1_SAMPLES=""
POP2_SAMPLES=""
OUTDIR=""
WIN=50000
STEP=10000
MAF=0.01
MAX_MISSING=0.8
MIN_SNPS_PER_WINDOW=10
THREADS=$(nproc)
KEEP_INTERMEDIATE=false

#  Help 
usage() {
  cat <<EOF
Usage: $0 --vcf <merged.vcf[.gz]> --pop1 <samples1.txt> --pop2 <samples2.txt> --outdir <dir> [options]

Required:
  --vcf        Merged VCF file (bgzipped or plain)
  --pop1       File with sample IDs for population 1 (one per line)
  --pop2       File with sample IDs for population 2 (one per line)
  --outdir     Output directory (will be created)

Optional:
  --win        Window size in bp            (default: ${WIN})
  --step       Window step in bp            (default: ${STEP})
  --maf        Minimum MAF                  (default: ${MAF})
  --max-miss   Maximum missing rate         (default: ${MAX_MISSING})
  --min-snps   Minimum SNPs per window      (default: ${MIN_SNPS_PER_WINDOW})
  --threads    Number of parallel threads   (default: ${THREADS})
  --keep-tmp   Keep intermediate per-chrom files
  -h, --help   Show this help

Dependencies:
  bcftools, vcftools, tabix, bgzip

Example:
  bash $0 \\
    --vcf merged.vcf.gz \\
    --pop1 pop1.samples.txt \\
    --pop2 pop2.samples.txt \\
    --outdir results/ \\
    --win 50000 --step 10000 \\
    --threads 16
EOF
  exit 0
}

# Parse args
while [[ $# -gt 0 ]]; do
  case $1 in
    --vcf)       MERGED_VCF="$2"; shift 2 ;;
    --pop1)      POP1_SAMPLES="$2"; shift 2 ;;
    --pop2)      POP2_SAMPLES="$2"; shift 2 ;;
    --outdir)    OUTDIR="$2"; shift 2 ;;
    --win)       WIN="$2"; shift 2 ;;
    --step)      STEP="$2"; shift 2 ;;
    --maf)       MAF="$2"; shift 2 ;;
    --max-miss)  MAX_MISSING="$2"; shift 2 ;;
    --min-snps)  MIN_SNPS_PER_WINDOW="$2"; shift 2 ;;
    --threads)   THREADS="$2"; shift 2 ;;
    --keep-tmp)  KEEP_INTERMEDIATE=true; shift ;;
    -h|--help)   usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Validate
if [[ -z "$MERGED_VCF" || -z "$POP1_SAMPLES" || -z "$POP2_SAMPLES" || -z "$OUTDIR" ]]; then
  echo "Error: --vcf, --pop1, --pop2 and --outdir are required."
  usage
fi

for f in "$MERGED_VCF" "$POP1_SAMPLES" "$POP2_SAMPLES"; do
  [[ -f "$f" ]] || { echo "Error: file not found: $f"; exit 1; }
done

command -v bcftools >/dev/null || { echo "Error: bcftools not found."; exit 1; }
command -v vcftools >/dev/null || { echo "Error: vcftools not found."; exit 1; }
command -v tabix    >/dev/null || { echo "Error: tabix not found."; exit 1; }

mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")
LOG_DIR="${OUTDIR}/logs"
mkdir -p "$LOG_DIR"

echo "Windowed FST pipeline"
echo "  VCF        : $MERGED_VCF"
echo "  Pop1       : $POP1_SAMPLES"
echo "  Pop2       : $POP2_SAMPLES"
echo "  Window/Step: ${WIN} / ${STEP} bp"
echo "  MAF        : ${MAF}"
echo "  Max missing: ${MAX_MISSING}"
echo "  Min SNPs   : ${MIN_SNPS_PER_WINDOW}"
echo "  Threads    : ${THREADS}"
echo "  Output     : $OUTDIR"

# Prepare bgzipped VCF
if [[ "$MERGED_VCF" == *.gz ]]; then
  VCF_GZ="$MERGED_VCF"
else
  VCF_GZ="${OUTDIR}/$(basename "$MERGED_VCF" .vcf).vcf.gz"
  if [[ ! -s "$VCF_GZ" ]]; then
    echo "[1/5] Compressing VCF with bgzip ..."
    bgzip -@ "$THREADS" -c "$MERGED_VCF" > "$VCF_GZ"
  fi
fi

if [[ ! -f "${VCF_GZ}.tbi" && ! -f "${VCF_GZ%.gz}.tbi" ]]; then
  echo "[1/5] Indexing VCF with tabix ..."
  tabix -p vcf "$VCF_GZ"
fi

# Validate sample lists
echo "[2/5] Validating sample lists ..."
vcf_samples="${OUTDIR}/all_vcf_samples.txt"
bcftools query -l "$VCF_GZ" > "$vcf_samples"

grep -F -x -v -f "$vcf_samples" "$POP1_SAMPLES" > "${LOG_DIR}/missing_pop1.txt" || true
grep -F -x -v -f "$vcf_samples" "$POP2_SAMPLES" > "${LOG_DIR}/missing_pop2.txt" || true

missing_pop1=$(wc -l < "${LOG_DIR}/missing_pop1.txt")
missing_pop2=$(wc -l < "${LOG_DIR}/missing_pop2.txt")

echo "  Pop1 samples in VCF: $(( $(wc -l < "$POP1_SAMPLES") - missing_pop1 ))"
echo "  Pop2 samples in VCF: $(( $(wc -l < "$POP2_SAMPLES") - missing_pop2 ))"

if [[ $missing_pop1 -gt 0 || $missing_pop2 -gt 0 ]]; then
  echo "  WARNING: $missing_pop1 pop1 and $missing_pop2 pop2 samples not in VCF."
  if [[ -t 0 ]]; then
    read -p "  Continue? (y/n): " -n 1 -r
    echo
    [[ $REPLY =~ ^[Yy]$ ]] || exit 1
  else
    echo "  Non-interactive session: continuing anyway."
  fi
fi

# Per-chromosome FST 
echo "[3/5] Computing FST per chromosome ..."
chromosomes=$(bcftools query -f '%CHROM\n' "$VCF_GZ" | sort -u)
if [[ -z "$chromosomes" ]]; then
  echo "Error: could not determine chromosome list from VCF."
  exit 1
fi

chr_count=$(echo "$chromosomes" | wc -l)
echo "  Found $chr_count chromosomes."

run_chr() {
  local chr="$1"
  local clean_chr
  clean_chr=$(echo "$chr" | sed 's/^chr//I')
  local prefix="${OUTDIR}/chr${clean_chr}.${WIN}_${STEP}"
  local fst_file="${prefix}.windowed.weir.fst"
  local log_file="${LOG_DIR}/chr${clean_chr}.log"

  if [[ -s "$fst_file" && $(wc -l < "$fst_file") -gt 1 ]]; then
    echo "  [skip] chr${clean_chr} already done"
    return 0
  fi

  vcftools --gzvcf "$VCF_GZ" \
           --weir-fst-pop "$POP1_SAMPLES" \
           --weir-fst-pop "$POP2_SAMPLES" \
           --chr "$chr" \
           --fst-window-size "$WIN" \
           --fst-window-step "$STEP" \
           --maf "$MAF" \
           --max-missing "$MAX_MISSING" \
           --out "$prefix" >"$log_file" 2>&1

  if [[ -s "$fst_file" ]]; then
    echo "  [done] chr${clean_chr}: $(( $(wc -l < "$fst_file") - 1 )) windows"
  else
    echo "  [fail] chr${clean_chr}: see $log_file"
  fi
}

export -f run_chr
export VCF_GZ POP1_SAMPLES POP2_SAMPLES OUTDIR LOG_DIR WIN STEP MAF MAX_MISSING

echo "$chromosomes" | xargs -P "$THREADS" -I {} bash -c 'run_chr "$@"' _ {}

# Merge results 
echo "[4/5] Merging results ..."
combined_file="${OUTDIR}/all_chromosomes.${WIN}_${STEP}.fst"
echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST" > "$combined_file"

for chr in $chromosomes; do
  clean_chr=$(echo "$chr" | sed 's/^chr//I')
  chr_file="${OUTDIR}/chr${clean_chr}.${WIN}_${STEP}.windowed.weir.fst"
  [[ -s "$chr_file" ]] && tail -n +2 "$chr_file"
done | sort -k1,1V -k2,2n >> "$combined_file"

total_windows=$(( $(wc -l < "$combined_file") - 1 ))
echo "  Total windows: $total_windows"

# Filter and top 1% 
echo "[5/5] Filtering and selecting top 1% ..."
filtered_file="${combined_file}.filter"
top1_file="${filtered_file}.top1pct"

awk -v min_snps="$MIN_SNPS_PER_WINDOW" 'BEGIN{OFS="\t"} NR==1{print; next} $4>=min_snps{print}' \
  "$combined_file" > "$filtered_file"

filtered=$(( $(wc -l < "$filtered_file") - 1 ))
echo "  Windows after filtering: $filtered"

if [[ $filtered -gt 0 ]]; then
  top=$(( filtered / 100 ))
  [[ $top -eq 0 ]] && top=1

  (head -1 "$filtered_file" && \
   tail -n +2 "$filtered_file" | sort --parallel="$THREADS" -k5,5gr) | \
   head -$((top + 1)) > "$top1_file"

  echo "  Top 1% windows: $top"
  echo ""
  echo "  Top 5 windows:"
  head -6 "$top1_file" | column -t
else
  echo "  No windows passed filtering."
fi

# Cleanup
if [[ "$KEEP_INTERMEDIATE" == false ]]; then
  rm -f "${OUTDIR}"/chr*.windowed.weir.fst
  rm -f "${OUTDIR}"/chr*.log
fi

echo ""
echo "Done. Results in: $OUTDIR"