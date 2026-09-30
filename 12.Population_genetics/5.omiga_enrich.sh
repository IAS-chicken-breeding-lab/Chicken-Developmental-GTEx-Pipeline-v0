#!/usr/bin/env bash

# Batch-run OmiGA enrichment for all tissues listed in a file.
#
# Usage:
#   bash run_omiga_enrich.sh \
#     --genotype /path/to/genotype_prefix \
#     --annot /path/to/annotation.bed \
#     --workdir /path/to/target_dir \
#     --outdir /path/to/output_dir \
#     --tissues /path/to/tissues.txt
#

set -euo pipefail

OMIGA="OmiGA"
GENOTYPE=""
ANNOT=""
WORKDIR=""
OUTDIR=""
TISSUES_FILE=""
PERMUTATIONS=1000
MAF_MATCH=0.02

usage() {
  cat <<EOF
Usage: $0 --genotype <prefix> --annot <bed> --workdir <dir> --outdir <dir> --tissues <file> [options]

Required:
  --genotype   PLINK binary prefix (without .bed/.bim/.fam)
  --annot      Annotation BED file
  --workdir    Directory containing <tissue>_bimid.txt files
  --outdir     Output directory
  --tissues    File with tissue names (one per line)

Optional:
  --omiga        Path to OmiGA executable (default: OmiGA)
  --permutations Number of permutations (default: 1000)
  --maf-match    MAF matching threshold (default: 0.02)
  -h, --help     Show this help
EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case $1 in
    --genotype)     GENOTYPE="$2"; shift 2 ;;
    --annot)        ANNOT="$2"; shift 2 ;;
    --workdir)      WORKDIR="$2"; shift 2 ;;
    --outdir)       OUTDIR="$2"; shift 2 ;;
    --tissues)      TISSUES_FILE="$2"; shift 2 ;;
    --omiga)        OMIGA="$2"; shift 2 ;;
    --permutations) PERMUTATIONS="$2"; shift 2 ;;
    --maf-match)    MAF_MATCH="$2"; shift 2 ;;
    -h|--help)      usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$GENOTYPE" || -z "$ANNOT" || -z "$WORKDIR" || -z "$OUTDIR" || -z "$TISSUES_FILE" ]]; then
  echo "Error: missing required arguments."
  usage
fi

if [[ ! -f "$TISSUES_FILE" ]]; then
  echo "Error: tissues file not found: $TISSUES_FILE"
  exit 1
fi

mapfile -t tissues < "$TISSUES_FILE"

echo "=== OmiGA enrichment batch ==="
echo "Genotype  : $GENOTYPE"
echo "Annotation: $ANNOT"
echo "Workdir   : $WORKDIR"
echo "Outdir    : $OUTDIR"
echo "Tissues   : ${tissues[*]}"
echo

for tissue in "${tissues[@]}"; do
  target_file="${WORKDIR}/${tissue}.txt"

  if [[ ! -f "$target_file" ]]; then
    echo "Warning: $target_file not found, skipping $tissue"
    continue
  fi

  echo "Processing: $tissue"

  bash omiga_enrich.sh \
    --tissue "$tissue" \
    --target "$target_file" \
    --genotype "$GENOTYPE" \
    --annot "$ANNOT" \
    --outdir "$OUTDIR" \
    --omiga "$OMIGA" \
    --permutations "$PERMUTATIONS" \
    --maf-match "$MAF_MATCH"

  echo
done

echo "All done. Results in: $OUTDIR"