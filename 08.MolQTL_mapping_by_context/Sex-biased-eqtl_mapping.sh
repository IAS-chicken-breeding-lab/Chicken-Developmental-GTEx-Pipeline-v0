#!/bin/bash
set -uo pipefail

# Sex-biased eQTL interaction mapping using OmiGA
#
# Usage:
# bash sb_eqtl_interaction_mapping.sh \
#   Geno_dir Pheno_dir Interaction_dir Output_dir OmiGA [threads]

GENO_DIR=$1
PHENO_DIR=$2
INTERACTION_DIR=$3
OUTPUT_DIR=$4
OMIGA=$5
THREADS=${6:-32}

TISSUES=(
  AF CT FZ GZ HC JG KC MC
  PZ SWM SZ XJ XQN XW XZ ZC
)

PERIODS=(A B C)

mkdir -p "$OUTPUT_DIR"

if [[ ! -x "$OMIGA" ]]; then
  echo "Error: OmiGA executable not found: $OMIGA"
  exit 1
fi


for PERIOD in "${PERIODS[@]}"; do

  for TISSUE in "${TISSUES[@]}"; do

    PREFIX="${PERIOD}_${TISSUE}"

    GENO_PREFIX="${GENO_DIR}/${PREFIX}_filtered"
    PHENO_FILE="${PHENO_DIR}/${PREFIX}_sbGene.bed"
    INTERACTION_FILE="${INTERACTION_DIR}/${PREFIX}_interaction.txt"

    OUTDIR="${OUTPUT_DIR}/${PERIOD}/${TISSUE}"
    mkdir -p "$OUTDIR"

    echo "Processing ${PREFIX} ..."

    # Check input files
    if [[ ! -f "${GENO_PREFIX}.bed" ||
          ! -f "${GENO_PREFIX}.bim" ||
          ! -f "${GENO_PREFIX}.fam" ]]; then

      echo "  Missing genotype files. Skip."
      continue
    fi

    if [[ ! -f "$PHENO_FILE" ]]; then
      echo "  Missing phenotype file: $PHENO_FILE"
      continue
    fi

    if [[ ! -f "$INTERACTION_FILE" ]]; then
      echo "  Missing interaction file: $INTERACTION_FILE"
      continue
    fi

    # OmiGA cis interaction mapping
    if "$OMIGA" \
      --mode cis_interaction \
      --genotype "$GENO_PREFIX" \
      --phenotype "$PHENO_FILE" \
      --interaction "$INTERACTION_FILE" \
      --dprop-pc-covar 0.001 \
      --rm-collinear-covar 0.95 \
      --prefix "$PREFIX" \
      --output-dir "$OUTDIR" \
      --thread "$THREADS"
    then
      echo "  Completed: ${PREFIX}"
    else
      echo "  Failed: ${PREFIX}"
    fi

  done
done

echo "All OmiGA interaction mapping tasks completed."