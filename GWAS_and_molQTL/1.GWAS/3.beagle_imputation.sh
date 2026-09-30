#!/usr/bin/env bash

set -euo pipefail

# Fill missing genotypes, perform reference-panel imputation, and create PLINK files.
# Usage: bash scripts/03_beagle_imputation.sh [config/config.sh]

CONFIG_FILE="${1:-config/config.sh}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: configuration file not found: ${CONFIG_FILE}" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "${CONFIG_FILE}"

required_vars=(
    JAVA_BIN PLINK_BIN BCFTOOLS_BIN SELF_BEAGLE_JAR REF_BEAGLE_JAR
    REFERENCE_PANEL_VCF QC_OUTPUT_PREFIX IMPUTATION_WORKDIR
    SELF_IMPUTED_PREFIX REF_IMPUTED_PREFIX IMPUTED_PLINK_PREFIX
    JAVA_MEMORY BEAGLE_NE BEAGLE_THREADS ANALYSIS_CHROMOSOMES
)

for var_name in "${required_vars[@]}"; do
    if [[ -z "${!var_name:-}" ]]; then
        echo "ERROR: ${var_name} is not set in ${CONFIG_FILE}" >&2
        exit 1
    fi
done

mkdir -p "${IMPUTATION_WORKDIR}"

INPUT_VCF="${QC_OUTPUT_PREFIX}-1_autosomes.vcf.gz"

echo "[1/5] Fill missing genotypes within the SNP-array dataset"
"${JAVA_BIN}" \
    "-Xmx${JAVA_MEMORY}" \
    -jar "${SELF_BEAGLE_JAR}" \
    gt="${INPUT_VCF}" \
    impute=true \
    out="${SELF_IMPUTED_PREFIX}"
"${BCFTOOLS_BIN}" index -f "${SELF_IMPUTED_PREFIX}.vcf.gz"

echo "[2/5] Perform reference-panel imputation"
"${JAVA_BIN}" \
    "-Xmx${JAVA_MEMORY}" \
    -jar "${REF_BEAGLE_JAR}" \
    gt="${SELF_IMPUTED_PREFIX}.vcf.gz" \
    ref="${REFERENCE_PANEL_VCF}" \
    impute=true \
    ne="${BEAGLE_NE}" \
    nthreads="${BEAGLE_THREADS}" \
    out="${REF_IMPUTED_PREFIX}"

echo "[3/5] Index the imputed VCF"
"${BCFTOOLS_BIN}" index -f "${REF_IMPUTED_PREFIX}.vcf.gz"

echo "[4/5] Convert the imputed VCF to PED/MAP"
"${PLINK_BIN}" \
    --allow-extra-chr \
    --chr-set 39 \
    --vcf "${REF_IMPUTED_PREFIX}.vcf.gz" \
    --chr "${ANALYSIS_CHROMOSOMES}" \
    --recode \
    --out "${IMPUTED_PLINK_PREFIX}"

echo "[5/5] Generate the final imputed PLINK binary dataset"
"${PLINK_BIN}" \
    --file "${IMPUTED_PLINK_PREFIX}" \
    --make-bed \
    --out "${IMPUTED_PLINK_PREFIX}" \
    --chr-set 39 \
    --allow-extra-chr \
    --chr "${ANALYSIS_CHROMOSOMES}"

echo "Beagle imputation completed: ${REF_IMPUTED_PREFIX}.vcf.gz"
echo "Imputed PLINK dataset: ${IMPUTED_PLINK_PREFIX}.bed/.bim/.fam"
