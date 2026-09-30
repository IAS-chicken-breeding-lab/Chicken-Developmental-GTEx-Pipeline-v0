#!/usr/bin/env bash

set -euo pipefail

# Perform PLINK quality control and prepare an autosomal VCF for Beagle.
# Usage: bash scripts/02_genotype_qc.sh [config/config.sh]

CONFIG_FILE="${1:-config/config.sh}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: configuration file not found: ${CONFIG_FILE}" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "${CONFIG_FILE}"

required_vars=(
    PLINK_BIN BCFTOOLS_BIN QC_WORKDIR QC_INPUT_PREFIX QC_OUTPUT_PREFIX
    ANALYSIS_CHROMOSOMES BCFTOOLS_REGIONS
)

for var_name in "${required_vars[@]}"; do
    if [[ -z "${!var_name:-}" ]]; then
        echo "ERROR: ${var_name} is not set in ${CONFIG_FILE}" >&2
        exit 1
    fi
done

mkdir -p "${QC_WORKDIR}"
cd "${QC_WORKDIR}"

QC_FILTERED_PREFIX="${QC_OUTPUT_PREFIX}-1"
SORTED_VCF="${QC_FILTERED_PREFIX}_sorted.vcf.gz"
AUTOSOMAL_VCF="${QC_FILTERED_PREFIX}_autosomes.vcf.gz"
CHR_COUNTS="${QC_WORKDIR}/chr_counts.txt"
EXCLUDE_CHR="${QC_WORKDIR}/exclude_chr.txt"

echo "[1/6] Retain biallelic A/C/G/T SNPs on chromosomes ${ANALYSIS_CHROMOSOMES}"
"${PLINK_BIN}" \
    --chr-set 39 \
    --allow-extra-chr \
    --bfile "${QC_INPUT_PREFIX}" \
    --chr "${ANALYSIS_CHROMOSOMES}" \
    --snps-only just-acgt \
    --recode tab \
    --allow-no-sex \
    --make-bed \
    --out "${QC_OUTPUT_PREFIX}"

echo "[2/6] Apply marker and sample quality filters"
# Original thresholds: geno=0.1, mind=0.1, maf=0.01, hwe=0.00001.
"${PLINK_BIN}" \
    --chr-set 39 \
    --allow-extra-chr \
    --bfile "${QC_OUTPUT_PREFIX}" \
    --chr "${ANALYSIS_CHROMOSOMES}" \
    --geno 0.1 \
    --mind 0.1 \
    --maf 0.01 \
    --hwe 0.00001 \
    --make-bed \
    --allow-no-sex \
    --out "${QC_FILTERED_PREFIX}"

echo "[3/6] Convert the QC dataset to VCF"
"${PLINK_BIN}" \
    --chr-set 39 \
    --allow-extra-chr \
    --bfile "${QC_FILTERED_PREFIX}" \
    --chr "${ANALYSIS_CHROMOSOMES}" \
    --recode vcf \
    --out "${QC_FILTERED_PREFIX}"

echo "[4/6] Sort, compress and index the VCF"
"${BCFTOOLS_BIN}" sort \
    "${QC_FILTERED_PREFIX}.vcf" \
    --output-file "${SORTED_VCF}" \
    --output-type z
"${BCFTOOLS_BIN}" index -f "${SORTED_VCF}"

echo "[5/6] Retain chromosomes ${ANALYSIS_CHROMOSOMES}"
"${BCFTOOLS_BIN}" view \
    --regions "${BCFTOOLS_REGIONS}" \
    "${SORTED_VCF}" \
    --output-type z \
    --output-file "${AUTOSOMAL_VCF}"
"${BCFTOOLS_BIN}" index -f "${AUTOSOMAL_VCF}"

echo "[6/6] Count variants by chromosome"
"${BCFTOOLS_BIN}" query -f '%CHROM\n' "${AUTOSOMAL_VCF}" \
    | sort \
    | uniq -c > "${CHR_COUNTS}"
awk '$1 < 3 {print $2}' "${CHR_COUNTS}" > "${EXCLUDE_CHR}"

echo "Genotype QC completed: ${AUTOSOMAL_VCF}"
