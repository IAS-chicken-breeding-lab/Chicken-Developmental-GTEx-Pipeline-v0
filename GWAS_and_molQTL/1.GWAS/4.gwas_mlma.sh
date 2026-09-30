#!/usr/bin/env bash

set -euo pipefail

# Construct a GRM, perform LD pruning and PCA, and run GCTA MLMA.
# Usage: bash scripts/05_gwas_mlma.sh [config/config.sh]

CONFIG_FILE="${1:-config/config.sh}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: configuration file not found: ${CONFIG_FILE}" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "${CONFIG_FILE}"

required_vars=(
    PLINK_BIN GCTA_BIN GWAS_WORKDIR GWAS_GENO_PREFIX GWAS_RAW_PHENO
    GWAS_CLEAN_PHENO GWAS_QCOVAR GWAS_GRM_PREFIX GWAS_PRUNE_PREFIX
    GWAS_PCA_PREFIX GWAS_RESULT_PREFIX N_TRAITS N_PCS THREADS
    ANALYSIS_CHROMOSOMES
)

for var_name in "${required_vars[@]}"; do
    if [[ -z "${!var_name:-}" ]]; then
        echo "ERROR: ${var_name} is not set in ${CONFIG_FILE}" >&2
        exit 1
    fi
done

mkdir -p "${GWAS_WORKDIR}"
cd "${GWAS_WORKDIR}"

echo "[1/4] Normalize phenotype delimiters and line endings"
# Replace visible ^I sequences with tabs and remove visible ^M at line ends.
# tr -d also removes actual carriage-return characters from CRLF files.
sed -e 's/\^I/\t/g' -e 's/\^M$//g' "${GWAS_RAW_PHENO}" \
    | tr -d '\r' > "${GWAS_CLEAN_PHENO}"

echo "[2/4] Construct the genomic relationship matrix"
"${GCTA_BIN}" \
    --bfile "${GWAS_GENO_PREFIX}" \
    --make-grm-alg 1 \
    --out "${GWAS_GRM_PREFIX}" \
    --autosome

echo "[3/4] Perform LD pruning and PCA"
"${PLINK_BIN}" \
    --allow-extra-chr \
    --chr-set 39 \
    --bfile "${GWAS_GENO_PREFIX}" \
    --chr "${ANALYSIS_CHROMOSOMES}" \
    --indep-pairwise 50 5 0.2 \
    --out "${GWAS_PRUNE_PREFIX}"

"${PLINK_BIN}" \
    --allow-extra-chr \
    --chr-set 39 \
    --bfile "${GWAS_GENO_PREFIX}" \
    --chr "${ANALYSIS_CHROMOSOMES}" \
    --extract "${GWAS_PRUNE_PREFIX}.prune.in" \
    --pca "${N_PCS}" header \
    --out "${GWAS_PCA_PREFIX}"

echo "[4/4] Run MLMA for ${N_TRAITS} phenotypes"
for ((trait_index = 1; trait_index <= N_TRAITS; trait_index++)); do
    echo "Running phenotype ${trait_index}/${N_TRAITS}"
    "${GCTA_BIN}" \
        --grm "${GWAS_GRM_PREFIX}" \
        --bfile "${GWAS_GENO_PREFIX}" \
        --pheno "${GWAS_CLEAN_PHENO}" \
        --mpheno "${trait_index}" \
        --qcovar "${GWAS_QCOVAR}" \
        --mlma \
        --thread-num "${THREADS}" \
        --out "${GWAS_RESULT_PREFIX}${trait_index}"
done

echo "GWAS completed in: ${GWAS_WORKDIR}"
