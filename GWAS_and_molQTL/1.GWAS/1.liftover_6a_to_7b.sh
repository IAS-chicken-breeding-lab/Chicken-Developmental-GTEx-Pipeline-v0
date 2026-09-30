#!/usr/bin/env bash

set -euo pipefail

# Convert chicken genotype coordinates from GRCg6a to GRCg7b.
# Usage: bash scripts/01_liftover_6a_to_7b.sh [config/config.sh]

CONFIG_FILE="${1:-config/config.sh}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: configuration file not found: ${CONFIG_FILE}" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "${CONFIG_FILE}"

required_vars=(
    PLINK_BIN LIFTOVER_BIN CHAIN_FILE LIFTOVER_WORKDIR
    LIFTOVER_INPUT_PREFIX LIFTOVER_OUTPUT_PREFIX ANALYSIS_CHROMOSOMES
)

for var_name in "${required_vars[@]}"; do
    if [[ -z "${!var_name:-}" ]]; then
        echo "ERROR: ${var_name} is not set in ${CONFIG_FILE}" >&2
        exit 1
    fi
done

mkdir -p "${LIFTOVER_WORKDIR}"
cd "${LIFTOVER_WORKDIR}"

INPUT_PREFIX="${LIFTOVER_INPUT_PREFIX}"
RENAMED_PREFIX="${INPUT_PREFIX}_renamed"
FILTERED_PREFIX="${INPUT_PREFIX}_mapped"

INITIAL_BED="${INPUT_PREFIX}_positions.bed"
INITIAL_LIFTED_BED="${INPUT_PREFIX}_lifted_round1.bed"
INITIAL_UNMAPPED_BED="${INPUT_PREFIX}_unmapped_round1.bed"
CLEAN_UNMAPPED_BED="${INPUT_PREFIX}_unmapped_round1.clean.bed"
EXCLUDE_IDS="${INPUT_PREFIX}_exclude_ids.txt"

FILTERED_BED="${FILTERED_PREFIX}_positions.bed"
FINAL_LIFTED_BED="${LIFTOVER_OUTPUT_PREFIX}_positions.bed"
FINAL_UNMAPPED_BED="${INPUT_PREFIX}_unmapped_round2.bed"

echo "[1/8] Convert PED/MAP to PLINK binary format"
"${PLINK_BIN}" \
    --file "${INPUT_PREFIX}" \
    --chr-set 39 \
    --allow-extra-chr \
    --chr "${ANALYSIS_CHROMOSOMES}" \
    --keep-allele-order \
    --make-bed \
    --out "${INPUT_PREFIX}"

echo "[2/8] Prepare coordinates for the first liftOver pass"
# This preserves the original coordinate rule: start=BP and end=BP+1.
awk 'BEGIN {OFS="\t"} {print $1, $4, $4 + 1}' \
    "${INPUT_PREFIX}.map" > "${INITIAL_BED}"

echo "[3/8] Run the first liftOver pass"
"${LIFTOVER_BIN}" \
    "${INITIAL_BED}" \
    "${CHAIN_FILE}" \
    "${INITIAL_LIFTED_BED}" \
    "${INITIAL_UNMAPPED_BED}"

echo "[4/8] Build the list of unmapped SNP identifiers"
awk '!/^#/' "${INITIAL_UNMAPPED_BED}" > "${CLEAN_UNMAPPED_BED}"
awk 'BEGIN {OFS="\t"} {print $1 "_" $2}' \
    "${CLEAN_UNMAPPED_BED}" > "${EXCLUDE_IDS}"

echo "[5/8] Rename SNP identifiers and remove unmapped loci"
awk 'BEGIN {OFS="\t"} {print $1, $1 "_" $4, $3, $4}' \
    "${INPUT_PREFIX}.map" > "${RENAMED_PREFIX}.map"
cp "${INPUT_PREFIX}.ped" "${RENAMED_PREFIX}.ped"

"${PLINK_BIN}" \
    --file "${RENAMED_PREFIX}" \
    --exclude "${EXCLUDE_IDS}" \
    --make-bed \
    --out "${FILTERED_PREFIX}" \
    --chr-set 39 \
    --allow-extra-chr \
    --chr "${ANALYSIS_CHROMOSOMES}"

"${PLINK_BIN}" \
    --bfile "${FILTERED_PREFIX}" \
    --recode \
    --out "${FILTERED_PREFIX}" \
    --chr-set 39 \
    --allow-extra-chr \
    --chr "${ANALYSIS_CHROMOSOMES}"

echo "[6/8] Report chromosome distributions"
awk '{count[$1]++} END {for (chr in count) print "chromosome", chr, count[chr], "markers"}' \
    "${FILTERED_PREFIX}.map" | sort -V
awk '{count[$1]++} END {for (chr in count) print "chromosome", chr, count[chr], "markers"}' \
    "${FILTERED_PREFIX}.bim" | sort -V

echo "[7/8] Run the second liftOver pass"
awk 'BEGIN {OFS="\t"} {print $1, $4, $4 + 1}' \
    "${FILTERED_PREFIX}.map" > "${FILTERED_BED}"

"${LIFTOVER_BIN}" \
    "${FILTERED_BED}" \
    "${CHAIN_FILE}" \
    "${FINAL_LIFTED_BED}" \
    "${FINAL_UNMAPPED_BED}"

echo "[8/8] Generate the final PLINK dataset"
# The third MAP column is kept as in the original workflow.
awk 'BEGIN {OFS="\t"} {print $1, $1 "_" $2, $2, $2}' \
    "${FINAL_LIFTED_BED}" > "${LIFTOVER_OUTPUT_PREFIX}.map"
cp "${FILTERED_PREFIX}.ped" "${LIFTOVER_OUTPUT_PREFIX}.ped"

"${PLINK_BIN}" \
    --file "${LIFTOVER_OUTPUT_PREFIX}" \
    --make-bed \
    --out "${LIFTOVER_OUTPUT_PREFIX}" \
    --chr-set 39 \
    --allow-extra-chr \
    --chr "${ANALYSIS_CHROMOSOMES}"

echo "Coordinate conversion completed: ${LIFTOVER_WORKDIR}/${LIFTOVER_OUTPUT_PREFIX}"
