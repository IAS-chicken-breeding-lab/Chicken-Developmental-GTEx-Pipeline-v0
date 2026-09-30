#!/bin/bash

CHR=$1
TIS=$2

GCTA="/home/cgbl/biosoft/gcta/gcta_1.93.2beta/gcta64"

CIS_GRM_DIR="gcta-h2/cis-h2/${CHR}_geno"
TRANS_GRM_DIR="gcta-h2/trans-h2/${CHR}"
PHENO_DIR="gcta-h2/cis-h2/${CHR}_gene_phenos"
QCOVAR_FILE="gcta-h2/cis-h2/${TIS}_qcovar.txt"
COVAR_FILE="gcta-h2/cis-h2/${TIS}_covar.txt"
OUTPUT_DIR="gcta-h2/h2-mgrm/${CHR}"

mkdir -p "${OUTPUT_DIR}"

for gene_dir in "${PHENO_DIR}"/*
do
    [ -d "${gene_dir}" ] || continue

    gene=$(basename "${gene_dir}")

    cis_grm="${CIS_GRM_DIR}/${gene}/${gene}_g1"
    trans_grm="${TRANS_GRM_DIR}/${gene}/${gene}_g1"
    pheno="${gene_dir}/${gene}_${CHR}.pheno.txt"
    mgrm="${OUTPUT_DIR}/${gene}_mgrm.txt"

    printf '%s\n%s\n' \
        "${trans_grm}" \
        "${cis_grm}" > "${mgrm}"

    "${GCTA}" \
        --reml \
        --mgrm "${mgrm}" \
        --pheno "${pheno}" \
        --qcovar "${QCOVAR_FILE}" \
        --covar "${COVAR_FILE}" \
        --thread-num 10 \
        --out "${OUTPUT_DIR}/${gene}-h2"
done
