#!/bin/bash

tis=$1

mkdir -p aFC_result/${tis}

for chr in {1..39}; do
    python3 calculate_effect_sizes.py \
        --vcf ${tis}/${tis}_geno_vcf.gz \
        --pheno ${tis}.bed.gz \
        --cov ${tis}_cov_.txt \
        --qtl ${tis}/${tis}_LMM/${tis}_LMM.cis_qtl.txt.gz \
        --chr ${chr} \
        --log_xform 0 \
        --output aFC_result/${tis}/${tis}.Chr${chr}.log2aFC.txt \
        --boot 100
done

awk 'FNR > 1 || NR == 1' \
    aFC_result/${tis}/${tis}.Chr{1..39}.log2aFC.txt |
    gzip > aFC_result/${tis}/${tis}.log2aFC.txt.gz

