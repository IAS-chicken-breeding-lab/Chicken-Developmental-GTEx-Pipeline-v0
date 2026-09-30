#!/bin/bash
tis=$1

# cis-QTL finemapping
echo 'cis-QTL finemapping'
python3 -m tensorqtl ${tis}/${tis}_geno ${tis}.bed.gz ${tis} --covariates ${tis}_cov_.txt \
    --fdr 0.05  --qvalue_lambda 0.85 --output_dir ${tis}/tensorqtl/susie/${tis}  \
    --mode cis_susie  --cis_output ${tis}/tensorqtl/nominal/${tis}/${tis}.cis_qtl.txt.gz