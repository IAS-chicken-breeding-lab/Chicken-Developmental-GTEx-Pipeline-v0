#!/bin/bash
tis=$1

# nominal cis-QTL mapping: nominal associations for all variant-phenotype pairs
echo 'nominal cis-QTL mapping'
python3 -m tensorqtl ${tis}/${tis}_geno ${tis}.bed.gz ${tis}/tensorqtl/nominal/${tis} \
    --covariates ${tis}_cov_.txt  --mode cis_nominal