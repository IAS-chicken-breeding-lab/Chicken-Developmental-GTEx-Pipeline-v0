#!/bin/bash

tis=$1  ### Adipose
threads=$2

omiga --mode her_est \
  --genotype  ${tis}/${tis}_geno \
  --phenotype ${tis}.bed.gz \
  --covariates ${tis}_cov_.txt \
  --prefix ${tis} \
  --threads $threads \
  --h2-model At+Ac \
  --output-dir ${tis}/${tis}_heritability