#!/bin/bash

tis=$1  ### Adipose
threads=$2

omiga --mode cis \
  --genotype  ${tis}/${tis}_geno \
  --phenotype ${tis}.bed.gz \
  --covariates ${tis}_cov_.txt \
  --prefix ${tis} \
  --threads $threads \
  --output-dir ${tis}/${tis}_LMM