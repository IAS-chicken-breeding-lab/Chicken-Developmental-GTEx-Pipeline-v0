#!/bin/bash

tis=$1  ### Adipose
threads=$2

omiga --genotype  ${tis}/${tis}_geno \
  --phenotype ${tis}.bed.gz \
  --dprop-pc-covar 0.001 \
  --prefix ${tis} \
  --output-dir ${tis}_cov_results \
  --rm-collinear-covar 0.95 \
  --threads $threads