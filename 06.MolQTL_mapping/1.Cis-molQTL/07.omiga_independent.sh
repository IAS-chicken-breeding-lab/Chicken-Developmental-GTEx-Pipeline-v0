#!/bin/bash

tis=$1  ### Adipose
threads=$2

omiga --mode cis_independent \
  --genotype  ${tis}/${tis}_geno \
  --phenotype ${tis}.bed.gz \
  --covariates ${tis}_cov_.txt \
  --prefix ${tis} \
  --threads $threads \
  --cis-file ${tis}/${tis}_LMM/${tis}_LMM.cis_qtl.txt.gz \
  --qtl-map-model a+A \
  --output-dir ${tis}/${tis}_independent