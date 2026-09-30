#!/bin/bash
tis=$1

# extract eGenes
zcat ${tis}/${tis}_LMM/${tis}_LMM.cis_qtl.txt.gz |
  awk 'NR==1 || $13<0.05' \
  > ${tis}/${tis}_LMM/${tis}.qval_g1_lt0.05.txt

