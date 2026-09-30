#!/bin/bash
tis=$1

# extract eGenes
zcat ${tis}/tensorqtl/nominal/${tis}/${tis}.cis_qtl.txt.gz |
  awk 'NR==1 || $18<0.05' \
  > ${tis}/tensorqtl/nominal/${tis}/${tis}.qval_lt0.05.txt
