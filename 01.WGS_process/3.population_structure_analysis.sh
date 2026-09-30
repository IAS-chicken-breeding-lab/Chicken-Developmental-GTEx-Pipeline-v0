#!/bin/bash
set -euo pipefail

# Population structure and genomic diversity analyses
# Usage:
# bash population_structure_analysis.sh WGSCommon filter.vcf output_dir threads

PLINK_PREFIX=$1
VCF=$2
OUT=$3
THREADS=${4:-32}

mkdir -p "$OUT"/{PCA,IBD,MAF,LD_decay}

# 1. PCA: LD pruning
plink \
  --bfile "$PLINK_PREFIX" \
  --indep-pairwise 25 5 0.2 \
  --allow-extra-chr \
  --chr 1-39 W Z \
  --chr-set 39 \
  --keep-allele-order \
  --out "$OUT/PCA/LDpruned"

plink \
  --bfile "$PLINK_PREFIX" \
  --extract "$OUT/PCA/LDpruned.prune.in" \
  --make-bed \
  --allow-extra-chr \
  --chr 1-39 W Z \
  --chr-set 39 \
  --keep-allele-order \
  --threads "$THREADS" \
  --out "$OUT/PCA/WGS_LDfilter"

plink \
  --bfile "$OUT/PCA/WGS_LDfilter" \
  --pca 20 \
  --allow-extra-chr \
  --chr 1-39 W Z \
  --chr-set 39 \
  --keep-allele-order \
  --threads "$THREADS" \
  --out "$OUT/PCA/WGS_PCA"

# 2. IBD / pairwise relatedness
plink \
  --bfile "$PLINK_PREFIX" \
  --genome \
  --allow-extra-chr \
  --chr 1-39 W Z \
  --chr-set 39 \
  --keep-allele-order \
  --threads "$THREADS" \
  --out "$OUT/IBD/IBD"

awk 'NR==1 {print "PI_HAT"} NR>1 {print $10}' \
  "$OUT/IBD/IBD.genome" \
  > "$OUT/IBD/PI_HAT.txt"

# 3. Minor allele frequency
plink \
  --bfile "$PLINK_PREFIX" \
  --freq \
  --allow-extra-chr \
  --chr 1-39 W Z \
  --chr-set 39 \
  --keep-allele-order \
  --out "$OUT/MAF/MAF_ALL"

# 4. LD decay
# Sample IDs must match those in the VCF header
awk '{print $1"_"$1}' "${PLINK_PREFIX}.fam" \
  > "$OUT/LD_decay/sample_names.txt"

PopLDdecay \
  -InVCF "$VCF" \
  -MaxDist 500k \
  -SubPop "$OUT/LD_decay/sample_names.txt" \
  -OutStat "$OUT/LD_decay/LD_decay"

echo "Population structure analysis completed."