#!/bin/bash
set -euo pipefail

# Usage:
# bash 02_joint_calling_qc.sh gvcf.list REF.fa OUTDIR THREADS

LIST=$1
REF=$2
OUT=$3
THREADS=${4:-16}

mkdir -p "$OUT"

# 1. Joint genotyping
sentieon driver \
  -r "$REF" \
  -t "$THREADS" \
  --algo GVCFtyper \
  "$OUT/cohort.raw.vcf.gz" \
  $(awk '{print "-v",$1}' "$LIST")

# 2. Select SNPs
gatk SelectVariants \
  -R "$REF" \
  -V "$OUT/cohort.raw.vcf.gz" \
  --select-type-to-include SNP \
  -O "$OUT/cohort.snp.vcf.gz"

# 3. Hard filtering
gatk VariantFiltration \
  -R "$REF" \
  -V "$OUT/cohort.snp.vcf.gz" \
  --filter-name "QD2" \
  --filter-expression "QD < 2.0" \
  --filter-name "FS60" \
  --filter-expression "FS > 60.0" \
  --filter-name "MQ40" \
  --filter-expression "MQ < 40.0" \
  --filter-name "MQRankSum" \
  --filter-expression "MQRankSum < -12.5" \
  --filter-name "ReadPosRankSum" \
  --filter-expression "ReadPosRankSum < -8.0" \
  --filter-name "SOR3" \
  --filter-expression "SOR > 3.0" \
  -O "$OUT/cohort.filtered.vcf.gz"

# 4. Retain PASS variants
gatk SelectVariants \
  -R "$REF" \
  -V "$OUT/cohort.filtered.vcf.gz" \
  --exclude-filtered true \
  -O "$OUT/cohort.pass.vcf.gz"

# 5. Individual missingness
vcftools \
  --gzvcf "$OUT/cohort.pass.vcf.gz" \
  --missing-indv \
  --out "$OUT/cohort"

awk 'NR > 1 && $5 > 0.10 {print $1}' \
  "$OUT/cohort.imiss" \
  > "$OUT/remove_individuals.txt"

# 6. Final SNP QC
vcftools \
  --gzvcf "$OUT/cohort.pass.vcf.gz" \
  --remove "$OUT/remove_individuals.txt" \
  --max-missing 0.90 \
  --maf 0.01 \
  --min-alleles 2 \
  --max-alleles 2 \
  --recode \
  --recode-INFO-all \
  --out "$OUT/cohort.final"

bgzip -f "$OUT/cohort.final.recode.vcf"
tabix -f -p vcf "$OUT/cohort.final.recode.vcf.gz"

echo "Final VCF: $OUT/cohort.final.recode.vcf.gz"