#!/bin/bash
set -euo pipefail

# IBS-based sample swap detection
#
# Usage:
# bash IBS_sample_swap_detection.sh \
#   SAMPLE R1.fastq.gz R2.fastq.gz \
#   reference.fa STAR_index DNA_genotype_prefix output_dir threads

SAMPLE=$1
R1=$2
R2=$3
REF=$4
STAR_INDEX=$5
DNA_PREFIX=$6
OUT=$7
THREADS=${8:-16}

mkdir -p "$OUT"/{01_fastp,02_STAR,03_RNA_VCF,04_merge,05_IBS}


# ============================================================
# 1. RNA-seq read QC
# ============================================================

fastp \
  -i "$R1" \
  -I "$R2" \
  -o "$OUT/01_fastp/${SAMPLE}_clean_1.fq.gz" \
  -O "$OUT/01_fastp/${SAMPLE}_clean_2.fq.gz" \
  -h "$OUT/01_fastp/${SAMPLE}.html" \
  -j "$OUT/01_fastp/${SAMPLE}.json" \
  -w "$THREADS"


# ============================================================
# 2. RNA-seq alignment
# ============================================================

STAR \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn \
    "$OUT/01_fastp/${SAMPLE}_clean_1.fq.gz" \
    "$OUT/01_fastp/${SAMPLE}_clean_2.fq.gz" \
  --readFilesCommand zcat \
  --runThreadN "$THREADS" \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix "$OUT/02_STAR/${SAMPLE}_"

BAM="$OUT/02_STAR/${SAMPLE}_Aligned.sortedByCoord.out.bam"

samtools index "$BAM"


# ============================================================
# 3. RNA genotype calling
# ============================================================

bcftools mpileup \
  -f "$REF" \
  -Ou \
  "$BAM" |
bcftools call \
  -mv \
  -Oz \
  -o "$OUT/03_RNA_VCF/${SAMPLE}.vcf.gz"

bcftools index \
  "$OUT/03_RNA_VCF/${SAMPLE}.vcf.gz"


# Convert RNA VCF to PLINK
plink \
  --vcf "$OUT/03_RNA_VCF/${SAMPLE}.vcf.gz" \
  --make-bed \
  --allow-extra-chr \
  --keep-allele-order \
  --out "$OUT/04_merge/${SAMPLE}_RNA"


# ============================================================
# 4. Identify shared RNA-DNA SNPs
# ============================================================

cut -f2 "$OUT/04_merge/${SAMPLE}_RNA.bim" \
  | sort > "$OUT/04_merge/RNA.snps"

cut -f2 "${DNA_PREFIX}.bim" \
  | sort > "$OUT/04_merge/DNA.snps"

comm -12 \
  "$OUT/04_merge/RNA.snps" \
  "$OUT/04_merge/DNA.snps" \
  > "$OUT/04_merge/shared.snps"


plink \
  --bfile "$OUT/04_merge/${SAMPLE}_RNA" \
  --extract "$OUT/04_merge/shared.snps" \
  --make-bed \
  --allow-extra-chr \
  --keep-allele-order \
  --out "$OUT/04_merge/${SAMPLE}_RNA_shared"


plink \
  --bfile "$DNA_PREFIX" \
  --extract "$OUT/04_merge/shared.snps" \
  --make-bed \
  --allow-extra-chr \
  --keep-allele-order \
  --out "$OUT/04_merge/DNA_shared"


# ============================================================
# 5. Merge RNA and DNA genotypes
# ============================================================

plink \
  --bfile "$OUT/04_merge/DNA_shared" \
  --bmerge \
    "$OUT/04_merge/${SAMPLE}_RNA_shared.bed" \
    "$OUT/04_merge/${SAMPLE}_RNA_shared.bim" \
    "$OUT/04_merge/${SAMPLE}_RNA_shared.fam" \
  --make-bed \
  --allow-extra-chr \
  --keep-allele-order \
  --out "$OUT/04_merge/RNA_DNA_merged"


# ============================================================
# 6. IBS / identity analysis
# ============================================================

plink \
  --bfile "$OUT/04_merge/RNA_DNA_merged" \
  --genome \
  --allow-extra-chr \
  --keep-allele-order \
  --out "$OUT/05_IBS/IBS"


# IBS distance matrix
plink \
  --bfile "$OUT/04_merge/RNA_DNA_merged" \
  --distance square 1-ibs \
  --allow-extra-chr \
  --keep-allele-order \
  --out "$OUT/05_IBS/IBS_distance"


# Extract major identity statistics
awk 'NR==1 {
    print "FID1 IID1 FID2 IID2 Z0 Z1 Z2 PI_HAT DST"
}
NR>1 {
    print $1,$2,$3,$4,$7,$8,$9,$10,$11
}' "$OUT/05_IBS/IBS.genome" \
  > "$OUT/05_IBS/IBS_summary.txt"

echo "IBS-based sample identity analysis completed."