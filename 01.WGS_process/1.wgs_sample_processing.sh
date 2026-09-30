#!/bin/bash
set -euo pipefail

# Usage:
# bash 01_sample_processing.sh SAMPLE R1.fq.gz R2.fq.gz REF.fa OUTDIR THREADS

ID=$1
R1=$2
R2=$3
REF=$4
OUT=$5
THREADS=${6:-16}

mkdir -p "$OUT"/{fastq,bam,gvcf,qc,tmp}

# 1. Quality control
fastp \
  -i "$R1" \
  -I "$R2" \
  -o "$OUT/fastq/${ID}_R1.clean.fq.gz" \
  -O "$OUT/fastq/${ID}_R2.clean.fq.gz" \
  -q 20 \
  -u 50 \
  -n 5 \
  -l 100 \
  --detect_adapter_for_pe \
  -w "$THREADS" \
  -j "$OUT/qc/${ID}.fastp.json" \
  -h "$OUT/qc/${ID}.fastp.html"

# 2. Alignment and sorting
sentieon bwa mem \
  -t "$THREADS" \
  -R "@RG\tID:${ID}\tSM:${ID}\tLB:${ID}\tPL:DNBSEQ" \
  "$REF" \
  "$OUT/fastq/${ID}_R1.clean.fq.gz" \
  "$OUT/fastq/${ID}_R2.clean.fq.gz" |
sentieon util sort \
  -r "$REF" \
  -t "$THREADS" \
  --sam2bam \
  -i - \
  -o "$OUT/bam/${ID}.sorted.bam"

# 3. Duplicate marking
sentieon driver \
  -t "$THREADS" \
  -i "$OUT/bam/${ID}.sorted.bam" \
  --algo LocusCollector \
  --fun score_info \
  "$OUT/tmp/${ID}.score.txt"

sentieon driver \
  -t "$THREADS" \
  -i "$OUT/bam/${ID}.sorted.bam" \
  --algo Dedup \
  --score_info "$OUT/tmp/${ID}.score.txt" \
  --metrics "$OUT/qc/${ID}.dedup.metrics.txt" \
  "$OUT/bam/${ID}.dedup.bam"

# 4. Single-sample GVCF calling
sentieon driver \
  -r "$REF" \
  -t "$THREADS" \
  -i "$OUT/bam/${ID}.dedup.bam" \
  --algo Haplotyper \
  --emit_mode gvcf \
  "$OUT/gvcf/${ID}.g.vcf.gz"

# 5. Alignment QC
samtools flagstat \
  -@ "$THREADS" \
  "$OUT/bam/${ID}.dedup.bam" \
  > "$OUT/qc/${ID}.flagstat.txt"

samtools stats \
  -@ "$THREADS" \
  "$OUT/bam/${ID}.dedup.bam" \
  > "$OUT/qc/${ID}.stats.txt"

echo "${ID} finished."