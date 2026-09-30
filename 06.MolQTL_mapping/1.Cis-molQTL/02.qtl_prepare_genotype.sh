#!/bin/bash

tis=$1

## 1. Genotype quality control
plink \
  --chr-set 39 \
  --allow-extra-chr \
  --chr 1-39 \
  --bfile WGS_genome \
  --mind 0.1 \
  --maf 0.01 \
  --geno 0.1 \
  --keep-allele-order \
  --make-bed \
  --out WGS_genome_qc

## 2. Extract individuals for eQTL analysis
plink \
  --bfile WGS_genome_qc \
  --keep ${tis}-sampleID.txt \
  --keep-allele-order \
  --chr-set 39 \
  --allow-extra-chr \
  --make-bed \
  --chr 1-39 \
  --out ${tis}/${tis}_geno
