#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 1 ]] || { echo "Usage: $0 CONFIG" >&2; exit 1; }
CONFIG=$1
[[ -f "$CONFIG" ]] || { echo "Missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

for tool in plink bedtools awk; do
  command -v "$tool" >/dev/null || { echo "Missing executable: $tool" >&2; exit 1; }
done

gene_map="$CROSSMAP_DIR/gene_mappability/gene_mappability.tsv"
snp_unique_bed="$MAPPABILITY_DIR/snp_mappability1.bed"
[[ -s "$gene_map" ]] || { echo "Missing gene mappability: $gene_map" >&2; exit 1; }
[[ -s "$snp_unique_bed" ]] || {
  echo "Missing $snp_unique_bed. Supply the project BED of uniquely mappable SNP positions." >&2
  exit 1
}
[[ -s "$REPEAT_BED" ]] || { echo "Missing repeat BED: $REPEAT_BED" >&2; exit 1; }

mkdir -p "$ANALYSIS_DIR/logs"
echo -e "Tissue\tN_SNPs\tN_Genes" > "$ANALYSIS_DIR/input_counts.tsv"
for tissue in "${TISSUES[@]}"; do
  out="$ANALYSIS_DIR/$tissue"
  mkdir -p "$out"/{geno,pheno,covar,logs,tmp}
  prefix="$GENOTYPE_DIR/$tissue"
  pheno="$PHENOTYPE_DIR/$tissue.bed"
  for ext in bed bim fam; do
    [[ -s "$prefix.$ext" ]] || { echo "Missing $prefix.$ext" >&2; exit 1; }
  done
  [[ -s "$pheno" ]] || { echo "Missing $pheno" >&2; exit 1; }

  awk 'BEGIN{OFS="\t"} {print $1,$4-1,$4,$2}' "$prefix.bim" > "$out/tmp/all_snps.bed"
  bedtools intersect -a "$out/tmp/all_snps.bed" -b "$snp_unique_bed" -wa -u | \
    bedtools intersect -a - -b "$REPEAT_BED" -wa -v > "$out/tmp/eligible_snps.bed"
  cut -f4 "$out/tmp/eligible_snps.bed" > "$out/tmp/eligible_snps.list"

  plink --bfile "$prefix" --chr-set "$AUTOSOME_NUM" \
    --extract "$out/tmp/eligible_snps.list" --make-bed --keep-allele-order \
    --out "$out/geno/$tissue" >"$out/logs/01_plink_filter.log" 2>&1

  awk -v threshold="$GENE_MAPPABILITY_THRESHOLD" '$2>=threshold {print $1}' \
    "$gene_map" > "$out/tmp/eligible_genes.list"
  awk 'NR==FNR {keep[$1]=1; next} FNR==1 || ($4 in keep)' \
    "$out/tmp/eligible_genes.list" "$pheno" > "$out/pheno/$tissue.bed"

  # Keep covariates in phenotype sample order. Base file: sample ID in col 1.
  head -1 "$BASE_COVARIATE_FILE" > "$out/covar/$tissue.base.tsv"
  head -1 "$out/pheno/$tissue.bed" | cut -f5- | tr '\t' '\n' > "$out/tmp/pheno_samples.list"
  awk 'NR==FNR {wanted[$1]=1; next} FNR>1 && ($1 in wanted)' \
    "$out/tmp/pheno_samples.list" "$BASE_COVARIATE_FILE" >> "$out/covar/$tissue.base.tsv"

  echo -e "$tissue\t$(wc -l < "$out/geno/$tissue.bim")\t$(($(wc -l < "$out/pheno/$tissue.bed")-1))" \
    >> "$ANALYSIS_DIR/input_counts.tsv"
done

echo "Prepared tissue inputs in $ANALYSIS_DIR"
