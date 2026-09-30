#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 1 ]] || { echo "Usage: $0 CONFIG" >&2; exit 1; }
CONFIG=$1
[[ -f "$CONFIG" ]] || { echo "Missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

required_r=(gtf_to_txt.R compute_mappability.R generate_ambiguous_kmers.R compute_cross_mappability.R)
for script in "${required_r[@]}"; do
  [[ -f "$CROSSMAP_TOOL_DIR/$script" ]] || {
    echo "Missing $CROSSMAP_TOOL_DIR/$script" >&2
    echo "Clone https://github.com/battle-lab/crossmap and update CROSSMAP_TOOL_DIR." >&2
    exit 1
  }
done
for tool in Rscript bowtie awk; do
  command -v "$tool" >/dev/null || { echo "Missing executable: $tool" >&2; exit 1; }
done

annot_dir="$CROSSMAP_DIR/annotation"
gene_map_dir="$CROSSMAP_DIR/gene_mappability"
kmer_dir="$CROSSMAP_DIR/ambiguous_kmers"
align_dir="$CROSSMAP_DIR/ambiguous_kmer_alignments"
pair_dir="$CROSSMAP_DIR/cross_mappability_parts"
log_dir="$CROSSMAP_DIR/logs"
mkdir -p "$annot_dir" "$gene_map_dir" "$kmer_dir" "$align_dir" "$pair_dir" "$log_dir"

annot_tmp="$annot_dir/annot.exon_utr.unfiltered.tsv"
annot_native="$annot_dir/annot.exon_three_prime_utr.tsv"
annot="$annot_dir/annot.exon_UTR.tsv"
gene_map="$gene_map_dir/gene_mappability.tsv"

Rscript "$CROSSMAP_TOOL_DIR/gtf_to_txt.R" \
  -g "$GTF_FILE" -f "exon,three_prime_utr" -o "$annot_tmp" \
  >"$log_dir/01_gtf_to_txt.log" 2>&1

# Normalize annotation chromosome names to the exact names used by the FASTA.
awk -v keep="$(IFS=,; echo "${CHROMOSOMES[*]}")" -v prefix="$ANNOTATION_CHROM_PREFIX" '
  BEGIN{FS=OFS="\t"; n=split(keep,a,","); for(i=1;i<=n;i++) ok[a[i]]=1}
  NR==1 {print; next}
  {if(prefix!="" && index($2,prefix)==1) $2=substr($2,length(prefix)+1); if($2 in ok) print}
' "$annot_tmp" > "$annot_native"
sed 's/five_prime_utr/UTR/g; s/three_prime_utr/UTR/g' "$annot_native" > "$annot"

Rscript "$CROSSMAP_TOOL_DIR/compute_mappability.R" \
  --annot "$annot_native" \
  --k_exon "$EXON_K" --k_utr "$UTR_K" \
  --kmap_exon "$MAPPABILITY_DIR/chicken_${EXON_K}_${MISMATCHES}.bed" \
  --kmap_utr "$MAPPABILITY_DIR/chicken_${UTR_K}_${MISMATCHES}.bed" \
  --verbose 1 --output "$gene_map" \
  >"$log_dir/02_gene_mappability.log" 2>&1

Rscript "$CROSSMAP_TOOL_DIR/generate_ambiguous_kmers.R" \
  --mappability "$gene_map" --genome "$MAPPABILITY_DIR/chromosomes" \
  --annot "$annot" --k_exon "$EXON_K" --k_utr "$UTR_K" \
  --kmap_exon "$MAPPABILITY_DIR/chicken_${EXON_K}_${MISMATCHES}.bed" \
  --kmap_utr "$MAPPABILITY_DIR/chicken_${UTR_K}_${MISMATCHES}.bed" \
  --th1 0 --th2 1 --dir_name_len 12 --verbose 1 --o "$kmer_dir" \
  >"$log_dir/03_ambiguous_kmers.log" 2>&1

find "$kmer_dir" -type f -name '*.kmer.txt' -print0 | while IFS= read -r -d '' f; do
  awk '{print ">" NR-1; print}' "$f" > "${f%.txt}.fa"
done

common_args=(--annot "$annot" --mappability "$gene_map" --kmer "$kmer_dir" \
  --align "$align_dir" --index "$BOWTIE_INDEX_PREFIX" --mismatch "$MISMATCHES" \
  --max_chr "$AUTOSOME_NUM" --max_gene 200 --dir_name_len 12 --verbose 1 --o "$pair_dir")

Rscript "$CROSSMAP_TOOL_DIR/compute_cross_mappability.R" "${common_args[@]}" \
  --n1 1 --n2 "$CROSSMAP_BATCH_SIZE" --initonly TRUE \
  >"$log_dir/04_crossmap_init.log" 2>&1

n_genes=$(wc -l < "$gene_map")
for ((n1=1; n1<=n_genes; n1+=CROSSMAP_BATCH_SIZE)); do
  n2=$((n1 + CROSSMAP_BATCH_SIZE - 1))
  Rscript "$CROSSMAP_TOOL_DIR/compute_cross_mappability.R" "${common_args[@]}" \
    --n1 "$n1" --n2 "$n2" --initonly FALSE \
    >"$log_dir/05_crossmap_${n1}_${n2}.log" 2>&1
done

find "$pair_dir" -type f -name '*.crossmap.txt' -print0 | sort -z | \
  xargs -0 cat > "$CROSSMAP_DIR/cross_mappability.tsv"

echo "Gene mappability: $gene_map"
echo "Cross-mappability: $CROSSMAP_DIR/cross_mappability.tsv"
