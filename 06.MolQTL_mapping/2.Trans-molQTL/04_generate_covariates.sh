#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 1 ]] || { echo "Usage: $0 CONFIG" >&2; exit 1; }
CONFIG=$1
[[ -f "$CONFIG" ]] || { echo "Missing config: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"
[[ -x "$OMIGA_BIN" ]] || { echo "OmiGA is not executable: $OMIGA_BIN" >&2; exit 1; }

for tissue in "${TISSUES[@]}"; do
  out="$ANALYSIS_DIR/$tissue"
  "$OMIGA_BIN" --genotype "$out/geno/$tissue" --phenotype "$out/pheno/$tissue.bed" \
    --dprop-pc-covar 0.001 --prefix "$tissue" --output-dir "$out/covar" \
    >"$out/logs/02_omiga_pca.log" 2>&1

  auto="$out/covar/$tissue.auto_covar"
  base="$out/covar/$tissue.base.tsv"
  combined="$out/covar/${tissue}_combined_covar.txt"
  [[ -s "$auto" && -s "$base" ]] || { echo "Missing covariate input for $tissue" >&2; exit 1; }

  awk 'BEGIN{FS=OFS="\t"}
    NR==FNR {if(NR>1){sex[$1]=$2; period[$1]=$3}; next}
    FNR==1 {for(i=2;i<=NF;i++) id[i]=$i}
    {print}
    END {
      printf "Sex"; for(i=2;i<=length(id)+1;i++) printf OFS ((id[i] in sex) ? sex[id[i]] : "NA"); print ""
      printf "Period"; for(i=2;i<=length(id)+1;i++) printf OFS ((id[i] in period) ? period[id[i]] : "NA"); print ""
    }' "$base" "$auto" > "$combined"

  if grep -Eq $'\tNA(\t|$)' "$combined"; then
    echo "Unmatched covariate sample in $tissue; inspect $combined" >&2
    exit 1
  fi
done

echo "Combined PCA, Sex and Period covariates written under $ANALYSIS_DIR"
