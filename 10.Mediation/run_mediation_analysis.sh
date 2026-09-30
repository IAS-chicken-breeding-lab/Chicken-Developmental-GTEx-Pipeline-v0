#!/bin/bash
set -uo pipefail
shopt -s nullglob

# Usage:
# bash run_mediation_analysis.sh \
#   input_dir covariate_dir MediationAnalysis.R output_dir \
#   mediator_type [periods] [tissues]
#
# mediator_type: eigengene | tf | hormone | cell
# periods: "A B C"
# tissues: "" for all tissues, or "AF SZ XJ"

INPUT_DIR=$1
COVAR_DIR=$2
R_SCRIPT=$3
OUTPUT_DIR=$4
MEDIATOR_TYPE=$5
PERIOD_STRING=${6:-"A B C"}
TISSUE_FILTER=${7:-""}

if [[ ! "$MEDIATOR_TYPE" =~ ^(eigengene|tf|hormone|cell)$ ]]; then
    echo "Mediator type must be eigengene, tf, hormone or cell."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

LOG_FILE="${OUTPUT_DIR}/mediation_$(date +%F_%H-%M-%S).log"

read -r -a PERIODS <<< "$PERIOD_STRING"

for PERIOD in "${PERIODS[@]}"; do

    FILES=("${INPUT_DIR}/${PERIOD}_"*_Mediation.txt)

    for GENOTYPE_FILE in "${FILES[@]}"; do

        [[ -f "$GENOTYPE_FILE" ]] || continue

        TISSUE=$(basename "$GENOTYPE_FILE" |
            sed -E "s/${PERIOD}_(.*)_Mediation\.txt/\1/")

        # Optional tissue filter
        if [[ -n "$TISSUE_FILTER" ]]; then
            if [[ ! " $TISSUE_FILTER " =~ " $TISSUE " ]]; then
                continue
            fi
        fi

        COV_FILE="${COVAR_DIR}/${PERIOD}/${TISSUE}/${PERIOD}_${TISSUE}.auto_covar"
        OUTDIR="${OUTPUT_DIR}/${PERIOD}/${TISSUE}/${MEDIATOR_TYPE}"

        [[ -f "$COV_FILE" ]] || {
            echo "Missing covariate file: $COV_FILE" | tee -a "$LOG_FILE"
            continue
        }

        mkdir -p "$OUTDIR"

        # --------------------------------------------------
        # Parse input columns
        # --------------------------------------------------

        IFS=$'\t' read -r -a COLUMNS < "$GENOTYPE_FILE"

        HORMONE_NAMES=("GHRH" "INS" "GH" "cortisol" "IGF" "ghrelin")

        SEX_INDEX=-1
        HORMONE_INDEX=-1
        EIGENGENE_INDEX=-1

        for i in "${!COLUMNS[@]}"; do

            COL="${COLUMNS[$i]}"

            [[ "$COL" == "Sex" ]] && SEX_INDEX=$i

            if [[ " ${HORMONE_NAMES[*]} " =~ " ${COL} " ]]; then
                if [[ $HORMONE_INDEX -eq -1 || $i -lt $HORMONE_INDEX ]]; then
                    HORMONE_INDEX=$i
                fi
            fi

            if [[ "$COL" =~ ^M([1-9]|[1-4][0-9]|50)$ ]]; then
                if [[ $EIGENGENE_INDEX -eq -1 || $i -lt $EIGENGENE_INDEX ]]; then
                    EIGENGENE_INDEX=$i
                fi
            fi
        done

        CELLS=()
        TFS=()
        HORMONES=()
        EIGENGENES=()
        GENES=()
        SNPS=()

        # Cells: between Sex and hormones
        if [[ $SEX_INDEX -ge 0 &&
              $HORMONE_INDEX -gt $((SEX_INDEX + 1)) ]]; then

            for ((i=SEX_INDEX+1; i<HORMONE_INDEX; i++)); do
                CELLS+=("${COLUMNS[$i]}")
            done
        fi

        # TFs: between hormone block and eigengenes
        if [[ $HORMONE_INDEX -ge 0 &&
              $EIGENGENE_INDEX -gt $((HORMONE_INDEX + 6)) ]]; then

            for ((i=HORMONE_INDEX+6; i<EIGENGENE_INDEX; i++)); do
                TFS+=("${COLUMNS[$i]}")
            done
        fi

        for COL in "${COLUMNS[@]}"; do

            if [[ " ${HORMONE_NAMES[*]} " =~ " ${COL} " ]]; then
                HORMONES+=("$COL")

            elif [[ "$COL" =~ ^M([1-9]|[1-4][0-9]|50)$ ]]; then
                EIGENGENES+=("$COL")

            elif [[ "$COL" =~ ^[0-9] ]]; then
                SNPS+=("$COL")

            elif [[ "$COL" != "ID" &&
                    "$COL" != "Sex" &&
                    ! " ${CELLS[*]} " =~ " ${COL} " &&
                    ! " ${TFS[*]} " =~ " ${COL} " ]]; then

                GENES+=("$COL")
            fi
        done

        case "$MEDIATOR_TYPE" in
            cell)
                MEDIATORS=("${CELLS[@]}")
                ;;
            tf)
                MEDIATORS=("${TFS[@]}")
                ;;
            hormone)
                MEDIATORS=("${HORMONES[@]}")
                ;;
            eigengene)
                MEDIATORS=("${EIGENGENES[@]}")
                ;;
        esac

        PAIR_COUNT=${#GENES[@]}

        if [[ ${#SNPS[@]} -lt $PAIR_COUNT ]]; then
            PAIR_COUNT=${#SNPS[@]}
        fi

        echo "${PERIOD} ${TISSUE}: ${PAIR_COUNT} gene-SNP pairs; ${#MEDIATORS[@]} mediators" |
            tee -a "$LOG_FILE"

        # --------------------------------------------------
        # Mediation analysis
        # --------------------------------------------------

        for MEDIATOR in "${MEDIATORS[@]}"; do

            # Skip cell types with all-zero abundance
            if [[ "$MEDIATOR_TYPE" == "cell" ]]; then

                COL_NUM=-1

                for i in "${!COLUMNS[@]}"; do
                    if [[ "${COLUMNS[$i]}" == "$MEDIATOR" ]]; then
                        COL_NUM=$((i + 1))
                        break
                    fi
                done

                if [[ $COL_NUM -gt 0 ]]; then

                    NON_ZERO=$(awk -F'\t' -v c="$COL_NUM" '
                        NR > 1 && $c != 0 && $c != "0" {
                            print 1
                            exit
                        }
                    ' "$GENOTYPE_FILE")

                    [[ "$NON_ZERO" == "1" ]] || continue
                fi
            fi

            for ((i=0; i<PAIR_COUNT; i++)); do

                GENE="${GENES[$i]}"
                SNP="${SNPS[$i]}"

                RESULT="${OUTDIR}/${PERIOD}_${TISSUE}_mediation_results_${GENE}_${MEDIATOR}_${SNP}.txt"

                [[ -f "$RESULT" ]] && continue

                echo "Running: $PERIOD $TISSUE $GENE $MEDIATOR $SNP" |
                    tee -a "$LOG_FILE"

                Rscript "$R_SCRIPT" \
                    --period "$PERIOD" \
                    --tissue "$TISSUE" \
                    --gene "$GENE" \
                    --"$MEDIATOR_TYPE" "$MEDIATOR" \
                    --snp "$SNP" \
                    --genotype_file "$GENOTYPE_FILE" \
                    --cov_file "$COV_FILE" \
                    --output_dir "$OUTDIR" \
                    >> "$LOG_FILE" 2>&1
            done
        done
    done
done

echo "Mediation analysis completed." | tee -a "$LOG_FILE"