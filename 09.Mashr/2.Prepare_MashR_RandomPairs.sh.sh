#!/bin/bash

tis_list=$1
subset_size=1000000

dir_output="tissue_random_pairs"

mkdir -p "${dir_output}"

rm -f "${dir_output}/nominal_combined_files.txt"

while read -r tis
do
    base_dir="${tis}/${tis}_LMM"

    for file in "${base_dir}/${tis}_LMM.cis_qtl_pairs."*.txt.gz
    do
        echo "${file}" >> "${dir_output}/nominal_combined_files.txt"
    done
done < "${tis_list}"

python3 combine_signif_pairs_tjy.py \
    "${dir_output}/nominal_combined_files.txt" \
    nominal_pairs \
    -o "${dir_output}"


while read -r tissue
do
{
    base_dir="${tissue}/${tissue}_LMM"

    rm -f "${dir_output}/${tissue}.nominal_files.txt"

    for file in "${base_dir}/${tissue}_LMM.cis_qtl_pairs."*.txt.gz
    do
        echo "${file}" >> "${dir_output}/${tissue}.nominal_files.txt"
    done

    python3 extract_pairs_tjy.py \
        "${dir_output}/${tissue}.nominal_files.txt" \
        "${dir_output}/nominal_pairs.combined_signifpairs.txt.gz" \
        "${tissue}_nominal_pairs" \
        -o "${dir_output}"
} &
done < "${tis_list}"

wait


rm -f "${dir_output}/nominal_pairs_files.txt"

for file in "${dir_output}"/*_nominal_pairs.extracted_pairs.txt.gz
do
    echo "${file}" >> "${dir_output}/nominal_pairs_files.txt"
done

base_dir=$(head -n 1 "${tis_list}")
base_dir="${base_dir}/${base_dir}_LMM"

python3 mashr_prepare_input.py \
    "${dir_output}/nominal_pairs_files.txt" \
    "nominal_pairs.${subset_size}_subset" \
    -o "${dir_output}" \
    --only_zscore \
    --dropna \
    --subset "${subset_size}" \
    --seed 9823

zcat "${dir_output}/nominal_pairs.${subset_size}_subset.MashR_input.txt.gz" |
    sed 's/.nominal_pairs_zval//g' |
    gzip > "${dir_output}/nominal_pairs.${subset_size}_subset.temp.txt.gz"

mv \
    "${dir_output}/nominal_pairs.${subset_size}_subset.temp.txt.gz" \
    "${dir_output}/nominal_pairs.${subset_size}_subset.MashR_input.txt.gz"
