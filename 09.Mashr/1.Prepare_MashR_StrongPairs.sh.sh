#!/bin/bash

tis_list=$1

dir_output="tissue-sharing-result"
base_dir="${tis}/${tis}_LMM"


mkdir -p "${dir_output}"

rm -f "${dir_output}/permutation_files.txt"

while read -r tis
do
    echo "${base_dir}/${tis}_LMM.cis_qtl.txt.gz" \
        >> "${dir_output}/permutation_files.txt"
done < "${tis_list}"

python3 combine_signif_pairs_tjy.py \
    "${dir_output}/permutation_files.txt" \
    strong_pairs \
    -o "${dir_output}"


while read -r tissue
do
{
    dir_nominal="${base_dir}"

    rm -f "${dir_output}/${tissue}.nominal_files.txt"

    for file in "${dir_nominal}/${tissue}_LMM.cis_qtl_pairs."*.txt.gz
    do
        echo "${file}" >> "${dir_output}/${tissue}.nominal_files.txt"
    done

    python3 extract_pairs_tjy.py \
        "${dir_output}/${tissue}.nominal_files.txt" \
        "${dir_output}/strong_pairs.combined_signifpairs.txt.gz" \
        "${tissue}" \
        -o "${dir_output}"
} &
done < "${tis_list}"

wait


rm -f "${dir_output}/strong_pairs_files.txt"

for file in "${dir_output}"/*.extracted_pairs.txt.gz
do
    echo "${file}" >> "${dir_output}/strong_pairs_files.txt"
done

python3 mashr_prepare_input.py \
    "${dir_output}/strong_pairs_files.txt" \
    strong_pairs \
    -o "${dir_output}" \
    --only_zscore

zcat "${dir_output}/strong_pairs.MashR_input.txt.gz" |
    sed 's/_zval//g' |
    gzip > "${dir_output}/strong_pairs.temp.txt.gz"

mv \
    "${dir_output}/strong_pairs.temp.txt.gz" \
    "${dir_output}/strong_pairs.MashR_input.txt.gz"