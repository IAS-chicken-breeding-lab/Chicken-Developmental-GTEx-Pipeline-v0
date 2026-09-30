#!/bin/bash

strong_file="tissue-sharing-result/strong_pairs.MashR_input.txt.gz"
random_file="tissue_random_pairs/nominal_pairs.1000000_subset.MashR_input.txt.gz"

export OPENBLAS_NUM_THREADS=1

Rscript run_MashR.R ${strong_file} ${random_file} 0 ./output_top_paris

Rscript run_MashR.R ${strong_file} ${random_file} 1 ./output_top_paris_across_all