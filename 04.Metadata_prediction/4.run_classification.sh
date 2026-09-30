#!/bin/bash
set -euo pipefail

# Usage:
# bash run_classification.sh work_path group_method num_splits num_repeats \
#   feature_method feature_number classif_method class_label output_folder \
#   input_X input_y input_X_test

work_path=$1
group_method=$2
num_splits=$3
num_repeats=$4
feature_method=$5
feature_number=$6
classif_method=$7
class_label=$8
folder_path=$9
input_X_file=${10}
input_y_file=${11}
input_X_test_file=${12:-NA}

cd "$work_path"

python3 classification.py \
  "$work_path" \
  "$group_method" \
  "$num_splits" \
  "$num_repeats" \
  "$feature_method" \
  "$feature_number" \
  "$classif_method" \
  "$class_label" \
  "$folder_path" \
  "$input_X_file" \
  "$input_y_file" \
  "$input_X_test_file"
