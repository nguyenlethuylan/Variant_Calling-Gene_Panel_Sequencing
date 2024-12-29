#!/bin/bash

# Add parameter handling
while getopts "i:o:a:r:w" opt; do
  case $opt in
    i) input=$OPTARG ;;
    o) output=$OPTARG ;;
    a) adapter=$OPTARG ;;
    r) references=$OPTARG ;;
    w) wetlab=$OPTARG ;;  # Ensure this is defined
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
  esac
done

# Ensure output directory for assembly exists
mkdir -p "$output/assembly"

# Assembly step for each dataset
for R1_file in "$output/trimmed/"*_1_trimmed_paired.fq.gz; do
  base_name=$(basename "$R1_file" _1_trimmed_paired.fq.gz)
  sample="$output/assembly/$base_name"
  mkdir -p "$sample"

  # Align using bwa
  bwa mem "$references/hg38.fa" \
    "$output/trimmed/${base_name}_1_trimmed_paired.fq.gz" \
    "$output/trimmed/${base_name}_2_trimmed_paired.fq.gz" \
    -o "$sample/${base_name}.sam" \
    2> "$sample/bwa_${base_name}_check.log"

done