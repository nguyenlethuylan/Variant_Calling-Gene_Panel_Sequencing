#!/bin/bash

while getopts "i:o:a:r:w:g:" opt; do
  case $opt in
    i) input=$OPTARG ;;
    o) output=$OPTARG ;;
    a) adapter=$OPTARG ;;
    r) references=$OPTARG ;;
    w) wetlab=$OPTARG ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
  esac
done

mkdir -p "$output/correct_base/coverage"
coverage="$output/correct_base/coverage"

for bam_file in "$output/bamfile/"*/*_sorted_RG_dedup.bam; do
  base_name=$(basename "$bam_file" _sorted_RG_dedup.bam)
  base_dir="$output/correct_base/$base_name"
  mkdir -p "$base_dir"

  # Create file coverage and target interval_list
  picard BedToIntervalList \
    I="$wetlab" \
    SD="$references/hg38.dict" \
    O="$coverage/target.interval_list"

  # For Exome sequencing or Gene panel Sequencing
  hs_metrics="$base_dir/${base_name}_collect_hs_metrics.txt"
  picard CollectHsMetrics \
    I="$bam_file" \
    O="$hs_metrics" \
    R="$references/hg38.fa" \
    BAIT_INTERVALS="$coverage/target.interval_list" \
    TARGET_INTERVALS="$coverage/target.interval_list"

  # Base Quality Score Recalibration (BQSR)
  recal_data="$base_dir/${base_name}_recal_data.table"
  gatk BaseRecalibrator \
    -I "$bam_file" \
    -R "$references/hg38.fa" \
    --known-sites "$references/All_20180418_fixed.vcf" \
    -O "$recal_data"

  gatk ApplyBQSR \
    -R "$references/hg38.fa" \
    -I "$bam_file" \
    --bqsr-recal-file "$recal_data" \
    -O "$base_dir/${base_name}_sorted_dedup_RG_recal.bam"
done
