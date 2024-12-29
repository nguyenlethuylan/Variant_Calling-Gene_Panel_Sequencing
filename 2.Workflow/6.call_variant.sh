#!/bin/bash

while getopts "i:o:a:r:w:" opt; do
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

# Create the variant directory if it doesn't exist
mkdir -p "$output/variant"

for variants in "$output/correct_base/"*/*_sorted_dedup_RG_recal.bam; do
  base_name=$(basename "$variants" _sorted_dedup_RG_recal.bam)
  variant_dir="$output/variant/$base_name"
  mkdir -p "$variant_dir"

  # Run GATK HaplotypeCaller
  gatk HaplotypeCaller \
    -R "$references/hg38.fa" \
    -I "$variants" \
    -O "$variant_dir/${base_name}_raw_variants.vcf" \
    2> "$variant_dir/gatk_check.log"

  # Run VariantFiltration on the generated VCF
  gatk VariantFiltration \
    -R "$references/hg38.fa" \
    -V "$variant_dir/${base_name}_raw_variants.vcf" \
    -O "$variant_dir/${base_name}_filtered_variants.vcf" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 4.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "FAIL" \
    2> "$variant_dir/${base_name}_filter_check.log"

  echo "Variant filtration completed for $base_name"
done
