#!/bin/bash

while getopts "i:o:a:r:w:" opt; do
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

# Create output directories
mkdir -p "$output/bamfile"
mkdir -p "$output/bamfile/bamTofastq"
mkdir -p "$output/bamfile/fastqc"

# BAM processing loop
for sam_file in "$output/assembly/"*/*.sam; do
  base_name=$(basename "$sam_file" .sam)
  base_folder=$(dirname "$sam_file")
  bam_dir="$output/bamfile/$base_name"
  mkdir -p "$bam_dir"

  # Convert SAM to BAM
  bam_file="$bam_dir/${base_name}.bam"
  samtools view -Sb "$sam_file" > "$bam_file"

  # Sort the BAM file
  sorted_bam="$bam_dir/${base_name}_sorted.bam"
  samtools sort "$bam_file" -o "$sorted_bam"

  # Index the sorted BAM file
  samtools index "$sorted_bam"

  # Add or replace read groups
  rg_bam="$bam_dir/${base_name}_sorted_RG.bam"
  picard AddOrReplaceReadGroups \
    I="$sorted_bam" \
    O="$rg_bam" \
    RGID=1 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM="$base_name"
    
  # Mark duplicates
  dedup_bam="$bam_dir/${base_name}_sorted_RG_dedup.bam"
  metrics_file="$bam_dir/${base_name}_sorted_RG_metrics.txt"
  picard MarkDuplicates \
    I="$rg_bam" \
    O="$dedup_bam" \
    M="$metrics_file"

  # Index the BAM with read groups
  samtools index "$dedup_bam"

  # Revert BAM to FASTQ
  bamTofastq=$output/bamfile/bamTofastq
  fastq_r1="$bamTofastq/${base_name}_sorted_RG_dedup_R1.fastq.gz"
  fastq_r2="$bamTofastq/${base_name}_sorted_RG_dedup_R2.fastq.gz"
  samtools fastq "$dedup_bam" \
    -1 >(gzip > "$fastq_r1") \
    -2 >(gzip > "$fastq_r2")

  # Run FastQC
  fastqc "$fastq_r1" "$fastq_r2" -o "$output/bamfile/fastqc/"

done

echo "Bam file processing complete."