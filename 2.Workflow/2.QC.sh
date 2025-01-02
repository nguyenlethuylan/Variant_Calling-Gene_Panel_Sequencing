#!/bin/bash

# Add parameter handling
while getopts "i:o:a:r:w:" opt; do
  case $opt in
    i) input=$OPTARG ;;
    o) output=$OPTARG ;;
    a) adapter=$OPTARG ;;
    w) wetlab=$OPTARG ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1;;
  esac
done


# Ensure output directories exist
mkdir -p "$output/raw_fastqc"
mkdir -p "$output/trimmed"

# Read the input CSV file, skipping the header, and run FastQC and Trimmomatic
tail -n +2 "$input" | while IFS=, read -r sample read1 read2; do
  echo "Processing sample: $sample"
  echo "Read1 path: $read1"
  echo "Read2 path: $read2"

  # Check and process Read1
  if [ -f "$read1" ]; then
    echo "Running FastQC for Read1..."
    fastqc "$read1" -o "$output/raw_fastqc"
  else
    echo "Error: File $read1 does not exist, skipping sample $sample."
    continue
  fi

  # Check and process Read2
  if [ -f "$read2" ]; then
    echo "Running FastQC for Read2..."
    fastqc "$read2" -o "$output/raw_fastqc"
  else
    echo "Error: File $read2 does not exist, skipping sample $sample."
    continue
  fi

  # Run Trimmomatic
  echo "Running Trimmomatic for $sample..."
  trimmomatic PE \
    -phred33 \
    -threads 64 \
    -trimlog "$output/trimmed/${sample}_sum.log" \
    "$read1" \
    "$read2" \
    "$output/trimmed/${sample}_1_trimmed_paired.fq.gz" \
    "$output/trimmed/${sample}_1_trimmed_unpaired.fq.gz" \
    "$output/trimmed/${sample}_2_trimmed_paired.fq.gz" \
    "$output/trimmed/${sample}_2_trimmed_unpaired.fq.gz" \
    ILLUMINACLIP:"$adapter:2:30:10:8:3:true" \
    CROP:140


  # Check Trimmomatic output
  if [ -f "$output/trimmed/${sample}_1_trimmed_paired.fq.gz" ] && [ -f "$output/trimmed/${sample}_2_trimmed_paired.fq.gz" ]; then
    echo "Trimmomatic completed successfully for $sample."
  else
    echo "Error: Trimmomatic failed for $sample. Check the logs."
    continue
  fi
done

